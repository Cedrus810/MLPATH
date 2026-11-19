"""Measurement (2026-11-15, user 「跑」): why do the fixed-rho transverse relaxations of the
bracket scan converge so rarely on b05? (PBS 28756: 24/388 and 39/724 scan points reached
transverse fmax <= 0.003 eV/A within ef.max_relax_steps = 200; every bracket no_bracket.)

Question: step cap (still decreasing at 200) or floor (fmax stalls above 0.003)?
Quantity: projected (transverse) fmax in eV/A per FIRE step, vs the threshold
ef.transverse_fmax_eV_A (0.003) and sigma_lambda of HF_V15 (0.0042, reaction-force jitter).
Also: does the sign of the reaction force at step 200 equal its sign at the end?

Re-generates the SC items of an existing run (its prepare.json copied in, SC nodes
re-minimized here), replays the chosen paths with run_paths, takes the same bracket frames
(bracket_frames on the first MAX_MAXIMA maxima) and relaxes each at its own rho with the
scan's settings (FIRE, maxstep 0.05, Newton restore) but --max-steps (default 2000),
logging fmax_projected and the reaction force every step. Not resumable (short).

    python bracket_convergence.py --run DIR --system b05 --seed S --shell R --sc-poses 32 \
        --ks 9,39,93 --out OUTDIR [--max-steps 2000]
"""

import argparse
import json
import shutil
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from run_mc_afir import MAX_MAXIMA, prepare  # noqa: E402  (sets sys.path)
from run_mc_afir_batch import sc_items  # noqa: E402
from run_d2 import make_atoms, make_config  # noqa: E402
from systems import HF_V15, HS_BRACKET_V15, SYSTEMS  # noqa: E402

from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.afir_batch import _evaluate, _Fire, bracket_frames, run_paths  # noqa: E402
from explore.afir_path import PathSettings  # noqa: E402
from explore.d2 import _ef_settings  # noqa: E402
from explore.ledger import Ledger  # noqa: E402
from explore.parallel import MACEBatchEvaluator  # noqa: E402
from explore.poses import sample_poses  # noqa: E402


def relax_logged(template, starts, cons, evaluate, ledger, *, fmax, max_steps, maxstep=0.05,
                 tol=1e-8, newton_iter=50):  # fmt: skip
    """afir_batch.minimize (constrained branch) with a per-step log of fmax_projected and
    the reaction force; a row leaves the batch at fmax_projected <= fmax."""
    x = np.stack([np.asarray(s, float) for s in starts])
    n = len(x)

    def restore(r):
        coord, q = cons[r]
        for _ in range(newton_iter):
            err = coord.value(x[r]) - q
            if abs(err) <= tol:
                return
            g = coord.gradient(x[r])
            x[r] -= err * g / np.vdot(g, g)

    for r in range(n):
        restore(r)
    fire = _Fire(n, x.shape[1:], maxstep)
    hist = [{"fmax": [], "rf": [], "E": []} for _ in range(n)]
    done_at = [None] * n
    active = list(range(n))
    with ledger.segment(task="measure", stage="relax_logged", attempt=0) as seg:
        for step in range(max_steps + 1):
            if not active:
                break
            e, f = _evaluate(evaluate, template, x[active], ledger)
            keep, fp_all = [], []
            for idx, r in enumerate(active):
                g = cons[r][0].gradient(x[r])
                rf = -float(np.vdot(f[idx], g) / np.vdot(g, g))
                fp = f[idx] + rf * g
                fm = float(np.linalg.norm(fp, axis=1).max())
                hist[r]["fmax"].append(fm)
                hist[r]["rf"].append(rf)
                hist[r]["E"].append(float(e[idx]))
                if fm <= fmax:
                    done_at[r] = step
                elif step < max_steps:
                    keep.append(r)
                    fp_all.append(fp)
            if keep:
                x[keep] += fire.step(keep, np.stack(fp_all))
                for r in keep:
                    restore(r)
            active = keep
        seg.outcome = f"{n} relaxations"
    return hist, done_at


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--run", required=True, help="existing SC run dir (prepare.json, wave_*.json)"
    )
    ap.add_argument("--system", required=True)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--shell", type=float, required=True)
    ap.add_argument("--sc-poses", type=int, default=32)
    ap.add_argument("--max-paths", type=int, default=4096)
    ap.add_argument("--gamma", type=float, default=106.9)
    ap.add_argument("--ks", required=True)
    ap.add_argument("--max-steps", type=int, default=2000)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    run, out = Path(args.run), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    shutil.copy(run / "prepare.json", out / "prepare.json")
    ks = [int(k) for k in args.ks.split(",")]
    recorded = {}
    for f in sorted(run.glob("wave_*.json")):
        for r in json.loads(f.read_text())["paths"]:
            if r["k"] in ks:
                recorded[r["k"]] = r

    t0 = time.monotonic()
    system = SYSTEMS[args.system]
    config = make_config(system)
    calc = polar1("cuda", "float32")
    atoms, _ = make_atoms(system, config)
    probe = atoms.copy()
    probe.calc = calc
    noise = 2 * float(np.spacing(np.float32(abs(float(probe.get_potential_energy())))))
    frags, confs, _ = prepare(out, system, config, atoms, calc)
    base = PathSettings(noise_eV=noise, gamma_kJ=args.gamma)
    settings = PathSettings(noise_eV=noise, gamma_kJ=args.gamma,
                            early_discard_steps=base.max_steps + 1, stop_after_maxima=2)  # fmt: skip
    ef = _ef_settings(HF_V15, HS_BRACKET_V15)
    poses = [p for p in sample_poses(atoms, frags, confs[0], confs[1], [args.shell],
                                     args.max_paths, args.seed) if p.status == "placed"]  # fmt: skip
    batch = MACEBatchEvaluator(calc)
    items, _ = sc_items(atoms, poses, config, batch, out,
                        SimpleNamespace(seed=args.seed, sc_poses=args.sc_poses, sc_coord="rho"))  # fmt: skip
    chosen = [items[k] for k in ks]
    for it in chosen:
        if it["label"] != recorded[it["k"]]["label"]:
            raise RuntimeError(
                f"k {it['k']}: label {it['label']} != {recorded[it['k']]['label']}"
            )

    ledger = Ledger(out / "measure.ledger.jsonl", arm="AFIR", replay=False)
    paths = run_paths(atoms, [it["start"] for it in chosen], [it["coord"] for it in chosen],
                      batch, config, settings, ledger, "measure|paths")  # fmt: skip
    rows, starts, cons = [], [], []
    replay = []
    for it, p in zip(chosen, paths):
        rec = recorded[it["k"]]
        replay.append({"k": it["k"], "label": it["label"], "status": p.status,
                       "recorded_status": rec["status"], "steps": len(p.frames),
                       "recorded_steps": rec["steps"], "maxima": p.maxima,
                       "recorded_maxima": rec["maxima"],
                       "bond_change": p.frames[-1]["edges"] != it["ref_edges"]})  # fmt: skip
        for j, top in enumerate(p.maxima[:MAX_MAXIMA]):
            back, fwd = bracket_frames(p.frames, top)
            for side, idx in (("back", back), ("fwd", fwd)):
                for pos, fi in enumerate(idx):
                    x = np.asarray(p.frames[fi]["positions"], float)
                    rows.append({"k": it["k"], "top": j, "side": side, "pos": pos, "frame": fi})
                    starts.append(x)
                    cons.append((it["coord"], it["coord"].value(x)))
    print(json.dumps(replay), flush=True)
    print(f"{len(rows)} scan points, max_steps {args.max_steps}", flush=True)
    t1 = time.monotonic()
    hist, done_at = relax_logged(atoms, starts, cons, batch, ledger,
                                 fmax=ef.transverse_fmax_eV_A, max_steps=args.max_steps)  # fmt: skip
    ledger.close()
    cap = ef.max_relax_steps
    pts = []
    for row, h, d, c in zip(rows, hist, done_at, cons):
        fm, rf = h["fmax"], h["rf"]
        at_cap = min(cap, len(fm) - 1)
        pts.append(row | {
            "q": c[1], "done_at": d, "n_steps": len(fm) - 1,
            "fmax_at_cap": fm[at_cap], "fmax_final": fm[-1], "fmax_min": min(fm),
            "rf_at_cap": rf[at_cap], "rf_final": rf[-1],
            "sign_same": bool(np.sign(rf[at_cap]) == np.sign(rf[-1])),
            "fmax_trace": fm[::10], "rf_trace": rf[::10],
        })  # fmt: skip
    dn = [p["done_at"] for p in pts if p["done_at"] is not None]
    summary = {
        "points": len(pts), "threshold": ef.transverse_fmax_eV_A, "cap": cap,
        "sigma_lambda": HF_V15.reaction_force_sigma, "max_steps": args.max_steps,
        "converged_by_cap": sum(d <= cap for d in dn), "converged_by_max": len(dn),
        "done_at_quantiles": np.quantile(dn, [0.1, 0.5, 0.9]).tolist() if dn else None,
        "unconverged_fmax_final": sorted(round(p["fmax_final"], 4) for p in pts
                                         if p["done_at"] is None),
        "sign_changed_after_cap": sum(not p["sign_same"] for p in pts),
        "forwards": ledger.spent, "relax_s": time.monotonic() - t1,
        "wall_s": time.monotonic() - t0,
    }  # fmt: skip
    (out / "points.json").write_text(json.dumps({"replay": replay, "points": pts}) + "\n")
    (out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
