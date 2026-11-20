"""E2 of docs/TS_LOCATOR_DESIGN_2026-11-19.md (2026-11-19, user「可以」to the E2 plan, then
"你先冒烟然后给我丢上来"): short string + fixed-direction climbing image on every path of
a finished SC-AFIR run whose HF kernel attempted a TS (path entries with "ts" in
wave_*.json / dwave_*.json), path by path against the HF kernel's own outcome.

  1. labels and HF statuses from RUN's waves; seed / shell / sc_poses / max_paths / gamma
     from RUN/manifest.json;
  2. replay of those paths exactly as n1_string_check.py (prepare.json copied, SC nodes
     re-minimized, only the named items built, pair / two-pair push, no early discard);
     frames cached in OUT/frames/<label>.npz and reused on a rerun;
  3. paths whose replay changed the bond graph: end minima (batched over all paths, fmax
     0.01), 20 nodes at equal arc length from the early low-energy frame, simplified string
     (batched over all paths), snapshots at --checkpoints (default 0, 2);
  4. from every snapshot: climbing image from the highest interior node along the fixed
     tangent x[i+1] - x[i-1] (n1_tangent_variants.climb_batch; fmax 0.01, ceiling +3 eV),
     all rows in one batch; converged images descend from x +/- 0.05 A tau (batched) and
     each side is labelled reactant / product / other against the two end minima;
  5. per chain: number of interior energy maxima with an uphill-to-downhill tangential force
     change (f_i < 0 < f_{i+1}), i.e. how many peaks a multi-peak path shows (no extra
     forward).
Writes OUT/e2.json, OUT/chains.npz, OUT/k<k>/wave_0000.json (converged images as
converged_candidate, for afir_validate.py). Resumable only through the frame cache.

    python e2_string_climb.py --run RUNDIR --system b04 --out OUT [--checkpoints 0,2]
"""

import argparse
import glob
import json
import shutil
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from run_mc_afir import prepare  # noqa: E402  (sets sys.path)
from run_mc_afir_batch import sc_items  # noqa: E402
from run_d2 import make_atoms, make_config  # noqa: E402
from systems import SYSTEMS  # noqa: E402
from n1_tangent_variants import climb_batch  # noqa: E402

from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.afir_batch import _evaluate, minimize, run_paths  # noqa: E402
from explore.afir_path import PathSettings  # noqa: E402
from explore.afir_string import reparametrize, resample, tangent  # noqa: E402
from explore.ledger import Ledger  # noqa: E402
from explore.parallel import MACEBatchEvaluator  # noqa: E402
from explore.poses import _jsonable, sample_poses  # noqa: E402


def hf_attempts(run):
    out = {}
    for f in sorted(set(glob.glob(str(run / "wave_*.json")))) + sorted(
        set(glob.glob(str(run / "dwave_*.json")))
    ):
        for p in json.load(open(f))["paths"]:
            if p.get("ts"):
                out[p["label"]] = {
                    "double": bool(p.get("double")),
                    "hf_status": [t.get("status") for t in p["ts"]],
                    "orig_first_energy_eV": (p.get("energies_eV") or [None])[0],
                }
    return out


def n_peaks(x, energies, forces):
    f = [np.vdot(forces[i], tangent(x, i)) for i in range(len(x))]
    c = 0
    for i in range(1, len(x) - 1):
        if energies[i] >= energies[i - 1] and energies[i] >= energies[i + 1]:
            if f[i - 1] < 0 < f[i + 1]:
                c += 1
    return c


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--system", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--checkpoints", default="0,2")
    ap.add_argument("--nodes", type=int, default=20)
    ap.add_argument("--fmax", type=float, default=0.01)
    ap.add_argument("--ci-steps", type=int, default=2000)
    ap.add_argument("--side-step", type=float, default=0.05)
    ap.add_argument("--limit", type=int, default=0, help="first N labels only (smoke)")
    ap.add_argument("--labels", default="", help="comma-separated subset of the labels")
    args = ap.parse_args()
    run, out = Path(args.run), Path(args.out)
    (out / "frames").mkdir(parents=True, exist_ok=True)
    shutil.copy(run / "prepare.json", out / "prepare.json")
    margs = json.load(open(run / "manifest.json"))["args"]
    hf = hf_attempts(run)
    labels = sorted(hf)
    if args.labels:
        want = set(args.labels.split(","))
        labels = [lb for lb in labels if lb in want]
    if args.limit:
        labels = labels[: args.limit]
    ks = sorted(int(k) for k in args.checkpoints.split(","))

    t0 = time.monotonic()
    system = SYSTEMS[args.system]
    config = make_config(system)
    calc = polar1("cuda", "float32")
    atoms, _ = make_atoms(system, config)
    batch = MACEBatchEvaluator(calc)
    safe = {lb: lb.replace("|", "_").replace("+", "_") for lb in labels}
    cache = {lb: out / "frames" / f"{safe[lb]}.npz" for lb in labels}
    frames_of, replay_info = {}, {}
    todo = [lb for lb in labels if not cache[lb].exists()]
    for lb in labels:
        if lb not in todo:
            z = np.load(cache[lb])
            frames_of[lb] = (z["positions"], z["energies"])
            replay_info[lb] = json.loads(str(z["info"]))
    if todo:
        probe = atoms.copy()
        probe.calc = calc
        noise = 2 * float(np.spacing(np.float32(abs(float(probe.get_potential_energy())))))
        frags, confs, _ = prepare(out, system, config, atoms, calc)
        settings = PathSettings(
            noise_eV=noise, gamma_kJ=margs["gamma"], early_discard_steps=1001
        )
        poses = [
            p
            for p in sample_poses(atoms, frags, confs[0], confs[1], [margs["shell"]],
                                  margs["max_paths"], margs["seed"])
            if p.status == "placed"
        ]  # fmt: skip
        doubles = []
        singles, _ = sc_items(
            atoms, poses, config, batch, out,
            SimpleNamespace(seed=margs["seed"], sc_poses=margs["sc_poses"], sc_coord="pair",
                            sc_double=1),
            doubles,
            double_filter=lambda lb: lb in todo,
        )  # fmt: skip
        by_label = {it["label"]: it for it in singles + doubles}
        # SC nodes are re-minimized in fp32 and deduplicated, so node numbering can differ
        # from the original run (b05 sA, 2026-11-19): unmatched labels are recorded, not run.
        missing = [lb for lb in todo if lb not in by_label]
        for lb in missing:
            replay_info[lb] = {"status": "not_rebuilt", "bond_change": False}
        todo = [lb for lb in todo if lb in by_label]
        chosen = [by_label[lb] for lb in todo]
        led = Ledger(out / "replay.ledger.jsonl", arm="E2", replay=False)
        paths = run_paths(atoms, [it["start"] for it in chosen], [it["coord"] for it in chosen],
                          batch, config, settings, led, "e2|replay")  # fmt: skip
        led.close()
        for lb, it, p in zip(todo, chosen, paths):
            fr = p.frames
            info = {"status": p.status, "steps": len(fr), "forwards": p.spent,
                    "bond_change": bool(fr) and fr[-1]["edges"] != it["ref_edges"],
                    "first_energy_eV": float(fr[0]["energy_eV"]) if fr else None}  # fmt: skip
            pos = np.stack([f["positions"] for f in fr])
            en = np.array([f["energy_eV"] for f in fr])
            np.savez(cache[lb], positions=pos, energies=en, info=json.dumps(info))
            frames_of[lb], replay_info[lb] = (pos, en), info
    t_replay = time.monotonic() - t0
    print(json.dumps({"replay_s": t_replay, "labels": len(labels),
                      "bond_change": sum(replay_info[lb]["bond_change"] for lb in labels)}),
          flush=True)  # fmt: skip

    ledger = Ledger(out / "e2.ledger.jsonl", arm="E2", replay=False)
    use = [lb for lb in labels if replay_info[lb]["bond_change"]]
    rec = {lb: {"label": lb, **hf[lb], "replay": replay_info[lb]} for lb in labels}
    for lb in labels:  # identity of the replayed path: start energy vs the original run's
        a, b = replay_info[lb].get("first_energy_eV"), hf[lb]["orig_first_energy_eV"]
        rec[lb]["first_energy_diff_eV"] = None if a is None or b is None else a - b
    for lb in labels:
        if lb not in use:
            rec[lb]["e2"] = "skipped: replay did not change the bond graph"
    if use:
        starts = []
        for lb in use:
            pos, _ = frames_of[lb]
            starts += [pos[0], pos[-1]]
        ends = minimize(atoms, starts, batch, config, ledger, "e2|ends",
                        fmax=args.fmax, max_steps=2000)  # fmt: skip
        xs = []
        for j, lb in enumerate(use):
            pos, en = frames_of[lb]
            low = int(np.argmin(en[: len(en) // 3 + 1]))
            mids = resample(list(pos[low:]), args.nodes)[1:-1]
            e_r, e_p = ends[2 * j], ends[2 * j + 1]
            xs.append(reparametrize(np.stack([e_r.positions] + list(mids) + [e_p.positions])))
            rec[lb].update(
                ends_converged=[bool(e_r.converged), bool(e_p.converged)],
                reaction_energy_eV=e_p.energy_eV - e_r.energy_eV,
                cost={"ends": e_r.forwards + e_p.forwards},
            )
        x = np.stack(xs)  # (P, N, n, 3)
        n_in = x.shape[1] - 2
        snaps = {}
        with ledger.segment(task="e2|string", stage="string", attempt=0) as seg:
            for it in range(0, ks[-1] + 1):
                inner = x[:, 1:-1].reshape((-1,) + x.shape[2:])
                e, f = _evaluate(batch, atoms, inner, ledger)
                e = np.asarray(e).reshape(len(use), n_in)
                f = np.asarray(f).reshape((len(use), n_in) + x.shape[2:])
                if it in ks:
                    snaps[it] = (x.copy(), e.copy(), f.copy())
                if it == ks[-1]:
                    break
                step = 0.05 * f
                size = np.linalg.norm(step, axis=3).max(axis=2)
                step *= np.minimum(1.0, 0.05 / np.where(size > 0, size, 1.0))[:, :, None, None]
                moved = x.copy()
                moved[:, 1:-1] += step
                x = np.stack([reparametrize(m) for m in moved])
            seg.outcome = f"{ks[-1]} iterations"
        rows, cstart, ctau = [], [], []
        chain_e, chain_x = {}, {}
        for j, lb in enumerate(use):
            e_r, e_p = ends[2 * j], ends[2 * j + 1]
            for k in ks:
                sx, se, sf = snaps[k]
                energies = np.concatenate([[e_r.energy_eV], se[j], [e_p.energy_eV]])
                forces = np.concatenate(
                    [np.zeros_like(sf[j][:1]), sf[j], np.zeros_like(sf[j][:1])]
                )
                top = int(np.argmax(energies[1:-1])) + 1
                rows.append((lb, k, top, n_peaks(sx[j], energies, forces)))
                cstart.append(sx[j][top])
                ctau.append(tangent(sx[j], top))
                chain_x[f"{safe[lb]}_k{k}"] = sx[j]
                chain_e[f"{safe[lb]}_k{k}"] = energies
        np.savez(out / "chains.npz", **{f"x_{q}": v for q, v in chain_x.items()},
                 **{f"e_{q}": v for q, v in chain_e.items()})  # fmt: skip
        ceil = []
        for lb, _, _, _ in rows:
            j = use.index(lb)
            ceil.append(ends[2 * j].energy_eV + 3.0)
        cis = climb_batch(atoms, cstart, ctau, batch, ledger, "e2|climb", ci_fmax=args.fmax,
                          ci_steps=args.ci_steps, ceiling_eV=np.array(ceil))  # fmt: skip
        conv = [r for r, c in enumerate(cis) if c["ci_status"] == "converged"]
        side_starts = []
        for r in conv:
            side_starts += [cis[r]["ci_positions"] - args.side_step * ctau[r],
                            cis[r]["ci_positions"] + args.side_step * ctau[r]]  # fmt: skip
        sides = (minimize(atoms, side_starts, batch, config, ledger, "e2|sides",
                          fmax=config.quench_fmax_eV_A, max_steps=2000)
                 if side_starts else [])  # fmt: skip
        waves = {}
        for r, ((lb, k, top, npk), c) in enumerate(zip(rows, cis)):
            j = use.index(lb)
            e_r, e_p = ends[2 * j], ends[2 * j + 1]
            arm = {"top": top, "n_peaks": npk, "ci_status": c["ci_status"],
                   "ci_steps": c["ci_steps"], "ci_rel_eV": c["ci_energy_eV"] - e_r.energy_eV,
                   "ci_fmax_eV_A": c["ci_fmax_eV_A"], "string": n_in * (k + 1),
                   "climb": c["ci_forwards"]}  # fmt: skip
            if r in conv:
                q = conv.index(r)
                labs, sc = [], 0
                for m in sides[2 * q : 2 * q + 2]:
                    sc += m.forwards
                    if m.edges == e_r.edges and abs(m.energy_eV - e_r.energy_eV) <= 0.05:
                        labs.append("reactant")
                    elif m.edges == e_p.edges and abs(m.energy_eV - e_p.energy_eV) <= 0.05:
                        labs.append("product")
                    else:
                        labs.append("other")
                arm.update(sides=labs, sides_forwards=sc)
                waves.setdefault(k, []).append({
                    "k": len(waves.get(k, [])), "label": lb, "new_events": [],
                    "ts": [{"status": "converged_candidate", "positions": c["ci_positions"],
                            "energy_eV": c["ci_energy_eV"], "link": {}}],
                })  # fmt: skip
            rec[lb].setdefault("arms", {})[f"k{k}"] = arm
            print(json.dumps({"label": lb, "k": k, **arm}), flush=True)
        for k in ks:
            d = out / f"k{k}"
            d.mkdir(parents=True, exist_ok=True)
            spent = sum(rec[lb]["cost"]["ends"] + rec[lb]["arms"][f"k{k}"]["string"]
                        + rec[lb]["arms"][f"k{k}"]["climb"] for lb in use)  # fmt: skip
            (d / "wave_0000.json").write_text(json.dumps(_jsonable(
                {"wave": 0, "paths": waves.get(k, []), "spent_total": spent}),
                default=str) + "\n")  # fmt: skip
    ledger.close()
    (out / "e2.json").write_text(json.dumps(_jsonable(
        {"args": vars(args), "manifest_args": margs, "replay_s": t_replay,
         "wall_s": time.monotonic() - t0, "paths": [rec[lb] for lb in labels]}),
        indent=1, default=str) + "\n")  # fmt: skip
    print(json.dumps({"wall_s": time.monotonic() - t0, "replay_s": t_replay}), flush=True)


if __name__ == "__main__":
    main()
