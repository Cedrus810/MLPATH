"""When is the path direction good enough to hand over to the climbing image? (2026-11-19;
user-pasted research route: "让 AFIR 路径提供方向，用很少的局部优化直接定位鞍点……最该研究的
问题是：方向究竟在第几步就已经够用了，能否提前交给爬坡？")

For each cached N1 DA path (frames_<label>.npz written by n1_string_check.py; no replay):
  1. end points: batched E/F minima of the first and last path frame (fmax 0.01);
     interior: the AFIR frames from the early low-energy frame to the end, resampled at
     equal arc length (20 nodes incl. ends), as the string start of n1_string_check.py;
  2. ONE simplified-string run (explore.afir_string: full-force step, max 0.05 A per node,
     equal-arc reparametrization), snapshots after k iterations, k in --checkpoints
     (k = 0: the AFIR path itself, no string iteration);
  3. from every snapshot independently: climbing image from the highest interior node along
     the fixed tangent of its two neighbours (afir_string.climb, fmax 0.01 eV/A, ceiling
     start minimum + 3 eV);
  4. converged images: descend from x +/- 0.05 A tau (batched E/F minimization, quench
     fmax) and label each side by bond graph and energy against the two end minima
     ("reactant" / "product" / "other"; a rough check made for this study, not the HC rule).
Cost of arm k = end minimizations + 18 * (k + 1) string evaluations up to the snapshot +
its own climb. Per arm k, OUT/k<k>/<seed>/wave_0000.json in the runner's wave format for
benchmarks/afir_mc/afir_validate.py. Not resumable (short).

    python n1_direction_study.py --out OUT --case SEED:LABEL:NPZ [--case ...]
"""

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import run_mc_afir  # noqa: E402,F401  (sets sys.path for run_d2 / systems)
from run_d2 import make_atoms, make_config  # noqa: E402
from systems import SYSTEMS  # noqa: E402

from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.afir_batch import _evaluate, minimize  # noqa: E402
from explore.afir_string import climb, reparametrize, resample, tangent  # noqa: E402
from explore.ledger import Ledger  # noqa: E402
from explore.parallel import MACEBatchEvaluator  # noqa: E402
from explore.poses import _jsonable  # noqa: E402
from explore.ts_ef import _align  # noqa: E402

FORMING = {"C1-C6": (0, 11), "C4-C7": (3, 12), "C1-C7": (0, 12), "C4-C6": (3, 11)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--case", action="append", required=True, help="SEEDTAG:LABEL:NPZ")
    ap.add_argument("--checkpoints", default="0,10,25,50,100,200,400")
    ap.add_argument("--nodes", type=int, default=20)
    ap.add_argument("--fmax", type=float, default=0.01)
    ap.add_argument("--ci-steps", type=int, default=2000)
    ap.add_argument("--side-step", type=float, default=0.05, help="A, along tau")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    ks = sorted(int(k) for k in args.checkpoints.split(","))
    t0 = time.monotonic()
    system = SYSTEMS["n1"]
    config = make_config(system)
    calc = polar1("cuda", "float32")
    atoms, _ = make_atoms(system, config)
    batch = MACEBatchEvaluator(calc)
    ledger = Ledger(out / "study.ledger.jsonl", arm="DIR", replay=False)
    report, waves = [], {}
    for case in args.case:
        tag, label, npz = case.split(":", 2)
        z = np.load(npz)
        pos, en = z["positions"], z["energies"]
        t1 = time.monotonic()
        s0 = ledger.spent
        ends = minimize(atoms, [pos[0], pos[-1]], batch, config, ledger, f"{label}|ends",
                        fmax=args.fmax, max_steps=2000)  # fmt: skip
        e0 = ends[0].energy_eV
        g_r, g_p = ends[0].edges, ends[1].edges
        ends_cost = ledger.spent - s0
        low = int(np.argmin(en[: len(en) // 3 + 1]))
        mids = resample(list(pos[low:]), args.nodes)[1:-1]
        x = reparametrize(np.stack([ends[0].positions] + list(mids) + [ends[1].positions]))
        snaps, string_cost = {}, 0
        with ledger.segment(task=f"{label}|string", stage="string", attempt=0) as seg:
            for it in range(0, ks[-1] + 1):
                e, f = _evaluate(batch, atoms, x[1:-1], ledger)
                string_cost += len(x) - 2
                if it in ks:
                    snaps[it] = (x.copy(), [e0] + [float(v) for v in e] + [ends[1].energy_eV],
                                 string_cost)  # fmt: skip
                if it == ks[-1]:
                    break
                step = 0.05 * f
                size = np.linalg.norm(step, axis=2).max(axis=1)
                step *= np.minimum(1.0, 0.05 / np.where(size > 0, size, 1.0))[:, None, None]
                moved = x.copy()
                moved[1:-1] += step
                x = reparametrize(moved)
            seg.outcome = f"{ks[-1]} iterations"
        for k in ks:
            xs, energies, scost = snaps[k]
            top = int(np.argmax(energies[1:-1])) + 1
            tau = tangent(xs, top)
            ci = climb(atoms, xs[top], tau, batch, ledger, f"{label}|k{k}", ci_fmax=args.fmax,
                       ci_steps=args.ci_steps, ceiling_eV=e0 + 3.0)  # fmt: skip
            rec = {"seed": tag, "label": label, "k": k, "top": top,
                   "peak_rel_eV": energies[top] - e0, "ci_status": ci["ci_status"],
                   "ci_steps": ci["ci_steps"], "ci_rel_eV": ci["ci_energy_eV"] - e0,
                   "ci_fmax_eV_A": ci["ci_fmax_eV_A"],
                   "start_to_ci_rmsd_A": float(np.sqrt(np.mean(np.sum(
                       (ci["ci_positions"] - _align(ci["ci_positions"], xs[top])) ** 2, axis=1)))),
                   "ci_forming_A": {n: round(float(np.linalg.norm(ci["ci_positions"][j]
                                                                   - ci["ci_positions"][i])), 3)
                                    for n, (i, j) in FORMING.items()},
                   "cost": {"ends": ends_cost, "string": scost, "climb": ci["ci_forwards"]},
                   "ci_trace": ci["ci_trace"]}  # fmt: skip
            if ci["ci_status"] == "converged":
                xc = ci["ci_positions"]
                c1 = ledger.spent
                sides = minimize(atoms, [xc - args.side_step * tau, xc + args.side_step * tau],
                                 batch, config, ledger, f"{label}|k{k}|sides",
                                 fmax=config.quench_fmax_eV_A, max_steps=2000)  # fmt: skip
                labels = []
                for m in sides:
                    if m.edges == g_r and abs(m.energy_eV - e0) <= 0.05:
                        labels.append("reactant")
                    elif m.edges == g_p and abs(m.energy_eV - ends[1].energy_eV) <= 0.05:
                        labels.append("product")
                    elif m.edges == g_r:
                        labels.append("reactant_graph")
                    elif m.edges == g_p:
                        labels.append("product_graph")
                    else:
                        labels.append("other")
                rec["sides"] = labels
                rec["sides_rel_eV"] = [m.energy_eV - e0 for m in sides]
                rec["cost"]["sides"] = ledger.spent - c1
                waves.setdefault(k, {}).setdefault(tag, []).append({
                    "k": len(waves.get(k, {}).get(tag, [])), "label": label, "new_events": [],
                    "ts": [{"status": "converged_candidate", "positions": xc,
                            "energy_eV": ci["ci_energy_eV"], "link": {}}],
                })  # fmt: skip
            rec["cost"]["total_locate"] = (rec["cost"]["ends"] + rec["cost"]["string"]
                                           + rec["cost"]["climb"])  # fmt: skip
            report.append(rec)
            print(json.dumps({q: v for q, v in rec.items() if q != "ci_trace"}), flush=True)
        print(json.dumps({"case": label, "wall_s": time.monotonic() - t1}), flush=True)
    ledger.close()
    for k, by_tag in waves.items():
        for tag, paths in by_tag.items():
            d = out / f"k{k}" / tag
            d.mkdir(parents=True, exist_ok=True)
            spent = sum(r["cost"]["total_locate"] for r in report
                        if r["k"] == k and r["seed"] == tag)  # fmt: skip
            (d / "wave_0000.json").write_text(json.dumps(_jsonable(
                {"wave": 0, "paths": paths, "spent_total": spent}), default=str) + "\n")  # fmt: skip
    (out / "study.json").write_text(json.dumps(_jsonable(
        {"args": vars(args), "wall_s": time.monotonic() - t0, "rows": report}), indent=1,
        default=str) + "\n")  # fmt: skip
    print(json.dumps({"wall_s": time.monotonic() - t0, "forwards": ledger.spent}), flush=True)


if __name__ == "__main__":
    main()
