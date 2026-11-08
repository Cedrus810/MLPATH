"""Read-only summary of midpoint_transverse_hessian.py output, joined with the cold-start
replay of the same failed frames (midpoint_cold_replay.py, group C, matched on run, task
and q). Writes OUT_JSON and prints the tables.

    python midpoint_transverse_summary.py OUT_JSON HESSIAN_DIR COLD_DIR
"""

import collections
import json
import sys
from pathlib import Path

import numpy as np


def q(v, ps=(0, 10, 50, 90, 100)):
    v = np.asarray(v, dtype=float)
    return {f"p{p}": float(np.percentile(v, p)) for p in ps} if v.size else None


def short(task):
    return task.split(":L0-v1:")[-1].split(":+1")[0]


def main():
    out, hdir, cdir = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    recs = [json.loads(p.read_text()) for p in sorted(hdir.glob("[FS]_[0-9]*.json"))]
    groups = {g: [r for r in recs if r["group"] == g] for g in ("F", "S")}

    cold_items = json.loads((cdir / "CD_items.json").read_text())["C"]
    cold = collections.defaultdict(list)
    for p in cdir.glob("C_*_r*.json"):
        r = json.loads(p.read_text())
        cold[(r["run"], r["task"], r["q"])].append(r["converged"])

    summary = {"n": {g: len(v) for g, v in groups.items()}, "groups": {}}
    for g, rs in groups.items():
        eig = np.array([r["transverse_eigs"][:5] for r in rs])
        mw = np.array([r["transverse_eigs_mass_weighted"][:3] for r in rs])
        entry = {
            "lowest_transverse": q(eig[:, 0]),
            "second_transverse": q(eig[:, 1]),
            "third_transverse": q(eig[:, 2]),
            "gap_second_over_lowest": q(eig[:, 1] / eig[:, 0]),
            "negative_lowest": int((eig[:, 0] < 0).sum()),
            "lowest_mass_weighted": q(mw[:, 0]),
            "rigid_only_lowest": q([r["rigid_only_lowest"][0] for r in rs]),
            "rigid_only_lowest_cos_grad_q": q(
                [r["rigid_only_lowest_cos_grad_q"][0] for r in rs]
            ),
            "fd_check_max": q([max(r["fd_check_deviations"]) for r in rs]),
            "fproj": q([r["fproj"] for r in rs]),
            "wall_s": q([r["wall_s"] for r in rs]),
        }
        if g == "F":
            entry |= {
                "wander_aligned_A": q([r["wander_norm_A"] for r in rs]),
                "wander_transverse_fraction": q([r["wander_transverse_fraction"] for r in rs]),
                "wander_fraction_lowest_1": q([r["wander_fraction_lowest"]["1"] for r in rs]),
                "wander_fraction_lowest_3": q([r["wander_fraction_lowest"]["3"] for r in rs]),
                "wander_fraction_lowest_5": q([r["wander_fraction_lowest"]["5"] for r in rs]),
            }
        by_intent = collections.defaultdict(list)
        for r in rs:
            by_intent[short(r["task"])].append(r["transverse_eigs"][0])
        entry["lowest_by_intent"] = {
            k: {"n": len(v), "p50": float(np.median(v))} for k, v in sorted(by_intent.items())
        }
        summary["groups"][g] = entry

    joined = []
    for r in groups["F"]:
        c = cold.get((r["run"], r["task"], r["q"]))
        if c is not None:
            joined.append((r["transverse_eigs"][0], sum(c), len(c)))
    by_k = collections.defaultdict(list)
    for lam, k, n in joined:
        by_k[f"{k}/{n}"].append(lam)
    summary["F_lowest_by_cold_outcome"] = {
        k: {"n": len(v), **q(v, (0, 50, 100))} for k, v in sorted(by_k.items())
    }
    summary["F_matched_to_cold"] = (
        f"{len(joined)} / {len(groups['F'])} (cold items {len(cold_items)})"
    )

    out.write_text(json.dumps(summary, indent=1) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
