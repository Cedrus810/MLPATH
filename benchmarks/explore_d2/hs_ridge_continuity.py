"""Does the L1 scan already show the geometric jump behind a false ridge? (TODO 9c'' R3; F4 of
TS_HESSIAN_FREE_DESIGN_RECORD_2026-10-24 §4). Read-only, from finished D2-HF journals; no model.

For every HS whose intent was promoted as `ridge_suspected`, the intent's L1 scan points are
read and the current `read_scan` ridge trigger is recomputed (peak index k). Written per HS:
the HS branch status and final status; the aligned (Kabsch) RMSD between consecutive scan
points at every step; the two steps adjacent to the peak (k-1 -> k, k -> k+1); the median of
the other steps; their ratio; the reaction forces at k-1, k, k+1; the energy rise and drop
around the peak. No threshold is applied.

    python hs_ridge_continuity.py OUT_JSON RUN_DIR [RUN_DIR ...]
"""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / "explore" / "src"), str(ROOT / "src")]

from explore.d2 import _points  # noqa: E402
from explore.probe import read_scan  # noqa: E402
from explore.ts_ef import _rmsd  # noqa: E402


def short(task):
    return task.split(":L0-v1:")[-1].split(":+1")[0].replace("f2.4:L0-v1|", "")


def scan(run_dir):
    run_dir = Path(run_dir)
    manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    noise = manifest["settings"]["noise_eV"]
    l1, hs = {}, {}
    with open(run_dir / "journal.jsonl", encoding="utf-8") as fh:
        for line in fh:
            if '"settle"' not in line or ('"L1|' not in line and '"HS|' not in line):
                continue
            ev = json.loads(line)
            if ev.get("event") != "settle":
                continue
            kind, intent = ev["task"].split("|", 1)
            if kind == "L1":
                l1[intent] = ev["payload"]
            elif kind == "HS" and ev["payload"].get("class") == "ridge_suspected":
                hs[intent] = (ev["status"], ev["payload"].get("branch_status"), ev["spent"])
    out = []
    for intent, (status, branch, spent) in hs.items():
        points = _points(l1[intent])
        ridge = [t for t in read_scan(points, noise) if t["kind"] == "ridge_suspected"]
        if not ridge:
            continue
        k = ridge[0]["peak"]
        steps = [_rmsd(a.positions, b.positions) for a, b in zip(points, points[1:])]
        adjacent = [steps[k - 1], steps[k]]
        others = [s for i, s in enumerate(steps) if i not in (k - 1, k)]
        base = float(np.median(others)) if others else None
        e = [p.energy_eV for p in points]
        out.append(
            {
                "run": run_dir.name,
                "intent": short("HS|" + intent),
                "status": status,
                "branch_status": branch,
                "spent": spent,
                "peak": k,
                "n_points": len(points),
                "steps_A": steps,
                "adjacent_A": adjacent,
                "others_median_A": base,
                "ratio": None if not base else max(adjacent) / base,
                "rf": [points[i].reaction_force for i in (k - 1, k, k + 1)],
                "rise_eV": e[k] - e[k - 1],
                "drop_eV": e[k] - e[k + 1],
                "same_edges": points[k].edges == points[k + 1].edges,
            }
        )
    return out


def main():
    out, runs = Path(sys.argv[1]), sys.argv[2:]
    rows = [r for run in runs for r in scan(run)]
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rows) + "\n", encoding="utf-8")
    print(len(rows), "ridge HS rows")


if __name__ == "__main__":
    main()
