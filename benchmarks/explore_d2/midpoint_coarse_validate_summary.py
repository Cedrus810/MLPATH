"""Read-only summary of midpoint_coarse_validate.py output: status transitions v3 -> v4,
cost, candidate agreement (energy difference, Kabsch RMSD), and every v4 coarse stage
(decided, or undecided with the reason visible from its two slopes). File names are
de-duplicated (183 NFS listing repeated entries, 2026-11-02).

    python midpoint_coarse_validate_summary.py OUT_JSON RESULT_DIR BOUND GATE
"""

import collections
import glob
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "explore" / "src"))

from explore.ts_ef import _rmsd  # noqa: E402


def short(task):
    return task.split(":L0-v1:")[-1].split(":+1")[0].replace("f2.4:L0-v1|", "")


def main():
    out, d = Path(sys.argv[1]), Path(sys.argv[2])
    bound, gate = float(sys.argv[3]), float(sys.argv[4])
    items = json.loads((d / "items.json").read_text())
    rec = {
        Path(f).name: json.loads(Path(f).read_text())
        for f in sorted(set(glob.glob(str(d / "*_v[34].json"))))
    }
    rows, transitions, cost = [], collections.Counter(), collections.Counter()
    coarse = {"decided": 0, "below_bound_or_mixed": [], "blocked_beyond_bound": []}
    for i, (group, task, run) in enumerate(items):
        a, b = rec[f"{group}_{i:02d}_v3.json"], rec[f"{group}_{i:02d}_v4.json"]
        transitions[f"{group}: {a['status']} -> {b['status']}"] += 1
        cost[f"{group}_v3"] += a["spent"]
        cost[f"{group}_v4"] += b["spent"]
        stages = [m for m in b["midpoints"] if m["origin"].endswith("_coarse")]
        for x, y in zip(stages[::2], stages[1::2]):
            if x["coarse_sign"]:
                coarse["decided"] += 1
                continue
            rx, ry = x["reaction_force"], y["reaction_force"]
            entry = {"item": i, "q": x["q"], "rf": [rx, ry], "diff": abs(rx - ry)}
            if abs(rx) > bound and abs(ry) > bound and rx * ry > 0:
                coarse["blocked_beyond_bound"].append(entry)
            else:
                coarse["below_bound_or_mixed"].append(entry)
        agree = None
        if a["candidate"] and b["candidate"]:
            agree = {
                "dE_eV": b["candidate"]["energy_eV"] - a["candidate"]["energy_eV"],
                "rmsd_A": _rmsd(
                    np.array(a["candidate"]["positions"]), np.array(b["candidate"]["positions"])
                ),
            }
        rows.append(
            {
                "item": i,
                "group": group,
                "intent": short(task),
                "run": run,
                "v3": {"status": a["status"], "spent": a["spent"]},
                "v4": {
                    "status": b["status"],
                    "spent": b["spent"],
                    "coarse_decisions": sum(1 for m in stages[::2] if m["coarse_sign"]),
                    "candidate_energy_eV": (b["candidate"] or {}).get("energy_eV"),
                },
                "candidate_agreement": agree,
            }
        )
    summary = {
        "items": collections.Counter(g for g, _, _ in items),
        "bound": bound,
        "gate": gate,
        "transitions": dict(transitions),
        "cost": dict(cost),
        "coarse_stages": coarse,
        "rows": rows,
    }
    out.write_text(json.dumps(summary, indent=1) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in summary.items() if k not in ("rows",)}, indent=1))


if __name__ == "__main__":
    main()
