"""Read-only summary of midpoint_fire_validate.py output: first-run status transitions
v4 -> v5, per-item outcomes over all repeats, candidate agreement on items where both
arms found one, where each remaining v5 failure stopped, and cost / wall time per group
and arm (first runs only). File names are de-duplicated (183 NFS listing repeated
entries, 2026-11-02).

    python midpoint_fire_validate_summary.py OUT_JSON RESULT_DIR
"""

import collections
import glob
import json
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "explore" / "src"))

from explore.ts_ef import _rmsd  # noqa: E402

NAME = re.compile(r"^(failed|control)_(\d+)_(v[45])(?:_r(\d))?\.json$")


def short(task):
    return task.split(":L0-v1:")[-1].split(":+1")[0].replace("f2.4:L0-v1|", "")


def main():
    out, d = Path(sys.argv[1]), Path(sys.argv[2])
    items = json.loads((d / "items.json").read_text())
    runs = collections.defaultdict(list)
    for name in sorted({Path(f).name for f in glob.glob(str(d / "*.json"))}):
        m = NAME.match(name)
        if m:
            runs[(int(m[2]), m[3])].append((int(m[4] or 0), json.loads((d / name).read_text())))
    for v in runs.values():
        v.sort(key=lambda x: x[0])
    transitions, cost, rows = collections.Counter(), collections.Counter(), []
    for i, (group, task, run) in enumerate(items):
        a = [r for _, r in runs[(i, "v4")]]
        b = [r for _, r in runs[(i, "v5")]]
        transitions[f"{group}: {a[0]['status']} -> {b[0]['status']}"] += 1
        for arm, r in (("v4", a[0]), ("v5", b[0])):
            cost[f"{group}_{arm}_forward"] += r["spent"]
            cost[f"{group}_{arm}_wall_s"] += r["wall_s"]
        agree = None
        if a[0]["candidate"] and b[0]["candidate"]:
            agree = {
                "dE_eV": b[0]["candidate"]["energy_eV"] - a[0]["candidate"]["energy_eV"],
                "rmsd_A": _rmsd(
                    np.array(a[0]["candidate"]["positions"]),
                    np.array(b[0]["candidate"]["positions"]),
                ),
            }
        stop = None
        if b[0]["candidate"] is None and b[0]["midpoints"]:
            last = b[0]["midpoints"][-1]
            stop = {k: last[k] for k in ("origin", "optimizer", "q", "reaction_force", "steps")}
            stop["fmax_projected"] = last["fmax_projected"]
        rows.append(
            {
                "item": i,
                "group": group,
                "intent": short(task),
                "run": run,
                "v4_status": [r["status"] for r in a],
                "v5_status": [r["status"] for r in b],
                "v5_candidate_energy_eV": [(r["candidate"] or {}).get("energy_eV") for r in b],
                "first_run_agreement": agree,
                "v5_first_run_stop": stop,
            }
        )
    summary = {
        "runs": sum(len(v) for v in runs.values()),
        "items": dict(collections.Counter(g for g, _, _ in items)),
        "recheck": json.loads((d / "recheck.json").read_text()),
        "first_run_transitions": dict(transitions),
        "first_run_cost": dict(cost),
        "rows": rows,
    }
    out.write_text(json.dumps(summary, indent=1) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in summary.items() if k != "rows"}, indent=1))


if __name__ == "__main__":
    main()
