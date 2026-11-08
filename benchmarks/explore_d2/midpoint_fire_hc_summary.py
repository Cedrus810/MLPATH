"""Read-only summary of midpoint_fire_hc.py output: HS and HC outcome per intent, side
status of unresolved connections, chemical-pair agreement with the reference new-E2
event, and cost. File names are de-duplicated (183 NFS listing repeated entries).

    python midpoint_fire_hc_summary.py OUT_JSON RESULT_DIR
"""

import collections
import glob
import json
import re
import sys
from pathlib import Path


def short(task):
    return task.split(":L0-v1:")[-1].split(":+1")[0].replace("f2.4:L0-v1|", "")


def main():
    out, d = Path(sys.argv[1]), Path(sys.argv[2])
    names = sorted(
        {Path(f).name for f in glob.glob(str(d / "failed_*.json"))}
        - {Path(f).name for f in glob.glob(str(d / "*.ledger.jsonl"))}
    )
    rows = []
    for name in names:
        if not re.match(r"^failed_\d+\.json$", name):
            continue
        r = json.loads((d / name).read_text(encoding="utf-8"))
        rows.append(
            {
                "file": name,
                "intent": short(r["task"]),
                "run": r["run"],
                "hs_status": r["hs_status"],
                "hs_spent": r["hs_spent"],
                "saddle_energy_eV": r.get("saddle_energy_eV"),
                "hc_status": r.get("hc_status"),
                "hc_spent": r.get("hc_spent"),
                "sides": [
                    {k: s[k] for k in ("side", "status", "end_energy_eV")}
                    for s in r.get("sides", [])
                ],
                "pair_matches_reference": r.get("pair_matches_reference"),
            }
        )
    summary = {
        "items": len(rows),
        "hs_status": dict(collections.Counter(r["hs_status"] for r in rows)),
        "hc_status": dict(collections.Counter(r["hc_status"] for r in rows if r["hc_status"])),
        "connected_pair_matches_reference": sum(
            1 for r in rows if r["hc_status"] == "connected" and r["pair_matches_reference"]
        ),
        "hs_spent": sum(r["hs_spent"] for r in rows),
        "hc_spent": sum(r["hc_spent"] or 0 for r in rows),
        "rows": rows,
    }
    out.write_text(json.dumps(summary, indent=1) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in summary.items() if k != "rows"}, indent=1))


if __name__ == "__main__":
    main()
