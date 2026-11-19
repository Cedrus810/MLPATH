"""Read-only brief of D2 run reports (report.json); run on the host that wrote them.

python v15_report_brief.py RUN_DIR [RUN_DIR ...]
"""

import json
import sys
from pathlib import Path

for d in map(Path, sys.argv[1:]):
    r = json.loads((d / "report.json").read_text())
    print(
        d.name,
        r["rule_version"],
        "ledger",
        r["ledger_total"],
        "unreconciled",
        len(r["unreconciled"]),
        json.dumps(r["accounts"]),
    )
    print("  promotion", r["promotion"], "skips", r["skips"])
    print("  same_graph", r.get("same_graph_ridge"))
    print("  geometry_noise", json.dumps(r.get("geometry_noise")))
    print("  hs", r["hs_status"], "hc", r["hc_status"], r["hf_forwards"])
    for e in r["hf_events"]:
        print(
            "   event",
            round(e["saddle_energy_eV"], 4),
            e["broken"],
            e["formed"],
            e["cumulative_spent"],
        )
