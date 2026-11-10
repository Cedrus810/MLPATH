"""HF-5 b04 run table (docs/PLAN_HF5_CROSS_SYSTEM_2026-11-08.md §7). Read-only: manifest.json
and report.json of each run; the journal is read only for HS settles (class, same-graph tag,
sample tag, status, spent) and HC settles.

    python hf5_b04_summary.py OUT_JSON RUN_DIR [RUN_DIR ...]
"""

import collections
import json
import sys
from pathlib import Path


def _settles(run):
    with open(run / "journal.jsonl", encoding="utf-8") as fh:
        for line in fh:
            if '"event": "settle"' not in line and '"event":"settle"' not in line:
                continue
            ev = json.loads(line)
            if ev.get("event") == "settle" and ev["task"][:3] in ("HS|", "HC|"):
                yield ev


def one(run):
    run = Path(run)
    m = json.loads((run / "manifest.json").read_text(encoding="utf-8"))
    r = json.loads((run / "report.json").read_text(encoding="utf-8"))
    hs, hc = {}, {}
    for ev in _settles(run):
        (hs if ev["task"].startswith("HS|") else hc)[ev["task"][3:]] = ev
    groups = collections.Counter()
    spent = collections.Counter()
    for iid, ev in hs.items():
        p = ev["payload"]
        key = (
            p.get("class"),
            "same_graph" if p.get("same_graph_ridge") else "bond_change",
            ev["status"],
            (hc[iid]["status"] if iid in hc else None),
        )
        groups[" / ".join(map(str, key))] += 1
        spent["same_graph" if p.get("same_graph_ridge") else "other"] += ev["spent"]
    ident = m["identity"]
    return {
        "run": run.name,
        "rule": m["rule_version"],
        "seed": m["settings"]["seed"],
        "rotate": (ident.get("args") or {}).get("rotate"),
        "magnitude_points": len(m["settings"]["magnitudes_A"]),
        "hs_bracket": m.get("hs_bracket"),
        "population": r.get("population"),
        "promotion": r.get("promotion"),
        "hs_status": r.get("hs_status"),
        "hc_status": r.get("hc_status"),
        "same_graph_ridge": {
            k: v
            for k, v in (r.get("same_graph_ridge") or {}).items()
            if k != "sample_candidates"
        }
        | {
            "sample_candidates": len(
                (r.get("same_graph_ridge") or {}).get("sample_candidates", [])
            )
        },
        "hs_by_class_graph_status_hc": dict(sorted(groups.items())),
        "hs_forwards_by_graph": dict(spent),
        "accounts_spent": {a: v["spent"] for a, v in r["accounts"].items()},
        "ledger_total": r.get("ledger_total"),
        "hf_forwards": r.get("hf_forwards"),
        "unreconciled": r.get("unreconciled"),
        "events": [
            {
                "saddle_energy_eV": e.get("saddle_energy_eV"),
                "broken": e.get("broken"),
                "formed": e.get("formed"),
                "cumulative_spent": e.get("cumulative_spent"),
                "class": e.get("class"),
            }
            for e in r.get("hf_events", [])
        ],
        "hosts": r.get("hosts"),
    }


def main():
    out = Path(sys.argv[1])
    rows = [one(p) for p in sys.argv[2:]]
    out.write_text(json.dumps(rows, indent=1) + "\n", encoding="utf-8")
    print(json.dumps(rows, indent=1))


if __name__ == "__main__":
    main()
