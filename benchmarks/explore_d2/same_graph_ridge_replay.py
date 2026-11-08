"""Offline replay of D2-HF-v1.2 promotion order on finished D2-HF runs
(docs/PLAN_HF_SAME_GRAPH_RIDGE_2026-11-07.md §4; user 2026-11-07 「可以」 to R3a-2).
Read-only, no model: every promoted intent is re-classified with `d2._classify` and marked
with `d2._same_graph_ridge` (the code the new rule runs); its recorded HS / HC outcome and
cost are attached. Reports, per run: promoted items by class, how many the new rule puts
last, their recorded HS forwards, and every candidate / HC-connected item among them.

    python same_graph_ridge_replay.py OUT_JSON RUN_DIR [RUN_DIR ...]
"""

import json
import sys
from collections import Counter
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent))

from midpoint_cold_replay import load  # noqa: E402  (sets sys.path)

from explore.d2 import _classify, _points, _same_graph_ridge, _task  # noqa: E402


def short(task):
    return task.split(":L0-v1:")[-1].split(":+1")[0].replace("f2.4:L0-v1|", "")


def replay(run_dir):
    run = load(Path(run_dir))
    view = SimpleNamespace(settled=run.settled, settings=run.settings)
    rows = []
    for row in run.rows.values():
        if _task("L1", row) not in run.settled:
            continue
        classified = _classify(view, row)
        if classified["class"] is None:
            continue
        points = _points(run.settled[_task("L1", row)]["payload"])
        hs = run.settled.get(_task("HS", row))
        hc = run.settled.get(_task("HC", row))
        rows.append(
            {
                "intent": short(_task("HS", row)),
                "class": classified["class"],
                "late": _same_graph_ridge(classified, points),
                "hs_status": None if hs is None else hs["status"],
                "hs_spent": 0 if hs is None else hs["spent"],
                "hc_status": None if hc is None else hc["status"],
                "hc_spent": 0 if hc is None else hc["spent"],
            }
        )
    late = [r for r in rows if r["late"]]
    return {
        "run": run.name,
        "promoted": dict(Counter(r["class"] for r in rows)),
        "late": len(late),
        "late_hs_status": dict(Counter(r["hs_status"] for r in late)),
        "hs_spent_all": sum(r["hs_spent"] for r in rows),
        "hs_spent_late": sum(r["hs_spent"] for r in late),
        "late_candidates": [r for r in late if r["hs_status"] == "converged_candidate"],
        "late_connected": [r for r in late if r["hc_status"] == "connected"],
        "candidates_all": sum(r["hs_status"] == "converged_candidate" for r in rows),
        "connected_all": sum(r["hc_status"] == "connected" for r in rows),
        "rows": rows,
    }


def main():
    out, runs = Path(sys.argv[1]), sys.argv[2:]
    results = [replay(r) for r in runs]
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(results, indent=1) + "\n", encoding="utf-8")
    for r in results:
        print(
            json.dumps(
                {k: r[k] for k in r if k not in ("rows", "late_candidates", "late_connected")}
                | {"late_candidates": len(r["late_candidates"])}
                | {"late_connected": len(r["late_connected"])}
            ),
            flush=True,
        )


if __name__ == "__main__":
    main()
