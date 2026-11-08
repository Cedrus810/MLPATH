"""Rebuild the HS bracket of every branch_bracketed item of finished D2-HF runs with the
shipped `_classify` / `_hf_bracket` (2026-10-28 fix: forward jump interval) and report
its q width next to the old hand-off (reactant point + L2 reverse point) and the HS
status the run recorded. Reads journals only; no model.

    python replay_hf_brackets.py OUT_JSON RUN_DIR [RUN_DIR ...]
"""

import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / "explore" / "src"), str(ROOT / "src")]

from explore.d2 import _classify, _hf_bracket, _intent, _points, _task  # noqa: E402


def _lines(path):
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.strip():
                yield json.loads(line)


def replay(run_dir):
    run_dir = Path(run_dir)
    m = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    settled = {}
    for ev in _lines(run_dir / "journal.jsonl"):
        if ev.get("event") != "settle":
            continue
        if ev["task"][:3] in ("HS|", "HC|"):
            ev["payload"] = {"error": ev["payload"].get("error")}
        settled[ev["task"]] = ev
    run = SimpleNamespace(
        settled=settled, settings=SimpleNamespace(noise_eV=m["settings"]["noise_eV"])
    )
    rows = []
    for row in _lines(run_dir / "population.jsonl"):
        if _task("L1", row) not in settled:
            continue
        c = _classify(run, row)
        if c["class"] != "branch_bracketed":
            continue
        points = _points(settled[_task("L1", row)]["payload"])
        a, b = _hf_bracket(c, points)
        product = np.asarray(c["info"]["product_side_positions"], float)
        old = abs(_intent(row).value(product) - points[c["info"]["reactant_side"]].q_target)
        hs = settled.get(_task("HS", row))
        rows.append(
            {
                "intent": row["intent_id"],
                "new_width_A": abs(b.q_target - a.q_target),
                "new_ends_edges_differ": a.edges != b.edges,
                "old_width_A": float(old),
                "hs_status": None if hs is None else hs["status"],
                "hs_spent": None if hs is None else hs["spent"],
            }
        )
    new = [r["new_width_A"] for r in rows]
    return {
        "run": str(run_dir),
        "rule_version": m["rule_version"],
        "branch_bracketed": len(rows),
        "new_width_min_A": min(new) if new else None,
        "new_width_max_A": max(new) if new else None,
        "new_ends_on_two_graphs": sum(r["new_ends_edges_differ"] for r in rows),
        "old_width_zero_or_eps": sum(r["old_width_A"] < 1e-12 for r in rows),
        "items": rows,
    }


def main():
    out, runs = Path(sys.argv[1]), sys.argv[2:]
    results = [replay(r) for r in runs]
    out.write_text(json.dumps(results, indent=1), encoding="utf-8")
    for r in results:
        print(
            r["run"],
            r["rule_version"],
            {k: r[k] for k in r if k not in ("run", "rule_version", "items")},
        )
        for it in r["items"]:
            if it["hs_status"] is not None:
                print(
                    "   ",
                    round(it["old_width_A"], 18),
                    "->",
                    round(it["new_width_A"], 4),
                    it["hs_status"],
                    it["hs_spent"],
                )


if __name__ == "__main__":
    main()
