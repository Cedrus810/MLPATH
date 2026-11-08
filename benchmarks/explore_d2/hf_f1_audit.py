"""Read-only F1/F6 diagnostic audit; never repairs historical journals.

python benchmarks/explore_d2/hf_f1_audit.py OUTPUT_JSON
"""

import json
import sys
from collections import Counter
from pathlib import Path

from f6_ab import l4


def audit(root):
    report = {}
    for block in ("sA", "sB"):
        old, _ = l4(root / "results/2026-10-23" / f"b05_fp32_{block}")
        new, _ = l4(root / "results/2026-10-24" / f"b05_v11_{block}")
        lost = [r for r in new.values() if r["status"] == "unresolved:saddle_lost_target_mode"]
        diagnostics = [r["payload"]["diagnostics"]["saddle_refine"] for r in lost]
        old_hits = {}
        for task, row in old.items():
            if row["status"] not in ("c2", "duplicate_c2"):
                continue
            now = new.get(task)
            old_hits[task] = {
                "old_status": row["status"],
                "old_event": row["payload"].get("event"),
                "new_status": None if now is None else now["status"],
                "new_refine": None
                if now is None
                else now["payload"].get("diagnostics", {}).get("saddle_refine"),
            }
        report[block] = {
            "lost_total": len(lost),
            "lost_at_step_zero": sum(d["steps"] == 0 for d in diagnostics),
            "null_initial_overlap_steps": dict(
                Counter(d["steps"] for d in diagnostics if d.get("target_overlap") is None)
            ),
            "lost_after_initial_selection": sum(d["steps"] > 0 for d in diagnostics),
            "old_hits": old_hits,
        }
    report["interpretation"] = {
        "initial_null": "All null target_overlap records exit at step 0; stale **bounds overwrote measured target_overlap and selected_mode_index. Historical values are not reconstructed.",
        "tracking": "After initial selection the anchor is the previous eigenvector, not the fixed scan tangent. Loss of overlap alone does not establish physical direction loss or methyl rotation.",
        "acceptance": "F6 alone did not recover the reactive E2 hit; the new E/F locator needs separate molecular validation.",
    }
    return report


if __name__ == "__main__":
    output = Path(sys.argv[1])
    output.parent.mkdir(parents=True, exist_ok=True)
    report = audit(Path(__file__).resolve().parent)
    output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    for block in ("sA", "sB"):
        print(block, {k: v for k, v in report[block].items() if k != "old_hits"})
