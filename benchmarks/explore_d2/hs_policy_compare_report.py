"""Summarise hs_policy_compare.py output (v3 halving retry vs v4 first-failure end, frozen
inputs, same machine). Read-only; run on the host that wrote the results.

Per item the two policies are paired. For each group ("failed": the recorded HS met a
non-converging trial; "control": it did not):
- forwards saved, v3 - v4, per item and summed;
- status pairs v3 -> v4, and whether v3 met a non-converging trial on this machine at all
  (if not, both policies ran the same code and any difference is non-reproducibility);
- candidates (converged_candidate) and continuous maxima under each policy, and whether a
  v3 candidate is kept by v4 (same item);
- control items as the same-machine noise baseline (identical code path).

    python hs_policy_compare_report.py RESULT_DIR OUT_JSON
"""

import json
import sys
from collections import Counter
from pathlib import Path


def main():
    res, out = Path(sys.argv[1]), Path(sys.argv[2])
    pairs = {}
    for f in sorted(res.glob("*_v[34].json")):
        rec = json.loads(f.read_text())
        key = f.name.rsplit("_", 1)[0]  # <run>_<idx>_<group>
        pairs.setdefault(key, {})[rec["policy"]] = rec
    rows = []
    for key, p in sorted(pairs.items()):
        if set(p) != {"v3", "v4"}:
            continue
        v3, v4 = p["v3"], p["v4"]
        rows.append(
            {
                "item": key,
                "block": key.split("_")[3],
                "group": key.rsplit("_", 1)[1],
                "v3_status": v3["status"],
                "v4_status": v4["status"],
                "v3_spent": v3["spent"],
                "v4_spent": v4["spent"],
                "saved": v3["spent"] - v4["spent"],
                "v3_failed_trials": sum(s["unconverged_trials"] for s in v3["sides"]),
                "v4_failed_trials": sum(s["unconverged_trials"] for s in v4["sides"]),
                "v3_branch": v3["branch_status"],
                "v4_branch": v4["branch_status"],
                "recorded": v3["recorded"],
            }
        )
    summary = {}
    for group in ("failed", "control"):
        g = [r for r in rows if r["group"] == group]
        hit = [r for r in g if r["v3_failed_trials"] > 0]
        summary[group] = {
            "items": len(g),
            "v3_met_failed_trial_here": len(hit),
            "v3_spent": sum(r["v3_spent"] for r in g),
            "v4_spent": sum(r["v4_spent"] for r in g),
            "saved_total": sum(r["saved"] for r in g),
            "saved_on_items_v3_hit": sum(r["saved"] for r in hit),
            "saved_on_items_v3_did_not_hit": sum(r["saved"] for r in g if r not in hit),
            "items_identical_spent": sum(r["saved"] == 0 for r in g),
            "status_pairs": dict(Counter(f"{r['v3_status']} -> {r['v4_status']}" for r in g)),
            "v3_candidates": sum(r["v3_status"] == "converged_candidate" for r in g),
            "v4_candidates": sum(r["v4_status"] == "converged_candidate" for r in g),
            "v3_continuous_maximum": sum(r["v3_branch"] == "continuous_maximum" for r in g),
            "v4_continuous_maximum": sum(r["v4_branch"] == "continuous_maximum" for r in g),
            "v3_candidate_lost_by_v4": [
                r["item"]
                for r in g
                if r["v3_status"] == "converged_candidate"
                and r["v4_status"] != "converged_candidate"
            ],
            "v4_candidate_not_in_v3": [
                r["item"]
                for r in g
                if r["v4_status"] == "converged_candidate"
                and r["v3_status"] != "converged_candidate"
            ],
        }
    out.write_text(json.dumps({"summary": summary, "items": rows}, indent=1) + "\n")
    print(json.dumps(summary, indent=1))
    for r in rows:
        print(
            f"{r['item'][12:]:28s} v3 {r['v3_status'][11:]:26s} {r['v3_spent']:5d} "
            f"(fail {r['v3_failed_trials']:2d}) | v4 {r['v4_status'][11:]:26s} "
            f"{r['v4_spent']:5d} (fail {r['v4_failed_trials']}) | saved {r['saved']:5d} "
            f"| rec {r['recorded']['status'][11:]} {r['recorded']['spent']}"
        )


if __name__ == "__main__":
    main()
