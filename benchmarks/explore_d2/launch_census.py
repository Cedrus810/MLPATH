"""Where does D2 launch, and what does each launch buy? (algorithm-defect census, 2026-11-14;
user: "发射的地方根本不在反应附近的 没用发射"). Read-only over finished runs; no model.

A launch = one L1 scan of one intent from one L0 pose. Per launch: its forward cost and its
furthest outcome -- event (its HC connected or hf_duplicate), candidate (HS converged),
searched (HS ran), promoted-not-searched (HS skipped / deferred), trigger-free (no ridge /
jump / audit), and the scan status (completed / censored_energy / ...). Grouped by what is
known before the launch: intent layer (bonded / near / far), scope (inter / intra),
element pair, start distance q0, and the pose.

    python launch_census.py OUT_JSON RUN_DIR [RUN_DIR ...]
"""

import collections
import json
import sys
from pathlib import Path

import numpy as np
from ase.data import chemical_symbols


def census(run):
    pop = [json.loads(x) for x in (run / "population.jsonl").read_text().splitlines() if x]
    poses = {}
    for line in (run / "poses.jsonl").read_text().splitlines():
        if line:
            r = json.loads(line)
            poses[r["pose_id"]] = r
    cost, status, q0 = collections.Counter(), {}, {}
    skipped = {}
    with open(run / "journal.jsonl", encoding="utf-8") as fh:
        for line in fh:
            ev = json.loads(line)
            t = ev.get("task", "")
            if ev["event"] in ("settle", "reconcile"):
                cost[t] += ev["spent"]
                if ev["event"] == "settle":
                    status[t] = ev["status"]
                    if t.startswith("L1|"):
                        q0[t[3:]] = ev["payload"].get("q0")
            elif ev["event"] == "skip":
                skipped[t] = ev["reason"]
    rows = []
    for r in pop:
        iid = r["intent_id"]
        l1 = "L1|" + iid
        if l1 not in status:
            continue
        nums = poses[r["pose_id"]]["numbers"]
        i, j = r["pairs"][0]
        pair = "-".join(sorted((chemical_symbols[nums[i]], chemical_symbols[nums[j]])))
        hs, hc = "HS|" + iid, "HC|" + iid
        if status.get(hc) in ("connected", "hf_duplicate"):
            out = "event"
        elif status.get(hs) == "converged_candidate":
            out = "candidate"
        elif hs in status:
            out = "searched"
        elif hs in skipped:
            out = "promoted_not_searched"
        else:
            out = "no_promotion"
        own = sum(cost[k + iid] for k in ("L1|", "L2|", "L3|"))
        rows.append(
            {
                "intent": iid,
                "pose": r["pose_id"],
                "layer": r["layer"],
                "scope": r["scope"].split(":")[0],
                "pair": pair,
                "sign": r["sign"],
                "q0": q0.get(iid),
                "l1_status": status[l1],
                "outcome": out,
                "search_cost": own,
                "refine_cost": cost[hs] + cost[hc],
            }
        )
    return rows


def table(rows, key):
    t = collections.defaultdict(lambda: collections.Counter())
    for r in rows:
        k = key(r)
        t[k]["launches"] += 1
        t[k]["search_cost"] += r["search_cost"]
        t[k]["refine_cost"] += r["refine_cost"]
        t[k][r["outcome"]] += 1
    return {str(k): dict(v) for k, v in sorted(t.items(), key=lambda kv: -kv[1]["search_cost"])}


def main():
    out = Path(sys.argv[1])
    report = {}
    for d in map(Path, sys.argv[2:]):
        rows = census(d)
        total = sum(r["search_cost"] for r in rows)
        productive = [r for r in rows if r["outcome"] in ("event", "candidate")]
        q_event = [r["q0"] for r in rows if r["outcome"] == "event" and r["q0"] is not None]
        report[d.name] = {
            "launches": len(rows),
            "search_cost": total,
            "outcome": table(rows, lambda r: r["outcome"]),
            "layer_scope": table(rows, lambda r: (r["layer"], r["scope"])),
            "pair": table(rows, lambda r: r["pair"]),
            "l1_status": table(rows, lambda r: r["l1_status"]),
            "q0_bins": table(
                rows,
                lambda r: None if r["q0"] is None else float(np.digitize(r["q0"], [1.5, 2.5, 3.5, 4.5, 6.0])),
            ),
            "productive_share_of_search_cost": sum(r["search_cost"] for r in productive) / total,
            "event_q0": sorted(q_event),
            "event_rows": [r for r in rows if r["outcome"] == "event"],
        }
        print(d.name, "launches", len(rows), "search", total,
              "productive share", round(report[d.name]["productive_share_of_search_cost"], 3))
        for k, v in report[d.name]["outcome"].items():
            print("   outcome", k, v)
        for k, v in report[d.name]["layer_scope"].items():
            print("   layer", k, v)
        for k, v in list(report[d.name]["pair"].items())[:8]:
            print("   pair", k, v)
        for k, v in report[d.name]["q0_bins"].items():
            print("   q0bin", k, v)
        print("   event q0", [round(x, 2) for x in q_event])
    out.write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
