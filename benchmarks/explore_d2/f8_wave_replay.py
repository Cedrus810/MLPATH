"""Offline estimate for F8 (design record §4 F8): when would the D2-HF-v1.5 runs have hit
their references if exploration and refinement alternated in waves? Read-only; no model.

Every task keeps the forward cost its run actually paid (journal settles); only the order of
payment changes. Wave w = the w-th of W equal slices of the L1 tasks in the run's own L1
order (the frozen population order, round-robin over poses). After a wave's L1 (and its L2),
the wave's first-group items (not second-round start-cluster members, not same-graph
ridges) pay HS then HC; GN is paid once, at the first wave holding a bond-changing ridge.
L3 releases, the second round and the same-graph group are paid after the last wave, as now
(so the L1 coverage is unchanged). A reference's hit cost is the cumulative spend when the
first HC that reached it (connected, or hf_duplicate of a connected event matched to it)
settles. Assumes outcomes do not change with order: duplicate labels, start clusters and fp32
noise can, so this is an estimate of the ordering effect only.

    python f8_wave_replay.py OUT_JSON
"""

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
NEW = ROOT / "benchmarks/explore_d2/results/2026-11-13"
WAVES = (1, 2, 4, 8, 16)


def load(run):
    order, cost, payload, status = [], {}, {}, {}
    with open(run / "journal.jsonl", encoding="utf-8") as fh:
        for line in fh:
            ev = json.loads(line)
            if ev["event"] in ("settle", "reconcile"):
                t = ev["task"]
                cost[t] = cost.get(t, 0) + ev["spent"]
                if ev["event"] == "settle":
                    order.append(t)
                    status[t] = ev["status"]
                    if t.startswith(("HS|", "HC|")):
                        p = ev["payload"]
                        payload[t] = {
                            k: p.get(k)
                            for k in ("same_graph_ridge", "start_cluster", "class", "duplicate_of")
                        }
    return order, cost, payload, status


def replay(run, hits_by_task, n_waves, first_round=1, l3_first=False):
    """l3_first: the current order (all L3 before any HS, as _promote does) -- a check that
    the replay reproduces the run's own first hits when n_waves = 1."""
    order, cost, payload, status = load(run)
    key = lambda t: t.split("|", 1)[1]  # noqa: E731
    l1 = [t for t in order if t.startswith("L1|")]
    wave_of = {key(t): min(i * n_waves // len(l1), n_waves - 1) for i, t in enumerate(l1)}
    first, late = {w: [] for w in range(n_waves)}, []
    for t in order:
        if not t.startswith("HS|"):
            continue
        p = payload[t]
        sc = p.get("start_cluster") or {}
        deferred = p.get("same_graph_ridge") or sc.get("rank", 0) >= first_round
        (late if deferred else first[wave_of[key(t)]]).append(t)
    gn = sum(c for t, c in cost.items() if t.startswith("GN|"))
    spent, hit = 0, {}
    gn_paid = gn == 0

    def pay_hs(t):
        nonlocal spent
        spent += cost[t]
        hc = "HC|" + key(t)
        if hc in cost:
            spent += cost[hc]
            ref = hits_by_task.get(hc)
            if ref and ref not in hit:
                hit[ref] = spent

    for w in range(n_waves):
        for t in order:
            if t[:3] in ("L1|", "L2|") and wave_of.get(key(t)) == w:
                spent += cost[t]
        if l3_first and w == n_waves - 1:
            spent += gn
            gn_paid = True
            for t in order:
                if t.startswith("L3|"):
                    spent += cost[t]
        bond_changing = [
            t for t in first[w] if payload[t].get("class") == "ridge_suspected"
        ]
        if not gn_paid and bond_changing:
            spent += gn
            gn_paid = True
        for t in first[w]:
            pay_hs(t)
    if not gn_paid:
        spent += gn
    if not l3_first:
        for t in order:
            if t.startswith("L3|"):
                spent += cost[t]
    for t in late:
        pay_hs(t)
    return hit, spent


def main():
    out = Path(sys.argv[1])
    val = json.loads((NEW / "v15_validate.json").read_text())
    res = {}
    for name in ("b04", "b05", "n1"):
        for tag in ("sA", "sB"):
            run = NEW / f"{name}_v15_{tag}"
            v = val[name][tag]
            # reference reached by each connected HC (order 1 and matched), then duplicates
            ref_of = {e["task"]: e["match"] for e in v["events"] if e["order"] == 1 and e["match"]}
            _, _, payload, status = load(run)
            hits = dict(ref_of)
            for t, p in payload.items():
                if t.startswith("HC|") and status[t] == "hf_duplicate":
                    if p.get("duplicate_of") in ref_of:
                        hits[t] = ref_of[p["duplicate_of"]]
            check, _ = replay(run, hits, 1, l3_first=True)
            row = {"B_star": v["budget"], "actual_first_hit": v["first_hit"],
                   "replay_current_order": check, "waves": {}}
            for n in WAVES:
                hit, total = replay(run, hits, n)
                row["waves"][n] = {"first_hit": hit, "total": total}
            res[f"{name}_{tag}"] = row
            print(name, tag, "B*", v["budget"], "actual", v["first_hit"], "replay", check)
            for n in WAVES:
                h = row["waves"][n]["first_hit"]
                print("  W", n, {k: (c, round(c / v["budget"], 3)) for k, c in sorted(h.items())},
                      "total", row["waves"][n]["total"])
    out.write_text(json.dumps(res, indent=1) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
