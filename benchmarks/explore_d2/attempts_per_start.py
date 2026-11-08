"""Follow-up to guess_separability.py: near-identical L4 starting guesses give different
certify outcomes. For leader clusters of the starting guess at radius r (journal order),
allow at most m certify attempts per cluster (the first m members) and report forwards
skipped and certified events kept; also the per-cluster success rate of those clusters that
contain an event. Radii are multiples of the measured b05 geometry jitter: repeated
constrained relaxations of one scan point differ by up to 0.0057 A
(results/2026-10-25/hf_ef_rf_jitter_root_b05_fp32_sA_10-13.json, tau 3e-3).

Reads the guess_<run>.json/.npz written by guess_separability.py; no model.

With --merge-symmetric, events whose saddles are within basin_rmsd_A / basin_energy_eV under
the pose-graph automorphisms (hf_event_symmetry.json in CENSUS_DIR) count as one event: the
run's own duplicate rule used the smaller saddle-graph group and kept them apart.

    python attempts_per_start.py OUT_JSON CENSUS_DIR [--merge-symmetric] RUN_NAME [RUN_NAME ...]
"""

import json
import sys
from pathlib import Path

import numpy as np

JITTER_A = 0.0057
FOUND = ("c2", "duplicate_c2", "connected", "hf_duplicate")


def leaders(d, r):
    lead, of = [], []
    for i in range(len(d)):
        k = next((k for k in lead if d[i, k] <= r), None)
        if k is None:
            lead.append(i)
            k = i
        of.append(k)
    return of


def symmetric_classes(census_dir, name):
    for r in json.loads((census_dir / "hf_event_symmetry.json").read_text()):
        if Path(r["run"]).name != name:
            continue
        ev, d = r["events"], r["pairwise_symmetric_rmsd_A"]
        rep = {}
        for i, e in enumerate(ev):
            rep[e["task"]] = next(
                (
                    rep[ev[j]["task"]]
                    for j in range(i)
                    if d[i][j] <= r["basin_rmsd_A"]
                    and abs(e["energy_eV"] - ev[j]["energy_eV"]) <= r["basin_energy_eV"]
                ),
                e["task"],
            )
        return rep
    raise KeyError(name)


def analyse(census_dir, name, merge=False):
    items = json.loads((census_dir / f"guess_{name}.json").read_text())[0]["items"]
    if merge:
        rep = symmetric_classes(census_dir, name)
        for it in items:
            if it["event"]:
                it["event"] = rep[it["event"]]
    z = np.load(census_dir / f"guess_{name}.npz")
    d, spent = z["rmsd"], z["spent"]
    events = {it["event"] for it in items if it["event"]}
    out = {"run": name, "merged_symmetric": merge, "total_forwards": float(spent.sum()), "events_total": len(events), "radii": []}
    for mult in (2, 4, 10):
        r = mult * JITTER_A
        of = leaders(d, r)
        clusters = {}
        for i, k in enumerate(of):
            clusters.setdefault(k, []).append(i)
        with_event = [m for m in clusters.values() if any(items[i]["event"] for i in m)]
        rows = []
        for m_cap in (1, 2, 3, 5, 10):
            kept = [i for m in clusters.values() for i in m[:m_cap]]
            skipped = sorted(set(range(len(items))) - set(kept))
            rows.append(
                {
                    "attempts_per_cluster": m_cap,
                    "items_run": len(kept),
                    "forwards_skipped_fraction": float(spent[skipped].sum() / spent.sum()),
                    "events_kept": len({items[i]["event"] for i in kept if items[i]["event"]}),
                }
            )
        out["radii"].append(
            {
                "radius_A": r,
                "radius_over_jitter": mult,
                "clusters": len(clusters),
                "cluster_size_max": max(len(m) for m in clusters.values()),
                "clusters_with_event": [
                    {
                        "size": len(m),
                        "found": sum(items[i]["status"] in FOUND for i in m),
                        "first_found_at": next(
                            (n + 1 for n, i in enumerate(m) if items[i]["status"] in FOUND), None
                        ),
                    }
                    for m in with_event
                ],
                "policies": rows,
            }
        )
    return out


def main():
    out_path, census_dir, rest = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3:]
    merge = "--merge-symmetric" in rest
    names = [n for n in rest if n != "--merge-symmetric"]
    results = [analyse(census_dir, n, merge) for n in names]
    out_path.write_text(json.dumps(results, indent=1), encoding="utf-8")
    for r in results:
        print("==", r["run"], "events", r["events_total"])
        for x in r["radii"]:
            print(f"  r={x['radius_A']:.4f} ({x['radius_over_jitter']}x) clusters {x['clusters']} max size {x['cluster_size_max']}")
            print("    event clusters (size, found, first found at):", [(c["size"], c["found"], c["first_found_at"]) for c in x["clusters_with_event"]])
            for p in x["policies"]:
                print(f"    m={p['attempts_per_cluster']}: run {p['items_run']}, skip {100 * p['forwards_skipped_fraction']:.0f}% fw, events {p['events_kept']}/{r['events_total']}")


if __name__ == "__main__":
    main()
