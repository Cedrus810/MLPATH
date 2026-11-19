"""Read-only (2026-11-15, user 「可以」): on an N1 SC run, which pushed pairs are the
Diels-Alder forming positions (diene termini C1 0 / C4 3 x dienophile C6 11 / C7 12), how far
their paths pushed rho, how they ended, and how the reactant-complex MINs place the diene
termini relative to the dienophile. No model calls.

    python n1_da_pairs.py RUN_DIR
"""

import collections
import glob
import json
import statistics as st
import sys

import numpy as np

DIENE, PHILE = (0, 3), (11, 12)
INTER = set(range(11))  # cyclopentadiene atoms 0..10; acrolein 11..18


def dist(x, i, j):
    return float(np.linalg.norm(np.asarray(x[i]) - np.asarray(x[j])))


def main(run):
    recs = []
    for f in sorted(glob.glob(f"{run}/wave_*.json")):
        recs += json.load(open(f))["paths"]
    kinds = collections.Counter()
    da, nodes = [], {}
    for r in recs:
        node, pair = r["label"].split("|")
        i, j = map(int, pair.split("-"))
        inter = (i in INTER) != (j in INTER)
        kinds["inter" if inter else "intra"] += 1
        st_min = (r.get("mins") or {}).get("start") or {}
        if st_min.get("positions") is not None and node not in nodes:
            x = st_min["positions"]
            nodes[node] = sorted(round(dist(x, a, b), 2) for a in DIENE for b in PHILE)
        if {i, j} & set(DIENE) and {i, j} & set(PHILE):
            en = (r.get("mins") or {}).get("end") or {}
            xe = en.get("positions")
            da.append({
                "label": r["label"], "status": r["status"], "steps": r["steps"],
                "a": r["a"], "b": r["b"],
                "rho0": round(r["rho_A"][0], 2), "rho_min": round(min(r["rho_A"]), 2),
                "rho_end": round(r["rho_A"][-1], 2),
                "dE_max_eV": round(max(r["energies_eV"]) - r["energies_eV"][0], 2),
                "end_min_CC": None if xe is None else
                    sorted(round(dist(xe, a, b), 2) for a in DIENE for b in PHILE),
                "bond_change": r["bond_change"],
            })  # fmt: skip
    print("paths", len(recs), dict(kinds))
    print("reactant-complex MINs:", len(nodes), "- sorted C1/C4 x C6/C7 distances (A):")
    for k, v in sorted(nodes.items(), key=lambda kv: kv[1][1]):
        print("  ", k, v)
    print("DA-position pairs pushed:", len(da))
    for d in da:
        print("  ", d)
    if da:
        print("rho_min median", st.median(d["rho_min"] for d in da),
              "status", collections.Counter(d["status"] for d in da))  # fmt: skip


if __name__ == "__main__":
    main(sys.argv[1])
