"""Read-only (2026-11-18, user 「可以」): why were the N1 s-cis DA saddles not found in PBS 29729?
Double-pair coverage, DA-position combinations (run / not run, outcome), and the acrolein
O9=C8-C7=C6 dihedral of every reactant complex and of the DA paths. No model calls.

    python n1_scis_check.py RUN_DIR [RUN_DIR ...]
"""

import json
import glob
import sys
import math
import numpy as np

DA = [{(0, 11), (3, 12)}, {(0, 12), (3, 11)}]


def dihed(x, a, b, c, d):
    x = np.asarray(x)
    b0, b1, b2 = x[a] - x[b], x[c] - x[b], x[d] - x[c]
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    return math.degrees(math.atan2(np.dot(np.cross(b1n, v), w), np.dot(v, w)))


for run in sys.argv[1:]:
    m = json.load(open(run + "/manifest.json"))
    recs = []
    for f in sorted(glob.glob(run + "/dwave_*.json")):
        recs += json.load(open(f))["paths"]
    ran = [r for r in recs if r["status"] != "budget_exhausted"]
    print(
        "==",
        run.split("/")[-1],
        "double planned",
        m["double_paths_planned"],
        "in wave files",
        len(recs),
        "ran (not budget_exhausted)",
        len(ran),
    )
    da = []
    for r in recs:
        p = r["label"].split("|")[1].split("+")
        prs = {tuple(sorted(map(int, q.split("-")))) for q in p}
        if prs in DA:
            da.append(r)
    print(
        "   DA-position combos in wave files:",
        len(da),
        " ran:",
        sum(r["status"] != "budget_exhausted" for r in da),
    )
    for r in da:
        if r["status"] == "budget_exhausted" or not r["rho_A"]:
            print("   ", r["label"], "budget_exhausted (not run)")
            continue
        mins = r.get("mins") or {}
        st = mins.get("start") or {}
        en = mins.get("end") or {}
        dst = (
            dihed(st["positions"], 14, 13, 12, 11) if st.get("positions") is not None else None
        )
        den = (
            dihed(en["positions"], 14, 13, 12, 11) if en.get("positions") is not None else None
        )
        sts = set(map(tuple, st.get("edges") or []))
        ens = set(map(tuple, en.get("edges") or []))
        ts = [
            (t.get("status"), round(t["energy_eV"], 3) if t.get("energy_eV") else None)
            for t in r.get("ts", [])
        ]
        print(
            "   ",
            r["label"],
            r["status"],
            r["steps"],
            "q",
            round(r["rho_A"][0], 2),
            "->",
            round(min(r["rho_A"]), 2),
            "O=C-C=C dihedral start",
            None if dst is None else round(dst),
            "end",
            None if den is None else round(den),
            "formed",
            sorted(ens - sts),
            "ts",
            ts,
        )

    nodes = {}
    for f in sorted(glob.glob(run + "/wave_*.json")):
        for r in json.load(open(f))["paths"]:
            n = r["label"].split("|")[0]
            st = (r.get("mins") or {}).get("start") or {}
            if n not in nodes and st.get("positions") is not None:
                nodes[n] = round(dihed(st["positions"], 14, 13, 12, 11))
    vals = sorted(nodes.values(), key=abs)
    print(
        "   acrolein O=C-C=C dihedral over",
        len(nodes),
        "reactant complexes (path-start MINs): |d| min",
        abs(vals[0]),
        "max",
        abs(vals[-1]),
        "; |d| < 90:",
        sum(abs(v) < 90 for v in vals),
    )
