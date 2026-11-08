"""Summarise the b05 QM chain (read-only, plain Python; run where the files are fresh).

Per stationary point: r2SCAN-3c final energy, lowest vibrational frequencies (imaginary
count), key distances; wB97M-V/def2-TZVPD single-point energies; relative energies in
kcal/mol. Atom indices: C0 central carbon, Cl1, O14, H15 (hydroxide H); the transferred H
and its beta carbon per TS: new_10-12 -> C10/H12, ctrl_6-9 -> C6/H9.

    python3 summarize.py
"""

import json
import math
import re
from pathlib import Path

HARTREE_KCAL = 627.509474
TS = {"new_10-12": (10, 12), "ctrl_6-9": (6, 9)}


def final_energy(path):
    vals = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", open(path).read())
    return float(vals[-1])


def frequencies(path):
    text = open(path).read()
    block = text.split("VIBRATIONAL FREQUENCIES")[-1].split("NORMAL MODES")[0]
    freqs = [float(x) for x in re.findall(r"^\s*\d+:\s+(-?\d+\.\d+) cm\*\*-1", block, re.M)]
    nonzero = [f for f in freqs if abs(f) > 1e-6]
    return nonzero


def xyz(path):
    lines = open(path).read().splitlines()
    n = int(lines[0])
    return [(l.split()[0], [float(v) for v in l.split()[1:4]]) for l in lines[2 : 2 + n]]


def dist(atoms, i, j):
    return math.dist(atoms[i][1], atoms[j][1])


out = {}
for name, (cb, h) in TS.items():
    rec = {}
    for tag, opt, sp in (
        ("ts", f"ts_{name}", f"wb_ts_{name}"),
        ("react", f"ep_{name}_react", f"wb_ep_{name}_react"),
        ("prod", f"ep_{name}_prod", f"wb_ep_{name}_prod"),
    ):
        # endpoints re-optimised to convergence by run_endpoint2.sh take precedence
        if tag != "ts" and Path(f"ep2_{name}_{tag}.out").exists():
            opt, sp = f"ep2_{name}_{tag}", f"wb_ep2_{name}_{tag}"
        converged = tag == "ts" or "OPTIMIZATION RUN DONE" in open(f"{opt}.out").read()
        geo = xyz(f"{opt}.xyz")
        fr = frequencies(f"{opt}.out")
        rec[tag] = {
            "source": opt,
            "opt_converged": converged,
            "E_r2scan": final_energy(f"{opt}.out"),
            "E_wb97mv": final_energy(f"{sp}.out"),
            "n_imag": sum(f < 0 for f in fr),
            "lowest_freqs": fr[:4],
            "C0-Cl1": dist(geo, 0, 1),
            "Cb-H": dist(geo, cb, h),
            "H-O14": dist(geo, h, 14),
            "C0-Cb": dist(geo, 0, cb),
            "O14-H15": dist(geo, 14, 15),
        }
    for level in ("E_r2scan", "E_wb97mv"):
        r = rec["react"][level]
        rec[f"barrier_{level}_kcal"] = (rec["ts"][level] - r) * HARTREE_KCAL
        rec[f"reaction_{level}_kcal"] = (rec["prod"][level] - r) * HARTREE_KCAL
    out[name] = rec

for level in ("E_r2scan", "E_wb97mv"):
    out[f"TS_new_minus_TS_ctrl_{level}_kcal"] = (
        out["new_10-12"]["ts"][level] - out["ctrl_6-9"]["ts"][level]
    ) * HARTREE_KCAL
    out[f"react_new_minus_react_ctrl_{level}_kcal"] = (
        out["new_10-12"]["react"][level] - out["ctrl_6-9"]["react"][level]
    ) * HARTREE_KCAL
    out[f"prod_new_minus_prod_ctrl_{level}_kcal"] = (
        out["new_10-12"]["prod"][level] - out["ctrl_6-9"]["prod"][level]
    ) * HARTREE_KCAL

open("summary.json", "w").write(json.dumps(out, indent=1) + "\n")
for name in TS:
    r = out[name]
    print(f"== {name}")
    for tag in ("react", "ts", "prod"):
        p = r[tag]
        print(
            f"  {tag:5s} {p['source']:22s} conv {p['opt_converged']!s:5s} imag {p['n_imag']} lowest {[round(f, 1) for f in p['lowest_freqs']]} | "
            f"C-Cl {p['C0-Cl1']:.2f} Cb-H {p['Cb-H']:.2f} H-O {p['H-O14']:.2f} "
            f"C0-Cb {p['C0-Cb']:.3f} O-H15 {p['O14-H15']:.2f}"
        )
    for level in ("E_r2scan", "E_wb97mv"):
        print(
            f"  {level}: barrier {r[f'barrier_{level}_kcal']:+.2f}  "
            f"reaction {r[f'reaction_{level}_kcal']:+.2f} kcal/mol"
        )
for k, v in out.items():
    if k.endswith("_kcal") and not isinstance(v, dict):
        print(f"{k}: {v:+.2f}")
