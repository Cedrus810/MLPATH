"""MACE-POLAR-1-M on the DFT geometries of the SI (Wang & Xue, JACS 2026, 148, 8723).

First look, no optimisation: at each M06-2X/def2-TZVP-SMD geometry, POLAR-1 (gas phase) gives
an energy, the largest atomic force, and the projected Hessian spectrum (prrs.runner
.hessian_spectrum: analytic Hessian checked against a finite difference, rigid-body modes
projected out). A DFT stationary point need not be stationary on the MLP surface, so the
max force says how far off it is and the imaginary modes are read with that in mind.

Relative energies use reactant = olefin + singlet carbene, matched by formula, and compare
POLAR-1 electronic energies against the SI's DFT electronic energies (Tables S2/S3, which
include SMD solvation; the MLP has none).

Usage: python benchmarks/carbene_ambimodal/evaluate.py --out <dir> [--only NAME ...]
"""

import argparse
import json
import re
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np
from ase.io import read

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT)]
from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402
from prrs.runner import _trivial_modes, _wavenumbers, hessian_spectrum  # noqa: E402
from prrs.validation import internal_eigenvalues  # noqa: E402

HARTREE_KCAL = 627.509474
EV_KCAL = 23.060548
TABLE_ROW = re.compile(
    r"^\s*(\S.*?)\s+(-\d+\.\d+)\s+(-\d+\.\d+)\s+(-\d+\.\d+)(?:\s+(-?\d+\.\d+))?\s*$"
)


def dft_table():
    """Name -> electronic energy (Hartree) and imaginary mode, from SI Tables S2 and S3."""
    text = (HERE / "si" / "ja5c13906_si.txt").read_text().split("\n")
    start = next(i for i, line in enumerate(text) if line.strip().startswith("Table S2:"))
    end = next(i for i, line in enumerate(text) if line.strip() == "Structural information")
    rows = {}
    for line in text[start:end]:
        m = TABLE_ROW.match(line)
        if m:
            rows[m.group(1).strip()] = dict(
                E_hartree=float(m.group(2)),
                nu_cm=None if m.group(5) is None else float(m.group(5)),
            )
    return rows


def short_name(name):
    """Olefins are titled '<long name> <code>' in the coordinates but '<code>' in the tables."""
    return name.split()[-1]


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--only", nargs="*", help="structure names; default all")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    index = json.loads((HERE / "structures" / "index.json").read_text())
    if args.only:
        index = [row for row in index if row["name"] in args.only]
    dft = dft_table()
    calc = polar1()
    config = SearchConfig()
    rows = []
    for row in index:
        atoms = read(HERE / row["file"])
        atoms.calc = calc
        t0 = time.time()
        energy = atoms.get_potential_energy()
        fmax = float(np.linalg.norm(atoms.get_forces(), axis=1).max())
        # drop the projected rigid-body zeros before reading signs
        trivial = len(_trivial_modes(atoms.positions, atoms.get_masses()))
        spectrum = internal_eigenvalues(hessian_spectrum(atoms, lambda: calc, config), trivial)
        imaginary = _wavenumbers(spectrum)
        key = row["name"] if row["name"] in dft else short_name(row["name"])
        out = dict(
            name=row["name"],
            section=row["section"],
            formula=row["formula"],
            spin=row["spin"],
            E_mlp_eV=energy,
            fmax_mlp_eV_A=fmax,
            imaginary_mlp_cm=imaginary,
            lowest_eigenvalues=[float(v) for v in spectrum[:3]],
            E_dft_hartree=dft.get(key, {}).get("E_hartree"),
            nu_dft_cm=dft.get(key, {}).get("nu_cm"),
            seconds=round(time.time() - t0, 1),
        )
        rows.append(out)
        print(
            f"{row['name']:34s} fmax {fmax:7.3f} eV/A  imag {[round(v) for v in imaginary]}"
            f"  dft nu {out['nu_dft_cm']}  ({out['seconds']} s)",
            flush=True,
        )

    # reactant = olefin + singlet carbene with the same total formula
    olefins = [r for r in rows if r["section"] == "olefins"]
    carbenes = [r for r in rows if r["section"] == "carbenes" and r["spin"] == 1]
    for r in rows:
        if r["section"] not in ("ts", "products"):
            continue
        target = _count(r["formula"])
        pairs = [
            (o, c)
            for o in olefins
            for c in carbenes
            if _count(o["formula"]) + _count(c["formula"]) == target
        ]
        pairs = [p for p in pairs if p[0]["E_dft_hartree"] and p[1]["E_dft_hartree"]]
        if not pairs:
            continue
        o, c = min(pairs, key=lambda p: p[0]["E_dft_hartree"] + p[1]["E_dft_hartree"])
        r["reactants"] = [o["name"], c["name"]]
        r["dE_mlp_kcal"] = (r["E_mlp_eV"] - o["E_mlp_eV"] - c["E_mlp_eV"]) * EV_KCAL
        if r["E_dft_hartree"] is not None:
            r["dE_dft_kcal"] = (
                r["E_dft_hartree"] - o["E_dft_hartree"] - c["E_dft_hartree"]
            ) * HARTREE_KCAL
    (args.out / "polar1_on_dft_geometries.json").write_text(json.dumps(rows, indent=1))
    print("\nrelative to olefin + carbene (electronic, kcal/mol): MLP gas vs DFT SMD")
    for r in rows:
        if "dE_mlp_kcal" in r:
            dft_value = r.get("dE_dft_kcal")
            print(
                f"{r['name']:24s} MLP {r['dE_mlp_kcal']:8.2f}   DFT "
                f"{'   n/a' if dft_value is None else f'{dft_value:8.2f}'}   ({' + '.join(r['reactants'])})"
            )


def _count(formula):
    counts = Counter()
    for symbol, number in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        counts[symbol] += int(number or 1)
    return counts


if __name__ == "__main__":
    main()
