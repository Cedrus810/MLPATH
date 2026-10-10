"""Read the POLAR-1 geomeTRIC re-optimisations of optimize_b08_b09.sh against the SI's DFT geometries.

Per structure: optimiser status, projected-Hessian imaginary modes at the optimised geometry
(prrs.runner.hessian_spectrum, rigid-body zeros dropped), Kabsch RMSD to the DFT starting
geometry (same atom order, no permutation search), and the carbene carbon's two closest
carbons -- the forming C-C bonds that tell a (2+1) from a (4+1) approach. Then energies
relative to the optimised olefin + carbene, and the RMSD between pairs of saddles that must
be different points (TS-CCl2 vs TS-CCl2-C4; the b09 saddles among themselves).

Usage: python benchmarks/carbene_ambimodal/analyze_opt.py <opt dir>
"""

import itertools
import json
import sys
from pathlib import Path

import numpy as np
from ase.io import read

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT)]
from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402
from prrs.runner import _trivial_modes, _wavenumbers, hessian_spectrum  # noqa: E402
from prrs.state import aligned_rmsd  # noqa: E402
from prrs.validation import internal_eigenvalues  # noqa: E402

EV_KCAL = 23.060548
REACTANTS = {
    "b08": ("1_2-bis_methylene_cycloheptane_7C", "CCl2"),
    "b09": ("1-Cyclohexene-1-carboxaldehyde_6O", "Dimethoxycarbene"),
}
SYSTEM = {"TS-CCl2": "b08", "TS-CCl2-C4": "b08"}


def carbene_carbon(atoms):
    """The carbon carrying both halogens / both methoxy oxygens."""
    symbols = atoms.get_chemical_symbols()
    partners = [i for i, s in enumerate(symbols) if s in ("Cl", "Br", "F", "I")]
    if not partners:  # dimethoxycarbene: the only carbon with two O at bonding distance
        partners = [i for i, s in enumerate(symbols) if s == "O"]
    carbons = [i for i, s in enumerate(symbols) if s == "C"]
    d = atoms.get_all_distances()
    return min(carbons, key=lambda c: sum(sorted(d[c, partners])[:2]))


def main():
    opt = Path(sys.argv[1])
    calc = polar1()
    config = SearchConfig()
    rows, geoms = {}, {}
    for directory in sorted(p for p in opt.iterdir() if p.is_dir()):
        name = directory.name
        result = json.loads((directory / "result.json").read_text())
        start = read(directory / "input.extxyz")
        final = read(directory / "geometry.extxyz")
        final.info.update(charge=0, spin=1)
        final.calc = calc
        energy = final.get_potential_energy()
        trivial = len(_trivial_modes(final.positions, final.get_masses()))
        spectrum = internal_eigenvalues(hessian_spectrum(final, lambda: calc, config), trivial)
        c = carbene_carbon(final)
        d = final.get_all_distances()
        others = sorted(
            (float(d[c, j]), j)
            for j, s in enumerate(final.get_chemical_symbols())
            if s in "CO" and j != c
        )
        rows[name] = dict(
            status=result.get("status"),
            E_eV=energy,
            fmax_eV_A=float(np.linalg.norm(final.get_forces(), axis=1).max()),
            imaginary_cm=_wavenumbers(spectrum),
            rmsd_to_dft_A=aligned_rmsd(start, final),
            carbene_carbon=c,
            closest_to_carbene=[(round(r, 3), j, final[j].symbol) for r, j in others[:5]],
        )
        geoms[name] = final
        print(
            f"{name:36s} {rows[name]['status']:10s} imag {[round(v) for v in rows[name]['imaginary_cm']]}"
            f"  fmax {rows[name]['fmax_eV_A']:.4f}  rmsd_dft {rows[name]['rmsd_to_dft_A']:.3f}"
            f"  C{c}... {rows[name]['closest_to_carbene']}",
            flush=True,
        )
    print("\nenergies relative to optimised olefin + carbene (POLAR-1, gas, kcal/mol)")
    for name, row in rows.items():
        system = SYSTEM.get(name, "b09" if ("OMe" in name or name.startswith("6O")) else None)
        if system is None or name in REACTANTS[system]:
            continue
        olefin, carbene = REACTANTS[system]
        if olefin in rows and carbene in rows:
            row["dE_kcal"] = (
                row["E_eV"] - rows[olefin]["E_eV"] - rows[carbene]["E_eV"]
            ) * EV_KCAL
            print(f"{name:20s} {system}  {row['dE_kcal']:8.2f}")
    print("\nRMSD between saddles that must be distinct points (A)")
    saddles = [n for n in rows if n.startswith("TS-")]
    for a, b in itertools.combinations(saddles, 2):
        if len(geoms[a]) == len(geoms[b]):
            print(f"{a:14s} {b:14s} {aligned_rmsd(geoms[a], geoms[b]):.3f}")
    (opt.parent / "opt_analysis.json").write_text(json.dumps(rows, indent=1))


if __name__ == "__main__":
    main()
