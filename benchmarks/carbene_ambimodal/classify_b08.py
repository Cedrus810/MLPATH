"""Which DFT saddle does each POLAR-1 b08 saddle correspond to?

The two DFT saddles are close (symmetry-aware RMSD 0.75 A): both attack C20=C16 at almost the
same distances and differ in the carbene's orientation -- C23-C17 is 3.73 A in TS-CCl2
(ambimodal, toward C4) and 4.23 A in TS-CCl2-C4 (away from C4). So "has one imaginary mode"
is not enough: every candidate is classified by the orientation markers and by the
symmetry-aware RMSD to *both* DFT saddles, and candidates are compared with each other.

Usage: python benchmarks/carbene_ambimodal/classify_b08.py <label>=<prrs optimize dir> ...
"""

import itertools
import json
import sys
from pathlib import Path

from ase.io import read

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT), str(HERE)]
from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402
from prrs.runner import _trivial_modes, _wavenumbers, hessian_spectrum  # noqa: E402
from prrs.validation import internal_eigenvalues  # noqa: E402
from scan_b08 import markers, rmsd_up_to_symmetry  # noqa: E402

EV_KCAL = 23.060548


def main():
    s = HERE / "structures" / "ts"
    ref = {"TS-CCl2": read(s / "TS-CCl2.xyz"), "TS-CCl2-C4": read(s / "TS-CCl2-C4.xyz")}
    for name, atoms in ref.items():
        m = markers(atoms)
        print(
            f"DFT {name:24s} C23-C20 {m['C23_C20']:.3f}  C23-C16 {m['C23_C16']:.3f}  C23-C17 {m['C23_C17']:.3f}"
        )
    calc = polar1()
    geoms, rows = {}, []
    for arg in sys.argv[1:]:
        label, directory = arg.split("=", 1)
        status = json.loads((Path(directory) / "result.json").read_text())["status"]
        atoms = read(Path(directory) / "geometry.extxyz")
        atoms.info.update(charge=0, spin=1)
        atoms.calc = calc
        trivial = len(_trivial_modes(atoms.positions, atoms.get_masses()))
        spectrum = internal_eigenvalues(
            hessian_spectrum(atoms, lambda: calc, SearchConfig()), trivial
        )
        row = dict(
            label=label,
            status=status,
            imaginary_cm=_wavenumbers(spectrum),
            E_eV=float(atoms.get_potential_energy()),
            **markers(atoms),
            rmsd_to_dft={k: rmsd_up_to_symmetry(v, atoms) for k, v in ref.items()},
        )
        row["closest_dft"] = min(row["rmsd_to_dft"], key=row["rmsd_to_dft"].get)
        rows.append(row)
        geoms[label] = atoms
        print(
            f"{label:28s} {status:10s} imag {[round(v) for v in row['imaginary_cm']]}"
            f"  C23-C20 {row['C23_C20']:.3f}  C23-C16 {row['C23_C16']:.3f}  C23-C17 {row['C23_C17']:.3f}"
            f"  RMSD to DFT TS-CCl2 {row['rmsd_to_dft']['TS-CCl2']:.3f} / TS-CCl2-C4 {row['rmsd_to_dft']['TS-CCl2-C4']:.3f}"
            f"  -> {row['closest_dft']}"
        )
    print("\nbetween candidates: symmetry-aware RMSD (A), energy difference (kcal/mol)")
    for (ka, a), (kb, b) in itertools.combinations(geoms.items(), 2):
        de = (b.get_potential_energy() - a.get_potential_energy()) * EV_KCAL
        print(f"  {ka:28s} {kb:28s} {rmsd_up_to_symmetry(a, b):.3f}   {de:+.3f}")
    print(json.dumps(rows))


if __name__ == "__main__":
    main()
