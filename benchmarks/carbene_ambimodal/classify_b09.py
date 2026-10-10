"""Which DFT saddle does each POLAR-1 b09 saddle correspond to?

b09 = dimethoxycarbene + 1-cyclohexene-1-carboxaldehyde. Four DFT saddles: TS-2OMe (ambimodal),
TS-2OMe-cis, TS-2OMe-C2, TS-2OMe-2 ((2+1) -> (4+1) rearrangement). Same discipline as
classify_b08.py: a candidate is named by its distance to *every* DFT saddle, not by its
starting point. RMSD is over heavy atoms only (methyl rotation permutes H and would inflate
an all-atom value), minimised over the mirror image; the DFT-DFT matrix gives the scale.

Markers, found from connectivity rather than fixed indices: the carbene carbon (the carbon
bonded to two methoxy oxygens), its distances to the two ring alkene carbons (Ca bears the
aldehyde, Cb is its neighbour) and to the aldehyde oxygen (the O-C bond of the (4+1) ring).

Usage: python benchmarks/carbene_ambimodal/classify_b09.py <label>=<prrs optimize dir> ...
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
from prrs.state import centred_kabsch_rmsd  # noqa: E402
from prrs.validation import internal_eigenvalues  # noqa: E402

EV_KCAL = 23.060548
DFT = ("TS-2OMe", "TS-2OMe-cis", "TS-2OMe-C2", "TS-2OMe-2")


def heavy_rmsd(a, b):
    heavy = [i for i, z in enumerate(a.numbers) if z > 1]
    x = a.positions[heavy] - a.positions[heavy].mean(axis=0)
    best = np.inf
    for sign in (1.0, -1.0):
        y = sign * b.positions[heavy]
        best = min(best, centred_kabsch_rmsd(x, y - y.mean(axis=0)))
    return float(best)


def atoms_of_interest(atoms):
    """Carbene C, ring alkene Ca (bearing CHO) and Cb, aldehyde O -- from the DFT TS-2OMe."""
    s = atoms.get_chemical_symbols()
    d = atoms.get_all_distances()
    oxygens = [i for i, x in enumerate(s) if x == "O"]
    carbons = [i for i, x in enumerate(s) if x == "C"]
    carbene = min(carbons, key=lambda c: sum(sorted(d[c, oxygens])[:2]))
    methoxy_o = sorted(oxygens, key=lambda o: d[carbene, o])[:2]
    aldehyde_o = next(o for o in oxygens if o not in methoxy_o)
    aldehyde_c = min((c for c in carbons if c != carbene), key=lambda c: d[c, aldehyde_o])
    ring_ca = min(
        (c for c in carbons if c not in (carbene, aldehyde_c)), key=lambda c: d[c, aldehyde_c]
    )
    hydrogens = [i for i, x in enumerate(s) if x == "H"]
    ring_cb = min(
        (c for c in carbons if c not in (carbene, aldehyde_c, ring_ca) and d[c, ring_ca] < 1.6),
        key=lambda c: sum(d[c, h] < 1.2 for h in hydrogens),  # the alkene CH, not a CH2
    )
    return dict(carbene=carbene, Ca=ring_ca, Cb=ring_cb, O_ald=aldehyde_o)


def markers(atoms, idx):
    d = atoms.get_all_distances()
    return {f"C:-{k}": float(d[idx["carbene"], idx[k]]) for k in ("Ca", "Cb", "O_ald")}


def main():
    s = HERE / "structures" / "ts"
    ref = {name: read(s / f"{name}.xyz") for name in DFT}
    idx = atoms_of_interest(ref["TS-2OMe"])
    print("atoms:", idx)
    for name, atoms in ref.items():
        m = markers(atoms, idx)
        print(f"DFT {name:14s} " + "  ".join(f"{k} {v:.3f}" for k, v in m.items()))
    print("DFT-DFT heavy-atom RMSD (A)")
    for (ka, a), (kb, b) in itertools.combinations(ref.items(), 2):
        print(f"  {ka:14s} {kb:14s} {heavy_rmsd(a, b):.3f}")
    calc = polar1()
    rows, geoms = [], {}
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
        rmsd = {k: heavy_rmsd(v, atoms) for k, v in ref.items()}
        row = dict(
            label=label,
            status=status,
            imaginary_cm=_wavenumbers(spectrum),
            E_eV=float(atoms.get_potential_energy()),
            markers=markers(atoms, idx),
            heavy_rmsd_to_dft=rmsd,
            closest_dft=min(rmsd, key=rmsd.get),
        )
        rows.append(row)
        geoms[label] = atoms
        print(
            f"{label:16s} {status:10s} imag {[round(v) for v in row['imaginary_cm']]}  "
            + "  ".join(f"{k} {v:.3f}" for k, v in row["markers"].items())
            + "  RMSD "
            + " / ".join(f"{k} {v:.3f}" for k, v in rmsd.items())
            + f"  -> {row['closest_dft']}"
        )
    print("\nbetween candidates: heavy-atom RMSD (A)")
    for (ka, a), (kb, b) in itertools.combinations(geoms.items(), 2):
        print(f"  {ka:16s} {kb:16s} {heavy_rmsd(a, b):.3f}")
    print(json.dumps(rows))


if __name__ == "__main__":
    main()
