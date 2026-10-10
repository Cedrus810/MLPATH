"""Starting structure for a b09 PRRS probe: the reactant-side complex below TS-2OMe on POLAR-1.

From the POLAR-1 TS-2OMe (results/<date>/opt/TS-2OMe/geometry.extxyz), step along the
imaginary mode in the direction that lengthens the forming carbene-C : Cb bond, then relax
with BFGS on POLAR-1. The result is the pre-reaction complex the paper's trajectories start
from ("vdW complex"). It is the complex that feeds TS-2OMe, so a search from it is biased
toward that channel -- acceptable for a first probe, to be revisited for a real benchmark.

Checks printed: carbene-C distances to Ca / Cb / aldehyde O, and the projected Hessian
(a minimum has no imaginary mode). Atom order is the SI's TS-2OMe order throughout
(carbene C = 2, Ca = 13, Cb = 15, aldehyde O = 1; see classify_b09.py).

Usage: python benchmarks/carbene_ambimodal/build_b09_input.py <TS geometry> <output extxyz>
"""

import sys
from pathlib import Path

import numpy as np
from ase.io import read, write
from ase.optimize import BFGS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT)]
from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402
from prrs.runner import _trivial_modes, _wavenumbers, hessian_spectrum  # noqa: E402
from prrs.validation import internal_eigenvalues  # noqa: E402

CARBENE, CA, CB, O_ALD = 2, 13, 15, 1
STEP_A = 0.15  # largest atomic displacement of the initial push along the mode


def spectrum(atoms, calc):
    trivial = len(_trivial_modes(atoms.positions, atoms.get_masses()))
    return internal_eigenvalues(hessian_spectrum(atoms, lambda: calc, SearchConfig()), trivial)


def main():
    ts = read(sys.argv[1])
    ts.info.update(charge=0, spin=1)
    calc = polar1()
    ts.calc = calc
    values, vectors = hessian_spectrum(ts, lambda: calc, SearchConfig(), return_vectors=True)
    mode = (vectors[:, 0] / np.sqrt(np.repeat(ts.get_masses(), 3))).reshape(-1, 3)
    mode *= STEP_A / np.linalg.norm(mode, axis=1).max()
    pushed = []
    for sign in (1.0, -1.0):
        trial = ts.copy()
        trial.positions += sign * mode
        pushed.append((trial.get_distance(CARBENE, CB), trial))
    atoms = max(pushed, key=lambda p: p[0])[1]  # the side where the forming bond lengthens
    atoms.info.update(charge=0, spin=1)
    atoms.calc = calc
    BFGS(atoms, logfile=None).run(fmax=0.005, steps=2000)
    d = atoms.get_all_distances()
    print(f"lowest TS eigenvalue {values[0]:.4f} eV/A^2/amu")
    print(
        f"complex: C:-Ca {d[CARBENE, CA]:.3f}  C:-Cb {d[CARBENE, CB]:.3f}  C:-O_ald {d[CARBENE, O_ALD]:.3f}"
    )
    print(f"max force {np.linalg.norm(atoms.get_forces(), axis=1).max():.4f} eV/A")
    print(f"imaginary modes {[round(v) for v in _wavenumbers(spectrum(atoms, calc))]}")
    print(f"E {atoms.get_potential_energy():.5f} eV  (TS {ts.get_potential_energy():.5f})")
    atoms.calc = None
    atoms.info = dict(charge=0, spin=1, name="b09 reactant complex (POLAR-1, below TS-2OMe)")
    write(sys.argv[2], atoms, format="extxyz")


if __name__ == "__main__":
    main()
