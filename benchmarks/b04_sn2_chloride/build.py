"""Build the b04 input: Cl- approaching CH3Cl from the backside, one fragment apart.

Deterministic and offline. The geometry here is a STARTING POINT, not an answer: the
search relaxes it and decides for itself what is a basin. In particular the separation
below is chosen to sit inside the proposal layer's reach, not at any structure taken
from the literature -- handing the search a pre-formed ion-dipole complex would be
handing it the product.

Why 3.6 A: `pair_cutoff_A = 4.0` is the largest separation at which the pair-approach
family will propose anything at all (prrs/perturbations.py, the non-bonded branch), and
the pair-hop filter does not apply across fragments because there is no bond path
between them. Closer than about 3.2 A and the starting point is already the complex.

Run:  python benchmarks/b04_sn2_chloride/build.py [out.extxyz]
"""

import sys
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io import write

# Total charge for the whole system. It is recorded here as well as in config.json
# because a structure file that does not carry its charge is a different chemical
# object read by a charge-aware model -- and MACE-OFF24 would return the same energy
# for all three charges, which is exactly the trap this benchmark exists to avoid.
TOTAL_CHARGE = -1
MULTIPLICITY = 1

# Neutral CH3Cl, near its own equilibrium; r(C-Cl) = 1.785 A is the starting guess the
# SN2 probe used and relaxed to 1.7816 A on both charge-aware models.
_METHYL_CHLORIDE = Atoms(
    "CClH3",
    positions=[
        [0.000, 0.000, 0.000],
        [0.000, 0.000, 1.785],
        [1.028, 0.000, -0.373],
        [-0.514, 0.890, -0.373],
        [-0.514, -0.890, -0.373],
    ],
)

# Backside, along -z, opposite the leaving chlorine. Placed on the C3 axis because any
# other placement breaks the symmetry the D3h analysis in the probe relies on, and a
# broken symmetry is a different starting point, not a more general one.
APPROACH_A = 3.6


def build(separation_A=APPROACH_A):
    atoms = _METHYL_CHLORIDE.copy()
    atoms += Atoms("Cl", positions=[[0.0, 0.0, -float(separation_A)]])
    atoms.info["charge"] = TOTAL_CHARGE
    atoms.info["spin"] = MULTIPLICITY
    # No cell, no pbc: the first version is finite and non-periodic by scope
    # (PRRS_ENGINEERING_PLAN.md section 4).
    atoms.pbc = False
    return atoms


def main(argv):
    atoms = build()
    carbon, leaving, nucleophile = 0, 1, 5
    print(f"formula          {atoms.get_chemical_formula()}")
    print(f"total charge     {atoms.info['charge']}   multiplicity {atoms.info['spin']}")
    print(f"r(C-Cl) leaving  {atoms.get_distance(carbon, leaving):.4f} A")
    print(f"r(C...Cl-) nuc   {atoms.get_distance(carbon, nucleophile):.4f} A")
    print(f"fragments        2 (the search must be able to reach 1)")
    if np.linalg.norm(atoms.positions[nucleophile] - atoms.positions[carbon]) > 4.0:
        raise SystemExit("nucleophile is outside pair_cutoff_A = 4.0; no probe would be proposed")
    if len(argv) > 1:
        write(str(Path(argv[1])), atoms)
        print(f"written          {argv[1]}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
