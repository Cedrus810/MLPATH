import numpy as np
from ase import Atoms
from ase.build import molecule


def ethanol():
    return molecule("CH3CH2OH")


def chloride_methyl_chloride():
    # b04-like: CH3Cl plus a free Cl- on the back side, charge handled by the config
    return Atoms(
        "CClHHHCl",
        positions=[
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.785],
            [1.028, 0.0, -0.373],
            [-0.514, 0.89, -0.373],
            [-0.514, -0.89, -0.373],
            [0.0, 0.0, -3.6],
        ],
    )


def water_dimer():
    a = molecule("H2O")
    b = molecule("H2O")
    b.positions += [0.3, 0.2, 3.0]
    return a + b


def relabel(atoms, seed=0):
    """A random atom permutation; returns (permuted atoms, perm) with new[k] = old[perm[k]]."""
    perm = np.random.default_rng(seed).permutation(len(atoms))
    return atoms[perm], perm
