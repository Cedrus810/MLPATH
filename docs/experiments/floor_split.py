"""The trivial-mode floor conflates two things. Separate them.

Translations are exactly zero on any PES that depends on coordinate differences,
whatever the gradient. Rotations are zero only at a STATIONARY point: differentiating
E(R(s)x) = E(x) twice gives

    xdot^T H xdot = omega0^2 (g . x_perp)

so a rotational Rayleigh quotient is O(|g|) and says nothing about model quality until
that term is removed. Prediction tested against the measured floor here.
"""
import sys; sys.path.insert(0, "/home/ruigengji/MLPATH/src")
import numpy as np
from ase.io import read
from ase.optimize import LBFGS
from prrs.calculators import mace_factory
from prrs.config import SearchConfig
from prrs.runner import _trivial_modes, _mass_weighted_hessian

factory = lambda: mace_factory("/home/kasuga/.cache/mace/MACE-OFF24_medium.model",
                               device="cuda", default_dtype="float64")

def report(atoms, label):
    masses = atoms.get_masses(); pos = atoms.positions
    g = -atoms.get_forces()
    cfg = SearchConfig(quench_fmax_eV_A=0.001)
    Hm, basis, prov = _mass_weighted_hessian(atoms, factory, cfg, "analytic")
    root = np.sqrt(np.repeat(masses, 3))
    com = pos - np.average(pos, axis=0, weights=masses)
    J = sum(np.dot(r, r) * np.eye(3) - np.outer(r, r) for r in com)
    print(f"\n{label}: fmax = {np.abs(atoms.get_forces()).max():.3e} eV/A   "
          f"source = {prov['source']}")
    print(f"{'mode':>6} {'measured':>13} {'predicted from g':>18} {'residual':>13}")
    for i, mode in enumerate(basis):
        xdot = (mode / root).reshape(-1, 3)
        translation = np.linalg.norm(xdot - xdot.mean(axis=0)) < 1e-8 * np.linalg.norm(xdot)
        if translation:
            predicted = 0.0
        else:
            L = sum(np.cross(r, v) for r, v in zip(com, xdot))
            omega = np.linalg.pinv(J) @ L
            speed = np.linalg.norm(omega)
            axis = omega / speed
            perp = com - np.outer(com @ axis, axis)
            predicted = speed ** 2 * float(np.sum(g * perp)) / np.sum(mode ** 2)
        measured = float(mode @ Hm @ mode)
        print(f"{i:>6} {measured:>13.3e} {predicted:>18.3e} {measured - predicted:>13.3e}"
              f"  {'trans' if translation else 'rot'}")

atoms = read("/home/ruigengji/MLPATH/runs/e2e_ethanol/input.extxyz")
atoms.calc = factory()
report(atoms, "P0 source structure as stored")
LBFGS(atoms, logfile=None).run(fmax=1e-4, steps=500)
report(atoms, "same, relaxed to fmax 1e-4")
