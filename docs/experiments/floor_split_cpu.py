"""Validate the gradient prediction for rotational Rayleigh quotients on CPU.

E(R(s)x) = E(x) for all s, so d2/ds2 gives  xdot^T H xdot = omega^2 (g . x_perp).
Translations have no such term. Tested at non-stationary geometries where the
contamination is large, and at the minimum where it must vanish.
"""
import sys; sys.path.insert(0, "/home/ruigengji/MLPATH/src"); sys.path.insert(0, "/home/ruigengji/MLPATH/tests")
import numpy as np
from ase import Atoms
from prrs.calculators import double_well_factory
from prrs.config import SearchConfig
from prrs.runner import _mass_weighted_hessian

def predict(mode, root, com, gradient, inertia):
    xdot = (mode / root).reshape(-1, 3)
    if np.linalg.norm(xdot - xdot.mean(axis=0)) < 1e-8 * max(np.linalg.norm(xdot), 1e-30):
        return 0.0, "trans"
    L = sum(np.cross(r, v) for r, v in zip(com, xdot))
    omega = np.linalg.pinv(inertia) @ L
    speed = np.linalg.norm(omega)
    if speed < 1e-12:
        return 0.0, "trans"
    axis = omega / speed
    perp = com - np.outer(com @ axis, axis)
    return speed ** 2 * float(np.sum(gradient * perp)) / float(np.sum(mode ** 2)), "rot"

cfg = SearchConfig(quench_fmax_eV_A=0.001)
for label, distance in [("non-stationary r=1.5", 1.5), ("non-stationary r=2.0", 2.0),
                        ("minimum r=1.2", 1.2), ("barrier top r=1.8", 1.8)]:
    atoms = Atoms("CH", positions=[[-distance/2, 0, 0], [distance/2, 0, 0]])
    atoms.calc = double_well_factory()
    masses = atoms.get_masses(); root = np.sqrt(np.repeat(masses, 3))
    gradient = -atoms.get_forces()
    com = atoms.positions - np.average(atoms.positions, axis=0, weights=masses)
    inertia = sum(np.dot(r, r)*np.eye(3) - np.outer(r, r) for r in com)
    Hm, basis, _ = _mass_weighted_hessian(atoms, double_well_factory, cfg, "fd")
    print(f"\n{label}   fmax = {np.abs(gradient).max():.4f} eV/A")
    print(f"{'kind':>6} {'measured':>13} {'predicted':>13} {'residual':>13}")
    for mode in basis:
        p, kind = predict(mode, root, com, gradient, inertia)
        m = float(mode @ Hm @ mode)
        print(f"{kind:>6} {m:>13.3e} {p:>13.3e} {m-p:>13.3e}")
