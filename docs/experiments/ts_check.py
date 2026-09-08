"""Does kappa_a along the reported unstable mode converge to lambda_1 as a -> 0?"""
import sys; sys.path.insert(0, "/home/ruigengji/MLPATH/src")
import numpy as np
from ase.io import read
from ase import units
from prrs.calculators import mace_factory
from prrs.config import SearchConfig
from prrs.reliability import GuardedCalculator
from prrs.runner import curvature_spectrum, response_curvature, _mass_weighted_direction

factory = lambda: mace_factory("/home/kasuga/.cache/mace/MACE-OFF24_medium.model",
                               device="cuda", default_dtype="float64")
atoms = read("/tmp/claude-1000/-home-ruigengji-MLPATH/bb9f534a-d589-44d0-a5ee-9f722d565b5e/scratchpad/p0_run2/ts_candidates/ts0000.extxyz")
atoms.calc = factory()
cfg = SearchConfig(quench_fmax_eV_A=0.001)
lam, vec, floor, prov = curvature_spectrum(atoms, factory, cfg, source="analytic")
masses = atoms.get_masses()
print(f"fmax = {np.abs(atoms.get_forces()).max():.3e}   source = {prov['source']}")
print(f"lowest 8 eigenvalues: {np.array2string(lam[:8], precision=6)}")
lam1 = lam[0]
scale = units._hbar*1e10/np.sqrt(units._e*units._amu)
print(f"lambda_1 = {lam1:+.8f}   as icm = {-scale*np.sqrt(-lam1)/units.invcm:+.4f}")

v1 = vec[:, 0]
cartesian = (v1.reshape(-1, 3) / np.sqrt(masses[:, None]))
cartesian /= np.linalg.norm(cartesian)
u = _mass_weighted_direction(masses, cartesian)
print(f"round trip |u . v1| = {abs(float(u.ravel() @ v1)):.12f}   (1.0 means exact)")

guard = GuardedCalculator(factory(), cfg, None)
probe = atoms.copy(); probe.calc = guard
E0 = float(guard.get_potential_energy(probe))
ref = probe.positions.copy()
disp = u / np.sqrt(masses[:, None])
print(f"\n{'a':>10} {'kappa_E':>13} {'kappa_E - lambda_1':>20}")
for a in (1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1):
    es = []
    for s in (1, -1):
        probe.set_positions(ref + s*a*disp)
        es.append(float(guard.get_potential_energy(probe)))
    probe.set_positions(ref)
    kap = (es[0]+es[1]-2*E0)/a**2
    print(f"{a:>10.0e} {kap:>13.6f} {kap-lam1:>20.3e}")

for band in [(0.001, 0.01), (0.01, 0.1)]:
    c = SearchConfig(quench_fmax_eV_A=0.001, response_band_A=band, response_samples=5)
    r = response_curvature(probe, guard, cartesian, c)
    print(f"band {band}: kappa = {r['kappa']:+.6f} +- {r['sigma']:.2e}  "
          f"c = {r['anharmonicity']:+.3f}  error vs lambda_1 = {r['kappa']-lam1:+.3e}")
