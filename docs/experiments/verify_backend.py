"""Does the new analytic path actually engage on MACE, and what does it cost?"""
import sys, time; sys.path.insert(0, "/home/ruigengji/MLPATH/src")
import numpy as np
from ase.build import molecule
from ase.optimize import LBFGS
from prrs.calculators import mace_factory
from prrs.config import SearchConfig
from prrs.runner import curvature_spectrum, response_curvature
from prrs.reliability import GuardedCalculator

MODEL = "/home/kasuga/.cache/mace/MACE-OFF24_medium.model"
factory = lambda: mace_factory(MODEL, device="cuda", default_dtype="float64")
atoms = molecule("CH3CH2OH"); atoms.calc = factory()
LBFGS(atoms, logfile=None).run(fmax=1e-4, steps=500)
print(f"minimum fmax = {np.abs(atoms.get_forces()).max():.2e} eV/A")

for source in ("analytic", "fd"):
    cfg = SearchConfig(quench_fmax_eV_A=0.001, hessian_source=source)
    t = time.perf_counter()
    lam, vec, floor, prov = curvature_spectrum(atoms, factory, cfg, source=source)
    dt = time.perf_counter() - t
    worst = max(abs(v) for v in floor)
    gate = cfg.trivial_floor_fraction * cfg.minimum_check_eigenvalue_tol
    print(f"\nsource={source:<9} {dt:7.2f} s   provenance={ {k:v for k,v in prov.items() if k!='meaning'} }")
    print(f"  trivial floor worst = {worst:.3e}   gate threshold = {gate:.3e}   "
          f"{'PASS' if worst <= gate else 'FAIL (closed)'}")
    print(f"  6 lowest physical-ish eigenvalues: {np.array2string(lam[6:12], precision=5)}")

cfg = SearchConfig(quench_fmax_eV_A=0.001, response_band_A=(0.01, 0.1), response_samples=5)
lam, vec, _, _ = curvature_spectrum(atoms, factory, cfg, source="analytic")
masses = atoms.get_masses()
mode = vec[:, 6].reshape(-1, 3) / np.sqrt(masses[:, None])  # cartesian displacement
guard = GuardedCalculator(factory(), cfg, None)
probe = atoms.copy(); probe.calc = guard
rep = response_curvature(probe, guard, mode, cfg)
print(f"\nlowest physical mode: lambda = {lam[6]:+.6f}")
print(f"  band regression kappa = {rep['kappa']:+.6f} +- {rep['sigma']:.2e}   "
      f"c = {rep['anharmonicity']:+.3f}")
print(f"  significant={rep['significant']}  sign_stable={rep['sign_stable']}  "
      f"residual_rms={rep['residual_rms']:.2e}")
print(f"  kappa - lambda = {rep['kappa']-lam[6]:+.3e}")
