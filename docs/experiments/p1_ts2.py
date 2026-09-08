"""The symmetric ridge, found by constraining d(O1-H) - d(O2-H) = 0 and relaxing.

For a symmetric double well the constrained minimum on that surface is the saddle. The
naive Cartesian midpoint is not, because the O...O contraction is part of the reaction
coordinate and averaging two structures does not contract anything.
"""
import sys; sys.path.insert(0, "/home/ruigengji/MLPATH/src")
import numpy as np
from ase.io import read, write
from ase.optimize import BFGS
from ase.constraints import FixInternals
from prrs.calculators import mace_factory
from prrs.config import SearchConfig
from prrs.reliability import GuardedCalculator
from prrs.runner import curvature_spectrum, _wavenumbers, response_curvature, floor_residual

factory = lambda: mace_factory("/home/kasuga/.cache/mace/MACE-OFF24_medium.model",
                               device="cuda", default_dtype="float64")
cfg = SearchConfig(quench_fmax_eV_A=0.002)
A = read("p1_run/chemical/c0000_m0000.extxyz"); A.calc = factory()
eA = A.get_potential_energy()

ts = A.copy(); ts.calc = factory()
# start the proton halfway along the O...O line so the constraint has something to hold
o1, o2 = ts.positions[0], ts.positions[4]
ts.positions[5] = 0.5 * (o1 + o2)
combo = [0.0, [[0, 5, 1.0], [4, 5, -1.0]]]
ts.set_constraint(FixInternals(bondcombos=[combo]))
BFGS(ts, logfile=None).run(fmax=2e-3, steps=500)
ts.set_constraint()
p = ts.positions
print(f"constrained relax: O1-H={np.linalg.norm(p[0]-p[5]):.4f}  O2-H={np.linalg.norm(p[4]-p[5]):.4f}"
      f"  O1-O2={np.linalg.norm(p[0]-p[4]):.4f}")
eTS = ts.get_potential_energy()
print(f"  E={eTS:.6f}   barrier={(eTS-eA)*1000:.2f} meV   fmax(unconstrained)={np.abs(ts.get_forces()).max():.4f} eV/A")

lam, vec, floor, prov = curvature_spectrum(ts, factory, cfg, source="analytic")
neg = [x for x in lam if x < -cfg.minimum_check_eigenvalue_tol]
print(f"\nnegatives={len(neg)}  {[round(x,5) for x in neg]}   icm {[round(x,1) for x in _wavenumbers(neg)]}")
print(f"  floor residual {floor_residual(floor):.2e}   6 lowest {np.array2string(lam[:8], precision=5)}")
if len(neg) == 1:
    masses = ts.get_masses()
    mode = vec[:, 0].reshape(-1, 3) / np.sqrt(masses[:, None]); mode /= np.linalg.norm(mode)
    guard = GuardedCalculator(factory(), cfg, eTS)
    probe = ts.copy(); probe.calc = guard
    rep = response_curvature(probe, guard, mode, cfg)
    print(f"  kappa={rep['kappa']:+.6f} +- {rep['sigma']:.2e} c={rep['anharmonicity']:+.3f} "
          f"sign_stable={rep['sign_stable']} significant={rep['significant']}")
    weights = sorted(((float(np.linalg.norm(mode[i])), i) for i in range(len(ts))), reverse=True)
    print(f"  unstable mode weight by atom (top 3): {[(i, round(w,3)) for w, i in weights[:3]]}"
          f"   (atom 5 is the transferring proton)")
    out = ts.copy(); out.calc = None; write("p1_ts.extxyz", out)
