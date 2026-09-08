"""Where does the FD Hessian's trivial-mode floor come from?

Displacing atom i by delta in a transverse direction changes the bond length only at
second order (chord vs arc), so the finite-difference transverse element is not zero but
~ E''(r) delta^2 / (2 r^2). Prediction: the fd floor scales as step^2 and the analytic
floor does not scale at all.
"""
import sys; sys.path.insert(0, "/home/ruigengji/MLPATH/src")
import numpy as np
from ase import Atoms
from prrs.calculators import double_well_factory
from prrs.config import SearchConfig
from prrs.runner import curvature_spectrum

def pair(d):
    a = Atoms("C2", positions=[[-d/2,0,0],[d/2,0,0]]); a.calc = double_well_factory(); return a

print(f"{'step (A)':>10} {'fd floor':>14} {'floor/step^2':>14} {'prediction':>14}")
r0, m = 1.2, 12.011
second = (1.6/0.36)*(3*1.0-1)
for step in (0.005, 0.01, 0.02, 0.04):
    cfg = SearchConfig(quench_fmax_eV_A=0.001, minimum_check_step_A=step)
    _,_,floor,_ = curvature_spectrum(pair(r0), double_well_factory, cfg, source="fd")
    worst = max(abs(v) for v in floor)
    pred = 2*second*step**2/(2*r0**2)/m
    print(f"{step:>10.3f} {worst:>14.3e} {worst/step**2:>14.4f} {pred:>14.3e}")
print()
print("ethanol-scale estimate with step=0.01, E''~50 eV/A^2, r~1.0 A, m=1 (H):")
print(f"  floor ~ {50*0.01**2/(2*1.0**2)/1.0:.2e} eV/(A^2 amu)   vs eigenvalue_tol 1e-3")
