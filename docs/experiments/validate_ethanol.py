"""Independent frequency cross-check on the real system: ethanol minimum and both TSs."""
import sys; sys.path.insert(0, "/home/ruigengji/MLPATH/src")
import numpy as np
from ase.io import read
from ase.optimize import LBFGS
from prrs.calculators import mace_factory
from prrs.config import SearchConfig
from prrs.runner import curvature_spectrum, _wavenumbers
from prrs import validation as v
np.set_printoptions(precision=6, suppress=True, linewidth=170)

MODEL = "/home/kasuga/.cache/mace/MACE-OFF24_medium.model"
factory = lambda: mace_factory(MODEL, device="cuda", default_dtype="float64")
cfg = SearchConfig(quench_fmax_eV_A=0.001)
RUN = "/tmp/claude-1000/-home-ruigengji-MLPATH/bb9f534a-d589-44d0-a5ee-9f722d565b5e/scratchpad/p0_run3"

def check(atoms, label):
    ours_all, _, floor, prov = curvature_spectrum(atoms, factory, cfg, source="analytic")
    energies, trivial, _ = v.ase_energies_from_our_hessian(atoms, factory, cfg, "analytic")
    ours = v.internal_eigenvalues(ours_all, trivial)
    same = v.internal_eigenvalues(v.eigenvalues_from_energies(energies), trivial)
    print(f"\n=== {label} ===  fmax={np.abs(atoms.get_forces()).max():.2e}  "
          f"source={prov['source']}  trivial modes={trivial}")
    print(f"  ours   6 lowest: {ours[:6]}")
    print(f"  ASE (our Hessian) max|d| = {np.abs(ours - same).max():.3e}   "
          f"signs identical: {np.sign(ours).tolist() == np.sign(same).tolist()}")
    own, trivial2, prov2 = v.ase_vibrations_energies(atoms, factory, cfg)
    theirs = v.internal_eigenvalues(v.eigenvalues_from_energies(own), trivial2)
    print(f"  ASE own FD Hessian (delta={prov2['delta_A']}, {prov2['force_evaluations']} force evals)")
    print(f"     6 lowest: {theirs[:6]}")
    print(f"     max|d| vs our analytic = {np.abs(ours - theirs).max():.3e}")
    neg_ours = [x for x in ours if x < -1e-3]
    neg_theirs = [x for x in theirs if x < -1e-3]
    print(f"  negative modes: ours {len(neg_ours)}  ASE {len(neg_theirs)}")
    if neg_ours:
        print(f"     icm ours {_wavenumbers(neg_ours)}   icm ASE {_wavenumbers(neg_theirs)}")

atoms = read(f"{RUN}/chemical/c0000/m0000.extxyz") if False else None
minimum = read("/home/ruigengji/MLPATH/runs/e2e_ethanol/input.extxyz")
minimum.calc = factory()
LBFGS(minimum, logfile=None).run(fmax=1e-4, steps=500)
check(minimum, "ethanol minimum")
import glob
for path in sorted(glob.glob(f"{RUN}/ts_candidates/*.extxyz")):
    ts = read(path); ts.calc = factory()
    check(ts, f"TS {path.split('/')[-1]}")
