"""Does min-mode following reach the saddle from a dividing-surface seed?

Round 3 bracketed the surface to [0.25, 0.275] A on the O-H stretch and every trial there
still quenched into a basin. If climbing from those frames reaches lambda_1 = -37.3, the
remaining gap is plumbing, not method.
"""
import sys; sys.path.insert(0, "/home/ruigengji/MLPATH/src")
import numpy as np
from ase.io import read
from prrs.calculators import mace_factory
from prrs.config import SearchConfig
from prrs.runner import follow_min_mode, curvature_spectrum
factory = lambda: mace_factory("/home/kasuga/.cache/mace/MACE-OFF24_medium.model",
                               device="cuda", default_dtype="float64")
cfg = SearchConfig(quench_fmax_eV_A=0.002, min_mode_fmax_eV_A=0.01,
                   min_mode_step_A=0.05, min_mode_max_steps=200, min_mode_refresh=5)
reference = -7274.549375
for tag, path in (("0.250 elastic", "p1_run3/trials/t000025/raw_endpoint.extxyz"),
                  ("0.275 crossing", "p1_run3/trials/t000027/raw_endpoint.extxyz"),
                  ("0.350 transient", "p1_run3/trials/t000022/raw_endpoint.extxyz")):
    try:
        seed = read(path)
    except Exception as e:
        print(f"{tag}: no file ({e})"); continue
    seed.calc = factory()
    print(f"\n=== seed {tag} ===  E-Emin={(seed.get_potential_energy()-reference)*1000:.1f} meV  "
          f"fmax={np.abs(seed.get_forces()).max():.3f}")
    saddle, diag = follow_min_mode(seed, factory, cfg)
    if saddle is None:
        print(f"   FAILED: {diag.get('reason')}  steps={diag.get('steps')} "
              f"order={diag.get('saddle_order')} negs={diag.get('negative_eigenvalues')}")
        continue
    saddle.calc = factory()
    print(f"   saddle reached in {diag['steps']} steps  fmax={diag['fmax_eV_A']:.5f}")
    print(f"   order={diag['saddle_order']}  lambda={diag['negative_eigenvalues']}  "
          f"icm={[round(x,1) for x in diag['imaginary_wavenumbers_icm']]}")
    print(f"   barrier={(saddle.get_potential_energy()-reference)*1000:.2f} meV  "
          f"floor={diag['trivial_mode_floor_worst_residual']:.1e}")
    p = saddle.positions
    print(f"   O1-H={np.linalg.norm(p[0]-p[5]):.4f}  O2-H={np.linalg.norm(p[4]-p[5]):.4f}  "
          f"O1-O2={np.linalg.norm(p[0]-p[4]):.4f}")
