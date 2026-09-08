"""Which bonds can actually bear E/Z, measured instead of guessed from a bond length.

The key currently calls a bond "locked" when it is shorter than 0.93*(r_a+r_b). The P1 run
showed the enol C-OH sitting 0.0002 A from that line. But a C-OH torsion is a
CONFORMATIONAL coordinate -- the run's own microstates include the rotated-out OH rotamer
0.36 eV up -- so recording its orientation as configurational parity violates the rule the
key already states for itself. Stiffness says which is which; a bond length cannot.
"""
import sys; sys.path.insert(0, "/home/ruigengji/MLPATH/src")
import numpy as np
from ase.io import read
from ase.data import covalent_radii
from prrs.calculators import mace_factory
from prrs.chemistry import canonical_labels
from prrs.config import SearchConfig
from prrs.perturbations import rotatable_torsions
from prrs.reliability import GuardedCalculator
from prrs.runner import torsion_response
from prrs.state import encode

factory = lambda: mace_factory("/home/kasuga/.cache/mace/MACE-OFF24_medium.model",
                               device="cuda", default_dtype="float64")
cfg = SearchConfig(quench_fmax_eV_A=0.002)
names = {0:"O1",1:"C1",2:"C2",3:"C3",4:"O2",5:"Hb",6:"H1",7:"H2",8:"H3"}
for path in ("p1_run/chemical/c0000_m0000.extxyz", "p1_run/chemical/c0001_m0000.extxyz"):
    atoms = read(path)
    guard = GuardedCalculator(factory(), cfg, None)
    atoms.calc = guard
    graph = encode(atoms, cfg.bond_scale)
    labels = canonical_labels(atoms.numbers, graph.edges, graph.index)
    genuine, rotors = rotatable_torsions(graph, labels)
    print(f"\n=== {path.split('/')[-1]} ===")
    print(f"{'torsion':<22} {'bond':<10} {'length':>8} {'0.93 limit':>11} {'k (eV/rad^2)':>14} {'tier':>8}")
    for t in genuine + rotors:
        i, j, k, l = t.indices
        d = np.linalg.norm(atoms.positions[j] - atoms.positions[k])
        lim = 0.93 * (covalent_radii[atoms.numbers[j]] + covalent_radii[atoms.numbers[k]])
        _, curvature = torsion_response(atoms, guard, t, cfg.soft_mode_step_rad)
        tier = ("stiff" if curvature > cfg.soft_mode_curvature_max_eV_rad2
                else "free" if abs(curvature) < cfg.soft_mode_curvature_floor_eV_rad2
                else "negative" if curvature < 0 else "soft")
        bond = f"{names[j]}-{names[k]}"
        mark = "  <- below limit" if d < lim else ""
        print(f"{str(t.indices):<22} {bond:<10} {d:>8.4f} {lim:>11.4f} {curvature:>14.3f} {tier:>8}{mark}")
