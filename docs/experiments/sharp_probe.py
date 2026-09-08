import sys; sys.path.insert(0,'/home/ruigengji/MLPATH/src')
import numpy as np
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from prrs.calculators import DoubleWell
from prrs.config import SearchConfig
from prrs.runner import curvature_spectrum

class SharpWell(DoubleWell):
    """Double well with a narrow bump on the ridge: same minima, a far sharper saddle."""
    HEIGHT, WIDTH = 0.6, 0.12
    def calculate(self, atoms=None, properties=("energy","forces"), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        v = atoms.positions[1]-atoms.positions[0]; r = float(np.linalg.norm(v))
        u = (r-1.8)/self.WIDTH
        bump = self.HEIGHT*np.exp(-u*u)
        d = -2*u*bump/self.WIDTH
        self.results["energy"] += float(bump)
        f = d*v/r
        self.results["forces"] = self.results["forces"] + np.array([f, -f])

cfg = SearchConfig(quench_fmax_eV_A=0.001)
for r in (1.2, 1.8, 2.4):
    a = Atoms("C2", positions=[[-r/2,0,0],[r/2,0,0]]); a.calc = SharpWell()
    lam,_,_,_ = curvature_spectrum(a, SharpWell, cfg, source="fd")
    print(f"r={r}  E={a.get_potential_energy():+.4f}  fmax={np.abs(a.get_forces()).max():.4f}  "
          f"lowest={lam[0]:+.4f}  highest={lam[-1]:+.4f}")
