import sys, json; sys.path.insert(0, "src")
import numpy as np
from ase import Atoms
from ase.constraints import FixInternals
from ase.optimize import FIRE
from prrs.calculators import mace_factory
from prrs.state import encode

S = "/tmp/claude-1000/-home-ruigengji-MLPATH/f74de3f6-4a7e-4dbb-aed5-a3d30f6fe292/scratchpad/"
meta = open(S+"rda_meta.txt").read().split("\n")
b1, b2 = eval(meta[1]), eval(meta[2])
atoms = Atoms("C6H10", positions=np.load(S+"rda_pos.npy"))
calc = mace_factory("models/MACE-OFF24_medium.model", device="cuda", default_dtype="float64")
atoms.calc = calc
E_ref = float(atoms.get_potential_energy())

def set_length(a, pair, target):
    """Move both atoms of `pair` symmetrically along their axis to the target length."""
    i, j = pair
    v = a.positions[j] - a.positions[i]
    r = np.linalg.norm(v)
    shift = 0.5 * (target - r) * v / r
    a.positions[i] -= shift
    a.positions[j] += shift

grid = np.round(np.arange(1.55, 3.201, 0.05), 4)
rows = []
for d in grid:
    set_length(atoms, b1, float(d))
    set_length(atoms, b2, float(d))
    atoms.set_constraint(FixInternals(bonds=[[float(d), list(b1)], [float(d), list(b2)]]))
    ok = FIRE(atoms, logfile=None).run(fmax=2e-3, steps=1500)
    atoms.set_constraint()
    E = float(atoms.get_potential_energy())
    r1 = float(np.linalg.norm(atoms.positions[b1[0]]-atoms.positions[b1[1]]))
    r2 = float(np.linalg.norm(atoms.positions[b2[0]]-atoms.positions[b2[1]]))
    g = encode(atoms)
    rows.append(dict(d=float(d), E=E, rel_meV=(E-E_ref)*1e3, converged=bool(ok),
                     fmax=float(np.abs(atoms.get_forces()).max()), r1=r1, r2=r2,
                     n_edges=len(g.edges)))
    print(f"d={d:.3f}  r=({r1:.3f},{r2:.3f})  E-E0={rows[-1]['rel_meV']:9.1f} meV  "
          f"fmax={rows[-1]['fmax']:.2e} conv={int(ok)} edges={rows[-1]['n_edges']}", flush=True)
json.dump(dict(E_ref=E_ref, b1=list(b1), b2=list(b2), rows=rows), open(S+"rda_scan.json","w"), indent=1)
