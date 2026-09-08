import sys; sys.path.insert(0, "/home/ruigengji/MLPATH/src")
import json
import numpy as np
from ase import Atoms
from ase.io import write
from ase.optimize import LBFGS
from prrs.calculators import mace_factory
from prrs.chemistry import chemical_key
from prrs.state import encode

O1, O2 = np.array([-1.10, 1.15, 0.0]), np.array([1.24, 1.28, 0.0])
SK = {"O1": O1, "C1": np.array([-1.19,-0.05,0.0]), "C2": np.array([0.00,-0.80,0.0]),
      "C3": np.array([1.19,-0.05,0.0]), "O2": O2,
      "H1": np.array([-2.15,-0.55,0.0]), "H2": np.array([0.00,-1.88,0.0]),
      "H3": np.array([2.15,-0.55,0.0])}
pos = dict(SK); pos["Hb"] = O2 + 0.98*(O1-O2)/np.linalg.norm(O1-O2)
order = ["O1","C1","C2","C3","O2","Hb","H1","H2","H3"]
atoms = Atoms("OCCCOHHHH", positions=[pos[k] for k in order])
atoms.calc = mace_factory("/home/kasuga/.cache/mace/MACE-OFF24_medium.model",
                          device="cuda", default_dtype="float64")
LBFGS(atoms, logfile=None).run(fmax=2e-3, steps=800)
print(f"relaxed fmax={np.abs(atoms.get_forces()).max():.2e}  E={atoms.get_potential_energy():.6f}")
k = chemical_key(atoms)
print("key:", k["key"][:16], "frags:", k["fragments"])
edges = sorted(tuple(sorted(e)) for e in encode(atoms).edges)
print("bonds:", edges)
print("O1..Hb =", np.linalg.norm(atoms.positions[0]-atoms.positions[5]),
      " O2-Hb =", np.linalg.norm(atoms.positions[4]-atoms.positions[5]),
      " O1..O2 =", np.linalg.norm(atoms.positions[0]-atoms.positions[4]))
out = atoms.copy(); out.calc = None
write("p1_input.extxyz", out)

base = json.load(open("p0_config.json"))
base.update(seed=17, max_trials=60,
            families=["stretch", "compress", "bend", "kick", "torsion", "torsion_kick"],
            geometry_amplitudes_A=[0.15, 0.35, 0.5],
            chemical_trials_per_node=120)
json.dump(base, open("p1_config.json", "w"), indent=2)
print("config families:", base["families"], "amplitudes:", base["geometry_amplitudes_A"])
