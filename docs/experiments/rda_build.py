import sys; sys.path.insert(0, "src")
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from ase import Atoms
from ase.optimize import FIRE
from prrs.calculators import mace_factory
from prrs.state import encode

mol = Chem.AddHs(Chem.MolFromSmiles("C1=CCCCC1"))
AllChem.EmbedMolecule(mol, randomSeed=0xC0FFEE)
AllChem.MMFFOptimizeMolecule(mol)
pos = mol.GetConformer().GetPositions()
sym = [a.GetSymbol() for a in mol.GetAtoms()]
atoms = Atoms(sym, positions=pos)
calc = mace_factory("models/MACE-OFF24_medium.model", device="cuda", default_dtype="float64")
atoms.calc = calc
FIRE(atoms, logfile=None).run(fmax=1e-3, steps=2000)
print("E =", atoms.get_potential_energy(), " fmax =", np.abs(atoms.get_forces()).max())
print("formula", atoms.get_chemical_formula())

# the ring, and the C=C inside it
graph = encode(atoms)
edges = {tuple(sorted(e)) for e in graph.edges}
carbons = [i for i, s in enumerate(sym) if s == "C"]
cc = [e for e in edges if e[0] in carbons and e[1] in carbons]
lengths = {e: float(np.linalg.norm(atoms.positions[e[0]] - atoms.positions[e[1]])) for e in cc}
for e, d in sorted(lengths.items(), key=lambda kv: kv[1]):
    print(f"  C-C {e} {d:.4f}")
double = min(lengths, key=lengths.get)
# ring order starting from the double bond
adj = {c: [o for e in cc for o in e if c in e and o != c] for c in carbons}
ring = [double[0], double[1]]
while len(ring) < 6:
    ring.append(next(n for n in adj[ring[-1]] if n != ring[-2]))
print("ring order (C1=C2,C3..C6):", ring)
# retro-DA breaks C3-C4 and C5-C6 in that ordering
b1 = tuple(sorted((ring[2], ring[3])))
b2 = tuple(sorted((ring[4], ring[5])))
print("breaking bonds:", b1, lengths[b1], "|", b2, lengths[b2])
np.save("/tmp/claude-1000/-home-ruigengji-MLPATH/f74de3f6-4a7e-4dbb-aed5-a3d30f6fe292/scratchpad/rda_pos.npy", atoms.positions)
open("/tmp/claude-1000/-home-ruigengji-MLPATH/f74de3f6-4a7e-4dbb-aed5-a3d30f6fe292/scratchpad/rda_meta.txt","w").write(
    f"{''.join(sym)}\n{b1}\n{b2}\n{ring}\n{atoms.get_potential_energy()}\n")
