"""The same C-C bond-breaking smoothness probe, on MACE-POLAR-1-M.

`rda_bondbreak_smoothness.md` answered the go/no-go on MACE-OFF24. The user has since
frozen P0-P3 on OFF24 and moved everything from P4 onward to POLAR-1, so that answer no
longer covers the model the retro-DA would actually run on. Smoothness in the bond-breaking
region is a property of the potential, not of the reaction, so it does not transfer.

POLAR-1 has a specific reason to be checked rather than assumed: P3 section 3.13 measured its
rigid-body floor residual at 2.2e-06 against OFF24's 3.2e-16 -- ten orders of magnitude
worse, consistent with its float32 weights. A floor that size is still far below the
eigenvalue tolerance, but this probe is about curvature in the one region where nothing
has ever been measured, so the model with the worse numerical floor is the one that has to
be looked at, not the one that can be waved through.

Every endpoint is relaxed by POLAR-1 itself. Handing it OFF24's geometry would measure
geometry transfer rather than the surface, which is the rule P3 section 5 already fixed.
"""
import sys, json; sys.path.insert(0, "src")
import numpy as np
from ase import Atoms
from ase.constraints import FixInternals
from ase.optimize import FIRE
from mace.calculators import MACECalculator
from prrs.state import encode

MODEL = "/home/ruigengji/MLP/mace/MACE-POLAR-1-M.model"
OUT = "docs/experiments/rda_polar1.json"
calc = MACECalculator(model_paths=MODEL, model_type="PolarMACE", device="cuda",
                      default_dtype="float64")

from rdkit import Chem
from rdkit.Chem import AllChem
mol = Chem.AddHs(Chem.MolFromSmiles("C1=CCCCC1"))
AllChem.EmbedMolecule(mol, randomSeed=0xC0FFEE)
AllChem.MMFFOptimizeMolecule(mol)
atoms = Atoms([a.GetSymbol() for a in mol.GetAtoms()],
              positions=mol.GetConformer().GetPositions())
atoms.calc = calc
FIRE(atoms, logfile=None).run(fmax=1e-3, steps=2000)
E_ref = float(atoms.get_potential_energy())
print(f"cyclohexene on POLAR-1: E = {E_ref:.6f} eV  fmax = "
      f"{np.abs(atoms.get_forces()).max():.2e}", flush=True)

graph = encode(atoms)
edges = {tuple(sorted(e)) for e in graph.edges}
carbons = [i for i, s in enumerate(atoms.symbols) if s == "C"]
cc = [e for e in edges if e[0] in carbons and e[1] in carbons]
length = {e: float(np.linalg.norm(atoms.positions[e[0]] - atoms.positions[e[1]])) for e in cc}
double = min(length, key=length.get)
adj = {c: [o for e in cc for o in e if c in e and o != c] for c in carbons}
ring = [double[0], double[1]]
while len(ring) < 6:
    ring.append(next(n for n in adj[ring[-1]] if n != ring[-2]))
b1, b2 = tuple(sorted((ring[2], ring[3]))), tuple(sorted((ring[4], ring[5])))
print(f"ring {ring}  C=C {double} {length[double]:.4f}  breaking {b1} {length[b1]:.4f} "
      f"/ {b2} {length[b2]:.4f}", flush=True)

def set_length(a, pair, target):
    i, j = pair
    v = a.positions[j] - a.positions[i]
    s = 0.5 * (target - np.linalg.norm(v)) * v / np.linalg.norm(v)
    a.positions[i] -= s; a.positions[j] += s

geoms, rows = {}, []
for d in np.round(np.arange(1.55, 2.401, 0.01), 4):
    set_length(atoms, b1, float(d)); set_length(atoms, b2, float(d))
    atoms.set_constraint(FixInternals(bonds=[[float(d), list(b1)], [float(d), list(b2)]]))
    ok = FIRE(atoms, logfile=None).run(fmax=2e-3, steps=2000)
    atoms.set_constraint()
    E = float(atoms.get_potential_energy())
    geoms[float(d)] = atoms.positions.copy()
    rows.append(dict(d=float(d), rel_meV=(E - E_ref) * 1e3, conv=bool(ok),
                     edges=len(encode(atoms).edges)))
    if abs(round(d, 2) * 100 % 5) < 1e-6:
        print(f"d={d:.2f}  E-E0={rows[-1]['rel_meV']:9.1f} meV  conv={int(ok)} "
              f"edges={rows[-1]['edges']}", flush=True)

e = np.array([r["rel_meV"] for r in rows])
steps = np.diff(e)
jump = int(np.argmax(np.abs(steps)))
print(f"\nlargest single step: d={rows[jump]['d']:.2f} -> {rows[jump+1]['d']:.2f}  "
      f"dE={steps[jump]:+.1f} meV", flush=True)
print(f"|second difference| excluding that step: max "
      f"{np.abs(np.delete(np.diff(e, 2), [max(jump-1, 0), jump])).max():.1f} meV", flush=True)

# The discriminator: a rigid line between the two relaxed geometries that bracket the
# largest step. No relaxation, so a kink here is the potential and not the optimizer.
probe = atoms.copy(); probe.calc = calc
A, B = geoms[rows[jump]["d"]], geoms[rows[jump + 1]["d"]]
line = []
for t in np.linspace(0.0, 1.0, 101):
    probe.set_positions((1 - t) * A + t * B)
    line.append(float(probe.get_potential_energy()))
line = np.array(line)
d2 = np.abs(np.diff(line, 2)) * 1e3
print(f"rigid line: E {(line.min()-E_ref)*1e3:.1f} .. {(line.max()-E_ref)*1e3:.1f} meV  "
      f"|second difference| max {d2.max():.4f} meV mean {d2.mean():.4f} meV  "
      f"worst at t={np.argmax(d2)/100:.2f}", flush=True)
json.dump(dict(model="MACE-POLAR-1-M", E_ref=E_ref, b1=list(b1), b2=list(b2),
               rows=rows, rigid_line_eV=[float(x) for x in line]), open(OUT, "w"), indent=1)
print("wrote", OUT, flush=True)
