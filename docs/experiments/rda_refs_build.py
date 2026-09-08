"""Freeze four geometries in the C-C bond-breaking region and read both models' curvature.

The smoothness probes (`rda_bondbreak_smoothness.md`) answer self-consistency: is the
surface smooth. They cannot answer accuracy -- `PRRS_STATUS` section 2.6c's standing
limitation, "a reference-free diagnostic can find inconsistency, never inaccuracy" --
and the bond-breaking region is exactly where nothing has ever been checked against
anything. There is now a calibrated ORCA protocol on the 40-core box, so a ruler exists
for the first time.

This produces what that comparison needs and nothing more:

  four geometries spanning the region, taken from the POLAR-1 relaxed walk because P4
  onward runs on POLAR-1. d = 1.55 is the control near equilibrium; 1.90, 2.05, 2.20 are
  in the breaking region, and 2.05 sits in the band where the constrained scan switched
  branches.

  BOTH models' analytic Hessian spectra AT THOSE SAME GEOMETRIES. Same geometry for both
  is the whole point: relaxing each model separately would compare two different points
  and measure geometry transfer instead of curvature.

The geometries are not stationary points, so by section 2.6g the three rotational
eigenvalues carry a gradient term and are not model error. They are reported, not
compared.
"""
import sys, json; sys.path.insert(0, "src")
import numpy as np
from ase import Atoms
from ase.io import write
from ase.constraints import FixInternals
from ase.optimize import FIRE
from mace.calculators import MACECalculator
from prrs.config import SearchConfig
from prrs.runner import curvature_spectrum, floor_residual
from prrs.state import encode
from rdkit import Chem
from rdkit.Chem import AllChem

OUT = "docs/experiments/rda_refs/"
POLAR = "/home/ruigengji/MLP/mace/MACE-POLAR-1-M.model"
OFF24 = "models/MACE-OFF24_medium.model"
WANT = [1.55, 1.90, 2.05, 2.20]

polar = lambda: MACECalculator(model_paths=POLAR, model_type="PolarMACE",
                               device="cuda", default_dtype="float64")
off24 = lambda: MACECalculator(model_paths=OFF24, model_type="MACE",
                               device="cuda", default_dtype="float64")
cfg = SearchConfig(closed_shell_only=False, quench_fmax_eV_A=0.002)

mol = Chem.AddHs(Chem.MolFromSmiles("C1=CCCCC1"))
AllChem.EmbedMolecule(mol, randomSeed=0xC0FFEE)
AllChem.MMFFOptimizeMolecule(mol)
atoms = Atoms([a.GetSymbol() for a in mol.GetAtoms()],
              positions=mol.GetConformer().GetPositions())
calc = polar()
atoms.calc = calc
FIRE(atoms, logfile=None).run(fmax=1e-3, steps=2000)

graph = encode(atoms)
cc = [e for e in ({tuple(sorted(x)) for x in graph.edges})
      if atoms.symbols[e[0]] == "C" and atoms.symbols[e[1]] == "C"]
length = {e: float(np.linalg.norm(atoms.positions[e[0]] - atoms.positions[e[1]])) for e in cc}
double = min(length, key=length.get)
adj = {c: [o for e in cc for o in e if c in e and o != c]
       for c in range(len(atoms)) if atoms.symbols[c] == "C"}
ring = [double[0], double[1]]
while len(ring) < 6:
    ring.append(next(n for n in adj[ring[-1]] if n != ring[-2]))
b1, b2 = tuple(sorted((ring[2], ring[3]))), tuple(sorted((ring[4], ring[5])))
print(f"breaking bonds {b1} {b2}   ring {ring}", flush=True)

def set_length(a, pair, target):
    i, j = pair
    v = a.positions[j] - a.positions[i]
    s = 0.5 * (target - np.linalg.norm(v)) * v / np.linalg.norm(v)
    a.positions[i] -= s; a.positions[j] += s

frozen = {}
for d in np.round(np.arange(1.55, 2.201, 0.01), 4):
    set_length(atoms, b1, float(d)); set_length(atoms, b2, float(d))
    atoms.set_constraint(FixInternals(bonds=[[float(d), list(b1)], [float(d), list(b2)]]))
    FIRE(atoms, logfile=None).run(fmax=2e-3, steps=2000)
    atoms.set_constraint()
    if any(abs(d - w) < 1e-9 for w in WANT):
        frozen[round(float(d), 2)] = atoms.positions.copy()
        print(f"  frozen d={d:.2f}", flush=True)

records = []
for d, positions in sorted(frozen.items()):
    probe = Atoms("C6H10", positions=positions)
    row = {"d": d, "n_edges": len(encode(probe).edges)}
    write(f"{OUT}rda_d{d:.2f}.xyz", probe)
    for name, factory in (("POLAR-1-M", polar), ("OFF24", off24)):
        probe.calc = factory()
        lam, vec, floor, prov = curvature_spectrum(probe, factory, cfg, source="analytic")
        internal = [float(x) for x in lam if abs(x) > 1e-4][:8]
        row[name] = {
            "E_eV": float(probe.get_potential_energy()),
            "fmax_eV_A": float(np.abs(probe.get_forces()).max()),
            "lowest_eigenvalues": [float(x) for x in lam[:9]],
            "n_negative": int(sum(1 for x in lam if x < -cfg.minimum_check_eigenvalue_tol)),
            "floor_residual": float(floor_residual(floor)),
            "source": prov["source"]}
        print(f"d={d:.2f} {name:10s} E={row[name]['E_eV']:14.6f} "
              f"fmax={row[name]['fmax_eV_A']:.2e} n_neg={row[name]['n_negative']} "
              f"floor={row[name]['floor_residual']:.2e} "
              f"lam[:4]={[round(x,4) for x in row[name]['lowest_eigenvalues'][:4]]}", flush=True)
    records.append(row)

json.dump({"breaking_bonds": [list(b1), list(b2)], "ring": ring, "records": records},
          open(f"{OUT}rda_refs.json", "w"), indent=1)
print("wrote", OUT, flush=True)
