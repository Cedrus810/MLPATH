"""Same-ordering endpoint construction for benzoylacetone. Capability probe, not preflight.

The first probe built A and B from two separate SMILES, so rdkit gave them different atom
orderings and classify_transition read the index mismatch as a double proton shuffle
(broken [[3,15],[5,16]], formed [[2,15],[3,16]]). That is a bookkeeping artefact. P2 avoided
it by building both endpoints from one skeleton and only moving the proton. This checks that
the same construction works here, and reports what the event key actually is.
"""
import sys, json
sys.path.insert(0, "/home/ruigengji/MLPATH/src")
import numpy as np
from ase import Atoms
from ase.optimize import LBFGS
from rdkit import Chem
from rdkit.Chem import AllChem
from prrs.calculators import mace_factory
from prrs.state import encode
from prrs.chemistry import chemical_key, classify_transition, reaction_event_key
from prrs.config import SearchConfig
from prrs.runner import curvature_spectrum, floor_residual

MODEL = "/home/kasuga/.cache/mace/MACE-OFF24_medium.model"
factory = lambda: mace_factory(MODEL, device="cuda", default_dtype="float64")
cfg = SearchConfig(quench_fmax_eV_A=0.002)

# A: chelated phenyl-side enol. Atom order is fixed by this one embedding and reused.
m = Chem.AddHs(Chem.MolFromSmiles(r"C\C(=O)\C=C(\O)c1ccccc1"))
ps = AllChem.ETKDGv3(); ps.randomSeed = 0xC0FFEE
assert AllChem.EmbedMolecule(m, ps) == 0
AllChem.MMFFOptimizeMolecule(m, maxIters=2000)
conf = m.GetConformer()
A = Atoms([a.GetSymbol() for a in m.GetAtoms()],
          positions=[list(conf.GetAtomPosition(i)) for i in range(m.GetNumAtoms())])

sym = list(A.symbols)
print("atom order:", " ".join(f"{i}{s}" for i, s in enumerate(sym)))

A.calc = factory()
LBFGS(A, logfile=None).run(fmax=2e-3, steps=2000)
gA = encode(A)
edgesA = sorted(tuple(sorted(e)) for e in gA.edges)

# Locate the enol proton and the two oxygens from A's own relaxed graph, not by assumption.
oxy = [i for i, s in enumerate(sym) if s == "O"]
hyd = [i for i, s in enumerate(sym) if s == "H"]
oh = [(o, h) for o in oxy for h in hyd if tuple(sorted((o, h))) in edgesA]
assert len(oh) == 1, oh
o_donor, h_shared = oh[0]
o_acceptor = [o for o in oxy if o != o_donor][0]
d_OO = np.linalg.norm(A.positions[o_donor] - A.positions[o_acceptor])
d_acc = np.linalg.norm(A.positions[o_acceptor] - A.positions[h_shared])
print(f"\nA relaxed: E={A.get_potential_energy():.6f} eV  fmax={np.abs(A.get_forces()).max():.2e}")
print(f"  donor O{o_donor} - H{h_shared}, acceptor O{o_acceptor}")
print(f"  r(O-O) = {d_OO:.3f} A   r(O_acceptor...H) = {d_acc:.3f} A   (chelate closed if < ~1.8)")

# B: same ordering, proton relocated onto the acceptor oxygen along the O...O line.
B = A.copy()
u = A.positions[o_donor] - A.positions[o_acceptor]
B.positions[h_shared] = A.positions[o_acceptor] + 0.98 * u / np.linalg.norm(u)
B.calc = factory()
LBFGS(B, logfile=None).run(fmax=2e-3, steps=2000)
print(f"\nB relaxed: E={B.get_potential_energy():.6f} eV  fmax={np.abs(B.get_forces()).max():.2e}")

report = {}
for tag, atoms in (("A", A), ("B", B)):
    key = chemical_key(atoms)
    lam, vec, floor, prov = curvature_spectrum(atoms, factory, cfg, source="analytic")
    neg = [float(x) for x in lam if x < -cfg.minimum_check_eigenvalue_tol]
    report[tag] = dict(E_eV=float(atoms.get_potential_energy()),
                       fmax=float(np.abs(atoms.get_forces()).max()),
                       key=key["key"][:16], graph_hash=key["graph_hash"],
                       fragments=key["fragments"],
                       stereo_unresolved=key["stereo_unresolved"],
                       n_locked_parity=len(key["locked_bond_parity"]),
                       n_minus=len(neg), lam6=[round(float(x), 6) for x in lam[:6]],
                       floor_residual=float(floor_residual(floor)),
                       n_edges=len(encode(atoms).edges))
    print(f"\n[{tag}] " + "  ".join(f"{k}={v}" for k, v in report[tag].items()))

cls = classify_transition(A, B, cfg.bond_scale, cfg.active_atoms)
ev_f = reaction_event_key(A, B, cfg.bond_scale, cfg.active_atoms)
ev_r = reaction_event_key(B, A, cfg.bond_scale, cfg.active_atoms)
print(f"\nclassify_transition(A,B) = {cls}")
print(f"event_key forward = {ev_f['key']}  broken={ev_f['broken']} formed={ev_f['formed']}")
print(f"event_key reverse = {ev_r['key']}  broken={ev_r['broken']} formed={ev_r['formed']}")
print(f"direction quotiented (keys equal) = {ev_f['key'] == ev_r['key']}")
print(f"keys differ = {report['A']['key'] != report['B']['key']}   "
      f"graph_hash differ = {report['A']['graph_hash'] != report['B']['graph_hash']}")
print(f"dE(B-A) = {(report['B']['E_eV'] - report['A']['E_eV'])*1e3:.3f} meV")

report["meta"] = dict(o_donor=o_donor, o_acceptor=o_acceptor, h_shared=h_shared,
                      r_OO_A=float(d_OO), event_forward=ev_f["key"],
                      event_reverse=ev_r["key"], classify=str(cls),
                      broken=ev_f["broken"], formed=ev_f["formed"],
                      atom_order="".join(sym))
json.dump(report, open("/home/ruigengji/MLPATH/runs/p3_ordering_check.json", "w"),
          indent=2, default=str)
from ase.io import write
write("/home/ruigengji/MLPATH/runs/p3_A_common_order.extxyz", A)
write("/home/ruigengji/MLPATH/runs/p3_B_common_order.extxyz", B)
print("\nwrote runs/p3_ordering_check.json and the two extxyz endpoints")
