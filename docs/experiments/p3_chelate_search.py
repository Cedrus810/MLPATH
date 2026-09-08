"""Chelated (Z) enol of benzoylacetone: endpoints, minima confirmation, conformer landscape.

Capability probe, zero search budget. NOT the preflight gate.

Two earlier probes were wrong about the substance and are kept in the record:
  p3_capability_probe   built A and B from separate SMILES, so different rdkit atom
                        orderings made classify_transition read the index mismatch as a
                        double proton shuffle.
  p3_ordering_check     fixed the ordering but used C\\C(=O)\\C=C(\\O)c1ccccc1, which rdkit
                        assigns STEREOE -- the OH points away from the carbonyl and
                        r(O...O) = 4.29 A over every one of its conformers. That enol
                        cannot chelate at all, so it is not the P3 reaction.

The chelating enol is the Z isomer, CC(=O)/C=C(\\O)c1ccccc1: 151 of 300 ETKDG conformers
come out with r(O...O) < 3.0 A, minimum 2.615 A.
"""
import sys, json
sys.path.insert(0, "/home/ruigengji/MLPATH/src")
import numpy as np
from ase import Atoms
from ase.io import write
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
SMI = r"C\C(=O)\C=C(/O)c1ccccc1"          # STEREOZ, the chelating enol
calc = factory()

mol = Chem.MolFromSmiles(SMI)
assert [str(b.GetStereo()) for b in mol.GetBonds()
        if b.GetStereo() != Chem.BondStereo.STEREONONE] == ["STEREOZ"]
m = Chem.AddHs(mol)
ps = AllChem.ETKDGv3(); ps.randomSeed = 0xC0FFEE
ids = AllChem.EmbedMultipleConfs(m, numConfs=300, params=ps)
AllChem.MMFFOptimizeMoleculeConfs(m, maxIters=4000)
sym = [a.GetSymbol() for a in m.GetAtoms()]
oxy = [i for i, s in enumerate(sym) if s == "O"]
hyd = [i for i, s in enumerate(sym) if s == "H"]
print(f"canonical SMILES {Chem.MolToSmiles(mol)}   atom order {''.join(sym)}   oxygens {oxy}")

rows = []
for cid in ids:
    c = m.GetConformer(cid)
    P = np.array([list(c.GetAtomPosition(i)) for i in range(m.GetNumAtoms())])
    rows.append((float(np.linalg.norm(P[oxy[0]] - P[oxy[1]])), int(cid), P))
rows.sort(key=lambda r: r[0])
closed = [r for r in rows if r[0] < 3.0]
print(f"MMFF conformers {len(rows)}  closed(r_OO<3.0) {len(closed)}  "
      f"r_OO min {rows[0][0]:.3f} max {rows[-1][0]:.3f}")

def measure(atoms, tag, spectrum=True):
    g = encode(atoms)
    edges = sorted(tuple(sorted(e)) for e in g.edges)
    oh = [(o, h) for o in oxy for h in hyd if tuple(sorted((o, h))) in edges]
    key = chemical_key(atoms)
    d = dict(tag=tag, E_eV=float(atoms.get_potential_energy()),
             fmax=float(np.abs(atoms.get_forces()).max()),
             r_OO=float(np.linalg.norm(atoms.positions[oxy[0]] - atoms.positions[oxy[1]])),
             n_OH=len(oh), n_edges=len(edges),
             key=key["key"][:16], graph_hash=key["graph_hash"],
             fragments=key["fragments"], stereo_unresolved=key["stereo_unresolved"],
             n_locked_parity=len(key["locked_bond_parity"]))
    if len(oh) == 1:
        donor, h_sh = oh[0]
        acc = [o for o in oxy if o != donor][0]
        d.update(donor=donor, acceptor=acc, h_shared=h_sh,
                 r_acceptorH=float(np.linalg.norm(atoms.positions[acc]
                                                  - atoms.positions[h_sh])))
    if spectrum:
        lam, _, floor, prov = curvature_spectrum(atoms, factory, cfg, source="analytic")
        neg = [float(x) for x in lam if x < -cfg.minimum_check_eigenvalue_tol]
        d.update(n_minus=len(neg), lam6=[round(float(x), 6) for x in lam[:6]],
                 floor_residual=float(floor_residual(floor)),
                 curvature_source=prov["source"])
    return d

# ---- A: tightest chelate, relaxed under MACE -------------------------------------
A = Atoms(sym, positions=closed[0][2]); A.calc = calc
LBFGS(A, logfile=None).run(fmax=2e-3, steps=2000)
mA = measure(A, "A_chelated_enol")
print("\n=== A (chelated phenyl-side enol) ===")
for k, v in mA.items(): print(f"  {k:18s} {v}")
print(f"  chelate survived MACE relaxation: {mA['r_OO'] < 3.0}  "
      f"(MMFF {closed[0][0]:.3f} -> MACE {mA['r_OO']:.3f} A)")

# ---- B: same atom ordering, proton moved across the chelate ----------------------
donor, acc, h_sh = mA["donor"], mA["acceptor"], mA["h_shared"]
B = A.copy()
u = A.positions[donor] - A.positions[acc]
B.positions[h_sh] = A.positions[acc] + 0.99 * u / np.linalg.norm(u)
B.calc = factory()
LBFGS(B, logfile=None).run(fmax=2e-3, steps=2000)
mB = measure(B, "B_chelated_enol")
print("\n=== B (methyl-side enol, proton moved across the chelate) ===")
for k, v in mB.items(): print(f"  {k:18s} {v}")

cls = classify_transition(A, B, cfg.bond_scale, cfg.active_atoms)
evf = reaction_event_key(A, B, cfg.bond_scale, cfg.active_atoms)
evr = reaction_event_key(B, A, cfg.bond_scale, cfg.active_atoms)
print(f"\nclassify_transition(A,B)   {cls}")
print(f"event forward              {evf['key']}  broken={evf['broken']} formed={evf['formed']}")
print(f"event reverse              {evr['key']}  broken={evr['broken']} formed={evr['formed']}")
print(f"direction quotiented       {evf['key'] == evr['key']}")
print(f"keys differ                {mA['key'] != mB['key']}")
print(f"graph_hash differ          {mA['graph_hash'] != mB['graph_hash']}")
print(f"fragments equal            {mA['fragments'] == mB['fragments']}")
print(f"dE(B-A)                    {(mB['E_eV'] - mA['E_eV'])*1e3:.3f} meV")

# ---- conformer landscape: how many distinct MACE minima share A's chemical_key ----
# This is a P3 scale number, not an acceptance result. Relaxes a spread of conformers and
# counts distinct energies; it bounds how many microstates the reservoir must hold.
print("\n=== conformer landscape (MACE-relaxed, A's key) ===")
sample = [closed[i] for i in np.linspace(0, len(closed) - 1, 16).astype(int)]
sample += [rows[-1], rows[-6]]
land = []
for r0, cid, P in sample:
    a = Atoms(sym, positions=P); a.calc = calc
    LBFGS(a, logfile=None).run(fmax=2e-3, steps=2000)
    mm = measure(a, f"conf{cid}", spectrum=False)
    land.append(mm)
    print(f"  conf{cid:<4d} r_OO {r0:.3f}->{mm['r_OO']:.3f}  E {mm['E_eV']:.6f}  "
          f"key {mm['key']}  nOH {mm['n_OH']}")
Es = np.array([x["E_eV"] for x in land])
uniq_E = len({round(e, 4) for e in Es})           # 0.1 meV resolution
keys = {x["key"] for x in land}
print(f"  {len(land)} relaxed  distinct energies (0.1 meV) {uniq_E}  distinct keys {len(keys)}")
print(f"  E span {(Es.max()-Es.min())*1e3:.1f} meV   lowest {Es.min():.6f} eV")
print(f"  A is the lowest of the sample: {abs(Es.min() - mA['E_eV']) < 1e-6}")

out = dict(smiles=SMI, canonical=Chem.MolToSmiles(mol), atom_order="".join(sym),
           oxygens=oxy, mmff_conformers=len(rows), mmff_closed=len(closed),
           A=mA, B=mB, classify=str(cls),
           event_forward=evf["key"], event_reverse=evr["key"],
           broken=evf["broken"], formed=evf["formed"],
           dE_B_minus_A_meV=(mB["E_eV"] - mA["E_eV"]) * 1e3,
           landscape=land, landscape_distinct_energies=uniq_E,
           landscape_distinct_keys=sorted(keys),
           landscape_span_meV=float((Es.max() - Es.min()) * 1e3))
json.dump(out, open("/home/ruigengji/MLPATH/runs/p3_chelate_search.json", "w"),
          indent=2, default=str)
write("/home/ruigengji/MLPATH/runs/p3_A_chelated.extxyz", A)
write("/home/ruigengji/MLPATH/runs/p3_B_chelated.extxyz", B)
print("\nwrote runs/p3_chelate_search.json and the two chelated endpoints")
