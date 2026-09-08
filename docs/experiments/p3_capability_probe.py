"""P3 capability probe. NOT the preflight gate, NOT a search run.

Purpose: measure whether the identity layer and the cost model can carry benzoylacetone
at all, BEFORE the P3 acceptance criteria are written. The project rule is that criteria
are frozen before the run; but a scale criterion cannot be written without knowing what
22 atoms cost, and an aromatic-ring criterion cannot be written without knowing whether
admissible_bond_orders resolves a Kekule pair. Both are measured here at zero search
budget. Everything printed is a capability number, not an acceptance result.
"""
import sys, time, json
sys.path.insert(0, "/home/ruigengji/MLPATH/src")
import numpy as np
from ase import Atoms
from ase.optimize import LBFGS
from rdkit import Chem
from rdkit.Chem import AllChem
from prrs.calculators import mace_factory
from prrs.state import encode, resolve_active
from prrs.chemistry import (chemical_key, classify_transition, reaction_event_key,
                            admissible_bond_orders, locked_edges)
from prrs.config import SearchConfig
from prrs.runner import curvature_spectrum, floor_residual

MODEL = "/home/kasuga/.cache/mace/MACE-OFF24_medium.model"
factory = lambda: mace_factory(MODEL, device="cuda", default_dtype="float64")
cfg = SearchConfig(quench_fmax_eV_A=0.002)

SMILES = {
    # chelated enols of benzoylacetone (1-phenyl-1,3-butanedione), C10H10O2
    "A_phenyl_enol": r"C\C(=O)\C=C(\O)c1ccccc1",   # OH on the phenyl-side carbon
    "B_methyl_enol": r"C\C(\O)=C\C(=O)c1ccccc1",   # OH on the methyl-side carbon
    "K_diketone":    r"CC(=O)CC(=O)c1ccccc1",      # the beta-diketone tautomer
}

def build(smi, seed=0xC0FFEE):
    m = Chem.AddHs(Chem.MolFromSmiles(smi))
    ps = AllChem.ETKDGv3(); ps.randomSeed = seed
    assert AllChem.EmbedMolecule(m, ps) == 0
    AllChem.MMFFOptimizeMolecule(m, maxIters=2000)
    c = m.GetConformer()
    return Atoms([a.GetSymbol() for a in m.GetAtoms()],
                 positions=[list(c.GetAtomPosition(i)) for i in range(m.GetNumAtoms())])

out = {}
relaxed = {}
for tag, smi in SMILES.items():
    atoms = build(smi); atoms.calc = factory()
    n = len(atoms)
    t0 = time.perf_counter(); atoms.get_potential_energy(); t_first = time.perf_counter() - t0
    t0 = time.perf_counter()
    for _ in range(20):
        atoms.calc.results.clear(); atoms.get_potential_energy(); atoms.get_forces()
    t_pt = (time.perf_counter() - t0) / 20

    t0 = time.perf_counter()
    opt = LBFGS(atoms, logfile=None); opt.run(fmax=2e-3, steps=2000)
    t_relax = time.perf_counter() - t0
    steps = opt.get_number_of_steps()

    key = chemical_key(atoms)
    graph = encode(atoms)
    edges = sorted(tuple(sorted(e)) for e in graph.edges)
    index = resolve_active(atoms, cfg.active_atoms)
    sols, info = admissible_bond_orders(atoms.numbers, graph.edges, index)
    n_sols = None if sols is None else len(sols)

    t0 = time.perf_counter()
    lam, vec, floor, prov = curvature_spectrum(atoms, factory, cfg, source="analytic")
    t_hess = time.perf_counter() - t0
    neg = [float(x) for x in lam if x < -cfg.minimum_check_eigenvalue_tol]

    relaxed[tag] = atoms
    out[tag] = dict(n_atoms=n, formula=str(atoms.symbols.get_chemical_formula()),
                    E_eV=float(atoms.get_potential_energy()),
                    fmax=float(np.abs(atoms.get_forces()).max()),
                    t_first_call_s=t_first, t_per_point_ms=t_pt*1e3,
                    t_relax_s=t_relax, relax_steps=int(steps),
                    t_hessian_s=t_hess, n_edges=len(edges),
                    n_admissible_bond_orders=n_sols, bond_order_info=info,
                    n_locked_parity=len(key["locked_bond_parity"]),
                    stereo_unresolved=key["stereo_unresolved"],
                    key=key["key"][:16], graph_hash=key["graph_hash"],
                    fragments=key["fragments"],
                    n_minus=len(neg), lambda_lowest=[round(float(x), 6) for x in lam[:6]],
                    floor_residual=float(floor_residual(floor)))
    print(f"\n=== {tag} ===")
    for k, v in out[tag].items():
        print(f"  {k:26s} {v}")

print("\n=== identity separation ===")
for x, y in (("A_phenyl_enol", "B_methyl_enol"), ("A_phenyl_enol", "K_diketone"),
             ("B_methyl_enol", "K_diketone")):
    kx, ky = chemical_key(relaxed[x]), chemical_key(relaxed[y])
    cls = classify_transition(relaxed[x], relaxed[y], cfg.bond_scale, cfg.active_atoms)
    ev = reaction_event_key(relaxed[x], relaxed[y], cfg.bond_scale, cfg.active_atoms)
    print(f"  {x} vs {y}")
    print(f"    key differs        {kx['key'] != ky['key']}")
    print(f"    graph_hash differs {kx['graph_hash'] != ky['graph_hash']}")
    print(f"    fragments          {kx['fragments']} / {ky['fragments']}")
    print(f"    classify           {cls if not isinstance(cls, dict) else {k: cls[k] for k in list(cls)[:6]}}")
    print(f"    event_key          {ev if not isinstance(ev, dict) else {k: ev[k] for k in list(ev)[:6]}}")
    dE = (relaxed[y].get_potential_energy() - relaxed[x].get_potential_energy())*1e3
    print(f"    dE(y-x)            {dE:.3f} meV")

json.dump(out, open("/home/ruigengji/MLPATH/runs/p3_capability_probe.json", "w"), indent=2, default=str)
print("\nwrote runs/p3_capability_probe.json")
