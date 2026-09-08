"""HEAVY 1 of 3. Full conformer landscape of the chelated benzoylacetone enol.

Why the full sweep: the local probe relaxed 18 of 300 conformers and they collapsed onto
4 distinct energies within 1.5 meV, all one chemical_key. That is a sample, not a count.
The reservoir capacity criterion for P3 (conformer_max_per_node) needs the real number of
distinct closed-chelate minima, and a Hessian on each -- the local probe skipped the
Hessian on the landscape members to save time, so none of those 18 is confirmed.

A Hessian, not `confirm_minimum`: this counts the negatives off `curvature_spectrum`
directly and none of the gates -- the rigid-body floor, the internal-mode recheck -- runs
here. An earlier version of this line said confirm_minimum and it was repeated in
docs/experiments/landscape_regression.py, where it turned into a claim that this job was
the one exercising the recheck gate. It is not; `ts` is, through `descend_saddle`.

Cost on the node (23.8 ms/point): ~150 relax (about 150 calls each) + 300 Hessians.
Estimate 300 x (150 x 0.024 + 4.5/2.2) s ~ 3.3 h. Sequential by design; see README.
"""
import sys, json, time
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "src"))
import numpy as np
from ase import Atoms
from ase.io import write
from ase.optimize import LBFGS
from rdkit import Chem
from rdkit.Chem import AllChem
from _env_guard import guard, shard_arg

SHARD, NSHARD, ARGV = shard_arg(sys.argv)
MODEL = ARGV[1] if len(ARGV) > 1 else str(HERE.parents[1] / "models" / "MACE-OFF24_medium.model")
OUT = Path(ARGV[2]) if len(ARGV) > 2 else HERE.parents[1] / "runs" / "p3_landscape_full"
env = guard(MODEL, label=f"landscape shard {SHARD}/{NSHARD}")
from prrs.calculators import mace_factory
from prrs.chemistry import chemical_key
from prrs.config import SearchConfig
from prrs.runner import curvature_spectrum, floor_residual
from prrs.state import encode

OUT.mkdir(parents=True, exist_ok=True)
factory = lambda: mace_factory(MODEL, device="cuda", default_dtype="float64")
cfg = SearchConfig(quench_fmax_eV_A=0.002)
calc = factory()
SMI = r"C\C(=O)\C=C(/O)c1ccccc1"          # STEREOZ; the E isomer cannot chelate at all

mol = Chem.MolFromSmiles(SMI)
assert [str(b.GetStereo()) for b in mol.GetBonds()
        if b.GetStereo() != Chem.BondStereo.STEREONONE] == ["STEREOZ"]
m = Chem.AddHs(mol)
ps = AllChem.ETKDGv3(); ps.randomSeed = 0xC0FFEE
ids = list(AllChem.EmbedMultipleConfs(m, numConfs=300, params=ps))
AllChem.MMFFOptimizeMoleculeConfs(m, maxIters=4000)
sym = [a.GetSymbol() for a in m.GetAtoms()]
OXY = [i for i, s in enumerate(sym) if s == "O"]
HYD = [i for i, s in enumerate(sym) if s == "H"]

# Conformers are independent, so the shard is a plain stride over the same ordered list.
# Stride rather than a contiguous block: the list is sorted by nothing in particular but
# ETKDG's later conformers tend to be the harder embeddings, and a block split would put
# all of them on one worker.
mine = [(n, cid) for n, cid in enumerate(ids) if n % NSHARD == SHARD]
print(f"shard {SHARD}/{NSHARD}: {len(mine)} of {len(ids)} conformers", flush=True)
SHARDJSON = OUT / f"landscape_shard{SHARD:02d}.json"

rows = []
t_start = time.time()
for k, (n, cid) in enumerate(mine):
    c = m.GetConformer(cid)
    P = np.array([list(c.GetAtomPosition(i)) for i in range(m.GetNumAtoms())])
    r_mmff = float(np.linalg.norm(P[OXY[0]] - P[OXY[1]]))
    atoms = Atoms(sym, positions=P); atoms.calc = calc
    t0 = time.time()
    opt = LBFGS(atoms, logfile=None); opt.run(fmax=2e-3, steps=2000)
    t_relax = time.time() - t0
    t0 = time.time()
    lam, _, floor, prov = curvature_spectrum(atoms, factory, cfg, source="analytic")
    t_hess = time.time() - t0
    neg = [float(x) for x in lam if x < -cfg.minimum_check_eigenvalue_tol]
    key = chemical_key(atoms)
    g = encode(atoms); edges = sorted(tuple(sorted(e)) for e in g.edges)
    oh = [(o, h) for o in OXY for h in HYD if tuple(sorted((o, h))) in edges]
    row = dict(conf=int(cid), r_OO_mmff=r_mmff,
               r_OO=float(np.linalg.norm(atoms.positions[OXY[0]] - atoms.positions[OXY[1]])),
               E_eV=float(atoms.get_potential_energy()),
               fmax=float(np.abs(atoms.get_forces()).max()),
               relax_steps=int(opt.get_number_of_steps()),
               n_minus=len(neg), lam1=float(lam[0]),
               floor_residual=float(floor_residual(floor)),
               n_OH=len(oh), n_edges=len(edges),
               key=key["key"][:16], graph_hash=key["graph_hash"],
               fragments=key["fragments"], stereo_unresolved=key["stereo_unresolved"],
               n_locked_parity=len(key["locked_bond_parity"]),
               t_relax_s=t_relax, t_hessian_s=t_hess)
    rows.append(row)
    write(str(OUT / f"conf{cid:04d}.extxyz"), atoms)
    eta = (time.time() - t_start) / (k + 1) * (len(mine) - k - 1) / 60
    print(f"[s{SHARD} {k+1:3d}/{len(mine)}] conf{cid:<4d} r_OO {r_mmff:.3f}->{row['r_OO']:.3f} "
          f"E {row['E_eV']:.6f} n- {row['n_minus']} key {row['key']} "
          f"({t_relax:.1f}+{t_hess:.1f}s, ETA {eta:.0f} min)", flush=True)
    json.dump(dict(env=env, smiles=SMI, atom_order="".join(sym),
                   shard=SHARD, nshard=NSHARD, rows=rows),
              open(SHARDJSON, "w"), indent=1, default=str)

# Per-shard summary only. The landscape question -- how many distinct confirmed minima --
# is global, so it is answered by p3_merge.py after every shard is in, never from one shard.
confirmed = [r for r in rows if r["n_minus"] == 0]
closed = [r for r in confirmed if r["r_OO"] < 3.0]
print(f"\n{len(rows)} relaxed, {len(confirmed)} confirmed minima, {len(closed)} closed")
if closed:
    Es = np.array([r["E_eV"] for r in closed])
    print(f"closed-chelate distinct energies at 0.1 meV: {len({round(e,4) for e in Es})}")
    print(f"closed-chelate distinct chemical_keys:       {len({r['key'] for r in closed})}")
    print(f"closed E span {(Es.max()-Es.min())*1e3:.2f} meV   lowest {Es.min():.6f} eV")
    op = [r for r in confirmed if r["r_OO"] >= 3.0]
    if op:
        print(f"open-form lowest {min(r['E_eV'] for r in op):.6f} eV  "
              f"chelate stabilisation {(min(r['E_eV'] for r in op)-Es.min())*1e3:.1f} meV")
print(f"shard {SHARD} wrote {SHARDJSON}; run p3_merge.py for the global count")
