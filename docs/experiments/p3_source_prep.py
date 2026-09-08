"""Which coordinate carries A's negative eigenvalue, and is the lower conformer a minimum?

p3_chelate_search left two things open:
  - A (tightest chelate) has n_minus = 1, lam_1 = -0.003183 eV/A^2, three times the
    minimum_check_eigenvalue_tol of 1e-3. Small next to P2's -0.067, but the P2 preflight
    gate P1' exists precisely to stop a saddle being used as a reactant.
  - eight sampled conformers relaxed to E = -14637.363922 eV, 1.50 meV BELOW A, so A is
    not the lowest closed conformer in the sample.
Both are answered by descending A along its soft negative mode, the P2 recipe: same
calculator, same analytic Hessian, same LBFGS tolerance, displacement normalised on the
largest single-atom move, and confirm_minimum on every endpoint rather than trusting fmax.
Zero search budget.
"""
import sys, json
sys.path.insert(0, "/home/ruigengji/MLPATH/src")
import numpy as np
from ase.io import read, write
from ase.optimize import LBFGS
from prrs.calculators import mace_factory
from prrs.chemistry import chemical_key, classify_transition, reaction_event_key
from prrs.config import SearchConfig
from prrs.runner import curvature_spectrum, floor_residual

MODEL = "/home/kasuga/.cache/mace/MACE-OFF24_medium.model"
factory = lambda: mace_factory(MODEL, device="cuda", default_dtype="float64")
cfg = SearchConfig(quench_fmax_eV_A=0.002)
calc = factory()

A = read("/home/ruigengji/MLPATH/runs/p3_A_chelated.extxyz")
A.calc = calc
sym = list(A.symbols)
OXY = [i for i, s in enumerate(sym) if s == "O"]
lam, vec, floor, prov = curvature_spectrum(A, factory, cfg, source="analytic")
neg_idx = [i for i, x in enumerate(lam) if x < -cfg.minimum_check_eigenvalue_tol]
print(f"A: lam[:8] = {[round(float(x),6) for x in lam[:8]]}")
print(f"   negative modes {neg_idx}  floor residual {floor_residual(floor):.2e}")

# Name the soft negative mode by where its amplitude sits, per atom.
# curvature_spectrum returns `vectors[:, order]`, so eigenvectors are COLUMNS, and they
# live in the mass-weighted basis. An earlier version of this script read vec[i] as mode i
# -- that is row i, a slice across all modes, and the atom distribution it printed was
# meaningless. Cartesian displacement needs the mass weighting undone: dx = v / sqrt(m).
v_mw = np.asarray(vec[:, neg_idx[0]]).reshape(len(A), 3)
v = v_mw / np.sqrt(A.get_masses())[:, None]
v /= np.linalg.norm(v)
amp = np.linalg.norm(v, axis=1)
order = np.argsort(amp)[::-1]
print("   negative-mode amplitude by atom (top 8):")
for i in order[:8]:
    print(f"     atom {i:2d} {sym[i]}  |v| = {amp[i]:.4f}")
# Groups are read off A's own graph, not hard-coded index ranges. An earlier version
# defined "the methyl" as atom 0 alone -- the carbon without its three hydrogens -- and so
# reported 0.000 for a mode that is almost entirely those three hydrogens.
from prrs.state import encode as _encode
_edges = sorted(tuple(sorted(e)) for e in _encode(A).edges)
def _nbrs(i):
    return [b if a == i else a for a, b in _edges if i in (a, b)]
ring = [i for i in range(len(A)) if sym[i] == "C" and
        sum(1 for j in _nbrs(i) if sym[j] == "C") >= 2 and
        len([j for j in _nbrs(i) if sym[j] == "H"]) <= 1 and i >= 6]
methyl_c = [i for i in range(len(A)) if sym[i] == "C"
            and len([j for j in _nbrs(i) if sym[j] == "H"]) == 3]
methyl = methyl_c + [j for i in methyl_c for j in _nbrs(i) if sym[j] == "H"]
total = float((amp**2).sum())
print(f"   fraction of |v|^2 on the phenyl ring {ring}: "
      f"{float((amp[ring]**2).sum())/total:.3f}")
print(f"   fraction on the methyl group {sorted(methyl)}: "
      f"{float((amp[sorted(methyl)]**2).sum())/total:.3f}")
print(f"   fraction on the shared proton (atom {[i for i in range(len(A)) if sym[i]=='H']}"
      f" -> H16): {float(amp[16]**2)/total:.3f}")

def confirm(atoms, tag):
    lam_, vec_, floor_, prov_ = curvature_spectrum(atoms, factory, cfg, source="analytic")
    neg = [float(x) for x in lam_ if x < -cfg.minimum_check_eigenvalue_tol]
    key = chemical_key(atoms)
    d = dict(tag=tag, E_eV=float(atoms.get_potential_energy()),
             fmax=float(np.abs(atoms.get_forces()).max()),
             n_minus=len(neg), lam6=[round(float(x), 6) for x in lam_[:6]],
             floor_residual=float(floor_residual(floor_)),
             r_OO=float(np.linalg.norm(atoms.positions[OXY[0]] - atoms.positions[OXY[1]])),
             key=key["key"][:16], graph_hash=key["graph_hash"],
             fragments=key["fragments"], stereo_unresolved=key["stereo_unresolved"],
             n_locked_parity=len(key["locked_bond_parity"]))
    print(f"\n[{tag}]")
    for k, val in d.items(): print(f"   {k:18s} {val}")
    return d

# Descend both ways along the negative mode, three amplitudes, largest-atom-move normalised.
step = v / np.abs(v).max()
endpoints = {}
for sign in (+1, -1):
    for a in (0.05, 0.10, 0.20):
        cand = A.copy()
        cand.positions += sign * a * step
        cand.calc = factory()
        LBFGS(cand, logfile=None).run(fmax=2e-3, steps=2000)
        tag = f"desc{'+' if sign>0 else '-'}{a:.2f}"
        endpoints[tag] = (confirm(cand, tag), cand)

# The lowest confirmed minimum on this side becomes the P3 source candidate.
mins = {t: (d, at) for t, (d, at) in endpoints.items() if d["n_minus"] == 0}
print(f"\nconfirmed minima among the six descents: {sorted(mins)}")
if mins:
    best_tag = min(mins, key=lambda t: mins[t][0]["E_eV"])
    dbest, abest = mins[best_tag]
    print(f"lowest confirmed minimum: {best_tag}  E = {dbest['E_eV']:.6f} eV  "
          f"r_OO = {dbest['r_OO']:.3f} A")
    print(f"  vs A (saddle)            E = {A.get_potential_energy():.6f} eV  "
          f"drop = {(dbest['E_eV'] - A.get_potential_energy())*1e3:.3f} meV")
    print(f"  matches the 1.50 meV lower conformer seen in the landscape sample: "
          f"{abs(dbest['E_eV'] - (-14637.363922)) < 5e-4}")
    write("/home/ruigengji/MLPATH/runs/p3_A_source.extxyz", abest)

    # Rebuild B from the confirmed source, same ordering, and re-measure the event.
    g = __import__("prrs.state", fromlist=["encode"]).encode(abest)
    edges = sorted(tuple(sorted(e)) for e in g.edges)
    hyd = [i for i, s in enumerate(sym) if s == "H"]
    oh = [(o, h) for o in OXY for h in hyd if tuple(sorted((o, h))) in edges]
    assert len(oh) == 1, oh
    donor, h_sh = oh[0]; acc = [o for o in OXY if o != donor][0]
    B = abest.copy()
    u = abest.positions[donor] - abest.positions[acc]
    B.positions[h_sh] = abest.positions[acc] + 0.99 * u / np.linalg.norm(u)
    B.calc = factory()
    LBFGS(B, logfile=None).run(fmax=2e-3, steps=2000)
    dB = confirm(B, "B_from_source")
    cls = classify_transition(abest, B, cfg.bond_scale, cfg.active_atoms)
    evf = reaction_event_key(abest, B, cfg.bond_scale, cfg.active_atoms)
    evr = reaction_event_key(B, abest, cfg.bond_scale, cfg.active_atoms)
    ok_p1 = dbest["n_minus"] == 0 and dB["n_minus"] == 0
    ok_p2 = dbest["key"] != dB["key"] and dbest["graph_hash"] != dB["graph_hash"]
    ok_p3 = (dbest["fragments"] == ["C10H10O2"] and dB["fragments"] == ["C10H10O2"]
             and dbest["stereo_unresolved"] is None and dB["stereo_unresolved"] is None)
    ok_p4 = (evf["key"] == evr["key"] and len(evf["broken"]) == 1
             and len(evf["formed"]) == 1)
    ok_p5 = dbest["r_OO"] < 3.0 and dB["r_OO"] < 3.0
    print("\n--- P2-style preflight on the confirmed pair ---")
    print(f"P1'  both n_minus == 0                          {ok_p1}")
    print(f"P2'  key and graph_hash differ                  {ok_p2}")
    print(f"P3'  both C10H10O2, valence-resolved            {ok_p3}")
    print(f"P4'  event key direction-invariant, one O-H     {ok_p4}")
    print(f"     classify={cls}  broken={evf['broken']} formed={evf['formed']}")
    print(f"     event_key={evf['key']}")
    print(f"P5'  chelate closed both ends                   {ok_p5}"
          f"  ({dbest['r_OO']:.3f} / {dB['r_OO']:.3f} A)")
    print(f"     dE(B-A) = {(dB['E_eV']-dbest['E_eV'])*1e3:.3f} meV")
    write("/home/ruigengji/MLPATH/runs/p3_B_source.extxyz", B)
    json.dump(dict(A_saddle_lam1=float(lam[neg_idx[0]]),
                   negative_mode_ring_fraction=float((amp[ring]**2).sum()/(amp**2).sum()),
                   negative_mode_top_atoms=[[int(i), sym[i], float(amp[i])] for i in order[:8]],
                   descents={t: d for t, (d, _) in endpoints.items()},
                   source=dbest, product=dB, classify=str(cls),
                   event_forward=evf["key"], event_reverse=evr["key"],
                   broken=evf["broken"], formed=evf["formed"],
                   dE_meV=(dB["E_eV"]-dbest["E_eV"])*1e3),
              open("/home/ruigengji/MLPATH/runs/p3_source_prep.json", "w"),
              indent=2, default=str)
    print("\nwrote runs/p3_source_prep.json, p3_A_source.extxyz, p3_B_source.extxyz")
