"""P2 preflight. Zero search budget: nothing is spent until these pass.

The benchmark's job is to check that a genuinely asymmetric rearrangement opens a second
chemical node, so the two tautomers must first be shown to be (a) the intended substance,
(b) real minima, and (c) distinguishable by the key. If any of that is untrue the run would
be testing the wrong thing.
"""
import sys; sys.path.insert(0, "/home/ruigengji/MLPATH/src")
import numpy as np
from ase import Atoms
from ase.io import write
from ase.optimize import LBFGS
from prrs.calculators import mace_factory
from prrs.chemistry import chemical_key, classify_transition, reaction_event_key
from prrs.config import SearchConfig
from prrs.runner import curvature_spectrum, floor_residual
from prrs.state import encode

factory = lambda: mace_factory("/home/kasuga/.cache/mace/MACE-OFF24_medium.model",
                               device="cuda", default_dtype="float64")
cfg = SearchConfig(quench_fmax_eV_A=0.002)
O1, O2 = np.array([-1.10, 1.15, 0.0]), np.array([1.24, 1.28, 0.0])
SK = {"O1": O1, "C1": np.array([-1.19,-0.05,0.0]), "C2": np.array([0.00,-0.80,0.0]),
      "C3": np.array([1.19,-0.05,0.0]), "O2": O2,
      "H1": np.array([-2.15,-0.55,0.0]), "H2": np.array([0.00,-1.88,0.0]),
      "H3": np.array([2.15,-0.55,0.0])}
NAMES = ["O1","C1","C2","C3","O2","Hb","H1","H2","Cm","Hm1","Hm2","Hm3"]

def build(on):
    host, other = (O1, O2) if on == "O1" else (O2, O1)
    pos = dict(SK); pos["Hb"] = host + 0.98*(other-host)/np.linalg.norm(other-host)
    axis = SK["H3"] - SK["C3"]; axis /= np.linalg.norm(axis)
    cm = SK["C3"] + 1.50*axis
    p1 = np.array([-axis[1], axis[0], 0.0]); p2 = np.array([0.0, 0.0, 1.0])
    hs = [cm + 1.09*(0.33*axis + 0.94*p1),
          cm + 1.09*(0.33*axis - 0.47*p1 + 0.82*p2),
          cm + 1.09*(0.33*axis - 0.47*p1 - 0.82*p2)]
    order = ["O1","C1","C2","C3","O2","Hb","H1","H2"]
    return Atoms("OCCCOHHH" + "C" + "HHH", positions=[pos[k] for k in order] + [cm] + hs)

results = {}
for tag, on in (("A (proton on O2, enol toward the ketone)", "O2"),
                ("B (proton on O1, enol toward the aldehyde)", "O1")):
    atoms = build(on); atoms.calc = factory()
    LBFGS(atoms, logfile=None).run(fmax=2e-3, steps=1000)
    graph = encode(atoms)
    bonds = sorted(tuple(sorted(e)) for e in graph.edges)
    key = chemical_key(atoms)
    lam, vec, floor, prov = curvature_spectrum(atoms, factory, cfg, source="analytic")
    negative = [x for x in lam if x < -cfg.minimum_check_eigenvalue_tol]
    print(f"\n=== {tag} ===")
    print(f"  fmax={np.abs(atoms.get_forces()).max():.2e}  E={atoms.get_potential_energy():.6f} eV")
    print(f"  key={key['key'][:12]}  graph_hash={key['graph_hash']}  frags={key['fragments']}")
    print(f"  stereo_unresolved={key['stereo_unresolved']}  locked_parity={len(key['locked_bond_parity'])}")
    print(f"  n_minus={len(negative)}  lowest={[round(float(x),5) for x in lam[:8]]}")
    print(f"  floor residual={floor_residual(floor):.2e}  source={prov['source']}")
    named = {f"{NAMES[a]}-{NAMES[b]}": round(float(np.linalg.norm(atoms.positions[a]-atoms.positions[b])),3)
             for a,b in bonds}
    print(f"  bonds: {named}")
    print(f"  O1..Hb={np.linalg.norm(atoms.positions[0]-atoms.positions[5]):.3f}  "
          f"O2..Hb={np.linalg.norm(atoms.positions[4]-atoms.positions[5]):.3f}  "
          f"O1..O2={np.linalg.norm(atoms.positions[0]-atoms.positions[4]):.3f}")
    results[on] = (atoms, key)

A, kA = results["O2"]; B, kB = results["O1"]
event = reaction_event_key(A, B)
print(f"\n=== 判定 ===")
print(f"  keys differ: {kA['key'] != kB['key']}   graph hashes differ: {kA['graph_hash'] != kB['graph_hash']}")
print(f"  broken={event['broken']}  formed={event['formed']}  event_key={event['key']}")
print(f"  classification: {classify_transition(kA['key'], kB['key'], event['broken'], event['formed'])}")
print(f"  dE(B-A) = {(B.get_potential_energy()-A.get_potential_energy())*1000:.2f} meV")
out = A.copy(); out.calc = None
write("p2_input.extxyz", out)
print("  wrote p2_input.extxyz (tautomer A as the search source)")
