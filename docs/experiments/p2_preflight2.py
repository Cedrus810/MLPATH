"""Full P1'-P4' re-run after source preparation. Not just n_minus: fixing A moved it, and a
conformational change can move the stereochemistry key, the event key and the energy gap
too, so every gate is evaluated again on the structures that will actually be used."""
import sys, json; sys.path.insert(0, "/home/ruigengji/MLPATH/src")
import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.optimize import LBFGS
from prrs.calculators import mace_factory
from prrs.chemistry import chemical_key, classify_transition, reaction_event_key
from prrs.config import SearchConfig
from prrs.runner import confirm_minimum, floor_residual
from prrs.state import encode

factory = lambda: mace_factory("/home/kasuga/.cache/mace/MACE-OFF24_medium.model",
                               device="cuda", default_dtype="float64")
CFG = SearchConfig(quench_fmax_eV_A=0.002)
NAMES = ["O1","C1","C2","C3","O2","Hb","H1","H2","Cm","Hm1","Hm2","Hm3"]
O1, O2 = np.array([-1.10, 1.15, 0.0]), np.array([1.24, 1.28, 0.0])
SK = {"O1": O1, "C1": np.array([-1.19,-0.05,0.0]), "C2": np.array([0.00,-0.80,0.0]),
      "C3": np.array([1.19,-0.05,0.0]), "O2": O2,
      "H1": np.array([-2.15,-0.55,0.0]), "H2": np.array([0.00,-1.88,0.0]),
      "H3": np.array([2.15,-0.55,0.0])}
def build(on):
    host, other = (O1, O2) if on == "O1" else (O2, O1)
    pos = dict(SK); pos["Hb"] = host + 0.98*(other-host)/np.linalg.norm(other-host)
    axis = SK["H3"] - SK["C3"]; axis /= np.linalg.norm(axis)
    cm = SK["C3"] + 1.50*axis
    p1 = np.array([-axis[1], axis[0], 0.0]); p2 = np.array([0.0, 0.0, 1.0])
    hs = [cm + 1.09*(0.33*axis + 0.94*p1), cm + 1.09*(0.33*axis - 0.47*p1 + 0.82*p2),
          cm + 1.09*(0.33*axis - 0.47*p1 - 0.82*p2)]
    order = ["O1","C1","C2","C3","O2","Hb","H1","H2"]
    return Atoms("OCCCOHHH" + "C" + "HHH", positions=[pos[k] for k in order] + [cm] + hs)

A = read("p2_input.extxyz"); A.calc = factory()          # prepared source
B = build("O1"); B.calc = factory()
LBFGS(B, logfile=None).run(fmax=2e-3, steps=1000)
# B may sit on a saddle for the same reason A did; check and repair it the same way.
okB, diagB = confirm_minimum(B, factory, CFG, seed=0)
if not okB and diagB.get("saddle_order") == 1:
    masses = B.get_masses()
    mode = np.asarray(diagB["unstable_mode"])
    mode = mode / np.linalg.norm(mode, axis=1).max()
    best = None
    for sign in (+1, -1):
        trial = B.copy(); trial.calc = factory()
        trial.set_positions(B.positions + sign * 0.10 * mode)
        LBFGS(trial, logfile=None).run(fmax=2e-3, steps=1000)
        ok, _ = confirm_minimum(trial, factory, CFG, seed=0)
        if ok and (best is None or trial.get_potential_energy() < best.get_potential_energy()):
            best = trial
    if best is not None:
        print(f"B repaired the same way: {B.get_potential_energy():.6f} -> {best.get_potential_energy():.6f}")
        B = best
    okB, diagB = confirm_minimum(B, factory, CFG, seed=0)

report = {}
for tag, atoms in (("A", A), ("B", B)):
    ok, diag = confirm_minimum(atoms, factory, CFG, seed=0)
    key = chemical_key(atoms)
    bonds = sorted(tuple(sorted(e)) for e in encode(atoms).edges)
    report[tag] = {"E": float(atoms.get_potential_energy()),
                   "fmax": float(np.abs(atoms.get_forces()).max()),
                   "confirmed_minimum": bool(ok), "saddle_order": diag.get("saddle_order"),
                   "key": key["key"][:12], "graph_hash": key["graph_hash"],
                   "fragments": key["fragments"], "stereo": key["stereo_unresolved"],
                   "floor": diag.get("trivial_mode_floor_worst_residual"),
                   "proton": [list(b) for b in bonds if 5 in b],
                   "O1_H": float(np.linalg.norm(atoms.positions[0]-atoms.positions[5])),
                   "O2_H": float(np.linalg.norm(atoms.positions[4]-atoms.positions[5]))}
    print(f"\n=== {tag} ===")
    for k, v in report[tag].items():
        print(f"  {k}: {v}")

fwd = reaction_event_key(A, B); rev = reaction_event_key(B, A)
gates = {
    "P1'  both are confirmed minima":
        report["A"]["confirmed_minimum"] and report["B"]["confirmed_minimum"],
    "P2'  keys differ and graph hashes differ":
        report["A"]["key"] != report["B"]["key"] and report["A"]["graph_hash"] != report["B"]["graph_hash"],
    "P3'  both C4H6O2 and valence-resolved":
        report["A"]["fragments"] == report["B"]["fragments"] == ["C4H6O2"]
        and report["A"]["stereo"] is None and report["B"]["stereo"] is None,
    "P4'  event key direction-free, one O-H each way":
        fwd["key"] == rev["key"] and len(fwd["broken"]) == 1 and len(fwd["formed"]) == 1,
}
print("\n=== preflight gates ===")
for name, passed in gates.items():
    print(f"  {'PASS' if passed else 'FAIL'}  {name}")
print(f"\n  broken={fwd['broken']} formed={fwd['formed']} event_key={fwd['key']}")
print(f"  classification: {classify_transition(report['A']['key'], report['B']['key'], fwd['broken'], fwd['formed'])}")
print(f"  dE(B-A) = {(report['B']['E']-report['A']['E'])*1000:.2f} meV")
print(f"\n  ALL GATES {'PASS' if all(gates.values()) else 'FAIL'}")
out = A.copy(); out.calc = None; write("p2_input.extxyz", out)
outB = B.copy(); outB.calc = None; write("p2_product_reference.extxyz", outB)
json.dump({"report": report, "gates": {k: bool(v) for k, v in gates.items()},
           "event": {"broken": fwd["broken"], "formed": fwd["formed"], "key": fwd["key"]},
           "dE_meV": (report["B"]["E"]-report["A"]["E"])*1000},
          open("p2_preflight2.json", "w"), indent=2)
