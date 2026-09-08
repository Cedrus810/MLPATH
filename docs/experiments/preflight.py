"""Preflight for the P1/P2 benchmarks: which cell of the classification table do they land in?

Zero search budget. Relax both tautomers, check the connectivity is what was intended, then
compare chemical_key (state identity) against the reaction event (path identity).
"""
import sys; sys.path.insert(0, "/home/ruigengji/MLPATH/src")
import numpy as np
from ase import Atoms
from ase.optimize import LBFGS
from prrs.calculators import mace_factory
from prrs.chemistry import chemical_key, reaction_event_key, classify_transition
from prrs.state import encode

factory = lambda: mace_factory("/home/kasuga/.cache/mace/MACE-OFF24_medium.model",
                               device="cuda", default_dtype="float64")

O1, O2 = np.array([-1.10, 1.15, 0.0]), np.array([1.24, 1.28, 0.0])
SK = {"O1": O1, "C1": np.array([-1.19,-0.05,0.0]), "C2": np.array([0.00,-0.80,0.0]),
      "C3": np.array([1.19,-0.05,0.0]), "O2": O2,
      "H1": np.array([-2.15,-0.55,0.0]), "H2": np.array([0.00,-1.88,0.0]),
      "H3": np.array([2.15,-0.55,0.0])}

def malonaldehyde(on):
    host, other = (O1, O2) if on == "O1" else (O2, O1)
    pos = dict(SK); pos["Hb"] = host + 0.98*(other-host)/np.linalg.norm(other-host)
    order = ["O1","C1","C2","C3","O2","Hb","H1","H2","H3"]
    return Atoms("OCCCOHHHH", positions=[pos[k] for k in order])

def oxobutanal(on):
    """Malonaldehyde with H3 replaced by a methyl: the two ends are no longer equivalent."""
    host, other = (O1, O2) if on == "O1" else (O2, O1)
    pos = dict(SK); pos["Hb"] = host + 0.98*(other-host)/np.linalg.norm(other-host)
    axis = (SK["H3"] - SK["C3"]); axis /= np.linalg.norm(axis)
    cm = SK["C3"] + 1.50*axis
    perp1 = np.array([-axis[1], axis[0], 0.0]); perp2 = np.array([0.0, 0.0, 1.0])
    hs = [cm + 1.09*( 0.33*axis + 0.94*perp1),
          cm + 1.09*( 0.33*axis - 0.47*perp1 + 0.82*perp2),
          cm + 1.09*( 0.33*axis - 0.47*perp1 - 0.82*perp2)]
    order = ["O1","C1","C2","C3","O2","Hb","H1","H2"]
    symbols = "OCCCOHHH" + "C" + "HHH"
    return Atoms(symbols, positions=[pos[k] for k in order] + [cm] + hs)

def relax(atoms, label):
    atoms.calc = factory()
    LBFGS(atoms, logfile=None).run(fmax=2e-3, steps=800)
    print(f"  {label}: fmax={np.abs(atoms.get_forces()).max():.2e}  E={atoms.get_potential_energy():.6f} eV")
    return atoms

def report(name, builder, proton_index, oxygens):
    print(f"\n=== {name} ===")
    A = relax(builder("O2"), "A (H on O2)")
    B = relax(builder("O1"), "B (H on O1)")
    for tag, m in (("A", A), ("B", B)):
        bonds = sorted(b for b in (tuple(sorted(e)) for e in encode(m).edges)
                       if proton_index in b)
        print(f"  {tag} proton bonds: {bonds}  (O1={oxygens[0]}, O2={oxygens[1]})")
    ka, kb = chemical_key(A), chemical_key(B)
    event = reaction_event_key(A, B)
    print(f"  graph_hash  A={ka['graph_hash']}  B={kb['graph_hash']}  same={ka['graph_hash']==kb['graph_hash']}")
    print(f"  chemical_key same: {ka['key'] == kb['key']}")
    print(f"  locked_bond_parity A={len(ka['locked_bond_parity'])} B={len(kb['locked_bond_parity'])}")
    print(f"  broken={event['broken']}  formed={event['formed']}")
    print(f"  event key fwd={event['key']}  rev={reaction_event_key(B, A)['key']}")
    print(f"  >>> classification: {classify_transition(ka['key'], kb['key'], event['broken'], event['formed'])}")
    print(f"  energy difference B - A = {(B.get_potential_energy()-A.get_potential_energy())*1000:.2f} meV")

report("P1 malonaldehyde (symmetric)", malonaldehyde, 5, (0, 4))
report("P2 3-oxobutanal enol (asymmetric)", oxobutanal, 5, (0, 4))
