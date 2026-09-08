"""Source preparation for P2: A came back with one negative mode, so it is a stationary
point but not a minimum. Descend both sides of that mode and see whether an in-domain
A minimum exists at all.

Same calculator, same Hessian source, same relaxation configuration throughout. The
displacement is normalized by the largest single-atom displacement and started small. Both
sides get confirm_minimum rather than being judged by the optimizer's fmax, which is what
produced the wrong verdict in the first place. A lower-energy endpoint is only usable as A
if it is still A -- same chemical key, same proton host, same tautomer -- and if both sides
land in B or elsewhere, the benchmark does not have an A and that is the finding.
"""
import sys, json; sys.path.insert(0, "/home/ruigengji/MLPATH/src")
import numpy as np
from ase import Atoms
from ase.io import write
from ase.optimize import LBFGS
from prrs.calculators import mace_factory
from prrs.chemistry import chemical_key
from prrs.config import SearchConfig
from prrs.runner import confirm_minimum, curvature_spectrum, floor_residual
from prrs.state import encode

factory = lambda: mace_factory("/home/kasuga/.cache/mace/MACE-OFF24_medium.model",
                               device="cuda", default_dtype="float64")
CFG = SearchConfig(quench_fmax_eV_A=0.002)
FMAX, STEPS = 2e-3, 1000
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

def described(atoms, label):
    key = chemical_key(atoms)
    bonds = sorted(tuple(sorted(e)) for e in encode(atoms).edges)
    host = [b for b in bonds if 5 in b]
    return {"label": label, "energy_eV": float(atoms.get_potential_energy()),
            "fmax": float(np.abs(atoms.get_forces()).max()),
            "key": key["key"][:12], "graph_hash": key["graph_hash"],
            "fragments": key["fragments"], "stereo_unresolved": key["stereo_unresolved"],
            "proton_bonds": [list(b) for b in host],
            "O1_H": float(np.linalg.norm(atoms.positions[0]-atoms.positions[5])),
            "O2_H": float(np.linalg.norm(atoms.positions[4]-atoms.positions[5]))}

provenance = {}
A = build("O2"); A.calc = factory()
LBFGS(A, logfile=None).run(fmax=FMAX, steps=STEPS)
lam, vec, floor, prov = curvature_spectrum(A, factory, CFG, source="analytic")
negative = [float(x) for x in lam if x < -CFG.minimum_check_eigenvalue_tol]
provenance["stationary_point"] = {**described(A, "A_as_optimised"),
                                  "n_minus": len(negative), "negative_eigenvalues": negative,
                                  "floor_residual": floor_residual(floor)}
print(f"A as optimised: n_minus={len(negative)} lambda={negative} "
      f"E={A.get_potential_energy():.6f} fmax={np.abs(A.get_forces()).max():.2e}")
masses = A.get_masses()
mode = vec[:, 0].reshape(-1, 3) / np.sqrt(masses[:, None])
mode = mode / np.linalg.norm(mode, axis=1).max()          # largest single-atom step is 1 A
provenance["unstable_mode_max_atom_component"] = mode.tolist()
write("p2_A_saddle.extxyz", A.copy())

results, origin = [], A.positions.copy()
for delta in (0.05, 0.10, 0.20):
    for sign in (+1, -1):
        trial = A.copy(); trial.calc = factory()
        trial.set_positions(origin + sign * delta * mode)
        LBFGS(trial, logfile=None).run(fmax=FMAX, steps=STEPS)
        confirmed, diag = confirm_minimum(trial, factory, CFG, seed=0)
        entry = {**described(trial, f"delta={sign*delta:+.2f}"),
                 "confirmed_minimum": bool(confirmed),
                 "saddle_order": diag.get("saddle_order"),
                 "lowest_eigenvalues": [float(x) for x in (diag.get("lowest_eigenvalues") or [])[:3]],
                 "floor_residual": diag.get("trivial_mode_floor_worst_residual")}
        results.append((entry, trial))
        print(f"  delta={sign*delta:+.2f}  E={entry['energy_eV']:.6f} "
              f"confirmed={confirmed} order={entry['saddle_order']} "
              f"key={entry['key']} proton={entry['proton_bonds']} "
              f"O1-H={entry['O1_H']:.3f} O2-H={entry['O2_H']:.3f}")
provenance["descents"] = [e for e, _ in results]

reference = chemical_key(A)["key"]
same_species = [(e, t) for e, t in results
                if e["confirmed_minimum"] and e["key"] == reference[:12]
                and e["proton_bonds"] == [[4, 5]]]
print(f"\nconfirmed minima that are still A: {len(same_species)}")
if same_species:
    entry, best = min(same_species, key=lambda pair: pair[0]["energy_eV"])
    provenance["chosen"] = entry
    out = best.copy(); out.calc = None
    write("p2_input.extxyz", out)
    print(f"chosen: {entry['label']} E={entry['energy_eV']:.6f} -> p2_input.extxyz")
else:
    provenance["chosen"] = None
    print("NO in-domain A minimum: the P2 benchmark as posed does not have its reactant.")
json.dump(provenance, open("p2_source_provenance.json", "w"), indent=2)
