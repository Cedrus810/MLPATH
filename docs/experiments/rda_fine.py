import sys, json; sys.path.insert(0, "src")
import numpy as np
from ase import Atoms
from ase.constraints import FixInternals
from ase.optimize import FIRE
from prrs.calculators import mace_factory
from prrs.state import encode

S = "/tmp/claude-1000/-home-ruigengji-MLPATH/f74de3f6-4a7e-4dbb-aed5-a3d30f6fe292/scratchpad/"
meta = open(S+"rda_meta.txt").read().split("\n")
b1, b2 = eval(meta[1]), eval(meta[2])
atoms = Atoms("C6H10", positions=np.load(S+"rda_pos.npy"))
calc = mace_factory("models/MACE-OFF24_medium.model", device="cuda", default_dtype="float64")
atoms.calc = calc
E_ref = float(atoms.get_potential_energy())

def set_length(a, pair, target):
    i, j = pair
    v = a.positions[j] - a.positions[i]; r = np.linalg.norm(v)
    s = 0.5 * (target - r) * v / r
    a.positions[i] -= s; a.positions[j] += s

# walk out to 1.90 the same way the coarse scan did, so the branch is the same one
for d in np.round(np.arange(1.55, 1.901, 0.05), 4):
    set_length(atoms, b1, float(d)); set_length(atoms, b2, float(d))
    atoms.set_constraint(FixInternals(bonds=[[float(d), list(b1)], [float(d), list(b2)]]))
    FIRE(atoms, logfile=None).run(fmax=2e-3, steps=1500); atoms.set_constraint()

geoms, rows = {}, []
for d in np.round(np.arange(1.90, 2.201, 0.01), 4):
    set_length(atoms, b1, float(d)); set_length(atoms, b2, float(d))
    atoms.set_constraint(FixInternals(bonds=[[float(d), list(b1)], [float(d), list(b2)]]))
    ok = FIRE(atoms, logfile=None).run(fmax=2e-3, steps=2000); atoms.set_constraint()
    E = float(atoms.get_potential_energy())
    geoms[float(d)] = atoms.positions.copy()
    rows.append(dict(d=float(d), rel_meV=(E-E_ref)*1e3, conv=bool(ok),
                     edges=len(encode(atoms).edges)))
    print(f"d={d:.2f}  E-E0={rows[-1]['rel_meV']:9.1f} meV  conv={int(ok)} "
          f"edges={rows[-1]['edges']}", flush=True)

e = np.array([r["rel_meV"] for r in rows])
print("\nsecond differences of the relaxed profile (meV, step 0.01 A):")
print(np.round(np.diff(e, 2), 1))
jump = int(np.argmax(np.abs(np.diff(e)))) if len(e) > 1 else 0
print(f"largest single step: d={rows[jump]['d']:.2f} -> {rows[jump+1]['d']:.2f}  "
      f"dE={e[jump+1]-e[jump]:+.1f} meV")

# rigid line between the two relaxed geometries that bracket the largest step: any kink
# here is the potential, not the optimizer switching branches.
A, B = geoms[rows[jump]["d"]], geoms[rows[jump+1]["d"]]
probe = atoms.copy(); probe.calc = calc
line = []
for t in np.linspace(0.0, 1.0, 101):
    probe.set_positions((1 - t) * A + t * B)
    line.append(float(probe.get_potential_energy()))
line = np.array(line)
print(f"\nrigid line between them, 101 points:")
print(f"  E range {(line.min()-E_ref)*1e3:.1f} .. {(line.max()-E_ref)*1e3:.1f} meV")
d2 = np.diff(line, 2) * 1e3
print(f"  |second difference| max {np.abs(d2).max():.4f} meV, mean {np.abs(d2).mean():.4f} meV")
print(f"  worst at t = {np.argmax(np.abs(d2))/100:.2f}")
json.dump(dict(E_ref=E_ref, rows=rows, rigid_line_eV=[float(x) for x in line]),
          open(S+"rda_fine.json","w"), indent=1)
