"""Why did p001's quench stop at fmax ~0.017 after 500 FIRE steps? Diagnostic, one structure.

Prints which atoms carry the residual force on the raw endpoint, then keeps minimizing the
same structure (FIRE, the quench's optimizer, and separately BFGS) and reports where each
lands: steps to fmax 0.01, energy drop, heavy-atom RMSD and the carbene markers.
Usage: python continue_quench.py <raw_endpoint.extxyz> <out dir>
"""
import sys
from pathlib import Path

import numpy as np
from ase.io import read, write
from ase.optimize import BFGS, FIRE

ROOT = Path(__file__).resolve().parents[5]
sys.path[:0] = [str(ROOT / "src"), str(ROOT), str(ROOT / "benchmarks/carbene_ambimodal")]
from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from classify_b09 import heavy_rmsd  # noqa: E402

start, out = read(sys.argv[1]), Path(sys.argv[2])
calc = polar1()
start.info.update(charge=0, spin=1)
start.calc = calc
f = start.get_forces()
norm = np.linalg.norm(f, axis=1)
s = start.get_chemical_symbols()
e_start = start.get_potential_energy()
print("raw endpoint E %.4f fmax %.4f" % (e_start, norm.max()))
for i in np.argsort(-norm)[:6]:
    print(f"  atom {i:2d} {s[i]:2s} |F| {norm[i]:.4f}")
for name, opt in (("FIRE", FIRE), ("BFGS", BFGS)):
    a = start.copy()
    a.info.update(charge=0, spin=1)
    a.calc = calc
    o = opt(a, logfile=None)
    ok = o.run(fmax=0.01, steps=3000)
    d = a.get_all_distances()
    # The calculator is shared: evaluate the end point last, so the results written with it
    # are its own (re-evaluating `start` here once wrote start's energy and forces into it).
    e_end = a.get_potential_energy()
    a.get_forces()
    print(f"{name}: converged {ok} steps {o.nsteps} E {e_end:.4f} "
          f"dE {e_end - e_start:+.4f} eV "
          f"RMSD {heavy_rmsd(start, a):.3f} A  C:-Ca {d[2,13]:.2f} C:-Cb {d[2,15]:.2f} "
          f"Ca-Cb {d[13,15]:.2f} C:-O_ald {d[2,1]:.2f}")
    write(out / f"{name.lower()}_end.extxyz", a)
