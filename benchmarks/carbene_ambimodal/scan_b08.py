"""Does TS-CCl2-C4 exist on MACE-POLAR-1-M? Constrained scans along the forming C23-C20 bond.

b08 = :CCl2 + 1,2-bis(methylene)cycloheptane. In the two DFT saddles the carbene C23 attacks the
same C20=C16 bond at nearly the same distances (C23-C20 2.18 / 2.31 A, C23-C16 2.72 / 2.74 A);
they differ in where the carbene points: C23-C17 (C17 = the other terminal CH2, "C4" in the
paper) is 3.73 A in TS-CCl2 (ambimodal) and 4.23 A in TS-CCl2-C4. A plain geomeTRIC TS
search from the DFT TS-CCl2-C4 slid to a loose structure (C23-C20 2.93 A, -4 cm-1).

Here: from each DFT saddle, fix C23-C20 at d = 3.0 ... 1.9 A (FixBondLength), relax the rest
with BFGS on POLAR-1 (gas, charge 0, singlet), and record the energy and the orientation
markers C23-C16, C23-C17 at each point. A maximum along the distal (TS-CCl2-C4) branch that
keeps C23-C17 large is a TS guess; a monotonic descent says there is none along this path.
Each point starts from the previous relaxed point (a relaxed scan, not independent).

Also: symmetry-aware RMSD between the DFT saddles and the POLAR-1 re-optimised saddles,
minimising over the Cl24/Cl25 swap and mirror image (the enantiomeric approach).

Usage: python benchmarks/carbene_ambimodal/scan_b08.py <out dir> <opt dir>
"""

import itertools
import json
import sys
from pathlib import Path

import numpy as np
from ase.constraints import FixBondLength
from ase.io import read, write
from ase.optimize import BFGS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT)]
from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from prrs.state import centred_kabsch_rmsd  # noqa: E402

CARBENE, ATTACKED, PARTNER, DISTAL = 23, 20, 16, 17
CL = (24, 25)
EV_KCAL = 23.060548
GRID = np.round(np.arange(3.0, 1.85, -0.1), 2)


def rmsd_up_to_symmetry(a, b):
    """Kabsch RMSD minimised over the Cl swap and the mirror image of b."""
    best = np.inf
    for swap, mirror in itertools.product((False, True), (False, True)):
        order = list(range(len(b)))
        if swap:
            order[CL[0]], order[CL[1]] = CL[1], CL[0]
        y = b.positions[order] * (-1.0 if mirror else 1.0)
        x = a.positions - a.positions.mean(axis=0)
        best = min(best, centred_kabsch_rmsd(x, y - y.mean(axis=0)))
    return float(best)


def markers(atoms):
    d = atoms.get_all_distances()
    return dict(
        C23_C20=float(d[CARBENE, ATTACKED]),
        C23_C16=float(d[CARBENE, PARTNER]),
        C23_C17=float(d[CARBENE, DISTAL]),
    )


def scan(start, calc, out, label):
    atoms = start.copy()
    atoms.info.update(charge=0, spin=1)
    rows = []
    for d in GRID:
        atoms.set_constraint()
        atoms.set_distance(CARBENE, ATTACKED, d, fix=1)  # move the carbene end only
        atoms.set_constraint(FixBondLength(CARBENE, ATTACKED))
        atoms.calc = calc
        opt = BFGS(atoms, logfile=None)
        converged = opt.run(fmax=0.01, steps=500)
        row = dict(
            d=float(d),
            E_eV=float(atoms.get_potential_energy()),
            converged=bool(converged),
            steps=opt.get_number_of_steps(),
            **markers(atoms),
        )
        rows.append(row)
        write(out / f"{label}.extxyz", atoms, append=True)
        print(
            f"{label:9s} d {d:.2f}  E {row['E_eV']:.5f}  C23-C16 {row['C23_C16']:.3f}"
            f"  C23-C17 {row['C23_C17']:.3f}  steps {row['steps']}  conv {converged}",
            flush=True,
        )
    e0 = rows[0]["E_eV"]
    for row in rows:
        row["dE_kcal_vs_first"] = (row["E_eV"] - e0) * EV_KCAL
    return rows


def main():
    out, opt = Path(sys.argv[1]), Path(sys.argv[2])
    out.mkdir(parents=True, exist_ok=True)
    s = HERE / "structures" / "ts"
    geoms = {
        "DFT TS-CCl2": read(s / "TS-CCl2.xyz"),
        "DFT TS-CCl2-C4": read(s / "TS-CCl2-C4.xyz"),
        "P1 opt from TS-CCl2": read(opt / "TS-CCl2" / "geometry.extxyz"),
        "P1 opt from TS-CCl2-C4": read(opt / "TS-CCl2-C4" / "geometry.extxyz"),
    }
    print("symmetry-aware RMSD (Cl swap, mirror), A")
    for (ka, a), (kb, b) in itertools.combinations(geoms.items(), 2):
        print(f"  {ka:24s} {kb:24s} {rmsd_up_to_symmetry(a, b):.3f}")
    calc = polar1()
    result = {
        "proximal": scan(geoms["DFT TS-CCl2"], calc, out, "proximal"),
        "distal": scan(geoms["DFT TS-CCl2-C4"], calc, out, "distal"),
    }
    print("\nrelaxed-scan energies, kcal/mol relative to each branch's d = 3.0 A point")
    for d_index, d in enumerate(GRID):
        p, q = result["proximal"][d_index], result["distal"][d_index]
        print(
            f"  d {d:.2f}   proximal {p['dE_kcal_vs_first']:7.2f} (C23-C17 {p['C23_C17']:.2f})"
            f"   distal {q['dE_kcal_vs_first']:7.2f} (C23-C17 {q['C23_C17']:.2f})"
            f"   distal-proximal {(q['E_eV'] - p['E_eV']) * EV_KCAL:6.2f}"
        )
    (out / "scan.json").write_text(json.dumps(result, indent=1))


if __name__ == "__main__":
    main()
