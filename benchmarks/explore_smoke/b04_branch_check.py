"""Smoke: is the 2.3 -> 2.2 A step of the 10-13 b04 scan a continuous ridge or a branch jump?

The 10-13 scan (b04_scan.py) read ``ridge_suspected`` at q = 2.3 A, but on 2.3 -> 2.2 the
energy change and the trapezoid integral of the reaction force disagreed by 0.045 eV
(every other step <= 0.007 eV), and no geometry was saved. This script repeats the same
constrained-relaxation scan and saves every point, then:

  forward    3.6 -> 1.9 A, 0.1 A steps (same settings as 10-13)
  reverse    from the 1.9 A forward point back to 3.6 A, 0.1 A steps
  fine_fwd   from the 2.4 A forward point to 2.1 A, 0.01 A steps
  fine_rev   from the 2.1 A fine_fwd point back to 2.4 A, 0.01 A steps

Per step it records dE, the trapezoid integral of dE*/dq, their difference, the aligned
RMSD between consecutive points and the C0-Cl(leaving) distance. A branch jump shows up
as a step whose mismatch does not shrink with the step size, a jump in the leaving-group
distance or RMSD, and forward/reverse curves that differ at the same q (hysteresis).
No thresholds are applied here; the numbers are reported for reading.
"""

import json
import sys
import time
from dataclasses import asdict
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io import read, write

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "explore" / "src"), str(ROOT)]

from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.ledger import Ledger, counted_factory  # noqa: E402
from explore.probe import Intent, ScanSettings, constrained_scan  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402
from prrs.search import _stamped_with_charge  # noqa: E402
from prrs.state import encode  # noqa: E402

STEP = 0.1
FINE = 0.01
Q_END = 1.9
FINE_FROM, FINE_TO = 2.4, 2.1
FMAX, MAX_STEPS, MAX_RISE = 0.01, 300, 3.0


def _kabsch_rmsd(a, b):
    a = a - a.mean(axis=0)
    b = b - b.mean(axis=0)
    u, _, vt = np.linalg.svd(a.T @ b)
    d = np.sign(np.linalg.det(u @ vt))
    rot = u @ np.diag([1.0, 1.0, d]) @ vt
    return float(np.sqrt(np.mean(np.sum((a @ rot - b) ** 2, axis=1))))


def _deltas(q_start, q_stop, step):
    n = int(round(abs(q_stop - q_start) / step))
    sign = 1.0 if q_stop > q_start else -1.0
    return tuple(round(sign * step * k, 10) for k in range(n + 1))


def _settings(deltas):
    return ScanSettings(deltas, FMAX, MAX_STEPS, MAX_RISE)


def _closest(points, q):
    return min(points, key=lambda p: abs(p.q_target - q))


def _summarize(name, result, template, carbon, leaving):
    pts = result.points
    rows = []
    for k, p in enumerate(pts):
        row = {k2: v for k2, v in asdict(p).items() if k2 not in ("positions", "forces_pes")}
        row["index"] = k
        row["r_C_Cl_leaving"] = float(
            np.linalg.norm(p.positions[leaving] - p.positions[carbon])
        )
        if k:
            a = pts[k - 1]
            dq = p.q_actual - a.q_actual
            de = p.energy_eV - a.energy_eV
            trapz = 0.5 * (a.reaction_force + p.reaction_force) * dq
            row["step"] = {
                "dq": dq,
                "dE": de,
                "dE_trapz": trapz,
                "mismatch_eV": de - trapz,
                "rmsd_aligned_A": _kabsch_rmsd(a.positions, p.positions),
                "max_atom_disp_A": float(
                    np.max(np.linalg.norm(p.positions - a.positions, axis=1))
                ),
            }
        rows.append(row)
    frames = []
    for p in pts:
        f = Atoms(template.get_chemical_symbols(), positions=p.positions)
        f.info.update(
            q_target=p.q_target,
            energy_eV=p.energy_eV,
            reaction_force=p.reaction_force,
            converged=p.converged,
        )
        frames.append(f)
    return {"name": name, "status": result.status, "q0": result.q0, "points": rows}, frames


def main(out_dir):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=False)
    config = SearchConfig(charge_sensitive=True, total_charge=-1, multiplicity=1)
    atoms = _stamped_with_charge(
        read(ROOT / "benchmarks/b04_sn2_chloride/input.extxyz"), config
    )
    symbols = atoms.get_chemical_symbols()
    graph = encode(atoms, config.bond_scale)
    bonded = {i for e in graph.edges for i in e}
    nucleophile = next(i for i, s in enumerate(symbols) if s == "Cl" and i not in bonded)
    leaving = next(i for i, s in enumerate(symbols) if s == "Cl" and i in bonded)
    carbon = next(i for i, s in enumerate(symbols) if s == "C")
    ledger = Ledger(out / "ledger.jsonl", arm="smoke", budget=60000)
    factory = counted_factory(polar1, ledger)
    intent = Intent(f"cmp_{nucleophile}_{carbon}", ((nucleophile, carbon),), (1.0,))
    q0 = intent.value(atoms.positions)

    record = {
        "atoms": {"nucleophile": nucleophile, "carbon": carbon, "leaving": leaving},
        "settings": {"step": STEP, "fine": FINE, "fmax": FMAX, "max_steps": MAX_STEPS},
        "scans": [],
    }

    def run(name, start_positions, deltas):
        start = atoms.copy()
        start.positions = start_positions.copy()
        t0 = time.time()
        res = constrained_scan(
            start, intent, _settings(deltas), factory, config, ledger, task=name
        )
        summary, frames = _summarize(name, res, atoms, carbon, leaving)
        summary["wall_s"] = time.time() - t0
        summary["spent_after"] = ledger.spent
        write(out / f"{name}.extxyz", frames)
        record["scans"].append(summary)
        (out / "branch_check.json").write_text(json.dumps(record, indent=2, default=str))
        print(
            f"{name}: {res.status}, {len(res.points)} points, {summary['wall_s']:.1f} s",
            flush=True,
        )
        return res

    try:
        fwd = run("forward", atoms.positions, _deltas(q0, Q_END, STEP))
        if not fwd.points:
            return
        last = fwd.points[-1]
        run("reverse", last.positions, _deltas(last.q_actual, q0, STEP))
        a = _closest(fwd.points, FINE_FROM)
        fine = run("fine_fwd", a.positions, _deltas(a.q_actual, FINE_TO, FINE))
        if fine.points:
            b = fine.points[-1]
            run("fine_rev", b.positions, _deltas(b.q_actual, a.q_actual, FINE))
    finally:
        ledger.close()
        record["spent"] = ledger.spent
        (out / "branch_check.json").write_text(json.dumps(record, indent=2, default=str))


if __name__ == "__main__":
    main(sys.argv[1])
