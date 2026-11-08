"""How much does the reaction force dE/dq of a converged constrained minimum jitter?

Basis for the `branch_force_disagreement` gate in explore/ts_ef.py
(|RF_left - RF_right| > 2 * reaction_force_tol with transverse fmax <= 0.003 eV/A),
which stopped the b04 bracket at its first midpoint (results/2026-10-25/hf_ef_bounded/b04.json:
geometries 5.6e-4 A apart, RF -0.057 vs -0.065 eV/A).

Same quantity as the gate, at the same q, same model (float32 POLAR-1), same relaxation
(constraint as in ts_ef, ASE BFGS maxstep 0.05):
  A. rigid-body floor: RF at the converged geometry after 10 random rigid motions
     (physically invariant, so the spread is the numerical floor of RF itself);
  B. repeats: 8 seeded perturbations (sigma 0.01 A) of the converged geometry, each
     relaxed at fixed q to tau = 3e-3 and to 1e-3 eV/A; spread of RF across repeats;
  C. tightening: the two saved continuations relaxed further to tau = 1e-3, 3e-4, 1e-4;
     RF difference and RMSD between them at each tau.

    python hf_ef_rf_jitter.py OUT_JSON [CASE [SOURCE_JSON]]

SOURCE_JSON (optional): a hf_jump_bounded.py record; the point is then the LAST bisection
midpoint of its hand-off to refine_ts_ef (near the reaction-force root) instead.

CASE is a name from hf_ef_bounded.CASES (default b04). The point is the first bisection
midpoint of that case in results/2026-10-25/hf_ef_bounded/; its left continuation is the
base geometry. C runs only when the two continuations are on one branch (b05 E2: they
are not, so C is skipped and recorded as such).
"""

import json
import sys
import time
from pathlib import Path

import numpy as np
from ase.optimize import BFGS
from scipy.spatial.transform import Rotation

sys.path.insert(0, str(Path(__file__).resolve().parent))

from hf_ef_bounded import CASES, frozen  # noqa: E402  (also sets sys.path)

from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.d2 import _intent  # noqa: E402
from explore.ledger import Ledger, view_factory  # noqa: E402
from explore.ts_ef import _Constraint, _rmsd  # noqa: E402

OUT = Path(__file__).resolve().parent / "results/2026-10-25"
MAX_STEPS = 400


def rf_of(atoms, intent):
    forces = atoms.get_forces(apply_constraint=False)
    g = intent.gradient(atoms.positions)
    slope = -float(np.vdot(forces, g) / np.vdot(g, g))
    return slope, float(np.linalg.norm(forces + slope * g, axis=1).max())


def relax(template, positions, intent, q, tau, factory):
    work = template.copy()
    work.set_constraint()
    work.positions = np.array(positions, float)
    work.set_constraint(_Constraint(intent, q, tol=1e-8))
    work.set_positions(work.positions)
    work.calc = factory()
    opt = BFGS(work, logfile=None, maxstep=0.05)
    converged = opt.run(fmax=tau, steps=MAX_STEPS)
    rf, fproj = rf_of(work, intent)
    return {
        "tau": tau,
        "converged": bool(converged),
        "steps": opt.nsteps,
        "reaction_force": rf,
        "fmax_projected": fproj,
        "energy_eV": float(work.get_potential_energy()),
        "positions": work.positions.copy(),
    }


def strip(r):
    return {k: v for k, v in r.items() if k != "positions"}


def main():
    out = Path(sys.argv[1])
    want = sys.argv[2] if len(sys.argv) > 2 else "b04"
    name, run_dir, suffix = next(c for c in CASES if c[0] == want)
    settings, config, row, pose, _, _ = frozen(run_dir, suffix)
    intent = _intent(row)
    if len(sys.argv) > 3:
        saved = json.loads(Path(sys.argv[3]).read_text())
        mids = [h for h in saved["result"]["refine"]["history"] if "branch_rmsd_A" in h]
        left = next(h for h in reversed(mids) if h["origin"] == "from_left")
        right = next(h for h in reversed(mids) if h["origin"] == "from_right")
    else:
        saved = json.loads((OUT / "hf_ef_bounded" / f"{name}.json").read_text())
        mids = [h for h in saved["search"]["history"] if "branch_rmsd_A" in h]
        left = next(h for h in mids if h["origin"] == "from_left")
        right = next(h for h in mids if h["origin"] == "from_right")
    q = left["q_target"]
    t0 = time.monotonic()
    calc = polar1("cuda", "float32")
    ledger = Ledger(out.with_suffix(".ledger.jsonl"), arm="HF-jitter")
    factory = view_factory(calc, ledger)
    rng = np.random.default_rng(20261025)
    report = {"case": name, "q": q, "saved_rf": [left["reaction_force"], right["reaction_force"]]}

    # A. rigid-body floor
    base = pose.copy()
    base.set_constraint()
    base.calc = factory()
    rigid = []
    x0 = np.asarray(left["positions"], float)
    for _ in range(10):
        rot = Rotation.random(random_state=rng.integers(2**31)).as_matrix()
        base.positions = (x0 - x0.mean(0)) @ rot.T + x0.mean(0) + rng.normal(0, 1.0, 3)
        rigid.append(rf_of(base, intent)[0])
    report["A_rigid_rf"] = {"values": rigid, "spread": max(rigid) - min(rigid), "std": float(np.std(rigid))}
    print("A", report["A_rigid_rf"]["spread"], flush=True)

    # B. perturbed repeats
    report["B_repeats"] = {}
    for tau in (3e-3, 1e-3):
        runs = []
        for _ in range(8):
            start = x0 + rng.normal(0, 0.01, x0.shape)
            runs.append(relax(pose, start, intent, q, tau, factory))
        rfs = [r["reaction_force"] for r in runs]
        report["B_repeats"][str(tau)] = {
            "runs": [strip(r) for r in runs],
            "rf_spread": max(rfs) - min(rfs),
            "rf_std": float(np.std(rfs)),
            "max_pair_rmsd_A": max(
                _rmsd(a["positions"], b["positions"]) for a in runs for b in runs
            ),
        }
        print("B", tau, report["B_repeats"][str(tau)]["rf_spread"], flush=True)

    # C. tightening the two saved continuations (only if they are one branch)
    report["C_tighten"] = []
    if left["branch_rmsd_A"] > 0.05:
        report["C_skipped"] = f"continuations {left['branch_rmsd_A']:.3f} A apart: two branches"
        report["forwards"] = ledger.spent
        report["wall_s"] = time.monotonic() - t0
        out.write_text(json.dumps(report, indent=2, default=float) + "\n")
        return
    a_pos, b_pos = x0, np.asarray(right["positions"], float)
    for tau in (1e-3, 3e-4, 1e-4):
        a = relax(pose, a_pos, intent, q, tau, factory)
        b = relax(pose, b_pos, intent, q, tau, factory)
        a_pos, b_pos = a["positions"], b["positions"]
        report["C_tighten"].append(
            {
                "tau": tau,
                "left": strip(a),
                "right": strip(b),
                "d_rf": abs(a["reaction_force"] - b["reaction_force"]),
                "rmsd_A": _rmsd(a["positions"], b["positions"]),
            }
        )
        print("C", tau, report["C_tighten"][-1]["d_rf"], a["converged"], b["converged"], flush=True)

    report["forwards"] = ledger.spent
    report["wall_s"] = time.monotonic() - t0
    out.write_text(json.dumps(report, indent=2, default=float) + "\n")


if __name__ == "__main__":
    main()
