"""Same-branch reaction-force jitter sigma and geometry noise for a new system, from its own
L1 scans (docs/PLAN_HF5_CROSS_SYSTEM_2026-11-08.md §5 n1). Generalises hf_ef_rf_jitter.py,
which needs a frozen 10-25 run; the measured quantities are the same:
  A. RF at a converged constrained point after 10 random rigid motions (floor of RF itself);
  B. 8 seeded perturbations (N(0, 0.01 A)) of that point, each relaxed at the same q with the
     ts_ef constraint and ASE BFGS maxstep 0.05 to tau = 3e-3 and to 1e-3 eV/A: std and
     range of RF (-> HFSettings.reaction_force_sigma) and the largest pairwise RMSD of the
     relaxed geometries (-> HFSettings.geometry_noise_A).
Points: L0 as run_d2.py (first placed pose); every distance intent of that pose is scanned
over the full magnitude grid; scans whose highest converged energy is at an interior point
AND whose own `edges` differ between that peak and the next converged point (a bond-changing
ridge, as `jump` / _same_graph_ridge; user 2026-11-09) are kept, N_SCANS of them in a seeded
order; on each, the peak point and the point before it. No threshold picks the points; which
statistic becomes sigma is decided afterwards. (First version, 11-09, took any interior peak:
on N1 that picked only loose H / O contacts -- see DONE_LOG.)

    python rf_jitter_scan.py OUT_JSON SYSTEM SEED SHELL N_SAMPLES N_FAR N_SCANS
"""

import json
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from ase.optimize import BFGS
from scipy.spatial.transform import Rotation

sys.path.insert(0, str(Path(__file__).resolve().parent))

import run_d2  # noqa: E402  (sets sys.path)

from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.intents import distance_intents  # noqa: E402
from explore.ledger import Ledger, view_factory  # noqa: E402
from explore.parallel import MACEBatchEvaluator, ScanJob, scan_many  # noqa: E402
from explore.streams import stream  # noqa: E402
from explore.ts_ef import _Constraint, _rmsd  # noqa: E402

MAX_STEPS = 400
REPEATS = 8
SIGMA_A = 0.01
TAUS = (3e-3, 1e-3)


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
        "converged": bool(converged),
        "steps": opt.nsteps,
        "reaction_force": rf,
        "fmax_projected": fproj,
        "positions": work.positions.copy(),
    }


def measure(pose, intent, x0, q, factory, rng):
    rec = {"q": q}
    base = pose.copy()
    base.set_constraint()
    base.calc = factory()
    rigid = []
    for _ in range(10):
        rot = Rotation.random(random_state=rng.integers(2**31)).as_matrix()
        base.positions = (x0 - x0.mean(0)) @ rot.T + x0.mean(0) + rng.normal(0, 1.0, 3)
        rigid.append(rf_of(base, intent)[0])
    rec["A_rigid"] = {"std": float(np.std(rigid)), "range": max(rigid) - min(rigid)}
    rec["B"] = {}
    for tau in TAUS:
        runs = [
            relax(pose, x0 + rng.normal(0, SIGMA_A, x0.shape), intent, q, tau, factory)
            for _ in range(REPEATS)
        ]
        ok = [r for r in runs if r["converged"]]
        rfs = [r["reaction_force"] for r in ok]
        rec["B"][str(tau)] = {
            "converged": len(ok),
            "rf_std": float(np.std(rfs)) if len(rfs) > 1 else None,
            "rf_range": (max(rfs) - min(rfs)) if len(rfs) > 1 else None,
            "rf_mean": float(np.mean(rfs)) if rfs else None,
            "max_pair_rmsd_A": max(
                (_rmsd(a["positions"], b["positions"]) for a in ok for b in ok), default=None
            ),
            "steps": [r["steps"] for r in runs],
        }
    return rec


def main():
    out = Path(sys.argv[1])
    name, seed = sys.argv[2], int(sys.argv[3])
    shell, n_samples, n_far, n_scans = (
        float(sys.argv[4]),
        int(sys.argv[5]),
        int(sys.argv[6]),
        int(sys.argv[7]),
    )
    system = run_d2.SYSTEMS[name]
    config = run_d2.make_config(system)
    t0 = time.monotonic()
    calc = polar1("cuda", "float32")
    ledger_path = out.with_suffix(".ledger.jsonl")
    atoms, _ = run_d2.make_atoms(system, config)
    args = SimpleNamespace(shell=[shell], n_samples=n_samples, seed=seed)
    poses, l0 = run_d2.l0_for(
        args, system, atoms, config, calc, out.with_suffix(".prepare.jsonl")
    )
    pose_id, pose = poses[0]
    specs = distance_intents(pose, pose_id, config, n_far, seed)[0]
    mags = run_d2.magnitudes(0.1)
    jobs = [
        ScanJob(
            pose,
            s.intent,
            run_d2.scan_settings(s.deltas(mags)),
            config,
            Ledger(ledger_path, "jitter-scan", replay=False),
            f"scan|{s.intent.intent_id}",
        )
        for s in specs
    ]
    scans, _ = scan_many(jobs, MACEBatchEvaluator(calc), 8)
    peaked = []
    for spec, sc in zip(specs, scans):
        pts = [p for p in (sc.result.points if sc.result else []) if p.converged]
        if len(pts) < 3:
            continue
        k = int(np.argmax([p.energy_eV for p in pts]))
        if 0 < k < len(pts) - 1:
            peaked.append((spec, pts, k))
    peaked.sort(key=lambda t: float(stream(seed, "jitter", t[0].intent.intent_id).random()))
    # user 2026-11-09: only bond-changing ridges -- the scan's own edges differ between the
    # peak and the next converged point (same definition as `jump` and _same_graph_ridge)
    changing = [t for t in peaked if t[1][t[2]].edges != t[1][t[2] + 1].edges]
    report = {
        "system": name,
        "seed": seed,
        "pose_id": pose_id,
        "l0": l0,
        "intents_scanned": len(specs),
        "scans_with_interior_peak": len(peaked),
        "peaked_scans": [
            {
                "intent_id": s.intent.intent_id,
                "bond_change_at_peak": pts[k].edges != pts[k + 1].edges,
                "peak_q": pts[k].q_target,
            }
            for s, pts, k in peaked
        ],
        "scans_with_bond_changing_peak": len(changing),
        "points": [],
    }
    ledger = Ledger(ledger_path, "jitter", replay=False)
    factory = view_factory(calc, ledger)
    rng = np.random.default_rng(seed)
    for spec, pts, k in changing[:n_scans]:
        for label, j in (("before_peak", k - 1), ("peak", k)):
            p = pts[j]
            rec = measure(
                pose, spec.intent, np.asarray(p.positions, float), p.q_target, factory, rng
            )
            rec |= {
                "intent_id": spec.intent.intent_id,
                "where": label,
                "scan_rf": p.reaction_force,
            }
            report["points"].append(rec)
            print(label, spec.intent.intent_id[-30:], json.dumps(rec["B"])[:200], flush=True)
    ledger.close()
    stats = {}
    for tau in TAUS:
        b = [r["B"][str(tau)] for r in report["points"]]
        stats[str(tau)] = {
            "rf_std_max": max(
                (x["rf_std"] for x in b if x["rf_std"] is not None), default=None
            ),
            "rf_std_median": float(
                np.median([x["rf_std"] for x in b if x["rf_std"] is not None])
            )
            if b
            else None,
            "max_pair_rmsd_A_max": max(
                (x["max_pair_rmsd_A"] for x in b if x["max_pair_rmsd_A"] is not None),
                default=None,
            ),
        }
    report["summary"] = stats
    report["wall_s"] = time.monotonic() - t0
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=1, default=float) + "\n", encoding="utf-8")
    print(json.dumps(stats, indent=1))


if __name__ == "__main__":
    main()
