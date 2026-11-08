"""Remaining fine-stage stops of HF-EF-bracket-v5 (TODO item 9c; handoff 2026-11-06 §5;
user 2026-11-04: "若仍卡在精阶段，就集中诊断剩余帧"). Offline diagnostic: ts_ef and the
runs are not changed; the Hessian below never enters a search.

Per item (index into the failed group of midpoint_coarse_validate.items), the whole HS is
rerun REPS times with v5 settings (as midpoint_fire_validate), and
  * every midpoint relaxation is recorded frame by frame (q, transverse / full force,
    reaction force, energy, positions);
  * the FIRE fine stage is instrumented per step: dt, alpha, velocity resets (P <= 0),
    |v|, fraction of v along the constraint normal at the new point, proposed step,
    maxstep clipping, and the fraction of the proposed step removed by the constraint
    correction (user 2026-11-06, item 2: velocity carried across a curved constraint).
For every fine-stage relaxation that did not converge, at its least-transverse-force
frame ("best") and its last frame:
  * NOISE repeated force evaluations at the same geometry (cache bypassed);
  * best frame only: the analytic Hessian projected on the transverse space (rigid body
    and grad q removed), the mean transverse force on its eigenvectors, the repeat noise
    along the softest mode, and NOISE evaluations at random transverse displacements of
    norm DISP_A with the Hessian-linear change removed (geometry-dependent roughness).

One JSON per (item, repeat); an existing file is skipped, so rerunning resumes.

    python midpoint_fine_diag.py OUT_DIR BOUND ITEMS_JSON SEED REPS ITEMS RESULTS_DIR

ITEMS_JSON is the item list of the 11-04 acceptance (midpoint_fire_validate items.json,
same intents and source runs); ITEMS is a comma list of failed-group indices there, e.g.
3,6,8,11,12. Only the source runs of those items are read, from RESULTS_DIR/<date>/<run>.
"""

import json
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT / "src"), str(ROOT / "explore" / "src")]
sys.path.insert(0, str(Path(__file__).resolve().parent))

from ase.calculators.calculator import all_changes  # noqa: E402
from ase.optimize import FIRE  # noqa: E402

import explore.ts_branch as ts_branch  # noqa: E402
import explore.ts_ef as ts_ef  # noqa: E402
from midpoint_coarse_validate import settings_for  # noqa: E402
from midpoint_cold_replay import context, load  # noqa: E402
from midpoint_transverse_hessian import spectrum  # noqa: E402

from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.d2 import _classify, _hf_bracket  # noqa: E402
from explore.d2 import _intent as intent_of  # noqa: E402
from explore.ledger import Ledger, view_factory  # noqa: E402
from prrs.runner import _hessian_providers  # noqa: E402

RELAX = ts_ef._relax
MIDPOINTS = ("from_left", "from_right", "from_left_coarse", "from_right_coarse")
NOISE = 16
DISP_A = 0.002
SOFT = 5
FIRE_LOGS = []


class LoggedFIRE(FIRE):
    """ASE FIRE, unchanged numerically; records its internal state after every step."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.diag = []
        FIRE_LOGS.append(self.diag)

    def step(self, f=None):
        opt = self.optimizable
        x0 = opt.get_x().copy()
        cons = opt.atoms.constraints[0]
        n0 = cons.gradient(opt.atoms.positions).ravel()
        n0 /= np.linalg.norm(n0)
        first, dt0 = self.vel is None, self.dt
        super().step(f)
        x1 = opt.get_x()
        dr = self.dt * self.vel
        size = float(np.linalg.norm(dr))
        clipped = size > self.maxstep
        if clipped:
            dr = dr * self.maxstep / size
        actual = x1 - x0
        n1 = cons.gradient(opt.atoms.positions).ravel()
        n1 /= np.linalg.norm(n1)
        vnorm = float(np.linalg.norm(self.vel))
        dnorm = max(float(np.linalg.norm(dr)), 1e-30)
        self.diag.append(
            {
                "dt": float(self.dt),
                "alpha": float(self.a),
                "reset": bool(not first and self.dt < dt0),
                "v_norm": vnorm,
                "v_normal_fraction": float(abs(self.vel @ n1) / max(vnorm, 1e-30)),
                "step_proposed_A": float(np.linalg.norm(dr)),
                "step_clipped": bool(clipped),
                "step_normal_fraction_old": float(abs(dr @ n0) / dnorm),
                "step_removed_fraction": float(np.linalg.norm(dr - actual) / dnorm),
                "step_actual_max_atom_A": float(
                    np.linalg.norm(actual.reshape(-1, 3), axis=1).max()
                ),
            }
        )


def frame_record(f):
    return {
        k: (f[k].tolist() if k == "positions" else f[k])
        for k in (
            "q_actual",
            "energy_eV",
            "reaction_force",
            "fmax_eV_A",
            "fmax_projected",
            "relax_step",
            "positions",
        )
    }


def rerun(run, task, settings, calc, out):
    row, points, ceiling = context(run, task)
    view = SimpleNamespace(settled=run.settled, settings=run.settings)
    classified = _classify(view, row)
    ledger = Ledger(out.with_suffix(".ledger.jsonl"), arm="HF", budget=run.hf.hs_cap)
    mids = []

    def watched(*args, **kwargs):
        origin, history = args[4], args[8]
        if origin not in MIDPOINTS:
            return RELAX(*args, **kwargs)
        n0, k0, stop = len(history), len(FIRE_LOGS), None
        try:
            return RELAX(*args, **kwargs)
        except ts_ef._Stop as exc:
            stop = str(exc)
            raise
        finally:
            fire = FIRE_LOGS[k0:]
            mids.append(
                {
                    "origin": origin,
                    "q": float(args[3]),
                    "target": args[10] if len(args) > 10 else kwargs.get("fmax"),
                    "optimizer": kwargs.get("optimizer", "bfgs"),
                    "converged": stop is None and len(history) > n0,
                    "stop": stop,
                    "_last": history[-1] if len(history) > n0 else None,
                    "frames": [frame_record(f) for f in history[n0:]],
                    "fire": fire[0] if fire else None,
                }
            )

    ts_ef._relax, ts_ef.FIRE = watched, LoggedFIRE
    t0 = time.monotonic()
    try:
        result = ts_branch.continue_branches_ef(
            run.poses[row["pose_id"]].copy(),
            intent_of(row),
            _hf_bracket(classified, points),
            view_factory(calc, ledger),
            run.config,
            ledger,
            task,
            branch_settings=ts_branch.BranchSettings(geometry_noise_A=run.hf.geometry_noise_A),
            settings=settings,
            ceiling_eV=ceiling,
        )
    finally:
        ts_ef._relax, ts_ef.FIRE = RELAX, FIRE
        ledger.close()
    for m in mids:  # refine_ts_ef sets coarse_sign on the returned frame afterwards
        last = m.pop("_last")
        m["coarse_sign"] = None if last is None else last.get("coarse_sign")
    refine = result.refine or {}
    return {
        "run": run.name,
        "task": task,
        "status": result.status,
        "spent": result.spent,
        "refine_spent": refine.get("spent"),
        "wall_s": time.monotonic() - t0,
        "candidate_energy_eV": refine.get("energy_eV"),
        "midpoints": mids,
    }


def evaluate(calc, atoms, positions):
    a = atoms.copy()
    a.positions = positions
    calc.calculate(a, ["energy", "forces"], all_changes)
    return float(calc.results["energy"]), np.asarray(calc.results["forces"], float).copy()


def transverse_basis(positions, grad_q):
    centered = positions - positions.mean(axis=0)
    rigid = []
    for axis in range(3):
        t = np.zeros_like(positions)
        t[:, axis] = 1.0
        rigid += [t.ravel(), np.cross(np.eye(3)[axis], centered).ravel()]
    return np.column_stack(rigid + [grad_q])


def project(forces, g):
    f = forces.ravel()
    return f - (f @ g) / (g @ g) * g


def atom_max(vec):
    return float(np.linalg.norm(vec.reshape(-1, 3), axis=1).max())


def noise_at(calc, atoms, positions, g, hessian=None, rng=None):
    rep = [evaluate(calc, atoms, positions) for _ in range(NOISE)]
    forces = np.array([f for _, f in rep])
    trans = np.array([project(f, g) for f in forces])
    out = {
        "energy_spread_eV": float(np.ptp([e for e, _ in rep])),
        "fproj_each": [atom_max(t) for t in trans],
        "fproj_of_mean": atom_max(trans.mean(axis=0)),
        "rf_each": [float(-(f.ravel() @ g) / (g @ g)) for f in forces],
        "fproj_noise_only": [atom_max(t - trans.mean(axis=0)) for t in trans],
        "component_std_max": float(forces.std(axis=0).max()),
        "component_std_median": float(np.median(forces.std(axis=0))),
    }
    if hessian is None:
        return out, trans
    n = len(positions)
    w, v = spectrum(hessian, transverse_basis(positions, g))
    coeff_mean = v.T @ trans.mean(axis=0)
    coeff_noise = (trans - trans.mean(axis=0)) @ v
    soft = v[:, 0]
    rough = []
    for _ in range(NOISE):
        d = v @ rng.normal(size=v.shape[1])
        d *= DISP_A / np.linalg.norm(d)
        _, f = evaluate(calc, atoms, positions + d.reshape(n, 3))
        rough.append(project(f, g) - (trans.mean(axis=0) - v @ (v.T @ (hessian @ d))))
    rough = np.array(rough)
    out |= {
        "transverse_eigs_lowest": w[:SOFT].tolist(),
        "mean_force_on_modes_lowest": coeff_mean[:SOFT].tolist(),
        "mean_force_norm": float(np.linalg.norm(coeff_mean)),
        "mean_force_fraction_soft": float(
            coeff_mean[0] ** 2 / max(coeff_mean @ coeff_mean, 1e-30)
        ),
        "mean_force_soft_only_atom_max": atom_max(soft * coeff_mean[0]),
        "mean_force_rest_atom_max": atom_max(trans.mean(axis=0) - soft * coeff_mean[0]),
        "noise_std_on_modes_lowest": coeff_noise[:, :SOFT].std(axis=0).tolist(),
        "noise_std_total": float(np.sqrt((coeff_noise**2).sum(axis=1).mean())),
        "soft_position_error_A": float(coeff_mean[0] / w[0]),
        "soft_noise_position_A": float(coeff_noise[:, 0].std() / w[0]),
        "roughness_fproj": [atom_max(r) for r in rough],
        "roughness_on_soft_std": float((rough @ soft).std()),
        "disp_A": DISP_A,
    }
    return out, trans


def hessian_at(calc, atoms, positions):
    a = atoms.copy()
    a.positions = positions
    a.calc = calc
    (provider,) = _hessian_providers(calc)
    scale = (
        float(getattr(provider, "energy_units_to_eV", 1.0))
        / float(getattr(provider, "length_units_to_A", 1.0)) ** 2
    )
    raw = scale * np.asarray(provider.get_hessian(a), dtype=float).reshape(3 * len(a), -1)
    return 0.5 * (raw + raw.T)


def diagnose(run, task, rec, calc, rng):
    row, _, _ = context(run, task)
    atoms = run.poses[row["pose_id"]].copy()
    atoms.set_constraint()
    intent = intent_of(row)
    out = []
    for j, m in enumerate(rec["midpoints"]):
        if m["optimizer"] != "fire" or m["converged"] or not m["frames"]:
            continue
        fr = m["frames"]
        best = int(np.argmin([f["fmax_projected"] for f in fr]))
        entry = {"midpoint": j, "origin": m["origin"], "best_frame": best}
        for label, k in (("best", best), ("last", len(fr) - 1)):
            pos = np.asarray(fr[k]["positions"])
            g = ts_ef._Constraint(intent, m["q"]).gradient(pos).ravel()
            hess = hessian_at(calc, atoms, pos) if label == "best" else None
            entry[label], _ = noise_at(calc, atoms, pos, g, hess, rng)
            entry[label]["recorded_fproj"] = fr[k]["fmax_projected"]
        out.append(entry)
    return out


def main():
    out_dir = Path(sys.argv[1])
    bound = float(sys.argv[2])
    picks = json.loads(Path(sys.argv[3]).read_text(encoding="utf-8"))
    seed, reps = map(int, sys.argv[4:6])
    wanted = [int(x) for x in sys.argv[6].split(",")]
    results = Path(sys.argv[7])
    failed = [p for p in picks if p[0] == "failed"]
    names = sorted({failed[i][2] for i in wanted})
    by_name = {r.name: r for r in (load(next(results.glob(f"*/{n}"))) for n in names)}
    out_dir.mkdir(parents=True, exist_ok=True)
    calc = polar1("cuda", "float32")
    rng = np.random.default_rng(seed)
    for i in wanted:
        _, task, name = failed[i]
        run = by_name[name]
        for k in range(reps):
            out = out_dir / f"item_{i:02d}_r{k}.json"
            if out.exists():
                continue
            rec = rerun(run, task, settings_for(run, "v5", bound), calc, out)
            t0 = time.monotonic()
            rec["diagnosis"] = diagnose(run, task, rec, calc, rng)
            rec["diagnosis_wall_s"] = time.monotonic() - t0
            rec["item"], rec["repeat"] = i, k
            out.write_text(json.dumps(rec) + "\n", encoding="utf-8")
            print(
                json.dumps(
                    {"file": out.name, "status": rec["status"], "spent": rec["spent"]}
                    | {"wall_s": round(rec["wall_s"], 1), "stuck": len(rec["diagnosis"])}
                ),
                flush=True,
            )


if __name__ == "__main__":
    main()
