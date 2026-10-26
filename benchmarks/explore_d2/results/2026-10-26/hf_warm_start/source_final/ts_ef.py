"""E/F-only stationary-candidate locator on a constrained, continuous branch.

HF-2/3 primitive; deliberately separate from C2 and the D2-v1.1 scheduler.
The full Hessian and endpoint connections are subsequent, independent checks.
Continuity checks below are bounded numerical diagnostics, not a proof that a
transverse minimum exists. No fallback to eigenvector following is performed.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
import time

import numpy as np
from ase.optimize import BFGS

from .ledger import BudgetExhausted
from .poses import _jsonable
from .probe import FixBondCombo


@dataclass(frozen=True)
class EFSettings:
    fmax_eV_A: float = 0.01
    transverse_fmax_eV_A: float = 0.003
    # Sign resolution for bracketing, not a separate stationarity tolerance.
    reaction_force_tol: float = 0.003
    constraint_tol: float = 1e-8
    min_bracket_width: float = 1e-7
    max_iterations: int = 40
    max_relax_steps: int = 200
    max_step_A: float = 0.05
    branch_rmsd_A: float = 0.05
    branch_step_A: float = 0.3
    # Measured 1-sigma jitter of the reaction force of a converged constrained minimum on
    # one branch (same q, same frame, independent re-relaxations), per system and model.
    # Two continuations then differ with sigma*sqrt(2); the branch-force gate sits at
    # three times that. None keeps the earlier 2 * reaction_force_tol rule, which on
    # fp32 POLAR-1 b04 lies below the measured same-branch spread
    # (benchmarks/explore_d2/results/2026-10-25/hf_ef_rf_jitter.json: sigma 0.0034 eV/A).
    reaction_force_sigma: float | None = None
    # Transverse force threshold for the two bisection-midpoint continuations only.
    # A relaxation that stops as soon as the projected force meets a loose threshold
    # keeps a start-dependent residual along soft transverse directions; near the root
    # that residual dominated the paired reaction-force difference (b05 q 1.4343: 0.027
    # at 3e-3, 0.0018 after both were tightened to 1e-3; DONE_LOG 2026-10-25 J1).
    # None uses transverse_fmax_eV_A. sigma must be measured at this threshold.
    midpoint_transverse_fmax_eV_A: float | None = None
    # Carry only same-branch curvature, never a cross-constraint BFGS secant pair.
    # Disable both options for the cold-start relaxation baseline.
    reuse_curvature: bool = True
    predict_tangent: bool = True

    def branch_force_limit(self):
        if self.reaction_force_sigma is None:
            return 2 * self.reaction_force_tol
        return 3 * np.sqrt(2) * self.reaction_force_sigma

    def __post_init__(self):
        for name, value in asdict(self).items():
            if name in ("reuse_curvature", "predict_tangent"):
                if not isinstance(value, bool):
                    raise ValueError(f"{name} must be a bool")
                continue
            if value is None and name in (
                "reaction_force_sigma",
                "midpoint_transverse_fmax_eV_A",
            ):
                continue
            if not np.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be finite and positive")
        for name in ("max_iterations", "max_relax_steps"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
                raise ValueError(f"{name} must be an integer")


@dataclass
class EFResult:
    status: str = "unresolved:not_started"
    positions: np.ndarray | None = None
    energy_eV: float | None = None
    forces_pes: np.ndarray | None = None
    direction_cartesian: np.ndarray | None = None
    rejected_positions: np.ndarray | None = None
    metrics: dict = field(default_factory=dict)
    history: list = field(default_factory=list)
    bracket: list = field(default_factory=list)
    provenance: dict = field(default_factory=dict)
    spent: int = 0
    wall_s: float = 0.0
    property_status: str = "not_run"
    connection_status: str = "not_run"

    def to_dict(self):
        """JSON-ready search record, also retained when validation cannot run."""
        return _jsonable(asdict(self))


class _Stop(RuntimeError):
    pass


class _Constraint(FixBondCombo):
    def gradient(self, positions):
        if not np.isfinite(positions).all() or any(
            np.linalg.norm(positions[j] - positions[i]) < 1e-10 for i, j in self.intent.pairs
        ):
            raise _Stop("degenerate_coordinate")
        g = self.intent.gradient(positions)
        if not np.isfinite(g).all() or np.vdot(g, g) < 1e-20:
            raise _Stop("degenerate_coordinate")
        return g

    def adjust_positions(self, atoms, new):
        for _ in range(self.max_iter):
            g = self.gradient(new)
            err = self.intent.value(new) - self.target
            if abs(err) <= self.tol:
                return
            new -= err * g / np.vdot(g, g)
        raise _Stop("constraint_not_converged")

    def adjust_forces(self, atoms, forces):
        g = self.gradient(atoms.positions)
        forces -= np.vdot(forces, g) / np.vdot(g, g) * g


def _align(reference, moving):
    a, b = reference - reference.mean(axis=0), moving - moving.mean(axis=0)
    u, _, vt = np.linalg.svd(b.T @ a)
    rot = u @ np.diag([1.0, 1.0, np.linalg.det(u @ vt)]) @ vt
    return b @ rot + reference.mean(axis=0)


def _rmsd(a, b):
    return float(np.sqrt(np.mean(np.sum((a - _align(a, b)) ** 2, axis=1))))


@dataclass
class _RelaxState:
    """Runtime hints for one branch; copy before a trial and commit only on acceptance.

    These are E/F-learned BFGS curvature and the last two accepted geometries, not
    a model Hessian or a force-convergence certificate. No optimizer secant history
    survives a change of q. States never appear in serialized search records.
    """

    hessian: np.ndarray | None = None
    points: tuple = ()  # (q_actual, positions), in one unchanged Cartesian frame

    def copy(self):
        return _RelaxState(
            None if self.hessian is None else self.hessian.copy(),
            tuple((q, x.copy()) for q, x in self.points),
        )


def _relax(
    atoms,
    intent,
    positions,
    q,
    origin,
    factory,
    config,
    cfg,
    history,
    ceiling_eV,
    fmax=None,
    *,
    state=None,
):
    """Constrained BFGS relaxation at q; every evaluated frame is appended to ``history``.

    Returns the first frame meeting ``fmax`` (or cfg.transverse_fmax_eV_A); raises
    _Stop on the energy ceiling, a constraint residual, or the step limit. Optional
    ``state`` is updated only on convergence; callers commit a trial copy only after
    branch checks. A changed input frame discards the hints.
    """
    from prrs.reliability import GuardedCalculator

    work = atoms.copy()
    work.set_constraint()
    work.positions = np.array(positions, dtype=float, copy=True)
    compatible = bool(
        state is not None
        and state.points
        and np.array_equal(state.points[-1][1], work.positions)
    )
    predicted_step = 0.0
    if compatible and cfg.predict_tangent and len(state.points) == 2:
        q_prev, previous = state.points[0]
        q_last, current = state.points[1]
        dq = q_last - q_prev
        if abs(dq) > cfg.constraint_tol and abs(q - q_last) > cfg.constraint_tol:
            # Predict from this branch alone, removing global motion from the
            # secant. Cap the largest atomic displacement by the BFGS step limit.
            delta = (q - q_last) / dq * (current - _align(current, previous))
            size = float(np.linalg.norm(delta, axis=1).max())
            if size > 0 and np.isfinite(delta).all():
                delta *= min(1.0, cfg.max_step_A / size)
                work.positions += delta
                predicted_step = float(np.linalg.norm(delta, axis=1).max())
    constraint = _Constraint(intent, q, tol=cfg.constraint_tol)
    work.set_constraint(constraint)
    work.set_positions(work.positions)
    work.calc = GuardedCalculator(factory(), config)
    opt = BFGS(work, logfile=None, maxstep=cfg.max_step_A)
    reused = compatible and cfg.reuse_curvature and state.hessian is not None
    if reused:
        normal = constraint.gradient(work.positions).ravel()
        normal /= np.linalg.norm(normal)
        nn = np.outer(normal, normal)
        tangent = np.eye(len(normal)) - nn
        opt.H0 = tangent @ state.hessian @ tangent + opt.alpha * nn
    target = cfg.transverse_fmax_eV_A if fmax is None else fmax
    last = None
    for _ in opt.irun(fmax=target, steps=cfg.max_relax_steps):
        energy = float(work.get_potential_energy())
        forces = work.get_forces(apply_constraint=False)
        g = constraint.gradient(work.positions)
        slope = -float(np.vdot(forces, g) / np.vdot(g, g))
        transverse = forces + slope * g
        last = {
            "q_target": float(q),
            "q_actual": intent.value(work.positions),
            "energy_eV": energy,
            "reaction_force": slope,
            "fmax_eV_A": float(np.linalg.norm(forces, axis=1).max()),
            "fmax_projected": float(np.linalg.norm(transverse, axis=1).max()),
            "constraint_residual": abs(intent.value(work.positions) - q),
            "positions": work.positions.copy(),
            "forces_pes": forces.copy(),
            "origin": origin,
            "relax_step": opt.nsteps,
            "curvature_reused": bool(reused),
            "predictor_step_A": predicted_step,
        }
        history.append(last)
        if ceiling_eV is not None and energy > ceiling_eV:
            raise _Stop("above_ceiling")
        if last["constraint_residual"] > cfg.constraint_tol:
            raise _Stop("constraint_not_converged")
        if last["fmax_projected"] <= target:
            if state is not None:
                learned = getattr(opt, "H", None)
                # Zero-step points have no H; retain an inherited H0 if present.
                learned = opt.H0 if learned is None and reused else learned
                state.hessian = None if learned is None else learned.copy()
                points = state.points if compatible else ()
                if points and abs(points[-1][0] - last["q_actual"]) <= cfg.constraint_tol:
                    points = points[:-1]
                state.points = (points + ((last["q_actual"], work.positions.copy()),))[-2:]
            return last
    raise _Stop("transverse_not_converged")


def refine_ts_ef(
    atoms,
    intent,
    bracket,
    factory,
    config,
    ledger,
    task,
    *,
    settings=None,
    ceiling_eV=None,
    attempt=0,
    _warm_starts=None,
):
    """Locate a candidate from two ScanPoints bracketing a local maximum in q.

    Endpoints are re-relaxed/re-evaluated: saved slopes and ``converged`` flags are
    never trusted. Each midpoint is reached separately from BOTH endpoints, with
    no geometry interpolation. Disagreeing branches stop unresolved. After branch
    agreement, either continuation can supply a candidate meeting the full-force,
    transverse-force and constraint tolerances. Only E/F is read; ASE BFGS learns
    a secant approximation from the relaxation trajectory.

    ``factory`` must be counted against ``ledger`` as for constrained_scan. Paid
    failures remain charged, with the last valid bracket and all evaluated frames.
    A converged candidate is not an index-one or C2 claim. This primitive does not
    implement D2 scheduling, crash recovery, or transverse stability certification.
    """
    from prrs.reliability import GateRejected

    cfg = settings or EFSettings()
    if np.any(atoms.pbc):
        raise ValueError("E/F branch locator currently supports nonperiodic molecules")
    if len(bracket) != 2:
        raise ValueError("exactly two bracket endpoints required")
    if ceiling_eV is not None and not np.isfinite(ceiling_eV):
        raise ValueError("ceiling_eV must be finite")
    ordered = sorted(
        zip(bracket, _warm_starts or (_RelaxState(), _RelaxState())),
        key=lambda item: item[0].q_target,
    )
    endpoints = [p for p, _ in ordered]
    left_state, right_state = [state.copy() for _, state in ordered]
    if not all(np.isfinite(p.q_target) for p in endpoints):
        raise ValueError("bracket q must be finite")
    if endpoints[0].q_target >= endpoints[1].q_target:
        raise ValueError("bracket must have positive width")
    for p in endpoints:
        if np.shape(p.positions) != atoms.positions.shape:
            raise ValueError("bracket atom shape mismatch")
        if not np.isfinite(p.positions).all():
            raise ValueError("bracket positions must be finite")
    result = EFResult(
        provenance={
            "algorithm": "HF-EF-bracket-v3",
            "candidate_rule": "branch_agreement_then_any_full_force_converged",
            "intent": asdict(intent),
            "settings": asdict(cfg),
            "branch_force_limit": float(cfg.branch_force_limit()),
            "numbers": atoms.numbers.tolist(),
            "atom_order": "input",
            "charge": config.total_charge,
            "multiplicity": config.multiplicity,
            "ceiling_eV": ceiling_eV,
            "input_bracket": [
                {"q_target": float(p.q_target), "positions": np.array(p.positions).copy()}
                for p in endpoints
            ],
        }
    )
    start, started = ledger.spent, time.monotonic()
    initial_tangent = (
        _align(endpoints[0].positions, endpoints[1].positions) - endpoints[0].positions
    )

    def relax(positions, q, origin, fmax=None, *, state=None):
        return _relax(
            atoms,
            intent,
            positions,
            q,
            origin,
            factory,
            config,
            cfg,
            result.history,
            ceiling_eV,
            fmax,
            state=state,
        )

    def stationary(p):
        # The unbiased Cartesian force includes the reaction-coordinate component.
        # A second |dE/dq| cutoff depends on q's normalization and can sit below
        # measured force jitter even when this structure meets the force contract.
        return (
            p["fmax_eV_A"] <= cfg.fmax_eV_A
            and p["fmax_projected"] <= cfg.transverse_fmax_eV_A
            and p["constraint_residual"] <= cfg.constraint_tol
        )

    def candidate(p):
        result.status = "converged_candidate"
        result.positions = p["positions"].copy()
        result.energy_eV = p["energy_eV"]
        result.forces_pes = p["forces_pes"].copy()
        result.metrics = {k: v for k, v in p.items() if k not in ("positions", "forces_pes")}
        left, right = result.bracket
        tangent = _align(p["positions"], right["positions"]) - _align(
            p["positions"], left["positions"]
        )
        norm = np.linalg.norm(tangent)
        if norm > 1e-12:
            result.direction_cartesian = tangent / norm
            # Compare in the original endpoint frame, removing only global motion.
            current = _align(endpoints[0].positions, right["positions"]) - _align(
                endpoints[0].positions, left["positions"]
            )
            denom = np.linalg.norm(initial_tangent) * np.linalg.norm(current)
            if denom > 1e-20:
                overlap = np.clip(abs(np.vdot(initial_tangent, current)) / denom, 0, 1)
                result.metrics["direction_drift_deg"] = float(np.degrees(np.arccos(overlap)))

    try:
        with ledger.segment(task=task, stage="ts_search_ef", attempt=attempt) as seg:
            try:
                left = relax(
                    endpoints[0].positions, endpoints[0].q_target, "left", state=left_state
                )
                result.bracket = [left]
                right = relax(
                    endpoints[1].positions, endpoints[1].q_target, "right", state=right_state
                )
                result.bracket = [left, right]
                # Strict sign and resolution requirements exclude minima, flat/noisy
                # profiles, and stationary endpoints without a resolved maximum bracket.
                if not (
                    left["reaction_force"] > cfg.reaction_force_tol
                    and right["reaction_force"] < -cfg.reaction_force_tol
                ):
                    raise _Stop("no_resolved_maximum_bracket")
                for iteration in range(cfg.max_iterations):
                    q = (left["q_actual"] + right["q_actual"]) / 2
                    mid = cfg.midpoint_transverse_fmax_eV_A
                    a_state, b_state = left_state.copy(), right_state.copy()
                    a = relax(left["positions"], q, "from_left", mid, state=a_state)
                    b = relax(right["positions"], q, "from_right", mid, state=b_state)
                    agreement = _rmsd(a["positions"], b["positions"])
                    travel = max(
                        _rmsd(a["positions"], left["positions"]),
                        _rmsd(b["positions"], right["positions"]),
                    )
                    a["branch_rmsd_A"] = b["branch_rmsd_A"] = agreement
                    a["branch_step_A"] = b["branch_step_A"] = travel
                    if agreement > cfg.branch_rmsd_A or travel > cfg.branch_step_A:
                        raise _Stop("branch_unresolved")
                    # Two continuations must agree also at the requested force
                    # resolution; geometry alone can hide a steep branch mismatch.
                    if (
                        abs(a["reaction_force"] - b["reaction_force"])
                        > cfg.branch_force_limit()
                    ):
                        raise _Stop("branch_force_disagreement")
                    # Agreement is a paired branch diagnostic; stationarity is a
                    # property of the individual returned structure. Never average
                    # forces/geometries or relax the full-force limit using sigma.
                    eligible = [p for p in (a, b) if stationary(p)]
                    if eligible:
                        candidate(min(eligible, key=lambda p: p["fmax_eV_A"]))
                        break
                    if right["q_actual"] - left["q_actual"] <= cfg.min_bracket_width:
                        raise _Stop("force_resolution_unresolved")
                    # Near-zero unresolved slopes must never replace a signed end.
                    if (
                        a["reaction_force"] > cfg.reaction_force_tol
                        and b["reaction_force"] > cfg.reaction_force_tol
                    ):
                        left = a
                        left_state = a_state
                    elif (
                        a["reaction_force"] < -cfg.reaction_force_tol
                        and b["reaction_force"] < -cfg.reaction_force_tol
                    ):
                        right = a
                        right_state = a_state
                    else:
                        raise _Stop("force_resolution_unresolved")
                    result.bracket = [left, right]
                else:
                    raise _Stop("iteration_limit")
            except _Stop as exc:
                result.status = f"unresolved:{exc}"
            except BudgetExhausted:
                result.status = "budget_exhausted"
            except GateRejected as exc:
                result.status = f"gate_rejected:{exc.code}"
                result.rejected_positions = exc.atoms.positions.copy()
            finally:
                seg.outcome = result.status
    finally:
        result.spent = ledger.spent - start
        result.wall_s = time.monotonic() - started
    return result
