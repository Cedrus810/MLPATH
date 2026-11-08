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
from ase.optimize import BFGS, FIRE

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
    # Two-stage midpoints (HF-EF-bracket-v4; user 2026-11-02,
    # docs/PLAN_HF_MIDPOINT_COARSE_2026-11-02.md). Both continuations first stop at the
    # coarse threshold. The bracket shrinks on that stage only when every paired check
    # passes, neither side is a candidate, and both slopes share a sign beyond the sign
    # bound; otherwise both continue to the midpoint threshold as in v3. On b05, 65 of
    # 71 midpoints that never met 1e-3 had one resolved sign at every frame <= 2e-3.
    # Both None keeps the single-stage rule.
    midpoint_coarse_transverse_fmax_eV_A: float | None = None
    midpoint_sign_bound_eV_A: float | None = None
    # Optimizer of the fine stage of two-stage midpoints: "bfgs" (v4) or "fire"
    # (HF-EF-bracket-v5; user 2026-11-04, docs/PLAN_HF_MIDPOINT_FINE_FIRE_2026-11-04.md).
    # In a cold same-frame comparison FIRE finished 33 of 38 near-root fine stages that
    # BFGS left unconverged (BFGS 3 of 38). FIRE starts from zero velocity on every call
    # and learns no curvature; the coarse stage and all other relaxations stay BFGS.
    midpoint_fine_optimizer: str = "bfgs"
    # Coarse-stage candidate (HF-EF-bracket-v6; user 2026-11-06,
    # docs/PLAN_HF_MIDPOINT_COARSE_CANDIDATE_2026-11-06.md). When the coarse pair passes
    # the geometry and paired reaction-force checks and either side already meets the
    # candidate contract, that side is the candidate; the fine stage, which only
    # sharpens slopes for bracketing, is skipped. In the 11-04 v5 acceptance 14% of the
    # failed-group HS forwards were fine stages run after such a coarse pair. v5 only.
    midpoint_coarse_candidate: bool = False
    # Endpoint candidate (HF-EF-bracket-v7; user 2026-11-07,
    # docs/PLAN_HF_ENDPOINT_CANDIDATE_2026-11-07.md). When the re-relaxed ends do not
    # resolve a maximum bracket, an end that meets the candidate contract while the other
    # end's slope points to a maximum on its side is the candidate (on b05 v4 sB all three
    # no_resolved_maximum_bracket stops had such an end). v6 only.
    endpoint_candidate: bool = False
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
            if name in (
                "reuse_curvature",
                "predict_tangent",
                "midpoint_coarse_candidate",
                "endpoint_candidate",
            ):
                if not isinstance(value, bool):
                    raise ValueError(f"{name} must be a bool")
                continue
            if name == "midpoint_fine_optimizer":
                if value not in ("bfgs", "fire"):
                    raise ValueError(f"{name} must be 'bfgs' or 'fire'")
                continue
            if value is None and name in (
                "reaction_force_sigma",
                "midpoint_transverse_fmax_eV_A",
                "midpoint_coarse_transverse_fmax_eV_A",
                "midpoint_sign_bound_eV_A",
            ):
                continue
            if not np.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be finite and positive")
        coarse, bound = self.midpoint_coarse_transverse_fmax_eV_A, self.midpoint_sign_bound_eV_A
        if (coarse is None) != (bound is None):
            raise ValueError("coarse midpoint threshold and sign bound go together")
        if coarse is not None:
            fine = self.midpoint_transverse_fmax_eV_A or self.transverse_fmax_eV_A
            if coarse <= fine:
                raise ValueError("coarse midpoint threshold must exceed the midpoint one")
            if bound < self.reaction_force_tol:
                raise ValueError("sign bound must not be below reaction_force_tol")
        elif self.midpoint_fine_optimizer != "bfgs":
            raise ValueError("a fine-stage optimizer needs the two-stage midpoint rule")
        if self.midpoint_coarse_candidate and self.midpoint_fine_optimizer != "fire":
            raise ValueError("the coarse-stage candidate rule needs v5 (two-stage, FIRE fine)")
        if self.endpoint_candidate and not self.midpoint_coarse_candidate:
            raise ValueError("the endpoint candidate rule needs v6 (coarse-stage candidate)")
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
    optimizer="bfgs",
):
    """Constrained BFGS (or FIRE) relaxation at q; every evaluated frame is appended
    to ``history``.

    Returns the first frame meeting ``fmax`` (or cfg.transverse_fmax_eV_A); raises
    _Stop on the energy ceiling, a constraint residual, or the step limit. Optional
    ``state`` is updated only on convergence; callers commit a trial copy only after
    branch checks. A changed input frame discards the hints. FIRE starts from zero
    velocity, uses no curvature hint, and on convergence clears the state's BFGS
    curvature, which no longer describes the accepted point; the geometry is kept.
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
    if optimizer == "fire":
        opt = FIRE(work, logfile=None, maxstep=cfg.max_step_A)
        reused = False
    else:
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
            "optimizer": optimizer,
        }
        history.append(last)
        if ceiling_eV is not None and energy > ceiling_eV:
            raise _Stop("above_ceiling")
        if last["constraint_residual"] > cfg.constraint_tol:
            raise _Stop("constraint_not_converged")
        if last["fmax_projected"] <= target:
            if state is not None:
                if optimizer == "fire":
                    state.hessian = None
                else:
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
            "algorithm": "HF-EF-bracket-v3"
            if cfg.midpoint_coarse_transverse_fmax_eV_A is None
            else "HF-EF-bracket-v7"
            if cfg.endpoint_candidate
            else "HF-EF-bracket-v6"
            if cfg.midpoint_coarse_candidate
            else "HF-EF-bracket-v5"
            if cfg.midpoint_fine_optimizer == "fire"
            else "HF-EF-bracket-v4",
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

    def relax(positions, q, origin, fmax=None, *, state=None, optimizer="bfgs"):
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
            optimizer=optimizer,
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

    def coarse_sign(a, b, left, right):
        """+1 / -1 shrinks the bracket at the coarse stage, 0 continues both sides to
        the midpoint threshold. The v3 checks all apply; failing one is not a stop here.
        A side that is already a candidate defers to the v3 candidate rule (v4 / v5), or,
        with midpoint_coarse_candidate (v6), marks the pair ``coarse_candidate`` once the
        geometry and paired-force checks have passed; the caller then accepts it."""
        if right["q_actual"] - left["q_actual"] <= cfg.min_bracket_width:
            return 0
        agreement = _rmsd(a["positions"], b["positions"])
        travel = max(
            _rmsd(a["positions"], left["positions"]),
            _rmsd(b["positions"], right["positions"]),
        )
        a["coarse_branch_rmsd_A"] = b["coarse_branch_rmsd_A"] = agreement
        a["coarse_branch_step_A"] = b["coarse_branch_step_A"] = travel
        if agreement > cfg.branch_rmsd_A or travel > cfg.branch_step_A:
            return 0
        if abs(a["reaction_force"] - b["reaction_force"]) > cfg.branch_force_limit():
            return 0
        if stationary(a) or stationary(b):
            if cfg.midpoint_coarse_candidate:
                a["coarse_candidate"] = b["coarse_candidate"] = True
            return 0
        bound = cfg.midpoint_sign_bound_eV_A
        if a["reaction_force"] > bound and b["reaction_force"] > bound:
            return 1
        if a["reaction_force"] < -bound and b["reaction_force"] < -bound:
            return -1
        return 0

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
                    # v7: an end that already meets the candidate contract, with the
                    # other end's slope pointing to a maximum on its side. Never a
                    # minimum along q: that end would have the opposite neighbour slope.
                    ends = []
                    if cfg.endpoint_candidate:
                        if (
                            stationary(right)
                            and left["reaction_force"] > cfg.reaction_force_tol
                        ):
                            ends.append(right)
                        if (
                            stationary(left)
                            and right["reaction_force"] < -cfg.reaction_force_tol
                        ):
                            ends.append(left)
                    if not ends:
                        raise _Stop("no_resolved_maximum_bracket")
                    best = min(ends, key=lambda p: p["fmax_eV_A"])
                    best["endpoint_candidate"] = True
                    candidate(best)
                    return result
                for iteration in range(cfg.max_iterations):
                    q = (left["q_actual"] + right["q_actual"]) / 2
                    mid = cfg.midpoint_transverse_fmax_eV_A
                    a_state, b_state = left_state.copy(), right_state.copy()
                    a_start, b_start = left["positions"], right["positions"]
                    if cfg.midpoint_coarse_transverse_fmax_eV_A is not None:
                        coarse = cfg.midpoint_coarse_transverse_fmax_eV_A
                        a = relax(a_start, q, "from_left_coarse", coarse, state=a_state)
                        b = relax(b_start, q, "from_right_coarse", coarse, state=b_state)
                        sign = coarse_sign(a, b, left, right)
                        a["coarse_sign"] = b["coarse_sign"] = sign
                        if a.get("coarse_candidate"):
                            eligible = [p for p in (a, b) if stationary(p)]
                            candidate(min(eligible, key=lambda p: p["fmax_eV_A"]))
                            break
                        if sign:
                            if sign > 0:
                                left, left_state = a, a_state
                            else:
                                right, right_state = a, a_state
                            result.bracket = [left, right]
                            continue
                        # Not decided: both continue from where they stopped.
                        a_start, b_start = a["positions"], b["positions"]
                    fine = cfg.midpoint_fine_optimizer
                    a = relax(a_start, q, "from_left", mid, state=a_state, optimizer=fine)
                    b = relax(b_start, q, "from_right", mid, state=b_state, optimizer=fine)
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
