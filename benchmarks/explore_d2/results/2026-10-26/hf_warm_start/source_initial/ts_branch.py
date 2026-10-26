"""Two-sided branch continuation for brackets whose ends lie on different branches.

docs/PLAN_HF_JUMP_2026-10-25.md J1; TS_HESSIAN_FREE_DESIGN_RECORD_2026-10-24 §5.2: shrink
the step and continue each branch separately; if no continuous locating interval forms,
keep both sides and return an explicit branch-unresolved result. Geometries of the two
branches are never averaged or interpolated, and an energy drop is never read as a TS.

Each side starts from its own end of the bracket and walks toward the other end by
constrained relaxation (the ``ts_ef`` relaxation, E/F only). A step is "on the branch"
when its aligned displacement stays within ``jump_factor`` times the larger of the
recent displacement rate times the step and the measured same-branch geometry noise.
A smooth branch moves proportionally to the step, a jump does not; so a suspected
jump is retried with the step halved, and only a jump that survives down to
``min_step_q`` terminates the branch. A branch on which the reaction force changes from
uphill to downhill before it terminates holds a continuous maximum; that same-branch
pair is handed to ``refine_ts_ef``. Otherwise the result is ``unresolved:hysteresis``
with both termination points and the geometries on either side of each jump.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from types import SimpleNamespace
import time

import numpy as np

from .ledger import BudgetExhausted
from .poses import _jsonable
from .ts_ef import EFSettings, _RelaxState, _relax, _rmsd, _Stop, refine_ts_ef


@dataclass(frozen=True)
class BranchSettings:
    # Measured same-branch geometry noise (max pairwise aligned RMSD of independent
    # re-relaxations at one q), per system and model; no default on purpose.
    geometry_noise_A: float
    step_q: float = 0.01
    min_step_q: float = 1e-3
    jump_factor: float = 5.0
    window: int = 3
    max_points: int = 400

    def __post_init__(self):
        for name, value in asdict(self).items():
            if not np.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be finite and positive")
        if self.min_step_q >= self.step_q:
            raise ValueError("min_step_q must be below step_q")
        if self.jump_factor <= 1:
            raise ValueError("jump_factor must exceed 1")
        for name in ("window", "max_points"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
                raise ValueError(f"{name} must be an integer")


@dataclass
class BranchResult:
    # status: final outcome (the hand-off's status when a continuous maximum was found);
    # branch_status: what the continuation itself established.
    status: str = "unresolved:not_started"
    branch_status: str = "unresolved:not_started"
    sides: list = field(default_factory=list)
    refine: dict | None = None
    history: list = field(default_factory=list)
    provenance: dict = field(default_factory=dict)
    spent: int = 0
    wall_s: float = 0.0

    def to_dict(self):
        return _jsonable(asdict(self))


def _brief(p):
    return {k: v for k, v in p.items() if k != "forces_pes"}


def _walk(sign, start, q_stop, relax, cfg, bcfg, state):
    """Continue one branch from ``start`` toward q_stop (sign +1: increasing q)."""
    tol = cfg.reaction_force_tol
    side = {"origin": start["origin"], "points": [_brief(start)], "status": "reached_end"}
    rates = []
    current, dq = start, bcfg.step_q
    for _ in range(bcfg.max_points):
        remaining = sign * (q_stop - current["q_actual"])
        if remaining <= 1e-12:
            return side, None
        step = min(dq, remaining)
        trial_state = state.copy()
        trial = relax(
            current["positions"],
            current["q_actual"] + sign * step,
            side["origin"],
            state=trial_state,
        )
        moved = _rmsd(trial["positions"], current["positions"])
        expected = np.median(rates[-bcfg.window :]) * step if rates else 0.0
        limit = bcfg.jump_factor * max(expected, bcfg.geometry_noise_A)
        if moved > limit:
            if step <= bcfg.min_step_q:
                side["status"] = "terminated"
                side["terminated_at_q"] = current["q_actual"]
                side["jump"] = {
                    "q": trial["q_actual"],
                    "step_q": step,
                    "rmsd_A": moved,
                    "limit_A": limit,
                    "dE_eV": trial["energy_eV"] - current["energy_eV"],
                    "reaction_force_after": trial["reaction_force"],
                    "positions_after": trial["positions"],
                }
                return side, None
            dq = step / 2
            continue
        rates.append(moved / step)
        side["points"].append(_brief(trial))
        before, after = current["reaction_force"], trial["reaction_force"]
        before_state = state
        state = trial_state
        current, dq = trial, min(bcfg.step_q, 2 * step)
        # Uphill-to-downhill along increasing q is a maximum; walking down in q the
        # same maximum appears as downhill-to-uphill. Unresolved slopes never count.
        if sign > 0 and before > tol and after < -tol:
            side["status"] = "maximum_bracketed"
            return side, ((side["points"][-2], side["points"][-1]), (before_state, state))
        if sign < 0 and before < -tol and after > tol:
            side["status"] = "maximum_bracketed"
            return side, ((side["points"][-1], side["points"][-2]), (state, before_state))
    side["status"] = "point_limit"
    return side, None


def continue_branches_ef(
    atoms,
    intent,
    bracket,
    factory,
    config,
    ledger,
    task,
    *,
    branch_settings,
    settings=None,
    ceiling_eV=None,
    attempt=0,
):
    """Locate a continuous maximum on either branch of a bracket, or report hysteresis.

    ``bracket`` holds two ScanPoints (any order). Paid evaluations stay charged on every
    exit. A ``converged_candidate`` from the hand-off is a search result only."""
    from prrs.reliability import GateRejected

    cfg = settings or EFSettings()
    bcfg = branch_settings
    if np.any(atoms.pbc):
        raise ValueError("branch continuation currently supports nonperiodic molecules")
    if len(bracket) != 2:
        raise ValueError("exactly two bracket endpoints required")
    low, high = sorted(bracket, key=lambda p: p.q_target)
    if not low.q_target < high.q_target:
        raise ValueError("bracket must have positive width")
    result = BranchResult(
        provenance={
            "algorithm": "HF-EF-branches-v2",
            "intent": asdict(intent),
            "settings": asdict(cfg),
            "branch_settings": asdict(bcfg),
            "numbers": atoms.numbers.tolist(),
            "ceiling_eV": ceiling_eV,
            "input_bracket": [
                {"q_target": float(p.q_target), "positions": np.array(p.positions).copy()}
                for p in (low, high)
            ],
        }
    )
    start, started = ledger.spent, time.monotonic()

    def relax(positions, q, origin, *, state):
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
            state=state,
        )

    found = None
    try:
        with ledger.segment(task=task, stage="ts_branches_ef", attempt=attempt) as seg:
            try:
                low_state, high_state = _RelaxState(), _RelaxState()
                ends = (
                    relax(low.positions, low.q_target, "low", state=low_state),
                    relax(high.positions, high.q_target, "high", state=high_state),
                )
                for sign, begin, stop, state in (
                    (1, ends[0], ends[1], low_state),
                    (-1, ends[1], ends[0], high_state),
                ):
                    side, found = _walk(sign, begin, stop["q_actual"], relax, cfg, bcfg, state)
                    result.sides.append(side)
                    if found is not None:
                        result.status = "continuous_maximum"
                        break
                else:
                    result.status = "unresolved:hysteresis"
            except _Stop as exc:
                result.status = f"unresolved:{exc}"
            except BudgetExhausted:
                result.status = "budget_exhausted"
            except GateRejected as exc:
                result.status = f"gate_rejected:{exc.code}"
            finally:
                seg.outcome = result.status
                result.branch_status = result.status
        if found is not None:
            points, states = found
            pair = [
                SimpleNamespace(q_target=p["q_actual"], positions=p["positions"])
                for p in points
            ]
            refined = refine_ts_ef(
                atoms,
                intent,
                pair,
                factory,
                config,
                ledger,
                task,
                settings=cfg,
                ceiling_eV=ceiling_eV,
                attempt=attempt,
                _warm_starts=states,
            )
            result.refine = refined.to_dict()
            result.status = refined.status
    finally:
        result.spent = ledger.spent - start
        result.wall_s = time.monotonic() - started
    return result
