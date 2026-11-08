"""Frozen copy of the HF-EF-branches-v3 branch walk (2026-10-28), for same-machine policy
comparisons only (hs_policy_compare.py). v3 halved the step after a non-converging trial
relaxation and ended the side only below min_step_q; v4 (explore.ts_branch, 2026-10-29)
ends the side at the first such failure. Everything else is the v4 code verbatim. Not
used by any production path."""

import numpy as np

from explore.ts_branch import _brief, _unconverged
from explore.ts_ef import _rmsd, _Stop


def _walk(sign, start, q_stop, relax, cfg, bcfg, state, history=None):
    """Continue one branch from ``start`` toward q_stop (sign +1: increasing q).
    ``history`` (the shared frame list ``relax`` appends to) only feeds diagnostics."""
    tol = cfg.reaction_force_tol
    side = {"origin": start["origin"], "points": [_brief(start)], "status": "reached_end"}
    rates = []
    current, dq = start, bcfg.step_q
    for _ in range(bcfg.max_points):
        remaining = sign * (q_stop - current["q_actual"])
        # The constraint holds q only to constraint_tol: a smaller "remaining" can never
        # be closed, and walking it re-relaxed the end point until max_points (b04 HF
        # smoke 2026-10-26: 378 repeated points at q 4.0895).
        if remaining <= cfg.constraint_tol:
            return side, None
        step = min(dq, remaining)
        trial_state = state.copy()
        target = current["q_actual"] + sign * step
        try:
            trial = relax(current["positions"], target, side["origin"], state=trial_state)
        except _Stop as exc:
            if str(exc) != "transverse_not_converged":
                raise
            # current and state stay accepted; the trial frames remain diagnostics only.
            failure = _unconverged(history, target, step, current)
            side.setdefault("unconverged_trials", []).append(
                {k: v for k, v in failure.items() if k != "positions"}
            )
            if step <= bcfg.min_step_q:
                side["status"] = "terminated"
                side["termination"] = "relax_unconverged"
                side["terminated_at_q"] = current["q_actual"]
                side["relax_unconverged"] = failure
                return side, None
            dq = step / 2
            continue
        moved = _rmsd(trial["positions"], current["positions"])
        expected = np.median(rates[-bcfg.window :]) * step if rates else 0.0
        limit = bcfg.jump_factor * max(expected, bcfg.geometry_noise_A)
        if moved > limit:
            if step <= bcfg.min_step_q:
                side["status"] = "terminated"
                side["termination"] = "jump"
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
