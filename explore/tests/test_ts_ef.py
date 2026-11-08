"""Analytic coupled-bond surface: search/force contracts, not molecular acceptance."""

import json
from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes

from explore.ledger import Ledger, counted_factory
from explore.probe import Intent
from explore.ts_ef import EFSettings, refine_ts_ef
from prrs.config import SearchConfig


class CoupledBonds(Calculator):
    implemented_properties = ["energy", "forces"]

    def __init__(self, calls, profile="maximum", branches=False, reaction_offset=0.0):
        super().__init__()
        self.calls, self.profile, self.branches = calls, profile, branches
        self.reaction_offset = reaction_offset

    def get_hessian(self, *args, **kwargs):
        raise AssertionError("Hessian entered E/F search")

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        self.calls.append(atoms.positions.copy())
        p = atoms.positions
        v, w = p[1] - p[0], p[2] - p[1]
        r, s = np.linalg.norm(v), np.linalg.norm(w)
        x = (r - 1.8) / 0.6
        if self.profile in ("maximum", "jitter"):
            energy, dr = 0.4 * (x * x - 1) ** 2, 1.6 * x * (x * x - 1) / 0.6
        elif self.profile == "minimum":
            energy, dr = x * x, 2 * x / 0.6
        else:
            energy, dr = 0.0, 0.0
        dr += self.reaction_offset
        if self.profile == "jitter":
            # Same-branch reaction-force jitter: alternates with every evaluation.
            dr += 2e-3 if len(self.calls) % 2 else -2e-3
        if self.profile == "noisy":
            # A deliberately unresolved, discontinuous normal force floor.
            dr = -0.01 if r >= 1.8 else 0.01
        if self.branches:
            y = s - 2.0
            energy += 10 * (y * y - 0.16) ** 2
            ds = 40 * y * (y * y - 0.16)
        else:
            # The second bond must move with q. Holding the input direction
            # fixed would miss the stationary point of this coupled surface.
            y = s - 2.0 - 0.4 * x
            energy += 3 * y * y
            ds = 6 * y
            dr -= ds * 0.4 / 0.6
        f = np.zeros_like(p)
        f[0] += dr * v / r
        f[1] -= dr * v / r
        f[1] += ds * w / s
        f[2] -= ds * w / s
        if self.profile == "invalid_forces":
            f[:] = np.nan
        self.results = {"energy": energy, "forces": f}


SETTINGS = EFSettings(fmax_eV_A=2e-4, transverse_fmax_eV_A=1e-5, reaction_force_tol=1e-4)
CONFIG = SearchConfig(closed_shell_only=False)
INTENT = Intent("coupled", ((0, 1),), (1.0,))


def frame(r, s=None):
    s = 2.0 + 0.4 * (r - 1.8) / 0.6 if s is None else s
    return Atoms("C3", positions=[[0, 0, 0], [r, 0, 0], [r + s, 0, 0]])


def endpoint(atoms, intent=INTENT):
    return SimpleNamespace(q_target=intent.value(atoms.positions), positions=atoms.positions)


def run(
    tmp_path,
    *,
    frames=None,
    settings=SETTINGS,
    budget=5000,
    intent=INTENT,
    ceiling=None,
    profile="maximum",
    branches=False,
    reaction_offsets=(),
):
    frames = frames or (frame(1.51, 1.9), frame(2.02, 2.1))
    calls = []
    ledger = Ledger(tmp_path / "ledger.jsonl", arm="HF", budget=budget)
    offsets = iter(reaction_offsets)
    factory = counted_factory(
        lambda: CoupledBonds(calls, profile, branches, next(offsets, 0.0)), ledger
    )
    result = refine_ts_ef(
        frames[0],
        intent,
        [endpoint(f, intent) for f in frames],
        factory,
        CONFIG,
        ledger,
        "test",
        settings=settings,
        ceiling_eV=ceiling,
    )
    ledger.close()
    assert result.spent == ledger.spent == len(calls)
    rows = [json.loads(line) for line in ledger.path.read_text().splitlines()]
    assert rows[-1]["outcome"] == result.status
    assert rows[-1]["stage"] == "ts_search_ef"
    assert not [r for r in rows if r.get("task") == "<unattributed>"]
    json.dumps(result.to_dict(), allow_nan=False)
    return result


def test_coupled_stationary_candidate_with_all_full_spectrum_routes_disabled(
    tmp_path, monkeypatch
):
    import prrs.runner as runner

    def forbidden(*args, **kwargs):
        raise AssertionError("full spectrum entered search")

    monkeypatch.setattr(runner, "curvature_spectrum", forbidden)
    monkeypatch.setattr(runner, "_mass_weighted_hessian", forbidden)
    result = run(tmp_path)
    assert result.status == "converged_candidate", result.status
    p = result.positions
    assert abs(np.linalg.norm(p[1] - p[0]) - 1.8) < 1e-4
    assert abs(np.linalg.norm(p[2] - p[1]) - 2.0) < 1e-4
    assert np.linalg.norm(result.forces_pes, axis=1).max() <= SETTINGS.fmax_eV_A
    assert result.property_status == result.connection_status == "not_run"
    assert np.linalg.norm(result.direction_cartesian) == pytest.approx(1)
    assert result.metrics["direction_drift_deg"] >= 0
    assert result.bracket and result.history and result.spent > 4


@pytest.mark.parametrize("reverse", [False, True])
def test_rigid_motion_and_scan_orientation(tmp_path, reverse):
    frames = [frame(1.51), frame(2.02)]
    for f in frames:
        f.rotate(61, (1, 2, 3))
        f.positions += [5, -3, 2]
    if reverse:
        frames.reverse()
    result = run(tmp_path, frames=frames)
    assert result.status == "converged_candidate", result.status
    assert abs(result.energy_eV - 0.4) < 1e-6


def test_negative_coordinate_definition(tmp_path):
    result = run(tmp_path, intent=Intent("negative", ((0, 1),), (-1.0,)))
    assert result.status == "converged_candidate", result.status


@pytest.mark.parametrize("profile", ["minimum", "flat"])
def test_minimum_and_flat_profile_are_not_maximum_brackets(tmp_path, profile):
    result = run(tmp_path, profile=profile)
    assert result.status == "unresolved:no_resolved_maximum_bracket"
    assert result.positions is None and len(result.bracket) == 2


def test_different_same_graph_branches_are_not_interpolated(tmp_path):
    result = run(tmp_path, frames=[frame(1.51, 1.6), frame(2.02, 2.4)], branches=True)
    assert result.status == "unresolved:branch_unresolved", result.status
    assert result.positions is None
    assert result.history[-1]["branch_rmsd_A"] > SETTINGS.branch_rmsd_A
    assert len(result.bracket) == 2


def test_budget_keeps_paid_frames_and_valid_bracket(tmp_path):
    result = run(tmp_path, frames=[frame(1.51), frame(2.02)], budget=3)
    assert result.status == "budget_exhausted"
    assert result.spent == 3 and result.history
    assert len(result.bracket) == 2


def test_reliability_failure_keeps_rejected_structure_and_cost(tmp_path):
    result = run(tmp_path, profile="invalid_forces")
    assert result.status == "gate_rejected:nonfinite"
    assert result.spent == 1
    assert result.rejected_positions is not None
    assert len(result.provenance["input_bracket"]) == 2


def test_ceiling_is_checked_during_relaxation(tmp_path):
    result = run(tmp_path, ceiling=0.01)
    assert result.status == "unresolved:above_ceiling"
    assert result.history and result.spent == 1


def test_unconverged_transverse_relaxation_does_not_supply_a_slope(tmp_path):
    result = run(
        tmp_path,
        frames=[frame(1.51, 2.8), frame(2.02)],
        settings=replace(SETTINGS, max_relax_steps=1),
    )
    assert result.status == "unresolved:transverse_not_converged"
    assert result.positions is None


def test_tiny_interval_alone_is_not_convergence(tmp_path):
    result = run(tmp_path, settings=replace(SETTINGS, min_bracket_width=1.0))
    assert result.status == "unresolved:force_resolution_unresolved"
    assert result.positions is None


def test_force_noise_floor_is_not_relabeled_convergence(tmp_path):
    result = run(tmp_path, profile="noisy")
    assert result.status == "unresolved:force_resolution_unresolved", result.status
    assert result.positions is None


def test_iteration_limit_keeps_a_signed_bracket(tmp_path):
    result = run(tmp_path, settings=replace(SETTINGS, max_iterations=1))
    assert result.status == "unresolved:iteration_limit"
    assert result.bracket[0]["reaction_force"] > 0 > result.bracket[1]["reaction_force"]


def test_degenerate_coordinate_has_explicit_outcome(tmp_path):
    intent = Intent("zero", ((0, 1), (0, 1)), (1.0, -1.0))
    frames = [endpoint(frame(1.51)), endpoint(frame(2.02))]
    ledger = Ledger(tmp_path / "zero.jsonl", arm="HF", budget=10)
    result = refine_ts_ef(frame(1.51), intent, frames, lambda: None, CONFIG, ledger, "zero")
    assert result.status == "unresolved:degenerate_coordinate"
    assert result.spent == 0


@pytest.mark.parametrize(
    "kwargs",
    [
        {"max_iterations": 0},
        {"max_relax_steps": 1.5},
        {"max_relax_steps": 1.0},
        {"fmax_eV_A": float("nan")},
    ],
)
def test_invalid_settings(kwargs):
    with pytest.raises(ValueError):
        EFSettings(**kwargs)


def test_branch_force_gate_scales_with_measured_jitter(tmp_path):
    # Continuations on one branch differ only by jitter of up to 4e-3 eV/A in the slope.
    # Against 2 * reaction_force_tol (2e-4) that is a "disagreement"; against
    # 3 * sqrt(2) * sigma with the jitter's own scale it is not.
    tight = run(tmp_path / "tight", profile="jitter")
    assert tight.status == "unresolved:branch_force_disagreement"
    measured = replace(SETTINGS, reaction_force_sigma=2e-3)
    loose = run(tmp_path / "loose", profile="jitter", settings=measured)
    assert loose.status != "unresolved:branch_force_disagreement"
    assert loose.provenance["branch_force_limit"] == pytest.approx(3 * np.sqrt(2) * 2e-3)


@pytest.mark.parametrize("offsets", [(0.005, 0.012), (0.012, 0.005)])
def test_one_force_converged_continuation_supplies_candidate_after_agreement(tmp_path, offsets):
    # At the symmetric midpoint, one continuation has a force above the full-force
    # tolerance. The other passes it but exceeds the old separate |RF| cutoff.
    settings = EFSettings(reaction_force_sigma=0.0034)
    result = run(
        tmp_path,
        frames=[frame(1.6), frame(2.0)],
        settings=settings,
        reaction_offsets=(0.0, 0.0, *offsets),
    )
    assert result.status == "converged_candidate", result.status
    assert settings.reaction_force_tol < abs(result.metrics["reaction_force"])
    assert result.metrics["fmax_eV_A"] <= settings.fmax_eV_A
    assert result.metrics["fmax_projected"] <= settings.transverse_fmax_eV_A
    assert result.metrics["constraint_residual"] <= settings.constraint_tol
    assert np.linalg.norm(result.forces_pes, axis=1).max() == pytest.approx(
        result.metrics["fmax_eV_A"]
    )
    assert result.metrics["origin"] == ("from_left" if offsets[0] == 0.005 else "from_right")
    last_pair = [p for p in result.history if "branch_rmsd_A" in p][-2:]
    assert max(p["fmax_eV_A"] for p in last_pair) > settings.fmax_eV_A
    assert abs(last_pair[0]["reaction_force"] - last_pair[1]["reaction_force"]) < (
        settings.branch_force_limit()
    )
    assert result.provenance["algorithm"] == "HF-EF-bracket-v3"
    assert (
        result.provenance["candidate_rule"] == "branch_agreement_then_any_full_force_converged"
    )
    assert result.property_status == result.connection_status == "not_run"


def test_force_converged_continuation_cannot_bypass_branch_force_disagreement(tmp_path):
    # Analogous to b04's eighth midpoint: a good force on one side cannot excuse
    # a much larger disagreement with its peer, even when the geometries agree.
    settings = EFSettings(reaction_force_sigma=0.0034)
    result = run(
        tmp_path,
        frames=[frame(1.6), frame(2.0)],
        settings=settings,
        reaction_offsets=(0.0, 0.0, 0.038, -0.0047),
    )
    assert result.status == "unresolved:branch_force_disagreement"
    assert result.positions is None
    assert result.history[-1]["fmax_eV_A"] <= settings.fmax_eV_A


def test_measured_sigma_does_not_loosen_full_force_tolerance(tmp_path):
    settings = EFSettings(reaction_force_sigma=0.0034, max_iterations=1)
    result = run(
        tmp_path,
        frames=[frame(1.6), frame(2.0)],
        settings=settings,
        reaction_offsets=(0.0, 0.0, 0.02, 0.024),
    )
    assert result.status == "unresolved:iteration_limit"
    assert result.positions is None
    assert all(p["fmax_eV_A"] > settings.fmax_eV_A for p in result.history)


def test_midpoint_threshold_applies_to_bisection_midpoints_only(tmp_path):
    tight = replace(SETTINGS, midpoint_transverse_fmax_eV_A=2e-6)
    result = run(tmp_path, settings=tight)
    assert result.status == "converged_candidate"
    finals = [h for h in result.history if "branch_rmsd_A" in h]
    assert finals and all(h["fmax_projected"] <= 2e-6 for h in finals)
    # Endpoints still stop at the general threshold.
    ends = [h for h in result.history if h["origin"] in ("left", "right")]
    assert max(h["fmax_projected"] for h in ends) > 2e-6


# Two-stage midpoints (HF-EF-bracket-v4; docs/PLAN_HF_MIDPOINT_COARSE_2026-11-02.md).
def coarse(bound, base=SETTINGS, **kwargs):
    return replace(
        base,
        midpoint_coarse_transverse_fmax_eV_A=1e-3,
        midpoint_sign_bound_eV_A=bound,
        **kwargs,
    )


def _coarse(history):
    return [h for h in history if h["origin"].endswith("_coarse") and "coarse_sign" in h]


@pytest.mark.parametrize(
    "kwargs",
    [
        # Default thresholds: midpoint falls back to 3e-3, reaction_force_tol is 3e-3.
        {"midpoint_coarse_transverse_fmax_eV_A": 4e-3},
        {"midpoint_sign_bound_eV_A": 5e-3},
        {"midpoint_coarse_transverse_fmax_eV_A": 3e-3, "midpoint_sign_bound_eV_A": 5e-3},
        {"midpoint_coarse_transverse_fmax_eV_A": 4e-3, "midpoint_sign_bound_eV_A": 1e-3},
    ],
)
def test_invalid_coarse_settings(kwargs):
    with pytest.raises(ValueError):
        EFSettings(**kwargs)


def test_unreachable_sign_bound_keeps_the_v3_outcome(tmp_path):
    old = run(tmp_path / "v3")
    new = run(tmp_path / "v4", settings=coarse(1e3))
    assert old.provenance["algorithm"] == "HF-EF-bracket-v3"
    assert new.provenance["algorithm"] == "HF-EF-bracket-v4"
    assert new.status == old.status == "converged_candidate"
    assert new.metrics["q_actual"] == pytest.approx(old.metrics["q_actual"], abs=1e-6)
    stages = _coarse(new.history)
    assert stages and all(h["coarse_sign"] == 0 for h in stages)
    # Every midpoint still ends at the fine threshold, after its coarse stage.
    fine = [h for h in new.history if "branch_rmsd_A" in h]
    assert len(fine) == len(stages)
    assert all(h["fmax_projected"] <= SETTINGS.transverse_fmax_eV_A for h in fine)


def test_resolved_sign_shrinks_on_the_coarse_stage(tmp_path):
    old = run(tmp_path / "v3")
    new = run(tmp_path / "v4", settings=coarse(1e-3))
    assert new.status == "converged_candidate"
    assert new.metrics["q_actual"] == pytest.approx(old.metrics["q_actual"], abs=1e-6)
    decided = [h for h in _coarse(new.history) if h["coarse_sign"] != 0]
    assert decided and all(abs(h["reaction_force"]) > 1e-3 for h in decided)
    assert all(h["fmax_projected"] <= 1e-3 for h in decided)
    # Decided midpoints skip the fine stage. Cost is not asserted here: on this
    # quadratic surface every relaxation takes a few steps, so the coarse stage's own
    # evaluations can outweigh what it saves; cost is measured on the model.
    fine_new = [h for h in new.history if "branch_rmsd_A" in h]
    fine_old = [h for h in old.history if "branch_rmsd_A" in h]
    assert len(fine_new) < len(fine_old)


def test_coarse_force_disagreement_continues_instead_of_stopping(tmp_path):
    # Calculator order: left end, right end, coarse a, coarse b, fine a, fine b, ...
    # The coarse pair disagrees by 0.1 eV/A; the fine pair agrees.
    settings = coarse(1e-3, reaction_force_sigma=2e-3)
    result = run(tmp_path, settings=settings, reaction_offsets=(0.0, 0.0, 0.05, -0.05))
    assert result.status != "unresolved:branch_force_disagreement"
    first = _coarse(result.history)[:2]
    assert abs(first[0]["reaction_force"] - first[1]["reaction_force"]) > (
        settings.branch_force_limit()
    )
    assert first[0]["coarse_sign"] == first[1]["coarse_sign"] == 0
    assert any(h["origin"] == "from_left" for h in result.history)


def test_coarse_candidate_defers_to_the_v3_candidate_rule(tmp_path):
    # Symmetric bracket: the first midpoint is the exact TS. Offsets put both coarse
    # slopes beyond the sign bound while each coarse frame already meets the
    # candidate contract, so the coarse stage must not shrink the bracket.
    settings = EFSettings(
        midpoint_transverse_fmax_eV_A=1e-3,
        midpoint_coarse_transverse_fmax_eV_A=2e-3,
        midpoint_sign_bound_eV_A=3e-3,
    )
    result = run(
        tmp_path,
        frames=[frame(1.6), frame(2.0)],
        settings=settings,
        reaction_offsets=(0.0, 0.0, 0.005, 0.005),
    )
    first = _coarse(result.history)[:2]
    assert all(h["reaction_force"] > 3e-3 for h in first)
    assert all(h["fmax_eV_A"] <= settings.fmax_eV_A for h in first)
    assert first[0]["coarse_sign"] == first[1]["coarse_sign"] == 0
    assert result.status == "converged_candidate"
    assert result.metrics["origin"] in ("from_left", "from_right")


def test_fine_failure_after_coarse_stage_is_still_transverse_not_converged(tmp_path):
    # Below the float64 force floor of this surface: the fine stage cannot finish.
    settings = coarse(1e3, midpoint_transverse_fmax_eV_A=1e-30, max_relax_steps=30)
    result = run(tmp_path, settings=settings)
    assert result.status == "unresolved:transverse_not_converged"
    assert _coarse(result.history)
    assert result.history[-1]["origin"] in ("from_left", "from_right")


# Fine-stage FIRE (HF-EF-bracket-v5; docs/PLAN_HF_MIDPOINT_FINE_FIRE_2026-11-04.md).
@pytest.mark.parametrize(
    "kwargs",
    [
        {"midpoint_fine_optimizer": "lbfgs"},
        # FIRE is a fine-stage option; without the two-stage rule there is no fine stage.
        {"midpoint_fine_optimizer": "fire"},
    ],
)
def test_invalid_fine_optimizer(kwargs):
    with pytest.raises(ValueError):
        EFSettings(**kwargs)


def test_fire_runs_only_the_fine_stage(tmp_path):
    old = run(tmp_path / "v3")
    new = run(tmp_path / "v5", settings=coarse(1e3, midpoint_fine_optimizer="fire"))
    assert new.provenance["algorithm"] == "HF-EF-bracket-v5"
    assert new.status == "converged_candidate"
    assert new.metrics["q_actual"] == pytest.approx(old.metrics["q_actual"], abs=1e-6)
    assert new.metrics["optimizer"] == "fire"
    fine = [h for h in new.history if h["origin"] in ("from_left", "from_right")]
    other = [h for h in new.history if h["origin"] not in ("from_left", "from_right")]
    assert fine and all(h["optimizer"] == "fire" for h in fine)
    assert other and all(h["optimizer"] == "bfgs" for h in other)
    assert not any(h["curvature_reused"] for h in fine)


def _state_case(tmp_path):
    from explore.ts_ef import _relax, _RelaxState

    calls = []
    ledger = Ledger(tmp_path / "ledger.jsonl", arm="HF", budget=5000)
    factory = counted_factory(lambda: CoupledBonds(calls), ledger)
    atoms = frame(1.75, 2.05)
    q = INTENT.value(atoms.positions)
    state = _RelaxState(np.eye(9) * 7.0, ((q, atoms.positions.copy()),))
    return _relax, atoms, q, factory, ledger, state


def test_fire_acceptance_keeps_geometry_and_clears_bfgs_curvature(tmp_path):
    _relax, atoms, q, factory, ledger, state = _state_case(tmp_path)
    history = []
    with ledger.segment(task="t", stage="s"):
        last = _relax(
            atoms,
            INTENT,
            atoms.positions,
            q,
            "from_left",
            factory,
            CONFIG,
            SETTINGS,
            history,
            None,
            state=state,
            optimizer="fire",
        )
    assert state.hessian is None
    assert np.array_equal(state.points[-1][1], last["positions"])
    assert all(h["optimizer"] == "fire" and not h["curvature_reused"] for h in history)


def test_failed_fire_trial_leaves_the_state_untouched(tmp_path):
    from explore.ts_ef import _Stop

    _relax, atoms, q, factory, ledger, state = _state_case(tmp_path)
    before = state.copy()
    cfg = replace(SETTINGS, max_relax_steps=2)
    with ledger.segment(task="t", stage="s"), pytest.raises(_Stop):
        _relax(
            atoms,
            INTENT,
            atoms.positions,
            q,
            "from_left",
            factory,
            CONFIG,
            cfg,
            [],
            None,
            1e-30,
            state=state,
            optimizer="fire",
        )
    assert np.array_equal(state.hessian, before.hessian)
    assert len(state.points) == len(before.points)
    assert all(np.array_equal(a[1], b[1]) for a, b in zip(state.points, before.points))


# Coarse-stage candidate (HF-EF-bracket-v6; docs/PLAN_HF_MIDPOINT_COARSE_CANDIDATE_2026-11-06.md).
V5 = EFSettings(
    midpoint_transverse_fmax_eV_A=1e-3,
    midpoint_coarse_transverse_fmax_eV_A=2e-3,
    midpoint_sign_bound_eV_A=3e-3,
    midpoint_fine_optimizer="fire",
)
V6 = replace(V5, midpoint_coarse_candidate=True)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"midpoint_coarse_candidate": 1},
        # No two-stage rule at all, and two-stage with a BFGS fine stage (v4).
        {"midpoint_coarse_candidate": True},
        {
            "midpoint_coarse_candidate": True,
            "midpoint_coarse_transverse_fmax_eV_A": 4e-3,
            "midpoint_sign_bound_eV_A": 5e-3,
        },
    ],
)
def test_invalid_coarse_candidate_settings(kwargs):
    with pytest.raises(ValueError):
        EFSettings(**kwargs)


def test_coarse_pair_with_a_candidate_ends_the_search_on_the_coarse_stage(tmp_path):
    # Same input as test_coarse_candidate_defers_to_the_v3_candidate_rule: both coarse
    # frames meet the candidate contract and agree, so v6 stops there; v5 does not.
    kwargs = dict(frames=[frame(1.6), frame(2.0)], reaction_offsets=(0.0, 0.0, 0.005, 0.005))
    old = run(tmp_path / "v5", settings=V5, **kwargs)
    new = run(tmp_path / "v6", settings=V6, **kwargs)
    assert old.provenance["algorithm"] == "HF-EF-bracket-v5"
    assert new.provenance["algorithm"] == "HF-EF-bracket-v6"
    assert any(h["origin"] in ("from_left", "from_right") for h in old.history)
    assert new.status == "converged_candidate"
    assert new.metrics["origin"] in ("from_left_coarse", "from_right_coarse")
    assert new.metrics["coarse_candidate"] and new.metrics["coarse_sign"] == 0
    assert not any(h["origin"] in ("from_left", "from_right") for h in new.history)
    assert new.metrics["fmax_eV_A"] <= V6.fmax_eV_A
    assert new.metrics["fmax_projected"] <= V6.transverse_fmax_eV_A
    assert new.metrics["constraint_residual"] <= V6.constraint_tol
    assert new.spent < old.spent
    assert new.direction_cartesian is not None


def test_coarse_candidate_still_needs_the_paired_force_check(tmp_path):
    # The coarse pair disagrees beyond the gate while the left frame meets the
    # candidate contract: no coarse candidate, the fine stage runs as in v5.
    settings = replace(V6, reaction_force_sigma=1e-3)
    result = run(
        tmp_path,
        frames=[frame(1.6), frame(2.0)],
        settings=settings,
        reaction_offsets=(0.0, 0.0, 0.002, 0.009),
    )
    first = _coarse(result.history)[:2]
    assert first[0]["fmax_eV_A"] <= settings.fmax_eV_A
    assert abs(first[0]["reaction_force"] - first[1]["reaction_force"]) > (
        settings.branch_force_limit()
    )
    assert not any(h.get("coarse_candidate") for h in first)
    assert any(h["origin"] in ("from_left", "from_right") for h in result.history)


# Endpoint candidate (HF-EF-bracket-v7; docs/PLAN_HF_ENDPOINT_CANDIDATE_2026-11-07.md).
V7 = replace(V6, endpoint_candidate=True)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"endpoint_candidate": 1},
        {"endpoint_candidate": True},
        {
            "endpoint_candidate": True,
            "midpoint_transverse_fmax_eV_A": 1e-3,
            "midpoint_coarse_transverse_fmax_eV_A": 2e-3,
            "midpoint_sign_bound_eV_A": 3e-3,
            "midpoint_fine_optimizer": "fire",
        },
    ],
)
def test_invalid_endpoint_candidate_settings(kwargs):
    with pytest.raises(ValueError):
        EFSettings(**kwargs)


def test_end_at_the_maximum_is_the_candidate(tmp_path):
    # The right end sits on the TS: its slope is unresolved, the left end is uphill.
    kwargs = dict(frames=[frame(1.6), frame(1.8)])
    old = run(tmp_path / "v6", settings=V6, **kwargs)
    new = run(tmp_path / "v7", settings=V7, **kwargs)
    assert old.status == "unresolved:no_resolved_maximum_bracket"
    assert new.provenance["algorithm"] == "HF-EF-bracket-v7"
    assert new.status == "converged_candidate"
    assert new.metrics["origin"] == "right" and new.metrics["endpoint_candidate"]
    assert new.metrics["fmax_eV_A"] <= V7.fmax_eV_A
    assert new.metrics["fmax_projected"] <= V7.transverse_fmax_eV_A
    assert new.metrics["constraint_residual"] <= V7.constraint_tol
    assert not any(h["origin"].startswith("from_") for h in new.history)
    assert new.spent == old.spent


def test_end_at_a_minimum_along_q_is_not_a_candidate(tmp_path):
    # Minimum profile: the right end is a stationary minimum, the left end is downhill.
    result = run(tmp_path, frames=[frame(1.6), frame(1.8)], settings=V7, profile="minimum")
    left = [h for h in result.history if h["origin"] == "left"][-1]
    right = [h for h in result.history if h["origin"] == "right"][-1]
    assert left["reaction_force"] < -V7.reaction_force_tol
    assert right["fmax_eV_A"] <= V7.fmax_eV_A
    assert result.status == "unresolved:no_resolved_maximum_bracket"
    assert result.positions is None
