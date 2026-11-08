"""Paid E/F savings and branch-state ownership, using invariant analytic surfaces."""

from dataclasses import replace
from contextlib import closing

import numpy as np
import pytest

from explore.ledger import BudgetExhausted, Ledger, counted_factory
from explore.ts_branch import BranchSettings, _walk
from explore.ts_ef import EFSettings, _RelaxState, _relax, _Stop
import explore.ts_ef as ts_ef

from test_ts_branch import CONFIG, INTENT, SETTINGS, Surface, frame, run


def test_curvature_and_predictor_save_paid_calls_at_the_same_stationary_point(tmp_path):
    results = {}
    for name, reuse, predict in (
        ("cold", False, False),
        ("curvature", True, False),
        ("predictor", False, True),
        ("warm", True, True),
    ):
        settings = replace(SETTINGS, reuse_curvature=reuse, predict_tangent=predict)
        result = run(
            tmp_path / name,
            "smooth",
            [frame(-0.47, -0.141), frame(0.5, 0.15)],
            settings=settings,
        )
        results[name] = result
        assert result.status == "converged_candidate", (name, result.status)
        pos = np.asarray(result.refine["positions"])
        assert np.linalg.norm(pos[1] - pos[0]) == pytest.approx(1.8, abs=2e-4)
        assert np.linalg.norm(pos[2] - pos[1]) == pytest.approx(2.0, abs=2e-4)
        assert result.refine["energy_eV"] == pytest.approx(0.4, abs=1e-7)
        assert np.linalg.norm(result.refine["forces_pes"], axis=1).max() <= settings.fmax_eV_A
        history = result.history + result.refine["history"]
        assert max(p["predictor_step_A"] for p in history) <= settings.max_step_A * (1 + 1e-12)
        if not reuse:
            assert not any(p["curvature_reused"] for p in history)
        if not predict:
            assert all(p["predictor_step_A"] == 0 for p in history)
        if reuse:
            # The accepted branch's curvature reaches the locator, including its
            # zero-step endpoint rechecks and both independent midpoint solvers.
            assert all(p["curvature_reused"] for p in result.refine["history"])
    for name in ("curvature", "predictor", "warm"):
        assert results[name].spent < results["cold"].spent
    assert any(p["predictor_step_A"] > 0 for p in results["warm"].history)


def test_cold_and_warm_continuation_both_report_real_hysteresis(tmp_path):
    results = []
    for reuse in (False, True):
        result = run(
            tmp_path / str(reuse),
            "fold",
            [frame(-1.2, 1.2416), frame(1.2, -1.2416)],
            settings=replace(SETTINGS, reuse_curvature=reuse, predict_tangent=reuse),
        )
        assert result.status == "unresolved:hysteresis"
        assert result.refine is None
        assert all(s["status"] == "terminated" for s in result.sides)
        results.append(result)
    for cold, warm in zip(results[0].sides, results[1].sides):
        assert abs(cold["terminated_at_q"] - warm["terminated_at_q"]) < 0.01


def test_rejected_step_does_not_seed_the_halved_retry_or_the_other_endpoint():
    start_atoms = frame(-0.5, -0.15)
    start = dict(
        origin="low", q_actual=1.5, positions=start_atoms.positions, reaction_force=1.0
    )
    state = _RelaxState(np.eye(9), ((1.5, start_atoms.positions.copy()),))
    incoming = []

    def trial(positions, q, origin, *, state):
        incoming.append(state.copy())
        x = (q - 1.8) / 0.6
        pos = frame(x, 0.3 * x).positions
        if len(incoming) == 1:
            pos[2, 0] += 1.0  # A converged but wrong-branch trial; must be rejected.
        state.hessian[:] = 99
        state.points = (state.points[-1], (q, pos.copy()))
        return dict(origin=origin, q_actual=q, positions=pos, reaction_force=-1.0)

    side, found = _walk(
        1,
        start,
        1.6,
        trial,
        SETTINGS,
        BranchSettings(geometry_noise_A=1e-3),
        state,
    )
    assert side["status"] == "maximum_bracketed"
    assert len(incoming) == 2
    assert incoming[1].points[-1][0] == 1.5
    np.testing.assert_array_equal(incoming[1].hessian, np.eye(9))
    np.testing.assert_array_equal(state.hessian, np.eye(9))
    points, states = found
    assert points[1]["q_actual"] == pytest.approx(1.505)
    assert states[0].points[-1][0] == 1.5
    assert states[1].points[-1][0] == pytest.approx(1.505)
    assert not np.shares_memory(states[0].hessian, states[1].hessian)


@pytest.mark.parametrize("failure", ["budget", "ceiling"])
def test_failed_relaxation_keeps_paid_frames_and_does_not_commit_hints(tmp_path, failure):
    atoms = frame(-0.5, -0.15)
    state = _RelaxState(np.eye(9), ((1.5, atoms.positions.copy()),))
    before = state.copy()
    calls, history = [], []
    with closing(Ledger(tmp_path / "ledger.jsonl", "HF", budget=1)) as ledger:
        factory = counted_factory(lambda: Surface(calls, "smooth"), ledger)
        with ledger.segment(task="failure", stage="ts_search_ef"):
            with pytest.raises(BudgetExhausted if failure == "budget" else _Stop):
                _relax(
                    atoms,
                    INTENT,
                    atoms.positions,
                    1.6,
                    "low",
                    factory,
                    CONFIG,
                    SETTINGS,
                    history,
                    -1.0 if failure == "ceiling" else None,
                    state=state,
                )
    assert ledger.spent == len(calls) == len(history) == 1
    np.testing.assert_array_equal(state.hessian, before.hessian)
    assert len(state.points) == 1 and state.points[0][0] == before.points[0][0]
    np.testing.assert_array_equal(state.points[0][1], before.points[0][1])


def test_projected_curvature_survives_zero_step_but_not_a_changed_frame(tmp_path, monkeypatch):
    atoms = frame(-0.5, -0.15)
    state = _RelaxState(np.diag(np.arange(1.0, 10.0)), ((1.5, atoms.positions.copy()),))
    learned = state.hessian.copy()
    optimizers = []
    bfgs = ts_ef.BFGS

    def record(*args, **kwargs):
        opt = bfgs(*args, **kwargs)
        assert opt.pos0 is None and opt.forces0 is None
        optimizers.append(opt)
        return opt

    monkeypatch.setattr(ts_ef, "BFGS", record)
    calls, history = [], []
    with closing(Ledger(tmp_path / "ledger.jsonl", "HF", budget=10)) as ledger:
        factory = counted_factory(lambda: Surface(calls, "smooth"), ledger)
        with ledger.segment(task="zero", stage="ts_search_ef"):
            result = _relax(
                atoms,
                INTENT,
                atoms.positions,
                1.5,
                "low",
                factory,
                CONFIG,
                SETTINGS,
                history,
                None,
                state=state,
            )
            assert result["relax_step"] == 0 and result["curvature_reused"]
            assert state.hessian is not None
            n = INTENT.gradient(result["positions"]).ravel()
            n /= np.linalg.norm(n)
            tangent = np.eye(9) - np.outer(n, n)
            h0 = optimizers[0].H0
            np.testing.assert_allclose(h0 @ n, optimizers[0].alpha * n, atol=1e-12)
            np.testing.assert_allclose(
                tangent @ h0 @ tangent, tangent @ learned @ tangent, atol=1e-12
            )
            changed = result["positions"] + [2.0, -1.0, 3.0]
            result = _relax(
                atoms,
                INTENT,
                changed,
                1.5,
                "low",
                factory,
                CONFIG,
                SETTINGS,
                history,
                None,
                state=state,
            )
            assert result["relax_step"] == 0 and not result["curvature_reused"]
            assert result["predictor_step_A"] == 0
            assert state.hessian is None and len(state.points) == 1
    assert ledger.spent == len(calls) == 2


@pytest.mark.parametrize("kwargs", [{"reuse_curvature": 1}, {"predict_tangent": "yes"}])
def test_warm_start_switches_require_booleans(kwargs):
    with pytest.raises(ValueError, match="must be a bool"):
        EFSettings(**kwargs)
