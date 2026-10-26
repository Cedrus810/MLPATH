"""Analytic fold surface: two-sided continuation contracts, not molecular acceptance."""

import json

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes

from explore.ledger import Ledger, counted_factory
from explore.probe import Intent
from explore.ts_branch import BranchSettings, continue_branches_ef
from explore.ts_ef import EFSettings
from prrs.config import SearchConfig

CONFIG = SearchConfig(closed_shell_only=False)
INTENT = Intent("r01", ((0, 1),), (1.0,))
SETTINGS = EFSettings(fmax_eV_A=2e-4, transverse_fmax_eV_A=1e-5, reaction_force_tol=1e-4)
BRANCH = BranchSettings(geometry_noise_A=1e-4)
# fold: E = 0.4[(y^2-1)^2 + 2xy + x^2/2]; for fixed x the transverse minima in y fold
# at |x| = 4/(3 sqrt 3); the only saddle (0, 0) lies on neither constrained-minimum
# branch, and dE/dx = 0.4(2y + x)/0.6 keeps its sign on each branch up to the fold.
X_FOLD = 4 / (3 * np.sqrt(3))


class Surface(Calculator):
    implemented_properties = ["energy", "forces"]

    def __init__(self, calls, kind):
        super().__init__()
        self.calls, self.kind = calls, kind

    def get_hessian(self, *args, **kwargs):
        raise AssertionError("Hessian entered E/F search")

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        self.calls.append(atoms.positions.copy())
        p = atoms.positions
        v, w = p[1] - p[0], p[2] - p[1]
        r, s = np.linalg.norm(v), np.linalg.norm(w)
        x, y = (r - 1.8) / 0.6, (s - 2.0) / 0.4
        if self.kind == "fold":
            energy = 0.4 * ((y * y - 1) ** 2 + 2 * x * y + 0.5 * x * x)
            dx, dy = 0.4 * (2 * y + x), 0.4 * (4 * y * (y * y - 1) + 2 * x)
        else:  # smooth maximum at x = 0 on the single branch y = 0.3 x
            energy = 0.4 * (1 - x * x) + 3 * (y - 0.3 * x) ** 2
            dx, dy = -0.8 * x - 1.8 * (y - 0.3 * x), 6 * (y - 0.3 * x)
        dr, ds = dx / 0.6, dy / 0.4
        f = np.zeros_like(p)
        f[0] += dr * v / r
        f[1] -= dr * v / r
        f[1] += ds * w / s
        f[2] -= ds * w / s
        self.results = {"energy": energy, "forces": f}


def frame(x, y):
    r, s = 1.8 + 0.6 * x, 2.0 + 0.4 * y
    return Atoms("C3", positions=[[0, 0, 0], [r, 0, 0], [r + s, 0, 0]])


def point(atoms):
    from types import SimpleNamespace

    return SimpleNamespace(q_target=INTENT.value(atoms.positions), positions=atoms.positions)


def run(tmp_path, kind, frames, budget=50000, settings=SETTINGS):
    calls = []
    ledger = Ledger(tmp_path / "ledger.jsonl", arm="HF", budget=budget)
    factory = counted_factory(lambda: Surface(calls, kind), ledger)
    result = continue_branches_ef(
        frames[0],
        INTENT,
        [point(f) for f in frames],
        factory,
        CONFIG,
        ledger,
        "test",
        branch_settings=BRANCH,
        settings=settings,
    )
    ledger.close()
    assert result.spent == ledger.spent == len(calls)
    json.dumps(result.to_dict(), allow_nan=False)
    return result


def test_fold_reports_hysteresis_with_both_termination_points(tmp_path):
    # Low end on the y > 0 branch, high end on the y < 0 branch: signs bracket a
    # maximum in q, but no continuous branch carries it.
    result = run(tmp_path, "fold", [frame(-1.2, 1.2416), frame(1.2, -1.2416)])
    assert result.status == result.branch_status == "unresolved:hysteresis"
    assert result.refine is None
    low, high = result.sides
    assert low["status"] == high["status"] == "terminated"
    # a converged trial that jumps is still halved down to min_step_q, then ends the side
    assert low["termination"] == high["termination"] == "jump"
    r_up, r_down = 1.8 + 0.6 * X_FOLD, 1.8 - 0.6 * X_FOLD
    assert r_up - 0.03 <= low["terminated_at_q"] <= r_up + 1e-3
    assert r_down - 1e-3 <= high["terminated_at_q"] <= r_down + 0.03
    for side in (low, high):
        assert side["jump"]["rmsd_A"] > 100 * BRANCH.geometry_noise_A
        # The reaction force kept its sign up to the fold on each branch.
        signs = {np.sign(p["reaction_force"]) for p in side["points"]}
        assert len(signs) == 1


def test_smooth_branch_hands_a_same_branch_pair_to_the_locator(tmp_path):
    result = run(tmp_path, "smooth", [frame(-0.5, -0.15), frame(0.5, 0.15)])
    assert result.branch_status == "continuous_maximum"
    assert result.status == "converged_candidate"
    assert result.sides[-1]["status"] == "maximum_bracketed"
    pos = np.asarray(result.refine["positions"])
    assert np.linalg.norm(pos[1] - pos[0]) == pytest.approx(1.8, abs=1e-3)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"geometry_noise_A": 0.0},
        {"geometry_noise_A": 1e-3, "min_step_q": 0.02},
        {"geometry_noise_A": 1e-3, "jump_factor": 1.0},
        {"geometry_noise_A": 1e-3, "window": 2.5},
    ],
)
def test_invalid_branch_settings(kwargs):
    with pytest.raises(ValueError):
        BranchSettings(**kwargs)


def test_a_branch_without_maximum_stops_at_the_far_end(tmp_path):
    # Both ends uphill of the smooth maximum: each side walks the whole interval once.
    result = run(tmp_path, "smooth", [frame(-0.8, -0.24), frame(-0.2, -0.06)])
    assert result.status == "unresolved:hysteresis"
    span = 0.6 * 0.6 / BRANCH.step_q
    for side in result.sides:
        assert side["status"] == "reached_end"
        assert len(side["points"]) <= 2 * span  # halvings allowed; the bug gave 400


# --- trial relaxations that do not converge (docs/PLAN_HF_WALK_RETRY_2026-10-28.md) ---


def _fake_relax(history, calls, fails, rf=lambda q: 1.0, failed_rf=None):
    """Relaxation on a straight branch (atom 2 moves with q). ``fails(q_from, q)`` decides
    a non-converging trial: it appends a diagnostic frame (as _relax does), poisons the
    trial state it was given, and raises like _relax."""
    from explore.ts_ef import _Stop

    def relax(positions, q, origin, *, state):
        q_from = float(positions[1][0])
        calls.append({"q_from": q_from, "q": q, "hessian": state.hessian})
        pos = np.array([[0.0, 0, 0], [q, 0, 0], [q + 2.0, 0.5 * q, 0]])
        frame = {
            "q_target": q,
            "q_actual": q,
            "energy_eV": -q,
            "reaction_force": rf(q),
            "fmax_eV_A": 0.1,
            "fmax_projected": 0.05,
            "constraint_residual": 0.0,
            "positions": pos,
            "origin": origin,
            "relax_step": 200,
        }
        if fails(q_from, q):
            frame["failed"] = True
            if failed_rf is not None:
                frame["reaction_force"] = failed_rf
            history.append(frame)
            state.hessian = np.array([-1.0])  # must never reach the accepted state
            raise _Stop("transverse_not_converged")
        history.append(frame)
        state.hessian = np.array([q])
        return frame

    return relax


def _start(q):
    return {
        "q_actual": q,
        "q_target": q,
        "positions": np.array([[0.0, 0, 0], [q, 0, 0], [q + 2.0, 0.5 * q, 0]]),
        "reaction_force": 1.0,
        "energy_eV": -q,
        "origin": "low",
    }


def test_an_unconverged_trial_ends_the_side_without_a_retry():
    from explore.ts_branch import _walk
    from explore.ts_ef import _RelaxState

    history, calls = [], []
    # would converge at a smaller step; v4 does not retry (b05 rerun, 2026-10-29)
    relax = _fake_relax(history, calls, fails=lambda a, b: abs(b - a) > 0.003)
    bcfg = BranchSettings(geometry_noise_A=1e-3)
    side, found = _walk(
        +1, _start(1.0), 1.05, relax, EFSettings(), bcfg, _RelaxState(), history
    )
    assert found is None
    assert side["status"] == "terminated" and side["termination"] == "relax_unconverged"
    assert len(calls) == 1 and len(side["unconverged_trials"]) == 1
    assert side["unconverged_trials"][0]["step_q"] == bcfg.step_q
    assert side["terminated_at_q"] == 1.0 and len(side["points"]) == 1
    assert any(f.get("failed") for f in history)  # charged and kept as diagnostics
    assert not any(p.get("failed") for p in side["points"])


def test_the_side_ends_at_the_first_unconverged_trial_after_accepted_steps():
    from explore.ts_branch import _walk
    from explore.ts_ef import _RelaxState

    history, calls = [], []
    relax = _fake_relax(history, calls, fails=lambda a, b: b > 1.02)
    bcfg = BranchSettings(geometry_noise_A=1e-3)
    side, found = _walk(
        +1, _start(1.0), 1.05, relax, EFSettings(), bcfg, _RelaxState(), history
    )
    assert found is None
    assert side["status"] == "terminated" and side["termination"] == "relax_unconverged"
    assert len(side["unconverged_trials"]) == 1  # no halving after the failure
    assert len(side["points"]) >= 2 and side["terminated_at_q"] <= 1.02
    assert side["terminated_at_q"] == side["points"][-1]["q_actual"]
    # the accepted state never inherited the failed trial's curvature
    assert all(c["hessian"] is None or c["hessian"][0] >= 0 for c in calls)
    diag = side["relax_unconverged"]
    assert {"fmax_projected", "dE_eV", "rmsd_A", "positions"} <= set(diag)
    assert "jump" not in side


def test_a_failed_trial_never_brackets_a_maximum():
    from explore.ts_branch import _walk
    from explore.ts_ef import _RelaxState

    history, calls = [], []
    # accepted points stay uphill; the failing trials report a downhill force
    relax = _fake_relax(history, calls, fails=lambda a, b: b > 1.02, failed_rf=-1.0)
    side, found = _walk(
        +1,
        _start(1.0),
        1.05,
        relax,
        EFSettings(),
        BranchSettings(geometry_noise_A=1e-3),
        _RelaxState(),
        history,
    )
    assert found is None and side["status"] == "terminated"
    assert all(p["reaction_force"] > 0 for p in side["points"])


@pytest.mark.parametrize("reason", ["above_ceiling", "constraint_not_converged", "budget"])
def test_other_stops_still_end_the_search(reason):
    from explore.ledger import BudgetExhausted
    from explore.ts_branch import _walk
    from explore.ts_ef import _RelaxState, _Stop

    def relax(positions, q, origin, *, state):
        if reason == "budget":
            raise BudgetExhausted("test")
        raise _Stop(reason)

    expected = BudgetExhausted if reason == "budget" else _Stop
    with pytest.raises(expected):
        _walk(
            +1,
            _start(1.0),
            1.05,
            relax,
            EFSettings(),
            BranchSettings(geometry_noise_A=1e-3),
            _RelaxState(),
            [],
        )


def test_unconverged_walk_end_to_end_is_charged_and_not_called_hysteresis(tmp_path):
    # Both ends sit exactly on the smooth branch (no relaxation needed); trial steps need
    # more than one BFGS step at this tolerance, so every trial fails.
    tight = EFSettings(
        fmax_eV_A=2e-4,
        transverse_fmax_eV_A=1e-5,
        reaction_force_tol=1e-4,
        max_relax_steps=1,
        predict_tangent=False,
    )
    result = run(tmp_path, "smooth", [frame(-0.5, -0.15), frame(0.5, 0.15)], settings=tight)
    assert result.status == result.branch_status == "unresolved:walk_relax_unconverged"
    assert result.refine is None
    for side in result.sides:
        assert side["termination"] == "relax_unconverged"
        assert len(side["points"]) == 1  # only the accepted end
    assert result.provenance["algorithm"] == "HF-EF-branches-v4"
