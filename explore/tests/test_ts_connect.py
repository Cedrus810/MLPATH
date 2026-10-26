"""Analytic double basin: E/F connection-tracing contracts, not molecular acceptance."""

import json
from types import SimpleNamespace

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes

from explore.ledger import Ledger, counted_factory
from explore.probe import Intent
from explore.ts_connect import ConnectSettings, branch_starts, trace_sides_ef
from explore.ts_ef import EFSettings, refine_ts_ef
from prrs.config import SearchConfig

CONFIG = SearchConfig(closed_shell_only=False)
INTENT = Intent("r01", ((0, 1),), (1.0,))
SEARCH = EFSettings(fmax_eV_A=2e-4, transverse_fmax_eV_A=1e-5, reaction_force_tol=1e-4)
CONNECT = ConnectSettings(reaction_force_sigma=1e-3)


class DoubleBasin(Calculator):
    """E = 0.4 (x^2 - 1)^2 + 3 (y - 0.3 x)^2: minima x = -1 (r 1.2, C0-C1 bonded) and
    x = +1 (r 2.4, not bonded), index-one saddle at x = 0 (r 1.8)."""

    implemented_properties = ["energy", "forces"]

    def __init__(self, calls):
        super().__init__()
        self.calls = calls

    def get_hessian(self, *args, **kwargs):
        raise AssertionError("Hessian entered E/F connection tracing")

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        self.calls.append(atoms.positions.copy())
        p = atoms.positions
        v, w = p[1] - p[0], p[2] - p[1]
        r, s = np.linalg.norm(v), np.linalg.norm(w)
        x, y = (r - 1.8) / 0.6, (s - 2.0) / 0.4
        energy = 0.4 * (x * x - 1) ** 2 + 3 * (y - 0.3 * x) ** 2
        dx = 1.6 * x * (x * x - 1) - 1.8 * (y - 0.3 * x)
        dy = 6 * (y - 0.3 * x)
        dr, ds = dx / 0.6, dy / 0.4
        f = np.zeros_like(p)
        f[0] += dr * v / r
        f[1] -= dr * v / r
        f[1] += ds * w / s
        f[2] -= ds * w / s
        self.results = {"energy": energy, "forces": f}


def frame(x, y=None):
    y = 0.3 * x if y is None else y
    r, s = 1.8 + 0.6 * x, 2.0 + 0.4 * y
    return Atoms("C3", positions=[[0, 0, 0], [r, 0, 0], [r + s, 0, 0]])


def point(atoms):
    return SimpleNamespace(q_target=INTENT.value(atoms.positions), positions=atoms.positions)


def setup(tmp_path, budget=200000):
    calls = []
    ledger = Ledger(tmp_path / "ledger.jsonl", arm="HF", budget=budget)
    factory = counted_factory(lambda: DoubleBasin(calls), ledger)
    return calls, ledger, factory


def search(tmp_path):
    calls, ledger, factory = setup(tmp_path)
    found = refine_ts_ef(
        frame(-0.5),
        INTENT,
        [point(frame(-0.5)), point(frame(0.5))],
        factory,
        CONFIG,
        ledger,
        "search",
        settings=SEARCH,
    )
    assert found.status == "converged_candidate"
    return found, calls, ledger, factory


def trace(found, starts, direction, calls, ledger, factory, settings=CONNECT):
    before = ledger.spent
    result = trace_sides_ef(
        frame(0.0),
        found.positions,
        found.energy_eV,
        starts,
        direction,
        factory,
        CONFIG,
        ledger,
        "connect",
        settings=settings,
    )
    assert result.spent == ledger.spent - before
    json.dumps(result.to_dict(), allow_nan=False)
    return result


def bond(atoms_positions):
    p = np.asarray(atoms_positions)
    return np.linalg.norm(p[1] - p[0])


def test_branch_frames_start_both_sides_and_reach_both_minima(tmp_path):
    found, calls, ledger, factory = search(tmp_path)
    q_star = found.metrics["q_actual"]
    starts = branch_starts(
        found.history,
        q_star,
        CONNECT.reaction_force_sigma,
        descent_fmax=CONFIG.quench_fmax_eV_A,
    )
    assert starts[0]["q_actual"] < q_star < starts[1]["q_actual"]
    assert all(s["fmax_eV_A"] > CONFIG.quench_fmax_eV_A for s in starts)
    result = trace(found, starts, found.direction_cartesian, calls, ledger, factory)
    ledger.close()
    assert ledger.spent == len(calls)
    assert result.status == "connected"
    # One bond (C0-C1) breaks from the low to the high basin and nothing forms.
    assert result.label == "dissociation"
    assert result.fingerprint["broken"] == [[0, 1]] and not result.fingerprint["formed"]
    low, high = result.sides
    assert low["start"]["source"] == high["start"]["source"] == "branch_frame"
    assert bond(low["descent"]["positions"]) == pytest.approx(1.2, abs=0.02)
    assert bond(high["descent"]["positions"]) == pytest.approx(2.4, abs=0.02)
    assert low["status"] == high["status"] == "ef_stable"


def test_missing_frames_fall_back_to_the_search_tangent(tmp_path):
    found, calls, ledger, factory = search(tmp_path)
    result = trace(found, (None, None), found.direction_cartesian, calls, ledger, factory)
    assert result.status == "connected"
    for side in result.sides:
        assert side["start"]["source"] == "tangent"
        assert side["start"]["tried"]


def test_saddle_endpoint_is_caught_by_the_perturbation_check(tmp_path):
    # Hand the saddle itself in as a "start" on both sides: the descent stops at once on
    # a zero-force point, and only the perturbed re-descents reveal it is not a minimum.
    found, calls, ledger, factory = search(tmp_path)
    saddle = {
        "q_actual": found.metrics["q_actual"],
        "reaction_force": 0.0,
        "energy_eV": found.energy_eV,
        "positions": found.positions,
    }
    settings = ConnectSettings(reaction_force_sigma=1e-3, perturb_A=0.05, n_perturb=4)
    result = trace(found, (saddle, saddle), None, calls, ledger, factory, settings)
    assert result.status == "unresolved:endpoint_unstable"
    assert all(s["stability"]["verdict"] == "ef_unstable" for s in result.sides)


def test_budget_exhaustion_keeps_paid_forwards(tmp_path):
    found, calls, _, _ = search(tmp_path)
    calls2 = []
    ledger = Ledger(tmp_path / "small.jsonl", arm="HF", budget=5)
    factory = counted_factory(lambda: DoubleBasin(calls2), ledger)
    starts = branch_starts(found.history, found.metrics["q_actual"], 1e-3)
    result = trace(found, starts, None, calls2, ledger, factory)
    ledger.close()
    assert result.status == "budget_exhausted"
    assert result.spent == ledger.spent == len(calls2) == 5


@pytest.mark.parametrize(
    "kwargs",
    [
        {"reaction_force_sigma": 0.0},
        {"reaction_force_sigma": 1e-3, "tangent_amplitudes_A": ()},
        {"reaction_force_sigma": 1e-3, "n_perturb": -1},
        {"reaction_force_sigma": 1e-3, "max_descent_steps": 0},
    ],
)
def test_invalid_connect_settings(kwargs):
    with pytest.raises(ValueError):
        ConnectSettings(**kwargs)


def test_start_already_meeting_the_descent_tolerance_is_skipped():
    frames = [
        {"q_actual": 1.79, "reaction_force": 0.02, "fmax_eV_A": 0.02, "positions": None},
        {"q_actual": 1.70, "reaction_force": 0.20, "fmax_eV_A": 0.20, "positions": None},
        {"q_actual": 1.81, "reaction_force": -0.25, "fmax_eV_A": 0.25, "positions": None},
    ]
    low, high = branch_starts(frames, 1.80, 1e-3, descent_fmax=0.03)
    assert low["q_actual"] == 1.70 and high["q_actual"] == 1.81
    low, _ = branch_starts(frames, 1.80, 1e-3)
    assert low["q_actual"] == 1.79
