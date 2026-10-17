import json

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.lj import LennardJones

from explore.ledger import Ledger, counted_factory
from explore.probe import (
    Intent,
    ScanPoint,
    ScanSettings,
    bracket_jump,
    constrained_scan,
    read_scan,
    release,
)
from prrs.calculators import double_well_factory
from prrs.config import SearchConfig
from prrs.io import snapshot
from prrs.reliability import GuardedCalculator

CONFIG = SearchConfig(closed_shell_only=False, quench_fmax_eV_A=0.001)


def dw_energy(r):
    x = (r - 1.8) / 0.6
    return 0.4 * (x * x - 1) ** 2


def dw_slope(r):
    x = (r - 1.8) / 0.6
    return 1.6 * x * (x * x - 1) / 0.6


def test_intent_gradient_matches_finite_difference():
    rng = np.random.default_rng(0)
    pos = rng.normal(size=(4, 3)) * 1.5
    intent = Intent("i", ((0, 1), (2, 3)), (1.0, -0.5))
    grad = intent.gradient(pos)
    h = 1e-6
    for a in range(4):
        for k in range(3):
            plus, minus = pos.copy(), pos.copy()
            plus[a, k] += h
            minus[a, k] -= h
            fd = (intent.value(plus) - intent.value(minus)) / (2 * h)
            assert abs(fd - grad[a, k]) < 1e-6


def scan_double_well(tmp_path, budget=1000, deltas=None):
    ledger = Ledger(tmp_path / "l.jsonl", arm="D2", budget=budget)
    atoms = Atoms("C2", positions=[[0, 0, 0], [1.2, 0, 0]])
    settings = ScanSettings(
        deltas_A=deltas or tuple(round(0.1 * k, 10) for k in range(11)),
        fmax_eV_A=1e-4,
        max_steps_per_point=50,
        max_rise_eV=10.0,
    )
    result = constrained_scan(
        atoms,
        Intent("cc", ((0, 1),), (1.0,)),
        settings,
        counted_factory(double_well_factory, ledger),
        CONFIG,
        ledger,
        task="dw",
    )
    return result, ledger


def test_double_well_scan_reproduces_the_analytic_profile(tmp_path):
    result, ledger = scan_double_well(tmp_path)
    assert result.status == "completed" and len(result.points) == 11
    for p in result.points:
        r = 1.2 + p.delta_A
        assert abs(p.q_actual - r) < 1e-6
        assert (
            abs((p.energy_eV - result.points[0].energy_eV) - (dw_energy(r) - dw_energy(1.2)))
            < 1e-8
        )
        assert abs(p.reaction_force - dw_slope(r)) < 1e-6
        assert p.forces_pes.shape == (2, 3)
        grad = Intent("cc", ((0, 1),), (1.0,)).gradient(p.positions)
        assert abs(-np.vdot(p.forces_pes, grad) / 2.0 - p.reaction_force) < 1e-12
    assert ledger.spent > 0


def test_unconstrained_degrees_of_freedom_relax(tmp_path):
    ledger = Ledger(tmp_path / "l.jsonl", arm="D2", budget=5000)
    atoms = Atoms("Ar3", positions=[[0, 0, 0], [1.12, 0, 0], [0.3, 1.6, 0]])
    start = atoms.positions[2].copy()
    settings = ScanSettings(
        deltas_A=(0.0, 0.1), fmax_eV_A=1e-3, max_steps_per_point=500, max_rise_eV=10.0
    )
    result = constrained_scan(
        atoms,
        Intent("ar", ((0, 1),), (1.0,)),
        settings,
        counted_factory(lambda: LennardJones(sigma=1.0, epsilon=0.01, rc=10.0), ledger),
        CONFIG,
        ledger,
        task="lj",
    )
    last = result.points[-1]
    assert last.converged and last.fmax_projected < 1e-3
    assert np.linalg.norm(last.positions[2] - start) > 0.01


def test_energy_cap_censors_and_keeps_points(tmp_path):
    ledger = Ledger(tmp_path / "l.jsonl", arm="D2", budget=1000)
    atoms = Atoms("C2", positions=[[0, 0, 0], [1.2, 0, 0]])
    settings = ScanSettings(
        tuple(round(0.1 * k, 10) for k in range(11)), 1e-4, 50, max_rise_eV=0.2
    )
    result = constrained_scan(
        atoms,
        Intent("cc", ((0, 1),), (1.0,)),
        settings,
        counted_factory(double_well_factory, ledger),
        CONFIG,
        ledger,
        task="cap",
    )
    assert result.status == "censored_energy"
    assert result.points[-1].energy_eV - result.points[0].energy_eV > 0.2
    assert all(p.energy_eV - result.points[0].energy_eV <= 0.2 for p in result.points[:-1])


def test_budget_exhaustion_mid_scan_keeps_points(tmp_path):
    result, ledger = scan_double_well(tmp_path, budget=3)
    assert result.status == "budget_exhausted"
    assert ledger.spent == 3


def test_overlap_gate_ends_scan_without_raising(tmp_path):
    deltas = tuple(round(-0.1 * k, 10) for k in range(12))  # compress 1.2 A towards 0
    result, _ = scan_double_well(tmp_path, deltas=deltas)
    assert result.status == "gate_rejected:overlap"
    assert len(result.points) >= 1


def test_negative_target_distance_is_infeasible_before_any_forward(tmp_path):
    result, ledger = scan_double_well(tmp_path, deltas=(0.0, -1.5))
    assert result.status == "infeasible" and result.points == [] and ledger.spent == 0


def test_first_delta_must_be_zero():
    with pytest.raises(ValueError):
        ScanSettings(
            deltas_A=(0.1, 0.2), fmax_eV_A=1e-3, max_steps_per_point=10, max_rise_eV=1.0
        )


def fake_point(delta, energy, rf, edges, converged=True):
    return ScanPoint(delta, 0.0, 0.0, energy, rf, 0.0, converged, 0, edges, np.zeros((2, 3)))


def test_double_well_scan_reads_one_ridge_at_the_barrier(tmp_path):
    result, _ = scan_double_well(tmp_path)
    triggers = read_scan(result.points, noise_eV=1e-3)
    ridges = [t for t in triggers if t["kind"] == "ridge_suspected"]
    assert len(ridges) == 1
    assert abs(result.points[ridges[0]["peak"]].q_actual - 1.8) < 1e-6
    assert ridges[0]["release"] == ridges[0]["peak"] + 1


def test_monotone_rise_reads_nothing():
    pts = [fake_point(0.1 * k, 0.1 * k, 1.0, [[0, 1]]) for k in range(6)]
    assert read_scan(pts, noise_eV=1e-3) == []


def test_peak_below_noise_is_not_a_ridge():
    es = [0.0, 0.0005, 0.0008, 0.0005, 0.0]
    rfs = [1.0, 1.0, 0.0, -1.0, -1.0]
    pts = [fake_point(0.1 * k, e, f, [[0, 1]]) for k, (e, f) in enumerate(zip(es, rfs))]
    assert read_scan(pts, noise_eV=1e-3) == []


def test_unconverged_point_before_the_peak_blocks_the_ridge():
    es, rfs = [0.0, 0.2, 0.4, 0.2, 0.0], [1.0, 1.0, 0.0, -1.0, -1.0]
    pts = [
        fake_point(0.1 * k, e, f, [[0, 1]], converged=(k != 1))
        for k, (e, f) in enumerate(zip(es, rfs))
    ]
    assert [t for t in read_scan(pts, 1e-3) if t["kind"] == "ridge_suspected"] == []


def test_graph_change_with_energy_drop_is_a_jump():
    pts = [
        fake_point(0.0, 0.0, 1.0, [[0, 1]]),
        fake_point(0.1, 0.5, 1.0, [[0, 1]]),
        fake_point(0.2, -0.3, -1.0, [[0, 1], [1, 2]]),
    ]
    assert {"kind": "jump", "before": 1, "after": 2} in read_scan(pts, noise_eV=1e-3)


def test_reverse_point_staying_on_the_new_graph_brackets_the_jump():
    pts = [fake_point(0.0, 0.0, 1.0, [[0, 1]]), fake_point(0.1, -0.3, -1.0, [[0, 1], [1, 2]])]
    jump = {"kind": "jump", "before": 0, "after": 1}
    stayed = fake_point(0.0, -0.1, -1.0, [[0, 1], [1, 2]])
    went_back = fake_point(0.0, 0.0, 1.0, [[0, 1]])
    assert bracket_jump(pts, jump, stayed)["kind"] == "branch_bracketed"
    assert bracket_jump(pts, jump, went_back) is None


def dw_source():
    atoms = Atoms("C2", positions=[[0, 0, 0], [1.2, 0, 0]])
    return snapshot(atoms, GuardedCalculator(double_well_factory(), CONFIG))


def test_release_past_the_ridge_lands_in_the_far_basin(tmp_path):
    ledger = Ledger(tmp_path / "l.jsonl", arm="D2", budget=5000)
    template = Atoms("C2", positions=np.zeros((2, 3)))
    out = release(
        np.array([[0, 0, 0], [2.0, 0, 0]], float),
        template,
        dw_source(),
        counted_factory(double_well_factory, ledger),
        CONFIG,
        ledger,
        task="rel",
        seed=0,
    )
    assert out["status"] == "c1_candidate"
    assert abs(out["endpoint"].get_distance(0, 1) - 2.4) < 1e-2
    assert out["label"] == "dissociation"
    stages = [json.loads(x)["stage"] for x in (tmp_path / "l.jsonl").read_text().splitlines()]
    assert "quench" in stages and "minimum_check" in stages
    assert "unattributed" not in stages
    assert out["seed"] == 0
    assert out["constraints_removed"] == []
    assert np.allclose(out["release_positions"], [[0, 0, 0], [2.0, 0, 0]])
    assert out["release_forces_pes"].shape == (2, 3)


def test_release_before_the_ridge_returns_to_source(tmp_path):
    ledger = Ledger(tmp_path / "l.jsonl", arm="D2", budget=5000)
    template = Atoms("C2", positions=np.zeros((2, 3)))
    out = release(
        np.array([[0, 0, 0], [1.6, 0, 0]], float),
        template,
        dw_source(),
        counted_factory(double_well_factory, ledger),
        CONFIG,
        ledger,
        task="rel",
        seed=0,
    )
    assert out["label"] == "elastic"
