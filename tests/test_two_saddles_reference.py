"""Independent geometry/derivative evidence for the analytic benchmark."""

import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "pars" / "src"))
from benchmarks.pars_two_saddles.run import (  # noqa: E402
    PROTOCOL,
    TwoSaddles,
    atoms_at,
    geometry,
    population,
    references,
)


def test_two_independent_index_one_saddles_connect_the_same_confirmed_minima():
    result = references(json.loads(PROTOCOL.read_text()))
    assert len(result) == 2
    assert result[0]["positions_A"][0][1] < 0 < result[1]["positions_A"][0][1]
    assert all({r["minimum"] for r in saddle["descents"]} == {-1, 1} for saddle in result)
    assert all(sum(value < -1e-4 for value in saddle["eigenvalues"]) == 1 for saddle in result)


def test_analytic_forces_match_independent_energy_differences():
    potential = json.loads(PROTOCOL.read_text())["potential"]
    positions = geometry(0.2, 0.5) + np.random.default_rng(42).normal(0, 0.03, (3, 3))
    forces = atoms_at(positions, TwoSaddles(potential)).get_forces()
    for j in range(positions.size):
        plus, minus = positions.copy(), positions.copy()
        plus.flat[j] += 1e-6
        minus.flat[j] -= 1e-6
        gradient = (
            atoms_at(plus, TwoSaddles(potential)).get_potential_energy()
            - atoms_at(minus, TwoSaddles(potential)).get_potential_energy()
        ) / 2e-6
        assert -gradient == pytest.approx(forces.flat[j], abs=1e-8)


def test_pool_has_geometric_contact_diversity_without_oracle_labels():
    candidates = population(json.loads(PROTOCOL.read_text()))
    assert len(candidates) == 12
    assert len({tuple(c.signatures) for c in candidates}) > 1
    assert all("saddle" not in c.provenance.extra for c in candidates)
    assert all(c.provenance.frame_stop == len(c.positions) for c in candidates)


def test_common_force_budget_is_enforced_even_during_endpoint_checks():
    from pars.config import ParsConfig
    from pars.refinement.ase_neb import AseNebBackend

    protocol = json.loads(PROTOCOL.read_text())
    candidate = population(protocol)[0]
    calculator = TwoSaddles(protocol["potential"], limit=2)
    result = AseNebBackend(calculator).refine(candidate, ParsConfig(**protocol["refinement"]))
    assert result.status.value == "failed"
    assert calculator.force_evaluations == 2


def test_connectivity_descent_uses_the_same_hard_force_budget():
    import dataclasses

    from pars.config import ParsConfig
    from pars.refinement.ase_neb import AseNebBackend

    protocol = json.loads(PROTOCOL.read_text())
    candidate = population(protocol)[0]
    config = ParsConfig(**protocol["refinement"])
    # Locate the first independent kick with a legitimate zero-step policy.
    preparation = TwoSaddles(protocol["potential"])
    stopped = AseNebBackend(preparation).refine(
        candidate, dataclasses.replace(config, refinement_descent_max_steps=0)
    )
    assert "descent accepted-step budget exhausted" in stopped.message
    limit = preparation.force_evaluations + 1
    calculator = TwoSaddles(protocol["potential"], limit=limit)
    result = AseNebBackend(calculator).refine(candidate, config)
    assert result.status.value == "failed"
    assert result.barrier_eV is None
    assert "ForceBudgetExceeded" in result.message
    assert calculator.force_evaluations == limit
    branch = result.provenance["connectivity_check"]["branches"][0]
    assert branch["status"] == "failed"
    assert branch["accepted_steps"] == 1


def test_finer_connectivity_replay_cannot_reset_the_hard_force_budget():
    from pars.config import ParsConfig
    from pars.refinement.ase_neb import AseNebBackend

    protocol = json.loads(PROTOCOL.read_text())
    config = ParsConfig(**protocol["refinement"])
    potential = protocol["potential"]
    y = min(np.roots([4, 3 * potential["b"], 2 * potential["a"] - 4]))
    mode = np.zeros((3, 3))
    mode[0, 0] = 1

    def inputs(calculator):
        return (
            atoms_at(geometry(0, y), calculator),
            mode,
            [atoms_at(geometry(-1), calculator), atoms_at(geometry(1), calculator)],
        )

    preparation = TwoSaddles(potential)
    coarse = {"branches": []}
    AseNebBackend(preparation)._check_connectivity_at_scale(
        *inputs(preparation), config, coarse
    )
    assert coarse["status"] == "matched"
    limit = preparation.force_evaluations + 1
    limited = TwoSaddles(potential, limit=limit)
    report = {"branches": []}
    AseNebBackend(limited)._check_connectivity(*inputs(limited), config, report)
    assert limited.force_evaluations == limit
    assert len(report["branches"]) == 2
    assert report["status"] == "failed"
    assert "ForceBudgetExceeded" in report["message"]
    assert report["resolution_check"]["status"] == "sensitive_or_unresolved"
