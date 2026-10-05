"""Independent downhill connectivity contracts on explicit analytic PESs."""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes

from pars.config import ParsConfig
from pars.refinement.ase_neb import AseNebBackend
from pars.refinement.descent import descend
from pars.schema import RefinementStatus

from test_refinement_contract import ThreeMinimaCalculator
from conftest import make_candidate

CONFIG = ParsConfig(refinement_remove_rigid_motion=False)


def atoms_at(position, calculator):
    atoms = Atoms("H", positions=[position])
    atoms.calc = calculator
    return atoms


@pytest.mark.parametrize("point", [(-0.8, 0.2, 0.1), (0.3, -0.1, 0.2)])
def test_three_minima_forces_match_independent_energy_differences(point):
    atoms = atoms_at(point, ThreeMinimaCalculator())
    force = atoms.get_forces()[0].copy()
    step = 1e-6
    for column in range(3):
        delta = step * np.eye(3)[column]
        plus = atoms_at(np.array(point) + delta, ThreeMinimaCalculator()).get_potential_energy()
        minus = atoms_at(
            np.array(point) - delta, ThreeMinimaCalculator()
        ).get_potential_energy()
        assert force[column] == pytest.approx(-(plus - minus) / (2 * step), abs=1e-8)


def check_three_minima(mode=1.0, config=CONFIG):
    calculator = ThreeMinimaCalculator()
    saddle = atoms_at((-1 / np.sqrt(3), 0, 0), calculator)
    endpoints = [atoms_at((-1, 0, 0), calculator), atoms_at((0, 0, 0), calculator)]
    report = {"branches": []}
    AseNebBackend(calculator)._check_connectivity(
        saddle, [[mode, 0, 0]], endpoints, config, report
    )
    return report


def test_negative_mode_sign_changes_labels_only():
    forward, reverse = check_three_minima(), check_three_minima(-1.0)
    assert forward["status"] == reverse["status"] == "matched"
    labels = [b["matched_endpoint"] for b in forward["branches"]]
    reversed_labels = [b["matched_endpoint"] for b in reverse["branches"]]
    assert labels == list(reversed(reversed_labels))
    assert forward["resolution_check"]["status"] == "stable_on_tested_scales"
    assert forward["continuous_connection_certified"] is False


def test_ambiguous_endpoint_pairing_is_refused():
    report = check_three_minima(
        config=dataclasses.replace(CONFIG, refinement_endpoint_identity_rmsd=2.0)
    )
    assert all(b["minimum_check"]["negative_modes"] == 0 for b in report["branches"])
    assert report["status"] == "failed"
    assert "uniquely" in report["message"]


def test_small_force_near_saddle_does_not_establish_a_minimum():
    report = check_three_minima(
        config=dataclasses.replace(CONFIG, refinement_descent_displacement=1e-4)
    )
    branch = report["branches"][0]
    assert branch["max_force_eV_A"] < CONFIG.refinement_descent_fmax
    assert branch["minimum_check"]["negative_modes"] == 1
    assert report["status"] == "failed"


class TiltedRing(Calculator):
    """One minimum and an index-one saddle; both signed descents meet below."""

    implemented_properties = ["energy", "forces"]

    def calculate(
        self, atoms=None, properties=("energy", "forces"), system_changes=all_changes
    ):
        super().calculate(atoms, properties, system_changes)
        x, y, z = self.atoms.positions[0]
        r = x * x + y * y - 1
        self.results = {
            "energy": float(r * r + 0.5 * y + 0.5 * z * z),
            "forces": np.array([[-4 * x * r, -4 * y * r - 0.5, -z]]),
        }


def test_both_signed_descents_to_one_basin_are_refused():
    roots = sorted(np.roots([4.0, 0.0, -4.0, 0.5]))
    calculator = TiltedRing()
    saddle = atoms_at((0, roots[-1], 0), calculator)
    endpoints = [atoms_at((0, roots[0], 0), calculator), atoms_at((-1, 0, 0), calculator)]
    report = {"branches": []}
    AseNebBackend(calculator)._check_connectivity(
        saddle, [[1, 0, 0]], endpoints, CONFIG, report
    )
    assert len(report["branches"]) == 2
    assert all(b["status"] == "minimum_checked" for b in report["branches"])
    assert all(b["endpoint_rmsd_A"]["reactant"] < 0.01 for b in report["branches"])
    assert report["status"] == "failed"


def test_energy_backtracking_rejects_inconsistent_uphill_forces():
    class Uphill(Calculator):
        implemented_properties = ["energy", "forces"]

        def calculate(
            self, atoms=None, properties=("energy", "forces"), system_changes=all_changes
        ):
            super().calculate(atoms, properties, system_changes)
            self.results = {
                "energy": float(self.atoms.positions[0, 0] ** 2),
                "forces": np.array([[2.0, 0.0, 0.0]]),
            }

    atoms = atoms_at((1, 0, 0), Uphill())
    report = {}
    config = dataclasses.replace(CONFIG, refinement_descent_max_backtracks=2)
    assert not descend(atoms, 2.0, config, report)
    assert report["accepted_steps"] == 0
    assert report["backtracks"] == 3
    assert np.array_equal(atoms.positions, [[1, 0, 0]])
    assert "backtracking budget" in report["message"]


@pytest.mark.parametrize("uncertainty", [None, 0.5])
def test_descent_uncertainty_gate_applies_outside_the_neb_band(uncertainty):
    class Uncertain(ThreeMinimaCalculator):
        def calculate(self, *args, **kwargs):
            super().calculate(*args, **kwargs)
            if uncertainty is not None:
                self.results["force_uncertainty_eV_A"] = uncertainty

    atoms = atoms_at((-0.6, 0, 0), Uncertain())
    report = {}
    config = dataclasses.replace(CONFIG, refinement_uncertainty_refusal=0.1)
    assert not descend(atoms, 4 / 27, config, report)
    assert report["status"] == "unresolved_mlp"


@pytest.mark.parametrize(
    "field,value",
    [
        ("refinement_descent_displacement", 0.0),
        ("refinement_descent_step", float("nan")),
        ("refinement_descent_max_steps", 1.5),
        ("refinement_descent_max_backtracks", True),
    ],
)
def test_invalid_descent_policy_is_refused(field, value):
    with pytest.raises(ValueError, match=field):
        ParsConfig(**{field: value})


class NarrowIntermediate(Calculator):
    """A 0.002 A wide well beyond the saddle, hidden inside the coarse kick."""

    implemented_properties = ["energy", "forces"]

    def calculate(
        self, atoms=None, properties=("energy", "forces"), system_changes=all_changes
    ):
        super().calculate(atoms, properties, system_changes)
        x, y, z = self.atoms.positions[0]
        u = (x - 0.015) / 0.002
        gaussian = np.exp(-u * u)
        self.results = {
            "energy": float((x * x - 1) ** 2 + 0.5 * (y * y + z * z) - 0.02 * gaussian),
            "forces": np.array([[-4 * x * (x * x - 1) - 20 * u * gaussian, -y, -z]]),
        }


@pytest.mark.parametrize("x", [0.013, 0.015, 0.017, 0.020])
def test_narrow_intermediate_forces_match_energy_differences(x):
    step = 1e-7
    plus = atoms_at((x + step, 0, 0), NarrowIntermediate()).get_potential_energy()
    minus = atoms_at((x - step, 0, 0), NarrowIntermediate()).get_potential_energy()
    force = atoms_at((x, 0, 0), NarrowIntermediate()).get_forces()[0, 0]
    assert force == pytest.approx(-(plus - minus) / (2 * step), abs=1e-6)


def test_smaller_resolution_detects_a_basin_skipped_by_the_coarse_kick():
    positions = np.zeros((9, 1, 3))
    positions[:, 0, 0] = np.linspace(-1, 1, 9)
    candidate = make_candidate(
        "narrow", "channel", np.zeros(9), [frozenset()] * 9, positions=positions
    )
    result = AseNebBackend(NarrowIntermediate()).refine(candidate, CONFIG)
    assert result.status is RefinementStatus.FAILED, result.message
    assert result.barrier_eV is None
    report = result.provenance["connectivity_check"]
    assert {b["matched_endpoint"] for b in report["branches"]} == {"reactant", "product"}
    fine = report["resolution_check"]["fine"]
    assert all(b["minimum_check"]["negative_modes"] == 0 for b in fine["branches"])
    assert min(abs(b["positions_A"][0][0] - 0.01500600110285) for b in fine["branches"]) < 1e-6
    assert report["resolution_check"]["status"] == "sensitive_or_unresolved"
    assert "finer descent resolution" in result.message


@pytest.mark.parametrize("factor", [0.0, 1.0, -0.25, float("nan")])
def test_resolution_policy_cannot_disable_finer_replay(factor):
    with pytest.raises(ValueError, match="refinement_descent_resolution_factor"):
        ParsConfig(refinement_descent_resolution_factor=factor)


class SpectatorThreeMinima(ThreeMinimaCalculator):
    """One reactive coordinate, with independent harmonic spectator atoms."""

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        Calculator.calculate(self, atoms, properties, system_changes)
        positions = self.atoms.positions
        fixed = np.zeros_like(positions)
        fixed[:, 0] = 3 * np.arange(len(positions))
        d = positions - fixed
        x = d[0, 0]
        forces = -d.copy()
        forces[0, 0] = -2 * x * (x * x - 1) * (3 * x * x - 1)
        energy = x * x * (x * x - 1) ** 2 + 0.5 * (np.sum(d[1:] ** 2) + np.sum(d[0, 1:] ** 2))
        self.results = {"energy": float(energy), "forces": forces}


def spectator_candidate(n, target=1):
    positions = np.zeros((9, n, 3))
    positions[:, :, 0] = 3 * np.arange(n)
    positions[:, 0, 0] = np.linspace(-1, target, 9)
    candidate = make_candidate(
        "spectators", "channel", np.zeros(9), [frozenset()] * 9, positions=positions
    )
    candidate.numbers = np.ones(n, dtype=int)
    return candidate


@pytest.mark.parametrize("n", [1, 500])
def test_spectators_cannot_turn_third_basin_into_a_valid_endpoint(n):
    candidate = spectator_candidate(n)
    result = AseNebBackend(SpectatorThreeMinima()).refine(
        candidate, dataclasses.replace(CONFIG, refinement_fmax_neb=0.001)
    )
    assert result.status is RefinementStatus.FAILED, result.message
    assert result.barrier_eV is None
    report = result.provenance["connectivity_check"]
    assert report["status"] == "failed"
    middle = min(report["branches"], key=lambda b: abs(b["positions_A"][0][0]))
    assert abs(middle["positions_A"][0][0]) < 0.001
    assert min(middle["endpoint_max_displacement_A"].values()) > 0.99


def test_spectators_do_not_collapse_distinct_valid_endpoints():
    # Whole-system RMSD is below 0.05 A, but the active atom changes by 1 A.
    result = AseNebBackend(SpectatorThreeMinima()).refine(
        spectator_candidate(500, target=0),
        dataclasses.replace(CONFIG, refinement_fmax_neb=0.001),
    )
    assert result.status is RefinementStatus.CI_NEB_CONVERGED, result.message
    assert result.provenance["connectivity_check"]["status"] == "matched"


def test_source_basin_witness_uses_the_same_local_identity_guard():
    candidate = spectator_candidate(500)
    wrong = candidate.positions[-1].copy()
    wrong[0, 0] = 0
    candidate.provenance = dataclasses.replace(
        candidate.provenance, extra={"endpoint_basin_references": {"product": wrong.tolist()}}
    )
    result = AseNebBackend(SpectatorThreeMinima()).refine(candidate, CONFIG)
    assert result.status is RefinementStatus.FAILED
    assert "basin identity mismatch" in result.message
    checked = result.provenance["endpoint_checks"]["product"]
    assert checked["basin_reference_rmsd_A"] < 0.05
    assert checked["basin_reference_max_displacement_A"] == pytest.approx(1)
