"""Physical and provenance regressions from the 2026-10-03 audit."""

import dataclasses
import json

import numpy as np
import pytest
from ase.calculators.calculator import Calculator, all_changes

from conftest import make_candidate
from pars.config import ParsConfig
from pars.refinement.ase_neb import AseNebBackend, _plain_band_resolved
from pars.schema import RefinementStatus
from test_cli import A, FACTORY_MODULE, _refine
from test_cli import factory_on_path as factory_on_path
from test_cli import selection_dir as selection_dir


class QuarticWell(Calculator):
    implemented_properties = ["energy", "forces"]

    def __init__(self, uncertainty=None):
        super().__init__()
        self.uncertainty = uncertainty

    def calculate(
        self, atoms=None, properties=("energy", "forces"), system_changes=all_changes
    ):
        super().calculate(atoms, properties, system_changes)
        x, y, z = atoms.positions[0]
        self.results = {
            "energy": 0.25 * x**2 * (x - 2) ** 2 + 0.5 * (y * y + z * z),
            "forces": np.array([[-x * (x - 2) * (x - 1), -y, -z]]),
        }
        if self.uncertainty is not None:
            self.results["force_uncertainty_eV_A"] = self.uncertainty


def quartic_candidate(start=0.0):
    positions = np.zeros((9, 1, 3))
    positions[:, 0, 0] = np.linspace(start, 2.0, 9)
    return make_candidate(
        "quartic", "channel", np.zeros(9), [frozenset()] * 9, positions=positions
    )


CONFIG = ParsConfig(refinement_remove_rigid_motion=False)


def test_zero_force_saddle_is_not_an_endpoint_minimum():
    result = AseNebBackend(QuarticWell()).refine(quartic_candidate(1.0), CONFIG)
    assert result.status is RefinementStatus.FAILED
    assert result.barrier_eV is None
    assert "negative curvature" in result.message


def test_endpoint_basin_witnesses_match_and_survive_refinement():
    candidate = quartic_candidate()
    candidate.provenance = dataclasses.replace(
        candidate.provenance,
        extra={"endpoint_basin_references": {"reactant": [[0, 0, 0]], "product": [[2, 0, 0]]}},
    )
    result = AseNebBackend(QuarticWell()).refine(candidate, CONFIG)
    assert result.status is RefinementStatus.CI_NEB_CONVERGED
    assert result.barrier_eV == pytest.approx(0.25)
    assert all(
        c["basin_identity"] == "matched" for c in result.provenance["endpoint_checks"].values()
    )
    assert result.provenance["connectivity_check"]["status"] == "matched"
    assert (
        result.provenance["connectivity_check"]["endpoint_identity_scope"]
        == "source_basin_witnesses"
    )


def test_true_minimum_in_the_wrong_endpoint_basin_is_refused():
    candidate = quartic_candidate()
    candidate.provenance = dataclasses.replace(
        candidate.provenance,
        extra={"endpoint_basin_references": {"reactant": [[0, 0, 0]], "product": [[0, 0, 0]]}},
    )
    result = AseNebBackend(QuarticWell()).refine(candidate, CONFIG)
    assert result.status is RefinementStatus.FAILED
    assert "basin identity mismatch" in result.message
    assert result.barrier_eV is None


def test_saddle_basin_witness_is_refused():
    candidate = quartic_candidate()
    candidate.provenance = dataclasses.replace(
        candidate.provenance,
        extra={"endpoint_basin_references": {"reactant": [[0, 0, 0]], "product": [[1, 0, 0]]}},
    )
    result = AseNebBackend(QuarticWell()).refine(candidate, CONFIG)
    assert result.status is RefinementStatus.FAILED
    assert "basin reference has negative curvature" in result.message


@pytest.mark.parametrize("value", [float("nan"), float("inf"), -float("inf")])
def test_nonfinite_refinement_threshold_is_rejected(value):
    with pytest.raises(ValueError, match="finite"):
        dataclasses.replace(CONFIG, refinement_uncertainty_refusal=value)


def test_plain_band_with_highest_energy_at_endpoint_is_unresolved():
    positions = quartic_candidate().positions
    assert not _plain_band_resolved(positions, np.linspace(1.0, 0.0, 9), 2.0)


@pytest.mark.parametrize("uncertainty", [float("nan"), float("inf"), -0.1])
def test_invalid_uncertainty_never_reports_success(uncertainty):
    result = AseNebBackend(QuarticWell(uncertainty)).refine(
        quartic_candidate(), dataclasses.replace(CONFIG, refinement_uncertainty_refusal=0.01)
    )
    assert result.status is RefinementStatus.FAILED
    assert "uncertainty" in result.message
    json.dumps(result.to_summary(), allow_nan=False)


def test_explicit_uncertainty_gate_requires_measurements():
    result = AseNebBackend(QuarticWell()).refine(
        quartic_candidate(), dataclasses.replace(CONFIG, refinement_uncertainty_refusal=0.01)
    )
    assert result.status is RefinementStatus.UNRESOLVED_MLP
    assert result.barrier_eV is None
    # the first unmeasured evaluation is an endpoint Hessian sample
    assert "uncertainty missing" in result.message


def test_partial_uncertainty_is_not_full_path_coverage():
    class Partial(QuarticWell):
        def calculate(self, *args, **kwargs):
            super().calculate(*args, **kwargs)
            # only the x=0.25 band image; Hessian samples sit within 0.001 of 0, 1, 2
            if 0.1 < self.atoms.positions[0, 0] < 0.5:
                self.results.pop("force_uncertainty_eV_A", None)

    result = AseNebBackend(Partial(0.005)).refine(
        quartic_candidate(), dataclasses.replace(CONFIG, refinement_uncertainty_refusal=0.01)
    )
    assert result.status is RefinementStatus.UNRESOLVED_MLP
    assert "partial" in result.message


def test_missing_descent_uncertainty_preserves_known_band_measurement():
    class BandOnly(QuarticWell):
        def calculate(self, *args, **kwargs):
            super().calculate(*args, **kwargs)
            # The stationary uniform band has x=0,.25,...,2 and Hessian samples
            # lie within 0.001 of it; the signed kick (0.05 away) lies outside
            # those measured points and must not inherit sigma.
            x = self.atoms.positions[0, 0]
            if abs(x - round(x * 4) / 4) > 0.002:
                self.results.pop("force_uncertainty_eV_A", None)

    result = AseNebBackend(BandOnly(0.005)).refine(
        quartic_candidate(), dataclasses.replace(CONFIG, refinement_uncertainty_refusal=0.01)
    )
    assert result.status is RefinementStatus.UNRESOLVED_MLP
    assert result.barrier_eV is None
    assert result.provenance["uncertainty_coverage"] == "complete"
    assert result.uncertainty_max_eV_A == pytest.approx(0.005)
    branch = result.provenance["connectivity_check"]["branches"][0]
    assert branch["uncertainty_coverage"] == "unavailable"


def test_missing_fine_kick_uncertainty_refuses_an_otherwise_matched_connection():
    class MissingFine(QuarticWell):
        def calculate(self, *args, **kwargs):
            super().calculate(*args, **kwargs)
            if abs(self.atoms.positions[0, 0] - 0.9875) < 1e-12:
                self.results.pop("force_uncertainty_eV_A", None)

    result = AseNebBackend(MissingFine(0.005)).refine(
        quartic_candidate(), dataclasses.replace(CONFIG, refinement_uncertainty_refusal=0.01)
    )
    assert result.status is RefinementStatus.UNRESOLVED_MLP
    assert result.barrier_eV is None
    report = result.provenance["connectivity_check"]
    assert len(report["branches"]) == 2
    assert report["resolution_check"]["fine"]["status"] == "unresolved_mlp"
    assert result.uncertainty_max_eV_A == pytest.approx(0.005)


def test_single_candidate_redo_groups_different_conditions(selection_dir, factory_on_path):
    assert _refine(selection_dir) == 0
    assert _refine(selection_dir, "--only", A, "--fmax-neb", "0.02", "--redo") == 0
    index = json.loads((selection_dir / "refined/refined_groups.json").read_text())
    assert index["homogeneous"] is False
    assert len(index["groups"]) == 2
    assert sorted(len(g["candidate_ids"]) for g in index["groups"]) == [1, 1]


def test_same_factory_replaced_checkpoint_refuses_resume(
    selection_dir, factory_on_path, tmp_path, monkeypatch
):
    import pars_test_factory

    checkpoint = tmp_path / "model.bin"
    checkpoint.write_bytes(b"original")
    monkeypatch.setattr(
        pars_test_factory,
        "calculator_metadata",
        lambda: {"model_paths": str(checkpoint)},
        raising=False,
    )
    assert _refine(selection_dir, "--only", A) == 0
    before = (selection_dir / "refined" / f"refined-{A}.json").read_bytes()
    checkpoint.write_bytes(b"replacement")
    with pytest.raises(ValueError, match="calculator"):
        _refine(selection_dir, "--only", A)
    assert (selection_dir / "refined" / f"refined-{A}.json").read_bytes() == before


def test_runtime_calculator_parameters_change_without_factory_rename(
    selection_dir, factory_on_path, monkeypatch
):
    import pars_test_factory

    original = pars_test_factory.make_calculator
    current_scale = {"scale": 1.0}

    def make():
        calculator = original()
        calculator.parameters["physical_scale"] = current_scale["scale"]
        return calculator

    monkeypatch.setattr(pars_test_factory, "make_calculator", make)
    assert _refine(selection_dir, "--only", A) == 0
    current_scale["scale"] = 2.0
    with pytest.raises(ValueError, match="calculator_identity"):
        _refine(selection_dir, "--only", A)


def test_cli_fixture_forces_match_energy_gradient():
    namespace = {}
    exec(FACTORY_MODULE, namespace)
    from ase import Atoms

    atoms = Atoms("HCO", positions=[[0.3, 0.2, 0.1], [1.8, 1.0, 0.2], [3.4, 0.3, 0.4]])
    atoms.calc = namespace["make_calculator"]()
    original = atoms.positions.copy()
    force = atoms.get_forces().copy()
    step = 1e-5
    for i in range(3):
        for j in range(3):
            offset = np.zeros_like(original)
            offset[i, j] = step
            atoms.positions = original + offset
            plus = atoms.get_potential_energy()
            atoms.positions = original - offset
            minus = atoms.get_potential_energy()
            assert force[i, j] == pytest.approx(-(plus - minus) / (2 * step), abs=1e-7)
