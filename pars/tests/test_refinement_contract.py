"""Refinement contract: candidate -> resampled chain -> MLP NEB -> RefinedPath
(design section 43, MVP-5), on an analytic potential, without prrs."""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest
from ase.calculators.calculator import Calculator, all_changes

from pars.config import ParsConfig
from pars.refinement.ase_neb import (
    AseNebBackend,
    _climbing_stage_accepted,
    _plain_band_resolved,
)
from pars.schema import RefinementStatus

from conftest import (
    ExplodingCalculator,
    NarrowBarrierCalculator,
    ProtonHopCalculator,
    TripleWellCalculator,
    UncertaintyCalculator,
    make_candidate,
)

D = 1.6

# the single-atom potentials here are external fields: rigid motion changes the
# energy, so removing it (the default, for isolated molecules) would be wrong
EXTERNAL_FIELD = ParsConfig(refinement_remove_rigid_motion=False)


def _well_candidate() -> object:
    from ase import Atoms
    from ase.optimize import LBFGS

    frames = 7
    calculator = TripleWellCalculator()

    def minimum(x: float, y: float) -> np.ndarray:
        atoms = Atoms(numbers=[1], positions=[[x, y, 0.0]])
        atoms.calc = calculator
        LBFGS(atoms, logfile=None).run(fmax=1e-6, steps=300)
        return atoms.positions.copy()

    # endpoints are the true minima, so the band is consistent from step zero
    # (mirrors a PRRS candidate whose ends already sit in basins)
    start, end = minimum(0.0, 0.0), minimum(D, 0.0)
    t = np.linspace(0.0, 1.0, frames)
    positions = (1.0 - t)[:, None, None] * start[None] + t[:, None, None] * end[None]
    positions[1:-1, 0, 2] += 0.02 * np.sin(np.pi * t[1:-1])
    energies = []
    for frame_positions in positions:
        atoms = Atoms(numbers=[1], positions=frame_positions)
        atoms.calc = calculator
        energies.append(atoms.get_potential_energy())
    return make_candidate(
        "pc-well-0",
        "rc-well",
        energies,
        [frozenset({(0, 1)})] * 2 + [frozenset()] * (frames - 2),
        positions=positions,
    )


@pytest.mark.parametrize("point", [(0.4, 0.3, 0.0), (1.1, -0.2, 0.1), (0.2, 1.3, -0.3)])
def test_triple_well_forces_are_minus_energy_gradient(point):
    # every NEB assertion below is meaningless if forces and energy disagree
    from ase import Atoms

    def energy(position):
        atoms = Atoms(numbers=[1], positions=[position])
        atoms.calc = TripleWellCalculator()
        return atoms.get_potential_energy()

    position = np.asarray(point)
    step = 1e-5
    finite_difference = -np.asarray(
        [
            (energy(position + step * axis) - energy(position - step * axis)) / (2 * step)
            for axis in np.eye(3)
        ]
    )
    atoms = Atoms(numbers=[1], positions=[position])
    atoms.calc = TripleWellCalculator()
    assert atoms.get_forces()[0] == pytest.approx(finite_difference, abs=1e-6)


def test_refined_path_contract():
    backend = AseNebBackend(TripleWellCalculator())
    config = dataclasses.replace(
        EXTERNAL_FIELD,
        refinement_images=7,
        refinement_fmax_neb=0.05,
        refinement_max_steps_neb=500,
    )
    refined = backend.refine(_well_candidate(), config)

    assert refined.status in (
        RefinementStatus.NEB_CONVERGED,
        RefinementStatus.CI_NEB_CONVERGED,
    )
    assert refined.positions.shape == (7, 1, 3)
    assert refined.energies.shape == (7,)
    assert refined.barrier_eV > 0.0  # the triple well separates its minima
    assert refined.neb_residual_eV_A <= config.refinement_fmax_neb + 1e-6
    assert refined.uncertainty_max_eV_A is None  # potential reports no uncertainty
    assert refined.climbing_image is True
    assert refined.provenance["candidate_id"] == "pc-well-0"


TRIPLE_WELL_BARRIER_EV = 2.11548  # E(saddle) - E(c0 minimum), from a root solve
# of the analytic gradient at (0.774, 0.0005); one negative Hessian eigenvalue


def test_climbing_image_reaches_the_saddle():
    # plain NEB images straddle the saddle and underestimate it (~2.100 eV
    # here); a correct climbing image RAISES the barrier onto the saddle and
    # must be accepted, not mistaken for a runaway
    backend = AseNebBackend(TripleWellCalculator())
    config = dataclasses.replace(
        EXTERNAL_FIELD, refinement_images=7, refinement_max_steps_neb=1000
    )
    refined = backend.refine(_well_candidate(), config)
    assert refined.status is RefinementStatus.CI_NEB_CONVERGED, refined.message
    assert refined.barrier_eV == pytest.approx(TRIPLE_WELL_BARRIER_EV, abs=0.005)
    assert refined.provenance["saddle_check"]["negative_modes"] == 1
    assert refined.provenance["saddle_check"]["connectivity"] == "matched"
    connection = refined.provenance["connectivity_check"]
    assert connection["endpoint_identity_scope"] == "refined_endpoints_only"
    assert {b["matched_endpoint"] for b in connection["branches"]} == {"reactant", "product"}
    for branch in connection["branches"]:
        assert np.all(np.diff([connection["saddle_energy_eV"], *branch["energies_eV"]]) < 0)
        assert branch["minimum_check"]["negative_modes"] == 0


class DegenerateSaddleCalculator(Calculator):
    implemented_properties = ["energy", "forces"]

    def __init__(self, flat=False):
        super().__init__()
        self.flat = flat

    def calculate(
        self, atoms=None, properties=("energy", "forces"), system_changes=all_changes
    ):
        super().calculate(atoms, properties, system_changes)
        x, y, z = self.atoms.positions[0]
        if self.flat:
            energy = (x**4 - 1) ** 2 + 0.5 * (y * y + z * z)
            forces = [-8 * x**3 * (x**4 - 1), -y, -z]
        else:
            energy = (x * x - 1) ** 2 + (y * y - 1) ** 2 + 0.5 * z * z
            forces = [-4 * x * (x * x - 1), -4 * y * (y * y - 1), -z]
        self.results = {"energy": float(energy), "forces": np.asarray([forces])}


class ThreeMinimaCalculator(Calculator):
    """Three minima at x=-1,0,+1 and two index-one saddles at +/-1/sqrt(3)."""

    implemented_properties = ["energy", "forces"]

    def calculate(
        self, atoms=None, properties=("energy", "forces"), system_changes=all_changes
    ):
        super().calculate(atoms, properties, system_changes)
        x, y, z = self.atoms.positions[0]
        self.results = {
            "energy": float(x * x * (x * x - 1) ** 2 + 0.5 * (y * y + z * z)),
            "forces": np.array([[-2 * x * (x * x - 1) * (3 * x * x - 1), -y, -z]]),
        }


def test_index_one_saddle_connecting_an_intermediate_basin_is_refused():
    candidate = _line_candidate((-1.0, 0.0, 0.0), (1.0, 0.0, 0.0), frames=11)
    config = dataclasses.replace(
        EXTERNAL_FIELD, refinement_images=11, refinement_fmax_neb=0.001
    )
    refined = AseNebBackend(ThreeMinimaCalculator()).refine(candidate, config)
    assert refined.status is RefinementStatus.FAILED, refined.message
    assert refined.barrier_eV is None
    assert refined.provenance["saddle_check"]["negative_modes"] == 1
    connection = refined.provenance["connectivity_check"]
    assert len(connection["branches"]) == 2
    assert all(b["minimum_check"]["negative_modes"] == 0 for b in connection["branches"])
    landed = [b["positions_A"][0][0] for b in connection["branches"]]
    assert min(abs(x) for x in landed) < 0.01  # the unrequested middle minimum
    assert "endpoint basins" in refined.message


def test_unfinished_descent_refuses_ci_success_and_preserves_partial_report():
    config = dataclasses.replace(EXTERNAL_FIELD, refinement_descent_max_steps=0)
    refined = AseNebBackend(TripleWellCalculator()).refine(_well_candidate(), config)
    assert refined.status is RefinementStatus.FAILED, refined.message
    assert refined.barrier_eV is None
    assert "descent accepted-step budget exhausted" in refined.message
    branch = refined.provenance["connectivity_check"]["branches"][0]
    assert len(branch["energies_eV"]) == 1
    assert branch["accepted_steps"] == 0
    assert refined.positions is not None


@pytest.mark.parametrize("flat,expected_modes", [(False, 2), (True, 0)])
def test_stationary_climbing_image_requires_one_resolved_negative_mode(flat, expected_modes):
    t = np.linspace(-1.0, 1.0, 9)
    positions = np.zeros((9, 1, 3))
    positions[:, 0, 0] = t
    if not flat:
        positions[:, 0, 1] = t
    candidate = make_candidate(
        "degenerate", "ch", np.zeros(9), [frozenset()] * 9, positions=positions
    )
    config = dataclasses.replace(EXTERNAL_FIELD, refinement_fmax_neb=0.001)
    refined = AseNebBackend(DegenerateSaddleCalculator(flat)).refine(candidate, config)
    assert refined.status is RefinementStatus.FAILED, refined.message
    assert refined.barrier_eV is None
    checked = refined.provenance["saddle_check"]
    assert checked["negative_modes"] == expected_modes
    assert "index-one saddle" in refined.message
    if not flat:
        assert checked["eigenvalues_eV_A2"] == pytest.approx([-4.0, -4.0, 1.0], abs=1e-4)
        assert refined.energies.max() == pytest.approx(2.0)
        # An explicit path through the mixed minimum (1,-1) has maximum 1 eV.
        # The stationary 2 eV diagonal is therefore a false useful barrier.
        detour = np.concatenate(
            [np.column_stack([t, -np.ones(9)]), np.column_stack([np.ones(9), t])]
        )
        energies = (detour[:, 0] ** 2 - 1) ** 2 + (detour[:, 1] ** 2 - 1) ** 2
        assert energies.max() == pytest.approx(1.0)


def _straight_band(images=7):
    band = np.zeros((images, 1, 3))
    band[:, 0, 0] = np.linspace(0.0, 1.5, images)
    return band


def test_climbing_stage_accepted_when_top_stays_bracketed():
    plain = _straight_band()
    ci = plain.copy()
    ci[3, 0, 0] += 0.1  # within the 0.25 A spacing to its plain neighbours
    energies = np.asarray([0.0, 1.0, 1.8, 2.1, 1.8, 1.0, 0.2])
    assert _climbing_stage_accepted(plain, ci, energies, converged=True)


def test_climbing_stage_rejected_when_top_runs_away():
    plain = _straight_band()
    ci = plain.copy()
    ci[3, 0, 1] += 2.0  # climbed far off the band, e.g. up a repulsive wall
    energies = np.asarray([0.0, 1.0, 1.8, 9.0, 1.8, 1.0, 0.2])
    assert not _climbing_stage_accepted(plain, ci, energies, converged=True)


def test_climbing_stage_rejected_when_unconverged_or_top_at_endpoint():
    plain = _straight_band()
    interior = np.asarray([0.0, 1.0, 1.8, 2.1, 1.8, 1.0, 0.2])
    assert not _climbing_stage_accepted(plain, plain, interior, converged=False)
    endpoint_top = np.asarray([0.0, 1.0, 1.8, 2.1, 1.8, 1.0, 3.0])
    assert not _climbing_stage_accepted(plain, plain, endpoint_top, converged=True)


def test_uncertainty_gate_reports_unresolved_not_success():
    backend = AseNebBackend(UncertaintyCalculator(TripleWellCalculator(), uncertainty=0.5))
    config = dataclasses.replace(
        EXTERNAL_FIELD,
        refinement_images=7,
        refinement_max_steps_neb=500,
        refinement_uncertainty_refusal=0.1,
    )
    refined = backend.refine(_well_candidate(), config)
    assert refined.status is RefinementStatus.UNRESOLVED_MLP
    assert refined.uncertainty_max_eV_A == pytest.approx(0.5)
    assert "not reporting success" in refined.message


def test_no_uncertainty_gate_without_threshold():
    backend = AseNebBackend(UncertaintyCalculator(TripleWellCalculator(), uncertainty=0.5))
    config = dataclasses.replace(
        EXTERNAL_FIELD, refinement_images=7, refinement_max_steps_neb=500
    )
    refined = backend.refine(_well_candidate(), config)
    assert refined.status in (
        RefinementStatus.NEB_CONVERGED,
        RefinementStatus.CI_NEB_CONVERGED,
    )


def test_backend_failure_is_reported_not_raised():
    backend = AseNebBackend(ExplodingCalculator())
    refined = backend.refine(_well_candidate(), EXTERNAL_FIELD)
    assert refined.status is RefinementStatus.FAILED
    assert "synthetic calculator failure" in refined.message
    assert refined.positions is None


def test_refined_path_summary_is_json_safe():
    backend = AseNebBackend(TripleWellCalculator())
    config = dataclasses.replace(
        EXTERNAL_FIELD, refinement_images=7, refinement_max_steps_neb=500
    )
    refined = backend.refine(_well_candidate(), config)
    summary = refined.to_summary()
    import json

    json.dumps(summary, allow_nan=False)
    assert summary["status"] == refined.status.value


def _line_candidate(start, end, frames=7):
    t = np.linspace(0.0, 1.0, frames)
    start, end = np.asarray(start, dtype=float), np.asarray(end, dtype=float)
    positions = ((1.0 - t)[:, None] * start + t[:, None] * end)[:, None, :]
    return make_candidate(
        "pc-line", "rc-well", np.zeros(frames), [frozenset()] * frames, positions=positions
    )


def test_endpoints_relaxing_into_one_basin_fail_instead_of_converging():
    # the last frame (0.35, 0) sits inside the c0 basin (saddle at x~0.774): it
    # relaxes back onto the reactant and a "converged" A->A band means nothing
    backend = AseNebBackend(TripleWellCalculator())
    refined = backend.refine(_line_candidate((0.0, 0.0, 0.0), (0.35, 0.0, 0.0)), EXTERNAL_FIELD)
    assert refined.status is RefinementStatus.FAILED
    assert "same minimum" in refined.message


def test_unconverged_endpoint_relaxation_fails():
    backend = AseNebBackend(TripleWellCalculator())
    config = dataclasses.replace(EXTERNAL_FIELD, refinement_max_steps_endpoints=1)
    refined = backend.refine(_line_candidate((0.25, 0.15, 0.0), (1.3, 0.1, 0.0)), config)
    assert refined.status is RefinementStatus.FAILED
    assert "endpoint relaxation" in refined.message


def test_refined_path_keeps_the_candidate_provenance():
    backend = AseNebBackend(TripleWellCalculator())
    config = dataclasses.replace(
        EXTERNAL_FIELD, refinement_images=7, refinement_max_steps_neb=500
    )
    refined = backend.refine(_well_candidate(), config)
    assert refined.provenance["source"]["attempt_id"] == "t-pc-well-0"
    assert refined.provenance["source"]["channel_id"] == "rc-well"


NARROW_BARRIER_EV = 0.41947  # dense 1D scan of NarrowBarrierCalculator


def _rattling_crossing_candidate():
    # like an MD response: most of the arc length is rattling (in y) inside the
    # basins and the crossing itself is short, so the uniform arc-length seed
    # puts no image near the saddle at x=0 (seed x: -1.03 ... -0.45, 0.51 ...)
    xa = np.linspace(-1.0342, -0.5, 9)
    xb = np.linspace(0.6, 1.0342, 7)
    ya = 0.25 * (-1.0) ** np.arange(9)
    yb = 0.25 * (-1.0) ** np.arange(7)
    ya[0] = yb[-1] = 0.0
    positions = np.zeros((16, 1, 3))
    positions[:, 0, 0] = np.concatenate([xa, xb])
    positions[:, 0, 1] = np.concatenate([ya, yb])
    return make_candidate(
        "pc-rattle", "rc-narrow", np.zeros(16), [frozenset()] * 16, positions=positions
    )


def test_weak_springs_leave_an_unresolved_band_that_is_not_reported_as_converged():
    # k=0.1 (the ASE default) "converges" with the saddle between two images and
    # a barrier of 0.29 eV for a 0.42 eV saddle: on the p1 runs this was 22/27
    config = dataclasses.replace(
        EXTERNAL_FIELD, refinement_spring_k=0.1, refinement_climbing_image=False
    )
    refined = AseNebBackend(NarrowBarrierCalculator()).refine(
        _rattling_crossing_candidate(), config
    )
    assert refined.status is RefinementStatus.FAILED
    assert refined.barrier_eV is None  # an unresolved band has no barrier to report
    assert "unresolved" in refined.message
    assert refined.positions is not None  # the band itself is kept for inspection


def test_default_springs_resolve_the_saddle_without_climbing():
    config = dataclasses.replace(EXTERNAL_FIELD, refinement_climbing_image=False)
    refined = AseNebBackend(NarrowBarrierCalculator()).refine(
        _rattling_crossing_candidate(), config
    )
    assert refined.status is RefinementStatus.NEB_CONVERGED, refined.message
    assert refined.provenance["connectivity_check"] == {
        "status": "not_checked",
        "reason": "plain_neb_band_estimate",
    }
    # plain NEB brackets the saddle from below
    assert NARROW_BARRIER_EV - 0.01 <= refined.barrier_eV <= NARROW_BARRIER_EV + 1e-4


def _band_from_gaps(gaps):
    positions = np.zeros((len(gaps) + 1, 1, 3))
    positions[:, 0, 0] = np.concatenate([[0.0], np.cumsum(gaps)])
    return positions


@pytest.mark.parametrize(
    "gaps, energies",
    [
        # p1_round7 t000029 and t000002 after plain NEB (image gaps in A, E - E0 in eV):
        # the top image sits next to the 0.5-0.7 A jump across the proton transfer
        (
            [0.092, 0.098, 0.145, 0.227, 0.415, 0.699, 0.516, 0.336],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.0, 0.0],
        ),
        (
            [0.736, 0.697, 0.463, 0.161, 0.158, 0.027, 0.114, 0.105],
            [0.0, 0.005, 0.001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ),
    ],
)
def test_real_flat_bands_are_unresolved(gaps, energies):
    assert not _plain_band_resolved(_band_from_gaps(gaps), np.asarray(energies), 2.0)


def test_uniform_band_is_resolved():
    energies = np.asarray([0.0, 0.1, 0.25, 0.38, 0.42, 0.37, 0.24, 0.1, 0.0])
    assert _plain_band_resolved(_band_from_gaps([0.25] * 8), energies, 2.0)


def test_coarse_ends_do_not_matter_when_the_top_is_finely_resolved():
    # relaxed endpoints far from their raw frames leave big end gaps; the
    # barrier is still resolved because the top image's neighbours are close
    gaps = [0.6, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.6]
    energies = np.asarray([0.0, 0.05, 0.2, 0.35, 0.42, 0.35, 0.2, 0.05, 0.0])
    assert _plain_band_resolved(_band_from_gaps(gaps), energies, 2.0)


def _tumbling_proton_hop_candidate():
    # the p1 t000002 failure: MD frames drift and tumble mid-trajectory (out and
    # back, so the endpoints nearly coincide in orientation); in Cartesian
    # coordinates NEB then spends the band on rigid motion and the proton hops
    # back and forth between images that all sit in a basin
    from scipy.spatial.transform import Rotation

    axis = np.asarray([0.3, 1.0, 0.5]) / np.linalg.norm([0.3, 1.0, 0.5])
    frames = []
    for index, q in enumerate(np.linspace(1.0342, -1.0342, 16)):
        out_and_back = np.sin(np.pi * index / 15)
        r_ah = (2.5 + q) / 2.0
        molecule = np.asarray([[-1.25, 0.0, 0.0], [1.25, 0.0, 0.0], [r_ah - 1.25, 0.0, 0.0]])
        rotation = Rotation.from_rotvec(axis * np.radians(120.0) * out_and_back).as_matrix()
        frames.append(molecule @ rotation.T + np.asarray([1.0, -0.6, 0.3]) * 1.5 * out_and_back)
    return make_candidate(
        "pc-hop", "rc-hop", np.zeros(16), [frozenset()] * 16, positions=np.asarray(frames)
    )


def test_rigid_motion_in_the_seed_does_not_fake_a_saddle():
    refined = AseNebBackend(ProtonHopCalculator()).refine(
        _tumbling_proton_hop_candidate(), ParsConfig()
    )
    assert refined.status is RefinementStatus.CI_NEB_CONVERGED, refined.message
    assert refined.barrier_eV == pytest.approx(NARROW_BARRIER_EV, abs=0.003)
