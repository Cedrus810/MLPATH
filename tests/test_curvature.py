"""Curvature as a measurement rather than a matrix: source, error floor, and scale.

Three things are pinned here. The six rigid-body modes have exactly zero curvature on any
invariant potential, so their measured size is the error floor of the whole calculation
and needs no reference data. An analytic second derivative from the backend must be
proven by execution, because the unit factors a calculator applies to energy and forces
are not necessarily applied to its Hessian and a wrong one is silent whenever the factors
happen to be 1. And a symmetric difference at amplitude a is not an approximation that
improves as a shrinks: it is the curvature of the potential low-passed at wavenumber
~2/a, which is a testable statement about the two available differences.
"""

import numpy as np
import pytest
from ase import Atoms
from ase.build import molecule
from ase.calculators.calculator import Calculator, all_changes
from prrs.calculators import DoubleWell, double_well_factory
from prrs.config import SearchConfig
from prrs.reliability import GateRejected, GuardedCalculator
from prrs.runner import (
    confirm_minimum,
    curvature_spectrum,
    floor_residual,
    hessian_spectrum,
    response_curvature,
    response_gradient,
)

CONFIG = dict(quench_fmax_eV_A=0.001)


def _pair(distance):
    atoms = Atoms("C2", positions=[[-distance / 2, 0, 0], [distance / 2, 0, 0]])
    atoms.calc = double_well_factory()
    return atoms


def _analytic_double_well_hessian(atoms):
    """d2E/dx2 of the double well in closed form, as a backend would supply it."""
    vector = atoms.positions[1] - atoms.positions[0]
    r = float(np.linalg.norm(vector))
    unit = vector / r
    x = (r - 1.8) / 0.6
    second = (1.6 / 0.36) * (3 * x * x - 1)
    first = 4 * 0.4 * x * (x * x - 1) / 0.6
    outer = np.outer(unit, unit)
    block = second * outer + (first / r) * (np.eye(3) - outer)
    hessian = np.zeros((6, 6))
    hessian[:3, :3] = hessian[3:, 3:] = block
    hessian[:3, 3:] = hessian[3:, :3] = -block
    return hessian


class WellWithHessian(DoubleWell):
    """The same physics, plus the analytic second derivative and a unit scale to abuse."""

    def __init__(self, energy_units_to_eV=1.0, length_units_to_A=1.0, corrupt=1.0):
        super().__init__()
        self.energy_units_to_eV = energy_units_to_eV
        self.length_units_to_A = length_units_to_A
        self._corrupt = corrupt

    def get_hessian(self, atoms):
        return self._corrupt * _analytic_double_well_hessian(atoms)


class TiltedWell(DoubleWell):
    """Double well plus a term in the absolute z coordinate: not translation invariant.

    Nothing about the forces looks wrong, and the potential is perfectly smooth. Only the
    rigid-body curvature betrays it, which is the point of measuring that floor.
    """

    def __init__(self, tilt=0.0):
        super().__init__()
        self.tilt = tilt

    def calculate(
        self, atoms=None, properties=("energy", "forces"), system_changes=all_changes
    ):
        super().calculate(atoms, properties, system_changes)
        z = atoms.positions[:, 2]
        self.results["energy"] += float(0.5 * self.tilt * np.sum(z**2))
        forces = self.results["forces"].copy()
        forces[:, 2] -= self.tilt * z
        self.results["forces"] = forces


def test_trivial_mode_floor_covers_every_rigid_body_mode_the_geometry_has():
    """A diatomic has five, not six: rotation about its own axis moves nothing."""
    config = SearchConfig(**CONFIG)
    _, _, floor, provenance = curvature_spectrum(_pair(1.2), double_well_factory, config)
    assert len(floor) == 5
    assert provenance["source"] == "fd"
    assert "analytic_unavailable" in provenance


def test_the_finite_difference_floor_is_the_step_own_chord_error_not_the_model():
    """Displacing an atom sideways by delta shortens no bond to first order, only to
    second, so the finite-difference transverse element is E'' delta^2 / (2 r^2) rather
    than zero -- the same chord-for-arc mistake the finite-amplitude probes avoid by
    construction. The double well is exactly invariant, so this floor is entirely the
    differentiation method's, and it must therefore scale as step^2 and vanish on the
    analytic path. Pinned because it is what decides whether a near-zero mode is readable:
    at a bond stiffness of 50 eV/A^2 on a hydrogen the same formula gives 2.5e-3, above
    the default eigenvalue tolerance.
    """
    stiffness = (1.6 / 0.36) * 2.0  # d2E/dr2 at the r=1.2 minimum
    distance, mass = 1.2, _pair(1.2).get_masses()[0]
    previous = None
    for step in (0.005, 0.01, 0.02):
        config = SearchConfig(minimum_check_step_A=step, **CONFIG)
        _, _, floor, _ = curvature_spectrum(
            _pair(distance), double_well_factory, config, source="fd"
        )
        worst = floor_residual(floor)
        predicted = stiffness * step**2 / (distance**2 * mass)
        assert worst == pytest.approx(predicted, rel=0.02), step
        if previous is not None:
            assert worst / previous == pytest.approx(4.0, rel=0.02)
        previous = worst

        _, _, analytic_floor, _ = curvature_spectrum(
            _pair(distance), WellWithHessian, config, source="analytic"
        )
        assert floor_residual(analytic_floor) < 1e-12


def test_trivial_mode_floor_exposes_a_non_invariant_potential_and_fails_closed():
    """A translation-dependent term is invisible in the forces and loud in the floor."""
    config = SearchConfig(**CONFIG)
    atoms = _pair(1.2)
    atoms.calc = TiltedWell(tilt=0.5)
    _, _, floor, _ = curvature_spectrum(atoms, lambda: TiltedWell(tilt=0.5), config)
    assert floor_residual(floor) > 1e-2

    confirmed, diagnostics = confirm_minimum(
        atoms, lambda: TiltedWell(tilt=0.5), config, seed=0
    )
    assert not confirmed
    assert diagnostics["reason"] == "trivial_mode_floor_above_tolerance"
    assert diagnostics["saddle_order"] is None
    # And the honest floor passes, so the gate is not simply always closed.
    clean = confirm_minimum(_pair(1.2), double_well_factory, config, seed=0)
    threshold = config.trivial_floor_fraction * config.minimum_check_eigenvalue_tol
    assert clean[0] and clean[1]["trivial_mode_floor_worst_residual"] < threshold


def test_analytic_hessian_agrees_with_finite_differences_and_is_recorded():
    config = SearchConfig(**CONFIG)
    factory = WellWithHessian
    for distance in (1.2, 1.8, 2.4):
        atoms = _pair(distance)
        atoms.calc = factory()
        fd = hessian_spectrum(atoms, factory, config, source="fd")
        analytic, _, _, provenance = curvature_spectrum(
            atoms, factory, config, source="analytic"
        )
        assert provenance["source"] == "analytic"
        assert provenance["finite_difference_deviation"] < config.analytic_hessian_check
        # The gap is the finite difference's own O(step^2) truncation, not disagreement.
        assert np.allclose(np.sort(analytic), np.sort(fd), atol=5e-4), distance


def test_a_wrong_unit_scale_on_the_analytic_hessian_is_caught_by_execution():
    """The trap: a calculator scales energy and forces but not its Hessian.

    With both factors at their default 1.0 nothing is wrong, so the failure only appears
    on a model that actually needs a conversion. Checking the method exists cannot see
    this; running it once and comparing can.
    """
    config = SearchConfig(**CONFIG)
    factory = lambda: WellWithHessian(corrupt=96.485)  # noqa: E731
    atoms = _pair(1.2)
    atoms.calc = factory()
    with pytest.raises(GateRejected) as caught:
        curvature_spectrum(atoms, factory, config, source="analytic")
    assert caught.value.code == "analytic_hessian_unavailable"

    # "auto" must not paper over it either: it falls back and says why.
    _, _, _, provenance = curvature_spectrum(atoms, factory, config, source="auto")
    assert provenance["source"] == "fd"
    assert "disagrees with a finite difference" in provenance["analytic_unavailable"]


def test_analytic_source_is_refused_when_the_backend_has_no_second_derivative():
    config = SearchConfig(**CONFIG)
    with pytest.raises(GateRejected) as caught:
        curvature_spectrum(_pair(1.2), double_well_factory, config, source="analytic")
    assert caught.value.code == "analytic_hessian_unavailable"
    with pytest.raises(ValueError):
        SearchConfig(hessian_source="analytic_maybe")


def test_band_regression_beats_the_finite_difference_hessian_it_replaces():
    """Extrapolating a band recovers the closed-form curvature; differencing once does not.

    The double well is quartic, so kappa(a) = kappa_0 + c a^2 holds exactly and the fit is
    exact too -- which is why the reported sigma collapses to rounding. That is the sigma
    behaving correctly, not a lucky case: it says the quadratic model explains every
    sample, and it would not say so on a band wide enough to need the next term.
    """
    config = SearchConfig(response_band_A=(0.005, 0.05), response_samples=6, **CONFIG)
    atoms = _pair(1.2)
    exact = ((1.6 / 0.36) * 2.0) / (atoms.get_masses()[0] / 2)
    guard = GuardedCalculator(double_well_factory(), config, None)
    probe = atoms.copy()
    probe.calc = guard
    along_bond = np.array([[-1.0, 0, 0], [1.0, 0, 0]])
    report = response_curvature(probe, guard, along_bond, config)
    finite_difference = hessian_spectrum(atoms, double_well_factory, config)[-1]

    assert report["kappa"] == pytest.approx(exact, rel=1e-9)
    assert abs(report["kappa"] - exact) < abs(finite_difference - exact)
    assert report["sigma"] < 1e-9
    assert report["exceeds_residual_scale"] and report["sign_stable"]
    assert len(report["samples"]) == config.response_samples
    assert "NOT an eigenvalue" in report["meaning"]


class FlankedTop(Calculator):
    """A maximum with minima on both sides: E = A cos(k u) + B cos(2 k u), u = r - 1.8.

    The measured shape of a torsional barrier top. Pushing far enough in either direction
    climbs back up past the flanking minima, so the symmetric second difference turns
    positive: the unstable direction is unstable at small amplitude and not at large.
    One number cannot describe that, which is what sign_stable is for.
    """

    implemented_properties = ["energy", "forces"]
    A, B, K = -1.0, 0.4, np.pi / 0.6

    def calculate(
        self, atoms=None, properties=("energy", "forces"), system_changes=all_changes
    ):
        super().calculate(atoms, properties, system_changes)
        vector = atoms.positions[1] - atoms.positions[0]
        r = float(np.linalg.norm(vector))
        u = r - 1.8
        energy = self.A * np.cos(self.K * u) + self.B * np.cos(2 * self.K * u)
        derivative = -self.A * self.K * np.sin(self.K * u) - 2 * self.B * self.K * np.sin(
            2 * self.K * u
        )
        force = derivative * vector / r
        self.results = {"energy": float(energy), "forces": np.array([force, -force])}


def test_a_band_that_reaches_past_the_flanking_minima_reports_itself_unstable_in_scale():
    narrow = SearchConfig(response_band_A=(0.02, 0.3), response_samples=5, **CONFIG)
    wide = SearchConfig(response_band_A=(0.05, 1.2), response_samples=8, **CONFIG)
    along_bond = np.array([[-1.0, 0, 0], [1.0, 0, 0]])
    reports = {}
    for config, expect_stable in ((narrow, True), (wide, False)):
        guard = GuardedCalculator(FlankedTop(), config, None)
        probe = Atoms("C2", positions=[[-0.9, 0, 0], [0.9, 0, 0]])
        probe.calc = guard
        report = response_curvature(probe, guard, along_bond, config)
        assert report["sign_stable"] is expect_stable, config.response_band_A
        if expect_stable:
            assert report["kappa"] < 0 and report["exceeds_residual_scale"]
        reports[expect_stable] = report
    # Cosines carry every order, so a band the quadratic model cannot cover is where the
    # standard error has to speak up rather than quietly reporting a biased number.
    assert reports[False]["sigma"] > 100 * reports[True]["sigma"]
    assert reports[False]["residual_rms"] > 100 * reports[True]["residual_rms"]


def test_the_energy_difference_filters_twice_as_gently_as_the_force_difference():
    """sinc^2(ka/2) against sinc(ka): the leading deviations differ by exactly two.

    This is the whole claim that a finite amplitude is a resolution rather than an error,
    reduced to something falsifiable. It also settles which difference to regress: only
    sinc^2 is non-negative, so only the energy form cannot invert a component's sign.
    """
    atoms = _pair(1.2)
    masses = atoms.get_masses()
    exact = ((1.6 / 0.36) * 2.0) / (masses[0] / 2)
    unit = np.array([[-1.0, 0, 0], [1.0, 0, 0]]) / np.sqrt(masses[:, None])
    unit /= np.linalg.norm(unit)  # mass-weighted unit direction
    reference = atoms.positions.copy()
    energy0 = atoms.get_potential_energy()

    for amplitude in (0.02, 0.05, 0.1):
        energies, gradients = {}, {}
        for sign in (1, -1):
            atoms.set_positions(reference + sign * amplitude * unit / np.sqrt(masses[:, None]))
            energies[sign] = atoms.get_potential_energy()
            gradients[sign] = -atoms.get_forces() / np.sqrt(masses[:, None])
        atoms.set_positions(reference)
        kappa_e = (energies[1] + energies[-1] - 2 * energy0) / amplitude**2
        kappa_f = float(np.sum((gradients[1] - gradients[-1]) * unit)) / (2 * amplitude)
        ratio = (kappa_f - exact) / (kappa_e - exact)
        assert ratio == pytest.approx(2.0, rel=0.05), amplitude


def test_the_response_gradient_is_the_derivative_of_the_response_curvature():
    """One pair of evaluations gives both the objective and its gradient.

    d/dv of the energy form is exactly twice the odd part of the gradient, so rotating a
    direction downhill in response curvature costs two force evaluations per step: the
    price of one Hessian-vector product, with no matrix and no spectrum.
    """
    config = SearchConfig(**CONFIG)
    atoms = _pair(1.5)
    masses = atoms.get_masses()
    guard = GuardedCalculator(double_well_factory(), config, None)
    probe = atoms.copy()
    probe.calc = guard
    amplitude = 0.05

    direction = np.array([[-1.0, 0.4, 0.0], [1.0, 0.2, -0.3]])
    weighted = direction / np.sqrt(masses[:, None])
    weighted /= np.linalg.norm(weighted)

    base_energy = float(guard.get_potential_energy(probe))

    def kappa_at(mass_weighted_direction):
        cartesian = mass_weighted_direction / np.sqrt(masses[:, None])
        energies = []
        for sign in (1, -1):
            probe.set_positions(atoms.positions + sign * amplitude * cartesian)
            energies.append(float(guard.get_potential_energy(probe)))
        probe.set_positions(atoms.positions)
        return (energies[0] + energies[1] - 2 * base_energy) / amplitude**2

    step = 1e-5
    numerical = np.zeros_like(weighted)
    for index in np.ndindex(weighted.shape):
        shifted = weighted.copy()
        shifted[index] += step
        plus = kappa_at(shifted)
        shifted = weighted.copy()
        shifted[index] -= step
        numerical[index] = (plus - kappa_at(shifted)) / (2 * step)

    analytic = response_gradient(probe, guard, direction, amplitude)
    assert np.allclose(analytic, numerical, atol=1e-3 * max(1.0, np.abs(numerical).max()))


def test_mass_weighting_maps_a_cartesian_direction_the_right_way():
    """Cartesian to mass-weighted is multiply by sqrt(m), and it must be, not divide.

    On a homonuclear pair the two differ only by a scale factor that normalization removes,
    so nothing notices. With unequal masses they pick different directions, and the
    Rayleigh quotient along the wrong one is simply a different number. Pinned with the
    heaviest-to-lightest ratio a molecule offers, carbon against hydrogen.
    """
    from prrs.runner import direction_stiffness

    config = SearchConfig(response_band_A=(0.002, 0.02), response_samples=5, **CONFIG)
    atoms = Atoms("CH", positions=[[-0.6, 0, 0], [0.6, 0, 0]])
    atoms.calc = double_well_factory()
    masses = atoms.get_masses()
    assert masses[0] / masses[1] > 10

    eigenvalues, vectors = hessian_spectrum(
        atoms, double_well_factory, config, return_vectors=True
    )
    mode = vectors[:, -1].reshape(-1, 3)
    cartesian = mode / np.sqrt(masses[:, None])  # the displacement it corresponds to

    assert direction_stiffness(atoms, double_well_factory, config, cartesian) == pytest.approx(
        eigenvalues[-1], rel=1e-3
    )
    # The wrong convention lands on a different direction, so it reports a different
    # curvature; if this ever stops differing, the test has stopped testing anything.
    wrong = cartesian / masses[:, None]
    assert direction_stiffness(atoms, double_well_factory, config, wrong) != pytest.approx(
        eigenvalues[-1], rel=1e-2
    )

    guard = GuardedCalculator(double_well_factory(), config, None)
    probe = atoms.copy()
    probe.calc = guard
    report = response_curvature(probe, guard, cartesian, config)
    assert report["kappa"] == pytest.approx(eigenvalues[-1], rel=1e-3)


def test_a_rotational_floor_is_mostly_the_residual_gradient_and_is_subtracted():
    """A structure that is not quite stationary shows a large rigid-body curvature.

    E(R(s)x) = E(x) differentiated twice gives xdot^T H xdot = omega^2 (g . x_perp), so
    the rotational quotient is proportional to the gradient and says nothing about the
    quality of the calculation until that term is removed. Gating on the raw number closes
    the gate on every geometry converged to an ordinary force tolerance, which is what it
    did before this was separated out. Translations have no such term at any geometry.
    """
    config = SearchConfig(**CONFIG)
    for distance in (1.5, 2.0):  # away from both minima
        atoms = Atoms("CH", positions=[[-distance / 2, 0, 0], [distance / 2, 0, 0]])
        atoms.calc = double_well_factory()
        assert np.abs(atoms.get_forces()).max() > 0.1
        _, _, floor, _ = curvature_spectrum(atoms, double_well_factory, config, source="fd")
        rotations = [entry for entry in floor if entry["kind"] == "rotation"]
        translations = [entry for entry in floor if entry["kind"] == "translation"]
        assert rotations and translations
        for entry in translations:
            assert entry["gradient_expectation"] == 0.0
            assert abs(entry["measured"]) < 1e-12
        for entry in rotations:
            assert abs(entry["measured"]) > 0.1  # gradient dominated
            assert entry["gradient_expectation"] == pytest.approx(entry["measured"], rel=1e-3)
            assert abs(entry["residual"]) < 1e-3 * abs(entry["measured"])
        # So the geometry passes the gate despite a raw floor a thousand times the
        # tolerance, which is the whole point of subtracting the expectation.
        assert (
            floor_residual(floor)
            < config.trivial_floor_fraction * config.minimum_check_eigenvalue_tol
        )


def test_the_probe_measures_curvature_instead_of_following_a_displacement():
    """A negative number is a stronger signal than "it slid away", and costs less.

    On a diatomic every direction with the rigid-body part removed is the bond direction,
    so the probe's lowest curvature must equal the one vibrational eigenvalue exactly.
    """
    config = SearchConfig(minimum_check="probe", **CONFIG)
    reduced = _pair(1.2).get_masses()[0] / 2
    for distance, expect_minimum in ((1.2, True), (2.4, True), (1.8, False)):
        atoms = _pair(distance)
        confirmed, diagnostics = confirm_minimum(atoms, double_well_factory, config, seed=0)
        assert confirmed is expect_minimum, distance
        assert diagnostics["method"] == "probe"
        assert len(diagnostics["probes"]) == config.minimum_check_probes
        x = (distance - 1.8) / 0.6
        expected = (1.6 / 0.36) * (3 * x * x - 1) / reduced  # closed form
        assert diagnostics["lowest_curvature"] == pytest.approx(expected, rel=1e-3)
        assert (diagnostics["negative_directions"] > 0) is (not expect_minimum)


def test_the_probe_finds_a_torsional_saddle_and_names_the_coordinate():
    """The failure the old probe had: a torsional barrier top it can slide off in two
    equivalent ways, where a random displacement may return to the top.

    The chain's torsion term is (V/2)(1 - cos 3 phi), so phi = 60 degrees is a maximum of
    known curvature -4.5 V, with every bond and angle already at its own minimum. Probing
    the molecule's own rotatable torsion finds it directly, and says which coordinate it
    was -- the reaction coordinate, with a chemically readable name and no eigenvector.
    """
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).parent))
    from test_convergence import Chain, _chain

    config = SearchConfig(minimum_check="probe", **CONFIG)
    top = _chain(60.0)
    confirmed, diagnostics = confirm_minimum(top, Chain, config, seed=0)
    assert not confirmed
    torsional = [entry for entry in diagnostics["probes"] if entry["source"] == "torsion"]
    assert len(torsional) == 1
    assert torsional[0]["negative"] and torsional[0]["curvature_eV_A2_amu"] < 0
    assert torsional[0]["indices"] == [0, 1, 2, 3]

    # Same coordinate at its minimum: positive, and the structure is confirmed.
    basin = _chain(0.0)
    confirmed, diagnostics = confirm_minimum(basin, Chain, config, seed=0)
    assert confirmed
    torsional = [entry for entry in diagnostics["probes"] if entry["source"] == "torsion"]
    assert torsional[0]["curvature_eV_A2_amu"] > 0
    assert diagnostics["negative_directions"] == 0


def test_the_probe_costs_two_force_evaluations_per_direction():
    calls = {"n": 0}

    class Counted(DoubleWell):
        def calculate(
            self, atoms=None, properties=("energy", "forces"), system_changes=all_changes
        ):
            calls["n"] += 1
            super().calculate(atoms, properties, system_changes)

    config = SearchConfig(minimum_check="probe", minimum_check_probes=3, **CONFIG)
    atoms = _pair(1.2)
    atoms.calc = Counted()
    atoms.get_potential_energy()
    calls["n"] = 0
    confirmed, diagnostics = confirm_minimum(atoms, Counted, config, seed=0)
    directions = len(diagnostics["probes"])
    assert directions == 3
    # One reference evaluation for the energy-rise gate, then two per direction.
    assert calls["n"] <= 2 * directions + 2


# ---------------------------------------------------------------------------
# A coordinate that carries no energy must not multiply states.
# ---------------------------------------------------------------------------


def _free_chain(phi_deg, barrier=0.0):
    """The four-atom chain with its torsion term switched off: phi is exactly flat."""
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).parent))
    from test_convergence import Chain, _chain

    atoms = _chain(phi_deg, barrier)
    atoms.calc = Chain(barrier)
    return atoms


def test_a_flat_torsion_is_quotiented_out_of_the_comparison():
    """Two geometries differing only along a barrierless coordinate are one microstate.

    Plain RMSD puts them far apart, so a genuinely free rotor would mint an unbounded
    number of spurious conformers. Rotating the comparison over that coordinate collapses
    them, and costs no force evaluations because the rotation is exact geometry.
    """
    from prrs.chemistry import automorphisms, free_aligned_symmetric_rmsd, symmetric_rmsd

    a, b = _free_chain(0.0), _free_chain(137.0)
    permutations_, _ = automorphisms(a)

    plain = symmetric_rmsd(a, b, permutations_)
    assert plain > 0.5, plain
    assert free_aligned_symmetric_rmsd(a, b, permutations_, free_bonds=()) == pytest.approx(
        plain
    )

    aligned = free_aligned_symmetric_rmsd(a, b, permutations_, free_bonds=((1, 2),), samples=72)
    assert aligned < 0.05, aligned
    assert aligned < 0.1 * plain


def test_only_a_measured_free_coordinate_is_quotiented_never_a_soft_one():
    """The verdict is energetic, so it is read off the quench rather than the graph.

    Merging a merely soft coordinate would destroy real conformers -- ethanol's anti and
    gauche differ along a torsion of 0.21 eV/rad^2 -- so nothing may infer freeness from
    connectivity. A structure with no quench report contributes no free coordinates at all.
    """
    from prrs.network import free_bonds_of
    from prrs.reliability import GuardedCalculator
    from prrs.runner import quench

    assert free_bonds_of(_free_chain(0.0)) == ()  # no report travels with it

    for barrier, expected in ((0.0, ((1, 2),)), (0.05, ())):
        atoms = _free_chain(30.0, barrier)
        config = SearchConfig(quench_fmax_eV_A=0.002, soft_polish_rounds=4)
        guard = GuardedCalculator(atoms.calc, config, None)
        atoms.calc = guard
        converged, diagnostics = quench(atoms, guard, config)
        assert converged, barrier
        tiers = [m["tier"] for m in diagnostics["rounds"][-1]["modes"]]
        atoms.info["quench"] = diagnostics
        assert free_bonds_of(atoms) == expected, (barrier, tiers)


def test_the_registry_admits_one_microstate_for_a_free_rotor_not_many():
    """The gap this closes: a free coordinate generating unbounded fake conformers."""
    from ase.calculators.singlepoint import SinglePointCalculator
    from prrs.network import Registry
    from prrs.reliability import GuardedCalculator
    from prrs.runner import quench

    def quenched(phi, barrier):
        atoms = _free_chain(phi, barrier)
        config = SearchConfig(quench_fmax_eV_A=0.002, soft_polish_rounds=4)
        guard = GuardedCalculator(atoms.calc, config, None)
        atoms.calc = guard
        converged, diagnostics = quench(atoms, guard, config)
        assert converged
        result = atoms.copy()
        result.calc = SinglePointCalculator(
            result, energy=atoms.get_potential_energy(), forces=atoms.get_forces()
        )
        result.info["quench"] = diagnostics
        return result

    # A three-fold torsion has minima at 0, 120 and 240 degrees, and the four starting
    # geometries relax into those three; a flat torsion has no minima at all, so all four
    # are one state.
    for barrier, expected_states in ((0.0, 1), (0.05, 3)):
        # A bare C4 chain has no neutral closed-shell valence assignment, so the
        # domain gate is meaningless here; this fixture is about the free rotor.
        registry = Registry(
            SearchConfig(
                conformer_max_per_node=8, chemical_max_nodes=4, closed_shell_only=False
            )
        )
        node = None
        for phi in (0.0, 40.0, 137.0, 231.0):
            admission = registry.admit(quenched(phi, barrier), node, "t0")
            if node is None:
                node = admission.node
        assert len(node.microstates) == expected_states, (barrier, len(node.microstates))
        if barrier == 0.0:
            assert node.microstates[0].free_bonds == ((1, 2),)


# ---------------------------------------------------------------------------
# The bottom direction without a spectrum: sphere descent on the response.
# ---------------------------------------------------------------------------


def test_sphere_descent_finds_the_softest_direction_without_a_hessian():
    """Descent on kappa_a^E must reach the lowest eigenvalue of the projected Hessian.

    At vanishing amplitude the response curvature is the Rayleigh quotient, so its minimum
    over the sphere is the lowest eigenvalue -- which is the check: the descent never builds
    or diagonalizes anything, and has to land on the same number.
    """
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).parent))
    from test_convergence import Chain, _chain
    from prrs.runner import curvature_spectrum, lowest_response_direction
    from prrs.reliability import GuardedCalculator

    # Default tolerance: the tangential fraction of the gradient below 1e-3.
    config = SearchConfig(
        response_band_A=(0.002, 0.02), response_descent_max_steps=200, quench_fmax_eV_A=0.002
    )
    for phi, expect_negative in ((0.0, False), (60.0, True)):
        atoms = _chain(phi)
        eigenvalues, *_ = curvature_spectrum(atoms, Chain, config, source="fd")
        internal = [v for v in eigenvalues if abs(v) > 1e-9]
        lowest = min(internal)
        assert bool(lowest < 0) is expect_negative, phi

        guard = GuardedCalculator(Chain(), config, None)
        probe = atoms.copy()
        probe.calc = guard
        rng = np.random.default_rng(3)
        found = lowest_response_direction(
            probe, guard, config, rng.normal(size=(len(atoms), 3)), amplitude=0.002
        )
        assert found["converged"], (phi, found["final_tangent_fraction"])
        assert found["curvature"] == pytest.approx(lowest, rel=2e-2), (phi, found["curvature"])
        assert bool(found["curvature"] < 0) is expect_negative
        assert found["rotation_from_initial_deg"] > 5.0  # it actually moved
        assert "NOT an eigenvector" in found["meaning"]


def test_sphere_descent_costs_two_force_evaluations_per_step():
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).parent))
    from test_convergence import Chain, _chain
    from prrs.runner import lowest_response_direction
    from prrs.reliability import GuardedCalculator

    calls = {"n": 0}

    class Counted(Chain):
        def calculate(
            self, atoms=None, properties=("energy", "forces"), system_changes=all_changes
        ):
            calls["n"] += 1
            super().calculate(atoms, properties, system_changes)

    config = SearchConfig(
        response_descent_max_steps=8, response_descent_tol=1e-12, quench_fmax_eV_A=0.002
    )
    atoms = _chain(60.0)
    guard = GuardedCalculator(Counted(), config, None)
    probe = atoms.copy()
    probe.calc = guard
    rng = np.random.default_rng(5)
    calls["n"] = 0
    found = lowest_response_direction(
        probe, guard, config, rng.normal(size=(len(atoms), 3)), amplitude=0.005
    )
    # Each pair of +-a v evaluations yields both the objective and its gradient, so the
    # bookkeeping the function reports has to be the real number of physical evaluations.
    assert calls["n"] == found["force_evaluations"] + 1  # plus the reference energy


def test_persistence_reports_the_band_and_refuses_to_rule_on_it():
    """The direction is tracked across amplitudes; no threshold is invented for it.

    On a chain whose torsion is the only soft coordinate the bottom direction barely moves,
    so this pins the reporting, not a verdict. What the tolerance should be cannot be
    settled on a molecule where nothing exercises it.
    """
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).parent))
    from test_convergence import Chain, _chain
    from prrs.runner import response_persistence
    from prrs.reliability import GuardedCalculator

    config = SearchConfig(
        response_band_A=(0.002, 0.05),
        response_samples=4,
        response_descent_max_steps=200,
        quench_fmax_eV_A=0.002,
    )
    atoms = _chain(60.0)
    guard = GuardedCalculator(Chain(), config, None)
    probe = atoms.copy()
    probe.calc = guard
    rng = np.random.default_rng(11)
    report = response_persistence(probe, guard, config, rng.normal(size=(len(atoms), 3)))

    assert len(report["band"]) == config.response_samples
    assert all(entry["converged"] for entry in report["band"])
    assert all(entry["curvature"] < 0 for entry in report["band"])
    assert report["sign_stable"] is True
    assert report["band"][0]["rotation_from_previous_deg"] is None
    assert all(entry["rotation_from_previous_deg"] is not None for entry in report["band"][1:])
    assert report["total_rotation_deg"] < 5.0
    assert "no persistence verdict" in report["meaning"]


def test_min_mode_following_climbs_to_the_saddle_a_quench_rolls_off():
    """Inverting the force along the lowest mode reaches a ridge that minimizing cannot.

    The four-atom chain's torsional maximum is a first-order saddle of known curvature.
    Started from a geometry displaced off it -- the shape a probe near the dividing surface
    produces -- a plain quench falls into a basin and min-mode following climbs back.
    """
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).parent))
    from test_convergence import Chain, _chain
    from prrs.runner import curvature_spectrum, follow_min_mode, quench
    from prrs.reliability import GuardedCalculator

    config = SearchConfig(quench_fmax_eV_A=0.002, min_mode_fmax_eV_A=0.002)
    displaced = _chain(48.0)  # off the phi = 60 ridge

    # A quench from here leaves the ridge entirely.
    rolled = displaced.copy()
    guard = GuardedCalculator(Chain(), config, None)
    rolled.calc = guard
    assert quench(rolled, guard, config)[0]
    eigenvalues, *_ = curvature_spectrum(rolled, Chain, config, source="fd")
    assert min(eigenvalues) > -config.minimum_check_eigenvalue_tol, "quench kept a saddle"

    saddle, diagnostics = follow_min_mode(displaced, Chain, config)
    assert saddle is not None, diagnostics
    assert diagnostics["saddle_order"] == 1
    assert diagnostics["fmax_eV_A"] < config.min_mode_fmax_eV_A
    assert diagnostics["negative_eigenvalues"][0] < 0
    assert diagnostics["imaginary_wavenumbers_icm"][0] < 0
    assert "unstable_mode" in diagnostics

    # It is the torsional maximum: (V/2)(1 - cos 3 phi) has curvature -4.5 V there.
    from prrs import internal

    phi = np.degrees(internal.coordinate(saddle.positions, "dihedral", (0, 1, 2, 3))) % 120
    assert min(abs(phi - 60), abs(phi - 60 + 120)) < 3.0, phi


def test_a_step_that_overshoots_the_ceiling_is_backtracked_not_fatal():
    """One overshoot is a step-length problem and must not be reported as a domain problem.

    Measured on malonaldehyde: the transient seed that had reached the proton-transfer
    saddle was refused because a fixed ascent step carried it over the line on the way,
    turning a recoverable overshoot into a failed search. The chain's torsional maximum is
    0.05 eV up, so a ceiling just above it plus a deliberately long step forces exactly that
    situation, and the walk still has to arrive.
    """
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).parent))
    from test_convergence import Chain, _chain
    from prrs.runner import follow_min_mode

    basin = _chain(0.0)
    basin.calc = Chain()
    floor = basin.get_potential_energy()
    config = SearchConfig(
        quench_fmax_eV_A=0.002,
        min_mode_fmax_eV_A=0.002,
        min_mode_step_A=0.30,  # long enough to overshoot
        min_mode_max_steps=300,
        quench_steps=400,
    )
    saddle, diagnostics = follow_min_mode(
        _chain(20.0), Chain, config, ceiling=floor + 0.06, reference=floor
    )
    assert saddle is not None, diagnostics
    assert diagnostics["saddle_order"] == 1
    assert diagnostics["backtracks"] > 0, "the test did not exercise the overshoot"
    assert diagnostics["ceiling_hit_step"] is not None
    assert diagnostics["max_rise_eV"] <= 0.06 + 1e-9
    saddle.calc = Chain()
    assert saddle.get_potential_energy() - floor == pytest.approx(0.05, abs=0.005)


def test_a_ceiling_cannot_be_walked_around_by_taking_smaller_steps():
    """The bound is on the region, so shortening the step must not get past it."""
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).parent))
    from test_convergence import Chain, _chain
    from prrs.runner import follow_min_mode

    basin = _chain(0.0)
    basin.calc = Chain()
    floor = basin.get_potential_energy()
    config = SearchConfig(
        quench_fmax_eV_A=0.002,
        min_mode_fmax_eV_A=0.002,
        min_mode_max_steps=300,
        quench_steps=400,
    )
    # Start below the ceiling so the verdict is about the ridge, not about the seed.
    saddle, diagnostics = follow_min_mode(
        _chain(20.0), Chain, config, ceiling=floor + 0.02, reference=floor
    )
    assert saddle is None
    assert diagnostics["reason"] in (
        "ceiling_blocked",
        "above_ceiling_at_saddle",
        "not_converged",
    )
    assert diagnostics["max_rise_eV"] <= 0.02 + 1e-9


def test_a_seed_already_above_the_ceiling_is_skipped_not_failed():
    """Its verdict is about the seed, so the next seed must not inherit it."""
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).parent))
    from test_convergence import Chain, _chain
    from prrs.runner import follow_min_mode

    basin = _chain(0.0)
    basin.calc = Chain()
    floor = basin.get_potential_energy()
    config = SearchConfig(quench_fmax_eV_A=0.002, min_mode_fmax_eV_A=0.002)
    saddle, diagnostics = follow_min_mode(
        _chain(60.0), Chain, config, ceiling=floor + 0.01, reference=floor
    )
    assert saddle is None
    assert diagnostics["reason"] == "seed_above_ceiling"
    assert diagnostics["steps"] == 0
    assert diagnostics["initial_rise_eV"] > 0.01


# ---------------------------------------------------------------------------
# Following a target coordinate is a different question from finding the softest saddle.
# ---------------------------------------------------------------------------


def _two_mode_chain():
    """The four-atom chain with a soft torsion and a much softer methyl-like rotor.

    Two ridges of very different curvature, so "the lowest mode" and "the mode I asked
    about" are different answers and choosing wrongly is visible.
    """
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).parent))
    from test_convergence import Chain, _chain

    return Chain, _chain


def test_a_target_direction_is_actually_used_to_choose_the_mode():
    """The bug this pins: `direction` was accepted and ignored, so the walk always climbed
    the globally softest mode. On 3-oxobutanal that returned eight rotor barriers at -65 to
    -80 cm^-1 while the proton transfer sits at -3186.

    Checked on the first step rather than the last: the mode indices are sorted by
    eigenvalue and reorder as the structure moves -- at a saddle the negative mode sorts to
    zero -- so the seed geometry is where "which mode was chosen" is a meaningful question.
    """
    Chain, build = _two_mode_chain()
    from prrs.runner import curvature_spectrum, follow_min_mode

    config = SearchConfig(
        quench_fmax_eV_A=0.002,
        min_mode_fmax_eV_A=0.002,
        min_mode_max_steps=60,
        quench_steps=400,
    )
    atoms = build(20.0)
    eigenvalues, vectors, _, _ = curvature_spectrum(atoms, Chain, config, source="fd")
    tolerance = config.minimum_check_eigenvalue_tol
    internal = [i for i, v in enumerate(eigenvalues) if abs(v) > tolerance]
    assert len(internal) >= 3, internal
    masses = atoms.get_masses()

    def first_step(diagnostics):
        history = diagnostics.get("history") or []
        return history[0] if history else None

    # With no target: the softest non-trivial mode, which is the behaviour that was there
    # before and must not change.
    _, unaimed = follow_min_mode(atoms, Chain, config)
    assert unaimed["target_supplied"] is False
    start = first_step(unaimed)
    if start is not None:
        assert start["followed_order"] == internal[0]

    # With a target: a different mode, and the overlap that justified choosing it.
    aimed_at = internal[2]
    target = vectors[:, aimed_at].reshape(-1, 3) / np.sqrt(masses[:, None])
    _, aimed = follow_min_mode(atoms, Chain, config, direction=target)
    assert aimed["target_supplied"] is True
    start = first_step(aimed)
    if start is not None:
        assert start["followed_order"] == aimed_at, (start, aimed_at)
        assert start["mode_overlap"] > 0.9
    else:
        # No step was taken, so it has to say why rather than report a result.
        assert aimed.get("reason"), aimed


def test_losing_the_target_mode_is_reported_rather_than_silently_switched():
    """If no eigenvector still resembles the coordinate, continuing would climb a different
    reaction. That has to be a stated outcome, not a quiet substitution."""
    Chain, build = _two_mode_chain()
    from prrs.runner import follow_min_mode

    config = SearchConfig(
        quench_fmax_eV_A=0.002,
        min_mode_fmax_eV_A=0.002,
        min_mode_overlap=0.999999,  # nothing can keep up
        min_mode_max_steps=50,
    )
    atoms = build(20.0)
    rng = np.random.default_rng(0)
    _, diagnostics = follow_min_mode(
        atoms, Chain, config, direction=rng.normal(size=(len(atoms), 3))
    )
    assert diagnostics["reason"] == "lost_target_mode"
    assert diagnostics["mode_overlap"] < config.min_mode_overlap


def test_a_line_search_that_accepts_nothing_is_stalled_not_converged():
    """Reporting a stall as convergence turns a failure into a result."""
    from prrs.runner import lowest_response_direction
    from prrs.reliability import GuardedCalculator

    Chain, build = _two_mode_chain()
    config = SearchConfig(
        quench_fmax_eV_A=0.002,
        response_descent_step=1e-12,
        response_descent_max_steps=3,
        response_descent_tol=1e-12,
    )
    atoms = build(60.0)
    guard = GuardedCalculator(Chain(), config, None)
    probe = atoms.copy()
    probe.calc = guard
    rng = np.random.default_rng(1)
    found = lowest_response_direction(
        probe, guard, config, rng.normal(size=(len(atoms), 3)), amplitude=0.01
    )
    assert found["stalled"] is True
    assert found["converged"] is False


def test_the_response_error_bar_is_not_called_a_significance():
    """sigma is the residual scale of a deterministic fit. It measures how well the a^2 model
    describes the samples and is blind to a systematic error such as probing the wrong
    direction -- measured on malonaldehyde it read 3.7e-06 while the value was wrong by 0.18.
    """
    from prrs.runner import response_curvature
    from prrs.reliability import GuardedCalculator

    config = SearchConfig(response_band_A=(0.005, 0.05), response_samples=6, **CONFIG)
    atoms = _pair(1.2)
    guard = GuardedCalculator(double_well_factory(), config, None)
    probe = atoms.copy()
    probe.calc = guard
    report = response_curvature(probe, guard, np.array([[-1.0, 0, 0], [1.0, 0, 0]]), config)
    assert "significant" not in report
    assert "exceeds_residual_scale" in report
    assert "not a calibrated" in report["sigma_meaning"]


def test_the_probe_mode_is_chosen_even_when_a_rotor_is_softer():
    """The exact P2 failure, in a system where the right answer is known.

    The chain with a nearly flat torsion makes that coordinate the softest thing in the
    molecule by two orders of magnitude, so climbing "the lowest mode" goes there -- which is
    what returned eight rotor barriers on 3-oxobutanal while the proton transfer sat at
    -3186. Pointing the walk at a stiff coordinate has to override that, and the claim is
    about which mode is chosen, so it is checked on the first step: whether the walk then
    converges is a separate question.
    """
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).parent))
    from test_convergence import Chain, _chain
    from prrs.runner import curvature_spectrum, follow_min_mode

    def factory():
        return Chain(barrier_eV=0.004)  # a nearly free torsion

    config = SearchConfig(
        quench_fmax_eV_A=0.002,
        min_mode_fmax_eV_A=0.005,
        min_mode_max_steps=40,
        min_mode_step_A=0.02,
    )
    atoms = _chain(20.0, 0.004)
    atoms.calc = factory()
    eigenvalues, vectors, floor, _ = curvature_spectrum(atoms, factory, config, source="fd")
    # Invariant to within the finite difference's own chord error at this step, which is
    # what the gate compares against; a non-invariant potential lands orders of magnitude
    # above this (the first version of this test used absolute-position anchors and read
    # 0.389, and the gate refused it).
    assert (
        floor_residual(floor)
        < config.trivial_floor_fraction * config.minimum_check_eigenvalue_tol
    )
    tolerance = config.minimum_check_eigenvalue_tol
    internal = [i for i, v in enumerate(eigenvalues) if abs(v) > tolerance]
    softest, stiff = internal[0], internal[-1]
    assert abs(eigenvalues[stiff]) > 100 * abs(eigenvalues[softest]), (
        float(eigenvalues[softest]),
        float(eigenvalues[stiff]),
    )

    masses = atoms.get_masses()
    tangent = vectors[:, stiff].reshape(-1, 3) / np.sqrt(masses[:, None])
    _, aimed = follow_min_mode(atoms, factory, config, direction=tangent)
    assert aimed["target_supplied"] is True
    first = (aimed.get("history") or [None])[0]
    assert first is not None, aimed
    assert first["followed_order"] == stiff, (first, softest, stiff)
    assert first["mode_overlap"] > 0.9
    assert aimed["target_overlap"] > 0.9

    # Without the tangent it goes to the softest mode, which is the wrong question here.
    _, unaimed = follow_min_mode(atoms, factory, config)
    assert unaimed["target_supplied"] is False
    assert (unaimed.get("history") or [{}])[0].get("followed_order") == softest


def test_a_tangent_with_no_internal_component_is_refused():
    """A pure translation is not a coordinate. Substituting the softest mode for it is how a
    proton transfer becomes a methyl rotor, so it has to be a stated refusal."""
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).parent))
    from test_convergence import Chain, _chain
    from prrs.runner import follow_min_mode

    config = SearchConfig(quench_fmax_eV_A=0.002, min_mode_fmax_eV_A=0.002)
    atoms = _chain(20.0)
    translation = np.tile(np.array([1.0, 0.0, 0.0]), (len(atoms), 1))
    saddle, diagnostics = follow_min_mode(atoms, Chain, config, direction=translation)
    assert saddle is None
    assert diagnostics["reason"] == "missing_target_direction"


def test_tracking_follows_one_branch_through_an_eigenvalue_crossing():
    """Mode indices are sorted by eigenvalue and reorder as the structure moves -- at a
    saddle the negative mode sorts to zero. Continuity has to come from the overlap, so the
    followed order is allowed to change while the overlap stays high."""
    Chain, build = _two_mode_chain()
    from prrs.runner import curvature_spectrum, follow_min_mode

    config = SearchConfig(
        quench_fmax_eV_A=0.002,
        min_mode_fmax_eV_A=0.002,
        min_mode_refresh=1,
        min_mode_max_steps=200,
        quench_steps=400,
    )
    atoms = build(20.0)
    eigenvalues, vectors, _, _ = curvature_spectrum(atoms, Chain, config, source="fd")
    tolerance = config.minimum_check_eigenvalue_tol
    internal = [i for i, v in enumerate(eigenvalues) if abs(v) > tolerance]
    masses = atoms.get_masses()
    tangent = vectors[:, internal[0]].reshape(-1, 3) / np.sqrt(masses[:, None])

    saddle, diagnostics = follow_min_mode(atoms, Chain, config, direction=tangent)
    history = diagnostics.get("history") or []
    assert history, diagnostics
    orders = {entry["followed_order"] for entry in history}
    assert len(orders) > 1, "the walk never crossed, so this test proves nothing"
    for entry in history:
        assert entry["mode_overlap"] >= config.min_mode_overlap, entry
    if saddle is not None:
        assert diagnostics["saddle_order"] == 1


def test_no_walk_ever_puts_a_nonfinite_number_in_the_record():
    """Ten of follow_min_mode's eleven returns are early ones, and they all carry `bounds`.

    max_rise_eV used to start at -inf and be mapped back to None only on the successful
    return, so every early return wrote a -inf into the record. `atomic_json` writes with
    allow_nan=False, so the write raised -- and the write that raised was publish(), which
    is also what the failure handler calls, so the run could not record that it had failed:
    status stayed "running", indistinguishable from being killed by a signal. Three seeds of
    the P2 T=300 scan died that way.

    This pins the property rather than the one path that happened to fire: whatever
    follow_min_mode returns must be writable.
    """
    from prrs.io import dumps_strict, nonfinite_paths
    from prrs.runner import follow_min_mode

    class Flat(Calculator):
        """No curvature anywhere, so mode selection finds no internal mode at step 0 --
        an early return taken before the energy block has touched `bounds`."""

        implemented_properties = ["energy", "forces"]

        def calculate(
            self, atoms=None, properties=("energy", "forces"), system_changes=all_changes
        ):
            super().calculate(atoms, properties, system_changes)
            self.results = {"energy": 0.0, "forces": np.zeros((len(atoms), 3))}

    atoms = Atoms("C2", positions=[[0.0, 0, 0], [1.5, 0, 0]])
    config = SearchConfig(closed_shell_only=False)
    structure, diagnostics = follow_min_mode(
        atoms, lambda: Flat(), config, reference=0.0, ceiling=1.0
    )
    assert structure is None
    assert diagnostics["reason"] == "no_internal_mode"
    assert nonfinite_paths(diagnostics) == [], nonfinite_paths(diagnostics)
    dumps_strict(diagnostics)  # must not raise
    assert diagnostics["max_rise_eV"] is None  # None, never -inf


def test_the_json_guard_names_the_field_it_rejected():
    """The message that cost a wrong diagnosis. "Out of range float values are not JSON
    compliant: -inf" named nothing, so the failure looked like an external kill."""
    from prrs.io import dumps_strict, nonfinite_paths

    record = {
        "ok": 1.0,
        "walk": {"max_rise_eV": float("-inf")},
        "history": [{"fmax": 0.1}, {"fmax": float("nan")}],
    }
    assert sorted(where for where, _ in nonfinite_paths(record)) == [
        "history[1].fmax",
        "walk.max_rise_eV",
    ]
    with pytest.raises(ValueError) as caught:
        dumps_strict(record)
    message = str(caught.value)
    assert "walk.max_rise_eV" in message and "history[1].fmax" in message
    assert "fix the source rather than relaxing this check" in message


def test_translation_is_projected_out_in_the_mass_weighted_metric():
    """A pure translation must leave nothing behind, on a heteronuclear molecule too.

    In mass-weighted coordinates the translation along axis a is (T_a)_ib = sqrt(m_i)
    delta_ab, not the uniform vector, so subtracting an unweighted column mean removes
    the wrong direction. On C2 the two coincide and nothing shows; on anything with more
    than one element a residual translation survives and is then read as travel. This
    feeds align_mode_with_travel, whose whole job is a sign, so the residual does not
    have to be large to flip a relay onto the wrong side of a saddle.
    """
    from prrs.runner import _projected_direction_from

    water = molecule("H2O")
    masses = water.get_masses()
    assert masses.max() / masses.min() > 10  # the case C2 cannot exercise

    for axis in np.eye(3):
        translation = np.tile(axis, (len(water), 1))
        assert _projected_direction_from(translation, masses) is None

    # An internal displacement must survive untouched when a translation rides along.
    internal_only = np.zeros((len(water), 3))
    internal_only[1] = [0.05, -0.02, 0.01]
    internal_only[2] = [-0.05, 0.02, -0.01]
    alone = _projected_direction_from(internal_only, masses)
    with_translation = _projected_direction_from(
        internal_only + np.tile([0.3, -0.7, 0.2], (len(water), 1)), masses
    )
    assert np.abs(alone - with_translation).max() < 1e-12

    # And what comes back really is orthogonal to every mass-weighted translation.
    root = np.sqrt(masses[:, None])
    for axis in np.eye(3):
        mode = root * np.tile(axis, (len(water), 1))
        mode = mode / np.linalg.norm(mode)
        assert abs(float(np.sum(alone * mode))) < 1e-12


def test_hessian_vector_product_restores_the_geometry_when_a_force_call_fails():
    """The guard refusing a displaced evaluation is its job, not an exceptional path.

    Without the restore the structure is left standing at R +/- h v, and every later
    measurement is silently taken on a geometry nobody asked for. _response_and_gradient
    already had the finally; this is the sibling that did not.
    """
    from prrs.runner import hessian_vector_product

    class FailsOnSecond(Calculator):
        implemented_properties = ["energy", "forces"]

        def __init__(self):
            super().__init__()
            self.calls = 0

        def calculate(
            self, atoms=None, properties=("energy", "forces"), system_changes=all_changes
        ):
            super().calculate(atoms, properties, system_changes)
            self.calls += 1
            if self.calls == 2:
                raise GateRejected("max_force", "probe refused", atoms)
            self.results = {"energy": 0.0, "forces": np.zeros((len(atoms), 3))}

    water = molecule("H2O")
    calculator = FailsOnSecond()
    water.calc = calculator
    guard = GuardedCalculator(calculator, SearchConfig(**CONFIG), 0.0)
    reference = water.positions.copy()
    direction = np.zeros((len(water), 3))
    direction[0] = [1.0, 0.0, 0.0]

    with pytest.raises(GateRejected):
        hessian_vector_product(water, guard, direction, 0.01)
    assert calculator.calls == 2  # it really did reach the second call
    assert np.array_equal(water.positions, reference)


def test_the_rigid_body_floor_passes_an_error_the_internal_recheck_catches():
    """The reason the floor gate is a necessary condition and not a sufficient one.

    An error confined to the internal subspace -- P_internal A P_internal -- has exactly
    zero residual on every rigid-body mode, so the floor gate cannot see it however large
    it is. Here the error is applied to the eigenvalue that decides the count, the floor
    is measured on the true geometry and passes, and the recheck is what notices. Passing
    D1 is evidence about invariance breaking; it was never evidence about curvature.
    """
    from prrs.runner import curvature_spectrum, internal_mode_recheck

    atoms = _pair(1.2)  # a real minimum of the well
    config = SearchConfig(closed_shell_only=False, **CONFIG)
    eigenvalues, vectors, floor, provenance = curvature_spectrum(
        atoms, double_well_factory, config
    )

    # The floor is clean: nothing here breaks translational or rotational invariance.
    assert (
        floor_residual(floor)
        <= config.trivial_floor_fraction * config.minimum_check_eigenvalue_tol
    )

    honest = internal_mode_recheck(
        atoms, double_well_factory, config, eigenvalues, vectors, provenance
    )
    assert honest["performed"] and honest["agree"]
    # The mode it checked is the internal one, not a rigid-body zero. Reading the boundary
    # off `len(negative)` put it on a translation, where the measured curvature is zero,
    # it agrees with the zero eigenvalue and the gate decides nothing -- which on a
    # minimum, having no negative modes, was every mode it looked at.
    assert [m["index"] for m in honest["modes"]] == [honest["boundary_mode_index"]]
    assert all(m["subspace"] == "internal" for m in honest["modes"])
    assert abs(honest["modes"][0]["eigenvalue"]) > 1.0

    # Now claim the deciding mode is unstable when the potential says it is not. The
    # eigenvectors are untouched, so the rigid-body residual is exactly what it was.
    corrupted = np.array(eigenvalues, dtype=float)
    lowest = int(np.argmin(np.abs(corrupted)))
    corrupted[lowest] = -50.0
    caught = internal_mode_recheck(
        atoms, double_well_factory, config, corrupted, vectors, provenance
    )
    assert caught["performed"] and not caught["agree"], caught
    assert any(not m["classification_agrees"] for m in caught["modes"])


def test_the_recheck_gate_separates_assembly_from_resolution():
    """A sign that moves with the step is the surface talking, not a defect.

    The gate compares the assembled eigenvalue against a direct product at the SAME step,
    which is different arithmetic at the same bandwidth. The step scan is reported and
    never gated, because a gaussian ridge narrower than the probe genuinely has a
    curvature whose sign depends on how wide the probe is.
    """
    from prrs.runner import curvature_spectrum, internal_mode_recheck

    atoms = _pair(1.2)
    config = SearchConfig(
        closed_shell_only=False, internal_mode_recheck_step_scale=2.0, **CONFIG
    )
    eigenvalues, vectors, _, provenance = curvature_spectrum(atoms, double_well_factory, config)
    report = internal_mode_recheck(
        atoms, double_well_factory, config, eigenvalues, vectors, provenance
    )
    assert report["gate_step_A"] == config.minimum_check_step_A
    for mode in report["modes"]:
        # The ladder brackets the Hessian's own step, so the extrapolation is an
        # interpolation of the trend rather than a reach past the smallest point.
        assert mode["steps_A"][1] == config.minimum_check_step_A
        assert mode["steps_A"][0] < mode["steps_A"][1] < mode["steps_A"][2]
        assert len(mode["curvature_by_step"]) == 3
    # On a smooth well nothing moves, so the flags are not simply always false.
    assert report["resolution_stable"]
    assert report["resolution_diagnosis"] == "stable"


def test_the_recheck_extrapolates_the_probe_truncation_away():
    """Why the gate compares an extrapolate rather than a single measurement.

    An analytic Hessian carries no truncation and a central difference carries O(h^2), so
    comparing them at one step charges the matrix for the probe's error. On a stiff mode
    that is harmless; on a soft one it is not. Measured on this project's MACE surface, a
    rotor mode whose analytic lambda is -0.00318 differences to -0.00498 at the default
    step -- a 56% overestimate, which near the tolerance is enough to move the
    classification and fail a correct Hessian closed.

    Fitting the h^2 law across a ladder removes it: here the single-step error is 3.4e-5
    and what is left after extrapolating is at machine precision.
    """
    from prrs.runner import GuardedCalculator, _curvature_along, curvature_spectrum

    def factory():
        return WellWithHessian()

    config = SearchConfig(
        closed_shell_only=False,
        hessian_source="analytic",
        internal_mode_recheck=False,
        **CONFIG,
    )
    for distance in (1.3, 1.7, 2.0):
        atoms = Atoms("C2", positions=[[-distance / 2, 0, 0], [distance / 2, 0, 0]])
        atoms.calc = factory()
        eigenvalues, vectors, _, _ = curvature_spectrum(
            atoms, factory, config, source="analytic"
        )
        index = int(np.argmax(np.abs(eigenvalues)))  # the one internal mode
        exact = float(eigenvalues[index])
        mode = vectors[:, index].reshape(-1, 3)

        guard = GuardedCalculator(factory(), config, atoms.get_potential_energy())
        probe = atoms.copy()
        probe.calc = guard
        base = config.minimum_check_step_A
        scale = config.internal_mode_recheck_step_scale
        ladder = np.array([base / scale, base, base * scale])
        measured = np.array([_curvature_along(probe, guard, mode, float(h)) for h in ladder])
        design = np.vstack([np.ones(3), ladder**2]).T
        (intercept, _), *_ = np.linalg.lstsq(design, measured, rcond=None)

        single_step = abs(measured[1] - exact)
        extrapolated = abs(intercept - exact)
        assert single_step > 1e-6, "this geometry has no truncation to remove"
        assert extrapolated < single_step / 1000, (
            f"r={distance}: extrapolation left {extrapolated:.3g}, "
            f"single step was {single_step:.3g}"
        )


def test_a_sign_that_moves_with_the_step_is_reported_not_gated():
    """A ridge narrower than the probe has a curvature whose sign is a bandwidth choice.

    On the bumpy-product surface the deciding quotient runs -30.5, -19.4, +16.8 across a
    four-fold step range. The assembly gate must still pass -- the matrix is right -- and
    the count must survive, with the movement reported instead.
    """
    from test_boundary_continuation import ProtonWithBumpyProduct, shared_proton_atoms

    def factory():
        return ProtonWithBumpyProduct()

    config = SearchConfig(closed_shell_only=False, quench_fmax_eV_A=0.001)
    atoms = shared_proton_atoms(ProtonWithBumpyProduct.bump_at).copy()
    atoms.calc = factory()
    _, diagnostics = confirm_minimum(atoms, factory, config, 0)

    assert diagnostics["saddle_order"] == 1
    report = diagnostics["internal_mode_recheck"]
    assert report["agree"], "the matrix is correct; the gate must not fire"
    assert not report["resolution_stable"]
    assert report["resolution_diagnosis"] == "classification_moves_with_step"
    signs = {np.sign(v) for v in report["modes"][0]["curvature_by_step"]}
    assert len(signs) > 1, "this surface is supposed to change sign across the ladder"


def test_confirm_minimum_records_that_the_floor_gate_claims_only_necessity():
    """The claim in the record has to match the claim the gate can support."""
    atoms = _pair(1.2)
    config = SearchConfig(closed_shell_only=False, **CONFIG)
    _, diagnostics = confirm_minimum(atoms, double_well_factory, config, 0)
    assert diagnostics["floor_gate"] == "necessary_condition_only"
    assert "says nothing about the internal subspace" in diagnostics["floor_gate_meaning"]
    assert diagnostics["internal_mode_recheck"]["performed"]


def test_the_chord_reproduces_the_hessian_quotient_away_from_a_stationary_point():
    """Why the probe map is a chord and not the arc the search actually walks.

    A retraction R(av) = R + a v + a^2 c_v / 2 gives kappa_0 = v^T H v + g . c_v. The
    chord has c_v = 0 exactly, so its limit is the Hessian quotient at ANY geometry --
    not only where the gradient vanishes. A curved map would be measuring something else
    wherever g != 0, which is most of the surface a search walks over.
    """
    from prrs.runner import _curvature_along, _mass_weighted_direction

    config = SearchConfig(closed_shell_only=False, **CONFIG)
    for distance in (1.35, 1.6, 2.1):  # all off the stationary points
        atoms = Atoms("C2", positions=[[-distance / 2, 0, 0], [distance / 2, 0, 0]])
        guard = GuardedCalculator(double_well_factory(), config, None)
        atoms.calc = guard
        gradient = -np.asarray(guard.get_forces(atoms))
        assert np.abs(gradient).max() > 1e-3, "this geometry is supposed to be off-minimum"

        cartesian = np.array([[-1.0, 0, 0], [1.0, 0, 0]])
        direction = _mass_weighted_direction(atoms.get_masses(), cartesian)
        quotient = _curvature_along(atoms, guard, direction, config.minimum_check_step_A)

        report = response_curvature(atoms, guard, cartesian, config)
        # kappa_0 is the a -> 0 extrapolate of the chord probe; it must land on v^T H v.
        assert report["kappa"] == pytest.approx(quotient, rel=2e-3), (
            f"r={distance}: kappa={report['kappa']:.6g} vs v^T H v={quotient:.6g}"
        )


def test_an_odd_second_order_term_in_the_probe_map_would_diverge_as_one_over_a():
    """The property that makes the symmetric difference work, stated as a measurement.

    The first-order term cancels only because the +a and -a endpoints carry the SAME
    second-order term. A chord has none, so it satisfies this trivially -- but any curved
    retraction added here later must be checked, and this is the check. With an odd
    component the numerator keeps a term in a . g . v and the quotient runs away as 1/a
    instead of settling.
    """
    config = SearchConfig(closed_shell_only=False, **CONFIG)
    atoms = Atoms("C2", positions=[[-0.8, 0, 0], [0.8, 0, 0]])
    guard = GuardedCalculator(double_well_factory(), config, None)
    atoms.calc = guard
    reference = atoms.positions.copy()
    reference_energy = float(guard.get_potential_energy(atoms))
    gradient = -np.asarray(guard.get_forces(atoms))
    assert np.abs(gradient).max() > 1e-3

    direction = np.array([[-1.0, 0, 0], [1.0, 0, 0]]) / np.sqrt(2)

    def quotient(amplitude, asymmetry):
        energies = []
        for sign in (1, -1):
            stretch = amplitude * (1.0 + asymmetry if sign > 0 else 1.0)
            atoms.set_positions(reference + sign * stretch * direction)
            energies.append(float(guard.get_potential_energy(atoms)))
        atoms.set_positions(reference)
        return (energies[0] + energies[1] - 2 * reference_energy) / amplitude**2

    amplitudes = np.array([0.02, 0.01, 0.005, 0.0025, 0.00125])
    even = np.array([quotient(a, 0.0) for a in amplitudes])
    odd = np.array([quotient(a, 0.05) for a in amplitudes])

    # Even: settles to a finite limit as the amplitude shrinks.
    assert np.ptp(even) < 0.05 * abs(even.mean())
    # Odd: runs away instead, and by the smallest amplitude has left the true value far
    # behind rather than converging on it.
    assert abs(odd[-1]) > 5 * abs(even.mean())

    # The difference isolates the term the cancellation was supposed to remove, and it is
    # exactly first order in 1/a: halving the amplitude doubles it. Measured here as
    # 2.13, 2.06, 2.03, 2.01 -- approaching 2 from above as the higher orders die out.
    divergent = odd - even
    ratios = divergent[1:] / divergent[:-1]
    assert np.all(ratios > 1.9) and np.all(ratios < 2.3), ratios
    assert ratios[-1] == pytest.approx(2.0, abs=0.1), ratios[-1]


def test_the_descent_residual_is_measured_in_the_subspace_it_can_move_in():
    """The walk must not be charged for a component it has deliberately excluded.

    The step is the tangential gradient with the rigid-body subspace projected out, so the
    convergence residual has to be the norm of that same projected vector. Measured before
    the projection it includes a part no step can ever remove, and the ratio it forms with
    an unprojected denominator cannot reach the tolerance however good the direction is --
    so the walk either runs to its step limit or reports a stall it is not in.
    """
    from prrs.runner import _projected_direction, _response_and_gradient, _without_trivial_modes

    config = SearchConfig(closed_shell_only=False, **CONFIG)
    atoms = _pair(1.7)
    guard = GuardedCalculator(double_well_factory(), config, None)
    atoms.calc = guard
    masses = atoms.get_masses()
    amplitude = float(config.response_band_A[0])

    direction = _projected_direction(
        atoms.positions, masses, np.array([[-1.0, 0, 0], [1.0, 0, 0]])
    )
    _, gradient = _response_and_gradient(
        atoms, guard, direction, amplitude, float(guard.get_potential_energy(atoms))
    )
    # Add a pure mass-weighted translation: a direction the projection removes entirely,
    # so it can change the residual only if the residual is measured before projecting.
    root = np.sqrt(masses[:, None])
    translation = root * np.tile([0.0, 1.0, 0.0], (len(atoms), 1))
    translation /= np.linalg.norm(translation)
    polluted = gradient + 5.0 * np.linalg.norm(gradient) * translation

    def tangential(g):
        return g - np.sum(g * direction) * direction

    before = float(np.linalg.norm(tangential(polluted)))
    after = float(
        np.linalg.norm(_without_trivial_modes(atoms.positions, masses, tangential(polluted)))
    )
    honest = float(
        np.linalg.norm(_without_trivial_modes(atoms.positions, masses, tangential(gradient)))
    )

    assert before > 3 * after, "the rigid-body part must dominate the naive residual"
    assert after == pytest.approx(honest, rel=1e-9), (
        "projecting must leave the physical residual untouched"
    )


def test_the_persistence_band_continues_instead_of_restarting():
    """ "Tracking" has to mean the next amplitude starts from the last answer.

    Restarting from the original seed at every amplitude makes
    rotation_from_previous_deg the distance between two independent solutions, which is
    exactly the quantity that cannot tell a direction that moved from a direction that
    jumped to another branch. The band is a path, and the record says so.
    """
    from prrs.runner import response_persistence

    config = SearchConfig(closed_shell_only=False, response_samples=4, **CONFIG)
    atoms = _pair(1.7)
    guard = GuardedCalculator(double_well_factory(), config, None)
    atoms.calc = guard

    report = response_persistence(atoms, guard, config, np.array([[-1.0, 0, 0], [1.0, 0, 0]]))
    assert report["warm_started"] is True
    assert "the result depends on the order" in report["warm_start_meaning"]
    assert len(report["band"]) == 4
    assert report["band"][0]["rotation_from_previous_deg"] is None
    for entry in report["band"][1:]:
        assert entry["rotation_from_previous_deg"] is not None
    # Drift from the original seed is reported alongside, so a steady walk is separable
    # from jitter of the same total size.
    for entry in report["band"]:
        assert entry["rotation_from_seed_deg"] is not None
        assert "stalled" in entry


def test_the_response_descent_reading_is_seeded_independently_of_the_hessian():
    """Corroboration needs the second reading not to have been told the answer.

    `response_curvature` is handed the Hessian's own eigenvector, so it refines that
    answer; a wrong Hessian makes both wrong together. This descent starts from a seeded
    random direction, so where it lands on the same direction the agreement means
    something. The overlap is the number, and it is reported rather than gated.
    """

    atoms = _pair(1.8)  # the double well's barrier top
    config = SearchConfig(
        closed_shell_only=False,
        response_descent_enabled=True,
        response_descent_max_steps=8,
        **CONFIG,
    )
    confirmed, diagnostics = confirm_minimum(atoms, double_well_factory, config, 0)
    assert not confirmed and diagnostics["saddle_order"] == 1, diagnostics

    reading = diagnostics["response_descent"]
    assert reading["seeded_from"] == "config.seed random direction"
    assert reading["overlap_with_unstable_mode"] is not None
    # One internal coordinate on a diatomic, so an honest descent has nowhere else to go.
    assert reading["overlap_with_unstable_mode"] > 0.99, reading
    assert reading["angle_from_unstable_mode_deg"] < 10.0
    assert reading["force_evaluations"] > 0


def test_the_independent_reading_is_off_unless_asked_for():
    """It costs evaluations per saddle and no recorded result was produced with it."""
    atoms = _pair(1.8)
    config = SearchConfig(closed_shell_only=False, **CONFIG)
    assert config.response_descent_enabled is False
    _, diagnostics = confirm_minimum(atoms, double_well_factory, config, 0)
    assert "response_descent" not in diagnostics
    assert "response" in diagnostics  # the band probe still runs


def test_the_recheck_boundary_mode_is_internal_and_not_a_rigid_body_zero():
    """Where the boundary of the count actually lives in a projected spectrum.

    `curvature_spectrum` returns all 3N eigenvalues of the projected matrix in ascending
    order, so the rigid-body modes are exact zeros sitting between the negative internal
    modes and the positive ones. Counting `len(negative)` therefore lands on the first
    rigid-body mode every time, never on the first internal one -- and a curvature
    measured along a translation is zero, agrees with the zero eigenvalue, and decides
    nothing. Checked on both verdicts, because on a minimum that WAS the whole gate.
    """
    from test_convergence import Chain, _chain

    from prrs.runner import curvature_spectrum, internal_mode_recheck

    config = SearchConfig(closed_shell_only=False, **CONFIG)
    tolerance = config.minimum_check_eigenvalue_tol
    for phi, expected_negative in ((0.0, 0), (60.0, 1)):
        atoms = _chain(phi)
        eigenvalues, vectors, _, provenance = curvature_spectrum(
            atoms, Chain, config, source="fd"
        )
        assert sum(1 for v in eigenvalues if v < -tolerance) == expected_negative
        report = internal_mode_recheck(atoms, Chain, config, eigenvalues, vectors, provenance)
        boundary = report["boundary_mode_index"]
        # It is the softest INTERNAL non-negative mode, which on a C4 chain sits at index
        # 6 or 7 -- after all six rigid-body zeros, not at index 0 or 1.
        assert boundary >= 6, (phi, boundary, list(eigenvalues))
        assert eigenvalues[boundary] > tolerance
        assert boundary in [mode["index"] for mode in report["modes"]]
        for mode in report["modes"]:
            if mode["eigenvalue"] >= -tolerance:
                assert mode["subspace"] == "internal", mode


def test_the_recheck_compares_like_with_like_for_each_hessian_source():
    """The reference and the probe have to carry the same truncation, or neither number
    is about the matrix.

    An analytic Hessian has none, so the gate must compare the h -> 0 extrapolate. An fd
    Hessian is itself a central difference at `step_A`, so its eigenvalue carries exactly
    the truncation the extrapolate removes and the gate must compare the rung at that same
    step. Fixing one direction and not the other is the same error twice.
    """
    from prrs.runner import curvature_spectrum, internal_mode_recheck

    config = SearchConfig(closed_shell_only=False, **CONFIG)
    atoms = _pair(1.2)
    eigenvalues, vectors, _, provenance = curvature_spectrum(
        atoms, double_well_factory, config, source="fd"
    )
    fd = internal_mode_recheck(
        atoms, double_well_factory, config, eigenvalues, vectors, provenance
    )
    assert fd["hessian_source"] == "fd"
    assert fd["gate_step_A"] == provenance["step_A"]
    for mode in fd["modes"]:
        assert mode["compared_against"] == "assembly_step"
        assert mode["reference_curvature"] == mode["curvature_by_step"][1]
        assert mode["steps_A"][1] == provenance["step_A"]

    def analytic_factory():
        return WellWithHessian()

    atoms = Atoms("C2", positions=[[-0.6, 0, 0], [0.6, 0, 0]])
    atoms.calc = analytic_factory()
    eigenvalues, vectors, _, provenance = curvature_spectrum(
        atoms, analytic_factory, config, source="analytic"
    )
    analytic = internal_mode_recheck(
        atoms, analytic_factory, config, eigenvalues, vectors, provenance
    )
    assert analytic["hessian_source"] == "analytic"
    for mode in analytic["modes"]:
        assert mode["compared_against"] == "extrapolate"
        assert mode["reference_curvature"] == mode["extrapolated_to_zero_step"]


def test_an_fd_boundary_mode_that_is_all_truncation_still_passes_the_gate():
    """The false positive that comparing an extrapolate against an fd matrix produces.

    On the bumpy-product surface the softest internal mode's ladder is pure h^2 with a
    zero intercept: the whole of the assembled +5e-3 IS the truncation of the differences
    it was built from. Extrapolating the probe to h -> 0 then reports a flat mode against
    a matrix that says positive, and the gate fails a Hessian that is exactly what it
    claims to be. Comparing at the assembly step compares the two arithmetics at one
    bandwidth, which is the only thing this gate was ever able to decide.

    This mode is only reachable at all because the boundary is now chosen in the internal
    subspace; before that the gate looked at a rigid-body zero and never saw it.
    """
    from test_boundary_continuation import ProtonWithBumpyProduct, shared_proton_atoms

    from prrs.runner import curvature_spectrum, internal_mode_recheck, quench

    def factory():
        return ProtonWithBumpyProduct()

    config = SearchConfig(closed_shell_only=False, **CONFIG)
    atoms = shared_proton_atoms(ProtonWithBumpyProduct.bump_at).copy()
    atoms.calc = factory()
    quench(atoms, GuardedCalculator(factory(), config, atoms.get_potential_energy()), config)
    eigenvalues, vectors, _, provenance = curvature_spectrum(
        atoms, factory, config, source=config.hessian_source
    )
    report = internal_mode_recheck(atoms, factory, config, eigenvalues, vectors, provenance)
    boundary = next(m for m in report["modes"] if m["index"] == report["boundary_mode_index"])

    # Pure truncation: the intercept is four decimal orders below the eigenvalue it is
    # meant to reproduce (5.9e-7 against 4.99e-3), and the ladder quadruples with the step.
    assert boundary["fits_h2_law"]
    assert abs(boundary["extrapolated_to_zero_step"]) < 1e-3 * boundary["eigenvalue"]
    ladder = boundary["curvature_by_step"]
    assert 3.5 < ladder[1] / ladder[0] < 4.5 and 3.5 < ladder[2] / ladder[1] < 4.5
    # So the extrapolate would have classified it flat against a positive eigenvalue...
    assert boundary["extrapolated_to_zero_step"] < config.minimum_check_eigenvalue_tol
    assert boundary["eigenvalue"] > config.minimum_check_eigenvalue_tol
    # ...and the same-step rung agrees with the matrix to a few percent.
    assert boundary["compared_against"] == "assembly_step"
    assert abs(boundary["reference_curvature"] / boundary["eigenvalue"] - 1.0) < 0.15
    assert boundary["classification_agrees"] and report["agree"]


def test_the_seed_drift_is_not_a_copy_of_the_step_rotation():
    """`rotation_from_seed_deg` has to be measured against the seed the band started at.

    Reading it off `lowest_response_direction`'s `rotation_from_initial_deg` did not do
    that once the band was warm-started: from the second amplitude onward that function's
    own "initial" IS the previous answer, so the field was a bit-identical copy of
    `rotation_from_previous_deg` and the drift was recorded nowhere.
    """
    from test_convergence import Chain, _chain

    from prrs.runner import response_persistence

    atoms = _chain(30.0)
    config = SearchConfig(closed_shell_only=False, response_samples=4, **CONFIG)
    guard = GuardedCalculator(Chain(0.05), config, atoms.get_potential_energy())
    initial = np.zeros((len(atoms), 3))
    initial[3] = [0.0, 0.0, 1.0]
    band = response_persistence(atoms, guard, config, initial)["band"]

    # This seed rotates tens of degrees on the first amplitude and then barely moves, so a
    # field that tracks the seed and a field that tracks the previous step cannot be
    # confused for one another here.
    assert band[0]["rotation_from_previous_deg"] is None
    assert band[0]["rotation_from_seed_deg"] > 10.0
    for entry in band[1:]:
        assert entry["rotation_from_seed_deg"] != entry["rotation_from_previous_deg"]
        assert entry["rotation_from_seed_deg"] > 10 * entry["rotation_from_previous_deg"]


def test_a_gate_refusal_on_the_final_hessian_is_a_return_and_not_an_exception():
    """Every other failure in `follow_min_mode` is a return carrying a `reason`.

    The final spectrum sat outside the walk's try, and it builds a GuardedCalculator and
    displaces along every axis just as the walk does. Nothing above catches GateRejected
    -- search.py calls follow_min_mode bare, and seek_saddle is called bare in turn -- so
    a refusal there took the whole run down instead of recording one dead climb.
    """
    from prrs.runner import follow_min_mode

    class HessianThatGoesBadAfterTheFirstLook(DoubleWell):
        """Correct analytic Hessian once, then one the finite-difference check rejects.

        `_analytic_hessian` returns no matrix when its seeded random HVPs disagree, and
        with hessian_source="analytic" that is a GateRejected -- raised, here, only on the
        walk's final spectrum.
        """

        looks = 0

        def get_hessian(self, atoms):
            type(self).looks += 1
            exact = _analytic_double_well_hessian(atoms)
            return exact if type(self).looks <= 1 else 3.0 * exact

    def factory():
        return HessianThatGoesBadAfterTheFirstLook()

    # r = 1.8 is the barrier top, so the walk is converged at step 0: one spectrum inside
    # the try, then the final one, which is the call that goes bad.
    atoms = Atoms("C2", positions=[[-0.9, 0, 0], [0.9, 0, 0]])
    atoms.calc = factory()
    config = SearchConfig(
        closed_shell_only=False, hessian_source="analytic", min_mode_refresh=1000, **CONFIG
    )
    structure, diagnostics = follow_min_mode(atoms, factory, config)
    assert HessianThatGoesBadAfterTheFirstLook.looks >= 2
    assert structure is None
    assert diagnostics["reason"] == "gate_analytic_hessian_unavailable", diagnostics
    assert diagnostics["stage"] == "final_spectrum"
    # The bounds every other early return carries are on this one too.
    assert "max_rise_eV" in diagnostics and "backtracks" in diagnostics


def test_a_failed_descent_keeps_the_polish_rounds_that_explain_it():
    """The failure path has to carry what the success path carries.

    `follow_unstable_mode` stored the quench diagnostics on the endpoint when it worked and
    returned `{"reason": ...}` when it did not, and `descend_saddle` then flattened them a
    second time -- so `soft_mode_moved_on_the_last_round`, which is every one of this
    project's six failed proton-transfer descents, arrived with no coordinate, no offset and
    no round history. A limit cycle between FIRE and the torsion polish and a budget one
    round short read identically through that hole.
    """
    from test_convergence import Chain, _chain

    from prrs.runner import polish_soft_modes, quench

    atoms = _chain(30.0)
    config = SearchConfig(closed_shell_only=False, quench_fmax_eV_A=0.001, soft_polish_rounds=1)
    guard = GuardedCalculator(Chain(0.05), config, atoms.get_potential_energy())
    converged, diagnostics = quench(atoms, guard, config)

    # Whatever the verdict, the rounds travel with it and each one says what it moved.
    assert diagnostics["rounds"], diagnostics
    for report in diagnostics["rounds"]:
        assert "movements" in report
        assert report["moved"] == bool(report["movements"])
        for movement in report["movements"]:
            assert len(movement["indices"]) == 4
            assert movement["kind"] in ("scan", "newton")
            assert isinstance(movement["offset_rad"], float)

    # And a polish that steps records the step rather than only the flag.
    report = polish_soft_modes(_chain(30.0), guard, config)
    assert report["moved"] == bool(report["movements"])


def test_a_recheck_that_could_not_run_fails_the_confirmation_closed():
    """The gate that certifies the count must not be skippable by refusing its probe.

    `internal_mode_recheck` has two ways to come back with `performed: False`. One is
    "no modes to check" -- no negative eigenvalue, no internal mode near the tolerance,
    nothing to certify -- and that is legitimately open. The other is a refused force
    evaluation, where there IS a classification to certify and the measurement could not
    be made. The caller used to test only `performed and not agree`, so the second case
    confirmed the count uncertified: a saddle order minted by a gate that never ran.
    That is the opposite direction from every other gate in `confirm_minimum`, including
    the rejection of `curvature_spectrum` itself a few lines above it.

    The refusal here is a real one -- the physical calculator returns a force above the
    configured limit and `GuardedCalculator` rejects it -- so the propagation path is the
    one the search actually uses, not a raise planted in the middle of it.
    """

    class Overforced(DoubleWell):
        """Physics the guard has to refuse: forces far above the configured ceiling."""

        def calculate(
            self, atoms=None, properties=("energy", "forces"), system_changes=all_changes
        ):
            super().calculate(atoms, properties, system_changes)
            self.results["forces"] = np.full((len(atoms), 3), 1e6)

    atoms = _pair(1.2)  # a real minimum of the well
    config = SearchConfig(closed_shell_only=False, **CONFIG)

    # Calibrated rather than assumed: whatever the spectrum consumes is left honest, and
    # the calculator the recheck builds for itself is the poisoned one. Hard-coding a
    # count here would silently start poisoning the assembly instead the day the Hessian
    # changes how many calculators it makes.
    consumed = {"n": 0}

    def counting_factory():
        consumed["n"] += 1
        return double_well_factory()

    curvature_spectrum(atoms, counting_factory, config)
    spectrum_calculators = consumed["n"]
    assert spectrum_calculators >= 1

    built = {"n": 0}

    def poisoned_factory():
        built["n"] += 1
        return double_well_factory() if built["n"] <= spectrum_calculators else Overforced()

    confirmed, diagnostics = confirm_minimum(atoms, poisoned_factory, config, 0)

    # The spectrum itself got through -- so this is the recheck failing, not the matrix.
    assert built["n"] > spectrum_calculators, "the recheck never built its calculator"
    report = diagnostics["internal_mode_recheck"]
    assert report["performed"] is False
    assert report["blocked"] is True
    assert report["reason"].startswith("gate_")

    assert confirmed is False, "a count no gate could certify must not be confirmed"
    assert diagnostics["reason"] == "internal_mode_recheck_blocked"
    assert diagnostics["saddle_order"] is None

    # And the open case stays open: the same surface with the gate able to run confirms,
    # so this fix cannot be passing by rejecting everything.
    ok, clean = confirm_minimum(atoms, double_well_factory, config, 0)
    assert ok and clean["internal_mode_recheck"]["performed"]


def test_the_two_sides_of_a_saddle_get_distinct_non_negative_streams():
    """The seed handed to each descent must be usable and must not be shared.

    The defect: `seed + 7919 * sign` goes negative on the sign = -1 side whenever the
    base seed is below 7919, and `numpy.random.default_rng` raises "expected non-negative
    integer" on it. It reaches that call only through `_probe_minimum`, so every frozen
    run missed it -- they all set minimum_check = "hessian" and `confirm_minimum` returns
    before the seed is read. Under minimum_check = "probe" it was unavoidable at the
    default cfg.seed = 42: the saddle search derives `cfg.seed + saddle_searches`, which
    never reaches 7919, so the negative side crashed every time.

    abs() is not the fix. At seed = 0 both sides map to 7919 and the two descents would
    silently share a stream, which is the one thing the offset exists to prevent.
    """
    from prrs.runner import side_seed

    for base in (0, 1, 42, 43, 7918, 7919, 7920, 10**9):
        plus, minus = side_seed(base, 1), side_seed(base, -1)
        assert plus >= 0 and minus >= 0, (base, plus, minus)
        assert plus != minus, f"both sides share a stream at seed={base}"
        # numpy is the actual consumer; make it prove it accepts them.
        np.random.default_rng(plus)
        np.random.default_rng(minus)

    # Deterministic: the same (seed, sign) must give the same stream every time.
    assert side_seed(17, 1) == side_seed(17, 1)
    assert side_seed(17, -1) == side_seed(17, -1)
