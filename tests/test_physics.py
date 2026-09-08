import numpy as np
import pytest
from ase import Atoms
from ase import units
from ase.md.verlet import VelocityVerlet
from prrs.calculators import DoubleWell, Committee, demo_atoms
from prrs.config import SearchConfig
from prrs.perturbations import displace_pair, kick_pair
from prrs.reliability import GuardedCalculator, GateRejected, PairPulse
from prrs.state import aligned_rmsd, encode, validate_atoms


@pytest.mark.parametrize("sign", [-1, 1])
def test_kick_adds_exact_energy_and_preserves_momentum_and_torque(sign):
    atoms = Atoms("HCO", positions=[[0.2, -0.3, 0.5], [1.3, 0.7, -0.2], [2.5, 1.5, 0]])
    atoms.set_momenta([[0.7, -0.3, 0.6], [-0.2, 0.8, 0.2], [0.1, 0.4, -0.1]])
    before_k = atoms.get_kinetic_energy()
    before_p = atoms.get_momenta().sum(axis=0)
    before_l = np.cross(atoms.positions, atoms.get_momenta()).sum(axis=0)
    kick_pair(atoms, (0, 1), 0.75, sign)
    assert atoms.get_kinetic_energy() - before_k == pytest.approx(0.75, abs=1e-12)
    np.testing.assert_allclose(atoms.get_momenta().sum(axis=0), before_p, atol=1e-12)
    np.testing.assert_allclose(
        np.cross(atoms.positions, atoms.get_momenta()).sum(axis=0), before_l, atol=1e-12
    )


def test_geometry_preserves_mass_weighted_com():
    atoms = Atoms("HC", positions=[[0, 0, 0], [1.0, 0, 0]])
    com = atoms.get_center_of_mass()
    displace_pair(atoms, (0, 1), 0.4)
    assert atoms.get_distance(0, 1) == pytest.approx(1.4)
    np.testing.assert_allclose(atoms.get_center_of_mass(), com, atol=1e-14)
    with pytest.raises(ValueError):
        displace_pair(atoms, (0, 1), -1.5)


@pytest.mark.parametrize("distance", [1.05, 1.4, 1.9, 2.6])
def test_analytic_forces_match_finite_differences(distance):
    atoms = demo_atoms()
    atoms.positions[1, 0] = atoms.positions[0, 0] + distance
    atoms.calc = DoubleWell()
    forces = atoms.get_forces()
    delta = 1e-6
    atoms.positions[1, 0] += delta
    plus = atoms.get_potential_energy()
    atoms.positions[1, 0] -= 2 * delta
    minus = atoms.get_potential_energy()
    assert forces[1, 0] == pytest.approx(-(plus - minus) / (2 * delta), abs=1e-7)
    np.testing.assert_allclose(forces.sum(axis=0), 0, atol=1e-14)


def test_nve_drift_converges_with_timestep():
    def drift(dt):
        atoms = demo_atoms()
        atoms.calc = DoubleWell()
        kick_pair(atoms, (0, 1), 0.2)
        initial = atoms.get_total_energy()
        dyn = VelocityVerlet(atoms, timestep=dt * units.fs, logfile=None)
        errors = []
        for _ in range(round(10 / dt)):
            dyn.run(1)
            errors.append(abs(atoms.get_total_energy() - initial))
        return max(errors)

    coarse, fine = drift(0.1), drift(0.05)
    assert fine < coarse * 0.35
    assert fine < 1e-5


def test_pulse_energy_gradient_and_work():
    atoms = demo_atoms()
    base = DoubleWell()
    atoms.calc = PairPulse(base, (0, 1), 1, 0.6, 1.2)
    forces = atoms.get_forces()
    atoms.positions[1, 0] += 1e-6
    plus = atoms.get_potential_energy()
    atoms.positions[1, 0] -= 2e-6
    minus = atoms.get_potential_energy()
    assert forces[1, 0] == pytest.approx(-(plus - minus) / 2e-6, abs=1e-7)
    atoms = demo_atoms()
    atoms.calc = PairPulse(base, (0, 1), 1, 0.6, 1.2)
    start = base.get_potential_energy(atoms) + atoms.get_kinetic_energy()
    initial_r = atoms.get_distance(0, 1)
    dyn = VelocityVerlet(atoms, timestep=0.05 * units.fs, logfile=None)
    dyn.run(100)
    work = 0.6 * (atoms.get_distance(0, 1) - initial_r)
    change = base.get_potential_energy(atoms) + atoms.get_kinetic_energy() - start
    assert change == pytest.approx(work, abs=1e-5)
    np.testing.assert_allclose(atoms.get_momenta().sum(axis=0), 0, atol=1e-12)


def test_rmsd_is_invariant_to_translation_and_proper_rotation():
    a = Atoms("HCO", positions=[[0, 0, 0], [1, 0, 0], [0.2, 1.3, 0.5]])
    b = a.copy()
    b.rotate(67, "z")
    b.translate([13, -7, 9])
    assert aligned_rmsd(a, b) < 1e-12


def test_graph_hysteresis_and_finite_scope():
    a = demo_atoms()
    cutoff = 1.2 * (0.76 + 0.76)
    a.positions[1, 0] = a.positions[0, 0] + cutoff * 1.02
    assert (0, 1) not in encode(a).edges
    assert (0, 1) in encode(a, previous={(0, 1)}).edges
    a.pbc = True
    with pytest.raises(ValueError, match="nonperiodic"):
        validate_atoms(a)


def test_uncertainty_is_not_silently_zero():
    a = demo_atoms()
    a.calc = GuardedCalculator(DoubleWell(), SearchConfig(max_uncertainty_eV_A=0.1))
    with pytest.raises(GateRejected, match="supplies none") as rejected:
        a.get_forces()
    assert rejected.value.code == "uncertainty_unavailable"
    a.calc = GuardedCalculator(
        Committee([DoubleWell(), DoubleWell()]), SearchConfig(max_uncertainty_eV_A=0.1)
    )
    assert np.isfinite(a.get_forces()).all()


@pytest.mark.parametrize(
    "options",
    [
        {"timestep_fs": 0},
        {"response_steps": 1.5},
        {"kick_energies_eV": [1, 0.5]},
        {"max_force_eV_A": float("nan")},
        {"families": ["banana"]},
        {"max_depth": -1},
        {"temperature_K": -1},
        {"seed": True},
    ],
)
def test_invalid_protocols_fail_before_running(options):
    with pytest.raises((ValueError, TypeError)):
        SearchConfig(**options)


def _scored(atoms, energy):
    from ase.calculators.singlepoint import SinglePointCalculator

    result = atoms.copy()
    result.calc = SinglePointCalculator(
        result, energy=energy, forces=np.zeros((len(result), 3))
    )
    return result


def _reactant_with_spectators():
    """Two reactive atoms plus far-away spectators standing in for an environment."""
    return Atoms(
        "C2ArArAr",
        positions=[[-0.6, 0, 0], [0.6, 0, 0], [0, 0, 9.0], [0, 0, 12.0], [0, 0, 15.0]],
    )


def test_active_subset_scopes_graph_and_uses_global_indices():
    atoms = _reactant_with_spectators()
    full = encode(atoms, active=None)
    scoped = encode(atoms, active=(0, 1))
    assert scoped.index.tolist() == [0, 1]
    assert scoped.edges == {(0, 1)}
    assert scoped.distances.shape == (2, 2)
    assert full.distances.shape == (5, 5)
    # Edge labels stay global so records remain interpretable against the structure.
    assert all(max(edge) < len(atoms) for edge in full.edges)


def test_environment_motion_does_not_create_a_new_basin():
    """A subset-scoped basin identity is what keeps environment noise out of the network."""
    from prrs.state import same_basin

    config = SearchConfig(active_atoms=(0, 1))
    a = _reactant_with_spectators()
    b = a.copy()
    b.positions[2:] += np.array([0.0, 0.0, 3.0])  # far beyond basin_rmsd_A
    assert aligned_rmsd(a, b) > config.basin_rmsd_A
    assert aligned_rmsd(a, b, config.active_atoms) == pytest.approx(0.0, abs=1e-12)
    assert same_basin(_scored(a, -1.0), _scored(b, -1.0), config)
    whole = SearchConfig()
    assert not same_basin(_scored(a, -1.0), _scored(b, -1.0), whole)


def test_probes_are_never_proposed_outside_the_active_subset():
    from prrs.perturbations import propose

    atoms = _reactant_with_spectators()
    config = SearchConfig(active_atoms=(0, 1), families=("stretch", "kick"))
    proposed = propose(atoms, config, seed=0)
    assert proposed
    for group in proposed:
        for probe in group:
            assert set(probe.pair) <= {0, 1}
    # Without a subset the spectators are fair game, which is the behaviour being scoped.
    unscoped = propose(atoms, SearchConfig(families=("stretch", "kick")), seed=0)
    assert any(set(p.pair) - {0, 1} for group in unscoped for p in group)


def test_active_subset_must_be_valid():
    with pytest.raises(ValueError):
        SearchConfig(active_atoms=(3,))
    with pytest.raises(ValueError):
        SearchConfig(active_atoms=(1, 1, 2))
    with pytest.raises(ValueError):
        validate_atoms(_reactant_with_spectators(), active=(0, 99))


def test_hessian_check_rejects_a_barrier_top_and_accepts_a_minimum():
    """The analytic double well has a stationary point at its barrier top, r=1.8 A.

    A converged quench cannot tell it from a basin; the projected Hessian can, and the
    sampled-probe fallback is only as good as the directions it happens to draw.
    """
    from ase.calculators.singlepoint import SinglePointCalculator
    from prrs.calculators import double_well_factory
    from prrs.runner import confirm_minimum, hessian_spectrum

    def scored(distance):
        atoms = Atoms("C2", positions=[[-distance / 2, 0, 0], [distance / 2, 0, 0]])
        atoms.calc = double_well_factory()
        result = atoms.copy()
        result.calc = SinglePointCalculator(
            result, energy=atoms.get_potential_energy(), forces=atoms.get_forces()
        )
        return result

    config = SearchConfig(quench_fmax_eV_A=0.001)
    # Two atoms leave a single vibrational degree of freedom after rigid-body projection.
    assert hessian_spectrum(scored(1.2), double_well_factory, config).size == 6

    for distance, expected in ((1.2, True), (2.4, True), (1.8, False)):
        confirmed, diagnostics = confirm_minimum(
            scored(distance), double_well_factory, config, seed=0
        )
        assert confirmed is expected, distance
        assert diagnostics["method"] == "hessian"
        if not expected:
            assert len(diagnostics["negative_eigenvalues"]) == 1
            assert diagnostics["imaginary_wavenumbers_icm"][0] < 0


def test_minimum_check_method_is_selectable_and_recorded():
    from ase.calculators.singlepoint import SinglePointCalculator
    from prrs.calculators import double_well_factory
    from prrs.runner import confirm_minimum

    atoms = Atoms("C2", positions=[[-0.9, 0, 0], [0.9, 0, 0]])
    atoms.calc = double_well_factory()
    top = atoms.copy()
    top.calc = SinglePointCalculator(
        top, energy=atoms.get_potential_energy(), forces=atoms.get_forces()
    )

    disabled = SearchConfig(minimum_check="none", quench_fmax_eV_A=0.001)
    confirmed, diagnostics = confirm_minimum(top, double_well_factory, disabled, seed=0)
    assert confirmed and diagnostics["method"] == "none"

    probed = SearchConfig(minimum_check="probe", quench_fmax_eV_A=0.001)
    _, diagnostics = confirm_minimum(top, double_well_factory, probed, seed=0)
    assert diagnostics["method"] == "probe"
    assert len(diagnostics["probes"]) == probed.minimum_check_probes

    with pytest.raises(ValueError):
        SearchConfig(minimum_check="hessian_maybe")


def _double_well_pair(distance):
    from prrs.calculators import double_well_factory

    atoms = Atoms("C2", positions=[[-distance / 2, 0, 0], [distance / 2, 0, 0]])
    atoms.calc = double_well_factory()
    return atoms


def test_hessian_eigenvalue_matches_the_analytic_double_well():
    """The one vibrational eigenvalue is k/mu, and both are known in closed form.

    E = 0.4 (x^2 - 1)^2 with x = (r - 1.8)/0.6, so d2E/dr2 = (1.6/0.36)(3x^2 - 1) and the
    mass-weighted eigenvalue is that over the reduced mass.
    """
    from prrs.calculators import double_well_factory
    from prrs.runner import hessian_spectrum

    config = SearchConfig(quench_fmax_eV_A=0.001)
    for distance in (1.2, 2.4, 1.8):
        x = (distance - 1.8) / 0.6
        stiffness = (1.6 / 0.36) * (3 * x * x - 1)
        reduced = _double_well_pair(distance).get_masses()[0] / 2
        spectrum = hessian_spectrum(_double_well_pair(distance), double_well_factory, config)

        nontrivial = spectrum[-1] if stiffness > 0 else spectrum[0]
        assert nontrivial == pytest.approx(stiffness / reduced, rel=1e-3), distance


@pytest.mark.parametrize("distance,expect_restoring", [(1.2, True), (2.4, True), (1.8, False)])
def test_response_to_a_small_push_has_the_sign_of_the_curvature(distance, expect_restoring):
    """Push along an eigenvector: a minimum pushes back, a saddle pushes further out.

    This is the whole method at infinitesimal amplitude. F(R + eps v) ~ -eps H v, so for
    an eigenvector the response is -eps lambda v: antiparallel where lambda > 0, parallel
    where lambda < 0. The finite-amplitude search is the same measurement continued.
    """
    from prrs.calculators import double_well_factory
    from prrs.runner import hessian_spectrum

    config = SearchConfig(quench_fmax_eV_A=0.001)
    atoms = _double_well_pair(distance)
    masses = atoms.get_masses()
    eigenvalues, vectors = hessian_spectrum(
        atoms, double_well_factory, config, return_vectors=True
    )
    index = 0 if not expect_restoring else len(eigenvalues) - 1
    value, mode = eigenvalues[index], vectors[:, index].reshape(-1, 3)
    assert bool(value > 0) is expect_restoring

    step = 1e-3
    pushed = atoms.copy()
    pushed.calc = double_well_factory()
    pushed.set_positions(atoms.positions + step * mode / np.sqrt(masses[:, None]))
    weighted_force = pushed.get_forces() / np.sqrt(masses[:, None])
    projection = float(np.sum(weighted_force * mode))
    assert projection == pytest.approx(-step * value, rel=0.02)
    assert bool(projection < 0) is expect_restoring


def test_direction_stiffness_costs_two_force_calls_and_agrees_with_the_hessian():
    from prrs.calculators import double_well_factory
    from prrs.runner import direction_stiffness, hessian_spectrum

    config = SearchConfig(quench_fmax_eV_A=0.001)
    atoms = _double_well_pair(1.2)
    spectrum = hessian_spectrum(atoms, double_well_factory, config)
    along_bond = np.array([[-1.0, 0, 0], [1.0, 0, 0]])
    assert direction_stiffness(atoms, double_well_factory, config, along_bond) == pytest.approx(
        spectrum[-1], rel=1e-3
    )
    saddle = _double_well_pair(1.8)
    assert direction_stiffness(saddle, double_well_factory, config, along_bond) < 0
