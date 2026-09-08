"""A torsion is free only if it is flat AND stationary.

Measured on 3-oxobutanal round three: the product minimum was reported 9.97 meV above the
same minimum reached with a tight tolerance, and the whole difference sat in one methyl
torsion the polish had labelled "free" -- curvature 0.0041 eV/rad^2 against a 0.01 floor,
initial gradient -0.014 eV/rad, and 0.71 rad from its minimum. Because that same verdict
feeds `free_bonds_of`, the coordinate was also quotiented out of microstate comparison. A
flat coordinate has no gradient; a small local curvature only means the Newton scale is
unusable near an inflection of the rotor barrier.
"""

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from prrs import internal
from prrs.config import SearchConfig
from prrs.network import free_bonds_of
from prrs.reliability import GuardedCalculator
from prrs.runner import polish_soft_modes


class Rotor(Calculator):
    """Ethane-like skeleton whose central torsion carries a threefold barrier.

    Everything except the torsion is held by stiff harmonic terms so the polish has exactly
    one soft coordinate to decide about. `amplitude` sets the barrier height and `phase`
    where its minima sit, so the same potential can be made genuinely flat or flat-looking
    at one particular angle.
    """

    implemented_properties = ["energy", "forces"]

    def __init__(self, amplitude=0.02, phase=0.0, **kwargs):
        super().__init__(**kwargs)
        self.amplitude = amplitude
        self.phase = phase

    def calculate(
        self, atoms=None, properties=("energy", "forces"), system_changes=all_changes
    ):
        super().calculate(atoms, properties, system_changes)
        positions = atoms.positions
        energy, forces = 0.0, np.zeros_like(positions)
        step = 1e-5
        for (i, j), r0 in self.bonds(atoms).items():
            vector = positions[j] - positions[i]
            r = np.linalg.norm(vector)
            energy += 20.0 * (r - r0) ** 2
            # forces[] accumulates dV/dx and is negated once at the end. Getting the
            # two halves of that convention to disagree is invisible whenever the bonds
            # start at equilibrium, which they do here -- so it is spelled out.
            gradient = 40.0 * (r - r0) * vector / r
            forces[i] -= gradient
            forces[j] += gradient
        angle = internal.coordinate(positions, "dihedral", (2, 0, 1, 4))
        energy += self.amplitude * (1 - np.cos(3 * (angle - self.phase)))
        # Torsional force by central difference on the exact coordinate, which keeps the
        # test potential's forces consistent with its energy without hand-derived chain
        # rules that would be their own source of error.
        for atom in range(len(atoms)):
            for axis in range(3):
                shifted = positions.copy()
                shifted[atom, axis] += step
                plus = internal.coordinate(shifted, "dihedral", (2, 0, 1, 4))
                shifted[atom, axis] -= 2 * step
                minus = internal.coordinate(shifted, "dihedral", (2, 0, 1, 4))
                derivative = np.arctan2(np.sin(plus - minus), np.cos(plus - minus)) / (2 * step)
                forces[atom, axis] += (
                    3 * self.amplitude * np.sin(3 * (angle - self.phase)) * derivative
                )
        self.results = {"energy": float(energy), "forces": -forces}

    @staticmethod
    def bonds(atoms):
        pairs = {(0, 1): 1.5, (0, 2): 1.1, (0, 3): 1.1, (1, 4): 1.1, (1, 5): 1.1}
        return pairs


def rotor_atoms(dihedral):
    """Two CH2 groups on a central bond, with the requested torsion between them."""
    positions = [
        [0.0, 0.0, 0.0],
        [1.5, 0.0, 0.0],
        [-0.4, 1.03, 0.0],
        [-0.4, -1.03, 0.0],
        [1.9, 1.03, 0.0],
        [1.9, -1.03, 0.0],
    ]
    atoms = Atoms("C2H4", positions=positions)
    internal.displace(
        atoms,
        "dihedral",
        (2, 0, 1, 4),
        dihedral - internal.coordinate(atoms.positions, "dihedral", (2, 0, 1, 4)),
    )
    return atoms


@pytest.fixture
def config():
    return SearchConfig(
        closed_shell_only=False, soft_polish_steps=10, soft_mode_curvature_floor_eV_rad2=0.01
    )


def polish(atoms, calculator, config):
    guard = GuardedCalculator(calculator, config)
    atoms.calc = guard
    return polish_soft_modes(atoms, guard, config)


def central_mode(report):
    return next(mode for mode in report["modes"] if sorted(mode["indices"][1:3]) == [0, 1])


def test_flat_local_curvature_with_a_live_gradient_is_not_free(config):
    """The failure mode itself. At an inflection of the threefold barrier the second
    derivative passes through zero while the first is at its largest."""
    inflection = np.pi / 6  # cos(3x) inflection: k = 0, |g| maximal
    atoms = rotor_atoms(inflection)
    report = polish(atoms, Rotor(amplitude=0.05), config)
    mode = central_mode(report)
    assert abs(mode["curvature_eV_rad2"]) < config.soft_mode_curvature_floor_eV_rad2
    assert abs(mode["initial_gradient_eV_rad"]) > mode["free_gradient_bound_eV_rad"]
    assert mode["tier"] == "flat_biased", "a coordinate with a gradient is not flat"


def test_a_genuinely_flat_rotor_is_still_free(config):
    """The fix must not turn every soft rotor into work. With no barrier at all the
    gradient is zero to numerical precision and the old verdict is the right one."""
    atoms = rotor_atoms(np.pi / 6)
    report = polish(atoms, Rotor(amplitude=0.0), config)
    mode = central_mode(report)
    assert mode["tier"] == "free"
    assert abs(mode["initial_gradient_eV_rad"]) <= mode["free_gradient_bound_eV_rad"]


def test_the_scan_finds_the_rotor_minimum_the_newton_step_could_not(config):
    """-0.014/0.0041 asked for 3.4 radians; the coordinate was 0.71 away. Scanning its own
    exact path over one symmetry period does not need a local model to be valid."""
    inflection = np.pi / 6
    atoms = rotor_atoms(inflection)
    calculator = Rotor(amplitude=0.05)
    before = atoms.copy()
    before.calc = GuardedCalculator(calculator, config)
    energy_before = before.get_potential_energy()
    report = polish(atoms, calculator, config)
    mode = central_mode(report)
    assert mode["scan_points"] == config.soft_mode_scan_points
    assert mode["scan_energy_gain_eV"] > 0.0
    assert atoms.get_potential_energy() < energy_before
    # The threefold minima sit at 0, +-2pi/3. Landing on one of them is the whole point.
    final = internal.coordinate(atoms.positions, "dihedral", (2, 0, 1, 4))
    distance = min(
        abs((final - centre + np.pi) % (2 * np.pi) - np.pi)
        for centre in (0.0, 2 * np.pi / 3, -2 * np.pi / 3)
    )
    assert distance < 0.1, f"left {distance:.3f} rad from the nearest rotor minimum"


def test_no_step_exceeds_half_the_coordinate_period(config):
    """The cap comes from the coordinate's periodicity, not from a trust radius. Without it
    the floor-clamped Newton step is as unusable as the raw one: on this rotor -0.15/0.01
    asks for 15 radians, which wraps past several symmetry copies of where it started. The
    period is the graph's, not the potential's -- the CH2 ends make it twofold."""
    from dataclasses import replace

    quiet = replace(config, soft_mode_scan_points=0)
    atoms = rotor_atoms(np.pi / 6)
    report = polish(atoms, Rotor(amplitude=0.05), quiet)
    mode = central_mode(report)
    assert mode["tier"] == "flat_biased"
    assert "scan_offset_rad" not in mode  # no scan was run
    assert abs(mode["last_offset_rad"]) <= np.pi / mode["symmetry_order"] + 1e-12
    # Capped stepping alone still improves the coordinate; the scan is what makes reaching
    # a minimum reliable rather than iteration-count dependent.
    assert abs(mode["final_gradient_eV_rad"]) < abs(mode["initial_gradient_eV_rad"])


def test_a_flat_biased_torsion_is_not_quotiented_out_of_microstate_comparison(config):
    """`free_bonds_of` reads this tier. Marking a biased coordinate free removed it from
    the RMSD comparison as well, so two structures differing only there looked identical."""
    atoms = rotor_atoms(np.pi / 6)
    report = polish(atoms, Rotor(amplitude=0.05), config)
    atoms.info["quench"] = {"rounds": [report]}
    assert (0, 1) not in free_bonds_of(atoms)
    flat = rotor_atoms(np.pi / 6)
    flat_report = polish(flat, Rotor(amplitude=0.0), config)
    flat.info["quench"] = {"rounds": [flat_report]}
    assert (0, 1) in free_bonds_of(flat)
