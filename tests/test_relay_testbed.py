"""The connect-mode testbed (PLAN_2026-09-28.md T0): what makes it a usable ground truth.

Each test is one property the draft testbed lacked and a later step depends on: the
mainline must accept its minima (G0 runs the real search on it), the eight placements
must be eight chemical states (d and the queue read them), and a hop's new bond must
not form inside the old well (the progress set of L1 would otherwise be entered by a
thermal fluctuation rather than by crossing the barrier).
"""

import itertools

import numpy as np
import pytest
from ase.optimize import LBFGS

from prrs.calculators import (
    RELAY_HEAVY,
    ProtonRelay,
    proton_relay_atoms,
    proton_relay_factory,
)
from prrs.config import SearchConfig
from prrs.network import Registry
from prrs.product import product_reference
from prrs.runner import confirm_minimum
from prrs.state import encode

PLACEMENTS = list(itertools.product((0, 1), repeat=3))
CONFIG = SearchConfig(closed_shell_only=False)


def _minimum(sides):
    atoms = proton_relay_atoms(sides)
    atoms.calc = ProtonRelay()
    LBFGS(atoms, logfile=None).run(fmax=1e-4, steps=2000)
    return atoms


def _expected_edges(sides):
    n = len(RELAY_HEAVY)
    return {(k + side, n + k) for k, side in enumerate(sides)}


def test_forces_are_the_exact_gradient():
    atoms = proton_relay_atoms((0, 1, 0))
    atoms.positions += 0.05 * np.random.default_rng(0).standard_normal(atoms.positions.shape)
    atoms.calc = ProtonRelay()
    forces = atoms.get_forces()
    numeric = np.zeros_like(forces)
    for i, c in itertools.product(range(len(atoms)), range(3)):
        energies = []
        for step in (1e-5, -1e-5):
            probe = atoms.copy()
            probe.positions[i, c] += step
            probe.calc = ProtonRelay()
            energies.append(probe.get_potential_energy())
        numeric[i, c] = -(energies[0] - energies[1]) / 2e-5
    assert np.abs(forces - numeric).max() < 1e-7


def test_energy_is_invariant_under_rigid_motion():
    atoms = proton_relay_atoms((1, 0, 1))
    atoms.calc = ProtonRelay()
    energy = atoms.get_potential_energy()
    moved = atoms.copy()
    moved.rotate(37, (1, 2, 3), center="COM")
    moved.translate((1.0, -2.0, 0.5))
    moved.calc = ProtonRelay()
    assert abs(moved.get_potential_energy() - energy) < 1e-12
    assert np.abs(atoms.get_forces().sum(axis=0)).max() < 1e-12


def test_every_placement_is_a_confirmed_minimum_and_its_own_chemical_state():
    registry = Registry(CONFIG)
    keys = set()
    for sides in PLACEMENTS:
        atoms = _minimum(sides)
        assert {tuple(sorted(e)) for e in encode(atoms, 1.2).edges} == _expected_edges(sides)
        ok, info = confirm_minimum(atoms, proton_relay_factory, CONFIG, 0)
        assert ok, (sides, info.get("reason"))
        keys.add(registry.key_of(atoms)["key"])
    assert len(keys) == len(PLACEMENTS)


@pytest.mark.parametrize("sides", [(0, 0, 0), (1, 1, 1)])
def test_the_only_zero_modes_are_rigid_motion(sides):
    """Six, not seven: without the held dihedral the X0-X1-X2-X3 torsion is exactly
    flat. confirm_minimum tolerates a flat torsion (the test above stays green without
    it), so this is the test that holds the term."""
    atoms = _minimum(sides)
    atoms.calc = ProtonRelay()
    LBFGS(atoms, logfile=None).run(fmax=1e-6, steps=5000)
    base = atoms.positions.copy()
    size = base.size
    hessian = np.zeros((size, size))
    for i in range(size):
        gradients = []
        for step in (1e-4, -1e-4):
            probe = atoms.copy()
            probe.positions = base.copy()
            probe.positions[i // 3, i % 3] += step
            probe.calc = ProtonRelay()
            gradients.append(-probe.get_forces().ravel())
        hessian[:, i] = (gradients[0] - gradients[1]) / 2e-4
    eigenvalues = np.linalg.eigvalsh(0.5 * (hessian + hessian.T))
    assert int(np.sum(np.abs(eigenvalues) < 1e-3)) == 6
    assert eigenvalues.min() > -1e-3


def test_distance_to_the_product_counts_the_hops_left():
    reference = product_reference(
        proton_relay_atoms((0, 0, 0)), _minimum((1, 1, 1)), Registry(CONFIG)
    )
    for sides in PLACEMENTS:
        assert reference.distance(_minimum(sides)) == 2 * (3 - sum(sides))


@pytest.mark.parametrize("gap", range(3))
def test_a_hops_new_bond_forms_only_past_the_barrier_top(gap):
    """Slide proton `gap` from its left well to its right one along the gap axis; the
    acceptor bond must first appear at or after the energy maximum, never before it."""
    left = _minimum((0, 0, 0))
    n = len(RELAY_HEAVY)
    h = n + gap
    axis = left.positions[gap + 1] - left.positions[gap]
    axis /= np.linalg.norm(axis)
    start = float(np.dot(left.positions[h] - left.positions[gap], axis))
    stop = (
        float(np.linalg.norm(left.positions[gap + 1] - left.positions[gap]))
        - (ProtonRelay().radii[gap + 1])
    )
    energies, bonded_right = [], []
    for s in np.linspace(start, stop, 201):
        probe = left.copy()
        probe.positions[h] = left.positions[gap] + s * axis
        probe.calc = ProtonRelay()
        energies.append(probe.get_potential_energy())
        edges = encode(probe, 1.2).edges
        bonded_right.append(tuple(sorted((gap + 1, h))) in edges)
    top = int(np.argmax(energies))
    assert 0 < top < 200, "the path must have an interior barrier"
    first = bonded_right.index(True)
    assert first >= top, (first, top)
    assert energies[top] - energies[0] > 0.3  # a barrier, not a fluctuation
