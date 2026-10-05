"""T-detour (PLAN_2026-09-28.md T0): the gated wire, where the cheap way to P raises d first.

The property the upper layer's no-pruning guard is tested against (plan section 1 item
4) is that the lowest-barrier route from R to P passes a state FARTHER from the product
than the one before it. These tests hold the testbed to that on its own potential --
rigid scans, so barriers are upper bounds -- and to the properties the plain wire is
held to: exact forces, rigid-motion invariance, sixteen confirmed minima that are
sixteen chemical states, six zero modes.
"""

import heapq
import itertools

import numpy as np
import pytest
from ase.optimize import LBFGS

from prrs.calculators import (
    RELAY_GATE_ATOM,
    RELAY_HEAVY,
    ProtonRelayGated,
    _relay_well_radius,
    proton_relay_gated_atoms,
    proton_relay_gated_factory,
)
from prrs.config import SearchConfig
from prrs.network import Registry
from prrs.product import product_reference
from prrs.runner import confirm_minimum

CONFIG = SearchConfig(closed_shell_only=False)
STATES = [
    (sides, home) for sides in itertools.product((0, 1), repeat=3) for home in (True, False)
]
RADIUS = {s: _relay_well_radius(s) for s in (*RELAY_HEAVY, RELAY_GATE_ATOM)}
CL, BR, GATEKEEPER = 3, 7, 8


def _minimum(sides, home):
    atoms = proton_relay_gated_atoms(sides, home)
    atoms.calc = ProtonRelayGated()
    LBFGS(atoms, logfile=None).run(fmax=1e-4, steps=3000)
    return atoms


@pytest.fixture(scope="module")
def minima():
    return {state: _minimum(*state) for state in STATES}


def _rigid_barrier(atoms, h, a, b, ra, rb):
    """Slide atom h from a's well to b's along a-b with everything else frozen."""
    axis = atoms.positions[b] - atoms.positions[a]
    length = np.linalg.norm(axis)
    axis /= length
    energies = []
    for s in np.linspace(ra, length - rb, 201):
        probe = atoms.copy()
        probe.positions[h] = atoms.positions[a] + s * axis
        probe.calc = ProtonRelayGated()
        energies.append(probe.get_potential_energy())
    return max(energies) - energies[0]


def _barriers(minima):
    edges = {}
    for (sides, home), atoms in minima.items():
        for k in range(3):
            if sides[k] == 0:
                to = (tuple(1 if i == k else v for i, v in enumerate(sides)), home)
                edges[(sides, home), to] = _rigid_barrier(
                    atoms, 4 + k, k, k + 1, RADIUS[RELAY_HEAVY[k]], RADIUS[RELAY_HEAVY[k + 1]]
                )
        a, b = (CL, BR) if home else (BR, CL)
        edges[(sides, home), (sides, not home)] = _rigid_barrier(
            atoms,
            GATEKEEPER,
            a,
            b,
            RADIUS["Cl" if home else "Br"],
            RADIUS["Br" if home else "Cl"],
        )
    return edges


def test_forces_are_the_exact_gradient_and_the_energy_is_invariant():
    atoms = proton_relay_gated_atoms((0, 1, 0), True)
    atoms.positions += 0.05 * np.random.default_rng(0).standard_normal(atoms.positions.shape)
    atoms.calc = ProtonRelayGated()
    forces = atoms.get_forces()
    numeric = np.zeros_like(forces)
    for i, c in itertools.product(range(len(atoms)), range(3)):
        energies = []
        for step in (1e-5, -1e-5):
            probe = atoms.copy()
            probe.positions[i, c] += step
            probe.calc = ProtonRelayGated()
            energies.append(probe.get_potential_energy())
        numeric[i, c] = -(energies[0] - energies[1]) / 2e-5
    assert np.abs(forces - numeric).max() < 1e-7
    moved = atoms.copy()
    moved.rotate(41, (2, -1, 1), center="COM")
    moved.translate((0.3, 1.1, -2.0))
    moved.calc = ProtonRelayGated()
    assert abs(moved.get_potential_energy() - atoms.get_potential_energy()) < 1e-12


def test_sixteen_confirmed_minima_sixteen_chemical_states(minima):
    registry = Registry(CONFIG)
    keys = set()
    for state, atoms in minima.items():
        ok, info = confirm_minimum(atoms, proton_relay_gated_factory, CONFIG, 0)
        assert ok, (state, info.get("reason"))
        keys.add(registry.key_of(atoms)["key"])
    assert len(keys) == len(STATES)


@pytest.mark.parametrize("state", [((0, 0, 0), True), ((1, 1, 1), True)])
def test_the_only_zero_modes_are_rigid_motion(state):
    atoms = _minimum(*state)
    atoms.calc = ProtonRelayGated()
    LBFGS(atoms, logfile=None).run(fmax=1e-6, steps=5000)
    base = atoms.positions.copy()
    hessian = np.zeros((base.size, base.size))
    for i in range(base.size):
        gradients = []
        for step in (1e-4, -1e-4):
            probe = atoms.copy()
            probe.positions = base.copy()
            probe.positions[i // 3, i % 3] += step
            probe.calc = ProtonRelayGated()
            gradients.append(-probe.get_forces().ravel())
        hessian[:, i] = (gradients[0] - gradients[1]) / 2e-4
    eigenvalues = np.linalg.eigvalsh(0.5 * (hessian + hessian.T))
    assert int(np.sum(np.abs(eigenvalues) < 1e-3)) == 6
    assert eigenvalues.min() > -1e-3


def test_the_gate_closes_proton_fives_hop_only_while_the_gatekeeper_is_home(minima):
    edges = _barriers(minima)
    gated = [v for ((s, home), (t, _)), v in edges.items() if s[1] == 0 and t[1] == 1 and home]
    open_ = [
        v for ((s, home), (t, _)), v in edges.items() if s[1] == 0 and t[1] == 1 and not home
    ]
    excursion = [v for (a, b), v in edges.items() if a[0] == b[0]]
    assert min(gated) > 2.0
    assert max(open_) < 1.0
    assert max(excursion) < 1.0


def test_the_cheapest_way_to_the_product_raises_d_on_the_way(minima):
    """The d-neutral first version is why the reactant's group is checked here: with the
    gatekeeper home on C it and proton 5 were exchangeable, and the rise vanished."""
    reference = product_reference(
        minima[((0, 0, 0), True)], minima[((1, 1, 1), True)], Registry(CONFIG)
    )
    assert reference.group_order == 1
    distance = {state: reference.distance(atoms) for state, atoms in minima.items()}
    edges = _barriers(minima)
    start, goal = ((0, 0, 0), True), ((1, 1, 1), True)
    best, previous, heap = {start: 0.0}, {}, [(0.0, start)]
    while heap:
        cost, state = heapq.heappop(heap)
        if cost > best[state]:
            continue
        for (a, b), barrier in edges.items():
            if a == state and max(cost, barrier) < best.get(b, np.inf):
                best[b], previous[b] = max(cost, barrier), a
                heapq.heappush(heap, (best[b], b))
    path = [goal]
    while path[-1] != start:
        path.append(previous[path[-1]])
    profile = [distance[state] for state in reversed(path)]
    assert any(later > earlier for earlier, later in zip(profile, profile[1:])), profile
    assert best[goal] < 1.0, "the detour must be cheaper than forcing the gate"


def test_the_excursion_is_uphill_in_energy_too(minima):
    """Home is the deeper well: every gatekeeper-out state lies above its home twin, so
    the detour costs energy as well as distance, and a search that settles by energy
    alone would never leave the gatekeeper out at the product."""
    for sides in itertools.product((0, 1), repeat=3):
        out = minima[(sides, False)].get_potential_energy()
        home = minima[(sides, True)].get_potential_energy()
        assert out > home + 0.05, sides
