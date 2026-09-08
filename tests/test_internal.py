"""Internal-coordinate perturbations: conventions, conservation and exactness."""

import numpy as np
import pytest
from ase import Atoms
from ase.build import molecule
from prrs import internal


@pytest.fixture
def ethanol():
    atoms = molecule("CH3CH2OH")
    rng = np.random.default_rng(3)
    atoms.set_momenta(rng.normal(scale=0.5, size=(len(atoms), 3)))
    return atoms


def test_coordinates_match_ase_conventions(ethanol):
    assert internal.coordinate(ethanol.positions, "bond", (0, 1)) == pytest.approx(
        ethanol.get_distance(0, 1)
    )
    assert np.degrees(
        internal.coordinate(ethanol.positions, "angle", (0, 1, 2))
    ) == pytest.approx(ethanol.get_angle(0, 1, 2))
    ase_dihedral = ethanol.get_dihedral(0, 1, 2, 3)
    ours = np.degrees(internal.coordinate(ethanol.positions, "dihedral", (0, 1, 2, 3))) % 360
    assert ours == pytest.approx(ase_dihedral % 360, abs=1e-6)


@pytest.mark.parametrize(
    "kind,indices", [("bond", (0, 1)), ("angle", (0, 1, 2)), ("dihedral", (0, 1, 2, 3))]
)
def test_gradient_is_invariant_to_rigid_motion(ethanol, kind, indices):
    """These two identities are what make the perturbations conserving by construction."""
    g = internal.gradient(ethanol.positions, kind, indices)
    assert np.abs(g.sum(axis=0)).max() < 1e-7
    torque = np.cross(ethanol.positions, g).sum(axis=0)
    assert np.abs(torque).max() < 1e-7


@pytest.mark.parametrize(
    "kind,indices,delta",
    [("bond", (0, 1), 0.15), ("angle", (0, 1, 2), 0.10), ("dihedral", (0, 1, 2, 3), 0.30)],
)
def test_displacement_hits_its_target_and_keeps_the_centre_of_mass(
    ethanol, kind, indices, delta
):
    before = ethanol.get_center_of_mass()
    _, achieved = internal.displace(ethanol, kind, indices, delta)
    assert achieved == pytest.approx(delta, rel=1e-6)
    assert np.abs(ethanol.get_center_of_mass() - before).max() < 1e-10


@pytest.mark.parametrize(
    "kind,indices", [("bond", (0, 1)), ("angle", (0, 1, 2)), ("dihedral", (0, 1, 2, 3))]
)
@pytest.mark.parametrize("sign", [1, -1])
def test_kick_is_exact_and_conserves_momentum_and_angular_momentum(
    ethanol, kind, indices, sign
):
    positions = ethanol.positions.copy()
    kinetic = ethanol.get_kinetic_energy()
    momentum = ethanol.get_momenta().sum(axis=0)
    centre = ethanol.get_center_of_mass()
    angular = np.cross(ethanol.positions - centre, ethanol.get_momenta()).sum(axis=0)

    internal.kick(ethanol, kind, indices, 0.4, sign)

    assert ethanol.get_kinetic_energy() - kinetic == pytest.approx(0.4, rel=1e-9)
    assert np.abs(ethanol.get_momenta().sum(axis=0) - momentum).max() < 1e-10
    after = np.cross(ethanol.positions - centre, ethanol.get_momenta()).sum(axis=0)
    assert np.abs(after - angular).max() < 1e-9
    # A kick changes velocities only; the geometry must be untouched.
    assert np.abs(ethanol.positions - positions).max() == 0.0


def test_kick_drives_the_intended_coordinate(ethanol):
    """Sign convention: a positive dihedral kick must increase the dihedral."""
    from prrs.calculators import DoubleWell  # noqa: F401  (import kept local)

    ethanol.set_momenta(np.zeros((len(ethanol), 3)))
    before = internal.coordinate(ethanol.positions, "dihedral", (0, 1, 2, 3))
    internal.kick(ethanol, "dihedral", (0, 1, 2, 3), 0.2, 1)
    velocities = ethanol.get_velocities()
    rate = float(
        np.sum(internal.gradient(ethanol.positions, "dihedral", (0, 1, 2, 3)) * velocities)
    )
    assert rate > 0
    assert before == internal.coordinate(ethanol.positions, "dihedral", (0, 1, 2, 3))


def test_a_dihedral_inside_a_ring_works():
    """The two-rigid-group construction is undefined here; the gradient one is not."""
    ring = molecule("C6H6")
    carbons = [i for i, s in enumerate(ring.get_chemical_symbols()) if s == "C"]
    indices = tuple(carbons[:4])
    g = internal.gradient(ring.positions, "dihedral", indices)
    assert np.abs(g.sum(axis=0)).max() < 1e-7
    assert np.abs(np.cross(ring.positions, g).sum(axis=0)).max() < 1e-7
    before = ring.get_center_of_mass()
    _, achieved = internal.displace(ring, "dihedral", indices, 0.05)
    assert achieved == pytest.approx(0.05, rel=1e-6)
    assert np.abs(ring.get_center_of_mass() - before).max() < 1e-10


def test_degenerate_geometry_is_refused():
    linear = Atoms("H3", positions=[[0, 0, 0], [1, 0, 0], [2, 0, 0]])
    with pytest.raises(ValueError):
        internal.gradient(linear.positions, "angle", (0, 1, 2))
    with pytest.raises(ValueError):
        internal.coordinate(linear.positions, "bond", (0, 0))


def _worst_internal_change(before, after, edges):
    from prrs import internal as ic

    bonds = max(
        (abs(after.get_distance(a, b) - before.get_distance(a, b)) for a, b in edges),
        default=0.0,
    )
    adjacency = {}
    for a, b in edges:
        adjacency.setdefault(a, []).append(b)
        adjacency.setdefault(b, []).append(a)
    angles = 0.0
    for vertex, neighbours in adjacency.items():
        for i in range(len(neighbours)):
            for j in range(i + 1, len(neighbours)):
                triple = (neighbours[i], vertex, neighbours[j])
                angles = max(
                    angles,
                    abs(
                        ic.coordinate(after.positions, "angle", triple)
                        - ic.coordinate(before.positions, "angle", triple)
                    ),
                )
    return bonds, angles


def test_rigid_rotation_preserves_every_bond_and_angle_exactly():
    """The exact realization for a torsion about a bridge bond, at any amplitude."""
    from prrs.perturbations import Probe, apply
    from prrs.state import encode

    source = molecule("CH3CH2OH")
    edges = sorted(encode(source).edges)
    for amplitude in (0.05, 0.52, 2.09, 3.0):
        for sign in (1, -1):
            atoms = source.copy()
            delivered = apply(atoms, Probe("torsion", (0, 1, 2, 3), sign, amplitude, 0))
            assert delivered["realization"] == "rigid_rotation"
            assert delivered["target_error"] < 1e-12
            bonds, angles = _worst_internal_change(source, atoms, edges)
            assert bonds < 1e-12, (amplitude, sign, bonds)
            assert angles < 1e-12, (amplitude, sign, angles)


def test_ring_torsion_falls_back_to_continuation_with_bounded_damage():
    """No rigid fragment exists inside a ring, so the incremental scheme has to hold up."""
    from prrs.perturbations import Probe, apply, torsion_for
    from prrs.state import encode

    source = molecule("C6H6")
    carbons = [i for i, s in enumerate(source.get_chemical_symbols()) if s == "C"]
    indices = tuple(carbons[:4])
    torsion = torsion_for(source, indices)
    assert torsion is None or not torsion.moving  # a ring bond is not a bridge
    atoms = source.copy()
    delivered = apply(atoms, Probe("torsion", indices, 1, 0.25, 0))
    assert delivered["realization"] == "continuation"
    assert delivered["target_error"] < 1e-6
    bonds, _ = _worst_internal_change(source, atoms, sorted(encode(source).edges))
    assert bonds < 0.05, bonds


def _torn_by_one_long_step(monkeypatch, source, indices=(0, 1, 2, 3), delta=2.09):
    """Reproduce the defect the increment cap prevents: one step along the start tangent."""
    from prrs import internal as ic

    monkeypatch.setitem(ic.MAX_STEP, "dihedral", 10.0)
    atoms = source.copy()
    _, achieved = ic.displace(atoms, "dihedral", indices, delta)
    return atoms, achieved


def test_one_long_step_would_tear_the_molecule_and_the_metrics_say_so(monkeypatch):
    """Target fidelity alone reports success on a structure that has been destroyed.

    With the increment cap lifted, a 2.09 rad torsion is delivered along the tangent at
    the starting geometry -- the chord of the arc -- so the requested dihedral is reached
    while the hydroxyl bond is pulled apart.
    """
    from prrs.perturbations import collateral_metrics

    source = molecule("CH3CH2OH")
    atoms, achieved = _torn_by_one_long_step(monkeypatch, source)
    assert abs(achieved - 2.09) < 1e-6  # the named coordinate is exact
    metrics = collateral_metrics(source, atoms, "dihedral", (0, 1, 2, 3))
    assert metrics["max_bond_change_A"] > 1.0  # and the O-H bond is gone
    assert source.get_distance(2, 3) < 1.0 < atoms.get_distance(2, 3)

    # With the cap in place the same request is clean.
    clean = source.copy()
    from prrs import internal as ic

    ic.MAX_STEP.setdefault("dihedral", 0.05)
    from prrs.perturbations import Probe, apply

    delivered = apply(clean, Probe("torsion", (0, 1, 2, 3), 1, 2.09, 0))
    assert delivered["collateral"]["max_bond_change_A"] < 1e-9


def test_collateral_gate_blocks_damage_but_allows_an_intended_bond_break(monkeypatch):
    from prrs.config import SearchConfig
    from prrs.perturbations import collateral_metrics
    from prrs.runner import GateRejected, gate_collateral
    from prrs import internal as ic

    config = SearchConfig()
    source = molecule("CH3CH2OH")

    torn, _ = _torn_by_one_long_step(monkeypatch, source)
    delivered = {"collateral": collateral_metrics(source, torn, "dihedral", (0, 1, 2, 3))}
    with pytest.raises(GateRejected) as raised:
        gate_collateral(delivered, config, torn)
    assert raised.value.code == "collateral_bond"

    # Stretching the bond a stretch probe targets is the point, not collateral damage.
    stretched = source.copy()
    ic.displace(stretched, "bond", (2, 3), 1.5)
    delivered = {"collateral": collateral_metrics(source, stretched, "bond", (2, 3))}
    metrics = delivered["collateral"]
    assert metrics["fragments_after"] == metrics["fragments_before"]
    gate_collateral(delivered, config, stretched)


def test_symmetry_order_and_fundamental_domain():
    from prrs.chemistry import canonical_labels
    from prrs.perturbations import rotatable_torsions
    from prrs.state import encode

    atoms = molecule("CH3CH2OH")
    graph = encode(atoms)
    labels = canonical_labels(atoms.numbers, graph.edges, graph.index)
    genuine, rotors = rotatable_torsions(graph, labels)
    assert [t.symmetry_order for t in genuine] == [1]
    assert [t.symmetry_order for t in rotors] == [3]
    assert rotors[0].domain == pytest.approx(2 * np.pi / 3)
    assert all(t.moving for t in genuine + rotors)  # both bonds are bridges

    ethane = molecule("C2H6")
    graph = encode(ethane)
    labels = canonical_labels(ethane.numbers, graph.edges, graph.index)
    genuine, rotors = rotatable_torsions(graph, labels)
    # Symmetric at both ends: the redundancy is lcm(3, 3) = 3.
    assert not genuine and [t.symmetry_order for t in rotors] == [3]


def test_rotor_amplitudes_avoid_symmetry_copies_of_zero():
    from prrs.config import SearchConfig
    from prrs.perturbations import propose

    atoms = molecule("CH3CH2OH")
    config = SearchConfig(
        families=("torsion",), max_directions=8, torsion_amplitudes_rad=(0.52, 1.05, 2.09)
    )
    for group in propose(atoms, config, 11):
        amplitudes = [p.amplitude for p in group]
        if group[0].indices == (7, 0, 1, 2):
            # 2.09 rad is a quarter of a degree from a C3 symmetry copy of zero.
            assert 2.09 not in amplitudes and amplitudes == [0.52, 1.05]
        else:
            assert amplitudes == [0.52, 1.05, 2.09]


def test_atom_indices_are_validated_not_coerced(ethanol):
    """int(i) accepts everything and complains about nothing.

    -1 is a perfectly good numpy index onto the last atom and 1.9 truncates to 1, so both
    return a real number for a triple nobody named. Nothing downstream can tell that the
    wrong atom moved, which is what makes silent coercion worse than a crash.
    """
    positions = ethanol.positions
    count = len(ethanol)
    for bad in (-1, 1.9, count, count + 5, "1", None, float("nan")):
        with pytest.raises(ValueError):
            internal.coordinate(positions, "bond", (0, bad))
        with pytest.raises(ValueError):
            internal.gradient(positions, "bond", (0, bad))
    with pytest.raises(ValueError):
        internal.displace(ethanol.copy(), "bond", (0, count), 0.01)
    with pytest.raises(ValueError):
        internal.rotate_fragment(positions, -1, 2, [3], 0.1)


def test_numpy_integers_are_accepted(ethanol):
    """Indices arrive from np.where and argmax on a bond graph, so np.int64 has to pass.

    A `type(i) is int` check would reject every index derived from an array, which is the
    ordinary way a caller finds the donor and acceptor of a transfer coordinate.
    """
    positions = ethanol.positions
    plain = internal.coordinate(positions, "bond", (0, 1))
    numpy_typed = internal.coordinate(positions, "bond", (np.int64(0), np.int64(1)))
    assert numpy_typed == plain
    assert internal.atom_index(np.int64(3), len(ethanol)) == 3
    assert isinstance(internal.atom_index(np.int64(3)), int)


def test_named_values_validates_its_own_indices(ethanol):
    """named_values converts indices itself, so it needs the same guard as _check.

    Its bond_difference branch cannot go through _check -- the four indices share the
    transferring atom by construction -- which is exactly how it came to have a second,
    unguarded conversion.
    """
    positions = ethanol.positions
    specs = [{"name": "q", "kind": "bond_difference", "indices": (2, 8, 0, 8)}]
    assert "q" in internal.named_values(positions, specs)
    for bad in ((0, -1), (0, 1.9), (0, len(ethanol))):
        with pytest.raises(ValueError):
            internal.named_values(positions, [{"name": "b", "kind": "bond", "indices": bad}])


def test_a_probe_index_is_validated_rather_than_truncated():
    """`Probe.__post_init__` used a bare int(), which is silent in both bad directions.

    1.9 becomes 1 and -1 becomes a numpy index onto the last atom, so the probe acts on an
    atom nobody named. The upper bound is not checkable here -- a Probe is built before it
    meets a molecule -- so the type and the sign are what this catches.
    """
    from prrs.perturbations import Probe

    with pytest.raises(ValueError):
        Probe("stretch", (0, 1.9), 1, 0.1, 0)
    with pytest.raises(ValueError):
        Probe("stretch", (0, -1), 1, 0.1, 0)
    with pytest.raises(ValueError):
        Probe("stretch", (0, True), 1, 0.1, 0)
    # numpy integers still pass: bond-graph indices arrive from np.where and argmax.
    assert Probe("stretch", (np.int64(0), np.int64(1)), 1, 0.1, 0).indices == (0, 1)
