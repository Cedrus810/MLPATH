"""Chemical identity must ignore conformation and atom labelling, but keep configuration."""

import numpy as np
import pytest
from ase import Atoms
from ase.build import molecule
from prrs.chemistry import chemical_key, canonical_labels, same_chemistry
from prrs.state import encode


def _ethanol_conformer(dihedral):
    atoms = molecule("CH3CH2OH")
    atoms.set_dihedral(
        0, 1, 2, 3, dihedral, mask=[1 if i == 3 else 0 for i in range(len(atoms))]
    )
    return atoms


def test_conformers_of_one_substance_share_a_chemical_key():
    """gauche+ and gauche- are the same substance; splitting them rebuilds the explosion."""
    anti = _ethanol_conformer(180)
    keys = {chemical_key(_ethanol_conformer(angle))["key"] for angle in (60, 120, 180, 300)}
    assert len(keys) == 1
    assert chemical_key(anti)["key"] in keys
    assert chemical_key(anti)["tetrahedral_parity"] == []  # no real stereocentre
    assert chemical_key(anti)["fragments"] == ["C2H6O"]


def test_a_broken_bond_changes_the_chemical_key():
    atoms = molecule("CH3CH2OH")
    dissociated = atoms.copy()
    dissociated.positions[3] += np.array([0.0, 0.0, 4.0])  # pull the hydroxyl H away
    assert not same_chemistry(atoms, dissociated)
    assert chemical_key(dissociated)["fragments"] == ["C2H5O", "H"]


def test_symmetry_equivalent_products_collapse_to_one_chemistry():
    """Which of three equivalent methyl hydrogens leaves must not create three products."""
    atoms = molecule("CH3CH2OH")
    labels = canonical_labels(atoms.numbers, encode(atoms).edges, np.arange(len(atoms)))
    methyl = [i for i in range(len(atoms)) if atoms.numbers[i] == 1 and labels[i] == labels[6]]
    assert len(methyl) >= 2, "expected equivalent methyl hydrogens"

    products, edge_sets = set(), set()
    for hydrogen in methyl[:2]:
        broken = atoms.copy()
        broken.positions[hydrogen] += 4.0 * np.array([0.3, 0.5, 0.8])
        products.add(chemical_key(broken)["key"])
        edge_sets.add(frozenset(encode(broken).edges))
    assert len(edge_sets) == 2, "the raw graphs really do differ by atom index"
    assert len(products) == 1, "but they are one chemistry"


def _chfclbr(mirror=False):
    directions = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], dtype=float)
    directions /= np.linalg.norm(directions, axis=1)[:, None]
    lengths = [1.10, 1.35, 1.77, 1.94]
    positions = [np.zeros(3)] + [d * r for d, r in zip(directions, lengths)]
    positions = np.array(positions)
    if mirror:
        positions[:, 2] *= -1
    return Atoms("CHFClBr", positions=positions)


def test_configurational_inversion_is_a_different_chemistry():
    left, right = _chfclbr(False), _chfclbr(True)
    assert len(chemical_key(left)["tetrahedral_parity"]) == 1
    assert (
        chemical_key(left)["tetrahedral_parity"][0][1]
        == -chemical_key(right)["tetrahedral_parity"][0][1]
    )
    assert not same_chemistry(left, right)
    # A pure rotation is not an inversion.
    rotated = left.copy()
    rotated.rotate(37, "xyz", rotate_cell=False)
    assert same_chemistry(left, rotated)


def _difluoroethene(cis):
    c1, c2 = np.array([0.0, 0, 0]), np.array([1.33, 0, 0])
    up1, down1 = c1 + [-0.54, 0.935, 0], c1 + [-0.54, -0.935, 0]
    up2, down2 = c2 + [0.54, 0.935, 0], c2 + [0.54, -0.935, 0]
    f2 = up2 if cis else down2
    h2 = down2 if cis else up2
    return Atoms("C2FHFH", positions=[c1, c2, up1, down1, f2, h2])


def test_cis_and_trans_across_a_locked_bond_differ():
    cis, trans = _difluoroethene(True), _difluoroethene(False)
    assert len(chemical_key(cis)["locked_bond_parity"]) == 1
    assert (
        chemical_key(cis)["locked_bond_parity"][0][1]
        == -chemical_key(trans)["locked_bond_parity"][0][1]
    )
    assert not same_chemistry(cis, trans)


def test_rotatable_single_bond_torsion_is_not_recorded_as_stereochemistry():
    """The locked-bond test must reject a long, freely rotating bond."""
    for angle in (60, 180, 300):
        assert chemical_key(_ethanol_conformer(angle))["locked_bond_parity"] == []


def test_charge_is_refused_when_the_potential_ignores_it():
    atoms = molecule("CH3CH2OH")
    with pytest.raises(ValueError):
        chemical_key(atoms, charge_sensitive=False, charge=-1)
    keyed = chemical_key(atoms, charge_sensitive=True, charge=-1, multiplicity=1)
    assert keyed["charge"] == -1
    assert keyed["key"] != chemical_key(atoms)["key"]


def _rotate_group(atoms, centre, axis_atom, group, turns):
    import math

    rotated = atoms.copy()
    axis = atoms.positions[axis_atom] - atoms.positions[centre]
    axis = axis / np.linalg.norm(axis)
    angle = 2 * math.pi * turns
    for atom in group:
        offset = atoms.positions[atom] - atoms.positions[centre]
        along = np.dot(offset, axis) * axis
        across = offset - along
        rotated.positions[atom] = (
            atoms.positions[centre]
            + along
            + across * math.cos(angle)
            + np.cross(axis, across) * math.sin(angle)
        )
    return rotated


def test_methyl_rotamers_are_not_new_conformers():
    """Fixed-index RMSD calls a 120 degree methyl turn a new conformer; it is not."""
    from prrs.chemistry import automorphisms, symmetric_rmsd
    from prrs.state import aligned_rmsd

    atoms = molecule("CH3CH2OH")
    permutations_, info = automorphisms(atoms)
    assert info["enumerated"] and sorted(info["orbit_sizes"]) == [2, 3]
    assert info["count"] == 12  # 3! methyl times 2! methylene

    methyl = [
        i
        for i in range(len(atoms))
        if atoms.numbers[i] == 1
        and np.linalg.norm(atoms.positions[i] - atoms.positions[0]) < 1.3
    ]
    assert len(methyl) == 3
    turned = _rotate_group(atoms, 0, 1, methyl, 1 / 3)

    assert aligned_rmsd(atoms, turned) > 0.5  # the naive measure is fooled
    assert symmetric_rmsd(atoms, turned, permutations_) < 0.02
    assert same_chemistry(atoms, turned)


def test_symmetric_rmsd_still_separates_real_conformers():
    """Symmetry must not be used to merge anti with gauche."""
    from prrs.chemistry import automorphisms, symmetric_rmsd

    anti, gauche = _ethanol_conformer(180), _ethanol_conformer(60)
    permutations_, _ = automorphisms(anti)
    assert symmetric_rmsd(anti, gauche, permutations_) > 0.15
    # gauche+ and gauche- are genuinely distinct under proper rotations only.
    assert symmetric_rmsd(_ethanol_conformer(60), _ethanol_conformer(300), permutations_) > 0.15


def test_automorphism_enumeration_refuses_to_guess_when_too_large():
    from prrs.chemistry import automorphisms

    atoms = molecule("CH3CH2OH")
    permutations_, info = automorphisms(atoms, limit=2)
    assert not info["enumerated"] and info["reason"] == "candidate_limit"
    assert len(permutations_) == 1  # identity only: over-count, never mis-merge


def _stereocentre(flatten=0.0):
    """Four distinct-colour neighbours, interpolated from tetrahedral towards coplanar.

    Four different elements means four different colours, which is exactly the condition
    for a real stereocentre. `flatten` walks the four directions from tetrahedral to a
    coplanar square, so the parity discriminant goes to zero at 1.0.

    The bond lengths are stretched beyond their physical values on purpose. They stay
    inside the C-X cutoffs, so the bond graph is four edges throughout, while the X-X
    distances stay outside every X-X cutoff -- with real lengths the flattest geometries
    bring two substituents into contact and a fifth edge appears, which would make this a
    test of connectivity rather than of parity.
    """
    lengths = np.array([1.25, 1.55, 2.10, 2.30])
    tetrahedral = np.array(
        [
            [0.0, 0.0, 1.0],
            [0.0, 2 * np.sqrt(2) / 3, -1 / 3],
            [-np.sqrt(6) / 3, -np.sqrt(2) / 3, -1 / 3],
            [np.sqrt(6) / 3, -np.sqrt(2) / 3, -1 / 3],
        ]
    )
    coplanar = np.array([[1.0, 0, 0], [0, 1.0, 0], [-1.0, 0, 0], [0, -1.0, 0]])
    directions = (1 - flatten) * tetrahedral + flatten * coplanar
    directions /= np.linalg.norm(directions, axis=1)[:, None]
    return Atoms("CHFClBr", positions=np.vstack([[0, 0, 0], directions * lengths[:, None]]))


def test_the_flattening_fixture_holds_the_bond_graph_still():
    """Otherwise the parity tests below would be measuring connectivity."""
    reference = frozenset(encode(_stereocentre(0.0), 1.2).edges)
    assert len(reference) == 4
    for flatten in (0.3, 0.6, 0.9, 0.99, 1.0):
        assert frozenset(encode(_stereocentre(flatten), 1.2).edges) == reference


def test_tetrahedral_parity_is_scale_free():
    """The discriminant must be a shape, not a size.

    `int(np.sign(det))` read the sign of a number whose scale was the product of three
    bond lengths, so no fixed threshold could serve a C-C skeleton and a C-H one at once.
    Scaling the molecule must not touch the recorded configuration; after normalising by
    those three lengths it cannot.
    """
    small = _stereocentre()
    large = small.copy()
    large.set_positions(large.positions * 3.0)
    key_small = chemical_key(small, bond_scale=1.2)
    key_large = chemical_key(large, bond_scale=3.6)
    assert key_small["tetrahedral_parity"], "this fixture is supposed to have a centre"
    assert key_small["tetrahedral_parity"] == key_large["tetrahedral_parity"]


def test_a_near_planar_centre_is_reported_unresolved_not_signed():
    """A sign decided by float noise mints macrostates out of thermal jiggle.

    Below the band every geometry gets the same stable token instead of a sign that turns
    over with the last digit. Outside the band nothing moves, which is why the frozen
    digests are untouched by this.
    """
    signed = {
        chemical_key(_stereocentre(f), bond_scale=1.2)["tetrahedral_parity"][0][1]
        for f in (0.0, 0.3, 0.6, 0.9)
    }
    assert signed == {-1}, signed

    for flatten in (0.99, 1.0):
        parity = chemical_key(_stereocentre(flatten), bond_scale=1.2)["tetrahedral_parity"]
        assert parity[0][1] == "planar", (flatten, parity)


def test_exactly_planar_is_the_same_token_as_the_rest_of_the_band():
    """np.sign(0.0) is 0, which the old code emitted as a third numeric parity that no
    configuration corresponds to. It has to join the band, not sit beside it."""
    flat = chemical_key(_stereocentre(1.0), bond_scale=1.2)["tetrahedral_parity"]
    near = chemical_key(_stereocentre(0.995), bond_scale=1.2)["tetrahedral_parity"]
    assert flat[0][1] == near[0][1] == "planar"


def test_a_clean_stereocentre_does_not_depend_on_the_tolerance():
    """A tolerance that could move a well-formed centre would rewrite frozen identities."""
    keys = {
        chemical_key(_stereocentre(0.0), bond_scale=1.2, parity_tolerance=tol)["key"]
        for tol in (0.0, 1e-6, 1e-3, 1e-2, 0.05)
    }
    assert len(keys) == 1
