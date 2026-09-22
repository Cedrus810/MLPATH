"""One net under every quotient in the codebase: relabel the atoms, nothing may move.

Six places quotient away a redundancy -- the chemical key, the symmetry-aware RMSD, the
free-aligned RMSD the conformer reservoir matches on, the colour dedup that picks one
representative per angle and torsion class, the reaction-event key, and the automorphism
enumeration all six of them lean on. Each one had its own relabelling test, written when
that one was built, testing that one at whatever level of detail it happened to need.
The gap between them is where the union-of-two-subgroups bug in `reaction_event_key`
lived for its whole life: every individual quotient had a test, and the property the
quotients share had none.

So this file tests the property, not the sites: permute the atom numbering, recompute,
and require the answer back bit for bit. It is a fuzzer rather than a fixture because the
failing case was a specific pair of incomparable automorphism groups that nobody thought
to write down -- REVIEW_2026-09-05.md:117 asked for exactly that test on 2026-09-05 and
it was never written, and a property test would have found it without anyone knowing to
look.

What is deliberately NOT asserted here, because it is order-dependent by construction and
an assertion would only be weakened later by someone who assumed it was a bug:

  * WHICH representative the colour dedup picks. `_angle_directions` and
    `_dihedral_directions` keep the first member of each colour class in enumeration
    order, so the representative is a function of the numbering. The QUOTIENT is what has
    to be invariant -- the set of classes and its size -- and that is what is checked.
  * WHICH directions `propose` selects. The per-family shuffle is seeded and runs over the
    enumeration order, so a renumbered structure draws a different sample. That is a
    sampling decision, not an identity.
  * The conformer reservoir's leader clustering. `_match_microstate` returns the first
    stored microstate within `basin_rmsd_A`, so its output depends on arrival order. The
    distance it clusters ON must be invariant, and that is `free_aligned_symmetric_rmsd`,
    which is checked.
"""

import itertools

import numpy as np
import pytest
from ase.build import molecule

from prrs.chemistry import (
    _generated_group,
    automorphisms,
    canonical_labels,
    chemical_key,
    free_aligned_symmetric_rmsd,
    reaction_event_key,
    symmetric_rmsd,
)
from prrs.config import SearchConfig
from prrs.perturbations import (
    _angle_directions,
    _bond_directions,
    _dihedral_directions,
)
from prrs.state import encode

# Small, fast, and between them they carry every feature the quotients key on: a genuine
# automorphism group (the methyls), a rotatable bond, a double bond, a ring-free skeleton
# with two heteroatoms, and one structure with no symmetry at all.
SUBJECTS = ("CH4", "CH3OH", "C2H6", "C2H4", "CH3CHO", "HCOOH", "CH3CH2OH", "CH3OCH3")
SEEDS = (1, 7, 13)


def _relabelled(atoms, seed):
    """A renumbering of the same structure, with the map back to the original indices.

    `shuffled[i]` is `atoms[order[i]]`, so `order` reads new index -> old index and is
    exactly the map needed to state a new-numbering result in the old numbering.
    """
    order = np.random.default_rng(seed).permutation(len(atoms))
    return atoms[order], {int(new): int(old) for new, old in enumerate(order)}


def _subject(name):
    atoms = molecule(name)
    atoms.center(vacuum=4.0)
    return atoms


def _jittered(atoms, seed, scale=0.03):
    """A nearby conformer: far enough that the RMSDs are not all zero, near enough that
    the bond graph is untouched, so the identity quotients still have to agree."""
    moved = atoms.copy()
    moved.positions += np.random.default_rng(seed).normal(scale=scale, size=(len(atoms), 3))
    return moved


@pytest.mark.parametrize("name", SUBJECTS)
@pytest.mark.parametrize("seed", SEEDS)
def test_colours_follow_the_atoms_they_belong_to(name, seed):
    """The colour of an atom is a property of the atom, not of the integer naming it."""
    atoms = _subject(name)
    shuffled, back = _relabelled(atoms, seed)
    plain = canonical_labels(atoms.numbers, encode(atoms).edges, encode(atoms).index)
    moved = canonical_labels(shuffled.numbers, encode(shuffled).edges, encode(shuffled).index)
    assert {new: moved[new] for new in moved} == {new: plain[back[new]] for new in moved}


@pytest.mark.parametrize("name", SUBJECTS)
@pytest.mark.parametrize("seed", SEEDS)
def test_chemical_key_is_blind_to_the_numbering(name, seed):
    atoms = _subject(name)
    shuffled, _ = _relabelled(atoms, seed)
    plain, moved = chemical_key(atoms), chemical_key(shuffled)
    assert plain["key"] == moved["key"], (name, seed)
    # Not only the digest: every component that feeds it, so a collision cannot hide a
    # component that did move.
    for field in ("graph", "fragments", "stereo"):
        if field in plain:
            assert plain[field] == moved[field], field


@pytest.mark.parametrize("name", SUBJECTS)
@pytest.mark.parametrize("seed", SEEDS)
def test_automorphism_group_size_is_blind_to_the_numbering(name, seed):
    atoms = _subject(name)
    shuffled, _ = _relabelled(atoms, seed)
    plain, plain_info = automorphisms(atoms)
    moved, moved_info = automorphisms(shuffled)
    assert plain_info["enumerated"] and moved_info["enumerated"]
    assert len(plain) == len(moved), (name, seed, len(plain), len(moved))
    assert sorted(plain_info["colour_class_sizes"]) == sorted(moved_info["colour_class_sizes"])


@pytest.mark.parametrize("name", SUBJECTS)
@pytest.mark.parametrize("seed", SEEDS)
def test_both_rmsds_are_blind_to_the_numbering(name, seed):
    """The reservoir clusters on these two numbers, so a numbering-dependent distance
    would mint or merge conformers according to the order atoms were written down."""
    atoms = _subject(name)
    other = _jittered(atoms, seed + 1000)
    shuffled, back = _relabelled(atoms, seed)
    other_shuffled = other[[back[i] for i in range(len(other))]]

    plain_perms, _ = automorphisms(atoms)
    moved_perms, _ = automorphisms(shuffled)
    plain = symmetric_rmsd(atoms, other, plain_perms)
    moved = symmetric_rmsd(shuffled, other_shuffled, moved_perms)
    assert plain > 0 and abs(plain - moved) < 1e-9, (name, seed, plain, moved)

    # The same for the free-aligned variant, with one rotatable bond declared free so the
    # extra rotational minimisation is actually exercised rather than skipped.
    bonds = sorted(encode(atoms).edges)
    free = (bonds[0],)
    free_moved = (
        tuple(sorted((back_index(back, bonds[0][0]), back_index(back, bonds[0][1])))),
    )
    plain_free = free_aligned_symmetric_rmsd(atoms, other, plain_perms, free_bonds=free)
    moved_free = free_aligned_symmetric_rmsd(
        shuffled, other_shuffled, moved_perms, free_bonds=free_moved
    )
    assert abs(plain_free - moved_free) < 1e-9, (name, seed, plain_free, moved_free)


def back_index(back, old):
    """The new index of an old one: `back` reads new -> old, so this inverts it."""
    return next(new for new, was in back.items() if was == old)


EVENT_SUBJECTS = ("CH3OH", "CH3CH2OH", "C2H6")


def _one_bond_broken(atoms, pull=1.2):
    """Pull one hydrogen off its neighbour: a single-bond break, chemically ordinary."""
    out = atoms.copy()
    graph = encode(atoms)
    hydrogen = next(i for i in range(len(atoms)) if atoms.numbers[i] == 1)
    partner = next(sum(e) - hydrogen for e in graph.edges if hydrogen in e)
    axis = out.positions[hydrogen] - out.positions[partner]
    out.positions[hydrogen] += pull * axis / np.linalg.norm(axis)
    assert encode(out).edges != graph.edges, "the probe structure has to differ"
    return out


def _renumbered_event(name, seed):
    before = _subject(name)
    after = _one_bond_broken(before)
    shuffled_before, back = _relabelled(before, seed)
    shuffled_after = after[[back[i] for i in range(len(after))]]
    return (
        reaction_event_key(before, after),
        reaction_event_key(shuffled_before, shuffled_after),
    )


@pytest.mark.parametrize("name", EVENT_SUBJECTS)
@pytest.mark.parametrize("seed", SEEDS)
def test_event_key_context_half_is_blind_to_the_numbering(name, seed):
    """The endpoint-graph half of the key is built from colours and does survive."""
    plain, moved = _renumbered_event(name, seed)
    assert plain["endpoint_graphs"] == moved["endpoint_graphs"], (name, seed)


@pytest.mark.parametrize("name", EVENT_SUBJECTS)
@pytest.mark.parametrize("seed", SEEDS)
def test_structural_key_is_blind_to_the_numbering(name, seed):
    """The quotient that was wrong twice, fixed 2026-09-20 (decision 1).

    Both endpoints are renumbered by the same permutation, which is what two builds of
    the same system produce: one structure, two orderings of its atoms.

    This was xfail(strict) while the fix was an open design decision, because the
    obvious fix -- canonical-order rendering inside `key` -- would have changed every
    event_key digest ever recorded, including the ones b04's frozen criteria quote.
    The taken route adds `structural_key`: the event as one edge-coloured
    (kept/broken/formed) superposition graph through `canonical_form`, numbering-free
    by construction. `key` itself REMAINS a within-numbering identity -- one run's
    structures share a numbering, so within-run channel dedupe is correct on it, and
    every recorded digest stands. What may no longer happen is quoting `key` across
    runs: that is `structural_key`'s job now, and this test is the marker coming off
    in the same change that introduced it (a strict xfail that passes is itself a
    failure).
    """
    plain, moved = _renumbered_event(name, seed)
    assert plain["structural_key"] == moved["structural_key"], (name, seed)
    assert plain["structural_form"] == moved["structural_form"], (name, seed)
    # the within-numbering fields are allowed to differ -- pinned so nobody reads
    # them as structural identities by accident
    assert plain["canonical"] != moved["canonical"] or plain["key"] == moved["key"]


def test_event_key_still_collapses_equivalent_atoms_within_one_numbering():
    """The job the key was built for, pinned so a numbering fix cannot quietly cost it.

    Breaking any one of ethane's six equivalent C-H bonds is one channel. That is what
    the automorphism quotient buys and it holds today; any future canonical-ordering
    change has to keep it.
    """
    atoms = _subject("C2H6")
    graph = encode(atoms)
    keys = set()
    for hydrogen in [i for i in range(len(atoms)) if atoms.numbers[i] == 1]:
        after = atoms.copy()
        partner = next(sum(e) - hydrogen for e in graph.edges if hydrogen in e)
        axis = after.positions[hydrogen] - after.positions[partner]
        after.positions[hydrogen] += 1.2 * axis / np.linalg.norm(axis)
        keys.add(reaction_event_key(atoms, after)["key"])
    assert len(keys) == 1, keys


@pytest.mark.parametrize("name", SUBJECTS)
@pytest.mark.parametrize("seed", SEEDS)
def test_direction_candidates_are_blind_to_the_numbering(name, seed):
    """Three enumerations, two different invariants, and the difference is the point.

    `_bond_directions` enumerates every pair, so the candidate SET must come back
    identical once its indices are read in the original numbering. The angle and torsion
    enumerations keep one representative per colour class, so the representative is a
    function of the numbering and only the quotient -- how many classes, and which colour
    signatures -- is invariant. Asserting the strict form on those two would be asserting
    a property the code does not have and does not want.
    """
    config = SearchConfig()
    atoms = _subject(name)
    shuffled, back = _relabelled(atoms, seed)
    families = tuple(config.families)

    def bonds(structure):
        graph = encode(structure, config.bond_scale, active=config.active_atoms)
        return _bond_directions(structure, graph, config, families)

    def rename(entry, mapping):
        family, indices, sign = entry[0], entry[1], entry[2]
        return (family, tuple(sorted(mapping(i) for i in indices)), sign)

    plain = {rename(e, lambda i: i) for e in bonds(atoms)}
    moved = {rename(e, lambda i: back[i]) for e in bonds(shuffled)}
    assert plain and plain == moved, (name, seed, plain ^ moved)

    # the two enumerators take different inputs now: angles dedupe on orbits alone,
    # dihedrals still want the colours for representative picking
    for enumerate_directions in (
        lambda graph, labels, orbits, families: _angle_directions(graph, orbits, families),
        _dihedral_directions,
    ):

        def classes(structure):
            from prrs.canonical import canonical_form

            graph = encode(structure, config.bond_scale, active=config.active_atoms)
            form = canonical_form(structure.numbers, graph.edges, graph.index)
            labels = form.colours
            orbits = form.orbits if form.info.get("enumerated") else labels
            return sorted(
                (e[0], tuple(sorted(labels[i] for i in e[1])), e[2])
                for e in enumerate_directions(graph, labels, orbits, families)
            )

        assert classes(atoms) == classes(shuffled), (name, seed, enumerate_directions)


def test_generated_group_is_invariant_under_conjugation():
    """The quotient the event key canonicalises on, fuzzed at the group level.

    Renumbering an event conjugates both endpoints' automorphism groups by the same
    permutation, so the group they generate must conjugate with them and keep its size.
    The union does not: it is not closed, and on the smallest incomparable pair it has 7
    elements where the group has 24. This is the property test that would have caught the
    original bug without anyone having to construct that pair on purpose.
    """
    rng = np.random.default_rng(20260913)
    n = 5
    index = list(range(n))
    checked = 0
    for _ in range(300):
        first = _random_group(rng, n)
        second = _random_group(rng, n)
        conjugator = tuple(rng.permutation(n))
        generated = _generated_group(index, first + second, 20000)
        conjugated = _generated_group(
            index,
            [_conjugate(g, conjugator) for g in first + second],
            20000,
        )
        assert len(generated) == len(conjugated)
        assert {_conjugate(g, conjugator) for g in generated} == set(conjugated)
        if len(set(first + second)) != len(generated):
            checked += 1
    assert checked > 0, "the fuzz never produced a case where the union is not a group"


def _conjugate(permutation, by):
    """by . permutation . by^-1, as a tuple over the same index range."""
    n = len(permutation)
    inverse = [0] * n
    for position, image in enumerate(by):
        inverse[image] = position
    return tuple(by[permutation[inverse[i]]] for i in range(n))


def _random_group(rng, n):
    """A subgroup of S_n, as the closure of one or two random permutations."""
    generators = [tuple(rng.permutation(n)) for _ in range(int(rng.integers(1, 3)))]
    group = {tuple(range(n))}
    frontier = list(group)
    while frontier:
        current = frontier.pop()
        for generator in generators:
            composed = tuple(generator[atom] for atom in current)
            if composed not in group:
                group.add(composed)
                frontier.append(composed)
    return sorted(group)


def test_the_fuzzer_can_fail():
    """A net nobody has seen catch anything is not known to be a net.

    Reintroduces the union-of-subgroups quotient on the pair that separates it and
    requires the invariance assertion to reject it, so a later refactor that quietly makes
    every assertion vacuous is visible here rather than three months downstream.
    """
    n, index = 4, list(range(4))
    star, path = [(0, 1), (0, 2), (0, 3)], [(1, 0), (0, 2), (2, 3)]

    def automorphisms_of(edges):
        target = {frozenset(e) for e in edges}
        return [
            p
            for p in itertools.permutations(range(n))
            if {frozenset((p[a], p[b])) for a, b in edges} == target
        ]

    generators = automorphisms_of(star) + automorphisms_of(path)
    union = sorted(set(generators))
    group = _generated_group(index, generators, 20000)
    assert len(union) == 7 and len(group) == 24

    conjugator = (2, 3, 0, 1)
    conjugated_union = sorted({_conjugate(g, conjugator) for g in union})
    assert conjugated_union != union, "the union must NOT survive conjugation"
    assert sorted({_conjugate(g, conjugator) for g in group}) == group


def test_structural_key_is_direction_free_and_discriminating():
    """The two properties the structural identity owes beyond numbering-freedom.

    Direction: reversing source and endpoint swaps broken and formed, which the
    smaller-certificate rule must absorb -- forward and reverse observations are two
    observations of one channel. Discrimination: two different reactions must not
    collide just because their element patterns match.
    """
    atoms = _subject("C2H6")
    graph = encode(atoms)
    hydrogens = [i for i in range(len(atoms)) if atoms.numbers[i] == 1]
    breaks = []
    for hydrogen in hydrogens[:2]:
        after = atoms.copy()
        partner = next(sum(e) - hydrogen for e in graph.edges if hydrogen in e)
        axis = after.positions[hydrogen] - after.positions[partner]
        after.positions[hydrogen] += 1.2 * axis / np.linalg.norm(axis)
        breaks.append(after)
    first = reaction_event_key(atoms, breaks[0])
    second = reaction_event_key(atoms, breaks[1])
    # equivalent hydrogens: same channel, and the structural key agrees
    assert first["structural_key"] == second["structural_key"]
    # direction-free: reversed arguments give the same structural identity
    reversed_event = reaction_event_key(breaks[0], atoms)
    assert reversed_event["structural_key"] == first["structural_key"]

    # Each orientation renders under whichever direction won the certificate min, so
    # the two forms are equal as an UNORDERED pair of sides -- asserting a per-side
    # correspondence would assume both orientations picked opposite winners, which
    # the min is free to contradict.
    def sides(form):
        return {repr(form["broken"]), repr(form["formed"])}

    assert sides(reversed_event["structural_form"]) == sides(first["structural_form"])
