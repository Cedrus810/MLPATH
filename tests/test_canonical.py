"""Acceptance for the canonical-form engine (PLAN_2026-09-20.md, section 1.1).

Five requirements, one per test below, in the order the plan's acceptance table lists
them. The anchor is the isomorphism fuzz: the certificate must be a function of the
labelled graph alone, not of the labelling it was handed -- every other guarantee here
(group completeness, degeneration, the wrapper contract) is a consequence of that one.

Two references are pinned against, deliberately, the OLD code rather than the new one:
the pre-change `canonical_labels` is inlined verbatim as `_wl_reference` (comparing the
degeneration against the new wrapper alone would be circular -- wrapper and engine now
share one implementation), and the automorphism set for ethanol is compared with the
recorded pre-change kernel digest from tests/test_automorphisms.py.
"""

import hashlib
import json
import sys
from collections import Counter
from itertools import permutations as all_permutations
from pathlib import Path

import numpy as np
import pytest
from ase.data import chemical_symbols
from ase.build import molecule
from ase.io import read

from prrs.canonical import canonical_form
from prrs.chemistry import BUDGET_REASONS, automorphisms, canonical_labels
from prrs.state import encode

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

# The pre-change kernel's answer for ethanol (tests/test_automorphisms.py ORACLE),
# recorded before that kernel was replaced and now carried across the IR rewrite.
ETHANOL_ORACLE = (12, "bd5aa54d2fc0d02c")


def fingerprint(permutations_):
    """Set identity of a permutation group, order-independent (same recipe as the
    automorphisms acceptance file, so the digests are comparable)."""
    rows = sorted(tuple(int(x) for x in p) for p in permutations_)
    assert len(rows) == len(set(rows)), "the same permutation was returned twice"
    blob = json.dumps(rows, separators=(",", ":")).encode()
    return len(rows), hashlib.sha256(blob).hexdigest()[:16]


def _wl_reference(numbers, edges, index):
    """`canonical_labels` as it stood before PLAN section 1.1 (chemistry.py:51-66),
    inlined verbatim. This is the byte-level target the degeneration is held to."""

    def digest(text):
        return hashlib.sha256(text.encode()).hexdigest()[:16]

    adjacency = {int(i): set() for i in index}
    for a, b in edges:
        adjacency[a].add(b)
        adjacency[b].add(a)
    labels = {int(i): chemical_symbols[int(numbers[i])] for i in index}
    for _ in range(len(index)):
        refined = {
            i: digest(labels[i] + "|" + ",".join(sorted(labels[j] for j in adjacency[i])))
            for i in labels
        }
        if Counter(refined.values()) == Counter(labels.values()):
            break
        labels = refined
    return labels


def _graph(atoms):
    """(numbers, edges, index) straight out of the project's own encoder."""
    graph = encode(atoms, 1.2)
    return atoms.numbers, graph.edges, graph.index


def _brute_force_automorphisms(numbers, edges, index, edge_colours=None):
    """Every permutation of `index` that preserves numbers, edges, and (if given) edge
    colours. Exhaustive over ALL permutations, so it shares nothing with the engine's
    colour-restricted search and can catch a partition that is too coarse."""
    atoms_in = [int(i) for i in index]
    normalized = {(min(a, b), max(a, b)) for a, b in edges}
    expected = set()
    for candidate in all_permutations(range(len(atoms_in))):
        image = {atom: atoms_in[candidate[k]] for k, atom in enumerate(atoms_in)}
        if any(numbers[image[a]] != numbers[a] for a in atoms_in):
            continue
        moved = {(min(image[a], image[b]), max(image[a], image[b])) for a, b in edges}
        if moved != normalized:
            continue
        if edge_colours is not None and any(
            edge_colours[(min(image[a], image[b]), max(image[a], image[b]))]
            != edge_colours[(min(a, b), max(a, b))]
            for a, b in edges
        ):
            continue
        expected.add(tuple(image[a] for a in atoms_in))
    return expected


def _assert_members_are_automorphisms(form, numbers, edges, index, edge_colours=None):
    """Belt and braces around the fuzz: every reported group member really is a
    bijection of the atom set preserving elements and (edge-coloured) edges. A group
    that grew by two leaves agreeing on rendered edges but not on chemistry would be
    caught here rather than trusted."""
    atoms_in = [int(i) for i in index]
    normalized = {(min(a, b), max(a, b)) for a, b in edges}
    for permutation in form.group:
        assert sorted(int(x) for x in permutation) == sorted(atoms_in)
        image = {atom: int(permutation[k]) for k, atom in enumerate(atoms_in)}
        assert all(numbers[image[a]] == numbers[a] for a in atoms_in)
        moved = {(min(image[a], image[b]), max(image[a], image[b])) for a, b in edges}
        assert moved == normalized
        if edge_colours is not None:
            assert all(
                edge_colours[(min(image[a], image[b]), max(image[a], image[b]))]
                == edge_colours[(min(a, b), max(a, b))]
                for a, b in edges
            )


# ---------------------------------------------------------------- the anchor
def test_certificate_is_byte_identical_under_one_thousand_random_relabelings():
    """Random graphs, random relabelings: same certificate, byte for byte.

    This is the plan's core assertion ("整件事的锚"). A canonical form that moved under
    relabelling would make every downstream identity -- chemical key, event key, conformer
    match -- a function of the input atom order, which is exactly the defect section 1
    exists to remove. Group size and membership validity are asserted alongside, since a
    certificate that is stable but bogus would otherwise pass this test alone.
    """
    rng = np.random.default_rng(20260920)
    for _ in range(1000):
        n = int(rng.integers(4, 10))
        numbers = np.array([int(rng.choice([1, 6, 8])) for _ in range(n)])
        edges = {(u, v) for u in range(n) for v in range(u + 1, n) if rng.random() < 0.4}
        first = canonical_form(numbers, edges, np.arange(n))
        assert first.info["enumerated"], first.info
        _assert_members_are_automorphisms(first, numbers, edges, np.arange(n))

        permutation = rng.permutation(n)
        image = {u: int(permutation[u]) for u in range(n)}
        relabelled_numbers = np.empty(n, dtype=int)
        for u in range(n):
            relabelled_numbers[image[u]] = numbers[u]
        relabelled_edges = {
            (min(image[u], image[v]), max(image[u], image[v])) for u, v in edges
        }
        second = canonical_form(relabelled_numbers, relabelled_edges, np.arange(n))
        assert second.info["enumerated"], second.info

        assert second.certificate == first.certificate
        assert len(second.group) == len(first.group)
        assert sorted(second.colours.values()) == sorted(first.colours.values())
        _assert_members_are_automorphisms(
            second, relabelled_numbers, relabelled_edges, np.arange(n)
        )


# ---------------------------------------------------------------- degeneration
def test_without_edge_colours_the_engine_degenerates_into_the_old_canonical_labels():
    """No edge colours => colours byte-identical to the pre-change canonical_labels.

    Run on three REAL graphs (ethanol, the b04 input) and two synthetic ones (a fully
    symmetric carbon ring, and twelve isolated hydrogens, which exercise the immediate
    stop: one pass, no split, the INITIAL labels returned untouched). The comparison is
    against the old code inlined at the top of this file, not against the new wrapper.
    """
    ring_numbers = np.array([6] * 6)
    ring_edges = {(u, (u + 1) % 6) for u in range(6)}
    star_numbers = np.array([1] * 12)
    cases = {
        "ethanol": _graph(molecule("CH3CH2OH")),
        "b04_input": _graph(read(ROOT / "benchmarks" / "b04_sn2_chloride" / "input.extxyz")),
        "carbon_ring": (ring_numbers, ring_edges, np.arange(6)),
        "twelve_isolated_hydrogens": (star_numbers, set(), np.arange(12)),
    }
    for name, (numbers, edges, index) in cases.items():
        reference = _wl_reference(numbers, edges, index)
        assert canonical_labels(numbers, edges, index) == reference, name
        assert canonical_form(numbers, edges, index).colours == reference, name


# ---------------------------------------------------------------- the group
@pytest.mark.parametrize("name", ["NH3", "CH3Cl", "H2O"])
def test_the_group_matches_a_brute_force_oracle(name):
    """|Aut| and every member, against exhaustive search: NH3 and CH3Cl have 6, water 2.

    The comparison is per-permutation in the old convention -- position-indexed arrays
    carrying global atom values -- not just a count, so a kernel that swapped the
    convention would fail here even with the right group order.
    """
    atoms = molecule(name)
    numbers, edges, index = _graph(atoms)
    expected = _brute_force_automorphisms(numbers, edges, index)
    found, info = automorphisms(atoms, 1.2, None, 20000)
    assert info["enumerated"], info
    assert len(expected) == {"NH3": 6, "CH3Cl": 6, "H2O": 2}[name]
    assert {tuple(int(x) for x in p) for p in found} == expected


def test_ethanol_keeps_the_recorded_group_and_a_tiny_limit_fails_closed():
    """Ethanol: the exact 12-permutation set the pre-change kernel returned, and the
    overflow contract: [identity] + reason, the string BUDGET_REASONS keys on."""
    atoms = molecule("CH3CH2OH")
    found, info = automorphisms(atoms, 1.2, None, 20000)
    assert info["enumerated"], info
    assert fingerprint(found) == ETHANOL_ORACLE

    starved, starved_info = automorphisms(atoms, 1.2, None, limit=2)
    assert starved_info["enumerated"] is False
    assert starved_info["reason"] == "candidate_limit"
    assert starved_info["reason"] in BUDGET_REASONS
    assert len(starved) == 1 and np.array_equal(starved[0], np.arange(len(atoms)))

    numbers, edges, index = _graph(atoms)
    form = canonical_form(numbers, edges, index, limit=2)
    assert form.info["reason"] == "candidate_limit"
    assert form.order == () and form.certificate == ()
    assert len(form.group) == 1 and np.array_equal(form.group[0], index)
    # the orbits reported alongside an uncertified group are the identity's own
    assert set(form.orbits.values()) == set(range(len(atoms)))


# ---------------------------------------------------------------- edge colours
def test_edge_colours_split_a_graph_the_vertex_colours_cannot():
    """One marked bond in a bare carbon ring: 1-WL without edge colours leaves all six
    carbons in one class forever (they are one orbit, |Aut| = 12); with the edge colour
    the refinement separates the marked bond's ends -- and the group collapses to the
    two maps that preserve the marking. Exhaustively checked, per permutation."""

    class Ring:
        def __init__(self):
            self.numbers = np.array([6] * 6)
            self.edges = {(u, (u + 1) % 6) for u in range(6)}
            self.index = np.arange(6)
            self.colours = {
                (min(a, b), max(a, b)): "double" if (a, b) == (0, 1) else "single"
                for a, b in self.edges
            }

    ring = Ring()
    plain = canonical_form(ring.numbers, ring.edges, ring.index)
    assert len(set(plain.colours.values())) == 1
    assert len(plain.group) == 12  # the full dihedral group, colours can not see the mark

    coloured = canonical_form(ring.numbers, ring.edges, ring.index, edge_colours=ring.colours)
    # stable classes: the marked bond's ends {0, 1}, then {2, 5} and {3, 4} by distance
    # from the mark -- measured, not assumed, which is why the assertion is on the
    # partition itself
    assert len(set(coloured.colours.values())) == 3
    assert coloured.colours[0] == coloured.colours[1] != coloured.colours[2]
    assert coloured.colours[2] == coloured.colours[5] != coloured.colours[3]
    assert coloured.info["colour_class_sizes"] == [2, 2, 2]
    expected = _brute_force_automorphisms(
        ring.numbers, ring.edges, ring.index, edge_colours=ring.colours
    )
    assert {tuple(int(x) for x in p) for p in coloured.group} == expected
    assert len(coloured.group) == 2  # identity + the reflection through the marked bond

    # both accepted spellings of the same colouring agree on the certificate
    via_callable = canonical_form(
        ring.numbers,
        ring.edges,
        ring.index,
        edge_colours=lambda a, b: ring.colours[(min(a, b), max(a, b))],
    )
    assert via_callable.certificate == coloured.certificate


def test_certificates_stay_invariant_under_relabeling_with_edge_colours():
    """The anchor's edge-coloured twin: random 3-way bond colourings (the kept/broken/
    formed alphabet the event key will use) survive random relabelings byte for byte."""
    rng = np.random.default_rng(20260921)
    for _ in range(250):
        n = int(rng.integers(4, 9))
        numbers = np.array([int(rng.choice([1, 6, 8])) for _ in range(n)])
        edges = {(u, v) for u in range(n) for v in range(u + 1, n) if rng.random() < 0.45}
        colours = {(u, v): str(rng.choice(["kept", "broken", "formed"])) for u, v in edges}
        first = canonical_form(numbers, edges, np.arange(n), edge_colours=colours)
        assert first.info["enumerated"], first.info
        _assert_members_are_automorphisms(first, numbers, edges, np.arange(n), colours)

        permutation = rng.permutation(n)
        image = {u: int(permutation[u]) for u in range(n)}
        relabelled_numbers = np.empty(n, dtype=int)
        for u in range(n):
            relabelled_numbers[image[u]] = numbers[u]
        relabelled_edges = {
            (min(image[u], image[v]), max(image[u], image[v])) for u, v in edges
        }
        relabelled_colours = {
            (min(image[u], image[v]), max(image[u], image[v])): colours[(u, v)]
            for u, v in edges
        }
        second = canonical_form(
            relabelled_numbers, relabelled_edges, np.arange(n), edge_colours=relabelled_colours
        )
        assert second.info["enumerated"], second.info
        assert second.certificate == first.certificate
        assert len(second.group) == len(first.group)
        _assert_members_are_automorphisms(
            second, relabelled_numbers, relabelled_edges, np.arange(n), relabelled_colours
        )


# ---------------------------------------------------------------- wrapper contract
def test_the_automorphisms_wrapper_keeps_the_old_contract_field_for_field():
    """automorphisms() is canonical_form(...).group now, but its return value and info
    field names are frozen by tests/test_automorphisms.py (17 tests) and consumed at
    network.py:311 and chemistry.symmetric_rmsd: the field SET is asserted exactly, on
    both the enumerated and the refused shape, and the wrapper is checked to agree with
    the engine it wraps."""
    atoms = molecule("CH3CH2OH")
    found, info = automorphisms(atoms, 1.2, None, 20000)
    assert set(info) == {
        "enumerated",
        "count",
        "complete_candidates",
        "colour_class_sizes",
        "search_nodes",
        "budget_unit",
        "limit",
    }
    assert "orbit_sizes" not in info
    assert isinstance(found, list) and all(isinstance(p, np.ndarray) for p in found)
    assert np.array_equal(found[0], np.arange(len(atoms)))
    assert fingerprint(found) == ETHANOL_ORACLE
    assert info["count"] == info["complete_candidates"] == len(found) == 12
    assert info["colour_class_sizes"] == sorted(info["colour_class_sizes"], reverse=True)
    assert info["budget_unit"] == "search_nodes"
    assert info["limit"] == 20000
    assert isinstance(info["search_nodes"], int) and info["search_nodes"] > 0

    starved, starved_info = automorphisms(atoms, 1.2, None, limit=2)
    assert set(starved_info) == {
        "enumerated",
        "reason",
        "colour_class_sizes",
        "search_nodes",
        "budget_unit",
        "limit",
    }
    assert starved_info["search_nodes"] > 2 - 1  # the budget is really what stopped it

    form = canonical_form(*_graph(atoms), limit=20000)
    assert len(form.group) == len(found)
    assert all(np.array_equal(a, b) for a, b in zip(form.group, found))
    assert form.info == info


def test_a_discrete_refinement_costs_one_search_node():
    """Sanity on the cheap path: four atoms, four distinct elements, so the INITIAL
    colouring is already discrete -- the root is the one leaf, the group is the identity
    at one search node. The budget exists for the symmetric cases; an asymmetric
    structure must never burn any of it."""
    numbers = np.array([6, 8, 1, 17])  # C-O-H plus a chlorine on the carbon
    edges = {(0, 1), (1, 2), (0, 3)}
    index = np.arange(4)
    form = canonical_form(numbers, edges, index)
    assert form.info["enumerated"] and form.info["search_nodes"] == 1
    assert len(form.group) == 1 and np.array_equal(form.group[0], index)
    expected = _brute_force_automorphisms(numbers, edges, index)
    assert {tuple(int(x) for x in p) for p in form.group} == expected


def test_the_canonical_labelling_itself_is_invariant_not_just_the_certificate():
    """The order must transform consistently, not merely agree on a certificate.

    The certificate being equal is necessary but weaker than what consumers need:
    `structural_form` renders bonds UNDER the canonical order (PLAN section 1.2), so
    the labelling itself has to be a function of the graph -- under a relabelling p,
    the canonical order must come back as exactly (p[old_order[k]]). A form whose
    order depends on which equal-certificate leaf exploration happened to visit
    first passes a certificate-only fuzz and quietly ships a labelling that is not
    a graph property.
    """
    import random

    from ase.build import molecule

    from prrs.state import encode

    rng = random.Random(20260920)
    cases = []
    for name in ("CH3CH2OH", "H2O", "CH4"):
        atoms = molecule(name)
        graph = encode(atoms)
        cases.append((tuple(atoms.numbers), sorted(tuple(sorted(e)) for e in graph.edges)))
    for t in range(8):
        n = rng.randint(4, 7)
        numbers = [6 if i < n // 2 else 1 for i in range(n)]
        edges = [(i, j) for i in range(n) for j in range(i + 1, n) if rng.random() < 0.5]
        cases.append((numbers, sorted(edges)))

    for numbers, edges in cases:
        n = len(numbers)
        index = tuple(range(n))
        base = canonical_form(numbers, edges, index)
        for _ in range(50):
            perm = list(range(n))
            rng.shuffle(perm)
            inv = [0] * n
            for a, new in enumerate(perm):
                inv[new] = a
            relabelled = canonical_form(
                [numbers[a] for a in inv],
                sorted(tuple(sorted((perm[a], perm[b]))) for a, b in edges),
                index,
            )
            assert relabelled.certificate == base.certificate
            # The labelling is a graph property UP TO Aut: two minimal-certificate
            # leaves differ by an automorphism, and which one exploration met first
            # must not leak into the answer. Pull the relabelled order back into old
            # atom names; it must be base_order composed with a group member.
            pulled = [inv[y] for y in relabelled.order]
            slot = {a: k for k, a in enumerate(range(n))}
            in_orbit = any(
                all(pulled[k] == int(g[slot[base.order[k]]]) for k in range(n))
                for g in base.group
            )
            assert in_orbit, (
                "the canonical labelling moved by more than an automorphism: "
                "exploration order leaked into the answer"
            )


def test_the_same_graph_twice_gives_byte_identical_output():
    """Reproducibility, the property every consumer actually leans on.

    Isomorphic graphs may label up to Aut -- that is theory, pinned one test up. But
    the SAME graph must return the same order, group, orbits and certificate on every
    call, or nothing downstream (keys, RMSD minima, structural forms) is a function
    of anything. This is the assertion an unseeded choice anywhere in the engine
    cannot survive.
    """
    from ase.build import molecule

    from prrs.state import encode

    atoms = molecule("CH3CH2OH")
    graph = encode(atoms)
    numbers, edges = tuple(atoms.numbers), sorted(tuple(sorted(e)) for e in graph.edges)
    index = tuple(range(len(numbers)))
    first = canonical_form(numbers, edges, index)
    for _ in range(20):
        again = canonical_form(numbers, edges, index)
        assert again.order == first.order
        assert again.certificate == first.certificate
        assert again.orbits == first.orbits
        assert all(np.array_equal(a, b) for a, b in zip(again.group, first.group)) and len(
            again.group
        ) == len(first.group)
