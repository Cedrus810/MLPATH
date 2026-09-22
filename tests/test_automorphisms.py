"""Acceptance for the automorphism enumeration kernel (user ruling, 2026-09-09).

The intervention was scoped to the kernel inside `automorphisms()`: point-wise
backtracking with an incremental edge/non-edge check, in place of the cartesian product
over colour-class permutations. `canonical_labels`, the identity definition, the event
key, the graph hash and the explicit permutation return value are all untouched, and the
group is still returned in full rather than as generators.

Six requirements were set with the ruling; each has a test below, in that order. The
oracle digests are the OLD implementation's answers, taken before it was replaced -- for
b05 at limit 10,000,000, which is what it needed to finish at all (19.8 s).
"""

import hashlib
import json
import sys
from itertools import permutations as all_permutations
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms
from ase.build import molecule
from ase.io import read

from prrs.chemistry import automorphisms, canonical_labels
from prrs.state import encode

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

# The old implementation's answers, recorded before the replacement. b05's row is the one
# it could only produce at limit 10,000,000; the others it produced at the default.
ORACLE = {
    "b05_reactants": (1296, "6fc11bfeff2a03d9"),
    "ethanol": (12, "bd5aa54d2fc0d02c"),
    "b04_input": (6, "8f03e36ace894433"),
}


def fingerprint(permutations_):
    """Set identity of a permutation group, order-independent."""
    rows = sorted(tuple(int(x) for x in p) for p in permutations_)
    assert len(rows) == len(set(rows)), "the same permutation was returned twice"
    blob = json.dumps(rows, separators=(",", ":")).encode()
    return len(rows), hashlib.sha256(blob).hexdigest()[:16]


# ---------------------------------------------------------------- requirement 1
def test_the_tertiary_substrate_enumerates_inside_the_default_budget():
    """b05's reactants: 1296 permutations at limit 20000, where the old kernel needed 10^7.

    Colour classes 3 and 9 gave 3! x 9! = 2,177,280 candidates before -- exactly
    computable, and exactly the wrong scaling -- for a group of order 1296.
    """
    from benchmarks.b05_tertiary_halide.build import build

    found, info = automorphisms(build(), 1.2, None, 20000)
    assert info["enumerated"], info
    assert len(found) == 1296
    assert info["count"] == info["complete_candidates"] == 1296
    assert info["search_nodes"] < 20000, info
    assert sorted(info["colour_class_sizes"], reverse=True) == [9, 3]


def test_the_neighbouring_grids_are_unchanged():
    """Products stay 144 and 288, b04 stays 6: this fix may not move anything else."""
    b04 = read(ROOT / "benchmarks" / "b04_sn2_chloride" / "input.extxyz")
    found, info = automorphisms(b04, 1.2, None, 20000)
    assert (len(found), info["enumerated"]) == (6, True)

    # From build.py, not from preflight.py: importing a gate script to get a structure
    # runs every gate in it, and its sys.exit -- which is how the first draft of this
    # test failed.
    from benchmarks.b05_tertiary_halide.build import (
        elimination_products,
        substitution_products,
    )

    for expected, structure in (
        (144, substitution_products()),
        (288, elimination_products()),
    ):
        found, info = automorphisms(structure, 1.2, None, 20000)
        assert info["enumerated"] and len(found) == expected, (expected, info)


# ---------------------------------------------------------------- requirement 2
def test_the_full_permutation_set_matches_the_old_implementation():
    """Set equality against the pre-replacement oracle, not just the group order.

    A kernel that dropped or duplicated members while keeping the count would pass an
    order-only check, and the group is consumed by minimising over it -- so what has to
    agree is the set.
    """
    from benchmarks.b05_tertiary_halide.build import build

    cases = {
        "b05_reactants": build(),
        "ethanol": molecule("CH3CH2OH"),
        "b04_input": read(ROOT / "benchmarks" / "b04_sn2_chloride" / "input.extxyz"),
    }
    for name, atoms in cases.items():
        found, info = automorphisms(atoms, 1.2, None, 20000)
        assert info["enumerated"], (name, info)
        assert fingerprint(found) == ORACLE[name], name


# ---------------------------------------------------------------- requirement 3
@pytest.mark.parametrize(
    "name", ["CH3CH2OH", "CH3OH", "C2H6", "CH3CHO", "CH3Cl", "H2O", "NH3", "CH4"]
)
def test_historical_structures_keep_their_group_and_their_consumers(name):
    """Same set, and the two consumers that read it give the same answers.

    The consumers matter more than the set: symmetric_rmsd minimises over the group and
    the event key canonicalises under it, so an unnoticed change here would move identity
    without moving any test that looks at the group directly.
    """
    from prrs.chemistry import symmetric_rmsd

    atoms = molecule(name)
    found, info = automorphisms(atoms, 1.2, None, 20000)
    assert info["enumerated"], info
    # every returned map is a bijection of the atom set and preserves the edge set
    graph = encode(atoms, 1.2)
    for permutation in found:
        assert sorted(int(x) for x in permutation) == sorted(range(len(atoms)))
        relabelled = {
            (
                min(int(permutation[a]), int(permutation[b])),
                max(int(permutation[a]), int(permutation[b])),
            )
            for a, b in graph.edges
        }
        assert relabelled == graph.edges
    # identity first, and RMSD against itself is zero under the group
    assert np.array_equal(found[0], np.arange(len(atoms)))
    assert symmetric_rmsd(atoms, atoms.copy(), found, None) == pytest.approx(0.0, abs=1e-12)


# ---------------------------------------------------------------- requirement 4
@pytest.mark.parametrize("size", [3, 4, 5])
def test_a_brute_force_oracle_on_small_graphs(size):
    """Exhaustive oracle: no misses, no duplicates, every map a bijection preserving edges.

    Brute force over ALL permutations of the atoms -- not just within colour classes --
    so the test does not inherit the kernel's own partition. A colour class that were too
    coarse would show up here as a missing automorphism.
    """
    rng = np.random.default_rng(11)
    for trial in range(6):
        positions = rng.uniform(0.0, 2.2, size=(size, 3))
        symbols = "H" * size if trial % 2 else ("C" + "H" * (size - 1))
        atoms = Atoms(symbols, positions=positions)
        graph = encode(atoms, 1.2)
        labels = canonical_labels(atoms.numbers, graph.edges, np.arange(size))

        expected = set()
        for candidate in all_permutations(range(size)):
            if any(labels[i] != labels[candidate[i]] for i in range(size)):
                continue  # a relabelling across colours is not an automorphism here
            mapped = {
                (min(candidate[a], candidate[b]), max(candidate[a], candidate[b]))
                for a, b in graph.edges
            }
            if mapped == graph.edges:
                expected.add(candidate)

        found, info = automorphisms(atoms, 1.2, None, 20000)
        assert info["enumerated"], info
        actual = {tuple(int(x) for x in p) for p in found}
        assert actual == expected, (trial, sorted(expected - actual), sorted(actual - expected))
        assert len(found) == len(actual), "duplicates returned"


# ---------------------------------------------------------------- requirement 5
def test_a_tiny_limit_still_fails_closed():
    """Overflow discards everything and returns the identity alone, with the old reason.

    `BUDGET_REASONS` keys on the string `candidate_limit`, and the event path turns it
    into IdentityUnavailable -- so the string is a contract, not a message.
    """
    from prrs.chemistry import BUDGET_REASONS

    atoms = molecule("CH3CH2OH")
    found, info = automorphisms(atoms, 1.2, None, limit=2)
    assert info["enumerated"] is False
    assert info["reason"] == "candidate_limit"
    assert info["reason"] in BUDGET_REASONS
    assert len(found) == 1 and np.array_equal(found[0], np.arange(len(atoms)))
    assert info["budget_unit"] == "search_nodes"


def test_a_hard_to_prune_symmetric_graph_still_fails_closed():
    """A big single colour class with no structure to prune on: refuse, do not grind.

    Twelve mutually non-bonded identical atoms have 12! = 479,001,600 automorphisms and
    nothing for the incremental check to cut on, because there are no edges to disagree
    about. The kernel must hit the budget and refuse rather than enumerate for a week --
    the pruning fix is not a claim that every group is now reachable.
    """
    positions = [[6.0 * i, 0.0, 0.0] for i in range(12)]
    atoms = Atoms("H" * 12, positions=positions)
    found, info = automorphisms(atoms, 1.2, None, 20000)
    assert info["enumerated"] is False, info
    assert info["reason"] == "candidate_limit"
    assert len(found) == 1
    assert info["search_nodes"] > 20000 - 1


# ---------------------------------------------------------------- requirement 6
def test_the_budget_reports_its_unit_and_stops_calling_classes_orbits():
    """`orbit_sizes` is gone: it named a refutable partition as a group-theoretic one.

    A colour class is what `canonical_labels` could not tell apart; the orbits of the
    automorphism group are only known once the group is enumerated. Both the enumerated
    and the refused shape carry the budget fields, because a refusal is the case where
    knowing how much was spent matters most.
    """
    atoms = molecule("CH3CH2OH")
    for limit, enumerated in ((20000, True), (2, False)):
        _, info = automorphisms(atoms, 1.2, None, limit)
        assert info["enumerated"] is enumerated
        assert "orbit_sizes" not in info
        assert info["colour_class_sizes"] == sorted(info["colour_class_sizes"], reverse=True)
        assert info["budget_unit"] == "search_nodes"
        assert info["limit"] == limit
        assert isinstance(info["search_nodes"], int)
        if enumerated:
            assert info["complete_candidates"] == info["count"]


# ------------------------------------- P1-4 (user decision 7, 2026-09-20)
def test_the_node_path_refuses_to_build_on_an_uncertified_group():
    """admit() no longer swallows an overflowed group as [identity] (P1-4).

    On overflow `automorphisms()` returns the identity alone with the reason in info,
    and the event path already refuses to canonicalise against that ("An incomplete
    group is never used"). The registry's new-node path used to build the node anyway,
    so conformer matching silently measured RMSD against a one-element "group" and a
    methyl rotation registered as a new conformer forever. The ruling: refuse the node,
    keep the run alive, leave a named reason the run's records can carry.
    """
    from ase.calculators.singlepoint import SinglePointCalculator

    from prrs.config import SearchConfig
    from prrs.network import Registry

    def ethanol():
        atoms = molecule("CH3CH2OH")
        atoms.calc = SinglePointCalculator(
            atoms, energy=-100.0, forces=np.zeros((len(atoms), 3))
        )
        return atoms

    starved = Registry(SearchConfig(automorphism_limit=2))
    admission = starved.admit(ethanol(), None, "t0")
    assert admission.outcome == "rejected"
    assert admission.reason == "automorphism_budget"
    assert admission.detail["reason"] == "candidate_limit"
    assert starved.nodes == []

    # The refusal is the budget's, not a blanket no: the same structure inside the
    # real budget is admitted as usual.
    fed = Registry(SearchConfig())
    assert fed.admit(ethanol(), None, "t0").outcome == "new_chemical_node"
