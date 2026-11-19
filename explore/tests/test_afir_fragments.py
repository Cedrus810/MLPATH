"""SC-AFIR fragment generation (explore.afir_fragments): perturbed structure, A / B groups,
candidate pairs with symmetry. Pure geometry."""

import numpy as np
import pytest
from ase import Atoms
from ase.build import molecule

from explore.afir_fragments import (
    candidate_pairs,
    perturbed_structure,
    ratios,
    sc_fragments,
)
from prrs.config import SearchConfig

CONFIG = SearchConfig(closed_shell_only=False)


def _b04():
    """Cl- (5) at the back of CH3Cl (0 C, 1 Cl, 2-4 H), along the C-Cl axis."""
    m = molecule("CH3Cl")
    c = next(i for i, z in enumerate(m.numbers) if z == 6)
    cl = next(i for i, z in enumerate(m.numbers) if z == 17)
    hs = [i for i, z in enumerate(m.numbers) if z == 1]
    order = [c, cl] + hs
    pos = m.positions[order]
    axis = (pos[0] - pos[1]) / np.linalg.norm(pos[0] - pos[1])
    pos = np.vstack([pos, pos[0] + 3.2 * axis])
    return Atoms("CClHHHCl", positions=pos)


def test_perturbed_structure_shortens_one_pair_and_keeps_the_rest():
    a = _b04()
    _, d0 = ratios(a.numbers, a.positions)
    x = perturbed_structure(a.numbers, a.positions, 0, 5, sign=+1)
    _, d = ratios(a.numbers, x)
    assert d[0, 5] == pytest.approx(0.8 * d0[0, 5], rel=1e-3)
    others = [(i, j) for i in range(6) for j in range(i + 1, 6) if (i, j) != (0, 5)]
    assert max(abs(d[i, j] - d0[i, j]) for i, j in others) < 0.5
    y = perturbed_structure(a.numbers, a.positions, 0, 5, sign=-1)
    assert ratios(a.numbers, y)[1][0, 5] == pytest.approx(1.2 * d0[0, 5], rel=1e-3)


def test_linear_triatomic_bends_under_push():
    """The papers' P-Q-R case: shortening P-R must not just shrink the line."""
    a = Atoms("OCO", positions=[[0, 0, 0], [1.16, 0, 0], [2.32, 0, 0]])
    x = perturbed_structure(a.numbers, a.positions, 0, 2, sign=+1)
    v1, v2 = x[0] - x[1], x[2] - x[1]
    cos = np.dot(v1, v2) / np.linalg.norm(v1) / np.linalg.norm(v2)
    assert cos > -0.99  # no longer linear (cos -1)


def test_sc_fragments_around_the_attacked_carbon():
    a = _b04()
    x = perturbed_structure(a.numbers, a.positions, 0, 5)
    A, B = sc_fragments(a.numbers, x, 0, 5)
    assert 0 in A and B == (5,)
    assert set(A) == {0, 1, 2, 3, 4}  # CH3Cl joins the carbon's side
    A2, B2 = sc_fragments(a.numbers, a.positions, 1, 0)  # bonded pair C-Cl
    assert 1 in A2 and 0 in B2 and not set(A2) & set(B2)


def test_candidate_pairs_respect_cutoff_and_symmetry():
    a = _b04()
    pairs = candidate_pairs(a, CONFIG)
    rr, _ = ratios(a.numbers, a.positions)
    assert all(rr[i, j] <= 3.0 for i, j in pairs)
    # the three H are equivalent: one C-H, one Cl-H, one H-H, one H-Cl(-) representative
    kinds = [tuple(sorted((a.numbers[i], a.numbers[j]))) for i, j in pairs]
    assert kinds.count((1, 6)) == 1 and kinds.count((1, 1)) == 1
    assert (0, 5) in pairs  # the SN2 pair


def test_two_pair_perturbation_is_one_fit_and_reduces_to_one_pair():
    """User 2026-11-15: both distances shortened 20 % in the same least-squares fit; one
    pair gives exactly perturbed_structure."""
    from explore.afir_fragments import perturbed_structure, perturbed_structure_pairs

    rng = np.random.default_rng(2)
    x = rng.normal(0.0, 1.6, (7, 3))
    z = np.full(7, 6)
    one = perturbed_structure(z, x, 0, 4, seed=3)
    assert np.array_equal(one, perturbed_structure_pairs(z, x, [(0, 4)], seed=3))
    two = perturbed_structure_pairs(z, x, [(0, 4), (2, 5)], seed=3)
    d0 = np.linalg.norm(x[:, None] - x[None, :], axis=-1)
    d = np.linalg.norm(two[:, None] - two[None, :], axis=-1)
    for k, m in ((0, 4), (2, 5)):
        assert d[k, m] == pytest.approx(0.8 * d0[k, m], rel=1e-3)
    flipped = perturbed_structure_pairs(z, x, [(2, 5), (0, 4)], seed=3)
    assert np.allclose(two, flipped, atol=1e-6)  # no order dependence


def test_two_contact_combos():
    """Non-bonded contacts (ratio <= 3.0), four different atoms, one per symmetry orbit of
    the pair of pairs; no ring-size filter."""
    from prrs.config import SearchConfig

    from explore.afir_batch import edges_of
    from explore.afir_fragments import ratios, two_contact_combos

    cfg = SearchConfig(closed_shell_only=False)
    # two H2 molecules side by side (a rectangle): bonds 0-1, 2-3; contacts 0-2, 1-3, 0-3, 1-2
    a = Atoms("H4", positions=[[0, 0, 0], [0.74, 0, 0], [0, 1.2, 0], [0.74, 1.2, 0]])
    combos = two_contact_combos(a, cfg)
    bonded = {tuple(e) for e in edges_of(a.numbers, a.positions, cfg.bond_scale)}
    rr, _ = ratios(a.numbers, a.positions)
    for p, q in combos:
        assert len({*p, *q}) == 4
        assert p not in bonded and q not in bonded
        assert rr[p] <= 3.0 and rr[q] <= 3.0
    # symmetry is the bond-graph automorphism group (as candidate_pairs): swapping 0 <-> 1
    # is a graph automorphism, so the parallel {0-2, 1-3} and crossed {0-3, 1-2} pairs of
    # contacts are one orbit even though they differ geometrically
    assert len(combos) == 1
