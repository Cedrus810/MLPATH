"""The direction budget is ranked by a spectrum that was already paid for.

Without one, every candidate's rank is 0, the stable sort in `propose` is a no-op, and the
order inside a family is the seeded shuffle -- a lottery. On b04 that lottery drew the one
direction of 28 that reacted first, with nothing putting it there; the probe in
`docs/experiments/spectral_proposal_probe.py` measured that the softest non-torsional
modes put the same direction first by construction, at 0.9651 against 0.6893 for the next
distinct one. These tests pin the mechanism, not that number: a calculator is not needed
to check that a ranking orders what it is given.
"""

import numpy as np
import pytest
from ase.build import molecule

from prrs.calculators import demo_atoms, double_well_factory
from prrs.config import SearchConfig
from prrs.perturbations import propose, spectral_overlap
from prrs.runner import soft_internal_modes
from prrs import internal


def _mass_weighted(atoms, kind, indices):
    gradient = internal.gradient(atoms.positions, kind, indices)
    vector = (gradient / np.sqrt(atoms.get_masses())[:, None]).reshape(-1)
    return vector / np.linalg.norm(vector)


def _order(atoms, config, modes):
    """The identity of every selected direction, in the order propose returns it."""
    report = {}
    list(propose(atoms, config, seed=17, report=report, soft_modes=modes))
    return report["selected_ids"], report["candidate_ids"]


def test_a_direction_that_is_the_soft_mode_is_ranked_first():
    """Hand the ranker one mode that IS a bond's own direction and it must lead."""
    atoms = molecule("CH3CH2OH")
    atoms.center(vacuum=4.0)
    config = SearchConfig(families=("stretch",), max_directions=64)
    target = (0, 1)
    modes = ((0.01, _mass_weighted(atoms, "bond", target)),)

    ranked, _ = _order(atoms, config, modes)
    assert ranked[0] == f"stretch:{target[0]}-{target[1]}:+", ranked[:4]

    # And the overlap that put it there is the self-overlap, so it is exactly 1.
    assert spectral_overlap(atoms, "bond", target, modes) == pytest.approx(1.0, abs=1e-12)


def test_the_ranking_is_what_changed_and_nothing_else():
    """Same structure, same seed: only the order moves, never the candidate set."""
    atoms = molecule("CH3CH2OH")
    atoms.center(vacuum=4.0)
    config = SearchConfig(families=("stretch", "compress", "bend"), max_directions=200)
    modes = ((0.01, _mass_weighted(atoms, "bond", (0, 1))),)

    plain, plain_pool = _order(atoms, config, ())
    ranked, ranked_pool = _order(atoms, config, modes)
    # The CANDIDATE set is a property of the structure and must not move. Which of them
    # get selected is exactly what a ranking is for, so that set is allowed to differ --
    # the quota takes a fraction of each family, ranked or not.
    assert set(plain_pool) == set(ranked_pool), set(plain_pool) ^ set(ranked_pool)
    assert plain != ranked, "a ranking that changes no order is not a ranking"


def test_no_spectrum_leaves_the_existing_order_untouched():
    """A microstate whose Hessian could not be had must fall back, not fail."""
    atoms = molecule("CH3OH")
    atoms.center(vacuum=4.0)
    config = SearchConfig(families=("stretch",), max_directions=64)
    assert _order(atoms, config, ()) == _order(atoms, config, ())


def test_soft_internal_modes_drops_the_rigid_body_block():
    """The bug that scored all 28 of b04's directions at exactly zero.

    `curvature_spectrum` projects the six trivial modes out rather than deleting them, so
    they stay in the spectrum as a null block at the SOFT end. Taking the softest modes
    literally selects translations and rotations, which are orthogonal to every internal
    coordinate by construction.
    """
    config = SearchConfig(minimum_check="hessian", hessian_source="fd")
    atoms = demo_atoms()
    atoms.calc = double_well_factory()
    modes, note = soft_internal_modes(atoms, double_well_factory, config)
    assert modes, note
    assert len(modes) <= config.spectral_rank_modes

    from prrs.runner import _trivial_modes

    trivial = _trivial_modes(atoms.positions, atoms.get_masses())
    for eigenvalue, mode in modes:
        assert float(np.linalg.norm(mode)) == pytest.approx(1.0, abs=1e-9)
        assert float(np.linalg.norm(trivial @ mode)) < 0.5, (eigenvalue, "still rigid-body")


def test_a_spectrum_that_cannot_be_had_returns_a_reason_not_an_exception():
    """A ranking is never worth failing a search for."""

    class Exploding:
        def __call__(self):
            raise RuntimeError("no calculator here")

    atoms = demo_atoms()
    atoms.calc = double_well_factory()
    modes, note = soft_internal_modes(atoms, Exploding(), SearchConfig())
    assert modes == () and "RuntimeError" in note, note
