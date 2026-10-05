"""Design section 31 invariants: permutation invariance and scale invariance."""

from __future__ import annotations

import random

import numpy as np

from pars.config import ParsConfig
from pars.frontier.archive import ParetoArchive

from conftest import TRIO_EARLY, TRIO_LATE, make_candidate


def _population(channel: str = "ch", n: int = 6) -> list:
    candidates = []
    for index in range(n):
        energies = np.asarray([-1.0 + 0.3 * index, 0.5 * index, 0.1 * index, -0.1 * index])
        signatures = (
            (TRIO_EARLY, TRIO_EARLY, TRIO_LATE, TRIO_LATE)
            if index % 2 == 0
            else (TRIO_EARLY, TRIO_LATE, TRIO_LATE, TRIO_LATE)
        )
        candidates.append(
            make_candidate(
                f"pc-{index:03d}", channel, energies, signatures, amplitude=0.1 * index
            )
        )
    return candidates


def _selection(candidates) -> tuple:
    archive = ParetoArchive("ch", candidates, ParsConfig())
    return (
        archive.frontier(),
        tuple(archive.representatives),
        {cid: t.value for cid, t in archive.tiers.items()},
        {cid: m for cid, m in archive.margins.items()},
    )


def test_frontier_is_invariant_under_candidate_order():
    candidates = _population()
    baseline = _selection(candidates)
    shuffled = list(candidates)
    random.Random(7).shuffle(shuffled)
    assert _selection(shuffled) == baseline


def test_scale_invariance_after_normalization():
    # multiplying every energy by 1000 (a unit change) must not alter selection
    candidates = _population()
    baseline = _selection(candidates)
    rescaled = []
    for candidate in candidates:
        scaled = make_candidate(
            candidate.candidate_id,
            candidate.channel_id,
            candidate.energies * 1000.0,
            [sorted(sig) for sig in candidate.signatures],
            amplitude=candidate.provenance.perturbation_amplitude,
        )
        rescaled.append(scaled)
    assert _selection(rescaled) == baseline


def test_channel_conditioned_archive_rejects_foreign_candidates():
    from pars.paths.candidate import build_candidates  # noqa: F401  (import sanity)

    candidates = _population()
    foreign = make_candidate("pc-foreign", "other", [0.0, 0.0, 0.0, 0.0], [TRIO_EARLY] * 4)
    try:
        ParetoArchive("ch", [*candidates, foreign], ParsConfig())
    except ValueError as error:
        assert "does not belong" in str(error)
    else:
        raise AssertionError("archive accepted a candidate from another channel")


def test_empty_archive_rejected():
    import pytest

    with pytest.raises(ValueError, match="empty archive"):
        ParetoArchive("ch", [], ParsConfig())
