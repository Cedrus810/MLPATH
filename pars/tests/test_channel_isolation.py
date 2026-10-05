"""Channel isolation: a superior path in A->B must not eliminate A->C
(design sections 4 and 31)."""

from __future__ import annotations

import numpy as np

from pars.config import ParsConfig
from pars.frontier.archive import ParetoArchive
from pars.schema import ParetoStatus

from conftest import TRIO_EARLY, TRIO_LATE, make_candidate


def _ab_population() -> list:
    """Channel A->B: many excellent low-barrier paths."""
    candidates = []
    for index in range(5):
        candidates.append(
            make_candidate(
                f"pc-ab-{index}",
                "rc-ab",
                np.asarray([-1.0, 0.2 + 0.05 * index, 0.1, -0.5]),
                [TRIO_EARLY, TRIO_EARLY, TRIO_LATE, TRIO_LATE],
            )
        )
    return candidates


def _ac_population() -> list:
    """Channel A->C: a single much worse path.

    Objectively dominated by every A->B path: lower barrier, lower curvature,
    whatever. Globally it would be pruned; per channel it must survive.
    """
    return [
        make_candidate(
            "pc-ac-0",
            "rc-ac",
            np.asarray([-1.0, 5.0, 4.0, 3.0]),  # huge excursion
            [TRIO_EARLY, TRIO_LATE, TRIO_EARLY, TRIO_LATE],  # and it recrosses
        )
    ]


def test_weak_channel_survives_channel_conditioning():
    ab = ParetoArchive("rc-ab", _ab_population(), ParsConfig())
    ac = ParetoArchive("rc-ac", _ac_population(), ParsConfig())

    assert ac.tiers["pc-ac-0"] is ParetoStatus.ACTIVE_FRONTIER
    assert ac.representatives == ["pc-ac-0"]
    assert len(ab.representatives) >= 1


def test_merged_archive_would_have_pruned_it():
    # the counterfactual that justifies channel conditioning: without the
    # channel split the single A->C path is not on the global frontier
    from pars.frontier.margin import assign_tiers

    merged_vectors = {}
    for candidate in _ab_population() + _ac_population():
        merged_vectors[candidate.candidate_id] = {
            "energy_excursion": candidate.energies.max() - candidate.energies[0],
            "recrossing": 0.0 if candidate.candidate_id.startswith("pc-ab") else 2.0,
        }
    config = ParsConfig()
    tiers, _ = assign_tiers(merged_vectors, config)
    assert tiers["pc-ac-0"] is not ParetoStatus.ACTIVE_FRONTIER


def test_archive_refuses_to_build_across_channels():
    import pytest

    candidates = _ab_population() + _ac_population()
    with pytest.raises(ValueError, match="does not belong"):
        ParetoArchive("rc-ab", candidates, ParsConfig())
