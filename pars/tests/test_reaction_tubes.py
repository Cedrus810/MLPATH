"""Reaction tubes (design sections 5, 15-16).

Same endpoints, same objectives, but different bond-change ordering = different
mechanism = different tube. The tube concept exists precisely so that
objective-space compression cannot merge physically distinct mechanisms.
"""

from __future__ import annotations

import pytest

from pars.config import ParsConfig
from pars.frontier.archive import ParetoArchive
from pars.paths.cluster import build_tubes
from pars.paths.fingerprint import frechet_path_distance

from conftest import TRIO_EARLY, TRIO_LATE, make_candidate

PROTON_FIRST = (TRIO_EARLY, TRIO_LATE, TRIO_EARLY, TRIO_LATE)  # different ordering
CLEAVAGE_FIRST = (TRIO_EARLY, TRIO_EARLY, TRIO_LATE, TRIO_LATE)


def test_different_mechanisms_form_different_tubes():
    candidates = [
        make_candidate("pc-mech-1", "ch", [0.0, 0.5, 0.1, 0.0], PROTON_FIRST),
        make_candidate("pc-mech-2", "ch", [0.0, 0.5, 0.1, 0.0], CLEAVAGE_FIRST),
    ]
    tubes = build_tubes(candidates, "ch", threshold=0.3)
    assert len(tubes) == 2
    assert {tube.member_ids[0] for tube in tubes} == {"pc-mech-1", "pc-mech-2"}


def test_near_identical_paths_share_a_tube():
    candidates = [
        make_candidate("pc-a", "ch", [0.0, 0.5, 0.1, 0.0], CLEAVAGE_FIRST),
        make_candidate("pc-b", "ch", [0.0, 0.5, 0.1, 0.0], CLEAVAGE_FIRST),
        make_candidate("pc-c", "ch", [0.0, 0.5, 0.1, 0.0], PROTON_FIRST),
    ]
    tubes = build_tubes(candidates, "ch", threshold=0.3)
    assert len(tubes) == 2
    members = {frozenset(tube.member_ids) for tube in tubes}
    assert {"pc-a", "pc-b"} in members
    assert {"pc-c"} in members


def test_frechet_distance_between_the_two_mechanisms_is_large():
    assert frechet_path_distance(PROTON_FIRST, CLEAVAGE_FIRST) >= 0.5


def test_missing_signatures_fail_closed():
    candidate = make_candidate("pc-x", "ch", [0.0, 0.1, 0.2], [TRIO_EARLY] * 3)
    candidate.signatures = None
    with pytest.raises(ValueError, match="signatures"):
        build_tubes([candidate], "ch", threshold=0.3)


def test_tube_alternates_are_retained_in_the_archive():
    """Failed refinement of one representative must not destroy the mechanism:
    the tube keeps its alternates (design section 16)."""
    candidates = [
        make_candidate(f"pc-tube-{i}", "ch", [0.0, 0.5 - 0.1 * i, 0.1, 0.0], PROTON_FIRST)
        for i in range(3)
    ]
    # a generous gamma keeps all three near-frontier candidates in the active
    # set, so the single tube holds the whole mechanism with its alternates
    archive = ParetoArchive("ch", candidates, ParsConfig(margin_gamma=2.5))
    assert len(archive.tubes) == 1
    tube = archive.tubes[0]
    assert len(tube.member_ids) == 3
    assert tube.representative_id is not None
    assert tube.representative_id in tube.member_ids
    alternates = set(tube.member_ids) - {tube.representative_id}
    assert alternates, "tube must retain alternates"


def test_threshold_extremes_fail_closed():
    candidate = make_candidate("pc-a", "ch", [0.0, 0.1, 0.2], [TRIO_EARLY] * 3)
    with pytest.raises(ValueError, match="threshold"):
        build_tubes([candidate], "ch", threshold=0.0)
