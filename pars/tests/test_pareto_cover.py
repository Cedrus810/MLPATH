"""Cover determinism and duplicate stability (design section 31)."""

from __future__ import annotations


from pars.config import ParsConfig
from pars.frontier.archive import ParetoArchive
from pars.frontier.cover import objective_space_cover

from conftest import TRIO_EARLY, TRIO_LATE, make_candidate


def _two_tradeoff_candidates() -> list:
    return [
        make_candidate(
            "pc-low-x",
            "ch",
            [0.0, 0.4, 0.2, 0.1],
            [TRIO_EARLY, TRIO_LATE, TRIO_LATE, TRIO_LATE],
        ),
        make_candidate(
            "pc-low-y",
            "ch",
            [0.4, 0.0, 0.2, 0.1],
            [TRIO_EARLY, TRIO_LATE, TRIO_LATE, TRIO_LATE],
        ),
    ]


def test_cover_is_deterministic():
    candidates = _two_tradeoff_candidates()
    first = ParetoArchive("ch", candidates, ParsConfig())
    second = ParetoArchive("ch", _two_tradeoff_candidates(), ParsConfig())
    assert first.representatives == second.representatives
    assert first.frontier() == second.frontier()


def test_duplicating_a_trajectory_100x_does_not_change_coverage():
    # design 31 duplicate stability: the representative set may pick any one
    # of the duplicates, but the number of representatives and the RAW
    # objective content they cover must not change. (Normalized values are
    # allowed to move: median/IQR are fitted per population, and a crowd of
    # duplicates legitimately changes the scale.)
    base = _two_tradeoff_candidates()
    single = ParetoArchive("ch", base, ParsConfig())
    baseline_vectors = {
        cid: dict(single.candidates[cid].objectives_raw) for cid in single.representatives
    }

    duplicated = list(base)
    for index in range(100):
        source = base[index % len(base)]
        duplicate = make_candidate(
            f"pc-dup-{index:03d}",
            source.channel_id,
            source.energies,
            [sorted(sig) for sig in source.signatures],
            amplitude=source.provenance.perturbation_amplitude,
        )
        duplicated.append(duplicate)
    crowd = ParetoArchive("ch", duplicated, ParsConfig())
    crowd_vectors = {
        cid: dict(crowd.candidates[cid].objectives_raw) for cid in crowd.representatives
    }

    assert len(crowd.representatives) == len(single.representatives)
    assert _content(crowd_vectors) == _content(baseline_vectors)


def _content(vectors: dict) -> list:
    def rounded(vec: dict) -> tuple:
        return tuple(sorted((k, None if v is None else round(v, 9)) for k, v in vec.items()))

    return sorted(rounded(vec) for vec in vectors.values())


def test_cover_respects_tolerance():
    vectors = {
        "a": {"x": 0.0, "y": 0.0},
        "b": {"x": 0.05, "y": 0.0},  # within tolerance of a
        "c": {"x": 1.0, "y": 1.0},  # a genuinely different trade-off
    }
    selected = objective_space_cover(vectors, ["a", "b", "c"], tolerance=0.10)
    assert selected[0] == "a"
    assert "b" not in selected
    assert "c" in selected


def test_diversity_preserving_selection_keeps_every_tube_represented():
    # a genuine trade-off (one path lower excursion, the other cleaner) keeps
    # both candidates non-dominated and their distinct mechanisms separate
    candidates = [
        make_candidate("pc-t1-a", "ch", [0.0, 0.5, 0.1, 0.0], [TRIO_EARLY, TRIO_LATE] * 2),
        make_candidate(
            "pc-t2-a",
            "ch",
            [0.0, 0.9, 0.1, 0.0],
            [TRIO_EARLY, TRIO_EARLY, TRIO_LATE, TRIO_LATE],
        ),
    ]
    archive = ParetoArchive("ch", candidates, ParsConfig())
    assert len(archive.tubes) == 2
    assert len(archive.representatives) == 2
    assert all(tube.representative_id for tube in archive.tubes)


def test_representatives_marked_on_candidates():
    candidates = _two_tradeoff_candidates()
    archive = ParetoArchive("ch", candidates, ParsConfig())
    for candidate_id in archive.representatives:
        assert archive.candidates[candidate_id].representative is True
        assert archive.candidates[candidate_id].tube_id is not None


def _twin_mechanisms() -> list:
    # identical objectives, different intermediates: two tubes, one objective point
    first = frozenset({(0, 1)})
    second = frozenset({(1, 2)})
    return [
        make_candidate(
            "pc-via-01", "ch", [0.0, 0.4, 0.2, 0.1], [TRIO_EARLY, first, TRIO_LATE, TRIO_LATE]
        ),
        make_candidate(
            "pc-via-12", "ch", [0.0, 0.4, 0.2, 0.1], [TRIO_EARLY, second, TRIO_LATE, TRIO_LATE]
        ),
    ]


def test_mechanism_diversity_keeps_one_representative_per_tube():
    archive = ParetoArchive("ch", _twin_mechanisms(), ParsConfig())
    assert len(archive.tubes) == 2
    assert archive.representatives == ["pc-via-01", "pc-via-12"]


def test_without_mechanism_diversity_the_cover_ignores_tubes():
    config = ParsConfig(preserve_mechanism_diversity=False)
    archive = ParetoArchive("ch", _twin_mechanisms(), config)
    assert len(archive.tubes) == 2  # tubes are still built and reported
    assert archive.representatives == ["pc-via-01"]
