"""A distinct mechanism must survive before objective-space compression."""

from pars.config import ParsConfig
from pars.frontier.archive import ParetoArchive
from pars.schema import ParetoStatus

from conftest import TRIO_EARLY, TRIO_LATE, make_candidate


def population():
    return [
        make_candidate(
            "cheap",
            "ch",
            [0.0, 0.1, 0.0, 0.0],
            [TRIO_EARLY, TRIO_LATE, TRIO_LATE, TRIO_LATE],
            amplitude=0.1,
        ),
        make_candidate(
            "rare",
            "ch",
            [0.0, 3.0, 2.0, 0.0],
            [TRIO_EARLY, TRIO_LATE, TRIO_EARLY, TRIO_LATE],
            amplitude=1.0,
        ),
        make_candidate(
            "rare-worse",
            "ch",
            [0.0, 4.0, 3.0, 0.0],
            [TRIO_EARLY, TRIO_LATE, TRIO_EARLY, TRIO_LATE],
            amplitude=2.0,
        ),
    ]


def test_dominated_distinct_mechanism_survives_and_its_alternates_remain():
    archive = ParetoArchive("ch", population(), ParsConfig(margin_gamma=0.0))
    assert archive.tiers["rare"] is ParetoStatus.ACTIVE_FRONTIER
    assert "rare" in archive.representatives
    assert "rare-worse" not in archive.representatives
    tube = next(t for t in archive.tubes if "rare" in t.member_ids)
    assert "rare-worse" in tube.member_ids
    assert archive.selection_dict()["dominance_scope"] == "reaction_tube"


def test_explicit_channel_pruning_retains_the_counterfactual():
    archive = ParetoArchive(
        "ch", population(), ParsConfig(margin_gamma=0.0, preserve_mechanism_diversity=False)
    )
    assert "rare" not in archive.representatives
    assert archive.selection_dict()["dominance_scope"] == "channel"


def test_limited_budget_visits_every_tube_before_an_alternate():
    # Budget ordering is independent of the number of selected alternates
    # and carries no refined barrier/oracle information.
    archive = ParetoArchive(
        "ch", population(), ParsConfig(margin_gamma=2.5, cover_objective_tolerance=0.0)
    )
    order = archive.refinement_order()
    assert set(order) == set(archive.representatives)
    membership = {cid: tube.tube_id for tube in archive.tubes for cid in tube.member_ids}
    assert len({membership[cid] for cid in order[: len(archive.tubes)]}) == len(archive.tubes)
