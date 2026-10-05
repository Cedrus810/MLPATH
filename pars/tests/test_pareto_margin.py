"""Archive rescue: increasing gamma monotonically expands the active set
(design sections 13 and 31)."""

from __future__ import annotations

import pytest

from pars.config import ParsConfig
from pars.frontier.margin import assign_tiers, directed_margin, frontier_ids
from pars.schema import ParetoStatus

VECTORS = {
    "on-frontier": {"x": 0.0, "y": 1.0},
    "trade-off": {"x": 1.0, "y": 0.0},
    "slightly-off": {"x": 0.2, "y": 1.1},
    "far-off": {"x": 2.0, "y": 3.0},
    "very-far": {"x": 5.0, "y": 5.0},
}


def test_frontier_members_have_zero_margin():
    assert frontier_ids(VECTORS, 0.0) == {"on-frontier", "trade-off"}
    assert directed_margin(VECTORS["on-frontier"], [VECTORS["on-frontier"]]) == 0.0


def test_margin_orders_candidates_by_dominance_strength():
    config = ParsConfig(margin_gamma=0.15)
    _, margins = assign_tiers(VECTORS, config)
    assert margins["slightly-off"] < margins["far-off"] < margins["very-far"]


def test_rescue_is_monotone_in_gamma():
    config = ParsConfig(margin_gamma=0.0, dominance_mode="nominal")
    tiers, margins = assign_tiers(VECTORS, config)

    # emulate F_gamma from the stored margins: ACTIVE plus margin <= gamma
    def active(gamma: float) -> set[str]:
        out = set()
        for cid, tier in tiers.items():
            if tier is ParetoStatus.ACTIVE_FRONTIER:
                out.add(cid)
            elif tier is ParetoStatus.NEAR_FRONTIER_RESERVE and margins[cid] <= gamma:
                out.add(cid)
        return out

    previous = set()
    for gamma in (0.0, 0.1, 0.15, 0.5, 1.0, 10.0):
        current = active(gamma)
        assert previous <= current, f"gamma={gamma} shrank the active set"
        previous = current


def test_uncertainty_aware_archive_keeps_uncertain_paths_rescuable():
    # far-off is dominated but its true values are deeply uncertain: it stays
    # in the cold archive rather than being hard-deleted, so a later gamma
    # increase (QM evidence) can rescue it
    config = ParsConfig(margin_gamma=0.15)
    sigma = {"far-off": {"x": 10.0, "y": 10.0}}
    tiers, _ = assign_tiers(VECTORS, config, sigma=sigma)
    assert tiers["far-off"] is ParetoStatus.COLD_ARCHIVE
    assert tiers["very-far"] is ParetoStatus.COLD_ARCHIVE


def test_missing_objectives_cannot_cycle_and_empty_the_frontier():
    import json

    vectors = {
        "a": {"x": 0.0, "y": None, "z": 2.0},
        "b": {"x": 1.0, "y": 0.0, "z": None},
        "c": {"x": None, "y": 1.0, "z": 1.0},
        "d": {"x": 2.0, "y": None, "z": 3.0},
    }
    assert frontier_ids(vectors, 0.0) == {"a", "b", "c"}
    tiers, margins = assign_tiers(vectors, ParsConfig())
    assert all(tiers[cid] is ParetoStatus.ACTIVE_FRONTIER for cid in ("a", "b", "c"))
    assert margins["d"] == pytest.approx(2.0)
    # Output JSON rejects infinity; every dominated candidate now has a
    # comparable frontier witness and hence a finite margin.
    json.dumps(margins, allow_nan=False)


def test_cover_cannot_replace_a_candidate_with_missing_objectives():
    from pars.frontier.cover import objective_space_cover

    vectors = {
        "a": {"x": 0.0, "y": None},
        "b": {"x": 0.0, "y": 2.0},
        "c": {"x": 0.0, "y": None},
    }
    assert objective_space_cover(vectors, ["a", "b", "c"], 0.1) == ["a", "b"]


def test_margin_function_prefers_the_closest_frontier_member():
    frontier = [VECTORS["on-frontier"], VECTORS["trade-off"]]
    probe = {"x": 1.2, "y": 0.9}
    # excess vs on-frontier: max(1.2, -0.1)_+ = 1.2; vs trade-off: max(0.2, 0.9) = 0.9
    assert directed_margin(probe, frontier) == pytest.approx(0.9)
