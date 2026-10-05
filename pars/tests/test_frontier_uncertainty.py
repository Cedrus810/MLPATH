"""Uncertainty-aware dominance and tiering (design sections 12-13)."""

from __future__ import annotations

import pytest

from pars.config import ParsConfig
from pars.frontier.dominance import certainly_dominates, dominates
from pars.frontier.margin import assign_tiers
from pars.schema import ParetoStatus

GOOD = {"x": 0.9, "y": 0.9}
BAD = {"x": 1.0, "y": 1.0}


def test_nominal_dominance_is_epsilon_dominance():
    assert dominates(GOOD, BAD)
    assert not dominates(BAD, GOOD)
    # an equal vector never strictly dominates itself
    assert dominates({"x": 1.0, "y": 1.0}, BAD, tolerance=0.05) is False


def test_dominance_requires_matching_availability_masks():
    partial_a = {"x": 0.5, "y": None}
    partial_b = {"x": 0.6, "y": 2.0}
    assert not dominates(partial_a, partial_b)
    assert dominates(partial_a, {"x": 0.6, "y": None})
    assert dominates(partial_a, {"x": 0.6})
    # no jointly available objective -> no dominance
    assert not dominates({"x": None}, {"x": None})


def test_certain_dominance_requires_interval_separation():
    # nominal: GOOD dominates BAD; but BAD's values are very uncertain, so the
    # separation is NOT certain
    assert not certainly_dominates(GOOD, BAD)
    assert certainly_dominates(
        GOOD, BAD, sigma_a={"x": 0.01, "y": 0.01}, sigma_b={"x": 0.01, "y": 0.01}
    )
    assert not certainly_dominates(GOOD, BAD, sigma_b={"x": 0.5, "y": 0.5})
    # with a's own uncertainty small and b's large, still uncertain
    assert not certainly_dominates(
        GOOD, BAD, sigma_a={"x": 0.02, "y": 0.02}, sigma_b={"x": 0.5, "y": 0.5}
    )
    # tolerance is extra required separation: it can only make a hard deletion
    # harder, never rescue an uncertain one into deletion
    assert certainly_dominates(GOOD, BAD, tolerance=0.05)
    assert not certainly_dominates(GOOD, BAD, tolerance=0.2)
    assert not certainly_dominates(GOOD, BAD, sigma_b={"x": 0.5, "y": 0.5}, tolerance=0.5)


def test_tolerance_never_empties_the_frontier():
    # within-tolerance trade-offs used to dominate each other mutually, leaving
    # no frontier and therefore no representative at all
    vectors = {
        "a": {"x": 0.0, "y": 0.05},
        "b": {"x": 0.05, "y": 0.0},
        "c": {"x": 0.5, "y": 0.5},
    }
    tiers, _ = assign_tiers(vectors, ParsConfig(dominance_tolerance=0.1))
    assert tiers["a"] is ParetoStatus.ACTIVE_FRONTIER
    assert tiers["b"] is ParetoStatus.ACTIVE_FRONTIER


def test_tolerance_widens_the_frontier():
    # b is worse than a by only 0.05 everywhere: dominated at tolerance 0,
    # indistinguishable (so kept on the frontier) at tolerance 0.1
    vectors = {"a": {"x": 0.0, "y": 0.0}, "b": {"x": 0.05, "y": 0.05}}
    strict, _ = assign_tiers(vectors, ParsConfig(dominance_tolerance=0.0))
    wide, _ = assign_tiers(vectors, ParsConfig(dominance_tolerance=0.1))
    assert strict["b"] is not ParetoStatus.ACTIVE_FRONTIER
    assert wide["b"] is ParetoStatus.ACTIVE_FRONTIER
    assert not dominates(vectors["a"], vectors["b"], tolerance=0.1)


def test_uncertain_dominated_path_archives_but_is_not_hard_deleted():
    vectors = {"good": GOOD, "bad": BAD}
    sigma = {"bad": {"x": 0.5, "y": 0.5}}
    # gamma below the nominal margin 0.1 so the uncertain path cannot hide in
    # the near-frontier reserve either
    config = ParsConfig(dominance_mode="uncertainty_aware", margin_gamma=0.05)
    tiers, _ = assign_tiers(vectors, config, sigma=sigma)
    assert tiers["good"] is ParetoStatus.ACTIVE_FRONTIER
    assert tiers["bad"] is ParetoStatus.COLD_ARCHIVE  # archived, recoverable

    # without any measured uncertainty and no tolerance, certainty cannot be
    # claimed at all: the path stays recoverable rather than hard-deleted
    tiers_no_sigma, _ = assign_tiers(vectors, config)
    assert tiers_no_sigma["bad"] is ParetoStatus.COLD_ARCHIVE


def test_measured_uncertainty_enables_hard_tiering():
    vectors = {"front": {"x": 0.0, "y": 0.0}, "near": {"x": 0.05, "y": 0.05}}
    config = ParsConfig(margin_gamma=10.0)
    sigma = {
        "front": {"x": 0.01, "y": 0.01},
        "near": {"x": 0.01, "y": 0.01},
    }  # complete measured intervals for this pair
    tiers, margins = assign_tiers(vectors, config, sigma=sigma)
    # interval separation is measured, so the hard tier wins even at large gamma
    assert tiers["near"] is ParetoStatus.CERTAINLY_DOMINATED
    assert margins["near"] == pytest.approx(0.05)


def test_tolerance_substitutes_for_missing_uncertainty():
    vectors = {"front": {"x": 0.0, "y": 0.0}, "near": {"x": 0.05, "y": 0.05}}
    config = ParsConfig(margin_gamma=10.0, dominance_tolerance=0.04)
    tiers, _ = assign_tiers(vectors, config)
    assert tiers["near"] is ParetoStatus.CERTAINLY_DOMINATED


def test_uncertain_near_frontier_path_stays_rescuable():
    vectors = {"front": {"x": 0.0, "y": 0.0}, "near": {"x": 0.05, "y": 0.05}}
    config = ParsConfig(margin_gamma=10.0)
    sigma = {"near": {"x": 2.0, "y": 2.0}}
    tiers, _ = assign_tiers(vectors, config, sigma=sigma)
    # near's true value might beat the frontier; it must remain recoverable
    assert tiers["near"] is ParetoStatus.NEAR_FRONTIER_RESERVE


def test_nominal_mode_is_strict_pareto_pruning():
    vectors = {"good": GOOD, "bad": BAD}
    config = ParsConfig(dominance_mode="nominal")
    sigma = {"bad": {"x": 0.5, "y": 0.5}}
    tiers, _ = assign_tiers(vectors, config, sigma=sigma)
    assert tiers["bad"] is ParetoStatus.CERTAINLY_DOMINATED


def test_archive_uses_measured_objective_uncertainty():
    # objective_uncertainty is in raw units; the archive must carry it into the
    # normalized tiering, where a measured, tiny sigma licenses hard deletion
    from pars.frontier.archive import ParetoArchive
    from pars.frontier.objectives import evaluate_candidate

    from conftest import TRIO_EARLY, TRIO_LATE, make_candidate

    def channel(with_sigma: bool):
        candidates = []
        for name, peak in (("pc-a", 0.2), ("pc-b", 0.5), ("pc-c", 1.0)):
            candidate = make_candidate(
                name, "ch", [0.0, peak, 0.1, 0.0], [TRIO_EARLY, TRIO_LATE, TRIO_LATE, TRIO_LATE]
            )
            if with_sigma:
                # This analytic fixture explicitly supplies exact (zero-width)
                # intervals for its deterministic geometric objectives.
                candidate.objective_uncertainty = {
                    key: 0.0
                    for key, value in evaluate_candidate(candidate, ParsConfig()).items()
                    if value is not None
                }
                candidate.objective_uncertainty["energy_excursion"] = 0.01
            candidates.append(candidate)
        return ParetoArchive("ch", candidates, ParsConfig(margin_gamma=10.0))

    assert channel(False).tiers["pc-c"] is ParetoStatus.NEAR_FRONTIER_RESERVE
    assert channel(True).tiers["pc-c"] is ParetoStatus.CERTAINLY_DOMINATED


@pytest.mark.parametrize("tolerance", [0.0, 0.05])
def test_unrelated_partial_sigma_does_not_hard_delete_missing_pair(tolerance):
    vectors = {
        "good": {"x": 0.0, "y": 0.0},
        "bad": {"x": 1.0, "y": 1.0},
        "noisy": {"x": 2.0, "y": 2.0},
    }
    config = ParsConfig(margin_gamma=10.0, dominance_tolerance=tolerance)
    before, _ = assign_tiers(vectors, config)
    after, _ = assign_tiers(vectors, config, {"noisy": {"x": 0.1}})
    assert after["bad"] is before["bad"]
    assert after["noisy"] is ParetoStatus.NEAR_FRONTIER_RESERVE


@pytest.mark.parametrize("tolerance", [0.0, 0.04])
def test_partial_intervals_cannot_establish_certain_dominance(tolerance):
    assert not certainly_dominates(
        GOOD,
        BAD,
        sigma_a={"x": 0.001, "y": 0.001},
        sigma_b={"x": 0.001},
        tolerance=tolerance,
    )
    assert not certainly_dominates(
        GOOD,
        BAD,
        sigma_a={"x": 0.001, "y": 0.001},
        tolerance=tolerance,
    )
