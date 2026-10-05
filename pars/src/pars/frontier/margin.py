"""Directed distance to the frontier and the adaptive margin (design 13).

Gamma(tau) = min over frontier members phi of max_k [R_k(tau) - R_k(phi)]_+
on normalized objectives. Gamma ~ 0 means on or near the frontier; larger
values mean more strongly dominated. The tier assignment follows design
section 13 with the refinement that CERTAINLY_DOMINATED (interval-separated
from the frontier) outranks a small nominal margin: hard exclusion wins over
rescue.
"""

from __future__ import annotations

import math

from ..config import ParsConfig
from ..schema import ParetoStatus
from .dominance import ObjectiveVector, certainly_dominates, comparable_keys, dominates

TierAssignment = dict[str, ParetoStatus]


def frontier_ids(vectors: dict[str, ObjectiveVector], tolerance: float) -> set[str]:
    """Non-dominated set under nominal epsilon-dominance."""
    frontier = set()
    for cid, vec in vectors.items():
        dominated = any(
            dominates(other, vec, tolerance)
            for other_id, other in vectors.items()
            if other_id != cid
        )
        if not dominated:
            frontier.add(cid)
    return frontier


def directed_margin(vec: ObjectiveVector, frontier_vectors: list[ObjectiveVector]) -> float:
    best = math.inf
    for phi in frontier_vectors:
        keys = comparable_keys(vec, phi)
        if not keys:
            continue
        excess = max(max(vec[k] - phi[k], 0.0) for k in keys)
        best = min(best, excess)
    return best


def assign_tiers(
    vectors: dict[str, ObjectiveVector],
    config: ParsConfig,
    sigma: dict[str, dict[str, float]] | None = None,
) -> tuple[TierAssignment, dict[str, float]]:
    """Return {candidate_id: tier} and {candidate_id: gamma}.

    Uncertainty semantics: claiming CERTAINLY_DOMINATED requires uncertainty
    information for every available objective of both compared candidates,
    or a configured tolerance policy when neither supplies interval data.
    Treating unknown noise as zero would fabricate certainty and turn every
    dominated path into a hard deletion, killing the adaptive margin; the
    design (sections 7/12) forbids exactly that.
    """
    tolerance = config.dominance_tolerance
    frontier = frontier_ids(vectors, tolerance)
    sigma = sigma or {}
    tiers: TierAssignment = {}
    margins: dict[str, float] = {}
    frontier_vectors = [vectors[cid] for cid in sorted(frontier)]
    for cid, vec in vectors.items():
        if cid in frontier:
            tiers[cid] = ParetoStatus.ACTIVE_FRONTIER
            margins[cid] = 0.0
            continue
        margin = directed_margin(vec, frontier_vectors)
        margins[cid] = margin
        if config.dominance_mode == "uncertainty_aware":
            certainly = any(
                certainly_dominates(
                    vectors[phi], vec, sigma.get(phi), sigma.get(cid), tolerance
                )
                for phi in frontier
            )
        else:
            # nominal mode IS strict Pareto pruning (design 38): dominated means
            # hard-excluded, with every dominated candidate dominated by some
            # frontier member in a finite set
            certainly = True
        if certainly:
            tiers[cid] = ParetoStatus.CERTAINLY_DOMINATED
        elif margin <= config.margin_gamma:
            tiers[cid] = ParetoStatus.NEAR_FRONTIER_RESERVE
        else:
            tiers[cid] = ParetoStatus.COLD_ARCHIVE
    return tiers, margins
