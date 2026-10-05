"""Objective-space Pareto cover and mechanism-aware representative selection
(design sections 14-16).

The cover is greedy and deterministic: candidates are visited in priority
order (frontier status, then margin, then id) and a candidate becomes a
representative unless it is already within the symmetric directed excess
tolerance of an earlier representative. Mechanism awareness enters by running
the cover per reaction tube: when preserve_mechanism_diversity is set, two
tubes with near-identical objectives still each keep a representative; when it
is not, a single cover runs over the whole active set.
"""

from __future__ import annotations

import math

from ..config import ParsConfig
from ..schema import ParetoStatus, PathCandidate, ReactionTube
from .dominance import ObjectiveVector, comparable_keys


def directed_excess(a: ObjectiveVector, b: ObjectiveVector) -> float:
    """max_k [a_k - b_k]_+ over identical nonempty availability masks.

    Different masks cannot certify coverage of the missing objectives.
    """
    keys = comparable_keys(a, b)
    if not keys:
        return math.inf
    return max(max(a[k] - b[k], 0.0) for k in keys)


def _symmetric_excess(a: ObjectiveVector, b: ObjectiveVector) -> float:
    return max(directed_excess(a, b), directed_excess(b, a))


def objective_space_cover(
    vectors: dict[str, ObjectiveVector],
    ordered_ids: list[str],
    tolerance: float,
) -> list[str]:
    selected: list[str] = []
    for cid in ordered_ids:
        if not any(
            _symmetric_excess(vectors[cid], vectors[chosen]) <= tolerance for chosen in selected
        ):
            selected.append(cid)
    return selected


def _priority(candidate_id: str, candidates: dict[str, PathCandidate]) -> tuple:
    candidate = candidates[candidate_id]
    status_rank = 0 if candidate.pareto_status is ParetoStatus.ACTIVE_FRONTIER else 1
    margin = candidate.margin_gamma if candidate.margin_gamma is not None else math.inf
    return (status_rank, margin, candidate_id)


def select_representatives(
    candidates: dict[str, PathCandidate],
    tubes: list[ReactionTube],
    vectors: dict[str, ObjectiveVector],
    config: ParsConfig,
) -> list[str]:
    active = [
        cid
        for cid in candidates
        if candidates[cid].pareto_status
        in (ParetoStatus.ACTIVE_FRONTIER, ParetoStatus.NEAR_FRONTIER_RESERVE)
    ]
    if not tubes or not config.preserve_mechanism_diversity:
        # one cover over the whole active set: tubes stay in the record, but two
        # mechanisms with near-identical objectives may share a representative
        ordered = sorted(active, key=lambda cid: _priority(cid, candidates))
        chosen = objective_space_cover(vectors, ordered, config.cover_objective_tolerance)
        for tube in tubes:
            tube.representative_id = next((c for c in chosen if c in tube.member_ids), None)
        return chosen
    representatives: list[str] = []
    for tube in sorted(tubes, key=lambda t: t.tube_id):
        # Tubes retain pruned alternates for later audit/recovery. Only active
        # members are eligible for the current refinement cover.
        ordered = sorted(
            (cid for cid in tube.member_ids if cid in active),
            key=lambda cid: _priority(cid, candidates),
        )
        if not ordered:
            continue
        chosen = objective_space_cover(vectors, ordered, config.cover_objective_tolerance)
        tube.representative_id = chosen[0]
        representatives.extend(chosen)
    return representatives
