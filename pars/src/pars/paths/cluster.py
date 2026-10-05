"""Reaction tubes: clusters of mechanistically equivalent candidates
(design section 16).

V1 uses deterministic leader clustering in candidate-id order: the first
candidate of a tube is its leader, later candidates join the first tube whose
leader is within the signature-distance threshold. Tubes are built per channel
over the near-frontier active set; every tube keeps its alternates so that a
failed refinement of one representative does not destroy the mechanism.
"""

from __future__ import annotations

from ..schema import PathCandidate, ReactionTube
from .fingerprint import frechet_path_distance


def build_tubes(
    candidates: list[PathCandidate], channel_id: str, threshold: float
) -> list[ReactionTube]:
    missing = [c.candidate_id for c in candidates if c.signatures is None]
    if missing:
        raise ValueError(
            f"tube clustering requires bond signatures; candidates without "
            f"signatures: {missing[:5]}"
        )
    if threshold <= 0 or threshold > 1:
        raise ValueError("tube threshold must be in (0, 1]")
    tubes: list[ReactionTube] = []
    leader_signatures: list[tuple] = []
    distance_cache: dict[tuple[str, str], float] = {}
    for candidate in sorted(candidates, key=lambda c: c.candidate_id):
        placed = False
        for index, tube in enumerate(tubes):
            pair = (tube.member_ids[0], candidate.candidate_id)
            if pair not in distance_cache:
                distance_cache[pair] = frechet_path_distance(
                    leader_signatures[index], candidate.signatures
                )
            if distance_cache[pair] <= threshold:
                tube.member_ids.append(candidate.candidate_id)
                placed = True
                break
        if not placed:
            tubes.append(
                ReactionTube(
                    tube_id=f"tube-{channel_id}-{len(tubes):04d}",
                    channel_id=channel_id,
                    member_ids=[candidate.candidate_id],
                )
            )
            leader_signatures.append(candidate.signatures)
    return tubes
