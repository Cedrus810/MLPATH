"""The channel-conditioned Pareto archive (design sections 4, 11, 13, 28).

One archive per reaction channel: a superior path in A->B must not erase
A->C, so Pareto pruning never crosses channel boundaries. The archive owns
the full V1 pipeline for one channel: objective evaluation, robust
normalization, tier assignment, tube construction over the active set, cover
selection, and the audit sample that will let the later QM feedback
experiment measure how often Pareto pruning discards unexpected successes.
"""

from __future__ import annotations

import random

from ..config import ParsConfig
from ..frontier.objectives import evaluate_candidate
from ..frontier.normalize import fit_normalization
from ..paths.cluster import build_tubes
from ..schema import ParetoStatus, PathCandidate, ReactionTube
from .cover import select_representatives
from .margin import assign_tiers


class ParetoArchive:
    def __init__(self, channel_id: str, candidates: list[PathCandidate], config: ParsConfig):
        if not candidates:
            raise ValueError(f"channel {channel_id}: empty archive")
        self.channel_id = channel_id
        self.config = config
        self.candidates: dict[str, PathCandidate] = {}
        self.objectives_raw: dict[str, dict] = {}
        self.vectors: dict[str, dict] = {}
        for candidate in candidates:
            if candidate.channel_id != channel_id:
                raise ValueError(
                    f"candidate {candidate.candidate_id} does not belong to "
                    f"channel {channel_id}"
                )
            raw = evaluate_candidate(candidate, config)
            candidate.objectives_raw = raw
            self.candidates[candidate.candidate_id] = candidate
            self.objectives_raw[candidate.candidate_id] = raw
        self.normalization = fit_normalization(self.objectives_raw, config)
        for cid, raw in self.objectives_raw.items():
            vector = self.normalization.normalize(raw)
            self.candidates[cid].objectives_normalized = vector
            self.vectors[cid] = vector
        # A channel may contain several mechanisms with the same endpoints.
        # Establish tubes BEFORE objective pruning: otherwise a dominated but
        # unique mechanism disappears before the diversity policy can protect it.
        all_tubes = (
            build_tubes(
                list(self.candidates.values()), channel_id, config.tube_distance_threshold
            )
            if config.preserve_mechanism_diversity
            else []
        )
        sigma = self._normalized_sigma()
        if all_tubes:
            tiers, margins = {}, {}
            for tube in all_tubes:
                ids = tube.member_ids
                local_tiers, local_margins = assign_tiers(
                    {cid: self.vectors[cid] for cid in ids},
                    config,
                    {cid: sigma[cid] for cid in ids if cid in sigma},
                )
                tiers.update(local_tiers)
                margins.update(local_margins)
        else:
            tiers, margins = assign_tiers(self.vectors, config, sigma)
        self.tiers: dict[str, ParetoStatus] = tiers
        self.margins: dict[str, float] = margins
        for cid, tier in tiers.items():
            self.candidates[cid].pareto_status = tier
            self.candidates[cid].margin_gamma = margins[cid]
        active = [
            self.candidates[cid]
            for cid in sorted(self.candidates)
            if tiers[cid] in (ParetoStatus.ACTIVE_FRONTIER, ParetoStatus.NEAR_FRONTIER_RESERVE)
        ]
        self.tubes: list[ReactionTube] = all_tubes or (
            build_tubes(active, channel_id, config.tube_distance_threshold) if active else []
        )
        for tube in self.tubes:
            for member in tube.member_ids:
                self.candidates[member].tube_id = tube.tube_id
        representatives = select_representatives(
            self.candidates, self.tubes, self.vectors, config
        )
        self.representatives: list[str] = sorted(representatives)
        for cid in self.representatives:
            self.candidates[cid].representative = True

    def _normalized_sigma(self) -> dict[str, dict[str, float]]:
        """Measured objective uncertainties, raw units -> normalized units.

        No V1 objective carries a measured uncertainty, so this is empty unless
        a caller set ``objective_uncertainty`` on the candidates; it is never
        filled in by guessing (design section 7).
        """
        sigma: dict[str, dict[str, float]] = {}
        for cid, candidate in self.candidates.items():
            entries = {
                name: float(value) / self.normalization.scale[name]
                for name, value in candidate.objective_uncertainty.items()
                if value is not None and name in self.normalization.scale
            }
            if entries:
                sigma[cid] = entries
        return sigma

    # -- queries --------------------------------------------------------------

    def refinement_order(self) -> list[str]:
        """Visit one representative per tube before spending on alternates.

        A limited refinement budget must not be consumed by one mechanism's
        id-sorted representatives while another tube has not been attempted.
        No refined energies, saddle labels or future costs enter this order.
        """
        queues = []
        selected = set(self.representatives)
        for tube in sorted(self.tubes, key=lambda t: t.tube_id):
            queue = sorted(selected.intersection(tube.member_ids))
            if queue:
                queues.append(queue)
        return [
            queue[i]
            for i in range(max(map(len, queues), default=0))
            for queue in queues
            if i < len(queue)
        ]

    def frontier(self) -> list[str]:
        return sorted(cid for cid, t in self.tiers.items() if t is ParetoStatus.ACTIVE_FRONTIER)

    def active_set(self, gamma: float | None = None) -> list[str]:
        """F_gamma of design section 13. Monotone in gamma: increasing gamma
        can only enlarge the set (CERTAINLY_DOMINATED is never rescued)."""
        limit = self.config.margin_gamma if gamma is None else gamma
        out = []
        for cid, tier in self.tiers.items():
            if tier is ParetoStatus.ACTIVE_FRONTIER:
                out.append(cid)
            elif tier is ParetoStatus.NEAR_FRONTIER_RESERVE and self.margins[cid] <= limit:
                out.append(cid)
        return sorted(out)

    def audit_sample(self) -> list[str]:
        """Deterministic sample of pruned candidates (design section 20):
        without an audit set the algorithm cannot estimate how often Pareto
        pruning discards unexpected QM successes."""
        pool = sorted(
            cid
            for cid, tier in self.tiers.items()
            if tier in (ParetoStatus.COLD_ARCHIVE, ParetoStatus.CERTAINLY_DOMINATED)
        )
        if self.config.audit_size <= 0 or not pool:
            return []
        rng = random.Random(self.config.audit_seed)
        picked = rng.sample(pool, min(self.config.audit_size, len(pool)))
        return sorted(picked)

    def selection_dict(self) -> dict:
        tiers_by_name: dict[str, list[str]] = {}
        for cid, tier in sorted(self.tiers.items()):
            tiers_by_name.setdefault(tier.value, []).append(cid)
        return {
            "channel_id": self.channel_id,
            "n_candidates": len(self.candidates),
            "gamma": self.config.margin_gamma,
            "dominance_scope": "reaction_tube"
            if self.config.preserve_mechanism_diversity
            else "channel",
            "tiers": tiers_by_name,
            "tubes": [tube.to_dict() for tube in self.tubes],
            "representatives": self.representatives,
            "refinement_order": self.refinement_order(),
            "audit": self.audit_sample(),
            "normalization": self.normalization.to_dict(),
            "objectives_raw": dict(self.objectives_raw),
            "objectives_normalized": dict(self.vectors),
        }
