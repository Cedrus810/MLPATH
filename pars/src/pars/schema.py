"""Core data objects for PARS.

The lifecycle follows the design document (section 6):

    DiscoveryTrace -> PathCandidate -> ReactionTube -> RefinedPath -> ValidatedPath

Everything here is deliberately plain: dataclasses, numpy arrays and dicts.
PARS never imports ``prrs``; the run directory is the contract.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum

import numpy as np


class ParetoStatus(str, Enum):
    """The four tiers of design section 13.

    Raising gamma can only promote a candidate (archive rescue); only
    CERTAINLY_DOMINATED is hard-excluded from rescue, and even those paths stay
    in the archive on disk -- uncertain dominated paths go to an archive, not
    to /dev/null (design section 12).
    """

    ACTIVE_FRONTIER = "active_frontier"
    NEAR_FRONTIER_RESERVE = "near_frontier_reserve"
    COLD_ARCHIVE = "cold_archive"
    CERTAINLY_DOMINATED = "certainly_dominated"


class RefinementStatus(str, Enum):
    NOT_REFINED = "not_refined"
    NEB_CONVERGED = "neb_converged"
    CI_NEB_CONVERGED = "ci_neb_converged"
    UNRESOLVED_MLP = "unresolved_mlp"
    FAILED = "failed"


class ValidationLevel(str, Enum):
    """Design section 6.5; these stages must never be conflated."""

    NONE = "none"
    MLP_ONLY = "mlp_only"
    QM_SINGLE_POINT_CHECKED = "qm_single_point_checked"
    QM_GRADIENT_CHECKED = "qm_gradient_checked"
    QM_TS_OPTIMIZED = "qm_ts_optimized"
    HESSIAN_VERIFIED = "hessian_verified"
    IRC_VERIFIED = "irc_verified"


class ElectronicStateStatus(str, Enum):
    """Design section 23. Ordinary PARS assumes a meaningful scalar PES; the
    default here is UNKNOWN so that multistate suspicion can only be claimed
    explicitly, never silently assumed away."""

    ORDINARY = "ordinary"
    UNKNOWN = "unknown"
    SUSPECTED_CROSSING = "suspected_crossing"
    MULTISTATE_REQUIRED = "multistate_required"


@dataclass(frozen=True)
class Provenance:
    """Complete traceability back to the original PRRS trial (design 6.2/27).

    No compressed or processed path may lose its provenance: every candidate
    keeps the run directory, the attempt id, the frame window of the original
    response trajectory, the perturbation that triggered it, and the model
    that produced it.
    """

    run_dir: str
    attempt_id: str
    channel_id: str | None
    event_key: str | None
    reactant_state_id: str | None
    product_state_id: str | None
    frame_start: int
    frame_stop: int  # exclusive; index into the original response trajectory
    perturbation_family: str | None = None
    perturbation_amplitude: float | None = None
    perturbation_unit: str | None = None
    seed: int | None = None
    implementation_sha256: str | None = None
    calculator: dict = field(default_factory=dict)
    # {"charge": int, "spin": multiplicity}: a PES input like the model itself.
    # None only for hand-built candidates that never came from a PRRS run.
    electronic_state: dict | None = None
    extra: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        return {
            "run_dir": self.run_dir,
            "attempt_id": self.attempt_id,
            "channel_id": self.channel_id,
            "event_key": self.event_key,
            "reactant_state_id": self.reactant_state_id,
            "product_state_id": self.product_state_id,
            "frame_start": self.frame_start,
            "frame_stop": self.frame_stop,
            "perturbation_family": self.perturbation_family,
            "perturbation_amplitude": self.perturbation_amplitude,
            "perturbation_unit": self.perturbation_unit,
            "seed": self.seed,
            "implementation_sha256": self.implementation_sha256,
            "calculator": dict(self.calculator),
            "electronic_state": (
                None if self.electronic_state is None else dict(self.electronic_state)
            ),
            "extra": dict(self.extra),
        }

    @classmethod
    def from_dict(cls, payload: dict) -> Provenance:
        known = {f for f in cls.__dataclass_fields__}
        missing = {"run_dir", "attempt_id", "frame_start", "frame_stop"} - set(payload)
        if missing:
            raise ValueError(f"provenance missing required keys: {sorted(missing)}")
        return cls(
            **{k: v for k, v in payload.items() if k in known},
        )


@dataclass(frozen=True)
class DiscoveryTrace:
    """Raw trajectory produced by PRRS; not yet an NEB path (design 6.1).

    Contains the forced perturbation frames, the free response, quench
    iterations and post-reaction relaxation. Phase labels come from the run's
    own observations, so trimming decisions downstream are authoritative
    rather than inferred.
    """

    trace_id: str
    numbers: np.ndarray  # (N,) atomic numbers
    positions: np.ndarray  # (T, N, 3) all recorded frames, Angstrom
    energies: np.ndarray  # (T,) potential energy, eV
    phases: tuple[str, ...]  # (T,) per-frame phase labels
    signatures: tuple[frozenset[tuple[int, int]], ...] | None  # (T,) bond sets
    uncertainty: np.ndarray | None  # (T,) eV/A; NaN where never measured
    endpoint_positions: np.ndarray | None  # (N, 3) product endpoint snapshot
    endpoint_energy: float | None
    provenance: Provenance

    def __post_init__(self) -> None:
        t = len(self.phases)
        if self.positions.shape[0] != t:
            raise ValueError(
                f"trace {self.trace_id}: {self.positions.shape[0]} frames but {t} phase labels"
            )
        if self.energies.shape[0] != t:
            raise ValueError(f"trace {self.trace_id}: energies must have {t} entries")


@dataclass
class PathCandidate:
    """A cleaned reaction segment A -> transition region -> B (design 6.2).

    The artificial forcing stage is excluded at construction; its parameters
    survive in provenance. Objectives, pareto status and refinement state
    accumulate as the candidate moves through the pipeline.
    """

    candidate_id: str
    channel_id: str
    trace_id: str
    numbers: np.ndarray
    positions: np.ndarray  # (T, N, 3) trimmed frames
    energies: np.ndarray  # (T,) eV
    uncertainty: np.ndarray | None  # (T,) eV/A; NaN where never measured
    signatures: tuple[frozenset[tuple[int, int]], ...] | None
    provenance: Provenance

    objectives_raw: dict[str, float | None] = field(default_factory=dict)
    objectives_normalized: dict[str, float | None] = field(default_factory=dict)
    objective_uncertainty: dict[str, float | None] = field(default_factory=dict)
    pareto_status: ParetoStatus | None = None
    margin_gamma: float | None = None
    tube_id: str | None = None
    representative: bool = False
    refinement_status: RefinementStatus = RefinementStatus.NOT_REFINED
    validation_level: ValidationLevel = ValidationLevel.NONE
    electronic_state_status: ElectronicStateStatus = ElectronicStateStatus.UNKNOWN

    def to_summary(self) -> dict:
        """JSON-safe record without frame arrays (frames stay in extxyz)."""
        return {
            "candidate_id": self.candidate_id,
            "channel_id": self.channel_id,
            "trace_id": self.trace_id,
            "n_frames": int(self.positions.shape[0]),
            "n_atoms": int(self.positions.shape[1]),
            "objectives_raw": dict(self.objectives_raw),
            "objectives_normalized": dict(self.objectives_normalized),
            "objective_uncertainty": dict(self.objective_uncertainty),
            "pareto_status": (
                self.pareto_status.value if self.pareto_status is not None else None
            ),
            "margin_gamma": self.margin_gamma,
            "tube_id": self.tube_id,
            "representative": self.representative,
            "refinement_status": self.refinement_status.value,
            "validation_level": self.validation_level.value,
            "electronic_state_status": self.electronic_state_status.value,
            "provenance": self.provenance.to_dict(),
        }


@dataclass
class ReactionTube:
    """Cluster of mechanistically equivalent PathCandidates (design 6.3/16)."""

    tube_id: str
    channel_id: str
    member_ids: list[str] = field(default_factory=list)
    representative_id: str | None = None

    def to_dict(self) -> dict:
        return {
            "tube_id": self.tube_id,
            "channel_id": self.channel_id,
            "member_ids": list(self.member_ids),
            "representative_id": self.representative_id,
        }


@dataclass
class RefinedPath:
    """Result of refining one candidate (design 6.4/18).

    Barrier values are meaningful on the MLP PES only. Missing quantities stay
    None; they are never fabricated (design section 7).
    """

    candidate_id: str
    backend: str
    status: RefinementStatus
    numbers: np.ndarray | None = None  # (N,) atomic numbers
    positions: np.ndarray | None = None  # (M, N, 3) refined chain
    energies: np.ndarray | None = None  # (M,) eV
    barrier_eV: float | None = None
    max_image_force_eV_A: float | None = None
    neb_residual_eV_A: float | None = None
    uncertainty_max_eV_A: float | None = None
    climbing_image: bool = False
    message: str = ""
    provenance: dict = field(default_factory=dict)

    def to_summary(self) -> dict:
        return {
            "candidate_id": self.candidate_id,
            "backend": self.backend,
            "status": self.status.value,
            "barrier_eV": self.barrier_eV,
            "max_image_force_eV_A": self.max_image_force_eV_A,
            "neb_residual_eV_A": self.neb_residual_eV_A,
            "uncertainty_max_eV_A": self.uncertainty_max_eV_A,
            "climbing_image": self.climbing_image,
            "message": self.message,
            "provenance": dict(self.provenance),
        }
