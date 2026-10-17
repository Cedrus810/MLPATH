"""Cheap pre-refinement objectives (design section 8).

All objectives are transformed so that smaller is better. A quantity that the
run artifacts do not carry is reported as None and never fabricated (design
section 7): None values are skipped by normalization and dominance, and the
selection record states them as unavailable.
"""

from __future__ import annotations

import numpy as np

from ..config import ParsConfig
from ..schema import PathCandidate
from ..paths.fingerprint import recrossing_count
from ..paths.resample import kabsch_align

OBJECTIVE_NAMES = (
    "energy_excursion",
    "trigger_severity",
    "mlp_uncertainty",
    "path_curvature",
    "recrossing",
    "endpoint_penalty",
)


def energy_excursion(energies: np.ndarray) -> float:
    """max_t E(t) - E_A over the candidate window (design 8.1).

    E_A is the first frame, i.e. the source-phase minimum. This is a cheap
    path descriptor, NOT a validated barrier.
    """
    if len(energies) == 0:
        raise ValueError("energy excursion requires frames")
    return float(np.max(energies) - energies[0])


def trigger_severity(candidate: PathCandidate) -> float | None:
    """Required perturbation amplitude (design 8.2). Not a free-energy barrier."""
    return candidate.provenance.perturbation_amplitude


def trigger_severity_key(candidate: PathCandidate) -> str:
    """Objective name the amplitude is reported under.

    Design 8.2 compares amplitudes only "under the specified perturbation
    family": 0.5 Angstrom of stretch and 0.2 eV of kick share no scale, and
    median/IQR normalization of one mixed column does not give them one.
    Keying the objective by family and unit gives candidates of different
    families different availability masks, so dominance and cover treat them
    as incomparable instead of ranking them by the bare number.
    """
    family = candidate.provenance.perturbation_family
    unit = candidate.provenance.perturbation_unit
    if family is None or unit is None:
        return "trigger_severity"
    return f"trigger_severity[{family}/{unit}]"


def mlp_uncertainty(uncertainty: np.ndarray | None, quantile: float) -> float | None:
    """Per-frame force uncertainty aggregated at a high quantile (design 8.3).

    The default quantile is 1.0 (the worst frame): a single unreliable
    transition-region frame must not hide behind a quiet mean.
    """
    if uncertainty is None:
        return None
    finite = uncertainty[np.isfinite(uncertainty)]
    if finite.size == 0:
        return None
    return float(np.quantile(finite, quantile))


def path_curvature(positions: np.ndarray, images: int) -> float:
    """Integral of |d2R/ds2|^2 ds along the path (design 8.4), in 1/Angstrom,
    measured at NEB image resolution.

    Each frame is Kabsch-aligned to the first frame first, so rigid drift of
    the molecule does not count as curvature. This is an alignment-aware
    approximation; symmetry-aware alignment stays authoritative on the PRRS
    side and can be injected later.

    The aligned path is resampled to ``images`` uniform arc-length points (the
    NEB image resolution; the refinement seed itself also keeps the source
    frame and is not aligned) before the finite difference. On raw MD
    frames the integral picks up every vibration the sampling resolves, so it
    grew with the number of frames and duplicated recrossing; at a fixed image
    count it measures the shape of the seed instead.
    """
    positions = np.asarray(positions, dtype=float)
    if positions.shape[0] < 3:
        raise ValueError("path curvature requires at least three frames")
    if images < 3:
        raise ValueError("path curvature requires at least three images")
    aligned = np.asarray([kabsch_align(frame, positions[0]) for frame in positions])
    flat = aligned.reshape(aligned.shape[0], -1)
    segments = np.linalg.norm(np.diff(flat, axis=0), axis=1)
    arc = np.concatenate([[0.0], np.cumsum(segments)])
    keep = np.concatenate([[True], segments > 1e-12])
    if int(keep.sum()) < 3:
        return 0.0
    arc, flat = arc[keep], flat[keep]
    grid = np.linspace(arc[0], arc[-1], images)
    sampled = np.column_stack([np.interp(grid, arc, flat[:, i]) for i in range(flat.shape[1])])
    spacing = grid[1] - grid[0]
    # d2R/ds2 ~ second / h^2, so the integral sum(|d2R/ds2|^2 h) is sum|second|^2 / h^3
    second = sampled[2:] - 2.0 * sampled[1:-1] + sampled[:-2]
    return float(np.sum(np.linalg.norm(second, axis=1) ** 2) / spacing**3)


def endpoint_penalty(candidate: PathCandidate) -> None:
    """Endpoint basin-confidence penalty (design 8.6).

    Not available in the run artifacts V1 consumes, so it stays None rather
    than being fabricated. The state classification that would feed this
    objective remains authoritative on the PRRS side.
    """
    return None


def recrossing(candidate: PathCandidate) -> float | None:
    """State crossings beyond the initial departure (design 8.5)."""
    if candidate.signatures is None:
        return None
    return float(recrossing_count(candidate.signatures))


def free_segment_start(candidate: PathCandidate) -> int:
    """Index of the first frame of the trailing contiguous run of kept frames.

    The candidate window is not contiguous: the forcing stage sits between the
    source frame and free_start (``provenance.extra["kept_frames"]``). The
    jump across that gap is the perturbation, not the shape of the path, so
    shape objectives are measured on the free segment after it. Without
    kept_frames the window is taken as contiguous.
    """
    kept = candidate.provenance.extra.get("kept_frames")
    if not kept:
        return 0
    start = len(kept) - 1
    while start > 0 and kept[start - 1] == kept[start] - 1:
        start -= 1
    return start


def evaluate_candidate(candidate: PathCandidate, config: ParsConfig) -> dict:
    segment = candidate.positions[free_segment_start(candidate) :]
    curvature = (
        path_curvature(segment, config.refinement_images) if segment.shape[0] >= 3 else None
    )
    severity_key = trigger_severity_key(candidate)
    # An amplitude without its family and unit has no comparison group.
    severity = trigger_severity(candidate) if severity_key != "trigger_severity" else None
    return {
        "energy_excursion": energy_excursion(candidate.energies),
        severity_key: severity,
        "mlp_uncertainty": mlp_uncertainty(candidate.uncertainty, config.uncertainty_quantile),
        "path_curvature": curvature,
        "recrossing": recrossing(candidate),
        "endpoint_penalty": endpoint_penalty(candidate),
    }
