"""Arc-length reparameterization (design section 17.4).

A direct upgrade of the seed-path export that used to live in the PRRS CLI:
identical arc-length logic, plus explicit rejection of zero-length segments.
The metric is plain Cartesian. Rigid motion is removed separately, by
``align_frames`` before resampling, and only for systems whose energy is
invariant under rigid motion (``refinement_remove_rigid_motion``): aligning
every frame onto the first keeps one globally consistent frame, but in an
external field or a fixed embedding it would move the system relative to it.
"""

from __future__ import annotations

import numpy as np

_ARC_EPSILON = 1e-12


def resample_positions(positions: np.ndarray, images: int) -> np.ndarray:
    """Redistribute ``positions`` (T, N, 3) onto ``images`` approximately
    uniform arc-length samples. The first and last frames are preserved
    exactly."""
    if images < 3:
        raise ValueError("At least three images required")
    positions = np.asarray(positions, dtype=float)
    if positions.ndim != 3 or positions.shape[0] < 2:
        raise ValueError("positions must be (T, N, 3) with T >= 2")
    flat = positions.reshape(positions.shape[0], -1)
    segments = np.linalg.norm(np.diff(flat, axis=0), axis=1)
    arc = np.concatenate([[0.0], np.cumsum(segments)])
    keep = np.concatenate([[True], segments > _ARC_EPSILON])
    if int(keep.sum()) < 2:
        raise ValueError("Trajectory has zero path length")
    arc, flat = arc[keep], flat[keep]
    targets = np.linspace(arc[0], arc[-1], images)
    interpolated = np.column_stack(
        [np.interp(targets, arc, flat[:, i]) for i in range(flat.shape[1])]
    )
    n_atoms = positions.shape[1]
    return np.ascontiguousarray(interpolated.reshape(images, n_atoms, 3))


def kabsch_align(frame: np.ndarray, reference: np.ndarray) -> np.ndarray:
    """Rigidly superimpose ``frame`` (N, 3) onto ``reference``; proper rotations only."""
    p = frame - frame.mean(axis=0)
    q = reference - reference.mean(axis=0)
    correlation = p.T @ q
    u, _, vt = np.linalg.svd(correlation)
    correction = np.eye(3)
    correction[-1, -1] = np.sign(np.linalg.det(u @ vt))
    # rows are atoms: p @ (U D V^T) maps p onto q (row-vector form of Kabsch)
    rotation = u @ correction @ vt
    return p @ rotation + reference.mean(axis=0)


def align_frames(positions: np.ndarray) -> np.ndarray:
    """Every frame of ``positions`` (T, N, 3) rigidly superimposed on the first."""
    positions = np.asarray(positions, dtype=float)
    return np.asarray([kabsch_align(frame, positions[0]) for frame in positions])
