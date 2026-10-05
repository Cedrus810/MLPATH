from __future__ import annotations

import numpy as np
import pytest

from pars.paths.resample import resample_positions


def _line(frames: int, n_atoms: int = 2) -> np.ndarray:
    return (
        np.linspace(0.0, 1.0, frames)[:, None, None]
        * np.ones((1, n_atoms, 3))
        * np.asarray([1.0, 2.0, 3.0])
    )


def test_endpoints_preserved_and_count_exact():
    positions = _line(9)
    chain = resample_positions(positions, 5)
    assert chain.shape == (5, 2, 3)
    assert np.allclose(chain[0], positions[0])
    assert np.allclose(chain[-1], positions[-1])


def test_images_must_be_at_least_three():
    with pytest.raises(ValueError, match="three images"):
        resample_positions(_line(6), 2)


def test_zero_path_length_rejected():
    frames = np.zeros((5, 2, 3))
    with pytest.raises(ValueError, match="zero path length"):
        resample_positions(frames, 4)


def test_duplicate_frames_removed_without_error():
    positions = _line(6)
    positions = np.concatenate([positions[:3], positions[2:]], axis=0)  # duplicated frame
    chain = resample_positions(positions, 4)
    assert chain.shape == (4, 2, 3)
    per_atom = np.linalg.norm(chain[-1][0] - chain[0][0])
    assert np.isclose(per_atom, np.sqrt(14.0))


def test_uniform_arc_length_distribution():
    # a path that lingers then rushes: images must come out uniformly spaced in s
    slow = np.linspace(0.0, 0.1, 6)[:, None, None] * np.ones((1, 1, 3))
    fast = np.linspace(0.2, 1.0, 6)[:, None, None] * np.ones((1, 1, 3))
    positions = np.concatenate([slow, fast], axis=0)
    chain = resample_positions(positions, 5)
    steps = np.linalg.norm(np.diff(chain[:, 0, :], axis=0), axis=1)
    assert np.ptp(steps) < 0.05 * steps.mean() + 1e-9
