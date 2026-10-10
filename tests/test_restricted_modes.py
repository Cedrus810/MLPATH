"""Soft modes belong to the restricted Hessian, not normalized projection scraps."""

import numpy as np
import pytest

from prrs.runner import _restricted_soft_modes


def test_an_almost_torsional_eigenvector_does_not_supply_a_fake_soft_mode():
    eps = 1e-6
    vectors = np.array(
        [[np.sqrt(1 - eps**2), -eps, 0], [eps, np.sqrt(1 - eps**2), 0], [0, 0, 1]]
    )
    values = np.array([0.01, 10.0, 1.0])
    modes = _restricted_soft_modes(values, vectors, np.array([[1.0, 0, 0]]), 2)
    # Normalizing the projected first eigenvector would choose y and falsely keep
    # its 0.01 eigenvalue. In the restricted Hessian z is softest and y costs ~10.
    assert modes[0][0] == pytest.approx(1.0)
    assert abs(float(modes[0][1] @ np.array([0, 0, 1]))) == pytest.approx(1.0)
    assert modes[1][0] == pytest.approx(10.0, abs=1e-9)


def test_dependent_exclusions_remove_only_their_actual_span():
    values = np.array([0.01, 1.0, 2.0])
    excluded = np.array([[1.0, 0, 0], [2.0, 0, 0]])
    modes = _restricted_soft_modes(values, np.eye(3), excluded, 4)
    assert len(modes) == 2
    assert [value for value, _ in modes] == pytest.approx([1.0, 2.0])
    assert np.array([mode for _, mode in modes]) @ np.array(
        [mode for _, mode in modes]
    ).T == pytest.approx(np.eye(2))
    assert np.max(np.abs(excluded @ np.array([mode for _, mode in modes]).T)) < 1e-12


def test_restricted_modes_are_eigenvectors_of_the_projected_operator():
    rng = np.random.default_rng(17)
    vectors, _ = np.linalg.qr(rng.normal(size=(6, 6)))
    values = np.arange(1.0, 7.0)
    excluded = rng.normal(size=(2, 6))
    modes = _restricted_soft_modes(values, vectors, excluded, 3)
    q, _ = np.linalg.qr(excluded.T)
    projection = np.eye(6) - q @ q.T
    hessian = (vectors * values) @ vectors.T
    assert len(modes) == 3
    for value, mode in modes:
        assert projection @ hessian @ projection @ mode == pytest.approx(
            value * mode, abs=1e-12
        )
        assert np.linalg.norm(mode) == pytest.approx(1.0)


def test_no_internal_dimension_left_returns_no_modes():
    assert _restricted_soft_modes(np.ones(3), np.eye(3), np.eye(3), 4) == ()
