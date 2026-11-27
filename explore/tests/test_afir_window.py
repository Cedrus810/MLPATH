"""afir_window.bond_windows: one window [j-1, j+2] per bond-graph change (2026-11-25)."""

import numpy as np

from explore.afir_window import bond_windows


def test_one_window_at_the_bond_change():
    # H2 pulled apart: bonded (0.74 A) for the first nodes, broken above 1.2 * 0.62 = 0.744 A
    numbers = np.array([1, 1])
    ds = [0.70, 0.72, 0.74, 0.76, 1.00, 1.20, 1.40]
    x = np.array([[[0.0, 0.0, 0.0], [d, 0.0, 0.0]] for d in ds])
    assert bond_windows(x, numbers, 1.2) == [(1, 4, 2)]  # cut 1.2 * 0.62 = 0.744 A


def test_no_change_no_window_and_ends_clipped():
    numbers = np.array([1, 1])
    x = np.array([[[0.0, 0.0, 0.0], [d, 0.0, 0.0]] for d in (0.7, 0.71, 0.72)])
    assert bond_windows(x, numbers, 1.2) == []
    y = np.array([[[0.0, 0.0, 0.0], [d, 0.0, 0.0]] for d in (0.7, 1.4, 1.5)])
    assert bond_windows(y, numbers, 1.2) == [(0, 2, 0)]
