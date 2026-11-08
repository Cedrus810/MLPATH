"""A refused initial target still carries the overlap just measured."""

import numpy as np
import pytest
from ase import Atoms

from prrs.calculators import double_well_factory
from prrs.config import SearchConfig
from prrs.runner import follow_min_mode


def test_initial_overlap_failure_keeps_measured_values(monkeypatch):
    import prrs.runner as runner

    # Split the bond tangent across two basis vectors; the 0.9 cutoff
    # deliberately refuses it at the first refresh, before any walk.
    monkeypatch.setattr(
        runner,
        "curvature_spectrum",
        lambda *a, **k: (np.ones(6), np.eye(6), {}, {}),
    )
    atoms = Atoms("C2", positions=[[0, 0, 0], [1.7, 0, 0]])
    saddle, info = follow_min_mode(
        atoms,
        double_well_factory,
        SearchConfig(closed_shell_only=False, min_mode_overlap=0.9),
        direction=np.array([[-1, 0, 0], [1, 0, 0]]),
    )
    assert saddle is None and info["reason"] == "lost_target_mode"
    assert info["steps"] == 0
    assert info["target_overlap"] == pytest.approx(1 / np.sqrt(2))
    assert info["target_overlap"] == info["mode_overlap"]
    assert info["selected_mode_index"] is not None
