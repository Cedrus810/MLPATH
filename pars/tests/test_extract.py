from __future__ import annotations

import numpy as np
import pytest

from pars.config import ParsConfig
from pars.paths.candidate import build_candidates
from pars.paths.ingest import ingest_run


def test_window_keeps_free_response_and_drops_forcing_and_quench(run_dir):
    traces = ingest_run(run_dir, ParsConfig())
    candidates = build_candidates(traces, ParsConfig())
    by_id = {c.candidate_id: c for c in candidates}

    candidate = by_id["pc-rc0000-t000001"]
    # source(0) + free_start(2) + free(3,4,5); perturbed(1), quench(6,7) dropped
    assert candidate.positions.shape[0] == 5
    assert candidate.energies.shape[0] == 5
    assert len(candidate.signatures) == 5
    assert candidate.provenance.frame_start == 0
    assert candidate.provenance.frame_stop == 6


def test_forcing_parameters_survive_in_provenance(run_dir):
    traces = ingest_run(run_dir, ParsConfig())
    candidate = build_candidates(traces, ParsConfig())[0]
    assert candidate.provenance.perturbation_amplitude == 0.6
    assert candidate.provenance.perturbation_family == "stretch"


def test_candidate_id_is_deterministic_per_channel_and_trial(run_dir):
    traces = ingest_run(run_dir, ParsConfig())
    again = build_candidates(ingest_run(run_dir, ParsConfig()), ParsConfig())
    assert [c.candidate_id for c in build_candidates(traces, ParsConfig())] == [
        c.candidate_id for c in again
    ]


def test_too_few_free_frames_fail_closed(run_dir):
    config = ParsConfig(min_path_frames=7)
    traces = ingest_run(run_dir, ParsConfig())
    with pytest.raises(ValueError, match="frames match"):
        build_candidates(traces, config)


def test_positions_are_a_copy_not_a_view(run_dir):
    traces = ingest_run(run_dir, ParsConfig())
    candidate = build_candidates(traces, ParsConfig())[0]
    before = candidate.positions.copy()
    candidate.positions[:] = 0.0
    traces2 = ingest_run(run_dir, ParsConfig())
    assert np.array_equal(build_candidates(traces2, ParsConfig())[0].positions, before)
