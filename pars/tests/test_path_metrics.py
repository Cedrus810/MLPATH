from __future__ import annotations

import numpy as np
import pytest

from pars.config import ParsConfig
from pars.frontier.objectives import (
    energy_excursion,
    evaluate_candidate,
    mlp_uncertainty,
    path_curvature,
    recrossing,
)
from pars.schema import Provenance

from conftest import TRIO_EARLY, TRIO_LATE, make_candidate


def test_energy_excursion_is_max_minus_source():
    energies = np.asarray([-1.0, 0.5, 2.0, -0.5])
    assert energy_excursion(energies) == pytest.approx(3.0)


def test_recrossing_counts_returns_to_source_basin():
    clean = [TRIO_EARLY, TRIO_EARLY, TRIO_LATE, TRIO_LATE]
    one_return_then_leave = [TRIO_EARLY, TRIO_LATE, TRIO_EARLY, TRIO_LATE, TRIO_LATE]
    two_returns = [
        TRIO_EARLY,
        TRIO_LATE,
        TRIO_EARLY,
        TRIO_LATE,
        TRIO_EARLY,
        TRIO_LATE,
    ]
    assert recrossing(make_candidate("c1", "ch", [0] * 4, clean)) == 0.0
    assert recrossing(make_candidate("c2", "ch", [0] * 5, one_return_then_leave)) == 1.0
    assert recrossing(make_candidate("c3", "ch", [0] * 6, two_returns)) == 2.0


def test_recrossing_none_without_signatures():
    candidate = make_candidate("c1", "ch", [0.0, 0.1, 0.2], [TRIO_EARLY] * 3)
    candidate.signatures = None
    assert recrossing(candidate) is None


def test_uncertainty_uses_worst_frame_by_default():
    uncertainty = np.asarray([np.nan, 0.02, 0.5, 0.03])
    assert mlp_uncertainty(uncertainty, 1.0) == pytest.approx(0.5)
    assert mlp_uncertainty(None, 1.0) is None
    assert mlp_uncertainty(np.asarray([np.nan]), 1.0) is None


def test_curvature_distinguishes_straight_from_bent():
    t = np.linspace(0.0, 1.0, 11)
    straight = np.zeros((11, 2, 3))
    straight[:, 0, 0] = t
    straight[:, 1, 0] = t + 1.0
    bent = straight.copy()
    # bend one atom against the other: the INTERNAL geometry now curves
    bend = np.sin(np.pi * t) * 0.8
    bent[:, 0, 1] = bend
    bent[:, 1, 1] = -bend
    assert path_curvature(bent, 9) > path_curvature(straight, 9)
    assert path_curvature(straight, 9) == pytest.approx(0.0)


def test_curvature_invariant_under_rigid_motion():
    # the same path carried by a rotating/translating molecule must not count
    rng = np.random.default_rng(0)
    t = np.linspace(0.0, 1.0, 9)
    base = np.zeros((9, 3, 3))
    base[:, 0, 0] = t
    base[:, 1, 1] = 1.0 - t
    base[:, 2, 2] = 0.5
    angle = np.pi / 3.0
    rotation = np.asarray(
        [
            [np.cos(angle), -np.sin(angle), 0.0],
            [np.sin(angle), np.cos(angle), 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    moved = base @ rotation.T + rng.uniform(-5, 5, size=3)
    assert path_curvature(moved, 9) == pytest.approx(path_curvature(base, 9), rel=1e-6)


def _rotation(axis, angle):
    axis = np.asarray(axis, dtype=float) / np.linalg.norm(axis)
    cross = np.asarray(
        [[0.0, -axis[2], axis[1]], [axis[2], 0.0, -axis[0]], [-axis[1], axis[0], 0.0]]
    )
    return np.eye(3) + np.sin(angle) * cross + (1.0 - np.cos(angle)) * cross @ cross


def test_curvature_of_a_tumbling_rigid_molecule_is_zero():
    # a frozen molecule that only tumbles has no internal path at all; a
    # constant rotation cannot catch a transposed alignment, a varying one can
    rng = np.random.default_rng(1)
    molecule = rng.normal(size=(5, 3))
    frames = np.asarray(
        [molecule @ _rotation((1.0, 2.0, 0.5), angle).T for angle in np.linspace(0.0, 1.2, 15)]
    )
    assert path_curvature(frames, 9) == pytest.approx(0.0, abs=1e-9)


def test_curvature_invariant_under_frame_dependent_rotation():
    t = np.linspace(0.0, 1.0, 11)
    base = np.zeros((11, 3, 3))
    base[:, 0, 0] = t
    base[:, 1, 1] = 1.0 - t
    base[:, 2, 2] = 0.5 + 0.3 * np.sin(np.pi * t)
    moved = np.asarray(
        [
            frame @ _rotation((0.3, -1.0, 2.0), 0.15 * index).T + 0.2 * index
            for index, frame in enumerate(base)
        ]
    )
    assert path_curvature(moved, 9) == pytest.approx(path_curvature(base, 9), rel=1e-6)


def test_evaluate_candidate_reports_unavailable_not_fabricated():
    candidate = make_candidate("c1", "ch", [0.0, 1.0, 0.2], [TRIO_EARLY, TRIO_LATE, TRIO_LATE])
    raw = evaluate_candidate(candidate, ParsConfig())
    assert raw["energy_excursion"] == pytest.approx(1.0)
    assert raw["trigger_severity[stretch/Angstrom]"] == 0.5
    assert raw["mlp_uncertainty"] is None  # not measured -> None, never a guess
    assert raw["endpoint_penalty"] is None  # design 8.6: PRRS stays the authority
    assert raw["path_curvature"] >= 0.0
    assert raw["recrossing"] == 0.0


def test_provenanceless_objectives_do_not_crash():
    candidate = make_candidate("c1", "ch", [0.0, 1.0, 0.2], [TRIO_EARLY, TRIO_LATE, TRIO_LATE])
    candidate.provenance = Provenance(
        run_dir="s",
        attempt_id="t",
        channel_id="ch",
        event_key=None,
        reactant_state_id=None,
        product_state_id=None,
        frame_start=0,
        frame_stop=3,
    )
    raw = evaluate_candidate(candidate, ParsConfig())
    assert raw["trigger_severity"] is None


def test_trigger_severity_compares_only_within_one_family():
    """Audit ALG-01: 0.5 Angstrom of stretch and 0.2 eV of kick share no
    scale, so the bare numbers must not decide dominance across families."""
    from dataclasses import replace

    from pars.frontier.dominance import dominates

    path = [0.0, 1.0, 0.2]
    signatures = [TRIO_EARLY, TRIO_LATE, TRIO_LATE]
    stretch = make_candidate("c1", "ch", path, signatures, amplitude=0.5)
    kick = make_candidate("c2", "ch", path, signatures, amplitude=0.2)
    kick.provenance = replace(
        kick.provenance, perturbation_family="kick", perturbation_unit="eV"
    )
    small_stretch = make_candidate("c3", "ch", path, signatures, amplitude=0.1)
    raw = {
        c.candidate_id: evaluate_candidate(c, ParsConfig())
        for c in (stretch, kick, small_stretch)
    }
    assert raw["c2"]["trigger_severity[kick/eV]"] == 0.2
    assert not dominates(raw["c2"], raw["c1"])
    assert not dominates(raw["c1"], raw["c2"])
    assert dominates(raw["c3"], raw["c1"])  # same family still ranks


def test_trigger_severity_without_family_or_unit_is_unavailable():
    from dataclasses import replace

    candidate = make_candidate("c1", "ch", [0.0, 1.0, 0.2], [TRIO_EARLY, TRIO_LATE, TRIO_LATE])
    candidate.provenance = replace(candidate.provenance, perturbation_unit=None)
    raw = evaluate_candidate(candidate, ParsConfig())
    assert raw["trigger_severity"] is None


def test_curvature_is_integral_of_squared_curvature_over_arc_length():
    # integral |d2R/ds2|^2 ds has units 1/length: scaling the whole path (and
    # molecule) by 2 halves it. The old /h normalization doubled it instead.
    t = np.linspace(0.0, 1.0, 11)
    path = np.zeros((11, 3, 3))
    path[:, 0, 0] = t
    path[:, 1, 1] = 1.0 - t
    path[:, 2, 2] = 0.5 + 0.3 * np.sin(np.pi * t)
    assert path_curvature(2.0 * path, 9) == pytest.approx(
        path_curvature(path, 9) / 2.0, rel=1e-9
    )


def test_curvature_objective_skips_the_forcing_gap():
    # kept_frames [0, 3, 4, ...]: frame 0 is the source, frames 1-2 of the
    # original trajectory (the forcing) are gone. The jump across that gap is
    # an artefact of the perturbation, not path shape: the free segment here is
    # straight, so the objective must be zero.
    import dataclasses

    # atoms 0/1 stretch symmetrically along x around atom 2 on the y axis, so the
    # optimal alignment is exactly the identity and the free segment is a
    # straight line; the source frame alone has atom 2 further out
    frames = 8
    half = np.concatenate([[1.0], np.linspace(1.0, 1.7, frames - 1)])
    positions = np.zeros((frames, 3, 3))
    positions[:, 0, 0] = -half
    positions[:, 1, 0] = half
    positions[:, 2, 1] = 1.0
    positions[0, 2, 1] = 1.6  # source geometry off the free segment's line
    candidate = make_candidate(
        "c-gap", "ch", np.zeros(frames), [TRIO_EARLY] * frames, positions=positions
    )
    candidate.provenance = dataclasses.replace(
        candidate.provenance, extra={"kept_frames": [0] + list(range(3, 3 + frames - 1))}
    )
    assert path_curvature(positions, 9) > 0.0  # the kink is real in the raw chain
    raw = evaluate_candidate(candidate, ParsConfig())
    assert raw["path_curvature"] == pytest.approx(0.0, abs=1e-12)


def _vibrating_bend(frames: int) -> np.ndarray:
    # one smooth bend with six 0.1 A vibrations riding on it, sampled at any
    # number of frames: the reaction path is the same, only the MD sampling is not
    t = np.linspace(0.0, 1.0, frames)
    positions = np.zeros((frames, 3, 3))
    positions[:, 0, 0] = -(1.0 + 0.6 * t)
    positions[:, 1, 0] = 1.0 + 0.6 * t
    positions[:, 2, 1] = 1.0 + 0.4 * np.sin(np.pi * t)
    positions[:, 2, 2] = 0.1 * np.sin(2 * np.pi * 6 * t)
    return positions


def test_curvature_objective_does_not_depend_on_frame_count():
    # measured on raw frames, curvature grew with the number of frames (rho
    # 0.86-0.92 on the p1 runs) and duplicated recrossing; at NEB image
    # resolution the same path must score the same however densely it was sampled
    def curvature(frames):
        candidate = make_candidate(
            f"c-{frames}",
            "ch",
            np.zeros(frames),
            [TRIO_EARLY] * frames,
            positions=_vibrating_bend(frames),
        )
        return evaluate_candidate(candidate, ParsConfig())["path_curvature"]

    assert curvature(12) == pytest.approx(curvature(40), rel=0.05)
