"""The committor shot machinery (PLAN item 2, decision 3 of 2026-09-20).

The protocol is frozen before any run: docs/P1_COMMITTOR_PROTOCOL.md. What is tested
here is the machinery, not the chemistry -- a T=0 config must be refused (fail
closed), shots must differ by seed and be reproducible, and the Wilson interval must
reproduce the number the project already quotes (16/16 -> lower bound 0.8064).
"""

import numpy as np
import pytest

from prrs.calculators import demo_atoms, double_well_factory
from prrs.config import SearchConfig
from prrs.perturbations import Probe
from prrs.runner import relax_source, shoot, wilson_interval


HOT = SearchConfig(
    families=("stretch", "compress"),
    geometry_amplitudes_A=(0.5,),
    timestep_fs=0.1,
    response_steps=25,
    quench_fmax_eV_A=1e-3,
    quench_steps=400,
    closed_shell_only=False,
    max_depth=1,
    temperature_K=300.0,
    committor_shots=4,
)


def _probe():
    return Probe(family="stretch", indices=(0, 1), sign=1, amplitude=0.5, seed=7)


def test_zero_temperature_is_refused_not_shot(tmp_path):
    frozen = HOT.__class__(**{**HOT.to_dict(), "temperature_K": 0.0})
    source = relax_source(demo_atoms(), double_well_factory, frozen)
    with pytest.raises(ValueError, match="temperature_K"):
        shoot(source, _probe(), double_well_factory, frozen, tmp_path, shots=2, seed0=1)


def test_shots_differ_by_seed_and_reproduce(tmp_path):
    source = relax_source(demo_atoms(), double_well_factory, HOT)
    first = shoot(source, _probe(), double_well_factory, HOT, tmp_path / "a", 4, 101)
    assert [o.record["probe"]["seed"] for o in first] == [101, 102, 103, 104]
    assert len({o.record["attempt_id"] for o in first}) == 4
    # reproducible: the same seeds give the same endpoint energies, other seeds do not
    again = shoot(source, _probe(), double_well_factory, HOT, tmp_path / "b", 4, 101)
    for one, two in zip(first, again):
        assert one.record["status"] == two.record["status"]
        if one.record["status"] == "completed":
            assert one.endpoint.get_positions() == pytest.approx(two.endpoint.get_positions())
    other = shoot(source, _probe(), double_well_factory, HOT, tmp_path / "c", 4, 900)
    seeds_matter = any(
        o1.record["status"] != o2.record["status"]
        or (
            o1.record["status"] == "completed"
            and not np.allclose(o1.endpoint.positions, o2.endpoint.positions)
        )
        for o1, o2 in zip(first, other)
    )
    assert seeds_matter, "T=300 shots with different seeds produced identical trajectories"


def test_wilson_matches_the_number_the_project_already_quotes():
    lo, hi = wilson_interval(16, 16)
    assert lo == pytest.approx(0.8064, abs=5e-5)
    assert hi == pytest.approx(1.0, abs=1e-9)
    lo, hi = wilson_interval(8, 16)
    # n=16 at p=0.5 spans roughly [0.28, 0.72]: that width is WHY the protocol says
    # N=16 only answers "near 0 or 1", never 0.5 vs 0.25
    assert 0.25 < lo < 0.35
    assert 0.65 < hi < 0.78
    with pytest.raises(ValueError):
        wilson_interval(2, 0)
    with pytest.raises(ValueError):
        wilson_interval(2, 1)


def test_committor_shots_defaults_off_and_validates():
    assert SearchConfig().committor_shots == 0
    assert HOT.committor_shots == 4
    with pytest.raises(ValueError):
        SearchConfig(committor_shots=-1)


def test_brackets_carry_a_committor_record_when_shots_are_enabled(tmp_path):
    """End to end on the analytic double well: the record's shape and its arithmetic.

    Each bracket gets {shots, reached_A, reached_B, other, p_B, wilson_95, ...}; the
    three counters add up to the shot count; other is never folded into p_B's
    denominator. The committor lands on the bracket (model report only), and with
    shots disabled nothing appears at all.
    """
    from prrs.search import ReactionSearch

    config = HOT.__class__(
        **{
            **HOT.to_dict(),
            "geometry_amplitudes_A": (0.1, 0.3, 0.7, 1.0),
            "response_steps": 20,
            "max_trials": 40,
            "refinement_steps": 2,
            "committor_shots": 3,
        }
    )
    network = ReactionSearch(double_well_factory, config).run(demo_atoms(), tmp_path / "run")
    brackets = [
        bracket for scan in network["amplitude_scans"] for bracket in scan["local_brackets"]
    ]
    assert brackets, "this potential brackets; the wiring needs one to exercise"
    for bracket in brackets:
        record = bracket["committor"]
        assert record["shots"] == 3
        assert record["reached_A"] + record["reached_B"] + record["other"] == 3
        assert record["ensemble"] == "maxwell-boltzmann"
        assert record["temperature_K"] == config.temperature_K
        terminated = record["reached_A"] + record["reached_B"]
        if terminated:
            assert record["p_B"] == pytest.approx(record["reached_B"] / terminated)
        else:
            assert record["p_B"] is None


def test_a_shot_that_terminates_is_judged_a_basin_not_other(tmp_path):
    """The hole the first production run fell into: A and B unreachable by construction.

    2026-09-22, P1 on kasuga02: every bracket reported A 0 / B 0 / other 16, and the
    48 endpoints were all quenched, converged, and sitting in a basin. Judging read
    `record["target_microstate"]` -- a key written only by the SCAN path, after
    registry.admit -- so no shot could ever be anything but "other". The shape test
    above stayed green throughout: its counters added up, and its p_B branch is
    guarded by `if terminated`, which was never true.

    So this asserts the one thing that test cannot: that a shot which terminates is
    attributed to a basin. Delete the registry.identify call in search.py and this
    goes red while every other committor test stays green.
    """
    from prrs.search import ReactionSearch

    config = HOT.__class__(
        **{
            **HOT.to_dict(),
            "geometry_amplitudes_A": (0.1, 0.3, 0.7, 1.0),
            "response_steps": 20,
            "max_trials": 40,
            "refinement_steps": 2,
            "committor_shots": 3,
        }
    )
    network = ReactionSearch(double_well_factory, config).run(demo_atoms(), tmp_path / "run")
    records = [
        bracket["committor"]
        for scan in network["amplitude_scans"]
        for bracket in scan["local_brackets"]
    ]
    assert records, "this potential brackets; the wiring needs one to exercise"
    for record in records:
        assert sum(record["landings"].values()) == record["shots"]
        named = {
            where: count
            for where, count in record["landings"].items()
            if not where.startswith("unidentified:")
        }
        assert named, (
            "not one shot could be attributed to a microstate -- the signature of "
            "judging reading a key nothing sets, not of a potential without basins; "
            f"landings were {record['landings']}"
        )
