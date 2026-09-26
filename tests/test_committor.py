"""The committor shot machinery (PLAN item 2, decision 3 of 2026-09-20).

The protocol is frozen before any run: docs/P1_COMMITTOR_PROTOCOL.md. What is tested
here is the machinery, not the chemistry -- a T=0 config must be refused (fail
closed), shots must differ by seed and be reproducible, and the Wilson interval must
reproduce the number the project already quotes (16/16 -> lower bound 0.8064).
"""

import json

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


def test_a_shot_starts_from_the_frame_it_is_given(tmp_path):
    # 2026-09-23: shots re-delivered the probe on top of the ridge frame (+0.4 eV on P1,
    # barrier 419 meV). A shot's first perturbed frame must sit at the source's
    # potential energy: thermal momenta only, no displacement, kick or pulse.
    import json

    source = relax_source(demo_atoms(), double_well_factory, HOT)
    for probe in (
        _probe(),
        Probe(family="kick", indices=(0, 1), sign=1, amplitude=0.5, seed=7),
    ):
        (outcome,) = shoot(
            source, probe, double_well_factory, HOT, tmp_path / probe.family, 1, 3
        )
        assert outcome.record["probe_delivery"] == "skipped:committor_shot"
        assert "probe_delivered" not in outcome.record
        rows = [
            json.loads(line)
            for line in (tmp_path / probe.family / "trials" / "shot0000" / "observations.jsonl")
            .read_text()
            .splitlines()
        ]
        assert rows[1]["phase"] == "perturbed"
        assert rows[1]["potential_eV"] == pytest.approx(rows[0]["potential_eV"], abs=1e-12)


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


# (0.1, 0.9) for the first-entry case: with the fixed-rule amplitudes, refinement ends
# the only window on a transient crossing (0.600), which names no B to enter; (0.1, 0.9)
# still refines twice and ends on a confirmed product (measured 2026-09-26).
@pytest.mark.parametrize(
    "entry_steps, amplitudes", [(0, (0.1, 0.3, 0.7, 1.0)), (10, (0.1, 0.9))]
)
def test_brackets_carry_a_committor_record_when_shots_are_enabled(
    tmp_path, entry_steps, amplitudes
):
    """End to end on the analytic double well: the record's shape and its arithmetic.

    Each bracket gets {shots, reached_A, reached_B, other, p_B_given_AB, wilson_95, ...}; the
    three counters add up to the shot count; other is never folded into p_B's
    denominator. The committor lands on the bracket (model report only), and with
    shots disabled nothing appears at all.
    """
    from prrs.search import ReactionSearch

    config = HOT.__class__(
        **{
            **HOT.to_dict(),
            "geometry_amplitudes_A": amplitudes,
            "response_steps": 20,
            "max_trials": 40,
            "refinement_steps": 2,
            "committor_shots": 3,
            "committor_entry_steps": entry_steps,
            "committor_max_steps": 200,
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
            assert record["p_B_given_AB"] == pytest.approx(record["reached_B"] / terminated)
        else:
            assert record["p_B_given_AB"] is None
        # bonded vs dissociated ends have different graphs, so entry applies wherever
        # the bracket names its B; a window whose high end never confirmed a product
        # has no B to enter and keeps the fixed rule, saying why
        if not entry_steps:
            expected = "fixed_t_free"
        elif bracket.get("target_microstate") is not None:
            expected = "first_entry"
        else:
            expected = "fixed_t_free:no_target_microstate"
        assert record["termination"] == expected
        if expected == "first_entry":
            assert record["entry_steps"] == entry_steps
            assert 0 <= record["entry_quench_mismatch"] <= terminated
        else:
            assert "entry_quench_mismatch" not in record
        assert record["p_B_conditioned_on"] == "reached_A or reached_B"
        assert record["other_fraction"] == pytest.approx(record["other"] / 3)
        assert record["max_nve_drift_eV_atom"] >= 0.0
    if entry_steps:
        used = [b for b in brackets if b["committor"]["termination"] == "first_entry"]
        assert used, "no bracket named its B; the first-entry path went unexercised"


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


def _barrier_frame():
    """C2 on the double-well barrier top (r = 1.8), where the bond graph is undecided."""
    atoms = demo_atoms()
    atoms.positions = [[-0.9, 0, 0], [0.9, 0, 0]]
    atoms.calc = double_well_factory()
    atoms.get_potential_energy()
    return atoms


def _end_graphs(config):
    from prrs.state import encode

    graphs = {}
    for side, r in (("A", 1.2), ("B", 2.4)):
        atoms = demo_atoms()
        atoms.positions = [[-r / 2, 0, 0], [r / 2, 0, 0]]
        graphs[side] = frozenset(encode(atoms, config.bond_scale).edges)
    assert graphs["A"] != graphs["B"]
    return graphs


def _offline_entry(directory, graphs, k):
    """The frozen rule applied afterwards to a full-length run, every step saved."""
    rows = [json.loads(line) for line in (directory / "observations.jsonl").open()]
    free = [r for r in rows if r["phase"] == "free"]
    assert [r["step"] for r in free] == list(range(1, len(free) + 1)), "every step saved"
    run, last = 0, None
    for row in free:
        edges = frozenset(tuple(e) for e in row["edges"])
        label = next((side for side, g in graphs.items() if edges == g), None)
        run = run + 1 if label is not None and label == last else int(label is not None)
        last = label
        if run >= k:
            return label, row["step"]
    return None, None


def _free_steps(directory):
    rows = [json.loads(line) for line in (directory / "observations.jsonl").open()]
    return max(r["step"] for r in rows if r["phase"] == "free")


@pytest.mark.parametrize("k", [5, 20])
def test_stopping_on_first_entry_matches_the_rule_applied_afterwards(tmp_path, k):
    """Protocol stage 2's requirement, on the analytic well: same seed, same verdict.

    One run stops at entry; the other runs the full length with every step saved and
    the rule is applied afterwards. Side and entry step must agree, and the early
    run's free segment must end exactly at that step.
    """
    config = HOT.__class__(**{**HOT.to_dict(), "response_steps": 400, "sample_interval": 1})
    graphs = _end_graphs(config)
    frame = _barrier_frame()
    entry = {"graphs": graphs, "k": k}
    (stopped,) = shoot(
        frame,
        _probe(),
        double_well_factory,
        config,
        tmp_path / "stop",
        1,
        31,
        first_entry=entry,
    )
    (full,) = shoot(frame, _probe(), double_well_factory, config, tmp_path / "full", 1, 31)
    side, step = _offline_entry(tmp_path / "full" / "trials" / "shot0000", graphs, k)
    assert side is not None, "the fixture must enter; otherwise this compares nothing"
    assert stopped.record["first_entry"]["side"] == side
    assert stopped.record["first_entry"]["step"] == step
    assert _free_steps(tmp_path / "stop" / "trials" / "shot0000") == step
    assert "first_entry" not in full.record, "off by default: the record is unchanged"


def test_a_shot_that_never_enters_runs_to_the_cap_uncommitted(tmp_path):
    config = HOT.__class__(**{**HOT.to_dict(), "response_steps": 60})
    entry = {"graphs": _end_graphs(config), "k": 10**6}
    (outcome,) = shoot(
        _barrier_frame(),
        _probe(),
        double_well_factory,
        config,
        tmp_path / "s",
        1,
        31,
        first_entry=entry,
    )
    assert outcome.record["first_entry"]["side"] is None
    assert outcome.record["first_entry"]["step"] is None
    assert _free_steps(tmp_path / "s" / "trials" / "shot0000") == 60
