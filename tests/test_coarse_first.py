"""A source's coarse coverage must not compete with its amplitude refinement."""

import json
from dataclasses import replace

import pytest

from prrs.calculators import demo_atoms, double_well_factory
from prrs.config import SearchConfig
from prrs.perturbations import propose
from prrs.search import ReactionSearch


def _protocol():
    config = SearchConfig(
        families=("stretch", "kick"),
        geometry_amplitudes_A=(0.1, 0.3, 0.7, 1.0),
        kick_energies_eV=(0.1, 0.5, 1.0),
        direction_mode="exhaustive",
        timestep_fs=0.1,
        response_steps=20,
        quench_fmax_eV_A=0.001,
        closed_shell_only=False,
        max_depth=1,
        max_trials=200,
        chemical_trials_per_node=200,
        conformer_trials_per_microstate=200,
        refinement_steps=2,
        saddle_search_enabled=False,
        boundary_continuation_enabled=False,
        irc_enabled=False,
    )
    coarse_count = sum(len(batch) for batch in propose(demo_atoms(), config, config.seed))
    return config, coarse_count


def _attempts(output):
    return [json.loads(line) for line in (output / "attempts.jsonl").read_text().splitlines()]


@pytest.mark.parametrize(
    "budget", ["max_trials", "chemical_trials_per_node", "conformer_trials_per_microstate"]
)
def test_exact_coarse_budget_covers_every_direction_before_refinement(tmp_path, budget):
    config, needed = _protocol()
    config = replace(config, **{budget: needed})
    network = ReactionSearch(double_well_factory, config).run(demo_atoms(), tmp_path)
    records = _attempts(tmp_path)
    assert len(records) == needed
    assert {record["stage"] for record in records} == {"coarse"}
    assert network["direction_ledger"]["directions_not_executed"] == 0
    assert network["direction_ledger"]["coverage_censored"] is False
    assert all(
        scan["proposal"]["schedule"]["amplitudes_not_executed"] == 0
        for scan in network["amplitude_scans"]
    )
    # This fixture really offers a return/departure window. A missing opportunity to
    # refine would make the stage assertion pass without exercising the competition.
    assert any(scan["local_brackets"] for scan in network["amplitude_scans"])


def test_spare_budget_refines_after_all_coarse_and_keeps_bracket_provenance(tmp_path):
    config, needed = _protocol()
    config = replace(config, max_trials=needed + 2)
    network = ReactionSearch(double_well_factory, config).run(demo_atoms(), tmp_path)
    records = _attempts(tmp_path)
    assert [record["stage"] for record in records[:needed]] == ["coarse"] * needed
    assert [record["stage"] for record in records[needed:]] == ["refine"] * 2
    by_id = {record["attempt_id"]: record for record in records}
    for scan in network["amplitude_scans"]:
        for bracket in scan["local_brackets"]:
            high = next(
                sample for sample in scan["samples"] if sample["amplitude"] == bracket["high"]
            )
            assert bracket["target_microstate"] == by_id[high["attempt"]]["target_microstate"]
        for sample in scan["samples"]:
            record = by_id[sample["attempt"]]
            assert record["source_microstate"] == scan["source_microstate"]
            if record["stage"] == "refine":
                assert all(parent in by_id for parent in record["refinement_parents"])


@pytest.mark.parametrize("kill_after", [5, 10, 11], ids=["coarse", "between_passes", "refine"])
def test_resume_in_both_passes_is_byte_identical(tmp_path, kill_after):
    import subprocess
    import sys
    from pathlib import Path

    config, needed = _protocol()
    config = replace(config, max_trials=needed + 2)
    full, killed = tmp_path / "full", tmp_path / "killed"
    ReactionSearch(double_well_factory, config).run(demo_atoms(), full)
    root = str(Path(__file__).resolve().parents[1])
    child = tmp_path / "child.py"
    child.write_text(
        f"""
import os, sys
sys.path.insert(0, {root!r} + '/src')
import prrs.search as search
from prrs.calculators import demo_atoms, double_well_factory
from prrs.config import SearchConfig
real = search.run_trial
count = 0
def counted(*args, **kwargs):
    global count
    count += 1
    if count > {kill_after}:
        os._exit(137)
    return real(*args, **kwargs)
search.run_trial = counted
search.ReactionSearch(double_well_factory, SearchConfig(**{config.to_dict()!r})).run(demo_atoms(), sys.argv[1])
"""
    )
    result = subprocess.run([sys.executable, str(child), str(killed)], capture_output=True)
    assert result.returncode == 137, result.stderr.decode()
    mid = _attempts(killed)
    assert len(mid) == kill_after
    assert mid[-1]["stage"] == ("refine" if kill_after > needed else "coarse")
    ReactionSearch(double_well_factory, config).run(demo_atoms(), killed, resume=True)
    for name in ["network.json", "attempts.jsonl"]:
        assert (killed / name).read_bytes() == (full / name).read_bytes(), name
