"""Resume acceptance (PLAN item 6): killed mid-run, resumed, byte-identical.

The only accepted form of proof (PLAN 6.3's own words): run to half, kill, resume,
and the result must equal an uninterrupted run's network.json and attempts.jsonl
byte for byte. The kill is a real one -- a child process calls os._exit(137) from
inside run_trial, so no exception handler runs and the network stays at status
"running", which is exactly the on-disk situation a power loss or SIGKILL leaves.
"""

import json
import subprocess
import sys
from pathlib import Path

import pytest

from prrs.calculators import demo_atoms, double_well_factory
from prrs.config import SearchConfig
from prrs.search import ReactionSearch

ROOT = Path(__file__).resolve().parents[1]


def _config():
    return SearchConfig(
        families=("stretch", "compress"),
        geometry_amplitudes_A=(0.1, 0.3, 0.7, 1.0),
        timestep_fs=0.1,
        response_steps=20,
        quench_fmax_eV_A=0.001,
        closed_shell_only=False,
        max_depth=2,
        max_trials=20,
        refinement_steps=2,
    )


CHILD = """
import os
import sys

sys.path.insert(0, {root!r})
sys.path.insert(0, {root!r} + "/src")

import prrs.search as search
from prrs.runner import run_trial as real_run_trial

limit = int(sys.argv[2])
seen = {{"n": 0}}


def counted(*args, **kwargs):
    seen["n"] += 1
    if seen["n"] > limit:
        os._exit(137)
    return real_run_trial(*args, **kwargs)


search.run_trial = counted

from prrs.calculators import demo_atoms, double_well_factory
from prrs.config import SearchConfig
from prrs.search import ReactionSearch

config = SearchConfig(
    families=("stretch", "compress"),
    geometry_amplitudes_A=(0.1, 0.3, 0.7, 1.0),
    timestep_fs=0.1,
    response_steps=20,
    quench_fmax_eV_A=0.001,
    closed_shell_only=False,
    max_depth=2,
    max_trials=20,
    refinement_steps=2,
)
ReactionSearch(double_well_factory, config).run(demo_atoms(), sys.argv[1])
"""


def _kill_after(tmp_path, trial_limit):
    script = tmp_path / "child.py"
    script.write_text(CHILD.format(root=str(ROOT)), encoding="utf-8")
    output = tmp_path / "killed"
    subprocess.run(
        [
            sys.executable,
            str(script),
            str(output),
            str(trial_limit),
        ],
        check=False,
        capture_output=True,
    )
    assert (output / "network.json").exists()
    return output


def test_killed_then_resumed_is_byte_identical_to_uninterrupted(tmp_path):
    full = ReactionSearch(double_well_factory, _config()).run(demo_atoms(), tmp_path / "full")
    assert full["status"] == "completed"

    killed = _kill_after(tmp_path, trial_limit=6)
    mid = json.loads((killed / "network.json").read_text(encoding="utf-8-sig"))
    assert mid["status"] == "running"
    assert mid.get("queue"), "the kill must happen with the queue still populated"

    resumed = ReactionSearch(double_well_factory, _config()).run(
        demo_atoms(), killed, resume=True
    )
    assert resumed["status"] == "completed"

    assert (killed / "network.json").read_bytes() == (
        tmp_path / "full" / "network.json"
    ).read_bytes()
    assert (killed / "attempts.jsonl").read_bytes() == (
        tmp_path / "full" / "attempts.jsonl"
    ).read_bytes()
    # the rollback is named, and lives OUTSIDE network.json so the comparison above
    # is even meaningful
    report = json.loads((killed / "resume.json").read_text(encoding="utf-8-sig"))
    assert "deleted_trial_dirs" in report and "dropped_attempt_rows" in report


def test_resume_refuses_a_finished_or_mismatched_run(tmp_path):
    finished = tmp_path / "finished"
    ReactionSearch(double_well_factory, _config()).run(demo_atoms(), finished)
    with pytest.raises(ValueError, match="running"):
        ReactionSearch(double_well_factory, _config()).run(demo_atoms(), finished, resume=True)

    killed = _kill_after(tmp_path, trial_limit=4)
    other = _config().__class__(**{**_config().to_dict(), "max_trials": 21})
    with pytest.raises(ValueError, match="config"):
        ReactionSearch(double_well_factory, other).run(demo_atoms(), killed, resume=True)

    fresh = tmp_path / "nowhere"
    with pytest.raises(FileNotFoundError):
        ReactionSearch(double_well_factory, _config()).run(demo_atoms(), fresh, resume=True)
