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
climb_limit = int(sys.argv[3]) if len(sys.argv) > 3 else None
seen = {{"n": 0, "climbs": 0}}


def counted(*args, **kwargs):
    seen["n"] += 1
    if seen["n"] > limit:
        os._exit(137)
    return real_run_trial(*args, **kwargs)


search.run_trial = counted

if climb_limit is not None:
    from prrs.search import follow_min_mode as real_follow_min_mode

    def counted_climb(*args, **kwargs):
        # killed INSIDE a climb: its paths/ directory exists, its scan is unpublished
        seen["climbs"] += 1
        if seen["climbs"] > climb_limit:
            os._exit(137)
        return real_follow_min_mode(*args, **kwargs)

    search.follow_min_mode = counted_climb

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


def _kill_after(tmp_path, trial_limit, climb_limit=None, name="killed"):
    script = tmp_path / "child.py"
    script.write_text(CHILD.format(root=str(ROOT)), encoding="utf-8")
    output = tmp_path / name
    subprocess.run(
        [
            sys.executable,
            str(script),
            str(output),
            str(trial_limit),
            *([] if climb_limit is None else [str(climb_limit)]),
        ],
        check=False,
        capture_output=True,
    )
    assert (output / "network.json").exists()
    return output


# 6 kills between microstates -- the only point the original acceptance exercised.
# 4 kills inside the first scan after its reaction was published, which is where the
# old reference-based rollback kept rows around holes and collided on t000002.
@pytest.mark.parametrize("trial_limit", [4, 6])
def test_killed_then_resumed_is_byte_identical_to_uninterrupted(tmp_path, trial_limit):
    full = ReactionSearch(double_well_factory, _config()).run(demo_atoms(), tmp_path / "full")
    assert full["status"] == "completed"

    killed = _kill_after(tmp_path, trial_limit=trial_limit)
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
    assert _tree(killed) == _tree(tmp_path / "full")


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


def _tree(root):
    """What the rollback is responsible for: the published record and the two
    directories written before their scan is published."""
    files = [root / "network.json", root / "attempts.jsonl"]
    for sub in ("paths", "ts_candidates"):
        files += [path for path in sorted((root / sub).rglob("*")) if path.is_file()]
    return {str(path.relative_to(root)): path.read_bytes() for path in files}


def test_a_kill_inside_a_climb_does_not_spend_the_saddle_budget(tmp_path):
    """The climb's directory is written before its scan is published. Counting it on
    resume spent budget on rolled-back work and shifted every later climb's number;
    the published record is the only count that survives a kill in that gap."""
    full = ReactionSearch(double_well_factory, _config()).run(demo_atoms(), tmp_path / "full")
    assert full["status"] == "completed"
    assert list((tmp_path / "full" / "paths").glob("climb*")), "the run must climb"

    killed = _kill_after(tmp_path, trial_limit=10**6, climb_limit=0, name="in_climb")
    mid = json.loads((killed / "network.json").read_text(encoding="utf-8-sig"))
    assert mid["status"] == "running"
    assert (killed / "paths" / "climb001").exists(), "the kill must leave the orphan"

    resumed = ReactionSearch(double_well_factory, _config()).run(
        demo_atoms(), killed, resume=True
    )
    assert resumed["status"] == "completed"
    report = json.loads((killed / "resume.json").read_text(encoding="utf-8-sig"))
    assert "paths/climb001" in report["deleted_path_entries"]
    assert _tree(killed) == _tree(tmp_path / "full")


def test_resume_refuses_a_different_input_structure(tmp_path):
    killed = _kill_after(tmp_path, trial_limit=4)
    moved = demo_atoms()
    moved.positions[0, 0] += 0.01
    before = sorted(str(p.relative_to(killed)) for p in killed.rglob("*"))
    with pytest.raises(ValueError, match="input_structure_sha256"):
        ReactionSearch(double_well_factory, _config()).run(moved, killed, resume=True)
    # refused before the rollback: nothing on disk was touched
    assert sorted(str(p.relative_to(killed)) for p in killed.rglob("*")) == before
