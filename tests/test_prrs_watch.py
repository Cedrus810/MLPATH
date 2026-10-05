"""PERF-11: the watcher re-reads a file only when it has changed.

scripts/prrs_watch.py polls every second for a whole observing session; parsing an
unchanged network.json, or a finished trial's result.json, on every frame costs in
proportion to the network and to the number of finished trials.
"""

import importlib.util
import json
import os
from pathlib import Path

import pytest

# overridable so a mutation check can point the same tests at a broken copy
SCRIPT = Path(
    os.environ.get("PRRS_WATCH_SCRIPT")
    or Path(__file__).resolve().parents[1] / "scripts" / "prrs_watch.py"
)


@pytest.fixture
def watch(monkeypatch):
    spec = importlib.util.spec_from_file_location("prrs_watch", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    reads = []
    original = module.load_json

    def counting(path):
        reads.append(os.path.basename(path))
        return original(path)

    monkeypatch.setattr(module, "load_json", counting)
    module.reads = reads
    return module


def publish(path, data):
    """The way search's atomic_json writes: a new file replaced into place."""
    temp = path.with_name(path.name + ".tmp")
    temp.write_text(json.dumps(data))
    os.replace(temp, path)


def observations(path, frames, phase="free"):
    with path.open("a") as handle:
        for k in frames:
            handle.write(json.dumps({"phase": phase, "step": k, "frame": k}) + "\n")


def test_unchanged_network_is_parsed_once(watch, tmp_path):
    publish(tmp_path / "network.json", {"status": "running", "attempts": 1})
    data = watch.RunData(str(tmp_path))
    for _ in range(3):
        data.poll()
    assert watch.reads.count("network.json") == 1
    publish(tmp_path / "network.json", {"status": "running", "attempts": 2})
    data.poll()
    assert watch.reads.count("network.json") == 2
    assert data.network["attempts"] == 2


def test_finished_trial_is_not_reread(watch, tmp_path):
    trial = tmp_path / "trials" / "t0000"
    trial.mkdir(parents=True)
    observations(trial / "observations.jsonl", range(3))
    data = watch.RunData(str(tmp_path))
    data.poll()
    publish(trial / "result.json", {"classification": "elastic"})
    for _ in range(4):
        data.poll()
    state = data.trials["t0000"]
    assert state.done and state.classification == "elastic" and state.n_frames == 3
    # the poll before it existed and the one that saw it; the finished trial costs a
    # stat from then on
    assert watch.reads.count("result.json") == 2


def test_replaced_trial_directory_is_read_again(watch, tmp_path):
    trial = tmp_path / "trials" / "t0000"
    trial.mkdir(parents=True)
    observations(trial / "observations.jsonl", range(3))
    publish(trial / "result.json", {"classification": "elastic"})
    data = watch.RunData(str(tmp_path))
    data.poll()
    data.poll()
    assert data.trials["t0000"].finished
    # a resume rolls the trial back and reruns it under the same id
    (trial / "result.json").unlink()
    fresh = trial / "observations.jsonl.new"
    observations(fresh, range(5), phase="quench")
    os.replace(fresh, trial / "observations.jsonl")
    data.poll()
    state = data.trials["t0000"]
    assert not state.done
    assert [f.phase for f in state.frames] == ["quench"] * 5
