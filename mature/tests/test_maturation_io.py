"""M1 io acceptance: same-config task tables rebuild identically, changed
protocols/checkpoints/geometry refuse to resume, half-written shot
directories are quarantined (not deleted), and completion counts from
commit markers only (plan M1/M4 acceptance; gate E4).
"""

import json
from dataclasses import replace

import pytest

from conftest import fake_sha, make_protocol, make_proposal, make_source_run
from mature.io import OutputConflict, OutputStore, read_json_strict, write_json_atomic
from mature.schema import ShotRecord, ShootingEstimate


@pytest.fixture
def batch(tmp_path):
    source = make_source_run(tmp_path / "source")
    return {
        "source": source,
        "out": tmp_path / "batch",
        "protocol": make_protocol(source["root"], source["files"]),
        "fingerprint": fake_sha("prrs-tree-at-batch-start"),
    }


def _init(batch):
    return OutputStore.initialize(
        batch["out"], batch["protocol"], prrs_implementation_sha256=batch["fingerprint"]
    )


# ---- lifecycle --------------------------------------------------------------


def test_initialize_creates_the_plan_layout_and_refuses_a_taken_directory(batch):
    _init(batch)
    for name in ("manifest.json", "protocol.json", "geometries", "shots", "commits"):
        assert (batch["out"] / name).exists()
    manifest = read_json_strict(batch["out"] / "manifest.json")
    assert manifest["protocol_sha256"] == batch["protocol"].sha256
    with pytest.raises(OutputConflict, match="non-empty"):
        _init(batch)


def test_initialize_verifies_source_hashes_before_writing_anything(batch):
    (batch["source"]["root"] / "network.json").write_text('{"chemical_nodes": ["tampered"]}\n')
    with pytest.raises(OutputConflict, match="content changed"):
        _init(batch)
    assert not batch["out"].exists() or not any(batch["out"].iterdir())


def test_resume_refuses_a_changed_protocol(batch):
    _init(batch)
    changed = make_protocol(
        batch["source"]["root"],
        batch["source"]["files"],
        dynamics={
            "temperature_K": 500.0,  # a different model temperature
            "timestep_fs": 0.25,
            "first_entry_k": 160,
            "free_step_cap": 1280,
        },
    )
    with pytest.raises(OutputConflict, match="protocol differs"):
        OutputStore.resume(batch["out"], changed)


def test_resume_refuses_changed_source_content(batch):
    store = _init(batch)
    (batch["source"]["root"] / "config.json").write_text('{"edited": true}\n')
    with pytest.raises(OutputConflict, match="content changed"):
        OutputStore.resume(batch["out"], batch["protocol"])
    # and the store object itself must not paper over it
    assert store.output_dir.exists()


def test_resume_refuses_changed_geometry_content(batch):
    store = _init(batch)
    proposal = make_proposal(batch["out"], "p1")
    store.add_proposal(proposal)
    (batch["out"] / proposal.geometry_file).write_text("tampered\n")
    with pytest.raises(OutputConflict, match="content changed"):
        OutputStore.resume(batch["out"], batch["protocol"])


def test_resume_refuses_a_changed_implementation_fingerprint(batch):
    _init(batch)
    with pytest.raises(OutputConflict, match="implementation fingerprint"):
        OutputStore.resume(
            batch["out"], batch["protocol"], prrs_implementation_sha256=fake_sha("other-tree")
        )


def test_resume_refuses_a_directory_that_is_not_a_batch(tmp_path):
    with pytest.raises(OutputConflict, match="no manifest"):
        OutputStore.resume(tmp_path / "empty", make_protocol(tmp_path, {"x": fake_sha("x")}))


# ---- proposals --------------------------------------------------------------


def test_add_proposal_is_idempotent_for_identical_content_and_refuses_a_clash(batch):
    store = _init(batch)
    proposal = make_proposal(batch["out"], "p1")
    store.add_proposal(proposal)
    store.add_proposal(make_proposal(batch["out"], "p1"))
    assert len(store.proposals()) == 1
    different = make_proposal(batch["out"], "p1")
    object.__setattr__(different, "n_atoms", 5)
    with pytest.raises(OutputConflict, match="id clash"):
        store.add_proposal(different)


def test_add_proposal_refuses_missing_or_mismatched_geometry(batch):
    store = _init(batch)
    proposal = make_proposal(batch["out"], "p1")
    (batch["out"] / proposal.geometry_file).unlink()
    with pytest.raises(OutputConflict, match="geometry missing"):
        store.add_proposal(proposal)
    (batch["out"] / proposal.geometry_file).write_text('{"n_atoms": 4, "name": "p1"}\n?\n')
    with pytest.raises(OutputConflict, match="hash mismatch"):
        store.add_proposal(proposal)


# ---- task table -------------------------------------------------------------


def test_task_table_rebuilds_identically_and_is_order_independent(batch):
    store = _init(batch)
    pa = make_proposal(batch["out"], "pa")
    pb = make_proposal(batch["out"], "pb")
    store.add_proposal(pa)
    store.add_proposal(pb)
    first = store.build_task_table(batch["protocol"], [pa, pb])
    again = store.build_task_table(batch["protocol"], [pb, pa])  # different call order
    assert first == again
    assert [t["seed"] for t in first] == [
        batch["protocol"].shooting.seed_base + i for i in range(len(first))
    ]
    assert [t["task_id"] for t in first] == sorted(t["task_id"] for t in first)


def test_task_table_refuses_a_conflicting_rebuild(batch):
    store = _init(batch)
    pa = make_proposal(batch["out"], "pa")
    store.add_proposal(pa)
    store.build_task_table(batch["protocol"], [pa])
    # same proposals, same shot counts, different seed table: the allocation
    # this directory was started with is the allocation it keeps
    other_shooting = replace(batch["protocol"].shooting, seed_base=999999)
    other = replace(batch["protocol"], shooting=other_shooting)
    with pytest.raises(OutputConflict, match="task table conflict"):
        store.build_task_table(other, [pa])


# ---- commit markers and completion ------------------------------------------


def _record(batch, proposal_id="p1", shot_id="s0000", outcome="A", **overrides):
    kwargs = dict(
        proposal_id=proposal_id,
        shot_id=shot_id,
        batch_role="explore",
        seed=42,
        protocol_sha256=batch["protocol"].sha256,
        implementation_sha256=fake_sha("impl"),
        outcome=outcome,
        trial_dir=f"shots/{proposal_id}/{shot_id}",
    )
    kwargs.update(overrides)
    if kwargs["outcome"] in ("A", "B") and not (
        kwargs.get("entry_side") or kwargs.get("landed_microstate")
    ):
        kwargs["entry_side"] = kwargs["outcome"]
        kwargs["entry_step"] = 200
    return ShotRecord(**kwargs)


def test_commit_then_rebuild_completion_from_markers(batch):
    store = _init(batch)
    proposal = make_proposal(batch["out"], "p1")
    store.add_proposal(proposal)
    tasks = store.build_task_table(batch["protocol"], [proposal])
    assert len(store.pending_tasks()) == len(tasks)
    store.commit(_record(batch))
    assert store.done_task_ids() == {"p1#s0000"}
    assert store.pending_tasks() == tasks[1:]


def test_recommit_identical_is_a_noop_and_a_content_clash_refuses(batch):
    store = _init(batch)
    store.commit(_record(batch))
    store.commit(_record(batch))  # same content: idempotent
    with pytest.raises(OutputConflict, match="commit conflict"):
        store.commit(_record(batch, outcome="B", entry_side="B", entry_step=10))


def test_half_written_shot_directory_is_quarantined_not_deleted(batch):
    store = _init(batch)
    stray = batch["out"] / "shots" / "p1" / "s0000"
    (stray / "trials").mkdir(parents=True)
    (stray / "result.json").write_text('{"partial": true}\n')
    destination = store.isolate_interrupted("p1", "s0000")
    assert destination is not None and "interrupted" in destination
    assert not stray.exists()
    assert (batch["out"] / "interrupted" / "p1").iterdir()
    # a committed shot is never quarantined
    store.commit(_record(batch, shot_id="s0001"))
    assert store.isolate_interrupted("p1", "s0001") is None
    assert (batch["out"] / "shots" / "p1" / "s0001").exists()


def test_resume_after_interrupt_redoes_only_the_unmarked_task(batch):
    store = _init(batch)
    proposal = make_proposal(batch["out"], "p1")
    store.add_proposal(proposal)
    tasks = store.build_task_table(batch["protocol"], [proposal])
    # a crash after writing the result but before the marker
    half = batch["out"] / "shots" / "p1" / "s0000"
    half.mkdir(parents=True)
    (half / "result.json").write_text("incomplete\n")
    reopened = OutputStore.resume(batch["out"], batch["protocol"])
    assert reopened.isolate_interrupted("p1", "s0000") is not None
    assert [t["task_id"] for t in reopened.pending_tasks()] == [
        t["task_id"] for t in tasks
    ]  # s0000 is redone; nothing was counted as done


def test_estimates_rebuild_from_markers_only(batch):
    store = _init(batch)
    assert store.estimates() == {}
    for i, outcome in enumerate(("A", "A", "A", "B")):
        store.commit(
            _record(
                batch,
                shot_id=f"s{i:04d}",
                outcome=outcome,
                entry_side=outcome,
                entry_step=200 + i,
            )
        )
    estimates = store.estimates()
    estimate = estimates["p1"]
    assert isinstance(estimate, ShootingEstimate)
    assert estimate.shots == 4
    assert estimate.counts["A"] == 3 and estimate.counts["B"] == 1
    assert estimate.conditional_denominator == 4


# ---- JSON strictness --------------------------------------------------------


def test_nan_is_refused_on_write_and_on_read(batch, tmp_path):
    with pytest.raises(ValueError, match="not JSON compliant"):
        write_json_atomic(tmp_path / "nan.json", {"x": float("nan")})
    poisoned = tmp_path / "poisoned.json"
    poisoned.write_text('{"x": NaN}\n', encoding="utf-8")
    with pytest.raises(ValueError, match="NaN"):
        read_json_strict(poisoned)


def test_written_json_is_byte_stable_for_identical_content(tmp_path):
    write_json_atomic(tmp_path / "a.json", {"b": 1, "a": [2, 3]})
    write_json_atomic(tmp_path / "a.json", {"a": [2, 3], "b": 1})
    assert (tmp_path / "a.json").read_text().count("{") == 1
    assert json.loads((tmp_path / "a.json").read_text()) == {"a": [2, 3], "b": 1}
