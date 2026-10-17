"""M2-M4 on a real (small) PRRS run: the proton-relay testbed, never a neural
checkpoint. The source run and its PARS selection are produced by the real
programs in tmp, once per module (about 15 s); nothing reads ``runs/``.

What is pinned (plan section 4 and gates E1-E4):
- candidates come from the source trajectory's free segment, all atoms kept,
  interpolation labelled as such, uncovered channels/tubes in the ledger;
- a batch shot is the machine's own shot: same geometry and seed give the same
  first entry and endpoint as calling ``prrs.runner.shoot`` directly, and no
  probe is delivered;
- judging keeps every excluded label apart from A/B;
- the explore stop rule and the budget act before a shot, never after;
- a killed shot is redone, a re-report never counts a shot twice.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from ase.io import read as ase_read

pytest.importorskip("prrs")

from mature import maturation  # noqa: E402
from mature.cli import main  # noqa: E402
from mature.io import OutputConflict, OutputStore, sha256_file  # noqa: E402
from mature.schema import MaturationProtocol, ShotRecord  # noqa: E402

from conftest import fake_sha  # noqa: E402

PARENT = Path(__file__).resolve().parents[2]
PARS_SRC = PARENT / "pars" / "src"


@pytest.fixture(scope="module")
def source(tmp_path_factory):
    if not PARS_SRC.is_dir():
        pytest.skip("standalone extraction: no sibling pars tree to select with")
    from prrs.calculators import proton_relay_atoms, proton_relay_factory
    from prrs.config import SearchConfig
    from prrs.search import ReactionSearch

    root = tmp_path_factory.mktemp("relay")
    config = SearchConfig(
        families=("stretch", "compress"),
        closed_shell_only=False,
        max_depth=1,
        max_trials=16,
        seed=7,
    )
    ReactionSearch(proton_relay_factory, config).run(
        proton_relay_atoms((0, 0, 0)), root / "run"
    )
    env = {**os.environ, "PYTHONPATH": str(PARS_SRC)}
    subprocess.run(
        [
            sys.executable,
            "-m",
            "pars",
            "select",
            "--run-dir",
            str(root / "run"),
            "--out",
            str(root / "selection"),
        ],
        check=True,
        env=env,
        capture_output=True,
    )
    return root


def _protocol(source, **changes):
    run, selection = source / "run", source / "selection"
    manifest = json.loads((run / "config.json").read_text())
    sides = {
        side: maturation.microstate_ref(run, state).to_dict()
        for side, state in (("A", "c0000/m0000"), ("B", "c0001/m0000"))
    }
    data = {
        "batch_role": "explore",
        "source": {
            "run_dir": str(run),
            "implementation_sha256": manifest["implementation_sha256"],
            "model_factory": "prrs.calculators:proton_relay_factory",
            "model_kwargs": {},
            "pars_selection_dir": str(selection),
            "pars_selection_sha256": sha256_file(selection / "pars_selection.json"),
            "files": {
                name: sha256_file(run / name)
                for name in ("config.json", "network.json", "attempts.jsonl")
            },
        },
        "regions": {"sides": sides},
        "dynamics": {
            "temperature_K": 300.0,
            "timestep_fs": manifest["config"]["timestep_fs"],
            "first_entry_k": 20,
            "free_step_cap": 400,
        },
        "candidates": {},
        "shooting": {
            "batch_role": "explore",
            "seed_base": 1000,
            "shots_per_candidate": 4,
            "drift_tolerance_eV_atom": 0.01,
        },
        "budget": {"max_shots": 200},
    }
    for key, value in changes.items():
        data[key] = {**data[key], **value} if isinstance(value, dict) else value
    return MaturationProtocol.from_dict(data)


# ---- M2 -----------------------------------------------------------------------


def test_candidates_come_from_the_free_segment_with_every_atom(source, tmp_path):
    protocol = _protocol(source)
    summary = maturation.prepare(protocol, tmp_path / "batch")
    store = OutputStore(tmp_path / "batch")
    proposals = store.proposals()
    assert len(proposals) == 5 and summary["tasks"] == 20
    reference = ase_read(str(source / "run" / protocol.regions.sides["A"].structure_path))
    for proposal in proposals:
        frames = ase_read(
            str(source / "run" / "trials" / proposal.source_attempt_id / "response.traj"), ":"
        )
        phases = [
            frames[i].info["phase"]
            for i in range(proposal.parent_frame_start, proposal.parent_frame_stop)
        ]
        assert set(phases) <= {"free_start", "free"}  # no forcing, no quench
        assert 0.0 < proposal.arc_fraction < 1.0  # endpoints are never proposed
        geometry = ase_read(str(tmp_path / "batch" / proposal.geometry_file))
        assert np.array_equal(geometry.numbers, reference.numbers)  # every H kept
        assert proposal.proposal_kind == "path_interpolation"
        frames_used = proposal.interpolation["frames"]
        assert all(
            proposal.parent_frame_start <= f < proposal.parent_frame_stop for f in frames_used
        )
        w = proposal.interpolation["weights"]
        mixed = (
            w[0] * frames[frames_used[0]].positions + w[1] * frames[frames_used[1]].positions
        )
        assert np.allclose(geometry.positions, mixed)
    ledger = summary["ledger"]["channels"]
    assert ledger["rc0000"]["status"] == "sampled"
    assert {ledger[c]["status"] for c in ledger if c != "rc0000"} == {"not_this_region"}


def test_without_interpolation_candidates_are_recorded_frames(source, tmp_path):
    protocol = _protocol(source, candidates={"allow_interpolation": False})
    maturation.prepare(protocol, tmp_path / "batch")
    for proposal in OutputStore(tmp_path / "batch").proposals():
        assert proposal.proposal_kind == "observed_frame"
        frame = ase_read(
            str(source / "run" / "trials" / proposal.source_attempt_id / "response.traj"),
            proposal.frame_index,
        )
        geometry = ase_read(str(tmp_path / "batch" / proposal.geometry_file))
        assert np.allclose(geometry.positions, frame.positions)


def test_region_edges_must_be_what_the_structure_encodes(source, tmp_path):
    protocol = _protocol(source)
    sides = protocol.regions.to_dict()["sides"]
    sides["B"]["edges"] = sides["B"]["edges"][:-1]
    tampered = _protocol(source, regions={"sides": sides})
    with pytest.raises(OutputConflict, match="does not encode"):
        maturation.prepare(tampered, tmp_path / "batch")
    assert not (tmp_path / "batch").exists()


def test_dry_run_writes_nothing(source, tmp_path):
    result = maturation.prepare(_protocol(source), tmp_path / "batch", dry_run=True)
    assert result["tasks"] == 20
    assert not (tmp_path / "batch").exists()


# ---- M5: explicit candidate lists ------------------------------------------------

ZERO = {
    "batch_role": "zero_measure",
    "shots_per_candidate": 2,
    "drift_tolerance_eV_atom": None,
}


def _listed_ids(source):
    return maturation.prepare(_protocol(source), "unused", dry_run=True)["proposals"]


def test_a_listed_candidate_is_the_only_one_prepared(source, tmp_path):
    chosen = _listed_ids(source)[2]
    protocol = _protocol(
        source,
        batch_role="zero_measure",
        shooting=ZERO,
        candidates={"proposal_ids": [chosen]},
    )
    summary = maturation.prepare(protocol, tmp_path / "batch")
    assert summary["proposals"] == [chosen] and summary["tasks"] == 2
    assert [p.name for p in (tmp_path / "batch" / "geometries").iterdir()] == [
        f"{chosen}.extxyz"
    ]
    (path,) = summary["ledger"]["channels"]["rc0000"]["paths"]
    skipped = [p for p in path["positions"] if "skipped" in p]
    assert len(skipped) == 4
    assert {p["skipped"] for p in skipped} == {"not in CandidateRule.proposal_ids"}


def test_a_listed_id_that_is_not_proposed_refuses_before_writing(source, tmp_path):
    protocol = _protocol(source, candidates={"proposal_ids": ["mp-rc0000-t9999-k0"]})
    with pytest.raises(OutputConflict, match="not proposed from this selection"):
        maturation.prepare(protocol, tmp_path / "batch")
    assert not (tmp_path / "batch").exists()


def _confirm(source, explore_dir, ids, **changes):
    return _protocol(
        source,
        batch_role="confirm",
        shooting={
            "batch_role": "confirm",
            "seed_base": 5000,
            "shots_per_candidate": None,
            "confirm_shots": 2,
            "selection_rule": "test: the listed explore candidates",
        },
        candidates={"from_batch": str(explore_dir), "proposal_ids": ids},
        **changes,
    )


def test_a_confirm_batch_shoots_the_explore_batchs_own_geometry(source, tmp_path):
    explore = tmp_path / "explore"
    maturation.prepare(_protocol(source), explore)
    original = {p.proposal_id: p for p in OutputStore(explore).proposals()}
    chosen = sorted(original)[3]
    out = tmp_path / "confirm"
    maturation.prepare(_confirm(source, explore, [chosen]), out)
    (proposal,) = OutputStore(out).proposals()
    assert proposal.batch_role == "confirm"
    assert (out / proposal.geometry_file).read_bytes() == (
        explore / original[chosen].geometry_file
    ).read_bytes()
    manifest = json.loads((explore / "manifest.json").read_text())
    assert proposal.origin_batch == {
        "dir": str(explore),
        "protocol_sha256": manifest["protocol_sha256"],
        "batch_role": "explore",
    }
    kept = {
        k: v for k, v in proposal.to_dict().items() if k not in ("batch_role", "origin_batch")
    }
    assert kept == {k: v for k, v in original[chosen].to_dict().items() if k != "batch_role"}
    maturation.shoot_batch(out, log=lambda _: None)
    (row,) = maturation.report(out)["candidates"]
    assert row["estimate"]["shots"] == 2
    handoff = json.loads((out / "refinement_handoff.json").read_text())
    assert all(c["origin_batch"]["dir"] == str(explore) for c in handoff["candidates"])


def test_a_confirm_batch_refuses_what_explore_did_not_hold_or_judged_elsewhere(
    source, tmp_path
):
    explore = tmp_path / "explore"
    maturation.prepare(_protocol(source), explore)
    chosen = OutputStore(explore).proposals()[0]
    with pytest.raises(OutputConflict, match="are not in"):
        maturation.prepare(_confirm(source, explore, ["mp-nope"]), tmp_path / "c1")
    sides = _protocol(source).regions.to_dict()["sides"]
    swapped = {"sides": {"A": sides["B"], "B": sides["A"]}}
    with pytest.raises(OutputConflict, match="another A/B region"):
        maturation.prepare(
            _confirm(source, explore, [chosen.proposal_id], regions=swapped), tmp_path / "c2"
        )
    (explore / chosen.geometry_file).write_text("tampered")
    with pytest.raises(OutputConflict, match="content changed"):
        maturation.prepare(_confirm(source, explore, [chosen.proposal_id]), tmp_path / "c3")
    assert not any((tmp_path / name).exists() for name in ("c1", "c2", "c3"))


# ---- M3 -----------------------------------------------------------------------


def test_a_batch_shot_is_the_machines_own_shot(source, tmp_path):
    """E2: same geometry, same seed -> the same first entry and endpoint as a
    direct prrs.runner.shoot; and the probe is never delivered."""
    from prrs.runner import shoot

    protocol = _protocol(source)
    maturation.prepare(protocol, tmp_path / "batch")
    maturation.shoot_batch(tmp_path / "batch", max_new_shots=1, log=lambda _: None)
    store = OutputStore(tmp_path / "batch")
    (record,) = store.committed_records()
    task = store.tasks()[0]
    proposal = next(p for p in store.proposals() if p.proposal_id == task["proposal_id"])
    trial = tmp_path / "batch" / record.trial_dir
    batch_result = json.loads((trial / "result.json").read_text())
    assert batch_result["probe_delivery"] == "skipped:committor_shot"

    context = maturation.shot_context(protocol, store)
    atoms = ase_read(str(tmp_path / "batch" / proposal.geometry_file))
    atoms.info = {}
    atoms.calc = context.counted.calculator
    (direct,) = shoot(
        atoms,
        maturation._probe_for(context, proposal, task["seed"]),
        lambda: context.counted.calculator,
        context.config,
        tmp_path / "direct",
        1,
        task["seed"],
        first_entry=context.first_entry,
    )
    assert direct.record["first_entry"] == batch_result["first_entry"]
    assert record.entry_step == batch_result["first_entry"]["step"]
    endpoint = ase_read(str(trial / "endpoint.extxyz"))
    assert np.allclose(direct.endpoint.positions, endpoint.positions)


def _trial(status="completed", side=None, step=None, converged=True, code=None):
    return {
        "status": status,
        "failure": {"code": code, "type": "X", "message": "m"} if code else None,
        "first_entry": {"side": side, "step": step, "k": 20, "time_fs": 1.0},
        "quench_converged": converged,
        "max_nve_drift_eV_atom": 0.001,
    }


SIDES = {"A": "c0000/m0000", "B": "c0001/m0000"}


@pytest.mark.parametrize(
    ("trial", "landed", "outcome", "consistent"),
    [
        (_trial(side="A", step=30), "c0000/m0000", "A", True),
        (_trial(side="B", step=30), "c0000/m0000", "B", False),  # counted, not re-judged
        (_trial("quench_failed", side="B", step=30, converged=False), None, "B", False),
        (_trial(), "c0002/m0000", "other", None),  # third basin
        (_trial(), None, "other", None),  # converged, unidentified
        (_trial(), "c0000/m0000", "unabsorbed", None),  # hovered, never k steps
        (_trial("quench_failed", converged=False), None, "unabsorbed", None),
        (_trial("rejected", code="nve_drift"), None, "gate_rejected", None),
        (_trial("failed", code="trial_error"), None, "failed", None),
    ],
)
def test_judging_keeps_every_excluded_label_apart(trial, landed, outcome, consistent):
    fields = maturation.judge(trial, landed, SIDES, 400)
    assert fields["outcome"] == outcome
    assert fields.get("entry_quench_consistent") is consistent


# ---- M4: stop rule, budget, resume, report --------------------------------------


def _record(protocol, shot, outcome):
    extra = (
        {"entry_side": outcome, "entry_step": 30}
        if outcome in ("A", "B")
        else {"outcome_reason": "x"}
    )
    return ShotRecord(
        proposal_id="p1",
        shot_id=f"s{shot:04d}",
        batch_role=protocol.batch_role,
        seed=shot,
        protocol_sha256=fake_sha("p"),
        implementation_sha256=fake_sha("i"),
        outcome=outcome,
        trial_dir="t",
        **extra,
    )


def test_the_explore_stop_rule_needs_the_interval_outside_the_band(source):
    protocol = _protocol(source, shooting={"shots_per_candidate": 16, "interim_look": 8})
    same_side = [_record(protocol, i, "A") for i in range(8)]
    assert maturation.stop_reason(protocol, same_side) == "interim_side"  # 0/8: upper 0.324
    one_across = same_side[:7] + [_record(protocol, 7, "B")]
    assert maturation.stop_reason(protocol, one_across) is None  # 1/8: upper 0.471
    assert maturation.stop_reason(protocol, same_side[:7]) is None  # look not complete
    excluded = [_record(protocol, 0, "unabsorbed"), _record(protocol, 1, "other")]
    assert maturation.stop_reason(protocol, excluded) == "excluded_share_already_failed"


def test_confirm_batches_take_every_planned_shot(source):
    protocol = _protocol(
        source,
        batch_role="confirm",
        shooting={
            "batch_role": "confirm",
            "shots_per_candidate": None,
            "confirm_shots": 16,
            "selection_rule": "frozen",
        },
    )
    excluded = [_record(protocol, i, "unabsorbed") for i in range(5)]
    assert maturation.stop_reason(protocol, excluded) is None


def test_the_budget_stops_before_a_shot_not_after(source, tmp_path):
    protocol = _protocol(source, budget={"max_shots": 3})
    maturation.prepare(protocol, tmp_path / "batch")
    result = maturation.shoot_batch(tmp_path / "batch", log=lambda _: None)
    assert result["state"] == "budget_exhausted"
    assert result["fired"] == 3
    assert result["budget"]["blocked_by"] == ["max_shots"]


def test_a_killed_shot_is_redone_and_a_rereport_counts_nothing_twice(source, tmp_path):
    out = tmp_path / "batch"
    maturation.prepare(_protocol(source), out)
    maturation.shoot_batch(out, max_new_shots=2, log=lambda _: None)
    store = OutputStore(out)
    task = store.pending_tasks()[0]
    half = out / "shots" / task["proposal_id"] / task["shot_id"] / "trials" / "shot0000"
    half.mkdir(parents=True)
    (half / "response.traj").write_bytes(b"partial")
    maturation.shoot_batch(out, max_new_shots=1, log=lambda _: None)
    assert (out / "interrupted" / task["proposal_id"]).is_dir()
    assert len(store.committed_records()) == 3
    first = maturation.report(out)
    second = maturation.report(out)
    assert first["candidates"] == second["candidates"]
    assert first["budget"]["used"]["shots"] == 3


def test_a_full_batch_reports_sides_and_hands_off_nothing_as_verified(source, tmp_path):
    out = tmp_path / "batch"
    maturation.prepare(_protocol(source), out)
    maturation.shoot_batch(out, log=lambda _: None)
    summary = maturation.report(out)
    labels = [row["diagnosis"]["label"] for row in summary["candidates"]]
    assert len(labels) == 5
    assert set(labels) <= {"A_side", "B_side", "transition_region_compatible", "unresolved"}
    for row in summary["candidates"]:
        counts = row["estimate"]["counts"]
        assert sum(counts.values()) == 4
        assert row["estimate"]["conditional_denominator"] == counts["A"] + counts["B"]
    handoff = json.loads((out / "refinement_handoff.json").read_text())
    assert "none of them is a verified saddle" in handoff["meaning"]
    lines = (out / "diagnostics.jsonl").read_text().splitlines()
    assert len(lines) == 5


def test_a_changed_prrs_tree_cannot_continue_the_batch(source, tmp_path, monkeypatch):
    out = tmp_path / "batch"
    maturation.prepare(_protocol(source), out)
    monkeypatch.setattr(maturation, "prrs_implementation_sha256", lambda: fake_sha("other"))
    with pytest.raises(OutputConflict, match="implementation fingerprint changed"):
        maturation.shoot_batch(out, log=lambda _: None)


def test_cli_runs_the_three_verbs(source, tmp_path, capsys):
    protocol_path = tmp_path / "protocol.json"
    protocol_path.write_text(json.dumps(_protocol(source).to_dict()))
    out = tmp_path / "batch"
    assert (
        main(["prepare", "--protocol", str(protocol_path), "--out", str(out), "--dry-run"]) == 0
    )
    assert not out.exists()
    assert main(["prepare", "--protocol", str(protocol_path), "--out", str(out)]) == 0
    assert main(["shoot", "--out", str(out), "--max-new-shots", "2"]) == 0
    assert main(["report", "--out", str(out)]) == 0
    assert "mature report" in capsys.readouterr().out
