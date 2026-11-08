import json
import socket

import numpy as np
import pytest
from ase.data import covalent_radii

from d2_bench import CONFIG, POSE, SPECS, Tick, backends, init, pose, settings, spec
from explore.d2 import init_run, open_run, run
from explore.journal import Journal
from explore.ledger import Ledger
from explore.probe import read_scan as real_read_scan
from twosaddle import SADDLE_Y, local_coords


def population(run_dir):
    lines = (run_dir / "population.jsonl").read_text(encoding="utf-8").splitlines()
    return [json.loads(line) for line in lines if line.strip()]


def two_pose_specs():
    pairs = (((1, 3), -1), ((0, 3), +1), ((0, 1), +1))
    return {p: [spec(pair, sign, pose=p) for pair, sign in pairs] for p in ("p0", "p1")}


def test_population_is_round_robin_keyed_and_quota_limited(tmp_path):
    poses = [("p0", pose()), ("p1", pose())]
    specs = two_pose_specs()
    cfg = settings(per_pose_quota=2, audit_fraction=0.5)
    init_run(tmp_path / "a", poses, specs, cfg, CONFIG, identity={})
    rows = population(tmp_path / "a")
    sampled = [r for r in rows if r["sampled"]]
    assert [r["order"][0] for r in sampled] == [0, 0, 1, 1]
    assert {r["pose_id"] for r in sampled[:2]} == {"p0", "p1"}
    assert sum(r["status"] == "not_sampled" for r in rows) == 2
    # appending an intent keeps every existing audit draw, pose key and relative order
    more = {p: v + [spec((0, 2), +1, pose=p)] for p, v in specs.items()}
    init_run(tmp_path / "b", poses, more, cfg, CONFIG, identity={})
    old = {r["intent_id"]: r for r in rows}
    for r in population(tmp_path / "b"):
        if r["intent_id"] in old:
            assert r["audit"] == old[r["intent_id"]]["audit"]
            assert r["order"][1] == old[r["intent_id"]]["order"][1]
    old_ids = [r["intent_id"] for r in rows if r["pose_id"] == "p0"]
    new_ids = [r["intent_id"] for r in population(tmp_path / "b") if r["pose_id"] == "p0"]
    assert [i for i in new_ids if i in old_ids] == old_ids


def test_domain_is_truncated_at_the_overlap_floor(tmp_path):
    # H-O1 starts at 3.2146 A; the floor is min_distance_ratio * (R_O + R_H)
    floor = CONFIG.min_distance_ratio * (covalent_radii[8] + covalent_radii[1])
    init(tmp_path / "r", specs=(spec((1, 3), -1),), magnitudes_A=(0.5, 1.0, 2.0, 3.0))
    (row,) = population(tmp_path / "r")
    q0 = 3.2146
    assert row["status"] == "ready" and row["domain_truncated"]
    assert row["deltas"] == [0.0, -0.5, -1.0, -2.0]
    assert q0 - 3.0 <= floor < q0 - 2.0


def test_a_domain_with_no_feasible_step_is_infeasible(tmp_path):
    init(tmp_path / "r", specs=(spec((1, 3), -1),), magnitudes_A=(3.0, 3.1))
    (row,) = population(tmp_path / "r")
    assert row["status"] == "infeasible" and row["deltas"] == [0.0]


def test_init_refuses_a_pose_without_electronic_state_and_an_existing_run(tmp_path):
    bare = pose()
    del bare.info["spin"]
    with pytest.raises(ValueError, match="spin"):
        init_run(tmp_path / "x", [(POSE, bare)], {POSE: list(SPECS)}, settings(), CONFIG, {})
    assert not (tmp_path / "x" / "manifest.json").exists()
    init(tmp_path / "r")
    with pytest.raises(FileExistsError):
        init(tmp_path / "r")


def test_settings_have_no_silent_values():
    with pytest.raises(ValueError):
        settings(audit_fraction=1.5)
    with pytest.raises(ValueError):
        settings(magnitudes_A=(0.2, 0.1))
    with pytest.raises(ValueError):
        settings(scan_cap=0)


def journal(run_dir):
    return Journal(run_dir / "journal.jsonl").events


def test_l1_scans_every_ready_intent_once_and_pays_exactly(tmp_path):
    init(tmp_path / "r")
    tick = Tick()
    rep = run(tmp_path / "r", backends(tick))
    settles = [e for e in journal(tmp_path / "r") if e["event"] == "settle"]
    l1 = [e for e in settles if e["task"].startswith("L1|")]
    assert len(l1) == 3 and all(e["attempt"] == 0 for e in l1)
    assert all(e["spent"] <= 600 for e in l1)
    assert rep["accounts"]["search"]["held"] == 0
    # D2 adds nothing to the count: model calls == ledger == accounts
    paid = rep["accounts"]["search"]["spent"] + rep["accounts"]["certify"]["spent"]
    assert tick.count == rep["ledger_total"] == paid


def test_coverage_truncation_is_recorded_not_hidden(tmp_path):
    init(tmp_path / "r", search_budget=600)
    rep = run(tmp_path / "r", backends(Tick()))
    events = journal(tmp_path / "r")
    assert sum(e["event"] == "settle" and e["task"].startswith("L1|") for e in events) == 1
    assert rep["skips"].get("L1_skipped:search_budget") == 2
    assert rep["accounts"]["search"]["spent"] <= 600


def test_interrupted_attempts_are_charged_and_retried(tmp_path):
    run_dir = tmp_path / "r"
    init(run_dir, specs=SPECS[:2])
    rows = population(run_dir)
    a, b = (f"L1|{r['intent_id']}" for r in rows[:2])
    j = Journal(run_dir / "journal.jsonl")
    for task in (a, b):
        j.append(
            {
                "event": "dispatch",
                "task": task,
                "attempt": 0,
                "account": "search",
                "cap": 600,
                "kind": "L1",
            }
        )
    led = Ledger(run_dir / "ledger.jsonl", arm="D2", replay=False)
    with led.segment(task=a, stage="probe"):
        for _ in range(37):
            led.charge()
    led._write({"event": "open", "task": b, "attempt": 0, "stage": "probe", "reserve": None})
    opened = open_run(run_dir)
    rec = {e["task"]: e for e in opened.journal.events if e["event"] == "reconcile"}
    assert (rec[a]["spent"], rec[a]["reconciled"]) == (37, True)
    assert (rec[b]["spent"], rec[b]["reconciled"]) == (600, False)
    assert opened.accounts["search"].spent == 637
    rep = run(run_dir, backends(Tick()))
    settled = {(e["task"], e["attempt"]) for e in journal(run_dir) if e["event"] == "settle"}
    assert (a, 1) in settled and (b, 1) in settled
    assert len([e for e in journal(run_dir) if e["event"] == "reconcile"]) == 2
    assert rep["unreconciled"][0]["task"] == b


def force_jump(monkeypatch):
    def fake(points, noise_eV):
        if len(points) >= 2:
            return [{"kind": "jump", "before": 0, "after": 1}]
        return real_read_scan(points, noise_eV)

    monkeypatch.setattr("explore.d2.read_scan", fake)


def test_a_jump_gets_one_reverse_scan_from_the_post_jump_point(tmp_path, monkeypatch):
    run_dir = tmp_path / "r"
    init(run_dir, specs=(SPECS[2],))
    force_jump(monkeypatch)
    rep = run(run_dir, backends(Tick()))
    iid = population(run_dir)[0]["intent_id"]
    settles = {e["task"]: e for e in journal(run_dir) if e["event"] == "settle"}
    l1, l2 = settles[f"L1|{iid}"], settles[f"L2|{iid}"]
    assert [t for t in settles if t.startswith("L2|")] == [f"L2|{iid}"]
    assert l2["cap"] == 200 and l2["attempt"] == 0
    assert abs(l2["payload"]["q0"] - l1["payload"]["points"][1]["q_actual"]) < 1e-9
    assert len(rep["l2_status"]) == 1


def test_reverse_scans_are_skipped_when_the_search_account_cannot_hold_them(
    tmp_path, monkeypatch
):
    run_dir = tmp_path / "r"
    init(run_dir, specs=(SPECS[2],), search_budget=600, reverse_cap=600)
    force_jump(monkeypatch)
    rep = run(run_dir, backends(Tick()))
    assert rep["skips"]["L2_skipped:search_budget"] == 1
    assert not [
        e for e in journal(run_dir) if e["event"] == "settle" and e["task"].startswith("L2|")
    ]


def test_resume_after_completion_spends_nothing_twice(tmp_path):
    run_dir = tmp_path / "r"
    init(run_dir)
    first = run(run_dir, backends(Tick()))
    n = len(journal(run_dir))
    tick = Tick()
    second = run(run_dir, backends(tick))
    assert tick.count == 0
    assert len(journal(run_dir)) == n
    assert second["accounts"] == first["accounts"]


def l4(run_dir):
    return [
        e for e in journal(run_dir) if e["event"] == "settle" and e["task"].startswith("L4|")
    ]


def test_ridges_are_released_and_certified_on_two_saddles(tmp_path):
    init(tmp_path / "r")
    rep = run(tmp_path / "r", backends(Tick()))
    assert rep["promotion"].get("ridge_suspected", 0) >= 1
    assert rep["deferred"].get("no_trigger", 0) >= 1  # O0-O1 stretch
    assert len(rep["c2"]) >= 1
    for e in l4(tmp_path / "r"):
        assert e["payload"]["class"] == "ridge_suspected"
        assert e["status"] in ("c2", "duplicate_c2")
    seen = []
    for e in l4(tmp_path / "r"):
        if e["status"] == "c2":
            _, y, _ = local_coords(e["payload"]["saddle_positions"])
            k = min(range(2), key=lambda i: abs(SADDLE_Y[i] - y))
            assert abs(SADDLE_Y[k] - y) < 1e-3
            assert k not in seen
            seen.append(k)
    c1 = rep["release_status"].get("c1_candidate", 0)
    assert c1 == rep["promotion"]["ridge_suspected"]


def test_audit_items_are_flagged_and_paid_by_d2(tmp_path):
    init(tmp_path / "r", specs=(SPECS[2],), audit_fraction=1.0)
    tick = Tick()
    rep = run(tmp_path / "r", backends(tick))
    assert rep["promotion"] == {"audit": 1}
    done = [e for e in journal(tmp_path / "r") if e["event"] == "settle"]
    tags = [e["payload"]["audit_release"] for e in done if e["task"][:2] in ("L3", "L4")]
    assert tags == [True, True]
    paid = rep["accounts"]["search"]["spent"] + rep["accounts"]["certify"]["spent"]
    assert tick.count == paid


def test_search_exhaustion_leaves_the_certification_reserve(tmp_path):
    init(tmp_path / "r", search_budget=1800, release_cap=1800)
    rep = run(tmp_path / "r", backends(Tick()))
    assert (
        rep["skips"].get("L3_skipped:search_budget", 0) == rep["promotion"]["ridge_suspected"]
    )
    assert len(l4(tmp_path / "r")) == rep["promotion"]["ridge_suspected"]


def test_certification_budget_exhaustion_is_a_skip_not_a_loss(tmp_path):
    init(tmp_path / "r", certify_budget=100)
    rep = run(tmp_path / "r", backends(Tick()))
    assert l4(tmp_path / "r") == []
    assert rep["skips"]["L4_skipped:certify_budget"] == rep["promotion"]["ridge_suspected"]
    assert rep["accounts"]["certify"]["spent"] == 0
    assert rep["release_status"].get("c1_candidate", 0) >= 1


def test_resume_after_completion_dispatches_nothing(tmp_path):
    init(tmp_path / "r")
    first = run(tmp_path / "r", backends(Tick()))
    n = len(journal(tmp_path / "r"))
    tick = Tick()
    second = run(tmp_path / "r", backends(tick))
    assert tick.count == 0 and len(journal(tmp_path / "r")) == n
    assert second["c2"] == first["c2"]


def test_journal_events_carry_host_and_time_and_the_report_lists_hosts(tmp_path):
    init(tmp_path / "r")
    rep = run(tmp_path / "r", backends(Tick()))
    events = journal(tmp_path / "r")
    assert events and all(e["host"] == socket.gethostname() and "t" in e for e in events)
    assert rep["hosts"] == [socket.gethostname()]
    manifest = json.loads((tmp_path / "r" / "manifest.json").read_text())
    assert manifest["host"] == socket.gethostname()


def test_report_counts_attempts_and_c1_labels(tmp_path):
    init(tmp_path / "r")
    rep = run(tmp_path / "r", backends(Tick()))
    assert rep["attempts"] == {}  # uninterrupted: every task dispatched once
    assert sum(rep["c1_by_label"].values()) == len(rep["c1"])


def test_f6_the_scan_tangent_reaches_the_refiner(tmp_path, monkeypatch):
    import explore.d2 as d2

    seen = []
    real = d2.certify

    def spy(*args, **kwargs):
        seen.append(kwargs.get("direction"))
        return real(*args, **kwargs)

    monkeypatch.setattr(d2, "certify", spy)
    init(tmp_path / "r", specs=(SPECS[0],))
    run(tmp_path / "r", backends(Tick()))
    assert seen and all(d is not None and np.asarray(d).shape == (4, 3) for d in seen)
    assert all(np.linalg.norm(d) > 0 for d in seen)
