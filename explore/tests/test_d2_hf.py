"""D2-HF-v1 on the TwoSaddle bench (docs/PLAN_HF4_D2_2026-10-26.md H2/H3): E/F search
and E/F connection replace L4; external counter == ledgers; kill->resume inside HS,
between HS and HC, and inside HC pays once and keeps the search geometry."""

import json
import os
import subprocess
import sys
from collections import Counter
from pathlib import Path

import pytest
from d2_bench import HF, Tick, backends, init
from twosaddle import SADDLE_Y, local_coords

from explore.d2 import HF_RULE_VERSION, HFSettings, run

ROOT = Path(__file__).resolve().parents[1]


def drive(run_dir, counter, crash_after, mode="hf", env=None):
    path = os.pathsep.join(str(ROOT / p) for p in ("src", "../src", "tests"))
    return subprocess.run(
        [
            sys.executable,
            str(ROOT / "tests" / "d2_driver.py"),
            str(run_dir),
            str(counter),
            str(crash_after),
            mode,
        ],
        env={**os.environ, "PYTHONPATH": path, **(env or {})},
        capture_output=True,
        text=True,
        timeout=900,
    )


def events(run_dir):
    return [json.loads(x) for x in (run_dir / "journal.jsonl").read_text().splitlines() if x]


def report(run_dir):
    return json.loads((run_dir / "report.json").read_text())


def paid(rep):
    return rep["accounts"]["search"]["spent"] + rep["accounts"]["refine"]["spent"]


@pytest.fixture(scope="module")
def reference(tmp_path_factory):
    base = tmp_path_factory.mktemp("hfref")
    done = drive(base / "run", base / "count", -1)
    assert done.returncode == 0, done.stderr
    return base / "run", (base / "count").stat().st_size


def test_hf_settings_have_no_silent_values():
    with pytest.raises(TypeError):
        HFSettings(reaction_force_sigma=1e-3, geometry_noise_A=1e-3)
    with pytest.raises(ValueError):
        HFSettings(1e-3, 0.0, 1e-3, 100, 100)
    with pytest.raises(ValueError):
        HFSettings(1e-3, 1e-3, 1e-3, 100.5, 100)


def test_hf_run_finds_and_connects_both_saddles_with_exact_accounting(reference):
    run_dir, external = reference
    rep = report(run_dir)
    assert rep["rule_version"] == HF_RULE_VERSION
    assert "certify_status" in rep and rep["certify_status"] == {}
    assert rep["hs_status"].get("converged_candidate", 0) >= 1
    assert rep["hc_status"].get("connected", 0) >= 1
    assert rep["property_status"] == "not_run"
    seen = set()
    for e in rep["hf_events"]:
        _, y, _ = local_coords(
            next(
                ev["payload"]["saddle_positions"]
                for ev in events(run_dir)
                if ev["event"] == "settle" and ev["task"] == e["task"]
            )
        )
        k = min(range(2), key=lambda i: abs(SADDLE_Y[i] - y))
        assert abs(SADDLE_Y[k] - y) < 1e-3
        seen.add(k)
    assert seen  # at least one of the two analytic saddles, located on E/F only
    assert external == rep["ledger_total"] == paid(rep)
    assert rep["unreconciled"] == []
    assert rep["hf_forwards"]["HS"] > 0 and rep["hf_forwards"]["HC"] > 0


def test_resume_after_completion_dispatches_nothing(reference, tmp_path):
    run_dir, _ = reference
    before = len(events(run_dir))
    rep = run(run_dir, backends(Tick()))
    assert len(events(run_dir)) == before
    assert rep["hc_status"] == report(run_dir)["hc_status"]


def crash_point(ref_dir, where):
    log = events(ref_dir)
    if where == "GN":
        gn = next(e for e in log if e["event"] == "settle" and e["task"].startswith("GN|"))
        return gn["cumulative_spent"] - gn["spent"] // 2
    hs = next(e for e in log if e["event"] == "settle" and e["task"].startswith("HS|"))
    hc = next(e for e in log if e["event"] == "settle" and e["task"].startswith("HC|"))
    if where == "HS":
        return hs["cumulative_spent"] - hs["spent"] // 2
    if where == "HC":
        return hc["cumulative_spent"] - hc["spent"] // 2
    # "between": no forward happens between an HS settle and the next HC dispatch, so the
    # earliest kill a forward counter can place after the HS settle is the first forward
    # of that intent's HC -- the search is on disk, the HC segment is open (charged cap)
    return hs["cumulative_spent"] + 1


@pytest.mark.parametrize("where", ["HS", "between", "HC"])
def test_kill_and_resume_pays_once_and_keeps_the_search_geometry(reference, tmp_path, where):
    _kill_and_resume(reference, tmp_path, where, "hf")


def _kill_and_resume(reference, tmp_path, where, mode, env=None):
    ref_dir, ref_external = reference
    run_dir, counter = tmp_path / "run", tmp_path / "count"
    first = drive(run_dir, counter, crash_point(ref_dir, where), mode, env)
    assert first.returncode == 137, first.stderr
    cut = events(run_dir)
    hs_before = {
        e["task"]: e["payload"]["refine"]["positions"]
        for e in cut
        if e["event"] == "settle"
        and e["task"].startswith("HS|")
        and e["status"] == "converged_candidate"
    }

    second = drive(run_dir, counter, -1, mode, env)
    assert second.returncode == 0, second.stderr
    log, rep, ref = events(run_dir), report(run_dir), report(ref_dir)

    # same science as the uninterrupted run
    assert rep["hs_status"] == ref["hs_status"]
    assert rep["hc_status"] == ref["hc_status"]
    assert sorted(round(e["saddle_energy_eV"], 6) for e in rep["hf_events"]) == sorted(
        round(e["saddle_energy_eV"], 6) for e in ref["hf_events"]
    )

    # a settled search is never redone: its geometry is the one written before the kill
    for task, pos in hs_before.items():
        settles = [e for e in log if e["event"] == "settle" and e["task"] == task]
        assert len(settles) == 1 and settles[0]["payload"]["refine"]["positions"] == pos

    # every task settles once; interrupted attempts are reconciled and retried
    assert max(Counter(e["task"] for e in log if e["event"] == "settle").values()) == 1
    reconciles = [e for e in log if e["event"] == "reconcile"]
    for r in reconciles:
        assert any(
            e["event"] == "settle" and e["task"] == r["task"] and e["attempt"] > r["attempt"]
            for e in log
        )
    assert reconciles

    # cost: nothing refunded, nothing paid twice
    wasted = sum(r["spent"] for r in reconciles)
    assert paid(rep) == paid(ref) + wasted
    external = counter.stat().st_size
    assert external <= paid(rep)
    assert external - ref_external <= wasted


def test_v11_runs_are_untouched_by_hf_mode(tmp_path):
    init(tmp_path / "old")
    rep = run(tmp_path / "old", backends(Tick()))
    assert rep["rule_version"] != HF_RULE_VERSION
    assert "hs_status" not in rep and "certify" in rep["accounts"]
    assert HF.hs_cap == 3000


def test_hf_duplicate_merges_a_symmetry_equivalent_saddle():
    """b05 E2 (DONE_LOG 2026-10-28): saddles transferring equivalent methyl H were kept as
    separate events because the duplicate rule took the symmetry of the saddle graph only."""
    from types import SimpleNamespace

    from prrs.config import SearchConfig
    from test_certify import _bridged_methyl

    from explore.d2 import _hf_duplicate

    minimum, a, b = _bridged_methyl()
    product = a.copy()
    product.positions[1] = (1.9, 0.0, -0.6)
    pair = ("k1", "k2")
    run = SimpleNamespace(
        config=SearchConfig(closed_shell_only=False),
        settled={
            "HC|first": {
                "status": "connected",
                "payload": {
                    "chemical_pair": list(pair),
                    "saddle_energy_eV": -1.0,
                    "saddle_positions": b.positions.tolist(),
                },
            }
        },
    )
    ends = (minimum.positions, product.positions)
    assert _hf_duplicate(run, minimum, a.positions, -1.0, pair, ends) == "HC|first"
    assert _hf_duplicate(run, minimum, a.positions, -1.0, ("k1", "k3"), ends) is None


def test_start_cluster_settings_and_freeze(tmp_path):
    from explore.d2 import StartCluster

    for bad in ((0.0, 1), (float("nan"), 1), (2.0, 0), (2.0, 1.5), (True, 1)):
        with pytest.raises(ValueError):
            StartCluster(*bad)
    with pytest.raises(ValueError):
        init(tmp_path / "no_hf", start_cluster=StartCluster(2.0, 1))
    init(tmp_path / "sc", hf=HF, start_cluster=StartCluster(2.0, 1))
    manifest = json.loads((tmp_path / "sc" / "manifest.json").read_text())
    assert manifest["rule_version"] == "D2-HF-v1.1"
    assert manifest["start_cluster"] == {"noise_multiple": 2.0, "first_round": 1}


def test_rmsd_lower_bound_never_exceeds_the_symmetric_rmsd():
    import numpy as np
    from prrs.chemistry import symmetric_rmsd
    from test_certify import _bridged_methyl

    from d2_bench import CONFIG

    from explore.certify import pair_symmetries
    from explore.d2 import _rmsd_lower_bound

    minimum, a, b = _bridged_methyl()
    group, _ = pair_symmetries(minimum, (a,), CONFIG)
    index = np.arange(len(a))
    rng = np.random.default_rng(3)
    for _ in range(50):
        c = a.copy()
        c.positions = a.positions + rng.normal(0, 0.3, a.positions.shape)
        for other in (b, minimum, c):
            lb = _rmsd_lower_bound(c, other, index)
            assert lb <= symmetric_rmsd(c, other, group, None) + 1e-12
    assert _rmsd_lower_bound(a, b, index) < 1e-12  # relabelled + rotated copy


def _one_cluster(tmp_path, **over):
    """All promoted items in one start cluster: radius 1e4 x 3.3e-4 A = 3.3 A."""
    from explore.d2 import StartCluster

    run_dir = tmp_path / "run"
    init(run_dir, hf=HF, start_cluster=StartCluster(1e4, 1), **over)
    tick = Tick()
    rep = run(run_dir, backends(tick))
    return run_dir, rep, tick.count


def test_start_cluster_with_ample_budget_runs_every_hs(tmp_path, reference):
    run_dir, rep, external = _one_cluster(tmp_path)
    ref = report(reference[0])
    assert rep["rule_version"] == "D2-HF-v1.1"
    sc = rep["start_clusters"]
    assert sc["clusters"] == 1 and sc["first_round"] == 1
    assert sc["first_round"] + sc["second_round"] == sum(ref["hs_status"].values())
    assert rep["hs_status"] == ref["hs_status"]  # order changed, nothing dropped
    assert not any(k.startswith("HS_") for k in rep["skips"])
    assert external == rep["ledger_total"] == paid(rep)
    tagged = [
        e["payload"]["start_cluster"]
        for e in events(run_dir)
        if e["event"] == "settle" and e["task"].startswith("HS|")
    ]
    assert sorted(t["rank"] for t in tagged) == list(range(len(tagged)))


def test_start_cluster_defers_the_rest_when_round_one_spends_the_budget(tmp_path):
    run_dir, rep, external = _one_cluster(tmp_path, certify_budget=HF.hs_cap)
    sc = rep["start_clusters"]
    assert sc["second_round"] >= 1
    hs_settled = [
        e for e in events(run_dir) if e["event"] == "settle" and e["task"].startswith("HS|")
    ]
    assert len(hs_settled) == 1 and hs_settled[0]["payload"]["start_cluster"]["rank"] == 0
    assert rep["skips"].get("HS_deferred:same_start") == sc["second_round"]
    assert "HS_skipped:refine_budget" not in rep["skips"]
    assert external == rep["ledger_total"] == paid(rep)


def _scan_record(q, energy, rf, edges, shift=0.0):
    positions = [[0.0, 0.0, 0.0], [q, 0.0, 0.0], [0.0, 2.0, shift]]
    return {
        "delta_A": q - 1.5,
        "q_target": q,
        "q_actual": q,
        "energy_eV": energy,
        "reaction_force": rf,
        "fmax_projected": 0.0,
        "converged": True,
        "steps": 1,
        "edges": edges,
        "positions": positions,
        "forces_pes": [[0.0] * 3] * 3,
    }


@pytest.mark.parametrize("direction", [+1, -1])
@pytest.mark.parametrize("reverse_offset", [0.0, 2.22e-16])
def test_bracketed_jump_hands_the_forward_jump_interval_to_hs(direction, reverse_offset):
    """b05 HF v1.1 (2026-10-28): 11 branch_bracketed HS failed `bracket must have positive
    width` because the L2 reverse point -- at the reactant point's q, or 2.2e-16 off --
    was used as the product end. The HS bracket is the forward scan's jump interval."""
    from types import SimpleNamespace

    from explore.d2 import _classify, _hf_bracket, _points, _task

    qs = [1.5 + direction * 0.1 * k for k in range(4)]
    reactant, product = [[0, 1]], [[0, 1], [1, 2]]
    l1 = [
        _scan_record(qs[0], 0.0, 1.0, reactant),
        _scan_record(qs[1], 0.2, 1.0, reactant),
        # same reaction-force sign across the jump: a jump, not a ridge
        _scan_record(qs[2], -0.5, 1.0, product, shift=0.4),
        _scan_record(qs[3], -0.6, 1.0, product, shift=0.4),
    ]
    back = [
        _scan_record(qs[2], -0.5, -1.0, product, shift=0.4),
        _scan_record(qs[1] + reverse_offset, -0.3, -1.0, product, shift=0.4),
    ]
    row = {"intent_id": "p0|x", "audit": False}
    run = SimpleNamespace(
        settled={
            _task("L1", row): {"status": "completed", "payload": {"points": l1}},
            _task("L2", row): {"status": "completed", "payload": {"points": back}},
        },
        settings=SimpleNamespace(noise_eV=1e-3),
    )
    classified = _classify(run, row)
    assert classified["class"] == "branch_bracketed"
    points = _points(run.settled[_task("L1", row)]["payload"])
    a, b = _hf_bracket(classified, points)
    assert (a.q_target, b.q_target) == (qs[1], qs[2])  # reactant end, product end
    assert a.edges == reactant and b.edges == product
    low, high = sorted((a, b), key=lambda p: p.q_target)
    assert high.q_target - low.q_target > 0.05  # one scan step, never the L2 point


# Same-graph ridges last (D2-HF-v1.2; docs/PLAN_HF_SAME_GRAPH_RIDGE_2026-11-07.md).
def test_same_graph_ridge_reads_the_scan_peak_bonds():
    from explore.d2 import _points, _same_graph_ridge

    graph, other = [[0, 1]], [[0, 1], [1, 2]]

    def points(after):
        return _points(
            {
                "points": [
                    _scan_record(1.5, 0.0, 1.0, graph),
                    _scan_record(1.6, 0.3, 0.1, graph),
                    _scan_record(1.7, 0.1, -1.0, after, shift=0.4),
                ]
            }
        )

    ridge = {"class": "ridge_suspected", "info": {"peak": 1}}
    assert _same_graph_ridge(ridge, points(graph))
    assert not _same_graph_ridge(ridge, points(other))
    jump = {"class": "branch_bracketed", "info": {"peak": 1}}
    assert not _same_graph_ridge(jump, points(graph))


def test_same_graph_ridge_settings_and_freeze(tmp_path):
    from explore.d2 import StartCluster

    with pytest.raises(ValueError):
        init(tmp_path / "no_sc", hf=HF, same_graph_ridge_last=True)
    with pytest.raises(ValueError):
        init(
            tmp_path / "bad", hf=HF, start_cluster=StartCluster(2.0, 1), same_graph_ridge_last=1
        )
    init(tmp_path / "sc", hf=HF, start_cluster=StartCluster(2.0, 1))
    assert "same_graph_ridge_last" not in json.loads(
        (tmp_path / "sc" / "manifest.json").read_text()
    )
    init(tmp_path / "sg", hf=HF, start_cluster=StartCluster(2.0, 1), same_graph_ridge_last=True)
    manifest = json.loads((tmp_path / "sg" / "manifest.json").read_text())
    assert manifest["rule_version"] == "D2-HF-v1.2"
    assert manifest["same_graph_ridge_last"] is True


def _first_promoted_is_same_graph(monkeypatch):
    """Mark the first promoted item (cluster rank 0 under D2-HF-v1.1) as a same-graph ridge."""
    import explore.d2 as d2

    calls = []

    def fake(classified, points):
        calls.append(classified)
        return len(calls) == 1

    monkeypatch.setattr(d2, "_same_graph_ridge", fake)
    return calls


def _sg_run(tmp_path, **over):
    from explore.d2 import StartCluster

    run_dir = tmp_path / "run"
    init(run_dir, hf=HF, start_cluster=StartCluster(1e4, 1), same_graph_ridge_last=True, **over)
    tick = Tick()
    rep = run(run_dir, backends(tick))
    return run_dir, rep, tick.count


def test_same_graph_ridge_runs_last_with_ample_budget(tmp_path, monkeypatch, reference):
    _first_promoted_is_same_graph(monkeypatch)
    run_dir, rep, external = _sg_run(tmp_path)
    ref = report(reference[0])
    assert rep["rule_version"] == "D2-HF-v1.2"
    assert rep["hs_status"] == ref["hs_status"]  # order changed, nothing dropped
    hs = [e for e in events(run_dir) if e["event"] == "settle" and e["task"].startswith("HS|")]
    flags = [e["payload"]["same_graph_ridge"] for e in hs]
    assert flags == [False] * (len(hs) - 1) + [True]  # dispatched after both rounds
    # The same-graph item no longer leads the start cluster.
    assert hs[0]["payload"]["start_cluster"]["rank"] == 0
    assert rep["same_graph_ridge"] == {"settled": 1, "deferred": 0}
    assert external == rep["ledger_total"] == paid(rep)


def test_same_graph_ridge_is_deferred_when_the_budget_is_spent(tmp_path, monkeypatch):
    _first_promoted_is_same_graph(monkeypatch)
    run_dir, rep, external = _sg_run(tmp_path, certify_budget=HF.hs_cap)
    hs = [e for e in events(run_dir) if e["event"] == "settle" and e["task"].startswith("HS|")]
    assert len(hs) == 1 and hs[0]["payload"]["same_graph_ridge"] is False
    assert rep["skips"].get("HS_deferred:same_graph_ridge") == 1
    assert "HS_deferred:same_start" not in rep["skips"]
    assert rep["same_graph_ridge"] == {"settled": 0, "deferred": 1}
    assert external == rep["ledger_total"] == paid(rep)


# HS bracket rule v7 in D2 (D2-HF-v1.3; docs/PLAN_HF_D2_V13_2026-11-07.md).
def test_hs_bracket_settings_and_freeze(tmp_path):
    from explore.d2 import HSBracket, StartCluster, _ef_settings

    for bad in ((0.0, 0.025), (2e-3, float("nan")), (True, 0.025)):
        with pytest.raises(ValueError):
            HSBracket(*bad)
    sc = StartCluster(2.0, 1)
    with pytest.raises(ValueError):  # needs D2-HF-v1.2
        init(tmp_path / "no_sg", hf=HF, start_cluster=sc, hs_bracket=HSBracket(2e-3, 0.025))
    with pytest.raises(ValueError):  # coarse threshold must exceed the midpoint one (1e-3)
        init(
            tmp_path / "loose",
            hf=HF,
            start_cluster=sc,
            same_graph_ridge_last=True,
            hs_bracket=HSBracket(5e-4, 0.025),
        )
    init(
        tmp_path / "v13",
        hf=HF,
        start_cluster=sc,
        same_graph_ridge_last=True,
        hs_bracket=HSBracket(2e-3, 0.025),
    )
    manifest = json.loads((tmp_path / "v13" / "manifest.json").read_text())
    assert manifest["rule_version"] == "D2-HF-v1.3"
    assert manifest["hs_bracket"] == {
        "midpoint_coarse_transverse_fmax_eV_A": 2e-3,
        "midpoint_sign_bound_eV_A": 0.025,
    }
    assert _ef_settings(HF) == _ef_settings(HF, None)  # earlier rules unchanged
    ef = _ef_settings(HF, HSBracket(2e-3, 0.025))
    assert ef.midpoint_fine_optimizer == "fire"
    assert ef.midpoint_coarse_candidate and ef.endpoint_candidate


def test_v13_run_searches_with_bracket_v7(tmp_path, reference):
    from explore.d2 import HSBracket, StartCluster

    run_dir = tmp_path / "run"
    init(
        run_dir,
        hf=HF,
        start_cluster=StartCluster(1e4, 1),
        same_graph_ridge_last=True,
        hs_bracket=HSBracket(2e-3, 0.025),
    )
    tick = Tick()
    rep = run(run_dir, backends(tick))
    ref = report(reference[0])
    assert rep["rule_version"] == "D2-HF-v1.3"
    assert rep["hs_status"] == ref["hs_status"]
    assert rep["hc_status"] == ref["hc_status"]
    hs = [e for e in events(run_dir) if e["event"] == "settle" and e["task"].startswith("HS|")]
    assert hs and all(
        e["payload"]["refine"]["provenance"]["algorithm"] == "HF-EF-bracket-v7" for e in hs
    )
    assert tick.count == rep["ledger_total"] == paid(rep)


# Same-graph ridge sample + parallel HS / HC (D2-HF-v1.4; docs/PLAN_HF_D2_V14_2026-11-07.md).
def _v13(run_dir, **over):
    from explore.d2 import HSBracket, StartCluster

    init(
        run_dir,
        hf=HF,
        start_cluster=StartCluster(1e4, 1),
        same_graph_ridge_last=True,
        hs_bracket=HSBracket(2e-3, 0.025),
        **over,
    )


def _settles(run_dir, prefix):
    return [
        e for e in events(run_dir) if e["event"] == "settle" and e["task"].startswith(prefix)
    ]


def _comparable(ev):
    """Per-task outcome of a settle event, independent of journal position."""
    p = ev["payload"]
    out = {"status": ev["status"], "spent": ev["spent"], "cap": ev["cap"]}
    if ev["task"].startswith("HS|") and p.get("refine") and p["refine"].get("positions"):
        out["positions"] = p["refine"]["positions"]
    if ev["task"].startswith("HC|"):
        out["pair"] = p.get("chemical_pair")
        out["duplicate_of"] = p.get("duplicate_of")
    return out


def test_v14_settings_and_freeze(tmp_path):
    from explore.d2 import HSBracket, StartCluster

    sc, br = StartCluster(2.0, 1), HSBracket(2e-3, 0.025)
    with pytest.raises(ValueError):  # the two v1.4 options go together
        init(
            tmp_path / "a",
            hf=HF,
            start_cluster=sc,
            same_graph_ridge_last=True,
            hs_bracket=br,
            same_graph_sample=True,
        )
    with pytest.raises(ValueError):  # needs v1.3
        init(
            tmp_path / "b",
            hf=HF,
            start_cluster=sc,
            same_graph_ridge_last=True,
            same_graph_sample=True,
            hf_parallel=True,
        )
    with pytest.raises(ValueError):
        _v13(tmp_path / "c", same_graph_sample=1, hf_parallel=1)
    _v13(tmp_path / "v13")
    assert "hf_parallel" not in json.loads((tmp_path / "v13" / "manifest.json").read_text())
    _v13(tmp_path / "v14", same_graph_sample=True, hf_parallel=True)
    manifest = json.loads((tmp_path / "v14" / "manifest.json").read_text())
    assert manifest["rule_version"] == "D2-HF-v1.4"
    assert manifest["same_graph_sample"] is True and manifest["hf_parallel"] is True


@pytest.mark.parametrize("max_active, wave_size", [(1, 1), (2, 4)])
def test_v14_parallel_waves_match_the_serial_rule(tmp_path, max_active, wave_size):
    """Exact bench evaluator: every HS / HC outcome equals the serial v1.3 one, whatever
    the wave size; with waves of one the journal order is the serial order too. Both
    bench ridges are same-graph, so the sample fraction is 1 (all of them run)."""
    old, new = tmp_path / "v13", tmp_path / "v14"
    _v13(old, max_active=max_active, wave_size=wave_size, audit_fraction=1.0)
    _v13(
        new,
        max_active=max_active,
        wave_size=wave_size,
        audit_fraction=1.0,
        same_graph_sample=True,
        hf_parallel=True,
    )
    a, b = Tick(), Tick()
    rep_old, rep_new = run(old, backends(a)), run(new, backends(b))
    assert rep_new["rule_version"] == "D2-HF-v1.4"
    for prefix in ("HS|", "HC|"):
        olds = {e["task"]: _comparable(e) for e in _settles(old, prefix)}
        news = {e["task"]: _comparable(e) for e in _settles(new, prefix)}
        assert olds and news == olds
    dispatched = [
        e["task"] for e in events(new) if e["event"] == "dispatch" and e["task"].startswith("H")
    ]
    assert [e["task"] for e in _settles(new, "H")] == dispatched  # settled in input order
    if wave_size == 1:
        order = [e["task"] for e in _settles(old, "H")]
        assert [e["task"] for e in _settles(new, "H")] == order
    assert rep_new["hs_status"] == rep_old["hs_status"]
    assert rep_new["hc_status"] == rep_old["hc_status"]
    assert b.count == rep_new["ledger_total"] == paid(rep_new) == paid(rep_old)


@pytest.mark.parametrize("fraction, sampled", [(0.0, False), (1.0, True)])
def test_v14_same_graph_sample(tmp_path, monkeypatch, fraction, sampled):
    _first_promoted_is_same_graph(monkeypatch)
    run_dir = tmp_path / "run"
    _v13(run_dir, audit_fraction=fraction, same_graph_sample=True, hf_parallel=True)
    tick = Tick()
    rep = run(run_dir, backends(tick))
    hs = _settles(run_dir, "HS|")
    marked = [e for e in hs if e["payload"]["same_graph_ridge"]]
    sg = rep["same_graph_ridge"]
    if sampled:
        assert len(marked) == 1 and marked[0]["payload"]["same_graph_sample"] is True
        assert hs[-1] is marked[0]  # still the last group
        assert sg["skipped"] == 0 and sg["settled"] == 1
    else:
        assert not marked and sg["skipped"] == 1 and sg["settled"] == 0
        assert rep["skips"]["HS_skipped:same_graph_ridge"] == 1
        skipped = [e["task"] for e in events(run_dir) if e["event"] == "skip"]
        hc_of_skipped = [t.replace("HS|", "HC|") for t in skipped if t.startswith("HS|")]
        assert not set(hc_of_skipped) & {e["task"] for e in _settles(run_dir, "HC|")}
    assert sg["sample_candidates"] == [
        e["task"] for e in marked if e["status"] == "converged_candidate"
    ]
    assert tick.count == rep["ledger_total"] == paid(rep)


def test_v14_group_order_holds_when_the_budget_runs_out(tmp_path, monkeypatch):
    """Budget for one HS: round one gets it; the sampled same-graph ridge is deferred."""
    _first_promoted_is_same_graph(monkeypatch)
    run_dir = tmp_path / "run"
    _v13(
        run_dir,
        audit_fraction=1.0,
        certify_budget=HF.hs_cap,
        same_graph_sample=True,
        hf_parallel=True,
    )
    tick = Tick()
    rep = run(run_dir, backends(tick))
    hs = _settles(run_dir, "HS|")
    assert hs and not any(e["payload"]["same_graph_ridge"] for e in hs)
    assert rep["skips"].get("HS_deferred:same_graph_ridge") == 1
    assert tick.count == rep["ledger_total"] == paid(rep)


@pytest.fixture(scope="module")
def reference_v14(tmp_path_factory):
    base = tmp_path_factory.mktemp("v14ref")
    done = drive(base / "run", base / "count", -1, "v14")
    assert done.returncode == 0, done.stderr
    return base / "run", (base / "count").stat().st_size


@pytest.mark.parametrize("where", ["HS", "HC"])
def test_v14_kill_inside_a_wave_and_resume_pays_once(reference_v14, tmp_path, where):
    assert report(reference_v14[0])["rule_version"] == "D2-HF-v1.4"
    _kill_and_resume(reference_v14, tmp_path, where, "v14")


# Three general rules (D2-HF-v1.5; docs/PLAN_HF_D2_V15_2026-11-13.md).
def test_same_graph_basin_compares_the_basins_around_the_ridge():
    """SN2-like: the bond change is spread over steps, peak and next point share a graph
    (old test: same-graph), but the descent ends in another graph (basin test: changed)."""
    from explore.d2 import _points, _same_graph_ridge

    a, mid, b = [[0, 1]], [[0, 2]], [[1, 2]]

    def pts(graphs, energies):
        return _points(
            {
                "points": [
                    _scan_record(1.5 + 0.1 * i, e, 1.0 - i, g)
                    for i, (g, e) in enumerate(zip(graphs, energies))
                ]
            }
        )

    ridge = {"class": "ridge_suspected", "info": {"peak": 2}}
    spread = pts([a, a, mid, mid, b], [0.0, 0.2, 0.4, 0.3, -0.5])
    assert _same_graph_ridge(ridge, spread)
    assert not _same_graph_ridge(ridge, spread, basin=True)
    conformational = pts([a, a, a, mid, a], [0.0, 0.2, 0.4, 0.3, -0.5])
    assert not _same_graph_ridge(ridge, conformational)  # neighbours differ
    assert _same_graph_ridge(ridge, conformational, basin=True)  # both basins in graph a
    # the lowest point of the converged run after the peak decides, not the last point
    dip = pts([a, a, a, b, a], [0.0, 0.2, 0.4, -0.6, -0.1])
    assert not _same_graph_ridge(ridge, dip, basin=True)


def test_geometry_noise_measure_validation():
    from d2_bench import gn_measure

    assert gn_measure().cap() == 2 * 2 * 3 * 60
    for bad in (
        dict(ridges=0),
        dict(repeats=1),
        dict(max_steps=1.5),
        dict(sigma_A=0.0),
        dict(fallback_A=float("nan")),
        dict(taus=()),
        dict(taus=(1e-3, -1e-3)),
    ):
        with pytest.raises(ValueError):
            gn_measure(**bad)


def test_v15_settings_and_freeze(tmp_path):
    from d2_bench import gn_measure

    with pytest.raises(ValueError):  # needs v1.4
        _v13(tmp_path / "a", geometry_noise_measure=gn_measure())
    with pytest.raises(ValueError):
        _v13(
            tmp_path / "b",
            same_graph_sample=True,
            hf_parallel=True,
            geometry_noise_measure={"ridges": 2},
        )
    _v13(tmp_path / "v14", same_graph_sample=True, hf_parallel=True)
    old = json.loads((tmp_path / "v14" / "manifest.json").read_text())
    assert old["rule_version"] == "D2-HF-v1.4" and "same_graph_basin" not in old
    _v13(
        tmp_path / "v15",
        same_graph_sample=True,
        hf_parallel=True,
        geometry_noise_measure=gn_measure(),
    )
    m = json.loads((tmp_path / "v15" / "manifest.json").read_text())
    assert m["rule_version"] == "D2-HF-v1.5"
    assert m["same_graph_basin"] is True and m["same_graph_fallback"] is True
    assert m["geometry_noise_measure"]["taus"] == [3e-3, 1e-3]


def _v15(run_dir, **over):
    from d2_bench import gn_measure

    _v13(
        run_dir,
        same_graph_sample=True,
        hf_parallel=True,
        geometry_noise_measure=gn_measure(),
        **over,
    )


def test_v15_measures_geometry_noise_on_bond_changing_ridges(tmp_path, monkeypatch):
    import explore.d2 as d2

    monkeypatch.setattr(d2, "_same_graph_ridge", lambda c, p, basin=False: False)
    run_dir = tmp_path / "run"
    _v15(run_dir)
    tick = Tick()
    rep = run(run_dir, backends(tick))
    gn = _settles(run_dir, "GN|")
    assert len(gn) == 2  # both bench ridges, at most ``ridges`` = 2
    assert all(e["account"] == "refine" and e["status"] == "completed" for e in gn)
    for e in gn:
        assert [(r["where"], r["tau"]) for r in e["payload"]["points"]] == [
            ("before_peak", 3e-3),
            ("before_peak", 1e-3),
            ("peak", 3e-3),
            ("peak", 1e-3),
        ]
    g = rep["geometry_noise"]
    values = [
        r["max_pair_rmsd_A"]
        for e in gn
        for r in e["payload"]["points"]
        if r["max_pair_rmsd_A"] is not None and r["max_pair_rmsd_A"] <= 0.05
    ]
    assert g["source"] == "measured" and g["geometry_noise_A"] == max(values)
    assert rep["start_clusters"]["radius_A"] == pytest.approx(1e4 * max(values))
    # GN tasks come before every HS, HS after all of them
    order = [e["task"][:2] for e in events(run_dir) if e["event"] == "settle"]
    assert order.index("HS") > max(i for i, t in enumerate(order) if t == "GN")
    assert tick.count == rep["ledger_total"] == paid(rep)
    # resume: nothing new, same value from the journal
    before = len(events(run_dir))
    again = run(run_dir, backends(Tick()))
    assert len(events(run_dir)) == before
    assert again["geometry_noise"] == rep["geometry_noise"]


def test_v15_geometry_noise_falls_back_without_budget(tmp_path, monkeypatch):
    import explore.d2 as d2
    from d2_bench import gn_measure

    monkeypatch.setattr(d2, "_same_graph_ridge", lambda c, p, basin=False: False)
    run_dir = tmp_path / "run"
    _v15(run_dir, certify_budget=gn_measure().cap() - 1)
    rep = run(run_dir, backends(Tick()))
    assert rep["skips"]["GN_skipped:refine_budget"] == 2
    assert rep["geometry_noise"]["source"] == "fallback"
    assert rep["geometry_noise"]["geometry_noise_A"] == 0.0104


def test_same_graph_fallback_reads_the_sampled_hc():
    from types import SimpleNamespace

    from explore.d2 import _same_graph_fallback

    rows = [({"intent_id": "a"},), ({"intent_id": "b"},)]
    for status, expect in (
        ("connected", True),
        ("hf_duplicate", True),
        ("same_basin", False),
        ("failed:tangent_no_resolved_drop", False),
    ):
        run_ = SimpleNamespace(settled={"HC|b": {"status": status}})
        assert _same_graph_fallback(run_, rows) is expect
    assert _same_graph_fallback(SimpleNamespace(settled={}), rows) is False


@pytest.mark.parametrize("triggered", [True, False])
def test_v15_unsampled_same_graph_ridges_wait_for_the_sample(tmp_path, monkeypatch, triggered):
    """Both bench ridges are same-graph; none sampled (fraction 0). Fallback on: they run
    last, tagged; off: skipped as in v1.4, and no HC for them."""
    import explore.d2 as d2

    monkeypatch.setattr(d2, "_same_graph_fallback", lambda run, sampled: triggered)
    run_dir = tmp_path / "run"
    _v15(run_dir, audit_fraction=0.0)
    tick = Tick()
    rep = run(run_dir, backends(tick))
    hs = _settles(run_dir, "HS|")
    sg = rep["same_graph_ridge"]
    if triggered:
        assert len(hs) == 2 and all(e["payload"]["same_graph_fallback"] for e in hs)
        assert all(e["payload"]["same_graph_sample"] is False for e in hs)
        assert sg["fallback"] == 2 and sg["skipped"] == 0
    else:
        assert not hs and sg["skipped"] == 2 and sg["fallback"] == 0
    assert tick.count == rep["ledger_total"] == paid(rep)


BOND_CHANGING = {"D2_BENCH_BOND_CHANGING": "1"}


@pytest.fixture(scope="module")
def reference_v15(tmp_path_factory):
    base = tmp_path_factory.mktemp("v15ref")
    done = drive(base / "run", base / "count", -1, "v15", BOND_CHANGING)
    assert done.returncode == 0, done.stderr
    return base / "run", (base / "count").stat().st_size


@pytest.mark.parametrize("where", ["GN", "HS"])
def test_v15_kill_and_resume_pays_once(reference_v15, tmp_path, where):
    rep = report(reference_v15[0])
    assert rep["rule_version"] == "D2-HF-v1.5"
    assert rep["geometry_noise"]["source"] == "measured"
    _kill_and_resume(reference_v15, tmp_path, where, "v15", BOND_CHANGING)
    assert report(tmp_path / "run")["geometry_noise"] == rep["geometry_noise"]
