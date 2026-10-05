from __future__ import annotations

import pytest

from pars.adapters.prrs_run import PrrsRun
from pars.config import ParsConfig
from pars.paths.ingest import ingest_run


def test_ingest_selects_completed_reactive_attempts(run_dir):
    traces = ingest_run(run_dir, ParsConfig())
    assert [t.trace_id for t in traces] == ["t000001", "t000002"]


def test_trace_carries_full_trajectory_and_phases(run_dir):
    trace = ingest_run(run_dir, ParsConfig())[0]
    assert len(trace.phases) == 8
    assert trace.phases[0] == "source"
    assert trace.phases[-1] == "quenched"
    assert trace.positions.shape == (8, 3, 3)
    assert trace.energies.shape == (8,)
    assert trace.signatures is not None
    assert len(trace.signatures) == 8


def test_trace_uncertainty_availability_mask(run_dir, run_dir_with_uncertainty):
    without = ingest_run(run_dir, ParsConfig())[1]
    assert without.uncertainty is None  # run had no committee

    with_unc = ingest_run(run_dir_with_uncertainty, ParsConfig())[1]
    assert with_unc.uncertainty is not None
    assert with_unc.uncertainty.shape == (8,)
    assert with_unc.uncertainty[0] == pytest.approx(0.01)


def test_provenance_is_complete(run_dir):
    run = PrrsRun(run_dir)
    trace = run.load_trace("t000001")
    provenance = trace.provenance
    assert provenance.channel_id == "rc0000"
    assert provenance.event_key == "aaaabbbbccccdddd"
    assert provenance.reactant_state_id == "c0000"
    assert provenance.product_state_id == "c0001"
    assert provenance.perturbation_family == "stretch"
    assert provenance.perturbation_amplitude == 0.6
    assert provenance.seed == 1235
    assert provenance.implementation_sha256 == "0" * 64
    assert provenance.frame_start == 0
    assert provenance.frame_stop == 8
    assert provenance.extra["classification"] == "reactive"


def test_old_format_run_rejected(tmp_path):
    import json

    from pars.config import ParsConfig

    run = tmp_path / "old"
    run.mkdir()
    (run / "network.json").write_text(json.dumps({"schema_version": 1}))
    (run / "config.json").write_text(json.dumps({"units": {}}))
    (run / "attempts.jsonl").write_text("")
    with pytest.raises(ValueError, match="reaction_channels"):
        ingest_run(run, ParsConfig())


def test_frame_row_mismatch_rejected(run_dir):

    run = PrrsRun(run_dir)
    observations = run.run_dir / "trials" / "t000001" / "observations.jsonl"
    rows = observations.read_text(encoding="utf-8").splitlines()
    observations.write_text("\n".join(rows[:-1]) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="refusing to guess"):
        run.load_trace("t000001")


def test_attempt_in_several_reactions_rejected(run_dir):
    import json

    network = json.loads((run_dir / "network.json").read_text())
    network["reactions"][1]["attempts"].append("t000001")
    (run_dir / "network.json").write_text(json.dumps(network))
    with pytest.raises(ValueError, match="multiple reactions"):
        PrrsRun(run_dir).load_trace("t000001")


@pytest.mark.parametrize("mutation", ["reverse", "duplicate", "missing", "boolean", "phase"])
def test_observation_identity_mismatch_is_rejected(run_dir, mutation):
    import json

    path = run_dir / "trials/t000001/observations.jsonl"
    rows = [json.loads(line) for line in path.read_text().splitlines()]
    if mutation == "reverse":
        rows.reverse()
    elif mutation == "duplicate":
        rows[1]["frame"] = 0
    elif mutation == "missing":
        del rows[0]["frame"]
    elif mutation == "boolean":
        rows[0]["frame"] = False
    else:
        rows[0]["phase"] = "quenched"
    path.write_text("".join(json.dumps(row) + "\n" for row in rows))
    with pytest.raises(ValueError, match="index mismatch|phase mismatch"):
        PrrsRun(run_dir).load_trace("t000001")


@pytest.mark.parametrize("artifact", ["trajectory", "endpoint"])
def test_import_rejects_species_sequence_changes(run_dir, artifact):
    from ase.io import read, write

    trial = run_dir / "trials/t000001"
    if artifact == "trajectory":
        path = trial / "response.traj"
        frames = read(str(path), index=":")
        frames[2].numbers[:] = frames[2].numbers[::-1]
    else:
        path = trial / "endpoint.extxyz"
        frames = [read(str(path))]
        frames[0].numbers[:] = frames[0].numbers[::-1]
    write(str(path), frames)
    with pytest.raises(ValueError, match="atomic mapping mismatch"):
        PrrsRun(run_dir).load_trace("t000001")


def test_ingest_rejects_disagreement_with_authoritative_result(run_dir):
    import json

    path = run_dir / "trials/t000001/result.json"
    result = json.loads(path.read_text())
    result["status"] = "failed"
    path.write_text(json.dumps(result))
    with pytest.raises(ValueError, match="ledger mismatch in status"):
        ingest_run(run_dir, ParsConfig())
