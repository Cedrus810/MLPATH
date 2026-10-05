"""Finite-horizon windows must not silently become static boundary claims."""

import json

import pytest

from prrs.calculators import demo_atoms, double_well_factory
from prrs.config import SearchConfig
from prrs.search import ReactionSearch
from prrs.window_audit import assess_horizons, replay_window_audit, response_windows


def test_matching_windows_do_not_hide_changed_product_identity():
    baseline = {"0.1": "returned", "0.2": "departed:B"}
    changed = {"0.1": "returned", "0.2": "departed:C"}
    assert response_windows(baseline) == response_windows(changed)
    result = assess_horizons({"10": {"labels": baseline}, "20": {"labels": changed}}, [10, 20])
    assert result["status"] == "sensitive"
    assert not result["static_basin_boundary"]


@pytest.mark.parametrize("other", [{}, {"0.1": "unresolved:quench_failed"}])
def test_missing_and_unresolved_never_certify_stability(other):
    result = assess_horizons({"10": {"labels": other}, "20": {"labels": other}}, [10, 20])
    assert not result["stable_on_tested_grid"]


def test_one_horizon_cannot_certify_itself():
    with pytest.raises(ValueError):
        assess_horizons({}, [10, 10])


@pytest.fixture
def source_run(tmp_path):
    config = SearchConfig(
        families=("stretch",),
        geometry_amplitudes_A=(0.1, 0.3),
        response_steps=10,
        timestep_fs=0.1,
        max_trials=2,
        max_depth=0,
        closed_shell_only=False,
        refinement_steps=0,
        quench_fmax_eV_A=0.001,
    )
    path = tmp_path / "run"
    network = ReactionSearch(double_well_factory, config).run(demo_atoms(), path)
    assert network["window_semantics"]["response_time_fs"] == 1.0
    assert all(
        s["window_semantics"]["static_basin_boundary"] is False
        for s in network["amplitude_scans"]
    )
    return path


def test_replay_records_provenance_and_leaves_original_unchanged(source_run, tmp_path):
    before = {
        p.relative_to(source_run): p.read_bytes() for p in source_run.rglob("*") if p.is_file()
    }
    report = replay_window_audit(
        source_run,
        tmp_path / "audit",
        double_well_factory,
        {"factory": "python_callable"},
        [20, 30],
    )
    assert report["baseline_reproduced"]
    assert report["source_unchanged"]
    assert report["horizons_steps"] == [10, 20, 30]
    assert report["status"] == "stable_on_tested_grid"
    assert before == {
        p.relative_to(source_run): p.read_bytes() for p in source_run.rglob("*") if p.is_file()
    }


def test_baseline_mismatch_runs_no_alternate_horizon(source_run, tmp_path):
    network = json.loads((source_run / "network.json").read_text())
    sample = network["amplitude_scans"][0]["samples"][0]
    records = [
        json.loads(line) for line in (source_run / "attempts.jsonl").read_text().splitlines()
    ]
    record = next(r for r in records if r["attempt_id"] == sample["attempt"])
    from ase.io import read, write

    endpoint = source_run / record["endpoint"]
    atoms = read(endpoint)
    atoms.positions[0, 0] += 0.2
    write(endpoint, atoms)
    output = tmp_path / "audit"
    with pytest.raises(ValueError, match="baseline endpoint mismatch"):
        replay_window_audit(
            source_run, output, double_well_factory, {"factory": "python_callable"}, [20]
        )
    report = json.loads((output / "report.json").read_text())
    assert not report["baseline_reproduced"]
    assert report["status"] == "failed"
    assert not list((output / "trials").glob("*h000020*"))
