"""Two record-layer fixes of 2026-09-26 (PLAN_2026-09-25 D1 and C1).

D1: a saddle's two-sided descent wrote both branches into one descent.extxyz through
two recorders that each counted frames from 0, so the second branch's first write
erased the first's geometry. C1: every scan sample now carries an outcome class and,
when unresolved, why -- recorded only, never routed on.
"""

import json

from ase.io import read

from prrs.calculators import demo_atoms, double_well_factory
from prrs.config import SearchConfig
from prrs.search import ReactionSearch, _sample_outcome

OUTCOMES = {"returned", "departed", "ridge", "transient", "unresolved"}


def _config():
    # the same small double-well search tests/test_resume.py kills and resumes; it
    # climbs to a saddle and descends both sides of it
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


def test_a_saddle_descent_keeps_the_geometry_of_both_sides(tmp_path):
    ReactionSearch(double_well_factory, _config()).run(demo_atoms(), tmp_path / "run")
    descents = sorted((tmp_path / "run" / "paths").glob("ts*"))
    assert descents, "the fixture must descend from a saddle; otherwise this checks nothing"
    for directory in descents:
        frames = read(str(directory / "descent.extxyz"), index=":")
        rows = [json.loads(line) for line in (directory / "descent.jsonl").open()]
        assert len(frames) == len(rows), "every logged row has its geometry"
        assert {f.info["phase"] for f in frames} == {"descent_plus", "descent_minus"}
        for row in rows:
            assert frames[row["frame"]].info["phase"] == row["phase"]
            assert row["phase"] == ("descent_plus" if row["sign"] > 0 else "descent_minus")


def test_each_outcome_class_and_its_reason():
    source = "c0000/m0000"
    completed = {"status": "completed", "admission": "known_microstate"}
    assert _sample_outcome({**completed, "target_microstate": source}, source) == (
        "returned",
        None,
    )
    assert _sample_outcome({**completed, "target_microstate": "c0001/m0000"}, source) == (
        "departed",
        None,
    )
    assert _sample_outcome(
        {"status": "completed", "target_microstate": None, "admission": "rejected"}, source
    ) == ("unresolved", "admission:rejected")
    assert _sample_outcome({"status": "ts_candidate"}, source) == ("ridge", None)
    assert _sample_outcome(
        {"status": "unconfirmed", "crossing": {"observed": True}}, source
    ) == ("transient", None)
    assert _sample_outcome(
        {"status": "quench_failed", "failure": {"code": "quench_nonconverged"}}, source
    ) == ("unresolved", "quench_failed:quench_nonconverged")
    assert _sample_outcome(
        {"status": "rejected", "failure": {"code": "nve_drift"}}, source
    ) == (
        "unresolved",
        "rejected:nve_drift",
    )


def test_scan_samples_carry_a_consistent_outcome(tmp_path):
    network = ReactionSearch(double_well_factory, _config()).run(demo_atoms(), tmp_path / "run")
    samples = 0
    for scan in network["amplitude_scans"]:
        for sample in scan["samples"]:
            samples += 1
            assert sample["outcome"] in OUTCOMES
            if sample["outcome"] == "returned":
                assert sample["target_microstate"] == scan["source_microstate"]
            if sample["outcome"] == "departed":
                assert sample["target_microstate"] not in (None, scan["source_microstate"])
            if sample["outcome"] == "unresolved":
                assert sample["outcome_reason"]
            else:
                assert sample["outcome_reason"] is None
    assert samples, "the fixture must sample; otherwise this checks nothing"
