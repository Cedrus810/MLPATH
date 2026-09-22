"""The schema anchor: SCHEMA_2 is a description, and runs/ is the thing described.

``src/prrs/schema.py`` states what the current writer puts in network.json. This test
holds that statement against every network.json already on disk:

* files written by older writers fail ``validate`` in recorded, era-explainable ways --
  that inventory is the DRIFT_LEDGER below, pinned so it can never grow silently;
* every other file -- i.e. anything written by the current writer family, starting
  with the b05 probes and ``b04_frozen1`` -- must validate clean.

Both directions are load-bearing. If ``validate`` were vacuous, the legacy files would
stop matching their recorded patterns and the ledger assertion would fail; if the
writer dropped or renamed a field, the clean set would fail. A green that cannot die
is not evidence, and this one can.
"""

import re
from pathlib import Path

import pytest

from prrs.calculators import demo_atoms, double_well_factory
from prrs.config import SearchConfig
from prrs.schema import read_network, validate
from prrs.search import ReactionSearch

RUNS = Path(__file__).resolve().parent.parent / "runs"

# Recorded 2026-09-20 against implementation c393da782a1350a1. Each entry is
# (files, deviation patterns) -- patterns with entry indices collapsed to [*].
# Every deviation below is an artefact of the writer era that produced the file
# (fields introduced by 5dc83ed4/fa15a829/89f41c80 and later), not an inconsistency:
# the older the run, the more fields it predates.
DRIFT_LEDGER = {
    "b04_probe/network.json knob_T300_s17/network.json knob_T300_s19/network.json "
    "p1_round7/network.json p1_round8/network.json p2_T300_seed17/network.json "
    "p2_T300_seed19/network.json p2_T300_seed23/network.json p2_T300_seed29/network.json "
    "p2_T300_seed31/network.json p2_T300_seed37/network.json p2_T300_seed41/network.json "
    "p2_T300_seed43/network.json p2_T300_seed47/network.json p2_T300_seed53/network.json "
    "p2_T300_seed59/network.json p2_T300_seed61/network.json p2_T300_seed67/network.json "
    "p2_T300_seed71/network.json p2_T300_seed73/network.json p2_T300_seed79/network.json "
    "p2_round4/network.json p2_round5/network.json p2_round6/network.json "
    "p2_round7/network.json p2_seed19/network.json p2_seed23/network.json "
    "p2_seed29/network.json p2_seed31/network.json p2_seed37/network.json "
    "p2_seed41/network.json p2_seed43/network.json": (
        "top: 'withheld_candidates' is not a list",
        "top: completed run missing 'direction_ledger'",
        "top: completed run missing 'direction_ledger_per_microstate'",
        "top: completed run missing 'termination_detail'",
    ),
    # The T=300 scan's killed-mid-run seeds: status running, so only the pre-5dc83ed4
    # gap shows.
    "p2_seed101/network.json p2_seed103/network.json p2_seed107/network.json "
    "p2_seed109/network.json p2_seed113/network.json p2_seed127/network.json "
    "p2_seed47/network.json p2_seed53/network.json p2_seed59/network.json "
    "p2_seed61/network.json p2_seed67/network.json p2_seed71/network.json "
    "p2_seed73/network.json p2_seed79/network.json p2_seed83/network.json "
    "p2_seed89/network.json p2_seed97/network.json": (
        "top: 'withheld_candidates' is not a list",
    ),
    "b04_probe2_openceiling/network.json b04_probe3_rigidbond/network.json": (
        "top: completed run missing 'direction_ledger'",
        "top: completed run missing 'direction_ledger_per_microstate'",
        "top: completed run missing 'termination_detail'",
    ),
    "demo/network.json smoke_ethanol/network.json": ("top: schema_version 1 != 2",),
    "b04_probe4_ledger/network.json": (
        "direction_ledger: missing microstates_searched",
        "top: completed run missing 'termination_detail'",
    ),
    "p1_round4/network.json p1_round5/network.json p1_round6/network.json "
    "p2_oxobutanal/network.json": (
        "amplitude_scans[*]: missing proposal, boundary_continuation",
        "top: 'withheld_candidates' is not a list",
        "top: completed run missing 'direction_ledger'",
        "top: completed run missing 'direction_ledger_per_microstate'",
        "top: completed run missing 'termination_detail'",
    ),
    "e2e_ethanol/network.json": (
        "amplitude_scans[*]: missing proposal, saddle_search, boundary_continuation",
        "chemical_nodes[*].microstates[*]: missing free_bonds",
        "top: 'reaction_channels' is not a list",
        "top: 'withheld_candidates' is not a list",
        "top: completed run missing 'direction_ledger'",
        "top: completed run missing 'direction_ledger_per_microstate'",
        "top: completed run missing 'termination_detail'",
        "ts_candidates[*]: missing found_by, curvature_source, "
        "trivial_mode_floor_worst_residual, response",
    ),
    "p0_after_batch2/network.json": (
        "amplitude_scans[*]: missing proposal, saddle_search, boundary_continuation",
        "top: 'reaction_channels' is not a list",
        "top: 'withheld_candidates' is not a list",
        "top: completed run missing 'direction_ledger'",
        "top: completed run missing 'direction_ledger_per_microstate'",
        "top: completed run missing 'termination_detail'",
        "ts_candidates[*]: missing found_by",
    ),
    "p0_after_curvature_batch/network.json": (
        "amplitude_scans[*]: missing proposal, saddle_search, boundary_continuation",
        "chemical_nodes[*].microstates[*]: missing free_bonds",
        "top: 'reaction_channels' is not a list",
        "top: 'withheld_candidates' is not a list",
        "top: completed run missing 'direction_ledger'",
        "top: completed run missing 'direction_ledger_per_microstate'",
        "top: completed run missing 'termination_detail'",
        "ts_candidates[*]: missing found_by",
    ),
    "p1_malonaldehyde/network.json": (
        "amplitude_scans[*]: missing proposal, saddle_search, boundary_continuation",
        "reactions[*]: missing mechanisms, continuations",
        "top: 'withheld_candidates' is not a list",
        "top: completed run missing 'direction_ledger'",
        "top: completed run missing 'direction_ledger_per_microstate'",
        "top: completed run missing 'termination_detail'",
    ),
    "p1_round3/network.json": (
        "amplitude_scans[*]: missing proposal, saddle_search, boundary_continuation",
        "top: 'withheld_candidates' is not a list",
        "top: completed run missing 'direction_ledger'",
        "top: completed run missing 'direction_ledger_per_microstate'",
        "top: completed run missing 'termination_detail'",
    ),
    "p2_round2/network.json": (
        "amplitude_scans[*]: missing boundary_continuation",
        "reactions[*]: missing mechanisms, continuations",
        "top: 'withheld_candidates' is not a list",
        "top: completed run missing 'direction_ledger'",
        "top: completed run missing 'direction_ledger_per_microstate'",
        "top: completed run missing 'termination_detail'",
    ),
    "p2_round3/network.json": (
        "amplitude_scans[*]: missing boundary_continuation",
        "top: 'withheld_candidates' is not a list",
        "top: completed run missing 'direction_ledger'",
        "top: completed run missing 'direction_ledger_per_microstate'",
        "top: completed run missing 'termination_detail'",
    ),
}

EXPECTED = {
    file: patterns for files, patterns in DRIFT_LEDGER.items() for file in files.split()
}


def patterns_of(network):
    return tuple(sorted({re.sub(r"\[\w+\]", "[*]", p) for p in validate(network)}))


def test_anchor_every_network_on_disk_matches_the_ledger_or_is_clean():
    seen = 0
    for path in sorted(RUNS.glob("*/network.json")) + sorted(RUNS.glob("*/*/network.json")):
        seen += 1
        relative = path.relative_to(RUNS).as_posix()
        problems = patterns_of(read_network(path))
        if relative in EXPECTED:
            assert problems == EXPECTED[relative], relative
        else:
            # The forward gate: anything written by the current writer family --
            # b04_frozen1 and the b05 probes today, every future run tomorrow --
            # must match the schema with zero deviations.
            assert problems == (), (relative, problems)
    assert seen > 60, "the runs/ corpus shrank; the anchor lost its teeth"


def test_fresh_run_is_clean_and_carries_no_warning(tmp_path):
    # The same small double-well protocol the search tests use: it reaches
    # completion, so the completed-lifecycle rules are actually exercised.
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
    search = ReactionSearch(double_well_factory, config)
    network = search.run(demo_atoms(), tmp_path / "run")
    assert validate(network) == []
    on_disk = read_network(tmp_path / "run" / "network.json")
    assert validate(on_disk) == []
    # A clean network writes no schema_warnings key: byte-for-byte the file the
    # pre-schema writer would have produced.
    assert "schema_warnings" not in on_disk


def test_version_gate_reports_and_stops():
    problems = validate({"schema_version": 1, "status": "running"})
    assert problems == ["top: schema_version 1 != 2"]


def test_completed_lifecycle_rules_are_real_demands(tmp_path):
    config = SearchConfig(
        families=("stretch", "compress"),
        geometry_amplitudes_A=(0.1, 0.3),
        timestep_fs=0.1,
        response_steps=10,
        quench_fmax_eV_A=0.001,
        closed_shell_only=False,
        max_depth=1,
        max_trials=6,
    )
    network = ReactionSearch(double_well_factory, config).run(demo_atoms(), tmp_path / "run")
    assert validate(network) == []
    broken = dict(network)
    del broken["termination"], broken["direction_ledger"]
    problems = validate(broken)
    assert any("termination" in p for p in problems)
    assert any("direction_ledger" in p for p in problems)
    # The same network mid-run (before the final update) is legitimate: a partial
    # snapshot must not be asked for completion fields.
    broken["status"] = "running"
    problems = validate(broken)
    assert not any("termination" in p or "direction_ledger" in p for p in problems)


def test_read_network_matches_plain_json(tmp_path):
    path = tmp_path / "n.json"
    path.write_text('{"schema_version": 2, "status": "running"}', encoding="utf-8")
    assert read_network(path)["status"] == "running"


def test_ledger_covers_only_existing_files():
    for relative in EXPECTED:
        assert (RUNS / relative).exists(), relative


@pytest.mark.parametrize("missing", ["reservoir", "conformer_transitions", "microstates"])
def test_chemical_node_fields_are_required(missing):
    network = {
        "schema_version": 2,
        "status": "running",
        "interpretation": "x",
        "attempts": 0,
        "chemical_nodes": [
            {
                "id": "c0000",
                "chemical_key": "k",
                "depth": 0,
                "discovered_by": "t0",
                "fragments": ["C"],
                "key_components": {},
                "trials_spent": 0,
                "reservoir": {},
                "conformer_transitions": [],
                "microstates": [],
            }
        ],
        "reactions": [],
        "reaction_channels": [],
        "ts_candidates": [],
        "ts_connected_channels": [],
        "amplitude_scans": [],
        "withheld_candidates": [],
    }
    assert validate(network) == []
    del network["chemical_nodes"][0][missing]
    (problem,) = validate(network)
    assert missing in problem
