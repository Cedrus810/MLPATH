"""M1 schema acceptance: unknown parameters, NaN, missing identity, hash
conflicts, same-graph regions and frozen units are refused; statistics
carry published values unchanged (plan section 4, M1 completion conditions;
section 3.3 data objects).
"""

import dataclasses
import math

from dataclasses import replace

import pytest

from conftest import fake_sha, make_protocol, make_sides, make_source_run
from mature.schema import (
    MaturationProtocol,
    MicrostateRef,
    PathProposal,
    RegionDefinition,
    SchemaError,
    ShotRecord,
    ShootingEstimate,
    wilson_interval,
)


@pytest.fixture
def source(tmp_path):
    return make_source_run(tmp_path)


@pytest.fixture
def protocol(source):
    return make_protocol(source["root"], source["files"])


# ---- construction and round-trip -------------------------------------------


def test_protocol_round_trips_through_dict_with_identical_hash(protocol):
    revived = MaturationProtocol.from_dict(protocol.to_dict())
    assert revived == protocol
    assert revived.sha256 == protocol.sha256


def test_protocol_is_frozen(protocol):
    with pytest.raises(dataclasses.FrozenInstanceError):
        protocol.batch_role = "confirm"


def test_region_graph_hashes_differ_by_construction(protocol):
    hashes = protocol.regions.graph_sha256s
    assert hashes["A"] != hashes["B"]


# ---- the refusals (M1 acceptance) -------------------------------------------


def test_unknown_field_is_refused_at_every_level(protocol):
    data = protocol.to_dict()
    data["unexpected"] = 1
    with pytest.raises(SchemaError, match="unknown fields"):
        MaturationProtocol.from_dict(data)
    data = protocol.to_dict()
    data["dynamics"]["temperature"] = 300.0  # plausible-looking wrong name
    with pytest.raises(SchemaError, match="unknown fields"):
        MaturationProtocol.from_dict(data)


def test_nan_and_infinity_are_refused(source):
    with pytest.raises(SchemaError, match="not finite"):
        make_protocol(source["root"], source["files"]).dynamics.__class__(
            temperature_K=float("nan"),
            timestep_fs=0.25,
            first_entry_k=160,
            free_step_cap=1280,
        )
    with pytest.raises(SchemaError, match="not finite"):
        make_protocol(
            source["root"], source["files"], budget={"max_wall_minutes": float("inf")}
        )


def test_missing_identity_is_refused():
    sides = make_sides()
    broken = sides["A"].to_dict()
    del broken["microstate_id"]
    with pytest.raises(SchemaError, match="missing required fields"):
        RegionDefinition.from_dict({"sides": {"A": broken, "B": sides["B"].to_dict()}})


def test_malformed_hash_is_refused(protocol):
    data = protocol.to_dict()
    data["source"]["implementation_sha256"] = "ebff4ec0"  # a real-looking short form
    with pytest.raises(SchemaError, match="64-hex"):
        MaturationProtocol.from_dict(data)


def test_same_graph_region_is_refused_at_definition_time():
    shared = MicrostateRef(
        node_id="c0000",
        microstate_id="c0000/m0000",
        chemical_key=fake_sha("k"),
        structure_path="chemical/c0000_m0000.extxyz",
        structure_sha256=fake_sha("s"),
        edges=[[0, 1], [1, 2]],
    )
    twin = MicrostateRef(
        node_id="c0000",
        microstate_id="c0000/m0001",
        chemical_key=fake_sha("k"),
        structure_path="chemical/c0000_m0001.extxyz",
        structure_sha256=fake_sha("s2"),
        edges=[[0, 1], [1, 2]],  # same graph, different conformer: NOT a region
    )
    with pytest.raises(SchemaError, match="unsupported_same_graph"):
        RegionDefinition({"A": shared, "B": twin})


def test_t0_dynamics_is_refused(source):
    # mirrors prrs.runner.shoot's own refusal: deterministic dynamics has no
    # thermal ensemble to shoot
    with pytest.raises(SchemaError, match="must be > 0"):
        make_protocol(
            source["root"],
            source["files"],
            dynamics={
                "temperature_K": 0,
                "timestep_fs": 0.25,
                "first_entry_k": 160,
                "free_step_cap": 1280,
            },
        )


def test_pulse_steps_are_refused_for_v1(source):
    with pytest.raises(SchemaError, match="pulse"):
        make_protocol(
            source["root"],
            source["files"],
            dynamics={
                "temperature_K": 300.0,
                "timestep_fs": 0.25,
                "first_entry_k": 160,
                "free_step_cap": 1280,
                "pulse_steps": 20,
            },
        )


def test_open_drift_tolerance_is_a_zero_measurement_privilege(source):
    make_protocol(source["root"], source["files"], batch_role="zero_measure")
    with pytest.raises(SchemaError, match="zero-measurement"):
        make_protocol(source["root"], source["files"]).shooting.__class__(
            batch_role="explore",
            seed_base=1,
            shots_per_candidate=16,
            drift_tolerance_eV_atom=None,
        )


def test_confirm_and_sidecheck_need_a_frozen_selection_rule(source):
    with pytest.raises(SchemaError, match="selection_rule"):
        make_protocol(
            source["root"],
            source["files"],
            batch_role="confirm",
            shooting={
                "batch_role": "confirm",
                "seed_base": 1,
                "confirm_shots": 100,
                "drift_tolerance_eV_atom": 0.01,
            },
        )


def test_budget_with_no_limits_is_refused(source):
    with pytest.raises(SchemaError, match="at least one limit"):
        make_protocol(source["root"], source["files"], budget={})


def test_units_are_pinned(source):
    with pytest.raises(SchemaError, match="units"):
        make_protocol(
            source["root"],
            source["files"],
            units={"energy": "kcal/mol", "length": "Bohr", "time": "ps", "temperature": "K"},
        )


def test_endpoint_arc_fractions_are_refused(source):
    with pytest.raises(SchemaError, match="interior positions"):
        make_protocol(
            source["root"],
            source["files"],
            candidates={"arc_fractions": [0.0, 1 / 6, 2 / 6, 3 / 6, 4 / 6, 5 / 6]},
        )


@pytest.mark.parametrize(
    ("candidates", "message"),
    [
        ({"proposal_ids": []}, "non-empty list"),
        ({"proposal_ids": ["p1", "p1"]}, "listed twice"),
        ({"from_batch": "explore"}, "name the candidates"),
    ],
)
def test_a_candidate_list_must_name_each_candidate_once(source, candidates, message):
    with pytest.raises(SchemaError, match=message):
        make_protocol(source["root"], source["files"], candidates=candidates)


def test_a_candidate_list_round_trips_and_its_absence_keeps_the_hash(source, protocol):
    listed = make_protocol(
        source["root"],
        source["files"],
        candidates={"from_batch": "explore", "proposal_ids": ["p1", "p2"]},
    )
    revived = MaturationProtocol.from_dict(listed.to_dict())
    assert revived.candidates.proposal_ids == ("p1", "p2")
    assert revived.sha256 == listed.sha256
    assert "proposal_ids" not in protocol.to_dict()["candidates"]


# ---- proposal provenance ----------------------------------------------------


def test_observed_frame_without_frame_index_is_refused(tmp_path):
    with pytest.raises(SchemaError, match="frame it came from"):
        make_proposal_kind(tmp_path, proposal_kind="observed_frame", frame_index=None)


def test_observed_frame_cannot_claim_interpolation(tmp_path):
    with pytest.raises(SchemaError, match="interpolation"):
        make_proposal_kind(
            tmp_path,
            proposal_kind="observed_frame",
            interpolation={"method": "linear", "weights": [1.0]},
        )


def test_interpolated_candidate_is_a_legal_proposal(tmp_path):
    proposal = make_proposal_kind(
        tmp_path,
        proposal_kind="path_interpolation",
        interpolation={"method": "linear", "weights": [0.5, 0.5]},
    )
    assert proposal.to_dict()["interpolation"]["method"] == "linear"


def test_endpoint_idpp_needs_its_method_named(tmp_path):
    with pytest.raises(SchemaError, match="idpp"):
        make_proposal_kind(
            tmp_path, proposal_kind="endpoint_idpp", interpolation={"method": "linear"}
        )


def test_proposal_without_region_is_refused(tmp_path):
    with pytest.raises(SchemaError, match="region"):
        make_proposal_kind(tmp_path, region_graph_sha256s={})


def make_proposal_kind(tmp_path, **overrides):

    kwargs = {
        "proposal_id": "p-test",
        "batch_role": "explore",
        "channel_id": "rc0000",
        "geometry_file": "geometries/p-test.extxyz",
        "geometry_sha256": fake_sha("geometry"),
        "n_atoms": 4,
        "proposal_kind": "observed_frame",
        "frame_index": 5,
        "region_graph_sha256s": {"A": fake_sha("ga"), "B": fake_sha("gb")},
    }
    kwargs.update(overrides)
    kwargs.pop("proposal_id", None)
    kwargs["proposal_id"] = "p-test"
    return PathProposal(**kwargs)


# ---- statistics -------------------------------------------------------------


def test_wilson_matches_the_published_confirm_batch_values():
    # runs/p1_first_entry_confirm/verdict.json, scan 0: B = 38 of A+B = 50;
    # tests never read runs/, so the published numbers are quoted here
    assert wilson_interval(38, 50) == pytest.approx([0.6258731623866494, 0.8570260860478336])


def test_wilson_at_8_of_16_is_the_feasibility_notes_interval():
    # GEN_COMPAS_FEASIBILITY_2026-10-03.md, corrected judgement 4:
    # "当 B 恰占一半时,N=16 的 Wilson 95% 区间约 [0.28,0.72]"
    low, high = wilson_interval(8, 16)
    assert round(low, 2) == 0.28
    assert round(high, 2) == 0.72


def test_wilson_with_no_denominator_is_empty():
    assert wilson_interval(0, 0) == []


def _record(outcome, proposal_id="p1", **overrides):
    kwargs = dict(
        proposal_id=proposal_id,
        shot_id="s0000",
        batch_role="explore",
        seed=1,
        protocol_sha256=fake_sha("protocol"),
        implementation_sha256=fake_sha("impl"),
        outcome=outcome,
        trial_dir="shots/p1/s0000",
    )
    kwargs.update(overrides)
    return ShotRecord(**kwargs)


def test_ab_outcome_without_evidence_is_refused():
    with pytest.raises(SchemaError, match="evidence"):
        _record("A")


def test_excluded_outcomes_need_a_reason():
    with pytest.raises(SchemaError, match="outcome_reason"):
        _record("unabsorbed")
    _record("unabsorbed", outcome_reason="cap reached without k consecutive steps")


def test_p_b_is_null_while_no_shot_reaches_a_side():
    records = [
        _record("unabsorbed", shot_id="s0000", outcome_reason="cap"),
        _record("other", shot_id="s0001", outcome_reason="third basin"),
    ]
    estimate = ShootingEstimate.from_records(records)
    assert estimate.conditional_denominator == 0
    assert estimate.p_B_given_AB is None
    assert estimate.wilson_95 is None
    assert estimate.band_class == "unresolved"
    assert estimate.to_dict()["counts"]["unabsorbed"] == 1


def test_band_classification_uses_the_interval_not_the_point():
    # 16/0: point estimate 0, but the claim is the interval -- [0, 0.194]
    # is entirely below the [0.4, 0.6] band, so it is A-side
    all_a = [
        _record("A", shot_id=f"s{i:04d}", entry_side="A", entry_step=200) for i in range(16)
    ]
    estimate = ShootingEstimate.from_records(all_a)
    assert estimate.p_B_given_AB == 0.0
    low, high = estimate.wilson_95
    assert high < 0.4
    assert estimate.band_class == "A_side"

    all_b = [
        _record("B", shot_id=f"s{i:04d}", entry_side="B", entry_step=900) for i in range(16)
    ]
    mirrored = ShootingEstimate.from_records(all_b)
    assert mirrored.wilson_95[0] > 0.6
    assert mirrored.band_class == "B_side"

    half = [
        _record("A", shot_id=f"s{i:04d}", entry_side="A", entry_step=200) for i in range(8)
    ] + [
        _record("B", shot_id=f"s{i:04d}", entry_side="B", entry_step=900) for i in range(8, 16)
    ]
    compatible = ShootingEstimate.from_records(half)
    assert compatible.p_B_given_AB == pytest.approx(0.5)
    assert compatible.band_class == "transition_region_compatible"

    # 12/4 straddles the band's lower edge: the interval [0.25, 0.60]-ish
    # overlaps [0.4, 0.6] without containing it -- still only "compatible"
    mixed = all_a[:12] + all_b[:4]
    assert ShootingEstimate.from_records(mixed).band_class == "transition_region_compatible"

    # 4/12 straddles the UPPER edge from inside: [0.46, 0.91]-ish. Until
    # 2026-10-11 this overlap fell through to "unresolved".
    upper = all_a[:4] + all_b[:12]
    estimate = ShootingEstimate.from_records(upper)
    assert 0.4 < estimate.wilson_95[0] < 0.6 < estimate.wilson_95[1]
    assert estimate.band_class == "transition_region_compatible"


def test_only_an_explore_batch_may_stop_a_candidate_early():
    from mature.schema import ShootingRule

    rule = ShootingRule(
        batch_role="explore",
        seed_base=0,
        shots_per_candidate=16,
        interim_look=8,
        drift_tolerance_eV_atom=0.01,
    )
    assert ShootingRule.from_dict(rule.to_dict()) == rule
    with pytest.raises(SchemaError, match="only an explore batch"):
        ShootingRule(
            batch_role="confirm",
            seed_base=0,
            confirm_shots=100,
            selection_rule="frozen",
            interim_look=8,
            drift_tolerance_eV_atom=0.01,
        )
    with pytest.raises(SchemaError, match="smaller than shots_per_candidate"):
        replace(rule, interim_look=16)


def test_estimates_spanning_two_proposals_are_refused():
    with pytest.raises(SchemaError, match="one proposal"):
        ShootingEstimate.from_records(
            [
                _record("A", entry_side="A", entry_step=200),
                _record("A", proposal_id="p2", entry_side="A", entry_step=200),
            ]
        )


def test_entry_quench_mismatch_is_counted_not_absorbed():
    records = [
        _record(
            "A",
            shot_id="s0000",
            entry_side="A",
            entry_step=200,
            landed_microstate="c0000/m0000",
            entry_quench_consistent=True,
        ),
        _record(
            "B",
            shot_id="s0001",
            entry_side="B",
            entry_step=300,
            landed_microstate="c0000/m0000",  # quench says A-side: mismatch
            entry_quench_consistent=False,
        ),
    ]
    estimate = ShootingEstimate.from_records(records)
    assert estimate.entry_quench_mismatches == 1
    # the raw entry-side counts survive; the mismatch invalidates interpretation,
    # it does not re-judge the shot (plan M3)
    assert estimate.counts["A"] == 1 and estimate.counts["B"] == 1
    assert math.isfinite(estimate.p_B_given_AB)
