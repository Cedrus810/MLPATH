"""The direction contract: what the proposal and scan orchestration owes the record.

Written before the scheduling slice is extracted from ReactionSearch.run, which is the
order the 2026-09-09 ruling sets: contract tests first, then move only that slice, never
a rewrite of the 1200-line function (docs/DECISIONS_PENDING.md, first section; the same
reason as REVIEW_2026-09-05.md 5.2 -- run() captures several mutable registries and
there is no orchestration-level net under it).

Two kinds of test live here and they are labelled:

    characterisation   pins behaviour that is already correct, so the extraction cannot
                       change it silently
    specification      xfail(strict), pins behaviour the ruling requires and the code
                       does not have yet. strict means it fails if it starts passing
                       while still marked, so it cannot rot into a lie

The distinction matters because the negative results this project reports depend on it:
a null is only a result if the ledger can say the direction was searched rather than
never selected, never executed, or refused at a gate.
"""

from dataclasses import replace

import pytest
from ase.build import molecule

from prrs.config import SearchConfig
from prrs.perturbations import propose


def _proposal(atoms, config, seed=42):
    report = {}
    batches = list(propose(atoms, config, seed, report=report))
    directions = [(b[0].family, tuple(b[0].indices), b[0].sign) for b in batches]
    return directions, report


@pytest.fixture
def ethanol():
    return molecule("CH3CH2OH")


# ---------------------------------------------------------------- characterisation


def test_tracked_coordinates_do_not_steer_the_search(ethanol):
    """Ruling point 1: tracked_coordinates audit coverage, they never select directions.

    Naming a coordinate is not expecting it (internal.named_values says so), and letting
    a named coordinate influence which directions are tried would feed the answer back
    into the search -- benchmarks/README.md rule 5. The test is byte-equality of the
    selection, not a weaker "still finds it": anything short of that is a channel through
    which the expectation could leak.
    """
    plain = SearchConfig()
    named = replace(
        plain,
        tracked_coordinates=(
            {"name": "r_O_H", "kind": "bond", "indices": (2, 3)},
            {"name": "q", "kind": "bond_difference", "indices": (0, 1, 2, 3)},
        ),
    )
    assert named.tracked_coordinates != plain.tracked_coordinates
    assert _proposal(ethanol, plain) == _proposal(ethanol, named)


def test_a_proposal_is_reproducible_and_without_replacement(ethanol):
    """The same (structure, config, seed) proposes the same directions, each once.

    Reproducibility is what makes a run's ledger auditable after the fact; no-replacement
    is what makes "candidates minus selected" meaningful. `propose` already removes each
    entry from its family pool, so this is characterisation -- it is here because the
    extraction must not lose it.
    """
    config = SearchConfig()
    first, _ = _proposal(ethanol, config)
    second, _ = _proposal(ethanol, config)
    assert first == second
    assert len(set(first)) == len(first), "a direction was proposed twice"


def test_the_report_reconciles_with_what_was_yielded(ethanol):
    """The ledger has to add up: per family, and against the batches actually produced."""
    config = SearchConfig()
    directions, report = _proposal(ethanol, config)
    families = report["families"]
    assert sum(f["selected"] for f in families.values()) == report["selected_total"]
    assert report["selected_total"] == len(directions)
    for family, entry in families.items():
        assert 0 <= entry["selected"] <= entry["quota"] <= entry["candidates"]
        assert sum(1 for f, _, _ in directions if f == family) == entry["selected"]
    assert report["cap_bound"] is (report["selected_total"] >= config.max_directions)


def test_coverage_is_truncated_even_when_the_budget_is_not(ethanol):
    """The fact the ruling rests on: the quota fraction censors coverage on its own.

    b04's probe3 is the measured case -- 28 candidates, 14 selected, `cap_bound false` on
    a six-atom system, so max_directions never bound and `direction_quota_fraction = 0.5`
    did the truncating. This is therefore NOT a large-molecule problem, and a null result
    cannot be read as "does not exist" at any system size. Asserted as the mechanism
    rather than as b04's numbers so it holds for whatever structure is passed in.
    """
    generous = replace(SearchConfig(), max_directions=10_000)
    _, report = _proposal(ethanol, generous)
    candidates = sum(f["candidates"] for f in report["families"].values())
    assert report["cap_bound"] is False, "this test is about the quota, not the cap"
    assert report["selected_total"] < candidates, (
        "with an unbounded cap the quota fraction still leaves candidates unsampled: "
        f"{report['selected_total']} of {candidates}"
    )
    for entry in report["families"].values():
        assert entry["quota"] <= entry["candidates"]


# ---------------------------------------------------------------- contract, implemented
# These two were specifications (xfail) for exactly as long as it took to write the
# ledger fields in perturbations.propose; the marks came off in the same change, because
# a strict xfail that passes is itself a failure and that is the point of strict.


def test_a_truncated_proposal_declares_its_coverage_censored(ethanol):
    """Ruling point 3: when the budget truncates coverage, the record must say so.

    Not a warning and not a number to be inferred by the reader: a status, because it
    decides whether a negative conclusion from this run is admissible at all.
    """
    _, report = _proposal(ethanol, SearchConfig())
    assert report["coverage_censored"] is True
    assert report["censored_by"] in ("direction_quota_fraction", "max_directions")


def test_every_candidate_has_a_stable_id(ethanol):
    """Ruling point 2: ids that survive a change in what gets selected.

    Without them the candidate -> selected -> executed -> refused chain cannot be
    written down, and every subdivision of `unanswered` is unattributable.
    """
    config = SearchConfig()
    _, wide = _proposal(ethanol, replace(config, max_directions=10_000))
    _, narrow = _proposal(ethanol, replace(config, max_directions=4))
    assert wide["candidate_ids"], "the report carries no candidate ids"
    assert set(narrow["candidate_ids"]) <= set(wide["candidate_ids"])
    assert set(narrow["selected_ids"]) <= set(narrow["candidate_ids"])


# ---------------------------------------------------------------- specification


def test_a_direction_selected_but_not_executed_is_not_counted_as_searched(tmp_path):
    """Ruling point 5: max_trials belongs to the same ledger as the direction budget.

    A run cut short by the trial budget leaves selected directions unexecuted. Nothing
    used to write them down, so a coordinate that was never tried was indistinguishable,
    in the record, from one that was tried and missed -- `proposed_not_executed` versus
    `searched_no_hit`, which the new criteria semantics treat as different answers.

    The ledger has TWO units and they must not be added together. An earlier draft of
    this test did exactly that -- `selected_total` counts DIRECTIONS while a scan's
    `samples` are the AMPLITUDES of one direction -- so it asserted an identity false by
    construction. Both are asserted below, separately, and the system is chosen to
    exercise both at once: three atoms give more than one direction, three amplitudes
    give a batch that a two-trial budget cuts in half.
    """
    from test_boundary_continuation import SharedProton, shared_proton_atoms

    from prrs.search import ReactionSearch

    config = SearchConfig(
        families=("stretch", "compress"),
        geometry_amplitudes_A=(0.1, 0.3, 0.6),
        timestep_fs=0.1,
        response_steps=10,
        closed_shell_only=False,
        max_depth=1,
        max_trials=2,
        refinement_steps=0,
        saddle_search_enabled=False,
    )
    network = ReactionSearch(SharedProton, config).run(
        shared_proton_atoms(-SharedProton.width), tmp_path / "cut"
    )
    assert network["termination"] == "trial_budget"

    # direction units: was this coordinate tried at all
    ledger = network["direction_ledger"]
    assert ledger["directions_selected"] == (
        ledger["directions_executed"] + ledger["directions_not_executed"]
    )
    assert ledger["directions_not_executed"] > 0, ledger
    assert ledger["not_executed_reasons"]["trial_budget"] > 0
    assert set(ledger["not_executed_ids"]).isdisjoint(ledger["executed_ids"])
    assert ledger["coverage_censored"] is True

    # amplitude units: how far along the one direction that ran did we get
    scan = network["amplitude_scans"][0]
    schedule = scan["proposal"]["schedule"]
    assert schedule["amplitudes_offered"] == (
        schedule["amplitudes_executed"] + schedule["amplitudes_not_executed"]
    )
    assert schedule["amplitudes_executed"] == len(scan["samples"])
    assert schedule["amplitudes_not_executed"] > 0
    assert schedule["halted_by"] == "trial_budget"


def test_a_run_that_dropped_nothing_says_so(tmp_path):
    """ "Nothing was dropped" is the only state in which a null is admissible.

    So it has to be a positive statement in the record rather than the absence of a
    complaint. The generous budget below exhausts the proposals instead of the trials,
    and the ledger must then show every selected direction executed -- while still
    reporting `coverage_censored` if the quota fraction left candidates unsampled, which
    is a different question and is not silently merged into this one.
    """
    network = _tiny_run(tmp_path, max_trials=500)
    ledger = network["direction_ledger"]
    assert network["termination"] == "configured_scope_exhausted"
    assert ledger["directions_not_executed"] == 0
    assert ledger["not_executed_ids"] == []
    assert ledger["directions_executed"] == ledger["directions_selected"]


# ------------------------------------------------- characterisation, orchestration
# The net that has to exist before the scheduling slice moves out of ReactionSearch.run.


def _tiny_run(tmp_path, **overrides):
    from prrs.calculators import DoubleWell, demo_atoms
    from prrs.search import ReactionSearch

    config = SearchConfig(
        **{
            "families": ("stretch", "compress"),
            "geometry_amplitudes_A": (0.1, 0.3),
            "timestep_fs": 0.1,
            "response_steps": 10,
            "closed_shell_only": False,
            "max_depth": 1,
            "max_trials": 12,
            "refinement_steps": 0,
            "saddle_search_enabled": False,
            **overrides,
        }
    )
    return ReactionSearch(DoubleWell, config).run(demo_atoms(), tmp_path / "run")


def test_every_executed_trial_belongs_to_exactly_one_scan(tmp_path):
    """A trial that ran is accounted for once, in the scan whose direction produced it.

    This is the property the ledger will be built on: if an attempt can appear in two
    scans, or in none, then "how many amplitudes of this direction were tried" has no
    answer and neither does `proposed_not_executed`.
    """
    network = _tiny_run(tmp_path)
    from collections import Counter

    seen = Counter(
        sample["attempt"] for scan in network["amplitude_scans"] for sample in scan["samples"]
    )
    assert seen, "the run produced no scans"
    assert max(seen.values()) == 1, f"an attempt appears in two scans: {seen.most_common(2)}"
    # Cross-check against what actually landed on disk: network["attempts"] is a COUNT,
    # the records are in attempts.jsonl, and every coarse one must appear in a scan.
    import json

    coarse = {
        json.loads(line)["attempt_id"]
        for line in (tmp_path / "run" / "attempts.jsonl").read_text().splitlines()
        if json.loads(line).get("stage") == "coarse"
    }
    assert coarse == set(seen), "a coarse attempt is missing from the scans"


def test_a_scan_reports_the_proposal_that_produced_it(tmp_path):
    """Each scan carries its own proposal budget, not the run's last one."""
    network = _tiny_run(tmp_path)
    for scan in network["amplitude_scans"]:
        budget = scan["proposal"]["budget"]
        assert budget["families"], scan["family"]
        assert scan["family"] in budget["families"]
        assert budget["selected_ids"], "the proposal ledger reached the scan record"


def test_the_termination_reason_distinguishes_budget_from_scope(tmp_path):
    """ "Ran out of trials" and "ran out of proposals" are different endings.

    Both are recorded today; the test pins them because the ledger's `not_executed`
    reasons have to agree with this field, and a run that ended on scope must not claim
    directions were dropped for budget.
    """
    generous = _tiny_run(tmp_path / "a", max_trials=500)
    assert generous["termination"] == "configured_scope_exhausted"
    tight = _tiny_run(tmp_path / "b", max_trials=2)
    assert tight["termination"] == "trial_budget"


# ------------------------------------------------------------------- the two modes


def test_exhaustive_mode_leaves_no_candidate_unsampled(ethanol):
    """Ruling point 4: `exhaustive` traverses every symmetry-reduced candidate.

    Symmetry-reduced is the whole enumeration here -- `propose` already collapses
    equivalent directions by canonical colour, so "every candidate" does not mean every
    atom pair. What it does mean is that `coverage_censored` comes back false, which is
    the only state in which the new criteria semantics admit a negative result.
    """
    sampled = SearchConfig()
    exhaustive = replace(sampled, direction_mode="exhaustive")

    _, sampled_report = _proposal(ethanol, sampled)
    directions, report = _proposal(ethanol, exhaustive)

    assert report["mode"] == "exhaustive"
    assert report["coverage_censored"] is False
    assert report["uncovered_ids"] == []
    assert report["cap_bound"] is False, "the cap may not decide coverage in this mode"
    assert len(directions) == len(report["candidate_ids"])
    assert report["selected_total"] > sampled_report["selected_total"]
    for family, entry in report["families"].items():
        assert entry["selected"] == entry["candidates"], family


def test_exhaustive_mode_refuses_to_start_when_the_budget_cannot_pay(tmp_path):
    """Ruling point 4: fail before starting, not halfway through.

    The message has to name both numbers, because the two ways out are different
    decisions: raise the budget, or accept sampled mode and the limit that comes with
    every null it produces.
    """
    from prrs.calculators import DoubleWell, demo_atoms
    from prrs.search import ReactionSearch

    config = SearchConfig(
        families=("stretch", "compress"),
        geometry_amplitudes_A=(0.1, 0.3, 0.6),
        direction_mode="exhaustive",
        max_trials=2,
        timestep_fs=0.1,
        response_steps=10,
        closed_shell_only=False,
        max_depth=1,
        saddle_search_enabled=False,
    )
    output = tmp_path / "refused"
    with pytest.raises(ValueError, match="exhaustive"):
        ReactionSearch(DoubleWell, config).run(demo_atoms(), output)

    # Refused before spending: no trial ran. The failure record is allowed to exist --
    # that is the 2026-09-05 fix -- but no attempt may.
    trials = output / "trials"
    assert not trials.exists() or not list(trials.iterdir())


def test_the_mode_is_recorded_so_a_null_can_be_read_correctly(ethanol):
    """A run's own record must say which claim it is entitled to make."""
    for mode in ("sampled", "exhaustive"):
        _, report = _proposal(ethanol, replace(SearchConfig(), direction_mode=mode))
        assert report["mode"] == mode
        assert "coverage_censored" in report
        assert (report["coverage_censored"] is False) == (mode == "exhaustive")


def test_the_run_level_ledger_closes_across_several_microstates(tmp_path):
    """The roll-up sums (microstate, direction) pairs; it does not unique direction ids.

    b04's probe4 is why this test exists. Two microstates were searched, each offered the
    same 14 direction ids, and the first roll-up took set unions: 28 selected against 17
    executed and 7 not executed -- an identity that does not close, in the one record
    whose job is to say whether a null is admissible. A direction searched on m0000 was
    not searched on m0001, so uniquing across microstates discards exactly what the
    ledger carries. The single-microstate test above could not catch it.
    """
    from test_boundary_continuation import SharedProton, shared_proton_atoms

    from prrs.search import ReactionSearch

    # The reactive settings are test_reaction_events.py's, which are known to transfer
    # the proton and so to admit a second microstate; collateral_gate is off there for a
    # stated reason -- this potential holds r01 + r12, so driving either bond necessarily
    # moves the other and the coupling IS the coordinate.
    config = SearchConfig(
        seed=11,
        families=("stretch", "compress"),
        geometry_amplitudes_A=(0.3, 0.6, 0.9),
        quench_fmax_eV_A=0.01,
        closed_shell_only=False,
        minimum_check="probe",
        temperature_K=0.0,
        collateral_gate=False,
        max_depth=1,
        max_trials=60,
        conformer_trials_per_microstate=4,  # forces directions to be dropped per state
        refinement_steps=0,
        saddle_search_enabled=False,
    )
    network = ReactionSearch(SharedProton, config).run(
        shared_proton_atoms(-0.4), tmp_path / "several"
    )
    ledger = network["direction_ledger"]
    per_state = network["direction_ledger_per_microstate"]
    assert ledger["microstates_searched"] == len(per_state) >= 2, per_state

    assert ledger["directions_selected"] == (
        ledger["directions_executed"] + ledger["directions_not_executed"]
    ), ledger
    for field in ("directions_selected", "directions_executed", "directions_not_executed"):
        assert ledger[field] == sum(entry[field] for entry in per_state.values()), field
    for entry in per_state.values():
        assert entry["directions_selected"] == (
            entry["directions_executed"] + entry["directions_not_executed"]
        ), entry

    # `termination` names the trial budget only, so whatever else dropped directions has
    # to be named separately -- otherwise "configured_scope_exhausted" reads as
    # "everything was tried", which is the claim this whole contract exists to refuse.
    detail = network["termination_detail"]
    assert detail["trials_spent"] <= detail["max_trials"]
    assert detail["directions_dropped_by"] == ledger["not_executed_reasons"]
    if ledger["directions_not_executed"]:
        assert detail["directions_dropped_by"], (
            "directions were dropped and the record does not say by what"
        )


def test_a_batch_cut_short_mid_way_names_the_budget_that_cut_it(tmp_path):
    """The reason has to be THIS refusal's, not the one that stopped the last direction.

    The line under test read `allowed, _ = registry.may_search(...)` and then reported
    the outer loop's `blocked`, so a direction whose amplitudes were cut off part-way by
    `conformer_trials_per_microstate` was filed under a generic `node_budget` -- or under
    whatever had been refused before it. b04's probe4 never reached that branch, so its
    ledger looked right; a ledger's job is to be right in the cases nobody has run yet.

    The setup makes the cut land INSIDE a batch: three amplitudes per direction against a
    per-microstate budget that is not a multiple of three.
    """
    from test_boundary_continuation import SharedProton, shared_proton_atoms

    from prrs.search import ReactionSearch

    config = SearchConfig(
        seed=11,
        families=("stretch", "compress"),
        geometry_amplitudes_A=(0.3, 0.6, 0.9),
        quench_fmax_eV_A=0.01,
        closed_shell_only=False,
        minimum_check="probe",
        temperature_K=0.0,
        collateral_gate=False,
        max_depth=1,
        max_trials=60,
        conformer_trials_per_microstate=2,  # 2 of a 3-amplitude batch, so the cut is mid-batch
        refinement_steps=0,
        saddle_search_enabled=False,
    )
    network = ReactionSearch(SharedProton, config).run(
        shared_proton_atoms(-0.4), tmp_path / "midbatch"
    )

    schedules = [
        scan["proposal"]["schedule"]
        for scan in network["amplitude_scans"]
        if scan["proposal"]["schedule"]["amplitudes_not_executed"] > 0
    ]
    assert schedules, "no batch was cut short; the budget above is not tight enough"
    for schedule in schedules:
        assert schedule["halted_by"] in (
            "conformer_trials_per_microstate",
            "chemical_trials_per_node",
            "trial_budget",
        ), schedule
    # The specific budget, not the generic word: `node_budget` is the fallback for a
    # refusal that gave no reason, and may_search always gives one here.
    assert any(
        schedule["halted_by"] == "conformer_trials_per_microstate" for schedule in schedules
    ), [s["halted_by"] for s in schedules]
    assert all(schedule["halted_by"] != "node_budget" for schedule in schedules)


def test_exhaustive_mode_checks_every_budget_that_can_truncate_it(tmp_path):
    """Not only max_trials: the per-microstate and per-node budgets cut the same run.

    b05's probe2 found this the way these things get found -- 24 directions selected, 9
    executed, 15 dropped by `conformer_trials_per_microstate`, all while max_trials still
    had 64 trials left. An exhaustive run wants 504 trials on that structure, so a raised
    max_trials alone would have been waved through the pre-check and then truncated at 32
    by a budget nobody looked at. "Refused before starting" has to mean all of them.
    """
    from prrs.calculators import DoubleWell, demo_atoms
    from prrs.search import ReactionSearch

    common = dict(
        families=("stretch", "compress"),
        geometry_amplitudes_A=(0.1, 0.3, 0.6),
        direction_mode="exhaustive",
        timestep_fs=0.1,
        response_steps=10,
        closed_shell_only=False,
        max_depth=1,
        saddle_search_enabled=False,
    )
    # Generous trial budget, tight per-microstate budget: the old check passed this.
    config = SearchConfig(
        max_trials=10_000,
        chemical_trials_per_node=10_000,
        conformer_trials_per_microstate=1,
        **common,
    )
    with pytest.raises(ValueError, match="conformer_trials_per_microstate"):
        ReactionSearch(DoubleWell, config).run(demo_atoms(), tmp_path / "per_state")

    # And the per-node budget, on its own.
    config = SearchConfig(
        max_trials=10_000,
        chemical_trials_per_node=1,
        conformer_trials_per_microstate=10_000,
        **common,
    )
    with pytest.raises(ValueError, match="chemical_trials_per_node"):
        ReactionSearch(DoubleWell, config).run(demo_atoms(), tmp_path / "per_node")

    # All three generous: it starts.
    config = SearchConfig(
        max_trials=10_000,
        chemical_trials_per_node=10_000,
        conformer_trials_per_microstate=10_000,
        **common,
    )
    network = ReactionSearch(DoubleWell, config).run(demo_atoms(), tmp_path / "ok")
    assert network["status"] == "completed"
    assert network["direction_ledger"]["directions_not_executed"] == 0
