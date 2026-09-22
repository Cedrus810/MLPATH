"""The explicit schema for ``network.json`` -- schema 2, described after the fact.

``search.py`` stamped ``schema_version: 2`` on every network since 2026-09-02 while
nothing anywhere defined what schema 2 IS: the on-disk shape was whatever the 16
``network[...]`` assignment sites happened to write, and it drifted with every writer
change. Two of the thirteen consumers are criteria (``b04/check.py``, ``b05/check.py``),
each hand-parsing the shape it needs. Measured over ``runs/`` at the time this module
was written: of 70 ``network.json`` files, 2 are schema 1 (``nodes``/``edges`` era),
3 schema-2 files predate ``reaction_channels``, and only 5 carry the direction ledger
-- so "schema 2" was a version stamp on an undefined object.

This module makes the schema a DESCRIPTION, not a wish: the tables below state what
the current writer puts on disk, entry field by entry, and ``tests/test_schema.py``
runs ``validate`` over every existing ``network.json`` to hold the description against
reality. Files written by older writers fail ``validate`` in recorded ways; their
deviations are pinned in that test as an explicit drift ledger, so a new run whose
network does not match the schema is a red test rather than a silent change.

Deliberately not jsonschema and not pydantic. The schema is one literal dict and the
validator is a handful of loops: anything heavier would be a dependency bought to
express less than the table already says, and the anchor test -- not the validator's
type algebra -- is what keeps the schema true.
"""

import json
from pathlib import Path

SCHEMA_VERSION = 2

# Keys the writer adds only in the final update of a run that reached its natural end.
# A ``partial``/``failed``/``running`` snapshot legitimately lacks all of them.
COMPLETED_DICTS = (
    "termination_detail",
    "outcome_counts",
    "classification_counts",
    "admission_counts",
    "budgets",
    "direction_ledger",
    "direction_ledger_per_microstate",
)

# Top-level keys that must be lists. Scalars are checked individually below --
# schema_version, attempts as ints; status, interpretation, termination as enums
# or presence -- because a single expected-type loop over mixed types is how the
# first draft of this validator called every network broken.
TOP_LEVEL_LISTS = (
    "chemical_nodes",
    "reactions",
    "reaction_channels",
    "ts_candidates",
    "ts_connected_channels",
    "amplitude_scans",
    "withheld_candidates",
)

STATUS_VALUES = frozenset({"running", "completed", "partial", "failed"})
TERMINATION_VALUES = frozenset({"trial_budget", "configured_scope_exhausted"})

# Per-collection required fields, exactly as the current writer builds each entry.
# Optional-by-code fields (ts_support, connects, transition_class, ...) stay optional
# here too: the schema describes what the writer always writes, not what it sometimes
# writes, so optionality below must trace to a branch in search.py, not to optimism.
SCHEMA_2 = {
    "chemical_nodes": {
        "required": (
            "id",
            "chemical_key",
            "depth",
            "discovered_by",
            "fragments",
            "key_components",
            "trials_spent",
            "reservoir",
            "microstates",
            "conformer_transitions",
        ),
        "meaning": (
            "one node per chemical_key; microstates are geometries within it, "
            "reservoir.evicted records what conformer budget cost"
        ),
    },
    "microstates (inside chemical_nodes)": {
        "required": (
            "id",
            "structure",
            "energy_eV",
            "discovered_by",
            "trials_spent",
            "torsions_rad",
            "free_bonds",
            "outcomes",
            "reservoir_score",
        ),
        "meaning": "Microstate.summary() verbatim; the reservoir's scoring unit",
    },
    "reactions": {
        "required": (
            "id",
            "source",
            "target",
            "observed_only",
            "attempts",
            "classifications",
            "best_by_family",
            "source_microstates",
            "mechanisms",
            "continuations",
        ),
        "meaning": (
            "directed observations only; ts_support is derived and therefore optional. "
            "continuations are not attempts (see record_reaction note [1])"
        ),
    },
    "reaction_channels": {
        "required": (
            "id",
            "event_key",
            "class",
            "ends",
            "self_loop",
            "broken",
            "formed",
            "symmetry_enumerated",
            "undirected",
            "observed_by",
            "ts_support",
            "meaning",
        ),
        "meaning": (
            "undirected reaction identity by event_key; self_loop = degenerate "
            "reaction, not a conformational change. event_key is identity WITHIN "
            "one atom numbering (see test_quotient_invariance) until structural_key"
        ),
    },
    "ts_candidates": {
        "required": (
            "id",
            "attempt",
            "structure",
            "reached_from",
            "saddle_order",
            "found_by",
            "imaginary_wavenumbers_icm",
            "curvature_source",
            "trivial_mode_floor_worst_residual",
            "response",
            "energy_eV",
            "negatively_curved_torsions",
            "status",
        ),
        "meaning": (
            "order-1 saddles; connects/connects_meaning only when the two-sided "
            "descent ran (irc_enabled), transition_class/reaction_channel only then"
        ),
    },
    "ts_connected_channels": {
        "required": (
            "ts",
            "undirected",
            "ends",
            "end_microstates",
            "transition_class",
            "meaning",
        ),
        "meaning": "undirected structural evidence from descending both saddle sides",
    },
    "amplitude_scans": {
        "required": (
            "source",
            "source_microstate",
            "family",
            "indices",
            "sign",
            "seed",
            "unit",
            "local_brackets",
            "samples",
            "nonmonotonic_return_observed",
            "proposal",
            "saddle_search",
            "boundary_continuation",
        ),
        "meaning": (
            "one scan per (microstate, direction); proposal holds the direction "
            "ledger for this scan's amplitudes, local_brackets the observed windows"
        ),
    },
    "withheld_candidates": {
        "required": (
            "source",
            "attempt",
            "reason",
            "detail",
            "confirmed_minimum",
            "structure",
            "energy_eV",
            "meaning",
        ),
        "meaning": (
            "minima that passed confirm_minimum and were kept out by a label; "
            "the retrieval path for a QM verification pass"
        ),
    },
    "direction_ledger": {
        "required": (
            "directions_selected",
            "directions_executed",
            "directions_not_executed",
            "microstates_searched",
            "not_executed_reasons",
            "not_executed_ids",
            "executed_ids",
            "uncovered_ids",
            "coverage_censored",
            "censored_by",
            "meaning",
        ),
        "meaning": (
            "(microstate, direction) pairs summed over microstates; ids lists are "
            "unions and NOT in that unit. direction_ledger_per_microstate holds the "
            "same fields per microstate"
        ),
    },
}

# Collections checked entry-by-entry; microstates are nested inside chemical_nodes
# entries and validated through the same table.
_COLLECTION_KEYS = tuple(
    key for key in SCHEMA_2 if " " not in key and key != "direction_ledger"
)


def validate(network):
    """Return the list of ways ``network`` deviates from schema 2, newest writer.

    Problems are strings, never exceptions: the validator's job is to describe a
    gap, and the caller decides what to do -- ``publish()`` records the list, tests
    and ``--strict-schema`` treat non-empty as failure.
    """
    problems = []
    if network.get("schema_version") != SCHEMA_VERSION:
        problems.append(
            f"top: schema_version {network.get('schema_version')!r} != {SCHEMA_VERSION}"
        )
        return problems
    for key in TOP_LEVEL_LISTS:
        if not isinstance(network.get(key), list):
            problems.append(f"top: '{key}' is not a list")
    if not isinstance(network.get("attempts"), int):
        problems.append("top: 'attempts' is not an int")
    status = network.get("status")
    if status not in STATUS_VALUES:
        problems.append(f"top: status {status!r} not in {sorted(STATUS_VALUES)}")
    if not network.get("interpretation"):
        problems.append("top: missing 'interpretation'")
    if status == "completed":
        if network.get("termination") not in TERMINATION_VALUES:
            problems.append(
                f"top: termination {network.get('termination')!r} not in "
                f"{sorted(TERMINATION_VALUES)}"
            )
        for key in COMPLETED_DICTS:
            if not isinstance(network.get(key), dict):
                problems.append(f"top: completed run missing '{key}'")

    for key in _COLLECTION_KEYS:
        required = SCHEMA_2[key]["required"]
        entries = network.get(key) or []
        for i, entry in enumerate(entries):
            missing = [f for f in required if f not in entry]
            if missing:
                problems.append(f"{key}[{i}]: missing {', '.join(missing)}")
    microstate_required = SCHEMA_2["microstates (inside chemical_nodes)"]["required"]
    for i, node in enumerate(network.get("chemical_nodes") or []):
        for j, state in enumerate(node.get("microstates") or []):
            missing = [f for f in microstate_required if f not in state]
            if missing:
                problems.append(
                    f"chemical_nodes[{i}].microstates[{j}]: missing {', '.join(missing)}"
                )
    ledger_required = SCHEMA_2["direction_ledger"]["required"]
    # microstates_searched is a run-level total; the per-microstate ledger() has no
    # such field, and demanding it there was the second bug of the first draft.
    per_state_required = tuple(f for f in ledger_required if f != "microstates_searched")
    ledger = network.get("direction_ledger")
    if isinstance(ledger, dict):
        missing = [f for f in ledger_required if f not in ledger]
        if missing:
            problems.append(f"direction_ledger: missing {', '.join(missing)}")
    for state_id, entry in (network.get("direction_ledger_per_microstate") or {}).items():
        if not isinstance(entry, dict):
            problems.append(f"direction_ledger_per_microstate[{state_id}]: not a dict")
            continue
        missing = [f for f in per_state_required if f not in entry]
        if missing:
            problems.append(
                f"direction_ledger_per_microstate[{state_id}]: missing {', '.join(missing)}"
            )
    return problems


def read_network(path):
    """Load a network.json the way every consumer should: schema checked on request.

    The existing consumers are left exactly as they are -- rewriting thirteen readers
    has no benefit. New code calls this so the path and the encoding are not retyped
    per script, and so ``validate`` is one import away.
    """
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))
