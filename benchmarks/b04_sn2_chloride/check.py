"""Mechanical evaluation of b04's criteria, one check per criterion, quoted.

Every criterion below is quoted from CRITERIA.md next to the code that evaluates it, so
drift between the document and the check is visible rather than silent. Nothing is
computed here that was not written down there first -- docs/experiments/p2_check.py
opens with that sentence and it is the whole discipline. The criteria text was written
on 2026-09-09 from the user's ruling; this file came after it, deliberately.

Exit codes say whether the RECORD could be adjudicated, never whether the chemistry
came out a particular way:

    0   adjudicated. The verdict is printed; `unanswered` is a verdict, not an error
    1   the record cannot be adjudicated -- a field the criteria name is missing or two
        fields contradict each other. That is a defect in the run or in this file
    2   called wrong (no run given)

    python benchmarks/b04_sn2_chloride/check.py runs/<the run>/network.json
"""

import json
import sys
from pathlib import Path

ADMITTED = ("existing_microstate", "new_microstate", "new_chemical_node")
CONFIRMED_ONLY_IF = (
    "these three outcomes are reachable only after descend_saddle's confirm_minimum "
    "passed (runner.py); a failure records not_a_minimum or quench_failed instead"
)

findings = []


def criterion(tag, quote, ok, detail):
    findings.append((tag, ok))
    print(f"  [{'PASS' if ok else 'fail'}] {tag}")
    for line in quote.strip().splitlines():
        print(f"         > {line.strip()}")
    for line in str(detail).strip().splitlines():
        print(f"         {line}")
    print()


def main(argv):
    if len(argv) < 2:
        print(__doc__.strip())
        return 2
    path = Path(argv[1])
    network = json.loads(path.read_text())
    # The run's configuration and provenance live in a sibling file, not in the network:
    # publish() rewrites network.json throughout a run, and config.json is written once.
    manifest = path.with_name("config.json")
    config = json.loads(manifest.read_text()) if manifest.exists() else {}
    settings = config.get("config", config)
    print(f"b04 check  ·  {argv[1]}")
    print(
        f"config: direction_mode={settings.get('direction_mode')}  "
        f"seed={settings.get('seed')}  "
        f"max_rise_eV={settings.get('saddle_search_max_rise_eV')}  "
        f"implementation={(config.get('implementation_sha256') or '?')[:16]}"
    )
    if settings.get("direction_mode") is None:
        print("        (this run predates direction_mode; it was sampled -- the field")
        print("         did not exist, so nothing here may be read as exhaustive)")
    print("=" * 78)

    channels = network.get("reaction_channels") or []
    ts_candidates = network.get("ts_candidates") or []

    # ---------------------------------------------------------------- observed_direct
    print("## observed_direct (CRITERIA.md section 2)\n")
    order_one = [t for t in ts_candidates if t.get("saddle_order") == 1]
    criterion(
        "D1  an order-1 saddle candidate exists",
        "存在 `saddle_order == 1` 的鞍点候选",
        bool(order_one),
        f"order-1: {[t['id'] for t in order_one]}   of {len(ts_candidates)} candidates",
    )

    two_sided = []
    for ts in order_one:
        sides = ts.get("connects") or []
        reached = [s.get("microstate") for s in sides if s.get("microstate")]
        if len(set(reached)) >= 2:
            two_sided.append((ts, sides, sorted(set(reached))))
    criterion(
        "D2  the two descents reach two DIFFERENT microstates",
        "两侧下降到达**两个不同的 microstate**。D2 用 microstate 而不是 chemical node，"
        "是这一格的关键：chemical_key 在自环两端相同，节点数扛不动这个证据",
        bool(two_sided),
        "\n".join(
            f"{ts['id']}: "
            + ", ".join(f"sign {s.get('sign'):+d} -> {s.get('microstate')}" for s in sides)
            for ts, sides, _ in two_sided
        )
        or "no order-1 saddle has two distinct descent endpoints",
    )

    confirmed = [
        (ts, sides)
        for ts, sides, _ in two_sided
        if all(s.get("status") in ADMITTED for s in sides)
    ]
    criterion(
        "D3  both endpoints are confirmed minima",
        "两侧端点都是**确认过的极小**：`connects[*].status` 属于 "
        "{existing_microstate, new_microstate, new_chemical_node}",
        bool(confirmed),
        "\n".join(
            f"{ts['id']}: " + ", ".join(str(s.get("status")) for s in sides)
            for ts, sides, _ in two_sided
        )
        + f"\n{CONFIRMED_ONLY_IF}",
    )

    supported = []
    for ts, sides in confirmed:
        for channel in channels:
            if ts["id"] in (channel.get("ts_support") or []):
                if channel.get("broken") and channel.get("formed") and channel.get("event_key"):
                    supported.append((ts, channel))
    criterion(
        "D4  the channel has a named bond-change event supported by that saddle",
        "`reaction_channels[*].broken` / `.formed` 非空，`event_key` 存在，"
        "且 `.ts_support` 含该 TS 的 id",
        bool(supported),
        "\n".join(
            f"{c['id']}: event_key={c['event_key']} broken={c['broken']} formed={c['formed']} "
            f"ts_support={c.get('ts_support')} self_loop={c.get('self_loop')}"
            for _, c in supported
        )
        or "no channel names this saddle in ts_support with a complete bond-change event",
    )

    direct = bool(supported)

    # ------------------------------------------------------------- observed_stepwise
    print("## observed_stepwise (CRITERIA.md section 3)\n")
    endpoints = {
        s.get("microstate")
        for ts in ts_candidates
        for s in (ts.get("connects") or [])
        if s.get("microstate")
    }
    all_microstates = {
        m["id"] for node in network.get("chemical_nodes", []) for m in node["microstates"]
    }
    intermediates = sorted(all_microstates - endpoints)
    criterion(
        "S1  a confirmed minimum on the channel that is not one of its endpoints",
        "通道上存在一个 `confirm_minimum` 通过的极小，且它**不是**该通道的任一端点",
        bool(intermediates),
        f"microstates {sorted(all_microstates)}\nendpoints   {sorted(endpoints)}\n"
        f"candidates  {intermediates or 'none'}",
    )
    stepwise = bool(intermediates)  # S2 is only asked if S1 holds
    if stepwise:
        criterion(
            "S2  evidence connecting it into this channel",
            "从它出发的鞍点或下降记录连到该通道的端点",
            False,
            "not evaluated: S1 produced candidates, so this needs a run that has them",
        )

    # ---------------------------------------------------- admissibility of a negative
    print("## whether a negative record is admissible (CRITERIA.md section 5)\n")
    per_state = network.get("direction_ledger_per_microstate") or {}
    ledger = network.get("direction_ledger") or {}
    if not ledger:
        print("  [fail] the run predates direction_ledger; no negative record from it is")
        print("         admissible, because coverage cannot be established. Re-run on a")
        print("         tree that writes it (implementation fa15a829 or later).\n")
        admissible = {}
        findings.append(("L0 direction_ledger present", False))
    else:
        if per_state and "microstates_searched" not in ledger:
            print("  note: this run's RUN-LEVEL roll-up predates the 2026-09-09 fix and")
            print("        unioned direction ids across microstates, so its three counts")
            print("        do not close. The per-microstate ledgers below are correct and")
            print("        are what section 5 adjudicates on; the id lists are unaffected.")
            print()
        # Section 5 names THREE conditions, not two. The third -- no gate refusal on that
        # microstate for that channel -- is not attributable per microstate in
        # network.json, which carries only the run-level outcome tally. So it is read
        # run-wide and conservatively: any refused trial anywhere censors every
        # microstate's negatives, because the record cannot say which one it belonged to.
        counts = network.get("outcome_counts")
        if counts is None:
            refused = {"outcome_counts": "absent"}
            print("  gate refusals: UNKNOWN -- this run has no outcome_counts, so whether")
            print("     a gate refused anything cannot be established. Absence of a")
            print("     refusal record is not evidence of no refusal, the same reading")
            print("     the no-ledger branch below refuses to make.")
            print()
        else:
            refused = {s: c for s, c in counts.items() if s != "completed"}
            if refused:
                print(f"  gate refusals in this run: {refused}")
                print("     CRITERIA.md 5 third condition: a microstate with a refusal on")
                print("     this channel admits no negative. network.json does not")
                print("     attribute refusals to microstates, so this censors all of them.")
                print()
        admissible = {}
        for state_id, entry in (per_state or {"(run)": ledger}).items():
            clean = (
                not entry.get("uncovered_ids")
                and not entry.get("not_executed_ids")
                and not refused
            )
            admissible[state_id] = clean
            print(
                f"  {state_id}: uncovered={len(entry.get('uncovered_ids') or [])} "
                f"not_executed={len(entry.get('not_executed_ids') or [])} "
                f"censored_by={entry.get('censored_by')} "
                f"refusals={'yes' if refused else 'no'} -> "
                f"{'negatives admissible' if clean else 'negatives NOT admissible'}"
            )
        print()

    # ---------------------------------------------------------------------- verdicts
    print("=" * 78)
    print("VERDICT, per CRITERIA.md section 1\n")
    if direct:
        print("  substitution channel: observed_direct")
        for ts, channel in supported:
            print(
                f"     {ts['id']} order 1, {channel['id']} "
                f"({'self_loop' if channel.get('self_loop') else 'directed'}), "
                f"ends {sorted({s.get('microstate') for s in (ts.get('connects') or [])})}"
            )
    else:
        print("  substitution channel: unanswered   (observed_direct not established)")

    if stepwise:
        print("  stepwise side: candidates exist -- see S2")
    elif not ledger:
        # This branch is here because the first draft of this file got it wrong, and got
        # it wrong in exactly the way the 2026-09-09 ruling exists to prevent: with no
        # ledger it fell through to `searched_no_hit`, the one reason that means "really
        # searched". Absence of a coverage record is not evidence of coverage.
        print("  stepwise side: NOT ADJUDICABLE")
        print("     the run predates direction_ledger, so coverage cannot be established")
        print("     and no reason from CRITERIA.md 4 applies. searched_no_hit in")
        print("     particular is NOT available: it would claim the basis was searched")
        print("     on the strength of the record not saying otherwise.")
    else:
        clean = [state for state, ok in admissible.items() if ok]
        reasons = []
        if ledger.get("uncovered_ids"):
            reasons.append("candidate_not_selected")
        if ledger.get("not_executed_ids"):
            reasons.append("proposed_not_executed")
        if not reasons and clean:
            reasons.append("searched_no_hit")
        print(
            f"  stepwise side: unanswered   reason(s): {', '.join(reasons) or 'undetermined'}"
        )
        if "searched_no_hit" in reasons:
            print(f"     admissible on: {clean}")
            print("     an admissible negative ON THOSE MICROSTATES, INSIDE the enumerated")
            print("     proposal basis. CRITERIA.md 5, verbatim: it is still NOT a")
            print("     statement about the potential energy surface -- the basis is")
            print("     finite and under 'sampled' it is truncated as well, so it supports")
            print("     no claim that anything does not exist.")
        else:
            print("     NOT searched_no_hit: the basis was censored, so nothing here is")
            print("     evidence that an intermediate does not exist (CRITERIA.md 1, 5)")

    print()
    print("  reminder (CRITERIA.md 1.1): even under direction_mode='exhaustive',")
    print("  searched_no_hit proves only 'no hit inside the enumerated basis'. Concerted")
    print("  still requires observed_direct.")

    broken = [tag for tag, ok in findings if not ok and tag.startswith("L")]
    return 1 if broken else 0


def _selftest():
    """Pin the one branch that already went wrong once.

    The first draft, with no ledger in the record, fell through to `searched_no_hit` --
    the single reason that means "really searched". No framework: the whole point is
    that this file can be run anywhere the criteria are, and a defect in the adjudicator
    is a defect in the criteria's enforcement.
    """
    import io
    import contextlib
    import tempfile

    base = {
        "ts_candidates": [
            {
                "id": "ts0000",
                "saddle_order": 1,
                "connects": [
                    {"sign": 1, "status": "existing_microstate", "microstate": "c0000/m0000"},
                    {"sign": -1, "status": "existing_microstate", "microstate": "c0000/m0001"},
                ],
            }
        ],
        "reaction_channels": [
            {
                "id": "rc0000",
                "event_key": "deadbeef",
                "broken": [[0, 1]],
                "formed": [[0, 5]],
                "ts_support": ["ts0000"],
                "self_loop": True,
            }
        ],
        "chemical_nodes": [{"microstates": [{"id": "c0000/m0000"}, {"id": "c0000/m0001"}]}],
    }
    cases = {
        "no ledger": ({}, ("NOT ADJUDICABLE",), ("searched_no_hit",), 1),
        "censored": (
            {"direction_ledger": {"uncovered_ids": ["stretch:0-1:+"], "not_executed_ids": []}},
            ("candidate_not_selected", "NOT searched_no_hit"),
            (),
            0,
        ),
        "clean": (
            {
                "direction_ledger": {"uncovered_ids": [], "not_executed_ids": []},
                "direction_ledger_per_microstate": {
                    "c0000/m0000": {"uncovered_ids": [], "not_executed_ids": []}
                },
                "outcome_counts": {"completed": 12},
            },
            ("searched_no_hit", "observed_direct"),
            ("NOT ADJUDICABLE",),
            0,
        ),
        # Section 5 has THREE conditions. A direction can be executed and still have every
        # one of its amplitudes refused by a gate, and then it was not searched either.
        "gate refused": (
            {
                "direction_ledger": {"uncovered_ids": [], "not_executed_ids": []},
                "direction_ledger_per_microstate": {
                    "c0000/m0000": {"uncovered_ids": [], "not_executed_ids": []}
                },
                "outcome_counts": {"completed": 10, "rejected": 2},
            },
            ("gate refusals in this run", "negatives NOT admissible"),
            ("searched_no_hit",),
            0,
        ),
        # And an absent refusal record is not a record of no refusals.
        "refusals unknown": (
            {
                "direction_ledger": {"uncovered_ids": [], "not_executed_ids": []},
                "direction_ledger_per_microstate": {
                    "c0000/m0000": {"uncovered_ids": [], "not_executed_ids": []}
                },
            },
            ("gate refusals: UNKNOWN",),
            ("searched_no_hit",),
            0,
        ),
    }
    for name, (extra, expected, forbidden, code) in cases.items():
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / "network.json"
            target.write_text(json.dumps({**base, **extra}))
            buffer = io.StringIO()
            with contextlib.redirect_stdout(buffer):
                actual = main(["check.py", str(target)])
            text = buffer.getvalue()
        for phrase in expected:
            assert phrase in text, f"{name}: expected {phrase!r}"
        # Only the VERDICT lines count. The reminder paragraph names searched_no_hit in
        # order to say it is unavailable, and a bare substring search reads that as the
        # verdict -- which is how the first version of this selftest failed a correct run.
        verdict = "\n".join(
            line for line in text.splitlines() if "reason(s)" in line or "channel:" in line
        )
        for phrase in forbidden:
            assert phrase not in verdict, f"{name}: {phrase!r} must not be the verdict"
        assert actual == code, f"{name}: exit {actual}, expected {code}"
        findings.clear()
    print("selftest: 5 cases pass -- no-ledger refuses to adjudicate, a censored basis")
    print("never reports searched_no_hit, a clean ledger does, and a gate refusal")
    print("(or an unknown one) takes searched_no_hit away again.")
    return 0


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        sys.exit(_selftest())
    sys.exit(main(sys.argv))
