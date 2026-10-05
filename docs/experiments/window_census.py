"""Branch C2 of PLAN_2026-09-25, offline part: which windows does the current rule miss?

search.py builds a bracket only where a return sits directly below a departure. PRRS-2
(REVIEW_2026-09-24 section 4) keeps every adjacent pair of different outcomes. Before
writing that, count on the runs already on disk what the difference would be. Pure
JSON, no dynamics: every scan's samples are sorted by amplitude and each adjacent pair
of different outcomes is classified.

Outcome of a sample:
  returned    completed, landed on the source microstate
  departed:X  completed, landed on microstate X (X != source)
  ridge       ts_candidate (the quench stopped on a saddle)
  transient   unconfirmed with an observed crossing
  unresolved  anything else (quench failed, gate rejection, unverified minimum, ...)

Window kinds:
  current     returned -> departed / ridge / transient   (what search.py brackets)
  reverse     departed / ridge / transient -> returned    (a return above a departure)
  switch      departed:X -> departed:Y, X != Y            (two different products)
  unresolved  either side unresolved

    python docs/experiments/window_census.py
"""

import collections
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "runs" / "window_census"


def outcome(sample, record, source):
    status = sample.get("status")
    if status == "completed":
        target = sample.get("target_microstate")
        if target == source:
            return "returned"
        return f"departed:{target}" if target else "unresolved"
    if status == "ts_candidate":
        return "ridge"
    if status == "unconfirmed" and (record.get("crossing") or {}).get("observed"):
        return "transient"
    return "unresolved"


def kind(low, high):
    left = low.split(":")[0]
    right = high.split(":")[0]
    if "unresolved" in (left, right):
        return "unresolved"
    if left == "returned":
        return "current"
    if right == "returned":
        return "reverse"
    if left == right == "departed":
        return "switch"
    return "other"  # e.g. ridge -> transient: both left the basin, differently


def main():
    totals = collections.Counter()
    runs = []
    for path in sorted((ROOT / "runs").glob("**/network.json")):
        network = json.loads(path.read_text(encoding="utf-8-sig"))
        scans = network.get("amplitude_scans") or []
        if not scans:
            continue
        records = {}
        attempts = path.parent / "attempts.jsonl"
        if attempts.exists():
            for line in attempts.open():
                if line.strip():
                    record = json.loads(line)
                    records[record["attempt_id"]] = record
        per_run = collections.Counter()
        examples = []
        for index, scan in enumerate(scans):
            source = scan.get("source_microstate")
            samples = sorted(scan.get("samples") or [], key=lambda s: s["amplitude"])
            pairs = [
                (s["amplitude"], outcome(s, records.get(s.get("attempt"), {}), source))
                for s in samples
            ]
            for (a, la), (b, lb) in zip(pairs, pairs[1:]):
                if la == lb:
                    continue
                k = kind(la, lb)
                per_run[k] += 1
                if k != "current" and len(examples) < 5:
                    examples.append({"scan": index, "window": [a, b], "outcome": [la, lb]})
            per_run["scans"] += 1
            per_run["interrupted_brackets"] += sum(
                1 for b in scan.get("local_brackets") or [] if b.get("interrupted_by")
            )
        totals.update(per_run)
        runs.append(
            {"run": str(path.parent.relative_to(ROOT)), **per_run, "examples": examples}
        )
        print(
            f"{str(path.parent.relative_to(ROOT)):55s} scans {per_run['scans']:3d}  "
            + "  ".join(
                f"{k} {per_run[k]}"
                for k in ("current", "reverse", "switch", "other", "unresolved")
            )
            + f"  interrupted {per_run['interrupted_brackets']}",
            flush=True,
        )
    print(
        "TOTAL "
        + "  ".join(
            f"{k} {totals[k]}"
            for k in (
                "scans",
                "current",
                "reverse",
                "switch",
                "other",
                "unresolved",
                "interrupted_brackets",
            )
        ),
        flush=True,
    )
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "census.json").write_text(
        json.dumps({"totals": totals, "runs": runs}, indent=1) + "\n"
    )
    print(f"wrote {OUT / 'census.json'}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
