"""Saddle-search rates with explicit, non-interchangeable denominators.

The audit that prompted this rewrite: the previous three rates read well and could not be
wrong, because each one silently excluded the cases that would have made it look bad. In
particular `requested_channel_coverage` put only channels that had already been observed in
its denominator, so "the target coordinate was never sampled at all" reported n/a rather
than 0 -- the single most important failure P2 exhibited was invisible in the metric meant
to catch it.

So every stage is reported with the population it is conditional on, and a pre-registered
truth set is scored separately against the whole run:

  considered            directions whose scan attempted a saddle search
  frame_materialised    a seed frame was actually read from disk
  climb_called          follow_min_mode ran
  ts_found              it returned an index-one saddle
  attribution_complete  both sides of the descent produced an admitted endpoint
  target_supported      the channel the scan requested is among the ones it supports

Each rate names its numerator and denominator. None of them is the run's reaction discovery
rate: that is a truth-set question and is reported as recall over a set fixed in advance.

Usage: saddle_metrics.py <run directory> [truth.json]

  truth.json: {"channels": [{"broken": [[4,5]], "formed": [[0,5]]}], ...}
              scored by whether any reaction_channel matches, so a coordinate that was
              never probed scores 0 rather than dropping out of the denominator.
"""
import json
import sys
from collections import defaultdict
from pathlib import Path


def load(run):
    network = json.loads((Path(run) / "network.json").read_text())
    attempts = []
    path = Path(run) / "attempts.jsonl"
    if path.exists():
        attempts = [json.loads(line) for line in path.read_text().splitlines()]
    return network, attempts


def _attribution(candidate, requested, network):
    ends = candidate.get("connects") or []
    complete = len(ends) == 2 and all(end.get("microstate") for end in ends)
    supports = [c["id"] for c in network["reaction_channels"]
                if candidate["id"] in c.get("ts_support", ())]
    if not complete or requested is None:
        return complete, None, supports
    return True, requested in supports, supports


def stages(network):
    by_id = {c["id"]: c for c in network["ts_candidates"]}
    rows = []
    for scan in network.get("amplitude_scans", []):
        report = scan.get("saddle_search") or {}
        if not report.get("attempted"):
            continue
        requested = report.get("channel")
        seeds = report.get("seeds", [])
        if not seeds:
            rows.append({"scan": f"{scan['family']}{scan['indices']}", "requested": requested,
                         "materialised": False, "climbed": False, "found": False,
                         "complete": False, "hit": None, "reason": report.get("reason")})
            continue
        for seed in seeds:
            candidate = by_id.get(seed.get("ts_candidate"))
            found = candidate is not None
            complete, hit, supports = _attribution(candidate, requested, network) \
                if found else (False, None, [])
            rows.append({"scan": f"{scan['family']}{scan['indices']}",
                         "seed_attempt": seed.get("seed_attempt"),
                         "transient": seed.get("transient"),
                         "amplitude": seed.get("amplitude"),
                         "frame": seed.get("seed_frame"),
                         "requested": requested,
                         "materialised": seed.get("seed_frame") is not None,
                         "climbed": seed.get("steps") is not None
                                    or seed.get("reason") is not None,
                         "found": found, "complete": complete, "hit": hit,
                         "supports": supports, "reason": seed.get("reason"),
                         "initial_rise_eV": seed.get("initial_rise_eV"),
                         "ts": seed.get("ts_candidate")})
    return rows


def truth_recall(network, truth):
    """Recall over a set of channels fixed before the run.

    Scored against the whole run, so a target that was never probed scores zero instead of
    leaving the denominator -- which is the failure the old coverage rate could not see.
    """
    observed = [(sorted(map(sorted, c["broken"])), sorted(map(sorted, c["formed"])),
                 c["id"]) for c in network["reaction_channels"]]
    found = []
    for wanted in truth.get("channels", []):
        broken = sorted(map(sorted, wanted["broken"]))
        formed = sorted(map(sorted, wanted["formed"]))
        match = next((cid for b, f, cid in observed
                      if (b, f) in ((broken, formed), (formed, broken))), None)
        found.append((wanted, match))
    return found


def main(run, truth_path=None):
    network, _ = load(run)
    rows = stages(network)

    def rate(name, numerator, denominator, conditional):
        value = "n/a" if not denominator else f"{numerator}/{denominator} = {numerator/denominator:.2f}"
        print(f"  {name:<26} {value:<18} | given {conditional}")

    considered = len(rows)
    materialised = [r for r in rows if r["materialised"]]
    climbed = [r for r in materialised if r["climbed"]]
    found = [r for r in climbed if r["found"]]
    complete = [r for r in found if r["complete"]]
    scored = [r for r in complete if r["hit"] is not None]
    hits = [r for r in scored if r["hit"]]

    print("stage rates (each conditional on the previous stage)")
    rate("frame_materialised", len(materialised), considered, "a saddle search was attempted")
    rate("climb_called", len(climbed), len(materialised), "a seed frame was read")
    rate("ts_found", len(found), len(climbed), "the climb ran")
    rate("attribution_complete", len(complete), len(found), "an index-one saddle was found")
    rate("target_scored", len(scored), len(complete),
         "attribution completed (null where the scan requested nothing)")
    rate("target_supported", len(hits), len(scored), "a target existed to score against")

    print("\nby seed kind (found / considered)")
    kinds = defaultdict(lambda: [0, 0])
    for row in rows:
        label = {True: "transient", False: "confirmed_or_boundary"}.get(row.get("transient"),
                                                                       "unknown")
        kinds[label][0] += int(row["found"])
        kinds[label][1] += 1
    for label, (hit, total) in sorted(kinds.items()):
        print(f"  {label:<24} {hit}/{total}")

    print("\nby amplitude (found / considered)")
    amplitudes = defaultdict(lambda: [0, 0])
    for row in rows:
        if row.get("amplitude") is None:
            continue
        amplitudes[row["amplitude"]][0] += int(row["found"])
        amplitudes[row["amplitude"]][1] += 1
    for amplitude, (hit, total) in sorted(amplitudes.items()):
        print(f"  {amplitude:<24} {hit}/{total}")

    if truth_path:
        truth = json.loads(Path(truth_path).read_text())
        found_truth = truth_recall(network, truth)
        matched = sum(1 for _, cid in found_truth if cid)
        print(f"\npre-registered truth set (scored against the whole run)")
        print(f"  recall                     {matched}/{len(found_truth)}"
              f" = {matched/len(found_truth):.2f}" if found_truth else "  recall n/a")
        for wanted, cid in found_truth:
            print(f"    broken={wanted['broken']} formed={wanted['formed']} -> "
                  f"{cid or 'NOT FOUND'}")

    print("\nper seed")
    for row in rows:
        verdict = ("no_seed" if not row["materialised"] else
                   "not_found" if not row["found"] else
                   "unattributed" if not row["complete"] else
                   "no_target" if row["hit"] is None else
                   "hit" if row["hit"] else "other_channel")
        print(f"  {row.get('seed_attempt')} ampl={row.get('amplitude')} "
              f"frame={row.get('frame')} -> {verdict}  ts={row.get('ts')} "
              f"requested={row['requested']} supports={row.get('supports')}"
              + (f"  reason={row['reason']}" if row.get("reason") else ""))


if __name__ == "__main__":
    main(*sys.argv[1:3])
