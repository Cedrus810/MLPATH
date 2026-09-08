"""Read the descent-recheck shards and say which of the three answers the data supports.

Prints one line per saddle and then the only thing step 3 has to deliver: whether the six
+1 descents agree, and if they do, which answer that is. No new computation, no new
tolerance, nothing read from anywhere but the shard json.
"""
import json, sys
from collections import Counter
from pathlib import Path

OUT = Path(sys.argv[1] if len(sys.argv) > 1 else "runs/p3_descent_recheck")
shards = sorted(OUT.glob("descent_recheck_shard*.json"))
if not shards:
    raise SystemExit(f"no shard files in {OUT}")

records, meta = [], None
for path in shards:
    payload = json.load(open(path))
    meta = meta or payload
    records.extend(payload["records"])
records.sort(key=lambda r: r["index"])

print(f"{len(records)} saddles from {len(shards)} shards   "
      f"quench_fmax={meta['quench_fmax_eV_A']}  steps={meta['quench_steps']}  sign=+1")
print(f"A {meta['key_A']}   B {meta['key_B']}\n")
head = f"{'saddle':22s} {'reason':34s} {'key':17s} {'min':5s} {'methyl':7s} {'motion':26s} verdict"
print(head); print("-" * len(head))
for r in records:
    if r.get("status") != "descended":
        print(f"{r['saddle_file']:22s} {r['status']}")
        continue
    key = r["final_key"] + (" (A)" if r["final_is_A"] else " (B)" if r["final_is_B"] else " (?)")
    frac = r.get("final_negative_mode_methyl_fraction")
    print(f"{r['saddle_file']:22s} {str(r['reason'])[:34]:34s} {key:17s} "
          f"{str(r['final_confirmed_minimum']):5s} "
          f"{('-' if frac is None else format(frac, '.3f')):7s} "
          f"{r['motion']['verdict']:26s} {r['verdict']}")

done = [r for r in records if r.get("status") == "descended"]
# Whether "order 0" is a statement about a stationary point depends on whether the
# geometry it was read on satisfies the force tolerance. Before the fix in runner.quench
# it did not; the reading below has to say which case this run is.
force_ok = bool(done) and all(r["final_fmax_eV_A"] <= meta["quench_fmax_eV_A"] for r in done)
print(f"\nendpoints within quench_fmax ({meta['quench_fmax_eV_A']}): {force_ok}  "
      f"fmax {min(r['final_fmax_eV_A'] for r in done):.2e} .. "
      f"{max(r['final_fmax_eV_A'] for r in done):.2e}" if done else "")
production = [r.get("production_path") for r in done if r.get("production_path")]
if production:
    reached = sum(1 for x in production if x.get("reached"))
    print(f"production path (descend_saddle): {reached}/{len(production)} returned an "
          f"endpoint; A: {sum(1 for x in production if x.get('is_A'))}; "
          f"admitted with unsettled soft modes: "
          f"{sum(1 for x in production if x.get('admitted_with_unsettled_soft_modes'))}")

verdicts = Counter(r.get("verdict") for r in done)
print("\nverdicts:", dict(verdicts))
if len(verdicts) == 1:
    verdict = next(iter(verdicts))
    reading = {
        "plus_side_reaches_B_too":
            "answer 3 -- both descents of these saddles land in B, so on the OFF24 surface "
            "they do not join A and B. A model boundary, not a search failure.",
        "plus_side_reaches_a_third_basin":
            "answer 3 -- the +1 side leaves for a basin that is neither A nor B.",
        "reached_A_and_confirms":
            "answer 2 -- the descent enters A's basin and the endpoint carries no negative "
            "curvature, so the archived refusal was the quench protocol, not the surface."
            + ("" if force_ok else
               " Note what this does NOT say: every endpoint here sits after an "
               "un-minimized polish displacement (fmax above quench_fmax), so it is "
               "evidence of the basin, not of a stationary point."),
        "reached_A_bookkeeping_refusal":
            "answer 2 -- the endpoint is in A's basin with no negative mode, refused only "
            "because the last polish round moved it.",
        "reached_A_rotor_holds_the_mode":
            "answer 1 -- the endpoint is in A's basin and the residual negative mode is the "
            "methyl rotor, which is the OFF24 softness §3.9 measured (0.00262 vs POLAR-1's "
            "0.00470). Model boundary.",
        "reached_A_but_not_a_minimum":
            "unresolved -- in A's basin, not a minimum, and the negative mode is not the "
            "rotor. Do not force this into one of the three answers.",
    }.get(verdict, "unrecognised verdict")
    print(f"\nall six agree: {verdict}\n  {reading}")
else:
    print("\nthe six do NOT agree; the saddles are not one phenomenon and step 3's "
          "three-way question has no single answer here.")
