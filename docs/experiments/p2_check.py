"""Evaluate the frozen P2 acceptance criteria against a finished run.

No criterion is computed here that was not written down in
docs/P2_OXOBUTANAL_ACCEPTANCE.md first. The wording of the criteria is quoted next to each
check so a drift between the doc and the code is visible rather than silent.

Usage: python p2_check.py runs/p2_round5/network.json
"""
import json
import re
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from prrs.chemistry import stereo_only_difference   # noqa: E402  the single implementation

net = json.load(open(sys.argv[1]))
results = []


def check(tag, ok, detail=""):
    results.append((tag, bool(ok), detail))


nodes = net["chemical_nodes"]
channels = net.get("reaction_channels", [])
edges = net["reactions"]
ts = net["ts_candidates"]

# --- A. network expansion ---
check("A1  chemical_nodes >= 2", len(nodes) >= 2, f"got {len(nodes)}")

# A2 is an EXISTENCE criterion: "two nodes exist whose keys differ and whose fragments are
# both ['C4H6O2']". With two nodes there is one such pair and A4's "those two nodes" is
# unambiguous. With three -- which happens as soon as the E/Z isomer is found -- there are
# three qualifying pairs and the frozen text does not say which. This loop used to keep
# whichever pair came last, i.e. (c0001, c0002), while the transfer channel joins
# (c0000, c0001): six criteria then failed (A4 A5 C4 D1 D2 and, separately, B1b) on four
# seeds whose runs had nothing wrong with them.
#
# So all qualifying pairs are kept, and A4 asks whether a reaction channel joins ANY of
# them. That is not a loosening: the arbitrary single choice was never the intent, and a
# criterion that depends on iteration order was not measuring anything. The pair that
# actually carries the transfer is then the one A5/C4/D1/D2 are evaluated against.
qualifying = []
for first in range(len(nodes)):
    for second in range(first + 1, len(nodes)):
        a, b = nodes[first], nodes[second]
        if (a["chemical_key"] != b["chemical_key"]
                and a["fragments"] == b["fragments"] == ["C4H6O2"]):
            qualifying.append((a, b))
pair = qualifying[0] if qualifying else None
check("A2  two nodes, keys differ and fragments identical (['C4H6O2'])", pair is not None,
      "" if pair is None else f"{pair[0]['id']} {pair[0]['chemical_key'][:8]} / "
                              f"{pair[1]['id']} {pair[1]['chemical_key'][:8]}  "
                              f"graphs {pair[0]['key_components']['graph_hash']} / "
                              f"{pair[1]['key_components']['graph_hash']}")

directed = [e for e in edges if e["source"] != e["target"]]
check("A3  len(reactions) >= 1 and the edge has source != target", bool(directed),
      f"reactions {len(edges)}, of which directed {len(directed)}"
      + ("" if not directed else
         "  mechanisms " + str(sorted({m for e in directed
                                       for m in e.get('mechanisms', [])}))))

target, ends = None, set()
for candidate_a, candidate_b in qualifying:
    candidate_ends = {candidate_a["id"], candidate_b["id"]}
    for channel in channels:
        if (channel["class"] == "reaction" and channel["self_loop"] is False
                and set(channel["ends"]) == candidate_ends):
            target, ends = channel, candidate_ends
check("A4  a reaction channel, not a self-loop, joining two such nodes",
      target is not None, "" if target is None else
      f"{target['id']} ends={target['ends']}  "
      f"({len(qualifying)} qualifying pair(s) considered)")


def is_oh(bond, structure_numbers):
    return sorted(structure_numbers[i] for i in bond) == [1, 8]


numbers = None
if pair is not None:
    numbers = [8, 6, 6, 6, 8, 1, 1, 1, 6, 1, 1, 1]      # from the P2 input, O C C C O H...
ok_a5, detail_a5 = False, "no target channel"
if target is not None and numbers is not None:
    broken, formed = target["broken"], target["formed"]
    if len(broken) == len(formed) == 1:
        shared = set(broken[0]) & set(formed[0])
        ok_a5 = (is_oh(broken[0], numbers) and is_oh(formed[0], numbers)
                 and len(shared) == 1 and numbers[shared.pop()] == 1
                 and sorted(broken[0]) != sorted(formed[0]))
    detail_a5 = f"broken={broken} formed={formed}"
check("A5  exactly one O-H broken and one O-H formed, same H", ok_a5, detail_a5)

# --- B. must not happen ---
by_id = {node["id"]: node for node in nodes}
empty_delta = [c for c in channels if not c["broken"] and not c["formed"]]
stereo_verdicts = {}
for channel in empty_delta:
    ends = [by_id.get(e) for e in channel["ends"]]
    if len(ends) == 2 and all(ends):
        # key_components, not the node dicts: the shared implementation in chemistry.py
        # takes the components themselves. Passing nodes made every field it looks for
        # absent, so it reported "no parity differs" about a genuine parity flip.
        stereo_verdicts[channel["id"]] = stereo_only_difference(
            ends[0]["key_components"], ends[1]["key_components"])
    else:
        stereo_verdicts[channel["id"]] = (False, "endpoints not both in the node list")

# REVISED 2026-09-02, see docs/P2_OXOBUTANAL_ACCEPTANCE.md "B1/B3 修订". The frozen text
# forbade every key_change_without_bond_change channel on the assumption that such a
# channel is always a defect. It is not: E/Z isomerisation changes configurational
# stereochemistry with no bond change, and a search told to explore freely is entitled to
# find it -- it is simply a higher-energy state, not an error. What the criterion should
# forbid is an empty-delta channel that is NOT explained by a parity flip.
bad_empty = [f"{cid}: {why}" for cid, (ok, why) in stereo_verdicts.items() if not ok]
# The class name the model now assigns to a real parity flip. The old name stays a defect.
LEGITIMATE_EMPTY_DELTA = "configurational_isomerisation"
# A run made before the fourth cell was split carries the old label on a channel that is
# genuinely a parity flip. That is a data-vintage fact, not a wrong classification, and the
# two must not read the same: one says "re-run this", the other says "fix the model".
misnamed, stale = [], []
for channel in empty_delta:
    if not stereo_verdicts[channel["id"]][0]:
        continue
    if channel["class"] == LEGITIMATE_EMPTY_DELTA:
        continue
    if channel["class"] == "key_change_without_bond_change":
        stale.append(channel["id"])
    else:
        misnamed.append(f"{channel['id']} class={channel['class']}")
check("B1a every empty-bond-delta channel is a parity flip, not an unexplained key change",
      not bad_empty, "; ".join(bad_empty) or
      ("; ".join(f"{cid}: {why}" for cid, (_, why) in stereo_verdicts.items())
       or "no empty-delta channels"))
# Same criterion, second half: a parity flip must also be NAMED as one by the model, not
# left under the defect name. Splitting this out keeps "is it real" and "does the model
# know it is real" from being one indivisible pass/fail.
flips = sum(1 for ok, _ in stereo_verdicts.values() if ok)
check("B1b a parity flip is classified configurational_isomerisation, not the defect name",
      not misnamed and not stale,
      "; ".join(misnamed) if misnamed else
      (f"{stale} still carry the pre-split label 'key_change_without_bond_change'; "
       "this run predates the classifier fix and must be re-run to be scored on B1b"
       if stale else f"{flips} parity flips, all named"))


def closed_shell(formula):
    counts = Counter()
    for element, number in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        if element:
            counts[element] += int(number or 1)
    electrons = {"H": 1, "C": 4, "N": 5, "O": 6, "F": 7, "Cl": 17, "Br": 35, "I": 53}
    return sum(electrons.get(e, 0) * n for e, n in counts.items()) % 2 == 0


odd = [(node["id"], f) for node in nodes for f in node["fragments"]
       if not closed_shell(f)]
check("B2  every microstate's fragments closed-shell", not odd, str(odd))

# REVISED with B1 and for the same reason: B3 as frozen was the same test in different
# words. Rotors cannot reach reaction_channels at all -- an unchanged key makes them
# conformational_transition, which record_channel drops -- so the only thing B3 could ever
# fire on was the E/Z case B1 now handles. What is still checked is that an empty-delta
# channel never joins two nodes whose keys are the SAME, which would be a genuine defect.
same_key_empty = [c["id"] for c in empty_delta
                  if len({by_id[e]["chemical_key"] for e in c["ends"] if e in by_id}) == 1
                  and not c["self_loop"]]
check("B3  no empty-delta channel between two nodes of the same key",
      not same_key_empty, str(same_key_empty) or f"{len(empty_delta)} empty-delta channels")

# --- C. curvature ---
check("C1  every TS saddle_order == 1",
      all(t["saddle_order"] == 1 for t in ts), f"{len(ts)} TS")
sources = Counter(t.get("curvature_source") for t in ts)
check("C2  every curvature_source analytic",
      set(sources) == {"analytic"} if ts else False, str(dict(sources)))
worst = max((t.get("trivial_mode_floor_worst_residual") or 0.0) for t in ts) if ts else None
check("C3  every trivial floor residual < 5e-4", worst is not None and worst < 5e-4,
      f"max {worst}")
supporting = [] if target is None else [t for t in ts if t["id"] in target["ts_support"]]
check("C4  supporting TS have response.sign_stable",
      bool(supporting) and all((t.get("response") or {}).get("sign_stable")
                               for t in supporting),
      f"{len(supporting)} supporting TS")

# --- D. TS and attribution ---
two_sided = []
for candidate in supporting:
    connects = candidate.get("connects") or []
    named = [end.get("chemical_node") for end in connects if end.get("microstate")]
    if len(named) == 2:
        two_sided.append((candidate["id"], named))
check("D1  at least one TS supports that channel by two-sided descent", bool(two_sided),
      str(two_sided))
check("D2  that TS's two descents land in two different chemical nodes",
      any(len(set(named)) == 2 for _, named in two_sided), str(two_sided))

tristate = []
for scan in net.get("amplitude_scans", []):
    for seed in (scan.get("saddle_search") or {}).get("seeds", []):
        if "supports_requested" in seed:
            tristate.append((seed.get("ts_candidate"), seed.get("attribution_complete"),
                             seed.get("supports_requested")))
violations = [row for row in tristate if row[1] is False and row[2] is not None]
check("D3  incomplete attribution records supports_requested = null, never false",
      not violations, f"{len(tristate)} scored seeds, {len(violations)} violations")

# --- integrity checks: NOT part of the frozen fifteen ---
# These do not score the run against the acceptance criteria; they check that the record is
# internally consistent. Reported separately and never folded into the count, because adding
# criteria after seeing data is how a frozen set stops meaning anything.
#
# The first one exists because P2 round five scored 15/15 while two of the three
# continuations on its directed edge came from rungs whose own verdict was
# "stopped_on_saddle, not evidence of arrival" -- a verdict-callback side effect leaked the
# endpoint out. The criteria could not see it; it was caught by reading the rungs by hand.
integrity = []


def integrity_check(tag, ok, detail=""):
    integrity.append((tag, bool(ok), detail))


walk_by_seed = {}
for scan in net.get("amplitude_scans", []):
    for entry in (scan.get("boundary_continuation") or {}).get("seeds", []):
        if entry.get("continuation") is not None:
            walk_by_seed.setdefault(entry["seed_attempt"], []).append(entry)

unsupported = []
for edge in edges:
    for continuation in edge.get("continuations") or []:
        attempt = continuation["seed_attempt"]
        entries = walk_by_seed.get(attempt, [])
        arrived = [e for e in entries
                   if (e["continuation"].get("outcome") or "").startswith("reached")]
        if not arrived:
            outcomes = [e["continuation"].get("outcome") for e in entries] or ["no walk"]
            unsupported.append(f"{edge['id']}<-{attempt} (walk outcome {outcomes})")
integrity_check("I1  every continuation on an edge comes from a walk that reached somewhere",
                not unsupported, "; ".join(unsupported) or
                f"{sum(len(e.get('continuations') or []) for e in edges)} continuations")

mismatched = []
for edge in edges:
    for attempt in edge["attempts"]:
        record = None
        path = Path(sys.argv[1]).parent / "trials" / attempt / "result.json"
        if path.exists():
            record = json.loads(path.read_text())
        if record and (record.get("source") != edge["source"]
                       or record.get("target") != edge["target"]):
            mismatched.append(f"{edge['id']}<-{attempt}")
integrity_check("I2  every attempt on an edge has a trial record naming the same target",
                not mismatched, "; ".join(mismatched) or
                f"{sum(len(e['attempts']) for e in edges)} attempts")

orphan_support = [t["id"] for t in ts
                  if t.get("reaction_channel")
                  and t["id"] not in
                  next((c["ts_support"] for c in channels
                        if c["id"] == t.get("reaction_channel")), [])]
integrity_check("I3  a TS naming a channel appears in that channel's ts_support",
                not orphan_support, str(orphan_support) or f"{len(ts)} TS")

width = max(len(tag) for tag, _, _ in results)
print("=" * 78)
for tag, ok, detail in results:
    print(f"{'PASS' if ok else 'FAIL'}  {tag:<{width}}  {detail}")
print("=" * 78)
print(f"{sum(1 for _, ok, _ in results if ok)}/{len(results)} criteria passed")
if integrity:
    iwidth = max(len(tag) for tag, _, _ in integrity)
    print("\n--- record integrity (NOT part of the frozen criteria) ---")
    for tag, ok, detail in integrity:
        print(f"{'ok  ' if ok else 'BAD '}  {tag:<{iwidth}}  {detail}")
    broken = [tag for tag, ok, _ in integrity if not ok]
    if broken:
        print(f"  {len(broken)} integrity problem(s): the criteria score above cannot be "
              "trusted until they are explained.")

print("\n--- observations ---")
for channel in channels:
    print(f"  channel {channel['id']} class={channel['class']} "
          f"self_loop={channel['self_loop']} ends={channel['ends']} "
          f"broken={channel['broken']} formed={channel['formed']} "
          f"observed_by={len(channel['observed_by'])} ts_support={channel['ts_support']}")
for edge in edges:
    print(f"  edge {edge['id']} {edge['source']} -> {edge['target']} "
          f"mechanisms={edge.get('mechanisms')} attempts={edge['attempts']} "
          f"continuations={edge.get('continuations')}")
for candidate in ts:
    response = candidate.get("response") or {}
    print(f"  {candidate['id']} {candidate['found_by']} E={candidate['energy_eV']:.5f} "
          f"icm={[round(v, 1) for v in (candidate['imaginary_wavenumbers_icm'] or [])]} "
          f"kappa={response.get('kappa')} sign_stable={response.get('sign_stable')} "
          f"channel={candidate.get('reaction_channel')}")
for node in nodes:
    energies = [round(m["energy_eV"], 6) for m in node["microstates"]]
    print(f"  node {node['id']} key={node['chemical_key'][:12]} "
          f"graph={node['key_components']['graph_hash']} "
          f"frags={node['fragments']} microstates={energies}")
if len(nodes) == 2:
    low = nodes[0]["microstates"][0]["energy_eV"]
    high = nodes[1]["microstates"][0]["energy_eV"]
    print(f"  dE(second - first) = {(high - low) * 1e3:+.3f} meV")
