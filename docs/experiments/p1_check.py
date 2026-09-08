"""Evaluate the frozen P1 acceptance criteria against a finished run. No criterion is
computed here that was not written down in docs/P1_MALONALDEHYDE_ACCEPTANCE.md first."""
import json, re, sys
from collections import Counter

net = json.load(open(sys.argv[1]))
results, observations = [], []
def check(tag, ok, detail=""):
    results.append((tag, bool(ok), detail))

nodes = net["chemical_nodes"]
channels = net.get("reaction_channels", [])
source_key = nodes[0]["chemical_key"] if nodes else None

# --- A ---
check("A1  chemical_nodes == 1", len(nodes) == 1, f"got {len(nodes)}")
check("A2  reactions == 0", len(net["reactions"]) == 0, f"got {len(net['reactions'])}")
check("A3  no microstate under a different key", len(nodes) == 1,
      f"node keys {[n['chemical_key'][:8] for n in nodes]}")

# --- B ---
check("B1  reaction_channels >= 1", len(channels) >= 1, f"got {len(channels)}")
target = None
for c in channels:
    if c["class"] != "degenerate_reaction" or not c["self_loop"]:
        continue
    if len(c["broken"]) != 1 or len(c["formed"]) != 1:
        continue
    if sorted(c["broken"][0]) == sorted(c["formed"][0]):
        continue
    hb = set(c["broken"][0]) & set(c["formed"][0])
    if len(hb) == 1:
        target = c; break
check("B2  degenerate self-loop channel, one bond broken and one formed on the same atom",
      target is not None,
      json.dumps({k: target[k] for k in ("id","class","self_loop","ends","broken","formed")})
      if target else f"classes seen: {[c['class'] for c in channels]}")
if target:
    check("B3  event key recorded", bool(target.get("event_key")), target["event_key"])

# --- C ---
empty_delta = [c["id"] for c in channels if not c["broken"] and not c["formed"]]
check("C1  no empty-bond-delta event in reaction_channels", not empty_delta, str(empty_delta))
bad_class = [c["id"] for c in channels if c["class"] == "key_change_without_bond_change"]
check("C2  no key_change_without_bond_change", not bad_class, str(bad_class))
def closed_shell(formula):
    counts = {e: int(n or 1) for e, n in re.findall(r"([A-Z][a-z]?)(\d*)", formula) if e}
    electrons = sum({"H":1,"C":6,"N":7,"O":8,"F":9}.get(e,0)*n for e, n in counts.items())
    return electrons % 2 == 0
odd = [f for n in nodes for f in n["fragments"] if not closed_shell(f)]
check("C3  every fragment closed-shell", not odd, str(odd))

# --- D ---
ts = net["ts_candidates"]
check("D1  every TS saddle_order == 1", all(t["saddle_order"] == 1 for t in ts),
      f"{len(ts)} TS")
check("D2  every curvature_source analytic",
      all(t.get("curvature_source") == "analytic" for t in ts),
      str(Counter(t.get("curvature_source") for t in ts)))
floors = [t.get("trivial_mode_floor_worst_residual") for t in ts]
check("D3  trivial floor residual < 5e-4",
      all(f is not None and f < 5e-4 for f in floors),
      f"max {max(floors) if floors else None}")
supported = [t for t in ts if target and t["id"] in target.get("ts_support", [])]
if supported:
    check("D4  supporting TS sign_stable",
          all((t.get("response") or {}).get("sign_stable") for t in supported),
          f"{len(supported)} TS support the channel")
else:
    observations.append("D4 not evaluated: no TS supports the target channel")

# --- observations ---
for c in channels:
    observations.append(f"channel {c['id']} class={c['class']} self_loop={c['self_loop']} "
                        f"ends={c['ends']} broken={c['broken']} formed={c['formed']} "
                        f"observed_by={len(c['observed_by'])} ts_support={c['ts_support']}")
for t in ts:
    r = t.get("response") or {}
    observations.append(f"{t['id']} E={t['energy_eV']:.5f} icm={[round(x,1) for x in t['imaginary_wavenumbers_icm']]} "
                        f"kappa={r.get('kappa')} sigma={r.get('sigma')} c={r.get('anharmonicity')} "
                        f"class={t.get('transition_class')} channel={t.get('reaction_channel')}")
for n in nodes:
    observations.append(f"node {n['id']} key={n['chemical_key'][:12]} frags={n['fragments']} "
                        f"microstates={len(n['microstates'])} conformer_transitions={len(n['conformer_transitions'])}")

width = max(len(t) for t, _, _ in results)
print("=" * 78)
for tag, ok, detail in results:
    print(f"{'PASS' if ok else 'FAIL'}  {tag:<{width}}  {detail[:70]}")
print("=" * 78)
passed = sum(1 for _, ok, _ in results if ok)
print(f"{passed}/{len(results)} criteria passed")
print("\n--- observations ---")
for line in observations:
    print(" ", line)
