"""A near-critical E/Z ladder: the test material a parity-tolerance fix actually needs.

Why the existing corpus cannot serve. Measured over all 300 benzoylacetone conformers plus
P1/P2's persisted structures:

    n_tetrahedral = 0 for every one of the 7 distinct chemical_keys in the golden sample
    -- there is not a single tetrahedral centre anywhere in P0-P3, so the "near-planar
    centre" half of a parity tolerance has nothing to act on here.
    the one locked bond (the enol C=C) sits at |cos| = 0.9993-0.9998 in every frame,
    where cos is the parity discriminant itself: chemistry._locked_bond_parities takes
    sign(dot(u_perp, v_perp)) of the two substituent components perpendicular to the bond
    axis. Flipping needs |cos| -> 0.

So the golden sample can show a parity fix BROKE nothing. It cannot show the fix WORKS,
because it contains no frame anywhere near the discriminant's zero. Those are different
claims and the sample only supports one.

This builds the other one: rotate one substituent about the locked C=C so |cos| sweeps
through zero, and record the discriminant, the parity, the graph hash and the key at each
step. Near |cos| = 0 the sign is decided by a vanishing quantity, which is exactly the
regime a scale-dependent tolerance is supposed to stabilise.

Geometry-only, no potential evaluation, no GPU. Rotation about the bond axis preserves
every bond length, and the run asserts graph_hash is constant across the ladder -- if the
graph moved, the ladder would be measuring connectivity rather than parity and would be
worthless for this.
"""
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
import numpy as np
from ase.io import read, write
from prrs.chemistry import chemical_key, locked_edges, canonical_labels
from prrs.state import encode, resolve_active

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "runs" / "p3_A_source_tight.extxyz"
OUT = ROOT / "runs" / "parity_critical_ladder"
OUT.mkdir(parents=True, exist_ok=True)


def locked_pair_and_substituents(atoms):
    g = encode(atoms)
    idx = resolve_active(atoms, None)
    locked, info = locked_edges(atoms.numbers, g.edges, idx)
    if not locked:
        raise SystemExit(f"no locked bond in {SRC.name}: {info}")
    labels = canonical_labels(atoms.numbers, g.edges, idx)
    adj = {int(i): [] for i in idx}
    for x, y in g.edges:
        adj[x].append(y); adj[y].append(x)
    for (p, q) in sorted(locked):
        picked = []
        for near, far in ((p, q), (q, p)):
            options = [n for n in adj[near] if n != far]
            if not options:
                break
            best = max(labels[n] for n in options)
            if sum(1 for n in options if labels[n] == best) != 1:
                break
            picked.append(next(n for n in options if labels[n] == best))
        if len(picked) == 2:
            return p, q, picked[0], picked[1], adj
    raise SystemExit("locked bond has no distinguishable substituent pair")


def discriminant(atoms, p, q, s0, s1):
    """The parity discriminant, normalised. sign() of the unnormalised form IS the parity."""
    ax = atoms.positions[q] - atoms.positions[p]
    ax = ax / np.linalg.norm(ax)
    u = atoms.positions[s0] - atoms.positions[p]
    v = atoms.positions[s1] - atoms.positions[q]
    u = u - np.dot(u, ax) * ax
    v = v - np.dot(v, ax) * ax
    nu, nv = np.linalg.norm(u), np.linalg.norm(v)
    return float(np.dot(u, v)), float(np.dot(u, v) / (nu * nv))


def side_of(adj, start, blocked):
    """Atoms reachable from `start` without crossing `blocked` -- the group that rotates."""
    seen, stack = {start}, [start]
    while stack:
        cur = stack.pop()
        for n in adj[cur]:
            if n == blocked or n in seen:
                continue
            seen.add(n); stack.append(n)
    return sorted(seen)


A = read(str(SRC))
P, Q, S0, S1, ADJ = locked_pair_and_substituents(A)
raw0, cos0 = discriminant(A, P, Q, S0, S1)
print(f"locked bond {A.symbols[P]}{P}={A.symbols[Q]}{Q}, "
      f"substituents {A.symbols[S0]}{S0} / {A.symbols[S1]}{S1}")
print(f"source discriminant raw {raw0:+.5f}  normalised {cos0:+.6f}")
moving = side_of(ADJ, Q, P)
moving = [i for i in moving if i != Q]
print(f"rotating the {len(moving)}-atom group beyond {A.symbols[Q]}{Q}: {moving}")

axis = A.positions[Q] - A.positions[P]
axis = axis / np.linalg.norm(axis)
origin = A.positions[Q].copy()
K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])

base_key = chemical_key(A)
rows = []
# ONLY near the two zero crossings. A full 360 sweep was tried first and had to be
# abandoned: rotating a 13-atom group that far brings atoms within bond_scale of each
# other, so graph_hash moved and the ladder measured connectivity as well as parity.
# The assertion below is what caught it. Near the crossings the displacement is small and
# the graph holds, which is the whole requirement for this to isolate parity.
WINDOWS = [(90.0, 6.0), (270.0, 6.0)]
degrees = sorted({round(c + d, 5) for c, half in WINDOWS
                  for d in np.concatenate([np.linspace(-half, half, 25),
                                           np.linspace(-0.25, 0.25, 51)])})
for deg in degrees:
    th = np.radians(deg)
    R = np.eye(3) + np.sin(th) * K + (1 - np.cos(th)) * (K @ K)
    cand = A.copy()
    for i in moving:
        cand.positions[i] = origin + R @ (A.positions[i] - origin)
    raw, cosv = discriminant(cand, P, Q, S0, S1)
    k = chemical_key(cand)
    rows.append(dict(deg=float(deg), discriminant_raw=raw, discriminant_cos=cosv,
                     key=k["key"][:16], graph_hash=k["graph_hash"],
                     locked_parity=k["locked_bond_parity"],
                     n_locked=len(k["locked_bond_parity"]),
                     stereo_unresolved=k["stereo_unresolved"],
                     fragments=k["fragments"]))

# Rungs whose graph moved are DROPPED, not kept with a warning: a frame that changed
# connectivity is not a parity test case, and leaving it in would let a graph change be
# scored as a parity change.
source_graph = base_key["graph_hash"]
dropped = [r for r in rows if r["graph_hash"] != source_graph]
rows = [r for r in rows if r["graph_hash"] == source_graph]
print(f"\n{len(rows)} usable rungs (graph == source); dropped {len(dropped)} whose graph moved")
if not rows:
    raise SystemExit("no rung kept its graph -- narrow the windows further")
graphs = {r["graph_hash"] for r in rows}
assert len(graphs) == 1, graphs
print("  graph constant across every kept rung -> the ladder isolates parity")

keys = {r["key"] for r in rows}
near = sorted(rows, key=lambda r: abs(r["discriminant_cos"]))[:6]
print(f"distinct chemical_key over the ladder: {len(keys)}")
print("\nclosest rungs to the discriminant's zero:")
for r in near:
    print(f"  {r['deg']:8.3f} deg  raw {r['discriminant_raw']:+.3e}  "
          f"cos {r['discriminant_cos']:+.3e}  key {r['key']}  n_locked {r['n_locked']}")
flips = [(rows[i]["deg"], rows[i + 1]["deg"]) for i in range(len(rows) - 1)
         if rows[i]["key"] != rows[i + 1]["key"]]
print(f"\nkey changes at {len(flips)} places: {flips[:6]}")
print("  Those are where a scale-dependent tolerance has to decide, and where the current")
print("  implementation decides on sign() of a vanishing number.")

json.dump(dict(source=str(SRC.relative_to(ROOT)), locked=[P, Q], substituents=[S0, S1],
               moving_group=moving, source_key=base_key["key"][:16],
               source_graph_hash=source_graph, source_discriminant_cos=cos0,
               dropped_graph_moved=len(dropped), rungs=rows),
          open(OUT / "ladder.json", "w"), indent=1, default=str)
for r, deg in zip(rows, [r["deg"] for r in rows]):
    pass
write(str(OUT / "ladder.extxyz"), [], format="extxyz")   # truncate
for r in rows:
    th = np.radians(r["deg"])
    R = np.eye(3) + np.sin(th) * K + (1 - np.cos(th)) * (K @ K)
    cand = A.copy()
    for i in moving:
        cand.positions[i] = origin + R @ (A.positions[i] - origin)
    cand.info = {"deg": r["deg"], "discriminant_cos": r["discriminant_cos"]}
    write(str(OUT / "ladder.extxyz"), cand, append=True)
print(f"\nwrote {OUT}/ladder.json and ladder.extxyz ({len(rows)} frames)")
