"""Measure colour-vs-orbit disagreement before changing any consumer (PLAN section 1.3).

WL colours are invariant but not complete: colour equality is a GUESS that two atoms
are equivalent, and the guess is wrong exactly when a true stereo centre's neighbours
share a colour yet no automorphism swaps them (P1-3), or when two graph-equivalent
bridge bonds are geometrically distinct in one fixed conformation (P1-7). The fix is
to consume certified orbits (canonical_form) instead of colours -- but whether that
changes ANY verdict on the recorded corpus is an empirical question, and decision 2
(2026-09-20) hinges on the answer: zero disagreement means the change is free
(digest stays zero-diff); nonzero means it is an identity tightening that must NOT
touch chemical_key now.

Four comparisons, and the fourth is not a variant of the third. Comparisons 1-3 all
ask one question -- does replacing a COLOUR by a certified GRAPH orbit move a verdict
-- and they are what decision 2 needed. But `propose` does not deduplicate on graph
orbits: it deduplicates on the GEOMETRIC subgroup (the automorphisms that also
preserve the current geometry), which is strictly smaller and therefore yields
strictly more candidates. That is the pair P1-7 is actually about, it is the pair
that shipped, and comparison 3 says nothing about it -- a zero there was quoted once
as though it did (`perturbations.propose`, corrected 2026-09-21). Comparison 4
measures it.

Usage:  python docs/experiments/orbit_vs_colour_probe.py
Prints a human summary and a final machine-readable JSON line.
"""

import glob
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

from ase.io import read  # noqa: E402

from prrs.chemistry import (  # noqa: E402
    automorphisms,
    canonical_labels,
    locked_edges,
)
from prrs.canonical import canonical_form  # noqa: E402
from prrs.perturbations import (  # noqa: E402
    FAMILIES,
    _adjacency,
    _angle_directions,
    _bridges,
    _dihedral_directions,
    geometric_automorphisms,
    orbits_from_group,
)
from prrs.state import encode  # noqa: E402


def corpus():
    seen = set()
    for pat in (
        "runs/p3_*.extxyz",
        "runs/p3_three_model/*.extxyz",
        "runs/p3_landscape_full/conf*.extxyz",
        "runs/p3_paths/*/*/*.extxyz",
        "runs/p1_round8/chemical/*.extxyz",
        "runs/p2_round7/chemical/*.extxyz",
    ):
        for f in sorted(glob.glob(str(ROOT / pat))):
            rel = str(Path(f).relative_to(ROOT))
            if rel not in seen:
                seen.add(rel)
                yield rel, f


def stabiliser_orbits(neighbours, group, slot, atom):
    """Orbits of `neighbours` under the subgroup fixing `atom` (Stab(atom))."""
    stab = [g for g in group if int(g[slot[atom]]) == atom]
    remaining, orbits = list(neighbours), []
    while remaining:
        seed = remaining.pop(0)
        orbit = {seed}
        frontier = [seed]
        while frontier:
            current = frontier.pop()
            for g in stab:
                image = int(g[slot[current]])
                if image not in orbit and image in remaining:
                    orbit.add(image)
                    remaining.remove(image)
                    frontier.append(image)
        orbits.append(sorted(orbit))
    return orbits


def analyse(atoms):
    """Four comparisons on one frame.

    Returns (skipped, stereo_diff, ez_diff, direction_diff, geometry) where `geometry`
    is None when skipped and otherwise (group_size, geometric_size, graph_candidates,
    geometric_candidates) -- counts, not a bool, because the fourth comparison's answer
    is a magnitude and reporting it as "differs / does not" is what hid it the first
    time.
    """
    graph = encode(atoms)
    index = graph.index
    labels = canonical_labels(atoms.numbers, graph.edges, index)
    adjacency = _adjacency(graph.edges, index)
    group, info = automorphisms(atoms, 1.2, None, 20000)
    if not info.get("enumerated"):
        return "unenumerated", None, None, None, None
    order = [int(a) for a in index]
    slot = {atom: k for k, atom in enumerate(order)}
    form = canonical_form(atoms.numbers, graph.edges, index)
    if not form.info.get("enumerated"):
        return "unenumerated", None, None, None, None
    orbits = form.orbits

    # 1. stereo centres: four distinct neighbour colours vs four distinct Stab-orbits
    stereo_diff = False
    for atom, neighbours in adjacency.items():
        if len(neighbours) != 4:
            continue
        colour_distinct = len({labels[n] for n in neighbours}) == 4
        stab_orbits = stabiliser_orbits(list(neighbours), group, slot, int(atom))
        orbit_distinct = len(stab_orbits) == 4
        if colour_distinct != orbit_distinct:
            stereo_diff = True

    # 2. E/Z substituent pick: colour max-label uniqueness vs orbit uniqueness
    ez_diff = False
    locked = locked_edges(atoms.numbers, graph.edges, index)
    for a, b in sorted(graph.edges):
        if (a, b) not in locked:
            continue
        for near, far in ((a, b), (b, a)):
            options = [n for n in adjacency[near] if n != far]
            if not options:
                continue
            best = max(labels[n] for n in options)
            colour_unique = sum(1 for n in options if labels[n] == best) == 1
            best_orbit = max(orbits[n] for n in options)
            orbit_unique = sum(1 for n in options if orbits[n] == best_orbit) == 1
            if colour_unique != orbit_unique:
                ez_diff = True

    # 3. direction dedup: colour keys vs orbit keys, candidate counts
    angle_colour, angle_orbit = set(), set()
    for vertex, neighbours in adjacency.items():
        for x in range(len(neighbours)):
            for y in range(x + 1, len(neighbours)):
                i, k = neighbours[x], neighbours[y]
                angle_colour.add((labels[vertex], tuple(sorted((labels[i], labels[k])))))
                angle_orbit.add((orbits[vertex], tuple(sorted((orbits[i], orbits[k])))))
    torsion_colour, torsion_orbit = set(), set()
    for j, k in sorted(_bridges(graph.edges, index)):
        left = [n for n in adjacency[j] if n != k]
        right = [n for n in adjacency[k] if n != j]
        if not left or not right:
            continue
        torsion_colour.add(tuple(sorted((labels[j], labels[k]))))
        torsion_orbit.add(tuple(sorted((orbits[j], orbits[k]))))
    direction_diff = len(angle_colour) != len(angle_orbit) or len(torsion_colour) != len(
        torsion_orbit
    )

    # 4. the pair that SHIPPED: graph orbits vs GEOMETRIC orbits (P1-7). Counted the
    #    way propose() counts, through the same two enumerators, so the number can be
    #    quoted in its comment without a second convention to reconcile.
    angle_families = tuple(f for f, spec in FAMILIES.items() if spec["kind"] == "angle")
    dihedral_families = tuple(f for f, spec in FAMILIES.items() if spec["kind"] == "dihedral")

    def candidates(basis):
        return len(list(_angle_directions(graph, basis, angle_families))) + len(
            list(_dihedral_directions(graph, labels, basis, dihedral_families))
        )

    subgroup = geometric_automorphisms(atoms, form.group, None)
    geometric = (
        len(form.group),
        len(subgroup),
        candidates(orbits),
        candidates(orbits_from_group(subgroup, index)),
    )
    return False, stereo_diff, ez_diff, direction_diff, geometric


def main():
    frames = stereo = ez = direction = unenum = 0
    shrunk = widened = 0
    graph_total = geometric_total = 0
    examples = {"stereo": [], "ez": [], "direction": [], "geometric": []}
    for rel, path in corpus():
        for i, atoms in enumerate(read(path, index=":")):
            frames += 1
            skipped, s, z, d, geometry = analyse(atoms)
            if skipped:
                unenum += 1
                continue
            if s:
                stereo += 1
                if len(examples["stereo"]) < 3:
                    examples["stereo"].append(f"{rel}#{i}")
            if z:
                ez += 1
                if len(examples["ez"]) < 3:
                    examples["ez"].append(f"{rel}#{i}")
            if d:
                direction += 1
                if len(examples["direction"]) < 3:
                    examples["direction"].append(f"{rel}#{i}")
            group_size, geometric_size, graph_n, geometric_n = geometry
            graph_total += graph_n
            geometric_total += geometric_n
            if geometric_size < group_size:
                shrunk += 1
            if geometric_n != graph_n:
                widened += 1
                if len(examples["geometric"]) < 3:
                    examples["geometric"].append(
                        f"{rel}#{i} |Aut|={group_size}->{geometric_size} "
                        f"cand {graph_n}->{geometric_n}"
                    )
    print(f"frames {frames} (unenumerated {unenum})")
    print(f"stereo-centre verdicts differ on   {stereo} frames  e.g. {examples['stereo']}")
    print(f"E/Z substituent verdicts differ on {ez} frames  e.g. {examples['ez']}")
    print(
        f"direction dedup counts differ on   {direction} frames  e.g. {examples['direction']}"
    )
    print(
        f"[4] geometry breaks a graph symmetry on {shrunk} frames; candidate count "
        f"differs on {widened} frames"
    )
    print(
        f"    angle+dihedral candidates summed over the corpus: "
        f"graph orbits {graph_total} -> geometric orbits {geometric_total}"
        + (
            f"  (+{100 * (geometric_total - graph_total) / graph_total:.1f}%)"
            if graph_total
            else ""
        )
    )
    for line in examples["geometric"]:
        print(f"    e.g. {line}")
    verdict = (
        "[1-3] zero disagreement: the colour -> graph-orbit change is free "
        "(digest stays zero-diff)"
        if stereo == 0 and ez == 0 and direction == 0
        else "[1-3] NONZERO disagreement: identity tightening -- decision 2 applies, "
        "chemical_key untouched"
    )
    print(verdict)
    print(
        "[4] the geometric basis is a REAL widening and comparisons 1-3 cannot see it: "
        "they compare colours with graph orbits, propose() deduplicates on the "
        "geometric subgroup"
    )
    print(
        json.dumps(
            {
                "frames": frames,
                "stereo_diff_frames": stereo,
                "ez_diff_frames": ez,
                "direction_diff_frames": direction,
                "geometry_shrinks_group_frames": shrunk,
                "geometric_widens_candidates_frames": widened,
                "candidates_graph_orbits": graph_total,
                "candidates_geometric_orbits": geometric_total,
            }
        )
    )


if __name__ == "__main__":
    main()
