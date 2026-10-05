"""b05 probe3: why 15 of 86 delivered probes tripped the collateral gate.

Task 1 of the 2026-09-09 ruling. Two families are reviewed SEPARATELY -- both fall back
to `continuation`, and that alone is not evidence of a shared root cause.

Reads the run's own records for what was delivered, then re-delivers on the same relaxed
source microstate so the numbers can be attributed atom by atom. No potential is
evaluated: delivery is geometry.
"""

import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from ase.io import read  # noqa: E402

from prrs.config import SearchConfig  # noqa: E402
from prrs.internal import coordinate  # noqa: E402
from prrs.perturbations import Probe, apply as apply_probe  # noqa: E402
from prrs.state import encode  # noqa: E402

RUN = ROOT / "runs" / "b05_probe3_budget"
SOURCE = RUN / "chemical" / "c0000_m0000.extxyz"


def main():
    text = (ROOT / "benchmarks" / "b05_tertiary_halide" / "config.json").read_text()
    cfg = SearchConfig(
        **json.loads(
            "\n".join(line for line in text.splitlines() if not line.lstrip().startswith("//"))
        )
    )
    source = read(SOURCE)
    graph = encode(source, cfg.bond_scale)
    symbols = source.get_chemical_symbols()

    refused = []
    for line in (RUN / "attempts.jsonl").read_text().splitlines():
        record = json.loads(line)
        if (record.get("failure") or {}).get("code") == "collateral_bond":
            refused.append(record)
    print(f"{len(refused)} refused deliveries in {RUN.name}\n")

    groups = defaultdict(list)
    for record in refused:
        probe = record["probe"]
        kind = "bend" if probe["family"] == "bend" else "pair"
        groups[kind].append(record)

    for kind, records in sorted(groups.items()):
        print("=" * 78)
        print(f"{kind}: {len(records)} refusals")
        print("=" * 78)
        for record in records:
            spec = record["probe"]
            probe = Probe(
                spec["family"], tuple(spec["indices"]), spec["sign"], spec["amplitude"], 0
            )
            moved = source.copy()
            delivered = apply_probe(moved, probe, cfg.bond_scale, cfg.active_atoms)
            collateral = delivered["collateral"]
            worst = collateral["worst_bond"]
            displacement = np.linalg.norm(moved.positions - source.positions, axis=1)
            movers = np.argsort(displacement)[::-1][:4]
            target = coordinate(
                source.positions,
                probe.kind,
                probe.indices,
            )
            print(
                f"\n  {spec['family']} {tuple(spec['indices'])} amp={spec['amplitude']} "
                f"({probe.unit})  realization={delivered['realization']}"
            )
            print(
                f"     target {target:.4f} -> requested {delivered['requested']:+.4f}, "
                f"achieved {delivered['achieved']:+.4f}, error {delivered['target_error']:.2e}"
            )
            print(
                f"     gate tripped on bond {worst} "
                f"({symbols[worst[0]]}{worst[0]}-{symbols[worst[1]]}{worst[1]}) "
                f"changed {collateral['max_bond_change_A']:.4f} A "
                f"(threshold {cfg.collateral_bond_A})"
            )
            print(
                "     largest displacements: "
                + ", ".join(
                    f"{symbols[i]}{i} {displacement[i]:.3f}"
                    for i in movers
                    if displacement[i] > 1e-6
                )
            )
            # is the tripped bond adjacent to the driven atoms, or elsewhere?
            driven = set(probe.indices)
            touching = driven & set(worst)
            print(
                f"     tripped bond shares atoms with the probe: "
                f"{sorted(touching) if touching else 'NO -- collateral is remote'}"
            )
    print("\nedges:", sorted(graph.edges))


def rigid_bend(positions, i, j, k, delta, adjacency):
    """The counterfactual for a bend: rotate the two side fragments about the vertex.

    Every bond length inside each fragment is preserved because a fragment moves as a
    rigid body, and |ji| and |jk| are preserved because the rotation axis passes through
    j. Only the angle, and torsions across the vertex, change. This is the exact analogue
    of internal.rotate_fragment for an angle, and it is the counterfactual the ruling asks
    for: if it reaches the same amplitude and passes the gate, the refusal was the
    delivery's fault, not the gate's.

    Returns None where no rigid split exists (i and k still connected without j).
    """
    positions = np.array(positions, dtype=float)

    def reach(start, blocked):
        seen, stack = {start}, [start]
        while stack:
            node = stack.pop()
            for neighbour in adjacency[node]:
                if neighbour == blocked or neighbour in seen:
                    continue
                seen.add(neighbour)
                stack.append(neighbour)
        return seen

    side_i, side_k = reach(i, j), reach(k, j)
    if k in side_i or side_i & side_k:
        return None
    u = positions[i] - positions[j]
    v = positions[k] - positions[j]
    axis = np.cross(u, v)
    norm = np.linalg.norm(axis)
    if norm < 1e-12:
        return None
    axis = axis / norm

    def rotate(members, angle):
        offsets = positions[list(members)] - positions[j]
        along = np.outer(offsets @ axis, axis)
        across = offsets - along
        perpendicular = np.cross(axis, across)
        positions[list(members)] = (
            positions[j] + along + across * np.cos(angle) + perpendicular * np.sin(angle)
        )

    # opening the angle means rotating the two sides apart by delta/2 each
    rotate(side_i, -delta / 2)
    rotate(side_k, +delta / 2)
    return positions


def counterfactual():
    """Task 2: can a bond-length-preserving delivery reach the same amplitude?"""
    from prrs.chemistry import canonical_labels  # noqa: F401 - kept for parity of imports
    from prrs.perturbations import collateral_metrics

    text = (ROOT / "benchmarks" / "b05_tertiary_halide" / "config.json").read_text()
    cfg = SearchConfig(
        **json.loads(
            "\n".join(line for line in text.splitlines() if not line.lstrip().startswith("//"))
        )
    )
    source = read(SOURCE)
    graph = encode(source, cfg.bond_scale)
    adjacency = {int(a): set() for a in graph.index}
    for a, b in graph.edges:
        adjacency[int(a)].add(int(b))
        adjacency[int(b)].add(int(a))

    print("\n" + "=" * 78)
    print("TASK 2a  bend: rigid two-sided rotation about the vertex")
    print("=" * 78)
    for indices, amplitudes in (((1, 0, 10), (0.17, 0.35)), ((10, 0, 6), (-0.17, -0.35))):
        for delta in amplitudes:
            moved = source.copy()
            positions = rigid_bend(source.positions, *indices, delta, adjacency)
            if positions is None:
                print(f"  bend {indices} {delta:+.2f}: no rigid split")
                continue
            moved.set_positions(positions)
            achieved = coordinate(moved.positions, "angle", indices) - coordinate(
                source.positions, "angle", indices
            )
            metrics = collateral_metrics(source, moved, "angle", indices, cfg.bond_scale)
            print(
                f"  bend {indices} requested {delta:+.3f} -> achieved {achieved:+.6f} rad "
                f"(error {abs(achieved - delta):.1e}); worst non-target bond "
                f"{metrics['max_bond_change_A']:.2e} A on {metrics['worst_bond']} "
                f"-> gate {'PASSES' if metrics['max_bond_change_A'] < cfg.collateral_bond_A else 'still refuses'}"
            )

    print("\n" + "=" * 78)
    print("TASK 2b  same-fragment H...H: what is reachable with bonds and angles fixed?")
    print("=" * 78)
    from prrs.internal import rotate_fragment

    for pair in ((5, 8), (9, 12), (8, 13), (8, 12), (9, 11)):
        a, b = pair
        carbon_a = next(iter(adjacency[a]))
        carbon_b = next(iter(adjacency[b]))
        if carbon_a == carbon_b:
            print(f"  {pair}: both on the same methyl -- rigid rotation cannot separate them")
            continue
        movers_a = [h for h in adjacency[carbon_a] if h != 0]
        movers_b = [h for h in adjacency[carbon_b] if h != 0]
        start = source.get_distance(a, b)
        reachable = []
        for phi in np.linspace(0, 2 * np.pi, 37):
            first = rotate_fragment(source.positions, 0, carbon_a, movers_a, phi)
            for psi in np.linspace(0, 2 * np.pi, 37):
                second = rotate_fragment(first, 0, carbon_b, movers_b, psi)
                reachable.append(float(np.linalg.norm(second[a] - second[b])))
        print(
            f"  {pair}: r = {start:.4f} A now; rigid methyl rotations reach "
            f"[{min(reachable):.4f}, {max(reachable):.4f}] A"
        )
        for requested in (0.3, 0.6, 1.0):
            target = start - requested
            verdict = (
                "reachable" if min(reachable) <= target <= max(reachable) else "NOT reachable"
            )
            print(f"       compress {requested} -> target {target:.4f} A: {verdict}")


if __name__ == "__main__":
    main()
    sys.exit(counterfactual())
