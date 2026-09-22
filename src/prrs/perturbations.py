"""Structured probes spanning the primitive internal-coordinate subspaces.

The proposal basis must cover stretch, bend and torsion, not only Cartesian pair
directions. Ethanol showed why: its only sub-eV motion is a torsion, so a search built
from pair displacements and pair kicks returns elastic at every amplitude while three
real basins sit a few meV away.

Directions are deduplicated by canonical atom colour. Three equivalent hydrogens on a
methyl group give one angle direction rather than three, and one torsion is proposed
per rotatable bond rather than one per dihedral about it, because all dihedrals about a
bond differ by a constant. Spending the direction budget on symmetry copies buys
nothing.
"""

from collections import deque
from dataclasses import asdict, dataclass
import math
import numpy as np
from . import internal
from .chemistry import canonical_labels, _bridges
from .state import encode

# mode: how the probe acts. sign: fixed by the family, or None when both are proposed.
FAMILIES = {
    "stretch": {
        "kind": "bond",
        "mode": "displace",
        "sign": 1,
        "amplitudes": "geometry_amplitudes_A",
        "bonded": True,
    },
    "compress": {
        "kind": "bond",
        "mode": "displace",
        "sign": -1,
        "amplitudes": "geometry_amplitudes_A",
        "bonded": False,
    },
    "kick": {
        "kind": "bond",
        "mode": "kick",
        "sign": None,
        "amplitudes": "kick_energies_eV",
        "bonded": None,
    },
    "pulse": {
        "kind": "bond",
        "mode": "pulse",
        "sign": None,
        "amplitudes": "pulse_forces_eV_A",
        "bonded": None,
    },
    "bend": {
        "kind": "angle",
        "mode": "displace",
        "sign": None,
        "amplitudes": "bend_amplitudes_rad",
        "bonded": None,
    },
    "bend_kick": {
        "kind": "angle",
        "mode": "kick",
        "sign": None,
        "amplitudes": "kick_energies_eV",
        "bonded": None,
    },
    "torsion": {
        "kind": "dihedral",
        "mode": "displace",
        "sign": None,
        "amplitudes": "torsion_amplitudes_rad",
        "bonded": None,
    },
    "torsion_kick": {
        "kind": "dihedral",
        "mode": "kick",
        "sign": None,
        "amplitudes": "kick_energies_eV",
        "bonded": None,
    },
}
UNITS = {
    "displace": {"bond": "Angstrom", "angle": "radian", "dihedral": "radian"},
    "kick": "eV",
    "pulse": "eV/Angstrom",
}


@dataclass(frozen=True)
class Probe:
    family: str
    indices: tuple[int, ...]
    sign: int
    amplitude: float
    seed: int
    # Why this direction was put in the front band, or None if it was not. Recorded on the
    # probe so a prior can never become hidden logic: whatever moved a direction forward is
    # in the trial record next to the result it produced.
    priority: int = 1
    priority_reason: dict = None

    def __post_init__(self):
        if self.family not in FAMILIES:
            raise ValueError(f"unknown perturbation family {self.family!r}")
        expected = internal.KINDS[FAMILIES[self.family]["kind"]]
        # internal.atom_index rather than int(): a bare int() turns 1.9 into 1 and -1 into
        # a numpy index onto the last atom, so a probe would silently act on an atom nobody
        # named. The upper bound cannot be checked here -- a Probe is built before it meets
        # a molecule -- so this catches the type and the sign only.
        index = tuple(internal.atom_index(i) for i in self.indices)
        if len(index) != expected or len(set(index)) != expected:
            raise ValueError(f"{self.family} needs {expected} distinct atom indices")
        object.__setattr__(self, "indices", index)

    @property
    def kind(self):
        return FAMILIES[self.family]["kind"]

    @property
    def mode(self):
        return FAMILIES[self.family]["mode"]

    @property
    def unit(self):
        mode = self.mode
        return UNITS[mode][self.kind] if mode == "displace" else UNITS[mode]

    @property
    def pair(self):
        if len(self.indices) != 2:
            raise AttributeError(f"{self.family} is not a pair probe")
        return self.indices

    def to_dict(self):
        return {
            **asdict(self),
            "indices": list(self.indices),
            "unit": self.unit,
            "kind": self.kind,
            "mode": self.mode,
        }


def pair_axis(atoms, pair):
    i, j = pair
    vector = atoms.positions[j] - atoms.positions[i]
    distance = np.linalg.norm(vector)
    if distance < 1e-12:
        raise ValueError("Pair axis undefined for coincident atoms")
    return vector / distance, float(distance)


def displace_pair(atoms, pair, delta_A):
    """Kept as the bond special case of the general internal-coordinate displacement."""
    internal.displace(atoms, "bond", tuple(pair), delta_A)


def kick_pair(atoms, pair, energy_eV, sign=1):
    return internal.kick(atoms, "bond", tuple(pair), energy_eV, sign)


def torsion_for(atoms, indices, bond_scale=1.2, active=None):
    """Look up the torsion descriptor for a dihedral at the current geometry."""
    from .canonical import canonical_form

    graph = encode(atoms, bond_scale, active=active)
    form = canonical_form(atoms.numbers, graph.edges, graph.index)
    labels = form.colours
    orbits = form.orbits if form.info.get("enumerated") else labels
    genuine, rotors = rotatable_torsions(graph, labels, orbits)
    target = tuple(indices)
    for torsion in genuine + rotors:
        if (
            torsion.indices == target
            or torsion.indices[1:3] == target[1:3]
            or torsion.indices[1:3] == target[2:0:-1]
        ):
            return torsion
    return None


def _graph_angles(graph):
    adjacency = _adjacency(graph.edges, graph.index)
    for vertex, neighbours in adjacency.items():
        for a in range(len(neighbours)):
            for b in range(a + 1, len(neighbours)):
                yield (neighbours[a], vertex, neighbours[b])


def collateral_metrics(before, after, kind, indices, bond_scale=1.2, active=None):
    """What a perturbation did to the coordinates it was not aiming at.

    Target-coordinate fidelity is not geometric fidelity. A torsion probe once reached
    its requested dihedral to one part in a million while stretching an O-H bond from
    0.957 to 2.346 Angstrom, and every check in place at the time passed. Acceptance has
    to be two-dimensional: the target coordinate, and the damage everywhere else.

    Bonds, clashes and fragment count are the "did it tear the system apart" measures and
    are meant to be gated. Non-target angles and torsions are reported but not gated:
    they couple to the target legitimately -- angles at a shared vertex are constrained to
    sum, so a real bend probe must move its neighbours.
    """
    from ase.data import covalent_radii

    graph = encode(before, bond_scale, active=active)
    target_bond = {indices[0], indices[1]} if kind == "bond" else None
    target_angle = tuple(indices) if kind == "angle" else None

    worst_bond, worst_bond_pair = 0.0, None
    for a, b in sorted(graph.edges):
        if target_bond == {a, b}:
            continue
        change = abs(after.get_distance(a, b) - before.get_distance(a, b))
        if change > worst_bond:
            worst_bond, worst_bond_pair = change, (a, b)

    worst_angle, worst_angle_triple = 0.0, None
    for triple in _graph_angles(graph):
        if target_angle in (triple, triple[::-1]):
            continue
        change = abs(
            internal.coordinate(after.positions, "angle", triple)
            - internal.coordinate(before.positions, "angle", triple)
        )
        if change > worst_angle:
            worst_angle, worst_angle_triple = change, triple

    from .canonical import canonical_form

    form = canonical_form(before.numbers, graph.edges, graph.index)
    labels = form.colours
    orbits = form.orbits if form.info.get("enumerated") else labels
    genuine, rotors = rotatable_torsions(graph, labels, orbits)
    worst_torsion, worst_torsion_indices = 0.0, None
    for torsion in genuine + rotors:
        if kind == "dihedral" and torsion.indices[1:3] in (
            tuple(indices[1:3]),
            tuple(indices[2:0:-1]),
        ):
            continue
        delta = internal.coordinate(
            after.positions, "dihedral", torsion.indices
        ) - internal.coordinate(before.positions, "dihedral", torsion.indices)
        change = abs((delta + np.pi) % (2 * np.pi) - np.pi)
        if change > worst_torsion:
            worst_torsion, worst_torsion_indices = change, torsion.indices

    radii = covalent_radii[after.numbers]
    distances = after.get_all_distances()
    np.fill_diagonal(distances, np.inf)
    ratio = distances / (radii[:, None] + radii[None, :])
    closest = int(np.argmin(ratio))
    return {
        "max_bond_change_A": worst_bond,
        "worst_bond": list(worst_bond_pair) if worst_bond_pair else None,
        "max_nontarget_angle_change_rad": worst_angle,
        "worst_angle": list(worst_angle_triple) if worst_angle_triple else None,
        "max_nontarget_torsion_change_rad": worst_torsion,
        "worst_torsion": list(worst_torsion_indices) if worst_torsion_indices else None,
        "closest_distance_ratio": float(ratio.min()),
        "closest_pair": [closest // len(after), closest % len(after)],
        "fragments_before": graph.components,
        "fragments_after": _collateral_fragments(
            graph, after, kind, indices, bond_scale, active
        ),
        "gated": ["max_bond_change_A", "closest_distance_ratio", "fragments_after"],
        "reported_only": ["max_nontarget_angle_change_rad", "max_nontarget_torsion_change_rad"],
    }


def _collateral_fragments(before_graph, after, kind, indices, bond_scale, active):
    """Fragment count after the probe, with the target bond held at its original state.

    A stretch probe that breaks the bond it is pulling has done exactly what it was asked
    to do; counting that as collateral damage would forbid the family outright. Only
    fragmentation the probe caused somewhere else is collateral.
    """
    from .state import Graph

    graph = encode(after, bond_scale, active=active)
    edges = set(graph.edges)
    if kind == "bond":
        pair = (min(indices[:2]), max(indices[:2]))
        if pair in before_graph.edges:
            edges.add(pair)
        else:
            edges.discard(pair)
    return Graph(edges, graph.weights, graph.distances, graph.index).components


def torsion_on_bond(graph, labels, bond, orbits=None):
    """The rotatable torsion whose central bond is `bond`, or None.

    A coordinate is identified by its bond, not by the four atoms: all dihedrals about one
    bond differ by a constant, and which representative `rotatable_torsions` picks depends
    only on the canonical colouring, so two structures of the same molecule agree on it.
    Looking a coordinate up this way is what lets one structure's measured stiffness be
    applied to another's geometry.
    """
    wanted = frozenset(bond)
    genuine, rotors = rotatable_torsions(graph, labels, orbits)
    for torsion in genuine + rotors:
        if frozenset(torsion.indices[1:3]) == wanted:
            return torsion
    return None


def step_torsion(atoms, torsion, delta, in_place=True):
    """Move along one torsion's own exact path, by delta radians.

    Uses the rigid rotation where the bond is a bridge, so no other internal coordinate
    changes at all, and the continuation scheme inside a ring. Having the coordinate's own
    path available is what lets its generalized force and stiffness be measured by
    differencing the energy along it, with no chain rule and no ambiguity about how the
    remaining coordinates were held.
    """
    target = atoms if in_place else atoms.copy()
    before = internal.coordinate(target.positions, "dihedral", torsion.indices)
    if torsion.moving:
        j, k = torsion.indices[1], torsion.indices[2]
        target.set_positions(
            internal.rotate_fragment(
                target.positions, j, k, torsion.moving, -delta, target.get_masses()
            )
        )
    else:
        internal.displace(target, "dihedral", torsion.indices, delta)
    after = internal.coordinate(target.positions, "dihedral", torsion.indices)
    return target, float((after - before + np.pi) % (2 * np.pi) - np.pi)


@dataclass(frozen=True)
class FragmentRotation:
    """A fragment's orientation relative to the rest of the system, as one coordinate.

    Same shape as `Torsion` on purpose, so the soft-mode polish can treat the two the
    same: an angle in radians, periodic, with a symmetry order bounding how much of the
    circle is distinct. `symmetry_order` is 1 here rather than measured -- a conservative
    cap of pi, which is what a rotation with no established symmetry is entitled to.
    """

    members: tuple
    axis: tuple
    label: str
    symmetry_order: int = 1

    @property
    def indices(self):
        """What the report keys on. Torsions report four atoms; this reports its fragment."""
        return tuple(self.members)

    @property
    def domain(self):
        return 2 * np.pi / self.symmetry_order

    def to_dict(self):
        return {
            "members": list(self.members),
            "axis": list(self.axis),
            "label": self.label,
            "symmetry_order": self.symmetry_order,
            "fundamental_domain_rad": self.domain,
            "realization": "rigid_rotation_about_fragment_com",
        }


def fragment_rotations(atoms, graph):
    """Each fragment's three rigid rotations, about its own principal axes.

    Only when there is more than one fragment: a single fragment's rotation is a trivial
    mode, already projected out of every spectrum and carrying no energy by construction.
    With two or more, the RELATIVE orientation is a real soft coordinate that no bond of
    the graph can reach -- see `internal.rotate_rigid` for the measurement that made this
    necessary.

    Principal axes rather than the lab frame, because the lab frame is not a property of
    the structure: the same geometry written down after a rotation would get a different
    set of coordinates and a different polish. Principal axes turn with the molecule, so
    the coordinates do too. A symmetric top makes two of them degenerate, and then the
    pair is an arbitrary basis of the same plane -- harmless for a relaxation, which only
    needs to span the space, and recorded rather than hidden by `degenerate_inertia`.

    A monatomic fragment (a bare halide) has no orientation and contributes nothing.
    """
    from .chemistry import _fragments

    fragments = [sorted(int(a) for a in part) for part in _fragments(graph.edges, graph.index)]
    if len(fragments) < 2:
        return ()
    fragments.sort()  # deterministic order; union-find hands them back in hash order
    masses = atoms.get_masses()
    modes = []
    for members in fragments:
        if len(members) < 2:
            continue
        block = atoms.positions[members]
        weights = masses[members]
        centre = np.average(block, axis=0, weights=weights)
        offsets = block - centre
        inertia = np.einsum(
            "i,ijk->jk",
            weights,
            np.eye(3) * (offsets**2).sum(axis=1)[:, None, None]
            - offsets[:, :, None] * offsets[:, None, :],
        )
        values, axes = np.linalg.eigh(inertia)
        for column in range(3):
            modes.append(
                FragmentRotation(
                    members=tuple(members),
                    axis=tuple(float(v) for v in axes[:, column]),
                    label=f"fragment{members[0]}_axis{column}",
                )
            )
    return tuple(modes)


def step_fragment_rotation(atoms, mode, delta, in_place=True):
    """Move along one fragment rotation's own exact path, by delta radians."""
    target = atoms if in_place else atoms.copy()
    target.set_positions(
        internal.rotate_rigid(
            target.positions, mode.members, mode.axis, delta, target.get_masses()
        )
    )
    return target, float(delta)


def apply(atoms, probe, bond_scale=1.2, active=None):
    """Realize a probe, preferring an exact realization, and report the collateral.

    Each coordinate is delivered by the exact rigid motion that moves it and nothing
    else, where one exists: a torsion about a bridge bond as a rotation of the fragment on
    one side, a bond as a translation of its two sides along the axis, an angle as
    opposite rotations of its two sides about the vertex. All three preserve every other
    bond length and bond angle to machine precision. Only where no rigid split exists --
    anything inside a ring, a non-bonded pair whose ends stay connected, or any other
    genuinely collective coordinate -- does delivery fall back to the continuation scheme:
    a large perturbation realized as a sequence of small, locally valid ones.
    """
    if probe.mode == "kick":
        before = atoms.get_kinetic_energy()
        internal.kick(atoms, probe.kind, probe.indices, probe.amplitude, probe.sign)
        return {
            "requested": probe.amplitude,
            "achieved": atoms.get_kinetic_energy() - before,
            "unit": probe.unit,
            "realization": "internal_coordinate_impulse",
        }
    if probe.mode != "displace":
        raise ValueError(f"{probe.family} is applied by the trial runner, not here")

    reference = atoms.copy()
    requested = probe.sign * probe.amplitude
    realization = "continuation"
    torsion = (
        torsion_for(atoms, probe.indices, bond_scale, active)
        if probe.kind == "dihedral"
        else None
    )
    if torsion is not None and torsion.moving:
        realization = "rigid_rotation"
        _, _, k, _ = torsion.indices
        j = torsion.indices[1]
        # A right-handed rotation of the l-side about j->k decreases the dihedral in the
        # convention coordinate() follows (ASE's), so the requested change is negated.
        # The paired test pins this rather than trusting the reasoning.
        atoms.set_positions(
            internal.rotate_fragment(
                atoms.positions, j, k, torsion.moving, -requested, atoms.get_masses()
            )
        )
        achieved = internal.coordinate(
            atoms.positions, "dihedral", probe.indices
        ) - internal.coordinate(reference.positions, "dihedral", probe.indices)
        achieved = float((achieved + np.pi) % (2 * np.pi) - np.pi)
    elif probe.kind == "angle" and (
        sides := _angle_sides(encode(atoms, bond_scale, active=active), *probe.indices)
    ):
        realization = "rigid_bend"
        atoms.set_positions(
            internal.rotate_about_vertex(
                atoms.positions, *probe.indices, *sides, requested, atoms.get_masses()
            )
        )
        achieved = internal.coordinate(
            atoms.positions, "angle", probe.indices
        ) - internal.coordinate(reference.positions, "angle", probe.indices)
    elif probe.kind == "bond" and (
        sides := _bond_sides(encode(atoms, bond_scale, active=active), *probe.indices)
    ):
        realization = "rigid_translation"
        atoms.set_positions(
            internal.translate_fragments(
                atoms.positions, *probe.indices, *sides, requested, atoms.get_masses()
            )
        )
        achieved = internal.coordinate(
            atoms.positions, "bond", probe.indices
        ) - internal.coordinate(reference.positions, "bond", probe.indices)
    else:
        _, achieved = internal.displace(atoms, probe.kind, probe.indices, requested)

    return {
        "requested": requested,
        "achieved": achieved,
        "unit": probe.unit,
        "realization": realization,
        "target_error": abs(achieved - requested),
        "collateral": collateral_metrics(
            reference, atoms, probe.kind, probe.indices, bond_scale, active
        ),
    }


def _adjacency(edges, index):
    table = {int(i): [] for i in index}
    for a, b in edges:
        table[a].append(b)
        table[b].append(a)
    return table


def _hop_counts(graph):
    """Bond-graph distance between every pair, for the pair-approach feasibility test."""
    adjacency = {int(i): set() for i in graph.index}
    for a, b in graph.edges:
        adjacency[a].add(b)
        adjacency[b].add(a)
    hops = {}
    for start in adjacency:
        seen = {start: 0}
        queue = deque([start])
        while queue:
            node = queue.popleft()
            for neighbour in adjacency[node]:
                if neighbour not in seen:
                    seen[neighbour] = seen[node] + 1
                    queue.append(neighbour)
        for other, distance in seen.items():
            hops[(min(start, other), max(start, other))] = distance
    return hops


ACCEPTORS = (7, 8, 9)  # N, O, F


def _hydrogen_bonded_donors(atoms, graph, config):
    """X-H bonds whose hydrogen already has an acceptor in reach, from the source geometry.

    A proton transfer's reaction coordinate is one particular bond, and nothing in the
    proposal distinguished it: directions inside a family were drawn uniformly, so on
    3-oxobutanal the O-H that carries the reaction had a 3-in-11 chance of being sampled at
    all and was not. It is however distinguished by the reactant's own geometry -- its
    hydrogen sits 1.69 A from a second oxygen at an angle of nearly 180 degrees -- and that
    is readable without knowing what the product is, which is the property that matters:
    the prior says "this bond is unusual", not "this reaction happens".

    Returns {(i, j): detail}. Used to order candidates, never to discard them.
    """
    numbers = atoms.numbers
    positions = atoms.positions
    donors = {}
    for a, b in graph.edges:
        for donor, hydrogen in ((a, b), (b, a)):
            if int(numbers[hydrogen]) != 1 or int(numbers[donor]) not in ACCEPTORS:
                continue
            best = None
            for acceptor in graph.index:
                acceptor = int(acceptor)
                if acceptor in (donor, hydrogen) or int(numbers[acceptor]) not in ACCEPTORS:
                    continue
                if (min(acceptor, hydrogen), max(acceptor, hydrogen)) in graph.edges:
                    continue
                reach = float(np.linalg.norm(positions[acceptor] - positions[hydrogen]))
                if reach > config.hbond_acceptor_A:
                    continue
                first = positions[donor] - positions[hydrogen]
                second = positions[acceptor] - positions[hydrogen]
                cosine = float(
                    np.dot(first, second) / (np.linalg.norm(first) * np.linalg.norm(second))
                )
                angle = float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))
                if angle < config.hbond_angle_deg:
                    continue
                if best is None or reach < best["h_acceptor_A"]:
                    best = {
                        "reason": "hydrogen_bond_donor",
                        "donor": donor,
                        "hydrogen": hydrogen,
                        "acceptor": acceptor,
                        "h_acceptor_A": reach,
                        "donor_h_acceptor_deg": angle,
                    }
            if best is not None:
                donors[(min(donor, hydrogen), max(donor, hydrogen))] = best
    return donors


def _bond_directions(atoms, graph, config, families):
    """Stretch and compress proposals along existing bonds.

    Notes
    -----
    [1] A pair two or three bonds apart is already held by the covalent skeleton, so
        compressing it is a bend or a torsion in disguise and the continuation has to drag
        the connecting bonds to deliver it. Measured on 3-oxobutanal, every pair-approach
        the collateral gate refused was three hops or fewer, and they were 45% of the trial
        budget. Refusing them here spends the budget on directions that can actually be
        delivered.
    """
    donors = _hydrogen_bonded_donors(atoms, graph, config)
    hops = _hop_counts(graph)
    for a in range(len(graph.index)):
        for b in range(a + 1, len(graph.index)):
            i, j = int(graph.index[a]), int(graph.index[b])
            bonded = (i, j) in graph.edges
            if not bonded:
                if graph.distances[a, b] > config.pair_cutoff_A:
                    continue
                # [1] a pair two or three bonds apart is a bend in disguise
                separation = hops.get((i, j))
                if separation is not None and separation < config.min_pair_hops:
                    continue
            detail = donors.get((i, j)) if bonded else None
            priority = 0 if detail is not None else 1
            for family in families:
                spec = FAMILIES[family]
                if spec["kind"] != "bond":
                    continue
                if spec["bonded"] is not None and spec["bonded"] != bonded:
                    continue
                if spec["sign"] is not None:
                    yield family, (i, j), spec["sign"], priority, 0, None, detail
                else:
                    yield family, (i, j), -1, priority, 0, None, detail
                    yield family, (i, j), 1, priority, 0, None, detail


def _angle_directions(graph, orbits, families):
    """One representative per symmetry class of angle, by certified orbits (P1-7)."""
    adjacency = _adjacency(graph.edges, graph.index)
    seen = set()
    for vertex, neighbours in adjacency.items():
        for a in range(len(neighbours)):
            for b in range(a + 1, len(neighbours)):
                i, k = neighbours[a], neighbours[b]
                key = (orbits[vertex], tuple(sorted((orbits[i], orbits[k]))))
                if key in seen:
                    continue
                seen.add(key)
                for family in families:
                    if FAMILIES[family]["kind"] != "angle":
                        continue
                    yield family, (i, vertex, k), -1, 1, 0, None, None
                    yield family, (i, vertex, k), 1, 1, 0, None, None


@dataclass(frozen=True)
class Torsion:
    """A rotatable-bond torsion together with the symmetry that makes part of it redundant.

    symmetry_order n means the conformational space of this coordinate is S^1 / C_n, so
    only the fundamental domain [0, 2*pi/n) needs searching. A methyl rotor has n = 3;
    ethane, symmetric at both ends, has n = lcm(3, 3) = 3.
    """

    indices: tuple[int, int, int, int]
    symmetry_order: int
    moving: tuple[int, ...] = ()  # rigid fragment across the bond, empty inside a ring

    @property
    def is_rotor(self):
        return self.symmetry_order > 1

    @property
    def domain(self):
        return 2 * np.pi / self.symmetry_order

    def to_dict(self):
        return {
            "indices": list(self.indices),
            "symmetry_order": self.symmetry_order,
            "fundamental_domain_rad": self.domain,
            "realization": "rigid_rotation" if self.moving else "continuation",
        }


def _end_symmetry(side, labels):
    """Rotation order contributed by one end: equivalent substituents make it symmetric."""
    if len(side) > 1 and len({labels[n] for n in side}) == 1:
        return len(side)
    return 1


def orbits_from_group(group, index):
    """Certified orbits of an (explicitly enumerated) group, ids by index scan."""
    order = [int(a) for a in index]
    slot = {atom: k for k, atom in enumerate(order)}
    orbits, next_id = {}, 0
    for atom in order:
        if atom in orbits:
            continue
        for member in sorted({int(g[slot[atom]]) for g in group}):
            orbits[member] = next_id
        next_id += 1
    return orbits


def geometric_automorphisms(atoms, group, active=None, tolerance=0.05):
    """The subgroup of `group` that also (approximately) preserves the current geometry.

    P1-7's whole point: two bridge bonds can be equivalent as a GRAPH while the fixed
    conformation breaks that equivalence, and deduplicating proposals by graph orbits
    merges them -- coverage lost without `coverage_censored` ever seeing it. The
    geometric subgroup is computed by one Kabsch fit per group element: the permutation
    preserves the geometry iff the atoms and their images superpose within `tolerance`
    RMSD. `tolerance` is a DEDUPLICATION threshold, not an identity threshold -- it
    decides how eager the proposal basis is to call two directions the same, never
    whether two structures are -- and it enters RATIONALE with this exact reading
    (PLAN section 1.3).
    """
    from .state import resolve_active

    index = resolve_active(atoms, active)
    positions = atoms.get_positions()
    order = [int(a) for a in index]
    slot = {atom: k for k, atom in enumerate(order)}
    keep = []
    for permutation in group:
        source = np.asarray([positions[a] for a in order])
        target = np.asarray([positions[int(permutation[slot[a]])] for a in order])
        source_centred = source - source.mean(axis=0)
        target_centred = target - target.mean(axis=0)
        covariance = source_centred.T @ target_centred
        u, _, vt = np.linalg.svd(covariance)
        # no reflection: a permutation that preserves a chiral geometry as a mirror
        # image is not a motion, and admitting it would merge enantiomer directions
        if np.linalg.det(u @ vt) < 0:
            u[:, -1] *= -1
        rotation = u @ vt
        fitted = source_centred @ rotation
        rmsd = float(np.sqrt(((fitted - target_centred) ** 2).sum(axis=1).mean()))
        if rmsd <= tolerance:
            keep.append(permutation)
    return keep


def _fragment_across(edges, index, j, k):
    """Atoms reachable from k without using the j-k bond, empty when the bond is in a ring."""
    adjacency = _adjacency(edges, index)
    seen, stack = {k}, [k]
    while stack:
        node = stack.pop()
        for neighbour in adjacency[node]:
            if (node, neighbour) in ((j, k), (k, j)):
                continue
            if neighbour not in seen:
                seen.add(neighbour)
                stack.append(neighbour)
    if j in seen:
        return ()  # not a bridge: no rigid fragment exists
    return tuple(sorted(seen - {k}))


def _bond_sides(graph, i, j):
    """The two rigid sides of a bond: what each end reaches without crossing it.

    None when the ends stay connected without the bond -- a ring member has no rigid
    split, exactly as a ring torsion has no fragment to rotate. Written out rather than
    reusing `_fragment_across`, which returns an empty tuple both for a ring bond and for
    a terminal atom: a terminal atom is a perfectly good rigid side (itself), so the two
    cases must not share a return value here.

    Atoms in neither side are spectator fragments and are left alone. The two ends being
    in different fragments to begin with is the normal case for an approach probe: each
    molecule is then its own side and moves rigidly.
    """
    adjacency = _adjacency(graph.edges, graph.index)

    def reach(start):
        seen, stack = {start}, [start]
        while stack:
            node = stack.pop()
            for neighbour in adjacency[node]:
                if (node, neighbour) in ((i, j), (j, i)):
                    continue
                if neighbour not in seen:
                    seen.add(neighbour)
                    stack.append(neighbour)
        return seen

    side_i = reach(i)
    if j in side_i:
        return None
    return tuple(sorted(side_i)), tuple(sorted(reach(j)))


def _angle_sides(graph, i, j, k):
    """The two rigid sides of an angle: what each arm reaches without passing the vertex.

    None when the arms stay connected without the vertex -- inside a ring there is no
    rigid split, exactly as for a ring bond or a ring torsion. Atoms in neither side are
    spectators and are left alone.
    """
    adjacency = _adjacency(graph.edges, graph.index)

    def reach(start):
        seen, stack = {start}, [start]
        while stack:
            node = stack.pop()
            for neighbour in adjacency[node]:
                if neighbour == j or neighbour in seen:
                    continue
                seen.add(neighbour)
                stack.append(neighbour)
        return seen

    side_i = reach(i)
    if k in side_i:
        return None
    side_k = reach(k)
    if side_i & side_k:
        return None
    return tuple(sorted(side_i)), tuple(sorted(side_k))


def rotatable_torsions(graph, labels, orbits=None):
    """One representative dihedral per rotatable bond.

    The dedup key is ORBITS when supplied and colours otherwise. On the recorded
    corpus the two agree on every frame (orbit_vs_colour_probe.py, 2026-09-20), so
    orbit keys change nothing there; they certify what colours only guessed (P1-7).
    The proposal side passes GEOMETRIC orbits -- a fixed conformation can break a
    graph symmetry, and merging geometrically distinct bridge bonds was invisible
    coverage loss. The measurement side passes graph orbits, keeping the coordinate
    set a property of the substance rather than of one conformer.

    A rotatable bond is a bridge -- so not a ring member, whose torsion the ring already
    fixes -- with substituents on both ends. All dihedrals about one bond differ by a
    constant, so one representative carries the whole degree of freedom; the substituents
    are chosen by highest canonical colour to make the choice deterministic. These same
    coordinates describe a molecule's conformational freedom, so the conformer reservoir
    measures diversity in them.

    Notes
    -----
    [1] A rotation that permutes equivalent substituents is a symmetry operation, so it can
        only reach conformers that are symmetry copies of the starting one. Such a bond
        consumed fifteen of ethanol's twenty-four trials and could not have found anything.
        It is kept, because a rotor may still couple to a reaction coordinate, but only its
        fundamental domain is searched and it is ranked last.
    """
    adjacency = _adjacency(graph.edges, graph.index)
    seen, genuine, rotors = set(), [], []
    for j, k in sorted(_bridges(graph.edges, graph.index)):
        left = [n for n in adjacency[j] if n != k]
        right = [n for n in adjacency[k] if n != j]
        if not left or not right:
            continue
        dedup = orbits if orbits is not None else labels
        key = tuple(sorted((dedup[j], dedup[k])))
        if key in seen:
            continue
        seen.add(key)
        # [1] a rotation that permutes equivalents is a symmetry. Deliberately the
        # GRAPH symmetry (colours), not the geometric subgroup: the rotor's period is
        # a property of the moving fragment, and measuring it against a global rigid
        # motion of the whole conformation breaks a methyl's C3 the way ethanol's
        # CH2/O side chain breaks it for every methyl everywhere -- tripling every
        # rotor's amplitude grid with trials that are symmetry copies by the
        # fragment's own internal symmetry. The geometric orbits above are for WHICH
        # coordinates exist; the period below is for how far one coordinate searches.
        order = _lcm(_end_symmetry(left, labels), _end_symmetry(right, labels))
        torsion = Torsion(
            indices=(
                max(left, key=lambda n: labels[n]),
                j,
                k,
                max(right, key=lambda n: labels[n]),
            ),
            symmetry_order=order,
            moving=_fragment_across(graph.edges, graph.index, j, k),
        )
        (rotors if torsion.is_rotor else genuine).append(torsion)
    return genuine, rotors


def _lcm(a, b):
    from math import gcd

    return a * b // gcd(a, b)


def _dihedral_directions(graph, labels, orbits, families):
    genuine, rotors = rotatable_torsions(graph, labels, orbits)
    for rank, group in ((0, genuine), (1, rotors)):
        for torsion in group:
            for family in families:
                if FAMILIES[family]["kind"] != "dihedral":
                    continue
                yield family, torsion.indices, -1, 1, rank, torsion, None
                yield family, torsion.indices, 1, 1, rank, torsion, None


def candidate_id(family, indices, sign):
    """A candidate direction's identity, stable across everything except the structure.

    It has to survive a change in what gets SELECTED -- that is the whole point. Without
    it the ledger the 2026-09-09 ruling requires (candidate -> selected -> executed ->
    refused) cannot be written, and every subdivision of `unanswered` is unattributable:
    "never selected" and "searched and missed" would look the same in the record.

    Derived from the coordinate the direction acts on, not from enumeration order, so
    widening the budget adds ids without renaming any. Atom indices are the structure's
    own; a relabelled structure is a different record and is not claimed to match.
    """
    return f"{family}:{'-'.join(str(int(i)) for i in indices)}:{'+' if sign >= 0 else '-'}"


def spectral_overlap(atoms, kind, indices, modes):
    """How much of a direction lies in the softest non-torsional modes, in [0, 1].

    `displace` steps along g/m in Cartesian coordinates -- the internal coordinate's own
    gradient over the masses -- so in mass-weighted coordinates q = sqrt(m) x that is
    g/sqrt(m). Using g, or g*sqrt(m), computes the overlap with a direction nothing
    travels in; that is the same mass-weighting slip the relay made once already
    (M^(1/2) g where M^-1 g was meant), so the conversion is written out rather than
    inherited. `modes` arrives mass-weighted and unit norm from `soft_internal_modes`.
    """
    gradient = internal.gradient(atoms.positions, kind, indices)
    vector = (gradient / np.sqrt(atoms.get_masses())[:, None]).reshape(-1)
    norm = float(np.linalg.norm(vector))
    if norm < 1e-30:
        return 0.0
    vector = vector / norm
    return max(abs(float(vector @ mode)) for _, mode in modes)


def propose(atoms, config, seed, report=None, soft_modes=()):
    """Amplitude scans, one list per direction, in a seeded deterministic order.

    `report`, if given, is filled with the candidate count, quota and selected count per
    family, plus whether the global cap bound the selection. A fixed total budget cannot
    hold a fixed hit rate while candidates grow with the molecule, so that limitation is
    written down where the run can be audited rather than left to be inferred.

    Notes
    -----
    [1] Stable sort by rank keeps the shuffle inside each band, so deprioritizing symmetry
        rotors survives randomisation instead of being undone by it. Sort by the
        within-family rank only. The priority band is element 3 and is handled by the
        two-pass allocation below, because it means something different: rank orders
        candidates inside a family, band decides which candidates get a front execution slot
        at all.
    [2] Allocate the direction budget per family, round robin. A pooled shuffle would
        allocate it by enumeration count instead: ethanol offers 35 atom pairs but only two
        rotatable bonds, so pair kicks would crowd out the torsion that is the one way to
        reach its other basins. Spanning the basis has to hold for the directions actually
        sampled, not only for the directions that could have been.
    [4] Without a spectrum every candidate's rank is 0, so the stable sort in [1] is a
        no-op and the order inside a family is the seeded shuffle alone. That is a lottery,
        and on b04 it happened to draw the C-Cl stretch first -- the one direction of 28
        that reacted -- with no mechanism putting it there. Given the softest non-torsional
        modes the rank becomes the overlap with them, which put that same direction first
        by construction at 0.9651 against 0.6893 for the next distinct one
        (`docs/experiments/spectral_proposal_probe.py`). The band in [1] is untouched: a
        hydrogen-bond donor still gets a front slot on its own evidence, and the spectrum
        only orders candidates inside a band. Absent modes changes nothing, which is what
        a microstate whose Hessian could not be had must fall back to -- a ranking is never
        worth failing a search for.
    [3] Search S^1 / C_n only. What disqualifies an amplitude is not being larger than the
        domain but being close to any symmetry copy of zero: 2.09 rad on a methyl rotor is
        119.75 degrees, a quarter of a degree from the identity, and a trial spent there
        cannot find anything. So the test is the distance to the nearest multiple of the
        domain.
    """
    graph = encode(atoms, config.bond_scale, active=config.active_atoms)
    labels = canonical_labels(atoms.numbers, graph.edges, graph.index)
    # The proposal basis deduplicates by GEOMETRIC orbits: a fixed conformation can
    # break a graph symmetry, and two bridge bonds a colour or graph orbit merges may
    # be geometrically distinct -- merging them was invisible coverage loss (P1-7).
    #
    # This is a LARGE widening, and the size of it is measured rather than assumed.
    # An earlier version of this comment cited the 2026-09-20 orbit probe's "zero
    # diff" for the opposite claim; that probe compared colours with GRAPH orbits and
    # says nothing about this pair. Comparison 4 was added to it on 2026-09-21 and
    # measures the pair that actually ships: over the same 2136-frame corpus geometry
    # breaks a graph symmetry on 2134 frames, and the angle+dihedral candidate count
    # goes 229,472 -> 331,840 (+44.6%). A relaxed conformer almost never keeps its
    # graph symmetry to within the 0.05 A deduplication tolerance, so "only where
    # geometry breaks symmetry" is, in practice, almost everywhere.
    #
    # The cost is real and lands on the budget, because the per-family quota is a
    # FRACTION of the candidate count: more candidates means more selected directions,
    # more trials, and coverage_censored true more often. That is the honest direction
    # -- over-count rather than merge on a guess -- but it is not free, and a run's
    # ledger is where it has to be read. b04's frozen grid is unaffected (|Aut| 6 -> 3,
    # bend 4 either way, candidates 28 / selected 14 exactly as frozen), so the frozen
    # record still recomputes; that is a property of a six-atom system, not a general
    # one.
    from .canonical import canonical_form

    form = canonical_form(atoms.numbers, graph.edges, graph.index)
    if form.info.get("enumerated"):
        subgroup = geometric_automorphisms(atoms, form.group, config.active_atoms)
        basis_orbits = orbits_from_group(subgroup, graph.index)
    else:
        # uncertified group: over-count rather than merge on a guess -- the identity
        # subgroup leaves every atom its own orbit, which is the same failure
        # direction the event key takes
        basis_orbits = {int(a): k for k, a in enumerate(graph.index)}
    families = tuple(config.families)
    directions = list(_bond_directions(atoms, graph, config, families))
    directions += list(_angle_directions(graph, basis_orbits, families))
    directions += list(_dihedral_directions(graph, labels, basis_orbits, families))
    # [4] rank by the spectrum when there is one
    if soft_modes:
        directions = [
            entry[:4]
            + (-spectral_overlap(atoms, FAMILIES[entry[0]]["kind"], entry[1], soft_modes),)
            + entry[5:]
            for entry in directions
        ]
    grouped = {}
    for entry in directions:
        grouped.setdefault(entry[0], []).append(entry)
    rng = np.random.default_rng(seed)
    for family in grouped:
        entries = grouped[family]
        order = rng.permutation(len(entries))
        shuffled = [entries[i] for i in order]
        # [1] stable sort keeps the shuffle inside each band
        grouped[family] = sorted(shuffled, key=lambda entry: entry[4])

    # [2] allocate the direction budget per family, round robin
    available = [family for family in families if grouped.get(family)]
    counts = {family: len(grouped[family]) for family in available}
    # The quota grows with the number of candidates: a constant budget cannot keep a
    # constant chance of sampling any given direction while the candidates grow with the
    # molecule. Ethanol offers eight bonds and 3-oxobutanal eleven, and the reaction
    # coordinate lost that draw.
    exhaustive = getattr(config, "direction_mode", "sampled") == "exhaustive"
    quotas = {
        family: counts[family]
        if exhaustive
        else max(
            config.min_directions_per_family,
            math.ceil(config.direction_quota_fraction * counts[family]),
        )
        for family in available
    }
    # In exhaustive mode the cap is not a cap: it may not silently decide coverage, and
    # the caller is required to have refused to start if the trial budget cannot pay for
    # the whole enumeration (search.ReactionSearch.run does that check).
    cap = sum(counts.values()) if exhaustive else config.max_directions
    selected, taken = [], {family: 0 for family in available}

    # Two passes. The priority band first, so a direction the reactant's own geometry marks
    # as unusual gets a front execution slot -- the trial budget is spent in the order this
    # list is returned, so being selected late is nearly the same as not being selected.
    for band in (0, 1):
        while len(selected) < cap:
            progressed = False
            for family in available:
                if len(selected) >= cap:
                    break
                if taken[family] >= quotas[family]:
                    continue
                entry = next((e for e in grouped[family] if e[3] == band), None)
                if entry is None:
                    continue
                grouped[family].remove(entry)
                selected.append(entry)
                taken[family] += 1
                progressed = True
            if not progressed:
                break

    if report is not None:
        candidate_ids = [candidate_id(e[0], e[1], e[2]) for e in directions]
        selected_ids = [candidate_id(e[0], e[1], e[2]) for e in selected]
        cap_bound = (not exhaustive) and len(selected) >= config.max_directions
        # Coverage is censored whenever anything was left unsampled, and WHICH mechanism
        # did it matters: the cap is a declared cost limit, the quota fraction censors
        # silently and at every system size -- b04's probe3 selected 14 of 28 candidates
        # with cap_bound false on six atoms. Either way a negative conclusion drawn from
        # this proposal is inadmissible (the 2026-09-09 ruling, point 3); the status is
        # recorded here so that judgement is made from the record rather than remembered.
        uncovered = [i for i in candidate_ids if i not in set(selected_ids)]
        report.update(
            families={
                family: {
                    "candidates": counts[family],
                    "quota": quotas[family],
                    "selected": taken[family],
                    "priority_candidates": sum(
                        1 for e in directions if e[0] == family and e[3] == 0
                    ),
                }
                for family in available
            },
            mode=("exhaustive" if exhaustive else "sampled"),
            max_directions=config.max_directions,
            effective_cap=cap,
            selected_total=len(selected),
            cap_bound=cap_bound,
            priority_selected=sum(1 for e in selected if e[3] == 0),
            candidate_ids=candidate_ids,
            selected_ids=selected_ids,
            uncovered_ids=uncovered,
            coverage_censored=bool(uncovered),
            censored_by=(
                None
                if not uncovered
                else "max_directions"
                if cap_bound
                else "direction_quota_fraction"
            ),
            meaning=(
                "candidate counts, per-family quotas and what was actually taken; "
                "cap_bound true means the global cap, not the quotas, decided; "
                "coverage_censored true means candidates were left unsampled, so a "
                "negative result from this proposal is not evidence of absence"
            ),
        )

    result = []
    for family, indices, sign, band, _rank, torsion, detail in selected:
        amplitudes = getattr(config, FAMILIES[family]["amplitudes"])
        if torsion is not None and torsion.is_rotor and FAMILIES[family]["mode"] == "displace":
            # [3] search S^1 / C_n only
            amplitudes = tuple(
                a
                for a in amplitudes
                if abs((a + torsion.domain / 2) % torsion.domain - torsion.domain / 2)
                > config.torsion_symmetry_margin_rad
            )
            if not amplitudes:
                continue
        direction_seed = int(rng.integers(0, 2**32))
        result.append(
            [
                Probe(
                    family,
                    indices,
                    sign,
                    float(a),
                    direction_seed,
                    priority=band,
                    priority_reason=detail,
                )
                for a in amplitudes
            ]
        )
    return result


def fragment_translations(atoms, graph, reference=None):
    """Three relative-translation directions per fragment pair (B'', PLAN item 5.2).

    The basis is decision 5 (2026-09-20): the COM line between the two heaviest
    fragments plus two deterministic orthonormal completions of it -- the coupled
    residual on b05's elimination endpoints is inter-fragment DISPLACEMENT, and this
    basis points at it directly. The completions are built by Gram-Schmidt against
    the coordinate axes in fixed order, so the basis is a function of the geometry
    and nothing else.

    INERT BY DESIGN: enumerate only. These become searchable/polishable coordinates
    when the three `fragment_translation_*` thresholds exist, and those must come
    from a measurement on the b05 endpoints (docs/experiments/
    b05_fragment_translation_measurement.py, needs the MACE checkpoint) -- not from
    converting the radian set, per `polish_soft_modes` footnote [1]. Until then no
    caller may wire them in; the function exists so the measurement has exactly the
    enumeration it will justify.
    """
    from .chemistry import _fragments

    fragments = [sorted(int(a) for a in part) for part in _fragments(graph.edges, graph.index)]
    if len(fragments) < 2:
        return ()
    fragments.sort()
    masses = atoms.get_masses()
    by_mass = sorted((m for m in fragments if len(m) >= 1), key=lambda m: -masses[m].sum())
    left, right = by_mass[0], by_mass[1]
    line = np.average(atoms.positions[right], axis=0, weights=masses[right]) - np.average(
        atoms.positions[left], axis=0, weights=masses[left]
    )
    norm = float(np.linalg.norm(line))
    if norm < 1e-12:
        return ()
    line = line / norm
    basis = [line]
    for axis in np.eye(3):
        candidate = axis - line * float(axis @ line)
        if np.linalg.norm(candidate) < 1e-9:
            continue
        basis.append(candidate / np.linalg.norm(candidate))
        if len(basis) == 3:
            break
    return tuple(
        {"members": list(left), "against": list(right), "direction": [float(x) for x in d]}
        for d in basis
    )
