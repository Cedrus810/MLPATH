"""Chemical identity of a structure, separated from its conformation.

The reaction network is a graph over chemical macrostates; conformers of one
macrostate live inside it. That split only works if chemical identity is computed in a
way that is blind to conformation and blind to atom labelling:

Blind to conformation. A signed torsion about a rotatable single bond distinguishes
gauche+ from gauche-, which are the same substance. Including it would split one
macrostate in two and reproduce the combinatorial explosion the hierarchy exists to
prevent. Only configurational stereochemistry -- parity that cannot change without
breaking a bond -- belongs in the key: tetrahedral parity at genuine stereocentres and
E/Z parity across bonds that are both short enough to be multiple and are bridges
rather than ring members.

Blind to atom labelling. Under the fixed atom mapping of this MVP, transferring any
one of three equivalent methyl hydrogens yields a different set of atom indices for
the same reaction. Keying on raw edge sets would record three products where there is
one, so the graph part of the key is a colour-refinement hash that is invariant under
graph automorphism. Refinement is invariant by construction but not a provably
complete canonical form; for molecular graphs at these sizes it separates every case
encountered, and a collision would merge two chemistries rather than invent one.

Charge and spin enter the key only when the potential actually responds to them.
MACE-OFF returns identical energies for total charge 0, -1 and +1, so declaring
charge_sensitive keeps charged chemistry out of the search space instead of letting it
in unnoticed.
"""
from collections import Counter
import hashlib
import json
import numpy as np
from ase.data import covalent_radii, chemical_symbols
from .state import encode, resolve_active


def _digest(text):
    return hashlib.sha256(text.encode()).hexdigest()[:16]


def canonical_labels(numbers, edges, index):
    """Automorphism-invariant colour per atom; equal colours mean equivalent atoms."""
    adjacency = {int(i): set() for i in index}
    for a, b in edges:
        adjacency[a].add(b)
        adjacency[b].add(a)
    labels = {int(i): chemical_symbols[int(numbers[i])] for i in index}
    for _ in range(len(index)):
        refined = {i: _digest(labels[i] + "|" + ",".join(sorted(labels[j]
                                                               for j in adjacency[i])))
                   for i in labels}
        if Counter(refined.values()) == Counter(labels.values()):
            break
        labels = refined
    return labels


def _bridges(edges, index):
    """Bonds whose removal disconnects the graph, so ring bonds are excluded."""
    adjacency = {int(i): [] for i in index}
    for a, b in edges:
        adjacency[a].append(b)
        adjacency[b].append(a)
    discovery, low, found, clock = {}, {}, set(), [0]

    def walk(node, parent):
        discovery[node] = low[node] = clock[0]
        clock[0] += 1
        for neighbour in adjacency[node]:
            if neighbour == parent:
                continue
            if neighbour not in discovery:
                walk(neighbour, node)
                low[node] = min(low[node], low[neighbour])
                if low[neighbour] > discovery[node]:
                    found.add((min(node, neighbour), max(node, neighbour)))
            else:
                low[node] = min(low[node], discovery[neighbour])

    for atom in adjacency:
        if atom not in discovery:
            walk(atom, None)
    return found


def _fragments(edges, index):
    parent = {int(i): int(i) for i in index}

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for a, b in edges:
        parent[root(a)] = root(b)
    groups = {}
    for atom in parent:
        groups.setdefault(root(atom), []).append(atom)
    return list(groups.values())


def _formula(numbers, atoms_in_fragment):
    from ase import Atoms
    return Atoms(numbers=[int(numbers[i]) for i in sorted(atoms_in_fragment)]
                 ).get_chemical_formula()


def _tetrahedral_parities(positions, numbers, edges, index, labels):
    """Parity at atoms with four neighbours that are all inequivalent.

    Four distinct colours is exactly the condition for a real stereocentre: if two
    neighbours share a colour, swapping them is an automorphism and there is no
    configuration to record. That is why a methyl carbon contributes nothing.
    """
    adjacency = {int(i): [] for i in index}
    for a, b in edges:
        adjacency[a].append(b)
        adjacency[b].append(a)
    parities = []
    for atom, neighbours in adjacency.items():
        if len(neighbours) != 4:
            continue
        colours = [labels[n] for n in neighbours]
        if len(set(colours)) != 4:
            continue
        ordered = [n for _, n in sorted(zip(colours, neighbours))]
        base = positions[ordered[0]]
        volume = np.linalg.det(np.array([positions[ordered[1]] - base,
                                         positions[ordered[2]] - base,
                                         positions[ordered[3]] - base]))
        parities.append([labels[atom], int(np.sign(volume))])
    return sorted(parities)


# Neutral closed-shell valences. Deliberately small: the first version has to be right on
# C/H/O, which covers every benchmark, plus the monovalent halogens, which are terminal and
# unambiguous in neutral organic molecules. Nitrogen, phosphorus and sulfur are left out on
# purpose -- their valences are not single-valued in that setting -- and anything outside
# this table is reported unresolved rather than guessed at.
VALENCE = {1: 1, 6: 4, 8: 2, 9: 1, 17: 1, 35: 1, 53: 1}


def admissible_bond_orders(numbers, edges, index, limit=20000):
    """Every integer bond-order assignment consistent with neutral closed-shell valences.

    A set, not one Lewis structure. Which edge of benzene is "the double bond" depends on
    which Kekule structure was picked, so any property read off a single assignment would
    be an artefact of that choice; asking instead what holds across all of them is well
    defined. This is also the only formulation that is conformation-independent, which
    chemical identity has to be -- measured on the malonaldehyde enol, the same C=C ranges
    over 0.66 to 4.59 eV/rad^2 depending on whether the intramolecular hydrogen bond is
    closed, so neither a bond length nor a stiffness can carry this.

    Returns (solutions, info). solutions is None when the valences cannot be satisfied,
    an element falls outside the supported set, or the enumeration exceeds the limit. The
    caller must then record the stereochemistry as unresolved rather than guess at it.
    """
    active = [int(i) for i in index]
    if any(int(numbers[i]) not in VALENCE for i in active):
        outside = sorted({chemical_symbols[int(numbers[i])] for i in active
                          if int(numbers[i]) not in VALENCE})
        return None, {"reason": "unsupported_elements", "elements": outside}
    incident = {i: [] for i in active}
    ordered = sorted(tuple(sorted(e)) for e in edges)
    for edge in ordered:
        for atom in edge:
            incident[atom].append(edge)
    remaining = {i: VALENCE[int(numbers[i])] - len(incident[i]) for i in active}
    if any(value < 0 for value in remaining.values()):
        return None, {"reason": "degree_exceeds_valence"}

    # Every edge starts at order one; what is left to place is the surplus valence.
    solutions, visited = [], 0

    def search(position, surplus, assignment):
        nonlocal visited
        visited += 1
        if visited > limit:
            raise RuntimeError("enumeration limit")
        if position == len(ordered):
            if all(value == 0 for value in surplus.values()):
                solutions.append(dict(assignment))
            return
        a, b = ordered[position]
        headroom = min(surplus[a], surplus[b], 2)
        for extra in range(headroom, -1, -1):
            surplus[a] -= extra
            surplus[b] -= extra
            assignment[(a, b)] = 1 + extra
            search(position + 1, surplus, assignment)
            surplus[a] += extra
            surplus[b] += extra
        assignment.pop((a, b), None)

    try:
        search(0, dict(remaining), {})
    except RuntimeError:
        return None, {"reason": "enumeration_limit", "limit": limit}
    if not solutions:
        return None, {"reason": "no_assignment_satisfies_valences"}
    return solutions, {"reason": None, "count": len(solutions)}


def locked_edges(numbers, edges, index, limit=20000):
    """Edges that are at least double in EVERY admissible assignment.

    Only those can carry E/Z: a bond that is double in one Kekule structure and single in
    another has no configuration to record, and a bond that is single in all of them cannot
    hold one. Returns (edges, info), with edges None when the assignment set is unresolved.
    """
    solutions, info = admissible_bond_orders(numbers, edges, index, limit)
    if solutions is None:
        return None, info
    always = frozenset(edge for edge in solutions[0]
                       if all(solution[edge] >= 2 for solution in solutions))
    return always, {**info, "locked": sorted(list(e) for e in always)}


def _locked_bond_parities(positions, numbers, edges, index, labels, locked):
    """E/Z parity across bonds that are double in every admissible bond-order assignment.

    Which bonds those are is decided by `locked_edges`, from the graph alone, because
    chemical identity has to be conformation-independent and neither of the two things
    tried before is. A bond length put a hard classifier inside an identity function and
    the malonaldehyde enol's C-OH sat 0.0002 A from it, so two structures differing by
    0.0003 A became two chemical nodes. Torsional stiffness moved that line rather than
    removing it: the same C=C measures 4.59 eV/rad^2 with the intramolecular hydrogen bond
    closed and 0.66 with it open, and those are conformers of one substance.

    An empty `locked` set therefore means "no bond can carry E/Z", and None is not passed
    here at all -- an unresolved valence is recorded in the key by the caller instead.
    """
    adjacency = {int(i): [] for i in index}
    for a, b in edges:
        adjacency[a].append(b)
        adjacency[b].append(a)
    parities = []
    for a, b in sorted(edges):
        if (a, b) not in locked:
            continue
        picked = []
        for near, far in ((a, b), (b, a)):
            options = [n for n in adjacency[near] if n != far]
            if not options:
                break
            best = max(labels[n] for n in options)
            if sum(1 for n in options if labels[n] == best) != 1:
                break            # substituents are equivalent: no E/Z to record
            picked.append(next(n for n in options if labels[n] == best))
        if len(picked) != 2:
            continue
        axis = positions[b] - positions[a]
        axis = axis / np.linalg.norm(axis)
        u = positions[picked[0]] - positions[a]
        v = positions[picked[1]] - positions[b]
        u = u - np.dot(u, axis) * axis
        v = v - np.dot(v, axis) * axis
        if min(np.linalg.norm(u), np.linalg.norm(v)) < 1e-6:
            continue
        parities.append([tuple(sorted((labels[a], labels[b]))), int(np.sign(np.dot(u, v)))])
    return sorted([[list(pair), parity] for pair, parity in parities])


def chemical_key(atoms, bond_scale=1.2, active=None, bond_order_limit=20000,
                 charge_sensitive=False, charge=None, multiplicity=None):
    """Identity of the chemical macrostate this structure belongs to.

    Returns the components as well as the digest so a record explains itself; two
    structures belong to the same macrostate exactly when their "key" agrees.
    """
    if not charge_sensitive and charge not in (None, 0):
        raise ValueError("charge_sensitive is false, so a nonzero total charge cannot be "
                         "represented; the potential does not respond to it")
    index = resolve_active(atoms, active)
    graph = encode(atoms, bond_scale, active=active)
    numbers = atoms.numbers
    labels = canonical_labels(numbers, graph.edges, index)
    locked, orders = locked_edges(numbers, graph.edges, index, bond_order_limit)
    components = {
        "graph_hash": _digest("|".join(sorted(labels.values())) + f"#{len(graph.edges)}"),
        "fragments": sorted(_formula(numbers, group)
                            for group in _fragments(graph.edges, index)),
        "tetrahedral_parity": _tetrahedral_parities(atoms.positions, numbers,
                                                    graph.edges, index, labels),
        "locked_bond_parity": _locked_bond_parities(atoms.positions, numbers, graph.edges,
                                                   index, labels, locked or frozenset()),
        # An unresolved valence is stated, never guessed at. It is part of the key, so an
        # unresolved substance is never compared against a resolved one on the quiet.
        "stereo_unresolved": None if locked is not None else orders["reason"],
        "charge_sensitive": bool(charge_sensitive),
    }
    if charge_sensitive:
        components["charge"] = charge
        components["multiplicity"] = multiplicity
    components["key"] = hashlib.sha256(
        json.dumps(components, sort_keys=True).encode()).hexdigest()
    components["note"] = ("conformation deliberately excluded; only configurational "
                          "stereochemistry contributes")
    return components


def same_chemistry(a, b, bond_scale=1.2, active=None, **kwargs):
    return (chemical_key(a, bond_scale, active, **kwargs)["key"]
            == chemical_key(b, bond_scale, active, **kwargs)["key"])


def _factorial(n):
    result = 1
    for value in range(2, n + 1):
        result *= value
    return result


def automorphisms(atoms, bond_scale=1.2, active=None, limit=20000):
    """Edge-preserving permutations of equivalent atoms, as index arrays.

    Conformer identity is measured on geometry with fixed atom indices, so rotating a
    methyl group by 120 degrees moves every hydrogen to an equivalent site and looks
    like a new conformer even though nothing physical changed. Minimising RMSD over
    this group removes those duplicates, which is the same automorphism argument the
    chemical key already uses, applied one level down.

    Candidates are permutations inside colour classes; only those that map the edge set
    onto itself are kept. The identity is always first. If the candidate count exceeds
    limit the group is not enumerated and only the identity is returned, so the caller
    over-counts conformers rather than silently mis-merging distinct ones.
    """
    from itertools import permutations, product
    index = resolve_active(atoms, active)
    graph = encode(atoms, bond_scale, active=active)
    labels = canonical_labels(atoms.numbers, graph.edges, index)
    classes = {}
    for atom in index:
        classes.setdefault(labels[int(atom)], []).append(int(atom))
    orbits = [group for group in classes.values() if len(group) > 1]
    total = 1
    for group in orbits:
        total *= _factorial(len(group))
        if total > limit:
            return [np.asarray(index)], {"enumerated": False, "reason": "candidate_limit",
                                         "orbit_sizes": [len(g) for g in orbits]}
    position = {int(atom): i for i, atom in enumerate(index)}
    found = []
    for choice in product(*[permutations(group) for group in orbits]):
        mapping = {int(atom): int(atom) for atom in index}
        for group, permuted in zip(orbits, choice):
            for source, target in zip(group, permuted):
                mapping[source] = target
        if {(min(mapping[a], mapping[b]), max(mapping[a], mapping[b]))
                for a, b in graph.edges} != graph.edges:
            continue
        found.append(np.asarray([mapping[int(atom)] for atom in index]))
    found.sort(key=lambda p: not np.array_equal(p, np.asarray(index)))
    return found, {"enumerated": True, "count": len(found),
                   "orbit_sizes": [len(g) for g in orbits]}



def _render_bonds(edges, mapping, numbers):
    """Bonds as element-tagged pairs under one relabelling, orientation-free."""
    rendered = []
    for a, b in edges:
        x, y = mapping[a], mapping[b]
        pair = sorted(((int(numbers[x]), x), (int(numbers[y]), y)))
        rendered.append(tuple(value for item in pair for value in item))
    return tuple(sorted(rendered))


def parity_map(components, field):
    """{locked bond (or centre) -> parity sign} from a node's key components."""
    out = {}
    for entry in components.get(field) or []:
        target, sign = entry
        out[tuple(sorted(target)) if isinstance(target, list) else target] = sign
    return out


def stereo_only_difference(first, second):
    """Do these two key-component sets differ ONLY in configurational stereochemistry?

    The single implementation of this test. It lives here, next to the classifier, because
    the acceptance checker needs the same question answered and having two copies is how the
    checker and the model come to disagree about what a class means.

    Two different things produce a key change with no bond change:

      real    E/Z isomerisation -- the locked-bond SET is identical and a parity sign is
              flipped. No bond changed because none had to; what changed is configurational
              stereochemistry, which is part of chemical identity by design.
      defect  the key changed for some other reason, e.g. conformation leaking into the key
              or a parity computed on a bond that is not actually locked.

    Returns (ok, reason).
    """
    if first is None or second is None:
        return False, "key components not supplied, so the key change cannot be explained"
    if first.get("graph_hash") != second.get("graph_hash"):
        return False, "graph hashes differ, so this is not a stereochemistry-only change"
    if first.get("fragments") != second.get("fragments"):
        return False, "fragments differ"
    flipped = []
    for field in ("locked_bond_parity", "tetrahedral_parity"):
        a, b = parity_map(first, field), parity_map(second, field)
        if set(a) != set(b):
            return False, (f"{field} is defined on different targets: "
                           f"{sorted(a)} vs {sorted(b)}")
        flipped += [target for target in a if a[target] != b[target]]
    if not flipped:
        return False, "no parity differs, so the key change is unexplained"
    # "key" is the digest of everything else, so of course it differs -- that the two keys
    # differ IS the premise of the question being asked, and treating it as counter-evidence
    # is self-refuting. It has to be ignored explicitly, because the two callers do not pass
    # the same shape: the acceptance checker reads `key_components` from network.json, where
    # publish() has already stripped "key", while the model passes `node.components`, where
    # it is still present. Measured on P2 T=300 seed 71: the checker said "parity flipped,
    # everything else identical" and the model said "key also differs" about the same pair.
    # Sharing one implementation is not the same as sharing one input convention.
    ignored = {"key", "locked_bond_parity", "tetrahedral_parity", "note", "automorphisms"}
    for field in set(first) | set(second):
        if field not in ignored and first.get(field) != second.get(field):
            return False, f"{field} also differs"
    return True, f"parity flipped on {flipped}, everything else identical"


def classify_transition(source_key, endpoint_key, broken, formed,
                        source_components=None, endpoint_components=None):
    """What happened between two minima, from the key and the bond delta together.

    Chemical-state identity and reaction-event identity are different questions. The key
    quotients away every redundancy of the state -- including the automorphism that makes
    a degenerate rearrangement end where it began -- and it must, or transferring one of
    three equivalent hydrogens would mint three products. But that same quotient cannot be
    what decides whether a reaction occurred, because it is blind to the path: a proton
    that moves from one oxygen to the other returns a structure with the identical key and
    a bond set that is not the one it started with.

    So both are needed. A vanishing bond delta is a conformational change however far the
    geometry moved; a non-vanishing one is a reaction whether or not the key changed, and
    when the key does not change the reaction is degenerate -- an automerization.

    The fourth cell -- key changed, no bond changed -- was originally named
    `key_change_without_bond_change` on the assumption that it is always a defect, because
    it is impossible if the key is a function of the bond graph alone. But the key is
    deliberately NOT that: it also carries configurational stereochemistry. So E/Z
    isomerisation lands in that cell legitimately, and it is a real event, just a
    higher-energy one. Measured on 3-oxobutanal, twice independently: the same locked bond,
    parity flipped from -1 to +1, graph hash unchanged, 365-404 meV above the reactant.

    Given both key-component sets the fourth cell is therefore split: a parity flip on an
    otherwise identical component set is `configurational_isomerisation`; anything else in
    that cell keeps the old name and remains a defect to be diagnosed. Without the
    components it cannot be told apart, so the conservative name is kept.
    """
    changed = bool(broken or formed)
    if source_key == endpoint_key:
        return "conformational_transition" if not changed else "degenerate_reaction"
    if changed:
        return "reaction"
    stereo, _ = stereo_only_difference(source_components, endpoint_components)
    return "configurational_isomerisation" if stereo else "key_change_without_bond_change"


def reaction_event_key(source, endpoint, bond_scale=1.2, active=None, limit=20000):
    """Identity of a reaction event, quotiented by symmetry and by direction.

    Storing raw atom indices would let equivalent atoms manufacture duplicate channels --
    transferring each of three equivalent hydrogens would look like three reactions -- so
    the broken and formed sets are canonicalized under the automorphism groups of the two
    endpoints, which is the same argument the chemical key already makes one level up.

    Direction is quotiented too: the forward and reverse observations of one channel are
    two observations, not two channels, which matches the existing treatment of a saddle's
    two endpoints as undirected structural evidence rather than as a pair of directed edges.

    If either automorphism group exceeds the enumeration limit only the identity is used,
    so the key over-counts channels rather than silently merging distinct ones -- the same
    failure direction the rest of the symmetry handling takes.
    """
    numbers = source.numbers
    ga = encode(source, bond_scale, active=active)
    gb = encode(endpoint, bond_scale, active=active)
    broken, formed = ga.edges - gb.edges, gb.edges - ga.edges
    index = resolve_active(source, active)

    relabellings, enumerated = [], True
    for structure in (source, endpoint):
        found, info = automorphisms(structure, bond_scale, active, limit)
        enumerated = enumerated and info.get("enumerated", False)
        for permutation in found:
            relabellings.append({int(a): int(b) for a, b in zip(index, permutation)})

    def canonical(first, second):
        return min((_render_bonds(first, mapping, numbers),
                    _render_bonds(second, mapping, numbers))
                   for mapping in relabellings)

    form = min(canonical(broken, formed), canonical(formed, broken))
    # The bond pattern alone is not a global identity: two different substances can break
    # and form the same elements at the same indices, and did in the malonaldehyde and
    # 3-oxobutanal preflights. Tagging the digest with the unordered pair of endpoint graph
    # hashes separates them without giving the key a direction.
    def graph_hash(structure):
        graph = encode(structure, bond_scale, active=active)
        labels = canonical_labels(structure.numbers, graph.edges,
                                  resolve_active(structure, active))
        return _digest("|".join(sorted(labels.values())) + f"#{len(graph.edges)}")

    context = tuple(sorted((graph_hash(source), graph_hash(endpoint))))
    form_with_context = (form, context)
    digest = hashlib.sha256(repr(form_with_context).encode()).hexdigest()[:16]
    return {"key": digest, "canonical": form, "endpoint_graphs": list(context),
            "broken": [list(e) for e in sorted(broken)],
            "formed": [list(e) for e in sorted(formed)],
            "symmetry_enumerated": enumerated,
            "relabellings": len(relabellings),
            "meaning": ("reaction-event identity: canonical under the endpoints' "
                        "automorphisms and undirected. NOT a chemical-state identity")}


def free_aligned_symmetric_rmsd(a, b, permutations_, active=None, free_bonds=(),
                                bond_scale=1.2, samples=24):
    """symmetric_rmsd, additionally minimized over rotations about coordinates that carry
    no energy.

    A torsion whose stiffness has fallen to zero connects the structures on either side of
    it by a path with no barrier, so they are one basin and one microstate however far
    apart their coordinates are. Plain RMSD says otherwise, and a genuinely free rotor
    therefore mints an unbounded number of spurious conformers. Quotienting the comparison
    by those coordinates is the same move the chemical key, the symmetry-aware RMSD and the
    rotor fundamental domain already make, applied to the one redundancy that is energetic
    rather than combinatorial. In Morse-Bott terms it is the neutral subspace: the
    transition-state criterion asks for exactly one unstable direction and says nothing
    about how many flat ones there are, so flat ones must not be allowed to multiply states.

    Which coordinates are free has to be measured, and only on the structure where it was
    measured is the answer known -- so `free_bonds` is supplied by the caller rather than
    guessed here. Passing a stiff or merely soft coordinate would wrongly merge two real
    conformers, which is why nothing infers this from the graph alone.

    Costs no force evaluations: rotation about a bridge bond is exact geometry.
    """
    from .perturbations import step_torsion, torsion_on_bond
    from .state import encode
    best = symmetric_rmsd(a, b, permutations_, active)
    if not free_bonds:
        return best
    graph = encode(b, bond_scale, active=active)
    labels = canonical_labels(b.numbers, graph.edges, graph.index)
    torsions = [t for t in (torsion_on_bond(graph, labels, bond) for bond in free_bonds)
                if t is not None]
    if not torsions:
        return best
    current = b
    grid = np.linspace(0.0, 2 * np.pi, samples, endpoint=False)[1:]
    for _ in range(2):
        improved = False
        for torsion in torsions:
            for delta in grid:
                rotated, _ = step_torsion(current, torsion, float(delta), in_place=False)
                distance = symmetric_rmsd(a, rotated, permutations_, active)
                if distance < best - 1e-12:
                    best, current, improved = distance, rotated, True
        if not improved:
            break
    return best


def symmetric_rmsd(a, b, permutations_, active=None):
    """Smallest Kabsch RMSD over a set of index permutations of b."""
    from .state import aligned_rmsd
    index = resolve_active(a, active)
    best = float("inf")
    for permutation in permutations_:
        relabelled = b.copy()
        positions = relabelled.positions.copy()
        positions[index] = b.positions[permutation]
        relabelled.set_positions(positions)
        best = min(best, aligned_rmsd(a, relabelled, active))
    return best
