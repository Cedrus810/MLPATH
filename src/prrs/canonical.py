"""A canonical form for labelled molecular graphs, by individualization-refinement.

`canonical_labels` (chemistry.py) is a 1-WL colour refinement: invariant under every
automorphism, cheap, and NOT a canonical form -- equal colours are a guess about orbits,
not a certified one, and that single gap is the hole behind three recorded symptoms
(PLAN_2026-09-20.md section 1: B / P1-3 / P1-7). This module replaces the guess with a
construction:

1. refine the element colours to stability (`refine`);
2. a stable colouring whose classes are all singletons is a LEAF -- one candidate
   labelling of the graph;
3. otherwise individualize each member of the lexicographically first non-singleton
   class -- a fixed rule on colour strings, never on enumeration order -- refine, recurse;
4. the canonical form is the leaf with the minimal certificate; its labelling is `order`;
5. two leaves with the same certificate are related by an automorphism (the map carrying
   the one leaf's labelling onto the other's), so one automorphism per equal-certificate
   leaf is the complete group, and orbits computed from THAT group are certified orbits.

No nauty-style pruning, deliberately (PLAN section 1.1): the largest group on record is
b05's |Aut| = 1296, a naive IR tree has as many leaves as the group has elements, and the
backtracking kernel this replaces enumerated that group in 3,889 search nodes / 41 ms --
the same order of magnitude. Pruning without a measured slow case has no object to work
on; the P0 fast path already paid for that lesson once.

Budget: `limit` bounds search nodes -- the vertices of the IR tree, one per
individualization plus the root -- and overflowing it is fail closed: the form is
reported as NOT computed (empty `order` and `certificate`, the identity alone in
`group`, `reason = "candidate_limit"`, the string `chemistry.BUDGET_REASONS` keys on,
so the event and registry paths refuse rather than canonicalise against it). An
incomplete group is never used.

Trust model, stated once: colours are 64-bit truncated SHA-256 digests, and everything
here -- like `canonical_labels` before it -- assumes distinct signature strings do not
collide. A collision would merge colour classes; the equal-certificate arguments above
would then certify a permutation that is not an automorphism. That is the identity
layer's pre-existing assumption, not a new one.
"""

from collections import Counter
from dataclasses import dataclass
import hashlib

import numpy as np
from ase.data import chemical_symbols


def _digest(text):
    # Must stay identical to chemistry._digest. chemistry.py imports this module (the
    # wrappers there are thin by mandate), so importing back is a cycle; two one-line
    # definitions with one job each is the cheaper side of that trade.
    return hashlib.sha256(text.encode()).hexdigest()[:16]


# Appended to a colour when its atom is individualized. Colours are otherwise element
# symbols or hex digests, so this cannot collide with one, and the marker is gone after
# the next refinement pass re-digests every colour anyway. It only has to make the
# individualized atom's colour unique at the moment of individualization.
_INDIVIDUALISED = "\x00individualised"

# The overflow reason. `chemistry.BUDGET_REASONS` keys on this exact string to decide
# that a budget ran out (as opposed to a molecule being unresolvable), and the network
# path turns it into a refusal -- it is a contract between modules, not a message.
_CANDIDATE_LIMIT = "candidate_limit"


def _edge_colour_lookup(edge_colours):
    """Normalize a caller's edge colouring into f(i, j) over ordered atom pairs.

    A callable is called with the ordered pair; a mapping is read with the ordered pair
    as key. A missing key raises KeyError -- an edge the caller forgot to colour must
    stop the run loudly, not quietly masquerade as some default colour.
    """
    if edge_colours is None:
        return None
    if callable(edge_colours):
        return lambda i, j: edge_colours(i, j) if i < j else edge_colours(j, i)
    return lambda i, j: edge_colours[(i, j) if i < j else (j, i)]


def _signature(node, labels, adjacency, edge_lookup):
    """What a node's neighbourhood looks like, as one string for the digest.

    Without edge colours this is byte-for-byte the signature `canonical_labels` has
    always hashed; that equality is what pins the degeneration below. With them, the
    (edge colour, neighbour colour) pairs are rendered with repr so the serialization
    stays unambiguous for any string colours.
    """
    if edge_lookup is None:
        return labels[node] + "|" + ",".join(sorted(labels[j] for j in adjacency[node]))
    return (
        labels[node]
        + "|"
        + repr(sorted((edge_lookup(node, j), labels[j]) for j in adjacency[node]))
    )


def _stabilise(labels, adjacency, edge_lookup):
    """Refine until a pass splits nothing; same loop, same stopping rule, as before."""
    for _ in range(len(labels)):
        refined = {i: _digest(_signature(i, labels, adjacency, edge_lookup)) for i in labels}
        if Counter(refined.values()) == Counter(labels.values()):
            break
        labels = refined
    return labels


def refine(numbers, edges, index, edge_colours=None):
    """Automorphism-invariant colour per atom; equal colours mean equivalent atoms.

    The engine's step 1, and the whole of `canonical_labels` (chemistry.py keeps that
    name as a thin wrapper -- network.py and perturbations.py still consume it). Without
    `edge_colours` the result is byte-identical to the pre-PLAN function, stop rule and
    signature included; tests/test_canonical.py pins that against the old code inlined.

    `edge_colours` refines across bonds as well as across atoms: the neighbour signature
    becomes the (edge colour, neighbour colour) pairs instead of the neighbour colours
    alone. It is either a callable f(i, j) or a mapping keyed by ordered atom pair
    (min(i, j), max(i, j)); every edge must be covered, and colours are coerced to str.
    This is what lets one graph carry two bond types (PLAN section 1.2 renders a reaction
    event as one edge-coloured graph with kept/broken/formed edges).
    """
    adjacency = {int(i): set() for i in index}
    for a, b in edges:
        adjacency[int(a)].add(int(b))
        adjacency[int(b)].add(int(a))
    labels = {int(i): chemical_symbols[int(numbers[i])] for i in index}
    return _stabilise(labels, adjacency, _edge_colour_lookup(edge_colours))


@dataclass(frozen=True)
class CanonicalForm:
    """The certified identity of one labelled graph. See `canonical_form`."""

    order: tuple  # canonical labelling: order[k] is the atom at canonical position k
    colours: dict  # the stable refinement BEFORE individualization (atom -> colour)
    orbits: dict  # atom -> orbit id, computed from the enumerated group, not from colours
    group: list  # automorphisms as position-indexed arrays, identity first, then sorted
    certificate: tuple  # (colours in canonical order, edge set rendered in that order)
    info: dict  # enumerated / reason / count / complete_candidates / budget fields


def canonical_form(numbers, edges, index, edge_colours=None, limit=20000):
    """The canonical form of one labelled graph, and its automorphism group.

    `numbers` maps atom index to atomic number, `edges` is an undirected edge set over
    those indices, `index` the sorted active subset (as `state.resolve_active` returns).
    Returns a `CanonicalForm`; see the module docstring for the construction.

    The certificate pairs the edge set rendered in the leaf's own labelling with that
    labelling's colour sequence. The rendered edges ALONE do not carry the plan's
    biconditional "isomorphic <=> certificate equal": water and hypofluorite both render
    the same three-position path while differing in elements, and within one search two
    leaves of one graph could order a carbon and a nitrogen into the same position,
    because digest order is content-order, not chemistry. The colour sequence closes both
    gaps: colours refine elements monotonically (distinct colours stay distinct through
    refinement, so a leaf colour determines the element), and edge colours ride along in
    the rendering.

    `limit` bounds search nodes -- vertices of the IR tree -- in `info["search_nodes"]`
    with `info["budget_unit"] == "search_nodes"`. On overflow nothing partial is
    returned as a form: `order` and `certificate` are empty, `group` is the identity
    alone, and `info` says `enumerated=False, reason="candidate_limit"`. The stable
    colours ARE still returned on overflow: refinement always completed, and a colour
    class is a reported guess, not a certified answer, so it was never the thing that
    had to fail closed.

    Determinism: the individualized class is the lexicographically first non-singleton
    colour, its members are tried in ascending index order, and the reported `order` is
    the first leaf in exploration order that achieves the minimal certificate. Same
    graph in, same everything out -- nothing depends on set or dict enumeration order.
    """
    atoms_in = [int(i) for i in index]
    slot = {atom: k for k, atom in enumerate(atoms_in)}
    adjacency = {k: set() for k in range(len(atoms_in))}
    resolve = _edge_colour_lookup(edge_colours)
    coloured = {}
    for a, b in edges:
        u, v = slot[int(a)], slot[int(b)]
        adjacency[u].add(v)
        adjacency[v].add(u)
        if resolve is not None:
            coloured[(u, v) if u < v else (v, u)] = str(resolve(int(a), int(b)))

    if resolve is not None:

        def edge_lookup(u, v):
            return coloured[(u, v) if u < v else (v, u)]

    else:
        edge_lookup = None

    ordered_edges = sorted((u, v) for u in adjacency for v in adjacency[u] if u < v)
    initial = {k: chemical_symbols[int(numbers[atom])] for k, atom in enumerate(atoms_in)}
    root = _stabilise(initial, adjacency, edge_lookup)
    root_classes = {}
    for node, colour in root.items():
        root_classes.setdefault(colour, []).append(node)
    budget = {
        "colour_class_sizes": sorted(
            (len(members) for members in root_classes.values() if len(members) > 1),
            reverse=True,
        ),
        "search_nodes": 0,
        "budget_unit": "search_nodes",
        "limit": limit,
    }

    leaves = []  # (certificate, labelling) in exploration order; bounded by the budget
    state = {"nodes": 0, "overflow": False}

    def render(labels):
        """One leaf: its certificate, and its labelling as a tuple of atom indices."""
        ordered = sorted(labels, key=lambda node: (labels[node], node))
        position = {node: k for k, node in enumerate(ordered)}
        if edge_lookup is None:
            rendered = tuple(
                sorted(
                    (min(position[u], position[v]), max(position[u], position[v]))
                    for u, v in ordered_edges
                )
            )
        else:
            rendered = tuple(
                sorted(
                    (
                        edge_lookup(u, v),
                        min(position[u], position[v]),
                        max(position[u], position[v]),
                    )
                    for u, v in ordered_edges
                )
            )
        certificate = (tuple(labels[node] for node in ordered), rendered)
        return certificate, tuple(atoms_in[node] for node in ordered)

    def search(labels):
        if state["overflow"]:
            return
        state["nodes"] += 1
        if state["nodes"] > limit:
            state["overflow"] = True
            return
        classes = {}
        for node, colour in labels.items():
            classes.setdefault(colour, []).append(node)
        if len(classes) == len(labels):
            leaves.append(render(labels))
            return
        # The plan's rule, verbatim: first non-singleton class by lexicographic colour.
        target = classes[min(c for c, members in classes.items() if len(members) > 1)]
        for node in sorted(target):
            child = dict(labels)
            child[node] = child[node] + _INDIVIDUALISED
            search(_stabilise(child, adjacency, edge_lookup))
            if state["overflow"]:
                return

    search(root)
    budget["search_nodes"] = state["nodes"]
    colours = {atoms_in[node]: colour for node, colour in root.items()}
    identity = np.asarray(atoms_in)

    if state["overflow"]:
        return CanonicalForm(
            order=(),
            colours=colours,
            orbits={atom: k for k, atom in enumerate(atoms_in)},
            group=[identity],
            certificate=(),
            info={"enumerated": False, "reason": _CANDIDATE_LIMIT, **budget},
        )
    if not leaves:
        # Unreachable: every individualization shortens the class it splits, so the tree
        # is finite and depth-bounded by len(index). Guarded anyway rather than guessed.
        raise RuntimeError("the refinement tree finished without reaching a leaf")

    best = min(certificate for certificate, _ in leaves)
    order = next(leaf_order for certificate, leaf_order in leaves if certificate == best)
    group = []
    for certificate, leaf_order in leaves:
        if certificate != best:
            continue
        # The automorphism carrying the canonical labelling onto this leaf's: it maps
        # order[k] to leaf_order[k], stored position-indexed with global atom values,
        # exactly the convention `symmetric_rmsd` and the event key consume.
        permutation = np.empty(len(atoms_in), dtype=int)
        for k, atom in enumerate(leaf_order):
            permutation[slot[order[k]]] = atom
        group.append(permutation)
    # Identity first, then lexicographic: the consumers minimise over the whole group, so
    # order carries no meaning -- which is exactly why it should not depend on the search.
    group.sort(key=lambda p: (not np.array_equal(p, identity), tuple(int(x) for x in p)))

    # Certified orbits: {g[a] for g in group} is a's whole orbit because group is the
    # complete group. Ids are assigned scanning atoms in index order, so they do not
    # depend on how the search happened to find the group members.
    orbits = {}
    next_orbit = 0
    for atom in atoms_in:
        if atom in orbits:
            continue
        for member in sorted({int(permutation[slot[atom]]) for permutation in group}):
            orbits[member] = next_orbit
        next_orbit += 1

    return CanonicalForm(
        order=order,
        colours=colours,
        orbits=orbits,
        group=group,
        certificate=best,
        info={
            "enumerated": True,
            "count": len(group),
            "complete_candidates": len(group),
            **budget,
        },
    )
