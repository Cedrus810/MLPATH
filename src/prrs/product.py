"""The product a connect-mode run was given, and every microstate's distance to it.

PLAN_2026-09-28.md item G1. This module RECORDS; nothing here routes. The product never
enters `chemical_key`, admission, a gate, or the queue -- a run given a product builds
exactly the network it would have built without one, plus the fields written from here.

Distance. With E(G) the labelled bond set of a geometry and E_P that of the product,

    d(G) = min over sigma in Aut(G_R) of |E(G) symmetric-difference sigma(E_P)|

where G_R is the reactant graph. The minimum is over the REACTANT's group, not the
product's: sigma in Aut(G_P) maps E_P onto itself by definition, so minimising over it
is the identity and was the draft's (and the plan's own K5's) mistake. What the minimum
has to absorb is the arbitrary choice made when the product was numbered -- which of
three equivalent methyl hydrogens the user wrote into the H-H bond -- and that choice
ranges over the atoms the reactant cannot tell apart, which is Aut(G_R). A trajectory
in which a different but equivalent hydrogen did the chemistry is at d = 0, not d = 4.

Bond orders are not in d (K5): the bond set is what `encode` perceives, the same object
the identity layer builds on. d is per MICROSTATE, not per chemical node: two
microstates of one node are isomorphic but may be labelled differently, and the
labelling is what d reads.

Fail closed, as the identity layer does: a product that does not share the reactant's
atom numbering is refused (no maximum-common-substructure matching, K4), and a reactant
group that does not certify inside `automorphism_limit` is refused rather than replaced
by the identity (decision 7).
"""

from dataclasses import dataclass

from .chemistry import BUDGET_REASONS, automorphisms
from .state import encode, resolve_active

DEFINITION = (
    "d = min over sigma in Aut(G_R) of |E(G) symmetric-difference sigma(E_P)|; "
    "E = labelled bond set as encode() perceives it, bond orders not included; "
    "recorded only -- the product routes nothing"
)


class ProductMismatch(ValueError):
    """The product cannot be compared with the reactant atom for atom."""


def _bond_set(structure, config):
    graph = encode(structure, config.bond_scale, active=config.active_atoms)
    return frozenset(frozenset(edge) for edge in graph.edges)


@dataclass(frozen=True)
class ProductReference:
    key: str
    images: tuple  # the distinct sigma(E_P), sigma in Aut(G_R)
    group_order: int
    bond_scale: float
    active_atoms: tuple

    def distance(self, structure):
        graph = encode(structure, self.bond_scale, active=self.active_atoms)
        edges = frozenset(frozenset(edge) for edge in graph.edges)
        return min(len(edges ^ image) for image in self.images)

    def key_difference(self, structure):
        """(broken, formed) from this geometry to the product: pairs to break, pairs to
        form. Taken against every image of E_P at the minimum distance and united, so a
        tie between equivalent numberings never hides one of them."""
        graph = encode(structure, self.bond_scale, active=self.active_atoms)
        edges = frozenset(frozenset(edge) for edge in graph.edges)
        best = min(len(edges ^ image) for image in self.images)
        broken, formed = set(), set()
        for image in self.images:
            if len(edges ^ image) == best:
                broken |= edges - image
                formed |= image - edges
        return broken, formed


def product_reference(reactant, product, registry):
    """Check the product against the reactant and fix everything d needs, once.

    `registry` supplies `key_of`, so "product-compatible" is decided by the same
    function that decides every node's identity -- a second, hand-assembled call to
    chemical_key would be a second identity waiting to drift from the first.
    """
    config = registry.config
    if len(product) != len(reactant):
        raise ProductMismatch(
            f"product has {len(product)} atoms, reactant {len(reactant)}; connect mode "
            "compares atom for atom and does no substructure matching"
        )
    if list(product.numbers) != list(reactant.numbers):
        first = next(
            i for i, (a, b) in enumerate(zip(reactant.numbers, product.numbers)) if a != b
        )
        raise ProductMismatch(
            f"product element sequence differs from the reactant's at atom {first}; "
            "the product must use the reactant's atom numbering"
        )
    group, info = automorphisms(
        reactant, config.bond_scale, config.active_atoms, config.automorphism_limit
    )
    if info.get("reason") in BUDGET_REASONS:
        raise ProductMismatch(
            "the reactant's automorphism group did not certify inside "
            f"automorphism_limit ({info['reason']}); d would be minimised over a group "
            "that was never enumerated"
        )
    index = [int(atom) for atom in resolve_active(reactant, config.active_atoms)]
    target = _bond_set(product, config)
    images = []
    for permutation in group:
        # permutation[k] is the image of atom index[k] (canonical.CanonicalForm.group);
        # the whole group is used, so the direction of that convention cannot matter.
        image_of = {atom: int(permutation[k]) for k, atom in enumerate(index)}
        image = frozenset(frozenset(image_of[a] for a in edge) for edge in target)
        if image not in images:
            images.append(image)
    return ProductReference(
        key=registry.key_of(product)["key"],
        images=tuple(images),
        group_order=len(group),
        bond_scale=config.bond_scale,
        active_atoms=config.active_atoms,
    )


PATHS_MEANING = (
    "shortest root -> product-compatible node paths, recorded only. 'observed' walks only "
    "directed edges a trial itself travelled; 'any_evidence' also walks continuations and "
    "the undirected edges of a saddle's two-sided descent. The views are never merged, and "
    "every edge lists its own attempts / continuations / ts, so a link carried by weaker "
    "evidence says so where it is. Not connected is an answer: 'closest' is the node(s) at "
    "the smallest d_to_product the run reached"
)


def product_paths(network):
    """The R -> P path of a connect run's network, in two evidence views (PLAN G4).

    Breadth-first from the root (the depth-0 node) to the nearest product-compatible
    node, neighbours in node-id order, so the same network always gives the same path.
    Reads the network and nothing else; refuses a network without exactly one root.
    """
    nodes = network.get("chemical_nodes") or []
    roots = [node["id"] for node in nodes if node.get("depth") == 0]
    if len(roots) != 1:
        raise ValueError(f"a product path needs exactly one root node, found {len(roots)}")
    root = roots[0]
    targets = {node["id"] for node in nodes if node.get("product_compatible")}

    attempts, continuations, ts = {}, {}, {}
    for edge in network.get("reactions") or []:
        pair = (edge["source"], edge["target"])
        attempts.setdefault(pair, []).extend(edge.get("attempts") or [])
        continuations.setdefault(pair, []).extend(edge.get("continuations") or [])
    for channel in network.get("ts_connected_channels") or []:
        ends = channel["ends"]
        if len(set(ends)) != 2:
            continue  # a saddle between two conformers of one node joins nothing here
        a, b = ends
        for pair in ((a, b), (b, a)):
            ts.setdefault(pair, []).append(channel["ts"])

    observed = {pair for pair, ids in attempts.items() if ids}
    any_evidence = observed | {pair for pair, ids in continuations.items() if ids} | set(ts)

    def evidence(a, b):
        return {
            "from": a,
            "to": b,
            "attempts": list(attempts.get((a, b), [])),
            "continuations": list(continuations.get((a, b), [])),
            "ts": sorted(ts.get((a, b), [])),
        }

    def shortest(edges):
        neighbours = {}
        for a, b in edges:
            neighbours.setdefault(a, set()).add(b)
        previous, frontier = {root: None}, [root]
        while frontier and not targets & set(frontier):
            following = []
            for a in frontier:
                for b in sorted(neighbours.get(a, ())):
                    if b not in previous:
                        previous[b] = a
                        following.append(b)
            frontier = following
        reached = sorted(targets & set(frontier))
        if not reached:
            return {"connected": False, "edges": []}
        chain, node = [], reached[0]
        while previous[node] is not None:
            chain.append(evidence(previous[node], node))
            node = previous[node]
        return {"connected": True, "target": reached[0], "edges": chain[::-1]}

    distances = {
        node["id"]: min(state["d_to_product"] for state in node["microstates"])
        for node in nodes
        if node.get("microstates")
    }
    best = min(distances.values()) if distances else None
    return {
        "observed": shortest(observed),
        "any_evidence": shortest(any_evidence),
        "closest": {
            "nodes": sorted(n for n, d in distances.items() if d == best),
            "d_to_product": best,
        },
        "meaning": PATHS_MEANING,
    }
