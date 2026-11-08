"""v1 intents (DISCOVERY_LAYER_DESIGN_2026-10-13 §3.1): distance primitives and their
binary combinations. An intent says which distances move and in which sign; it never
names a product, a reference bond change or a saddle direction. Units are Å throughout
(v1 has only distances, so §3.2's per-family scale is 1 Å and nothing is mixed)."""

from __future__ import annotations

import itertools
from dataclasses import dataclass

from .poses import RULE_VERSION, fragments
from .probe import Intent
from .streams import stream


@dataclass(frozen=True)
class IntentSpec:
    """An intent with its sign. `Intent` itself is sign-free: every scan must take its
    deltas from `deltas`, the single entry point that applies the sign (ruling R10)."""

    intent: Intent
    sign: int
    family: str
    layer: str
    scope: str

    def deltas(self, magnitudes):
        return (0.0, *[self.sign * abs(float(m)) for m in magnitudes])


def _ident(pose_id, family, pairs, coefs, sign):
    body = ";".join(f"{i}-{j}" for i, j in pairs)
    cs = ",".join(f"{c:+g}" for c in coefs)
    return f"{pose_id}|{family}:{body}:{cs}:{'+' if sign > 0 else '-'}:{RULE_VERSION}"


def distance_intents(atoms, pose_id, config, n_far, seed):
    from prrs.chemistry import canonical_labels
    from prrs.perturbations import _hop_counts
    from prrs.state import encode

    graph = encode(atoms, config.bond_scale)
    labels = canonical_labels(atoms.numbers, graph.edges, graph.index)
    hops = _hop_counts(graph)
    owner = {i: k for k, f in enumerate(fragments(atoms, config)) for i in f.atoms}
    near, far_classes = [], {}
    excluded = {"bend_in_disguise": 0}
    n = len(atoms)
    for i, j in itertools.combinations(range(n), 2):
        bonded = (i, j) in graph.edges
        scope = "inter" if owner[i] != owner[j] else f"intra:f{owner[i]}"
        if bonded:
            near.append(((i, j), "bonded", scope))
            continue
        sep = hops.get((i, j))
        if sep is not None and sep < config.min_pair_hops:
            excluded["bend_in_disguise"] += 1
            continue
        if atoms.get_distance(i, j) <= config.pair_cutoff_A:
            near.append(((i, j), "near", scope))
        else:
            key = tuple(sorted((str(labels[i]), str(labels[j])))) + (scope.split(":")[0],)
            far_classes.setdefault(key, []).append(((i, j), "far", scope))
    classes = sorted(far_classes)
    order = stream(seed, "far", pose_id).permutation(len(classes)) if classes else []
    chosen = [classes[k] for k in list(order)[:n_far]]
    rows = near + [m for key in chosen for m in far_classes[key]]
    specs = []
    for pair, layer, scope in rows:
        for sign in (-1, 1):
            specs.append(
                IntentSpec(
                    Intent(_ident(pose_id, "distance", (pair,), (1.0,), sign), (pair,), (1.0,)),
                    sign,
                    "distance",
                    layer,
                    scope,
                )
            )
    report = {
        "excluded": excluded,
        "bonded_pairs": sum(1 for r in near if r[1] == "bonded"),
        "near_pairs": sum(1 for r in near if r[1] == "near"),
        "far_pairs_total": sum(len(v) for v in far_classes.values()),
        "far_classes_total": len(classes),
        "far_classes_sampled": len(chosen),
    }
    return specs, report


def _adjacent(p, q, edges):
    if set(p) & set(q):
        return True
    return any((min(a, b), max(a, b)) in edges for a in p for b in q)


def pair_intents(singles, atoms, config, locality, pose_id):
    """Unordered pairs of distinct single-pair intents (bonded and near layers only),
    coefficients (1, 1) and (1, -1), both signs: 4 intents per combination. Pairs are
    put in index order so (1, -1) and (-1, 1) are not both generated."""
    from prrs.state import encode

    edges = encode(atoms, config.bond_scale).edges
    base = sorted({s.intent.pairs[0] for s in singles if s.layer in ("bonded", "near")})
    scope = {s.intent.pairs[0]: s.scope for s in singles}
    specs, counts = [], {"shared_atom": 0, "adjacent": 0, "all": 0}
    for p, q in itertools.combinations(base, 2):
        shared = bool(set(p) & set(q))
        adj = _adjacent(p, q, edges)
        counts["all"] += 1
        counts["adjacent"] += adj
        counts["shared_atom"] += shared
        keep = {"shared_atom": shared, "adjacent": adj, "all": True}[locality]
        if not keep:
            continue
        cat = (
            "inter"
            if "inter" in (scope[p], scope[q])
            else ("intra" if scope[p] == scope[q] else "cross_intra")
        )
        for coefs in ((1.0, 1.0), (1.0, -1.0)):
            for sign in (-1, 1):
                specs.append(
                    IntentSpec(
                        Intent(
                            _ident(pose_id, "distance_pair", (p, q), coefs, sign), (p, q), coefs
                        ),
                        sign,
                        "distance_pair",
                        locality,
                        cat,
                    )
                )
    return specs, {"combinations_by_locality": counts, "intents": len(specs)}
