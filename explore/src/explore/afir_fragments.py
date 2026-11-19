"""SC-AFIR fragment generation (design docs/NONHESSIAN_AFIR_DESIGN_2026-11-14.md §5.2; user
2026-11-14 chose it over raising gamma: "b   走起来").

From AFIR (Chem. Rec. 2016 §2.4; J. Comput. Chem. 2018), for an atom pair (k, l):
  1. perturbed structure: the k-l distance reduced by x % (push) or increased by x % (pull),
     x = 20; the other atoms placed by least squares so that the distance matrix stays as
     close as possible to the original one (no energy evaluation);
  2. on the perturbed structure: k in A, l in B; every atom connected to A (B) joins A (B),
     connected = r_ij / (R_i + R_j) < 1.5; done twice (first and second layer); finally
     i in A and j in B are both removed when r_ij / (R_i + R_j) < r_kl / (R_k + R_l);
  3. pairs with r_kl / (R_k + R_l) above the cut-off are skipped (user 2026-11-14: 3.0,
     the MLP does not see the papers' 5.0 / 8.0); symmetry-equivalent pairs are done once.
Pure geometry; R = ase covalent radii.
"""

from __future__ import annotations

import numpy as np
from ase.data import covalent_radii
from scipy.optimize import minimize as scipy_minimize

CONNECT_RATIO = 1.5
LAYERS = 2
PERTURB_FRACTION = 0.20
PAIR_RATIO_MAX = 3.0
_KL_WEIGHT = 1.0e4  # holds the k-l distance at its target in the least-squares fit


def ratios(numbers, positions):
    x = np.asarray(positions, float)
    d = np.linalg.norm(x[:, None] - x[None, :], axis=-1)
    r = covalent_radii[np.asarray(numbers)]
    return d / (r[:, None] + r[None, :]), d


def perturbed_structure(numbers, positions, k, l, sign=+1, fraction=PERTURB_FRACTION, seed=0):  # noqa: E741 (papers: pair (k, l))
    """Least-squares structure with d_kl scaled by (1 - fraction) (push, sign +1) or
    (1 + fraction) (pull, sign -1), every other distance kept as close as possible to the
    original. A 1e-3 A seeded jitter breaks exact symmetry (the papers' linear P-Q-R case:
    a symmetric start would stay linear)."""
    return perturbed_structure_pairs(numbers, positions, [(k, l)], sign, fraction, seed)


def perturbed_structure_pairs(numbers, positions, pairs, sign=+1, fraction=PERTURB_FRACTION,
                              seed=0):  # fmt: skip
    """``perturbed_structure`` for several pairs at once: every listed distance scaled in
    the SAME least-squares fit, so the start does not depend on an order of the pairs
    (two-pair SC paths, user 2026-11-15). One pair gives exactly perturbed_structure."""
    x0 = np.asarray(positions, float)
    n = len(x0)
    _, d0 = ratios(numbers, x0)
    target = d0.copy()
    w = np.ones((n, n))
    for k, l in pairs:  # noqa: E741
        target[k, l] = target[l, k] = d0[k, l] * (1.0 - sign * fraction)
        w[k, l] = w[l, k] = _KL_WEIGHT
    iu = np.triu_indices(n, 1)

    def f(flat):
        x = flat.reshape(n, 3)
        diff = x[:, None] - x[None, :]
        d = np.linalg.norm(diff, axis=-1)
        np.fill_diagonal(d, 1.0)
        err = (d - target) * w
        np.fill_diagonal(err, 0.0)
        val = float(np.sum(((d - target) ** 2 * w)[iu]))
        g = 2.0 * np.sum((err / d)[..., None] * diff, axis=1)
        return val, g.ravel()

    rng = np.random.default_rng(seed)
    start = x0 + rng.normal(0.0, 1e-3, x0.shape)
    res = scipy_minimize(f, start.ravel(), jac=True, method="L-BFGS-B",
                         options={"maxiter": 2000, "gtol": 1e-10})  # fmt: skip
    return res.x.reshape(n, 3)


def sc_fragments(numbers, positions, k, l):  # noqa: E741
    """Fragments A, B around (k, l) on (perturbed) ``positions`` (rules in module doc)."""
    rr, _ = ratios(numbers, positions)
    n = len(rr)
    connected = rr < CONNECT_RATIO
    np.fill_diagonal(connected, False)
    a, b = {int(k)}, {int(l)}
    for _ in range(LAYERS):
        grow_a = {j for i in a for j in np.nonzero(connected[i])[0]} - b
        grow_b = {j for i in b for j in np.nonzero(connected[i])[0]} - a - grow_a
        a |= {int(j) for j in grow_a}
        b |= {int(j) for j in grow_b}
    ref = rr[k, l]
    drop_a = {i for i in a for j in b if i != k and rr[i, j] < ref}
    drop_b = {j for i in a for j in b if j != l and rr[i, j] < ref}
    a -= drop_a
    b -= drop_b
    if not a or not b or a & b or k not in a or l not in b:
        raise ValueError(f"bad fragments for ({k}, {l}): A={sorted(a)} B={sorted(b)}")
    del n
    return tuple(sorted(a)), tuple(sorted(b))


def candidate_pairs(atoms, config, ratio_max=PAIR_RATIO_MAX):
    """Unordered atom pairs with r/(R_i+R_j) <= ratio_max, one per symmetry orbit under the
    structure's bond-graph automorphisms (prrs.chemistry.automorphisms)."""
    from prrs.chemistry import automorphisms

    rr, _ = ratios(atoms.numbers, atoms.positions)
    group, _ = automorphisms(
        atoms, config.bond_scale, config.active_atoms, config.automorphism_limit
    )
    group = [np.asarray(g) for g in group] or [np.arange(len(atoms))]
    seen, out = set(), []
    n = len(atoms)
    for i in range(n):
        for j in range(i + 1, n):
            if rr[i, j] > ratio_max:
                continue
            orbit = {tuple(sorted((int(g[i]), int(g[j])))) for g in group}
            key = min(orbit)
            if key in seen:
                continue
            seen.add(key)
            out.append((i, j))
    return out


def two_contact_combos(atoms, config, ratio_max=PAIR_RATIO_MAX):
    """Two-pair SC combinations (user 2026-11-15): two non-bonded contacts (k, l), (m, n),
    each with r/(R_i+R_j) <= ratio_max, on four different atoms; no ring-size or other
    chemical filter; one combination per symmetry orbit of the pair of pairs under the
    structure's bond-graph automorphisms. Bonds by the edges_of rule (bond_scale)."""
    from prrs.chemistry import automorphisms

    from .afir_batch import edges_of

    rr, _ = ratios(atoms.numbers, atoms.positions)
    bonded = {tuple(e) for e in edges_of(atoms.numbers, atoms.positions, config.bond_scale)}
    n = len(atoms)
    contacts = [(i, j) for i in range(n) for j in range(i + 1, n)
                if rr[i, j] <= ratio_max and (i, j) not in bonded]  # fmt: skip
    group, _ = automorphisms(
        atoms, config.bond_scale, config.active_atoms, config.automorphism_limit
    )
    group = [np.asarray(g) for g in group] or [np.arange(n)]

    def canon(p, q, g):
        a = tuple(sorted((int(g[p[0]]), int(g[p[1]]))))
        b = tuple(sorted((int(g[q[0]]), int(g[q[1]]))))
        return tuple(sorted((a, b)))

    seen, out = set(), []
    for x in range(len(contacts)):
        for y in range(x + 1, len(contacts)):
            p, q = contacts[x], contacts[y]
            if len({*p, *q}) < 4:
                continue
            key = min(canon(p, q, g) for g in group)
            if key in seen:
                continue
            seen.add(key)
            out.append((p, q))
    return out
