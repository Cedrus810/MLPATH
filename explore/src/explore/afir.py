"""AFIR function pieces for the non-Hessian AFIR line (docs/NONHESSIAN_AFIR_DESIGN_2026-11-14.md
§2; plan docs/PLAN_AFIR_MC_2026-11-14.md).

Taken from AFIR (Maeda et al., Chem. Rec. 2016, 16, 2232, eqs. 1-3 and 5; J. Comput. Chem.
2018, 39, 233):

    F(Q)  = E(Q) + rho * alpha * sum_{i in A, j in B} w_ij r_ij / sum w_ij
    w_ij  = ((R_i + R_j) / r_ij)^p,  p = 6, R = covalent radii
    alpha = gamma / ([2^(-1/6) - (1 + sqrt(1 + gamma/eps))^(-1/6)] R0)
    gamma = -R T ln(h / (t k_B T))
    R0 = 3.8164 A, eps = 1.0061 kJ/mol (Ar-Ar Lennard-Jones)

rho = +1 pushes A and B together, -1 pulls them apart. Eq. (3) is garbled in both text
extractions we could read; the form above follows the papers' wording (mean force between
two Ar atoms colliding with energy gamma) and matches the readable fragments. Check it
against the original PDF (design §12) before relying on alpha to more than its size.

``WeightedDistance`` exposes ``value`` / ``gradient`` / ``pairs`` like ``probe.Intent``, so
the HF E/F kernel (``ts_ef._Constraint``, ``refine_ts_ef``, ``continue_branches_ef``) can
hold it fixed without change. Only energies and forces are used; no Hessian anywhere.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from ase.calculators.calculator import Calculator, all_changes
from ase.data import covalent_radii

KJ_PER_EV = 96.48533212  # kJ/mol per eV
R0_A = 3.8164
EPS_KJ = 1.0061
GAS_R_KJ = 8.314462618e-3  # kJ/(mol K)
K_B = 1.380649e-23
H_PLANCK = 6.62607015e-34


def gamma_from_T_t(T_K, t_s):
    """Chem. Rec. eq. (5): the barrier crossed in time t at temperature T, kJ/mol."""
    return GAS_R_KJ * T_K * math.log(t_s * K_B * T_K / H_PLANCK)


def alpha_from_gamma(gamma_kJ):
    """Chem. Rec. eq. (3) (form to be checked against the PDF), kJ/mol/A."""
    if not gamma_kJ > 0:
        raise ValueError("gamma must be positive")
    turning = (1.0 + math.sqrt(1.0 + gamma_kJ / EPS_KJ)) ** (-1.0 / 6.0)
    return gamma_kJ / ((2.0 ** (-1.0 / 6.0) - turning) * R0_A)


def alpha_eV_A(gamma_kJ):
    return alpha_from_gamma(gamma_kJ) / KJ_PER_EV


@dataclass(frozen=True)
class WeightedDistance:
    """rho_AB(x) = sum w_ij r_ij / sum w_ij over i in A, j in B, w = ((R_i+R_j)/r)^p."""

    a: tuple[int, ...]
    b: tuple[int, ...]
    radii_sum: np.ndarray = field(repr=False, compare=False)
    p: float = 6.0
    intent_id: str = "afir"

    @classmethod
    def build(cls, numbers, a, b, p=6.0, intent_id="afir"):
        a, b = tuple(int(i) for i in a), tuple(int(j) for j in b)
        if not a or not b or set(a) & set(b):
            raise ValueError("A and B must be non-empty and disjoint")
        numbers = np.asarray(numbers)
        rs = (
            covalent_radii[numbers[list(a)]][:, None]
            + covalent_radii[numbers[list(b)]][None, :]
        )
        return cls(a, b, rs, float(p), intent_id)

    @property
    def pairs(self):
        return tuple((i, j) for i in self.a for j in self.b)

    def _parts(self, positions):
        x = np.asarray(positions, float)
        d = x[list(self.b)][None, :, :] - x[list(self.a)][:, None, :]  # (|A|, |B|, 3)
        r = np.linalg.norm(d, axis=-1)
        w = (self.radii_sum / r) ** self.p
        return d, r, w

    def value(self, positions) -> float:
        _, r, w = self._parts(positions)
        return float(np.sum(w * r) / np.sum(w))

    def gradient(self, positions) -> np.ndarray:
        d, r, w = self._parts(positions)
        s0 = np.sum(w)
        rho = np.sum(w * r) / s0
        # d(w r)/dr = (1 - p) w, dw/dr = -p w / r  ->  drho/dr = w ((1 - p) + p rho / r) / s0
        dr = w * ((1.0 - self.p) + self.p * rho / r) / s0
        u = d / r[..., None]  # unit vector from i (A) to j (B)
        grad = np.zeros_like(np.asarray(positions, float))
        contrib = dr[..., None] * u
        np.add.at(grad, list(self.b), contrib.sum(axis=0))
        np.subtract.at(grad, list(self.a), contrib.sum(axis=1))
        return grad


@dataclass(frozen=True)
class FrozenWeightedDistance:
    """q(x) = sum_ij c_ij r_ij(x), i in A, j in B, with the AFIR weights frozen at x0:
    c_ij = w_ij(x0) / sum_ab w_ab(x0), w = ((R_i+R_j)/r)^p.

    SC-AFIR bias coordinate decided by the user 2026-11-15 (design §2.1): x0 is the path's
    perturbed start, c is not updated along the path (the next path takes its own start),
    so dq/dr_ij = c_ij > 0 for every pair. With rho_AB the weights follow the geometry and
    drho/dr_ij ~ w_ij (6 rho / r_ij - 5): the closest pair takes the force and pairs with
    r_ij > 1.2 rho are pulled apart (N1, docs/AFIR_PUSH_OPTIONS_2026-11-15.md). Same
    interface as WeightedDistance (value / gradient / pairs / a / b / p)."""

    a: tuple[int, ...]
    b: tuple[int, ...]
    c: np.ndarray = field(repr=False, compare=False)  # (|A|, |B|), sums to 1
    p: float = 6.0
    intent_id: str = "afir-frozen"

    @classmethod
    def build(cls, numbers, a, b, x0, p=6.0, intent_id="afir-frozen"):
        rho = WeightedDistance.build(numbers, a, b, p=p)
        _, _, w = rho._parts(x0)
        return cls(rho.a, rho.b, w / np.sum(w), float(p), intent_id)

    @property
    def pairs(self):
        return tuple((i, j) for i in self.a for j in self.b)

    def _parts(self, positions):
        x = np.asarray(positions, float)
        d = x[list(self.b)][None, :, :] - x[list(self.a)][:, None, :]
        return d, np.linalg.norm(d, axis=-1)

    def value(self, positions) -> float:
        _, r = self._parts(positions)
        return float(np.sum(self.c * r))

    def gradient(self, positions) -> np.ndarray:
        d, r = self._parts(positions)
        contrib = (self.c / r)[..., None] * d  # c_ij * unit vector i -> j
        grad = np.zeros_like(np.asarray(positions, float))
        np.add.at(grad, list(self.b), contrib.sum(axis=0))
        np.subtract.at(grad, list(self.a), contrib.sum(axis=1))
        return grad


@dataclass(frozen=True)
class DisplacementCoordinate:
    """Path-displacement coordinate of design §4.2 for a branch jump, in distance space
    (user 2026-11-15 「可以」 after the Cartesian form failed on the separable fold bench:
    a linear Cartesian coordinate is rotation invariant to first order only, and the
    fixed-q relaxation escaped it by rotating the fragments 10-82 deg):

        q_d(x) = sum_{i<j} u_ij [r_ij(x) - r_ij(x_before)],  u = dr / |dr|,
        dr_ij  = r_ij(x_after) - r_ij(x_before)  over all atom pairs.

    Exactly translation and rotation invariant, no alignment needed; q_d(before) = 0,
    q_d(after) = |dr| (``length_A``). Same interface as WeightedDistance (value /
    gradient / pairs), so the HF kernel holds it fixed unchanged; every pair distance that
    changes across the jump enters the locating coordinate."""

    i: np.ndarray = field(repr=False, compare=False)
    j: np.ndarray = field(repr=False, compare=False)
    u: np.ndarray = field(repr=False, compare=False)
    r_ref: np.ndarray = field(repr=False, compare=False)
    length_A: float = 0.0
    intent_id: str = "afir-displacement"

    @classmethod
    def build(cls, x_before, x_after, intent_id="afir-displacement"):
        xb, xa = np.asarray(x_before, float), np.asarray(x_after, float)
        i, j = np.triu_indices(len(xb), 1)
        rb = np.linalg.norm(xb[j] - xb[i], axis=1)
        dr = np.linalg.norm(xa[j] - xa[i], axis=1) - rb
        n = float(np.linalg.norm(dr))
        if not n > 1e-8:
            raise ValueError("no internal displacement between the two structures")
        return cls(i, j, dr / n, rb, n, intent_id)

    @property
    def pairs(self):
        return tuple((int(a), int(b)) for a, b in zip(self.i, self.j))

    def value(self, positions) -> float:
        x = np.asarray(positions, float)
        r = np.linalg.norm(x[self.j] - x[self.i], axis=1)
        return float(np.dot(self.u, r - self.r_ref))

    def gradient(self, positions) -> np.ndarray:
        x = np.asarray(positions, float)
        d = x[self.j] - x[self.i]
        r = np.linalg.norm(d, axis=1)
        contrib = (self.u / r)[:, None] * d
        grad = np.zeros_like(x)
        np.add.at(grad, self.j, contrib)
        np.subtract.at(grad, self.i, contrib)
        return grad


def _pair_parts(pairs, positions):
    x = np.asarray(positions, float)
    i = np.array([p[0] for p in pairs])
    j = np.array([p[1] for p in pairs])
    d = x[j] - x[i]
    return i, j, d, np.linalg.norm(d, axis=1)


def _pair_gradient(i, j, d, r, coef, shape):
    contrib = (coef / r)[:, None] * d
    grad = np.zeros(shape)
    np.add.at(grad, j, contrib)
    np.subtract.at(grad, i, contrib)
    return grad


@dataclass(frozen=True)
class PairSumDistance:
    """sum_p w_p r_p over listed atom pairs with fixed weights. Two-pair SC push (user
    2026-11-15): pairs (k, l), (m, n) on four different atoms, w = 1/sqrt(2) each, so
    V_bias = alpha (r_kl + r_mn) / sqrt(2): each pair gets alpha / sqrt(2) and the force
    along the joint direction is alpha (gamma unchanged). Interface as WeightedDistance;
    ``a`` / ``b`` list the first / second atom of each pair (records only)."""

    pair_list: tuple[tuple[int, int], ...]
    w: np.ndarray = field(repr=False, compare=False)
    p: float = 1.0
    intent_id: str = "afir-pairsum"

    @classmethod
    def build(cls, pairs, weights=None, intent_id="afir-pairsum"):
        pairs = tuple((int(i), int(j)) for i, j in pairs)
        if not pairs or any(i == j for i, j in pairs):
            raise ValueError("pairs must be non-empty and of two different atoms")
        w = np.full(len(pairs), 1.0 / np.sqrt(len(pairs))) if weights is None else weights
        return cls(pairs, np.asarray(w, float), 1.0, intent_id)

    @property
    def a(self):
        return tuple(i for i, _ in self.pair_list)

    @property
    def b(self):
        return tuple(j for _, j in self.pair_list)

    @property
    def pairs(self):
        return self.pair_list

    def value(self, positions) -> float:
        _, _, _, r = _pair_parts(self.pair_list, positions)
        return float(np.dot(self.w, r))

    def gradient(self, positions) -> np.ndarray:
        i, j, d, r = _pair_parts(self.pair_list, positions)
        return _pair_gradient(i, j, d, r, self.w, np.shape(positions))


@dataclass(frozen=True)
class FrozenPairSet:
    """TS locating coordinate (user 2026-11-15): q(x) = sum_{(i,j) in P} c_ij r_ij(x) over
    a set P of unordered atom pairs, c_ij = w_ij(x_f) / sum w(x_f), w = ((R_i+R_j)/r)^6,
    frozen at a frame x_f (the bracket's peak). For one SC pair P = canon(A x B) (same
    value as FrozenWeightedDistance on A, B); for two pairs P = canon(A_kl x B_kl) union
    canon(A_mn x B_mn) -- pairs, not atom sets, so no A/B overlap and no cross pairs the
    two groupings did not have."""

    pair_list: tuple[tuple[int, int], ...]
    c: np.ndarray = field(repr=False, compare=False)
    p: float = 6.0
    intent_id: str = "afir-ts-pairset"

    @staticmethod
    def canon(pairs):
        """Unordered, de-duplicated, sorted atom pairs (i < j)."""
        out = {tuple(sorted((int(i), int(j)))) for i, j in pairs}
        if any(i == j for i, j in out):
            raise ValueError("a pair joins an atom to itself")
        return tuple(sorted(out))

    @classmethod
    def build(cls, numbers, pairs, x_frame, p=6.0, intent_id="afir-ts-pairset"):
        pairs = cls.canon(pairs)
        if not pairs:
            raise ValueError("empty pair set")
        numbers = np.asarray(numbers)
        _, _, _, r = _pair_parts(pairs, x_frame)
        rs = np.array(
            [covalent_radii[numbers[i]] + covalent_radii[numbers[j]] for i, j in pairs]
        )
        w = (rs / r) ** p
        return cls(pairs, w / w.sum(), float(p), intent_id)

    @property
    def a(self):
        return tuple(i for i, _ in self.pair_list)

    @property
    def b(self):
        return tuple(j for _, j in self.pair_list)

    @property
    def pairs(self):
        return self.pair_list

    def value(self, positions) -> float:
        _, _, _, r = _pair_parts(self.pair_list, positions)
        return float(np.dot(self.c, r))

    def gradient(self, positions) -> np.ndarray:
        i, j, d, r = _pair_parts(self.pair_list, positions)
        return _pair_gradient(i, j, d, r, self.c, np.shape(positions))


class AFIRCalculator(Calculator):
    """E + sign * alpha * rho_AB around an unbiased calculator. The unbiased energy and
    forces of the last call are kept in ``unbiased`` (eV, eV/A); every call costs exactly
    one unbiased evaluation, charged wherever ``base`` charges it."""

    implemented_properties = ["energy", "forces"]

    def __init__(self, base, coordinate, alpha_eV_A, sign=+1, **kwargs):
        super().__init__(**kwargs)
        if sign not in (+1, -1):
            raise ValueError("sign must be +1 (push) or -1 (pull)")
        self.base, self.coord = base, coordinate
        self.alpha, self.sign = float(alpha_eV_A), int(sign)
        self.unbiased = None

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        work = self.atoms.copy()
        work.calc = self.base
        e = float(work.get_potential_energy())
        f = np.asarray(work.get_forces(), float)
        x = self.atoms.positions
        k = self.sign * self.alpha
        rho = self.coord.value(x)
        self.unbiased = {"energy": e, "forces": f.copy(), "rho": rho}
        self.results = {"energy": e + k * rho, "forces": f - k * self.coord.gradient(x)}
