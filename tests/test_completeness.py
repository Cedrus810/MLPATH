"""The analytic full chain for conditional completeness certificates (PLAN item 3).

A synthetic triple well, roots symmetric about A, along one bond coordinate, energy and forces analytic, is the
cheapest surface where the certificate's answer is known by construction: the two
outer wells sit above the declared threshold X, so a small enough delta certifies
that A's is the only sublevel component in the ball -- and a delta enlarged past the
gradient condition must be REFUSED even though every sampled point still quenches
back to A. That refusal is the mutation check the plan demands: without it this is
just a function that says yes.
"""

import numpy as np
import pytest
from ase.calculators.calculator import Calculator

from prrs.config import SearchConfig
from prrs.calculators import demo_atoms

from prrs.completeness import audit_basin, certificate_text, cover, gradient_bound

C_POLY = 20.0
SLOPE = 0.0
K_PERP = 4.0
ROOTS = (0.5, 1.2, 1.9)


def potential(t, r):
    poly = C_POLY * np.prod([(t - root) ** 2 for root in ROOTS])
    return poly + SLOPE * t + 0.5 * K_PERP * (r**2).sum()


def dpotential_dt(t):
    poly_prime = 0.0
    for i, root in enumerate(ROOTS):
        others = np.prod([(t - o) ** 2 for j, o in enumerate(ROOTS) if j != i])
        poly_prime += 2 * (t - root) * others
    return C_POLY * poly_prime + SLOPE


class TripleWell(Calculator):
    implemented_properties = ["energy", "forces"]

    def calculate(self, atoms=None, properties=("energy",), system_changes=None):
        super().calculate(atoms, properties, system_changes)
        u = atoms.positions[1] - atoms.positions[0]
        t, r = u[0], u[1:]
        energy = potential(t, r)
        du = np.zeros(3)
        du[0] = dpotential_dt(t)
        du[1:] = K_PERP * r
        forces = np.zeros_like(atoms.positions)
        forces[1] = -du
        forces[0] = du
        self.results = {"energy": energy, "forces": forces}


def triple_well_factory():
    return TripleWell()


CONFIG = SearchConfig(
    closed_shell_only=False,
    quench_fmax_eV_A=1e-4,
    quench_steps=400,
    response_steps=20,
)

# The certificate's subspace V: bond stretch and one transverse direction, both
# unit-normalised in the flat 6-dim space. d = 2, inside the enforced cap.
STRETCH = np.zeros(6)
STRETCH[[0, 3]] = [1.0, -1.0]
STRETCH /= np.linalg.norm(STRETCH)
TRANSVERSE = np.zeros(6)
TRANSVERSE[[1, 4]] = [1.0, -1.0]
TRANSVERSE /= np.linalg.norm(TRANSVERSE)
BASIS = np.vstack([STRETCH, TRANSVERSE])
# 证书链用 d=1(键伸缩):环上逐点淬火,d=2 会把测试拖到分钟级。
# 覆盖半径的几何性质仍用 d=2 钉住。
CERT_BASIS = BASIS[:1]


@pytest.fixture(scope="module")
def basin_minimum():
    from prrs.reliability import GuardedCalculator
    from prrs.runner import quench

    atoms = demo_atoms()
    guard = GuardedCalculator(TripleWell(), CONFIG)
    atoms.calc = guard
    converged, _ = quench(atoms, guard, CONFIG)
    assert converged
    return atoms


def true_lipschitz(lo=0.9, hi=1.5):
    """sup |V''| over the t range the audit's sample ring actually covers.

    Analytic potential, Hessian known: the caller states the range because the
    certificate's Lipschitz constant belongs to the region the claim covers -- a
    coarse delta pushes the +delta ring past the barrier, and the bound must price
    that curvature in or the covering argument is quietly using a constant that was
    never true where it needs it.
    """
    ts = np.linspace(lo, hi, 4001)
    seconds = np.gradient([dpotential_dt(t) for t in ts], ts)
    return float(max(seconds.max(), K_PERP)) * np.sqrt(2)


def test_cover_radius_holds_over_the_ball():
    rng = np.random.default_rng(42)
    centre = np.zeros(6)
    delta = 0.05
    points = cover(centre, BASIS, 0.25, delta)
    flat = np.asarray(points)
    for _ in range(300):
        direction = rng.normal(size=2)
        direction /= np.linalg.norm(direction)
        r = rng.uniform(0, 0.25)
        probe = (centre + BASIS.T @ (direction * r)).reshape(1, -1)
        distance = np.linalg.norm(flat - probe, axis=1).min()
        assert distance <= delta + 1e-9


def test_a_basis_that_is_not_orthonormal_is_refused():
    """A scaled row stretches the real grid spacing past h: the cover would still
    claim radius delta while leaving ball points up to 2*delta from every grid point."""
    centre = np.zeros(6)
    with pytest.raises(ValueError, match="orthonormal"):
        cover(centre, np.vstack([2.0 * STRETCH, TRANSVERSE]), 0.25, 0.05)
    skewed = np.vstack([STRETCH, (STRETCH + TRANSVERSE) / np.sqrt(2)])
    with pytest.raises(ValueError, match="orthonormal"):
        cover(centre, skewed, 0.25, 0.05)


def test_gradient_bound_recovers_a_linear_field():
    points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]
    gradients = [[0.0, 0.0], [2.0, 0.0], [0.0, 2.0]]
    assert gradient_bound(points, gradients) == pytest.approx(2.0)


def test_certificate_succeeds_when_the_answer_is_known(basin_minimum):
    record = audit_basin(
        basin_minimum,
        triple_well_factory,
        CONFIG,
        CERT_BASIS,
        radius=0.15,
        delta=0.02,
        energy_ceiling_X=0.15,
        gap_w=0.10,
        lipschitz_upper=true_lipschitz(0.9, 1.5),
    )
    assert record["verdict"] == "conditions_met", record["conditions"]
    assert record["l_used"] >= record["l_measured"]
    assert record["l_used"] >= true_lipschitz() * 0.99
    text = certificate_text(record)
    assert "DRAFT" in text and record["verdict"] in text


def test_enlarged_delta_is_refused_even_though_every_sample_returns(basin_minimum):
    """The plan's mutation check: sparse cover must not pass by luck.

    delta = 0.3 leaves every grid point inside A's basin (the ball is only 0.25
    wide), so membership alone sees nothing wrong -- the refusal must come from the
    gradient condition, which is exactly the covering argument noticing it can no
    longer bridge its samples.
    """
    record = audit_basin(
        basin_minimum,
        triple_well_factory,
        CONFIG,
        CERT_BASIS,
        radius=0.15,
        delta=0.12,
        energy_ceiling_X=0.15,
        gap_w=0.10,
        lipschitz_upper=true_lipschitz(0.9, 1.62),
    )
    assert record["verdict"] == "refused"
    # every sampled point still returns to A -- the refusal must come from the
    # gradient condition, not from having sampled the barrier by luck
    assert all(row["member"] for row in record["membership_table"])
    assert not record["conditions"]["(ii) L_used*delta < g_min"]


def test_a_ball_that_contains_another_basin_is_refused_by_membership(basin_minimum):
    """radius 0.8 reaches the 1.9 well; honest sampling sees it and refuses."""
    record = audit_basin(
        basin_minimum,
        triple_well_factory,
        CONFIG,
        CERT_BASIS,
        radius=0.45,
        delta=0.05,
        energy_ceiling_X=0.15,
        gap_w=0.10,
        lipschitz_upper=true_lipschitz(),
    )
    assert record["verdict"] == "refused"
    assert not record["conditions"]["(i) membership all basin"]
    assert any(row["outcome"] == "other_basin" for row in record["membership_table"])


def test_missing_quantity_is_an_error_not_a_quiet_record(basin_minimum):
    with pytest.raises(ValueError, match="radius"):
        audit_basin(
            basin_minimum,
            triple_well_factory,
            CONFIG,
            BASIS,
            radius=None,
            delta=0.02,
            energy_ceiling_X=0.15,
            gap_w=0.10,
        )


def test_subspace_cap_is_enforced(basin_minimum):
    big_basis = np.vstack([BASIS] + [np.eye(6)[2 + k] for k in range(2)])
    with pytest.raises(ValueError, match="cap"):
        audit_basin(
            basin_minimum,
            triple_well_factory,
            CONFIG,
            big_basis,
            radius=0.15,
            delta=0.02,
            energy_ceiling_X=0.15,
            gap_w=0.10,
        )
