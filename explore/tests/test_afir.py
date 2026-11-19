"""Non-Hessian AFIR, phase 1 pieces (docs/PLAN_AFIR_MC_2026-11-14.md §4): the AFIR
coordinate and function, one AFIR path, and the project's own MIN / TS / connection on the
analytic TwoSaddle bench. CPU only, no MLP."""

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from twosaddle import SADDLE_Y, TwoSaddleInvariant, local_coords, two_saddle_atoms

from explore.afir import (
    AFIRCalculator,
    FrozenWeightedDistance,
    WeightedDistance,
    alpha_eV_A,
    alpha_from_gamma,
    gamma_from_T_t,
)
from explore.afir_local import link, path_minima, relax_min, ts_from_path
from explore.afir_path import PathSettings, run_afir_path
from explore.ledger import Ledger, view_factory
from explore.ts_branch import BranchSettings
from explore.ts_ef import EFSettings
from prrs.config import SearchConfig

CONFIG = SearchConfig(closed_shell_only=False)
# TwoSaddle bench numerics, as d2_bench.HF (measured on this bench 2026-10-26)
EF = EFSettings(reaction_force_sigma=7e-4, midpoint_transverse_fmax_eV_A=1e-3)
BRANCH = BranchSettings(geometry_noise_A=3.3e-4)


def _molecule(seed=0):
    rng = np.random.default_rng(seed)
    return Atoms("CCHOHN", positions=rng.normal(0.0, 1.5, (6, 3)))


def test_rho_gradient_matches_finite_differences_and_is_invariant():
    a = _molecule()
    rho = WeightedDistance.build(a.numbers, (0, 2, 3), (1, 4, 5))
    x = a.positions
    g = rho.gradient(x)
    h = 1e-6
    num = np.zeros_like(x)
    for i in range(len(x)):
        for k in range(3):
            xp, xm = x.copy(), x.copy()
            xp[i, k] += h
            xm[i, k] -= h
            num[i, k] = (rho.value(xp) - rho.value(xm)) / (2 * h)
    assert np.allclose(g, num, atol=1e-7)
    q, _ = np.linalg.qr(np.random.default_rng(1).normal(size=(3, 3)))
    moved = x @ q.T + np.array([0.3, -2.0, 1.0])
    assert rho.value(moved) == pytest.approx(rho.value(x), abs=1e-12)
    assert rho.pairs == ((0, 1), (0, 4), (0, 5), (2, 1), (2, 4), (2, 5), (3, 1), (3, 4), (3, 5))
    with pytest.raises(ValueError):
        WeightedDistance.build(a.numbers, (0, 1), (1, 2))


def test_rho_is_dominated_by_close_pairs():
    """w = ((R_i+R_j)/r)^6: a far pair barely moves rho (the point of the AFIR weight)."""
    a = Atoms("CCC", positions=[[0, 0, 0], [1.6, 0, 0], [20.0, 0, 0]])
    near = WeightedDistance.build(a.numbers, (0,), (1,)).value(a.positions)
    both = WeightedDistance.build(a.numbers, (0,), (1, 2)).value(a.positions)
    assert near == pytest.approx(1.6)
    assert both == pytest.approx(1.6, abs=1e-5)


def test_frozen_weights_coordinate():
    """Design §2.1 (user 2026-11-15): q = sum c_ij r_ij, c frozen at x0. Gradient matches
    finite differences, q is rigid-motion invariant, equals rho_AB at x0, and every pair
    derivative dq/dr_ij = c_ij > 0 at any geometry (rho_AB pulls pairs with r > 1.2 rho apart)."""
    a = _molecule()
    x0 = a.positions
    A, B = (0, 2, 3), (1, 4, 5)
    q = FrozenWeightedDistance.build(a.numbers, A, B, x0)
    rho = WeightedDistance.build(a.numbers, A, B)
    assert q.value(x0) == pytest.approx(rho.value(x0), abs=1e-12)
    assert q.c.sum() == pytest.approx(1.0) and (q.c > 0).all()
    x = x0 + np.random.default_rng(3).normal(0.0, 0.4, x0.shape)
    g = q.gradient(x)
    h = 1e-6
    num = np.zeros_like(x)
    for i in range(len(x)):
        for k in range(3):
            xp, xm = x.copy(), x.copy()
            xp[i, k] += h
            xm[i, k] -= h
            num[i, k] = (q.value(xp) - q.value(xm)) / (2 * h)
    assert np.allclose(g, num, atol=1e-7)
    rot, _ = np.linalg.qr(np.random.default_rng(1).normal(size=(3, 3)))
    assert q.value(x @ rot.T + 1.0) == pytest.approx(q.value(x), abs=1e-12)
    # contrast at a geometry where rho_AB pulls a far pair apart: q never does
    far = x.copy()
    far[B[0]] += 3.0 * (x[B[0]] - x[A[0]]) / np.linalg.norm(x[B[0]] - x[A[0]])
    _, r, w = rho._parts(far)
    drho = w * (6.0 * rho.value(far) / r - 5.0)
    assert (drho < 0).any() and (q.c > 0).all()
    assert q.pairs == rho.pairs


def test_alpha_and_gamma_reproduce_the_paper_numbers():
    # Chem. Rec. 2016: 298.15 K -> 93.3 (1 h), 101.2 (1 day), 106.9 (10 days) kJ/mol
    for t, gamma in ((3600, 93.3), (86400, 101.2), (864000, 106.9)):
        assert gamma_from_T_t(298.15, t) == pytest.approx(gamma, abs=0.05)
    assert alpha_from_gamma(100.0) == pytest.approx(118.8, abs=0.1)
    assert alpha_eV_A(100.0) == pytest.approx(1.231, abs=0.001)
    with pytest.raises(ValueError):
        alpha_from_gamma(0.0)


def test_afir_calculator_forces_are_the_gradient_of_the_biased_energy():
    a = two_saddle_atoms(-0.6, 0.3, 0.05)
    rho = WeightedDistance.build(a.numbers, (3,), (1, 2))
    calc = AFIRCalculator(TwoSaddleInvariant(), rho, 0.8, +1)
    a.calc = calc
    f = a.get_forces()
    e0 = TwoSaddleInvariant()
    b = a.copy()
    b.calc = e0
    assert calc.unbiased["energy"] == pytest.approx(b.get_potential_energy())
    assert a.get_potential_energy() == pytest.approx(
        b.get_potential_energy() + 0.8 * rho.value(a.positions)
    )
    h = 1e-5
    x = a.positions.copy()
    for i, k in ((3, 0), (3, 1), (1, 2), (0, 0)):
        xp, xm = x.copy(), x.copy()
        xp[i, k] += h
        xm[i, k] -= h
        a.positions = xp
        ep = a.get_potential_energy()
        a.positions = xm
        em = a.get_potential_energy()
        assert -(ep - em) / (2 * h) == pytest.approx(f[i, k], abs=1e-4)
    with pytest.raises(ValueError):
        AFIRCalculator(TwoSaddleInvariant(), rho, 0.8, 0)


def _counted(tmp_path, name, budget=None):
    ledger = Ledger(tmp_path / f"{name}.jsonl", arm="T", budget=budget, replay=False)
    return ledger, view_factory(TwoSaddleInvariant(), ledger)


def _path(tmp_path, **over):
    """Push H (atom 3) towards O1 (atom 1) from the left minimum (the D2 bench intent)."""
    a = two_saddle_atoms(-1.0, 0.0)
    a.info.update(charge=0, spin=1)
    rho = WeightedDistance.build(a.numbers, (3,), (1,))
    ledger, factory = _counted(tmp_path, "path", over.pop("budget", None))
    settings = PathSettings(noise_eV=1e-4, **over)
    res = run_afir_path(a, rho, factory, CONFIG, settings, ledger, "path|0")
    return a, rho, res, ledger, factory


def test_afir_path_crosses_the_barrier_and_records_every_step(tmp_path):
    _, rho, res, ledger, _ = _path(tmp_path)
    assert res.status in ("biased_minimum", "two_maxima")
    assert res.maxima, "the path passes an energy maximum"
    x_end, _, _ = local_coords(res.frames[-1]["positions"])
    assert x_end > 0.5  # product side of the bench
    assert len(res.frames) == ledger.spent  # one unbiased forward per step
    assert res.spent == ledger.spent
    assert all(abs(f["rho_A"] - rho.value(f["positions"])) < 1e-9 for f in res.frames)
    top = res.maxima[0]
    assert res.frames[top]["energy_eV"] > res.frames[0]["energy_eV"]


def test_afir_path_guards(tmp_path):
    _, _, res, _, _ = _path(tmp_path / "a", domain_min_ratio=50.0)
    assert res.status == "domain_exit" and len(res.frames) == 1
    _, _, res, ledger, _ = _path(tmp_path / "b", budget=15)
    assert res.status == "budget_exhausted" and ledger.spent <= 15
    _, _, res, _, _ = _path(tmp_path / "c", max_steps=3)
    assert res.status == "max_steps"
    with pytest.raises(ValueError):
        PathSettings(noise_eV=1e-4, sign=0)


@pytest.fixture
def located(tmp_path):
    a, rho, res, _, _ = _path(tmp_path)
    ledger, factory = _counted(tmp_path, "ts")
    status, points, branch = ts_from_path(
        a,
        rho,
        res.frames,
        res.maxima[0],
        factory,
        CONFIG,
        ledger,
        "ts|0",
        ef_settings=EF,
        branch_settings=BRANCH,
        ceiling_eV=res.frames[0]["energy_eV"] + 5.0,
    )
    return a, res, status, points, branch, ledger, tmp_path


def test_ts_from_the_path_top_is_a_bench_saddle(located):
    _, _, status, points, branch, ledger, _ = located
    assert status == "converged_candidate", status
    found = branch.refine
    x, y, _ = local_coords(found["positions"])
    assert abs(x) < 1e-2
    assert min(abs(y - s) for s in SADDLE_Y) < 1e-2
    f = TwoSaddleInvariant()
    b = two_saddle_atoms(0, 0)
    b.positions = np.asarray(found["positions"])
    b.calc = f
    assert np.linalg.norm(b.get_forces(), axis=1).max() <= 0.01
    assert any(p.reaction_force > 0 for p in points) and any(
        p.reaction_force < 0 for p in points
    )
    assert branch.spent > 0 and ledger.spent >= branch.spent


def test_link_connects_two_minima(located):
    a, _, status, _, branch, _, tmp_path = located
    assert status == "converged_candidate"
    ledger, factory = _counted(tmp_path, "link")
    conn = link(
        a, branch, factory, CONFIG, ledger, "link|0", ef_settings=EF, sigma=7e-4, seed=3
    )
    assert conn.status == "connected", conn.status
    xs = sorted(local_coords(sd["descent"]["positions"])[0] for sd in conn.sides)
    assert xs[0] < -0.9 and xs[1] > 0.9


def test_path_minima_and_min_confirmation(tmp_path):
    a, _, res, _, _ = _path(tmp_path)
    mins = path_minima(res.frames, 1e-4)
    assert mins[-1] == len(res.frames) - 1
    ledger, factory = _counted(tmp_path, "min")
    out = relax_min(
        a,
        res.frames[-1]["positions"],
        factory,
        CONFIG,
        ledger,
        "min|0",
        rng=np.random.default_rng(0),
    )
    assert out["status"] == "confirmed", out
    x, y, _ = local_coords(out["positions"])
    assert abs(x - 1.0) < 1e-2 and abs(y) < 1e-2


def test_displacement_coordinate():
    """Design §4.2 in distance space: q_d = sum u_ij (r_ij - r_ij(before)), u ~ dr over all
    pairs; q_d(before) = 0, q_d(after) = |dr|, exact rigid-motion invariance, gradient
    matches finite differences."""
    from explore.afir import DisplacementCoordinate

    rng = np.random.default_rng(5)
    before = rng.normal(0.0, 1.5, (6, 3))
    rot, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    after = (before + rng.normal(0.0, 0.2, before.shape)) @ rot.T + 2.0
    q = DisplacementCoordinate.build(before, after)
    assert q.value(before) == pytest.approx(0.0, abs=1e-12)
    assert q.value(after) == pytest.approx(q.length_A)
    assert np.linalg.norm(q.u) == pytest.approx(1.0)
    x = before + rng.normal(0.0, 0.1, before.shape)
    rot2, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    assert q.value(x @ rot2.T + np.array([0.7, -1.0, 3.0])) == pytest.approx(q.value(x))
    g = q.gradient(x)
    h = 1e-6
    for i, k in ((0, 0), (3, 2), (5, 1)):
        xp, xm = x.copy(), x.copy()
        xp[i, k] += h
        xm[i, k] -= h
        assert (q.value(xp) - q.value(xm)) / (2 * h) == pytest.approx(g[i, k])
    assert len(q.pairs) == 15
    with pytest.raises(ValueError):
        DisplacementCoordinate.build(before, before @ rot.T)


class _SeparableFold(Calculator):
    """The fold of test_ts_branch with x on the 0-1 distance and y on the 2-3 distance of
    two separated diatomics: the jump displacement (pure y) is orthogonal to x in
    Cartesian space, so q_d = y exactly."""

    implemented_properties = ["energy", "forces"]

    def __init__(self, calls):
        super().__init__()
        self.calls = calls

    def get_hessian(self, *args, **kwargs):
        raise AssertionError("Hessian entered E/F search")

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        self.calls.append(atoms.positions.copy())
        p = atoms.positions
        v, w = p[1] - p[0], p[3] - p[2]
        r, s = np.linalg.norm(v), np.linalg.norm(w)
        x, y = (r - 1.8) / 0.6, (s - 2.0) / 0.4
        energy = 0.4 * ((y * y - 1) ** 2 + 2 * x * y + 0.5 * x * x)
        dr, ds = 0.4 * (2 * y + x) / 0.6, 0.4 * (4 * y * (y * y - 1) + 2 * x) / 0.4
        f = np.zeros_like(p)
        f[0] += dr * v / r
        f[1] -= dr * v / r
        f[2] += ds * w / s
        f[3] -= ds * w / s
        self.results = {"energy": energy, "forces": f}


def _sep_frame(x, y):
    r, s = 1.8 + 0.6 * x, 2.0 + 0.4 * y
    return Atoms("C4", positions=[[0, 0, 0], [r, 0, 0], [0, 0, 6.0], [0, s, 6.0]])


def _fold_xy(p, separable):
    p = np.asarray(p)
    r = np.linalg.norm(p[1] - p[0])
    s = np.linalg.norm((p[3] - p[2]) if separable else (p[2] - p[1]))
    return (r - 1.8) / 0.6, (s - 2.0) / 0.4


@pytest.mark.parametrize(
    "separable",
    [
        pytest.param(
            True,
            marks=pytest.mark.xfail(
                strict=True,
                reason="2026-11-15: the jump displacement between the two independently "
                "walked branch structures also changes the cross-diatomic distances "
                "(u 0-3 -0.36, 1-2 +0.29 beside 2-3 -0.885); they cost no energy on this "
                "bench, so the fixed-q_d relaxation moves them and slides to the global "
                "minimum. The Cartesian form failed the same way through rotations.",
            ),
        ),
        False,
    ],
)
def test_jump_is_rebracketed_on_the_displacement_coordinate(tmp_path, separable):
    """Fold bench: the saddle (0, 0) lies on neither fixed-x branch, so the kernel on the
    push coordinate (x) ends in hysteresis, both sides jumping in y.
    separable: x and y are distances of two separated diatomics, the jump displacement is
      pure y in Cartesian space, the q_d branch passes the saddle and the same kernel
      locates it (design §4.2).
    collinear (test_ts_branch's C3): the Cartesian jump displacement (atom 2 along the
      chain, centred) overlaps the x stretch, q_d ~ 0.6 x + 0.8 y in bench units, and
      g.H^-1.g = 0.85 > 0 at the saddle: the q_d branch misses it too. Expected:
      branch_unresolved and no candidate -- never a false TS."""
    from test_ts_branch import BRANCH, INTENT, SETTINGS, Surface, frame, point

    from explore.afir_local import rebracket_jump
    from explore.ledger import counted_factory
    from explore.ts_branch import continue_branches_ef

    calls = []
    ledger = Ledger(tmp_path / "ledger.jsonl", arm="HF", budget=200000)
    if separable:
        factory = counted_factory(lambda: _SeparableFold(calls), ledger)
        frames = [_sep_frame(-1.2, 1.2416), _sep_frame(1.2, -1.2416)]
    else:
        factory = counted_factory(lambda: Surface(calls, "fold"), ledger)
        frames = [frame(-1.2, 1.2416), frame(1.2, -1.2416)]
    br = continue_branches_ef(frames[0], INTENT, [point(f) for f in frames], factory,
                              CONFIG, ledger, "t", branch_settings=BRANCH,
                              settings=SETTINGS)  # fmt: skip
    assert br.status == "unresolved:hysteresis"
    out = rebracket_jump(frames[0], br, factory, CONFIG, ledger, "t|qd",
                         ef_settings=SETTINGS, branch_settings=BRANCH)  # fmt: skip
    ledger.close()
    assert ledger.spent == len(calls)
    if not separable:
        assert out["status"] == "branch_unresolved"
        assert all(j["status"] != "bracketed" for j in out["jumps"])
        return
    assert out["status"] == "converged_candidate"
    x, y = _fold_xy(out["branch"].refine["positions"], separable)
    assert x == pytest.approx(0.0, abs=1e-2)
    assert y == pytest.approx(0.0, abs=1e-2)


def test_pair_sum_and_pair_set_coordinates():
    """Two-pair push alpha (r_kl + r_mn)/sqrt 2 and the TS pair-set coordinate (user
    2026-11-15): gradients match finite differences; each push pair gets 1/sqrt 2; the
    pair set over canon(A x B) equals FrozenWeightedDistance on A, B; the union of two
    pair sets has no duplicate or self pair."""
    from explore.afir import FrozenPairSet, PairSumDistance

    a = _molecule()
    x = a.positions
    push = PairSumDistance.build([(0, 1), (2, 4)])
    r = lambda i, j: np.linalg.norm(x[j] - x[i])  # noqa: E731
    assert push.value(x) == pytest.approx((r(0, 1) + r(2, 4)) / np.sqrt(2))
    A, B = (0, 2, 3), (1, 4, 5)
    pairs = [(i, j) for i in A for j in B]
    ps = FrozenPairSet.build(a.numbers, pairs, x)
    assert ps.value(x) == pytest.approx(
        FrozenWeightedDistance.build(a.numbers, A, B, x).value(x)
    )
    y = x + np.random.default_rng(4).normal(0.0, 0.3, x.shape)
    for coord in (push, ps):
        g = coord.gradient(y)
        h = 1e-6
        for i, k in ((0, 0), (2, 1), (4, 2)):
            yp, ym = y.copy(), y.copy()
            yp[i, k] += h
            ym[i, k] -= h
            assert (coord.value(yp) - coord.value(ym)) / (2 * h) == pytest.approx(
                g[i, k], abs=1e-7
            )
    union = FrozenPairSet.canon([(1, 0), (0, 1), (4, 2), (2, 4), (3, 5)])
    assert union == ((0, 1), (2, 4), (3, 5))
