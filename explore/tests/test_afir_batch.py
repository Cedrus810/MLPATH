"""Batched AFIR paths and relaxations (explore.afir_batch) on the TwoSaddle bench: every row
behaves like its own optimization, the counts close, and the batched bracket hands the HF
kernel a pair that locates the bench saddle. CPU only."""

import numpy as np
import pytest
from twosaddle import SADDLE_Y, TwoSaddleInvariant, local_coords, two_saddle_atoms

from explore.afir import FrozenWeightedDistance, WeightedDistance
from explore.afir_batch import (
    bracket_frames,
    confirm_minima,
    edges_of,
    expand_brackets,
    minimize,
    pick_bracket,
    run_paths,
    scan_points,
    with_coordinate,
)
from explore.afir_path import PathSettings
from explore.ledger import Ledger
from explore.ts_branch import BranchSettings, continue_branches_ef
from explore.ts_ef import EFSettings
from prrs.config import SearchConfig

CONFIG = SearchConfig(closed_shell_only=False)
EF = EFSettings(reaction_force_sigma=7e-4, midpoint_transverse_fmax_eV_A=1e-3)


class Batch:
    """Batched evaluator over the bench; counts structures and calls."""

    def __init__(self):
        self.structures = 0
        self.calls = 0

    def __call__(self, frames):
        self.calls += 1
        out = []
        for a in frames:
            self.structures += 1
            b = a.copy()
            b.calc = TwoSaddleInvariant()
            out.append((b.get_potential_energy(), b.get_forces()))
        return out


STARTS = [(-1.0, 0.0), (-1.0, 0.1), (-0.95, -0.05), (-1.02, 0.03)]


def _template():
    a = two_saddle_atoms(-1.0, 0.0)
    a.info.update(charge=0, spin=1)
    return a


def _paths(tmp_path, budget=None, **over):
    t = _template()
    starts = [two_saddle_atoms(x, y).positions for x, y in STARTS]
    coords = [WeightedDistance.build(t.numbers, (3,), (1,)) for _ in starts]
    ledger = Ledger(tmp_path / "paths.jsonl", arm="T", budget=budget, replay=False)
    ev = Batch()
    res = run_paths(t, starts, coords, ev, CONFIG, PathSettings(noise_eV=1e-4, **over), ledger,
                    "wave|0")  # fmt: skip
    return t, coords, res, ledger, ev


def test_edges_of_matches_the_encode_rule():
    from prrs.state import encode

    a = two_saddle_atoms(-1.0, 0.0)
    a.positions[3] = a.positions[1] + [0.6, 0.0, 0.0]  # H bonded to O1
    assert edges_of(a.numbers, a.positions, CONFIG.bond_scale) == sorted(
        [int(i), int(j)] for i, j in encode(a, CONFIG.bond_scale).edges
    )


def test_batched_paths_cross_the_barrier_and_counts_close(tmp_path):
    _, coords, res, ledger, ev = _paths(tmp_path)
    assert all(r.status in ("biased_minimum", "two_maxima") for r in res)
    assert all(local_coords(r.frames[-1]["positions"])[0] > 0.5 for r in res)
    assert all(r.maxima for r in res)
    assert sum(r.spent for r in res) == ledger.spent == ev.structures
    assert all(r.spent == len(r.frames) for r in res)
    assert ev.calls == max(r.spent for r in res)  # one batched call per step
    for r, c in zip(res, coords):
        assert all(abs(f["rho_A"] - c.value(f["positions"])) < 1e-9 for f in r.frames)


def test_batched_paths_respect_the_budget(tmp_path):
    _, _, res, ledger, ev = _paths(tmp_path, budget=30)
    assert ledger.spent <= 30 and ledger.spent == ev.structures
    assert any(r.status == "budget_exhausted" for r in res)


def test_minimize_free_and_at_fixed_coordinate(tmp_path):
    t = _template()
    ledger = Ledger(tmp_path / "m.jsonl", arm="T", replay=False)
    ev = Batch()
    starts = [two_saddle_atoms(-0.9, 0.1).positions, two_saddle_atoms(0.8, -0.1).positions]
    free = minimize(t, starts, ev, CONFIG, ledger, "m", fmax=1e-3, max_steps=2000)
    assert all(m.converged for m in free)
    xs = sorted(local_coords(m.positions)[0] for m in free)
    assert xs[0] == pytest.approx(-1.0, abs=1e-2) and xs[1] == pytest.approx(1.0, abs=1e-2)
    coord = WeightedDistance.build(t.numbers, (3,), (1,))
    x0 = two_saddle_atoms(-0.5, -0.3).positions
    q = coord.value(x0)
    held = minimize(t, [x0], ev, CONFIG, ledger, "c", fmax=1e-3, max_steps=2000,
                    maxstep=0.05, constraints=[(coord, q)])[0]  # fmt: skip
    assert held.converged and abs(held.q - q) < 1e-6
    assert held.fmax_projected <= 1e-3 and held.reaction_force is not None
    assert ledger.spent == ev.structures == sum(m.forwards for m in free) + held.forwards


def test_batched_bracket_locates_the_bench_saddle(tmp_path):
    t, coords, res, _, _ = _paths(tmp_path)
    ledger = Ledger(tmp_path / "b.jsonl", arm="T", replay=False)
    ev = Batch()
    chosen = [bracket_frames(r.frames, r.maxima[0]) for r in res]
    frames_list = [[r.frames[k] for k in b + f] for r, (b, f) in zip(res, chosen)]
    pts = scan_points(t, frames_list, coords, ev, CONFIG, ledger, "s",
                      fmax=EF.transverse_fmax_eV_A, max_steps=EF.max_relax_steps)  # fmt: skip
    assert ledger.spent == ev.structures
    found = 0
    for r, (b, f), p, coord in zip(res, chosen, pts, coords):
        pair = pick_bracket(p[: len(b)], p[len(b) :])
        assert pair is not None
        branch = continue_branches_ef(
            t.copy(), coord, pair, lambda: TwoSaddleInvariant(), CONFIG,
            Ledger(tmp_path / f"k{found}.jsonl", arm="T", replay=False), "k",
            branch_settings=BranchSettings(geometry_noise_A=3.3e-4), settings=EF,
            ceiling_eV=5.0,
        )  # fmt: skip
        assert branch.status == "converged_candidate", branch.status
        x, y, _ = local_coords(branch.refine["positions"])
        assert abs(x) < 1e-2 and min(abs(y - s) for s in SADDLE_Y) < 1e-2
        found += 1
    assert found == len(res)


def test_expanding_brackets_match_the_full_scan_with_fewer_relaxations(tmp_path):
    """User 2026-11-15: grow the bracket outward from the peak and stop at the first one.
    Same bracket as the full scan (pick_bracket), fewer relaxed points, every round of
    all searches in batched calls, and the kernel still locates the bench saddle."""
    t, coords, res, _, _ = _paths(tmp_path)
    chosen = [bracket_frames(r.frames, r.maxima[0]) for r in res]
    full_ledger = Ledger(tmp_path / "full.jsonl", arm="T", replay=False)
    frames_list = [[r.frames[k] for k in b + f] for r, (b, f) in zip(res, chosen)]
    full = scan_points(t, frames_list, coords, Batch(), CONFIG, full_ledger, "s",
                       fmax=EF.transverse_fmax_eV_A, max_steps=EF.max_relax_steps)  # fmt: skip
    ledger = Ledger(tmp_path / "e.jsonl", arm="T", replay=False)
    ev = Batch()
    items = [(r.frames, b, f, c) for r, (b, f), c in zip(res, chosen, coords)]
    out = expand_brackets(t, items, ev, CONFIG, ledger, "e",
                          fmax=EF.transverse_fmax_eV_A, max_steps=EF.max_relax_steps)  # fmt: skip
    assert ledger.spent == ev.structures
    relaxed = sum(len(pb) + len(pf) for pb, pf, _ in out)
    assert relaxed < sum(len(p) for p in full)
    for (pb, pf, pair), p, (b, f), coord in zip(out, full, chosen, coords):
        ref = pick_bracket(p[: len(b)], p[len(b) :])
        assert pair is not None and ref is not None
        assert [q.q_target for q in pair] == pytest.approx([q.q_target for q in ref])
        branch = continue_branches_ef(
            t.copy(), coord, pair, lambda: TwoSaddleInvariant(), CONFIG,
            Ledger(tmp_path / f"k{id(pair)}.jsonl", arm="T", replay=False), "k",
            branch_settings=BranchSettings(geometry_noise_A=3.3e-4), settings=EF,
            ceiling_eV=5.0,
        )  # fmt: skip
        assert branch.status == "converged_candidate", branch.status
        x, y, _ = local_coords(branch.refine["positions"])
        assert abs(x) < 1e-2 and min(abs(y - s) for s in SADDLE_Y) < 1e-2


def test_expanding_brackets_without_a_sign_change_exhaust_both_sides(tmp_path):
    t, coords, res, _, _ = _paths(tmp_path)
    r = res[0]
    # frames before the peak only on the uphill side: no + -> - pair can form
    back, _ = bracket_frames(r.frames, r.maxima[0])
    items = [(r.frames, back[-2:], [], coords[0])]
    ledger = Ledger(tmp_path / "x.jsonl", arm="T", replay=False)
    (pb, pf, pair), = expand_brackets(t, items, Batch(), CONFIG, ledger, "x",
                                      fmax=EF.transverse_fmax_eV_A,
                                      max_steps=EF.max_relax_steps)  # fmt: skip
    assert pair is None and len(pb) == 2 and pf == []


@pytest.mark.parametrize("b_group, expect", [((0, 2), "saddle"), ((0, 1, 2), "no_bracket")])
def test_decoupled_ts_coordinate(tmp_path, b_group, expect):
    """User 2026-11-15: SC pushes only the selected pair (here H-O1, V_bias = alpha r) and
    the TS bracket / kernel use a different multi-pair coordinate, q over H x B with
    weights frozen at the peak frame.
    B = {O0, O2}: q is monotone along the path; brackets on the re-labelled frames and the
      same kernel give the bench saddle.
    B = {O0, O1, O2}: H moves from O0 to O1, one distance grows as the other shrinks, and q
      (c = 0.51 / 0.47 / 0.03) has its minimum at the peak, so dE/dq < 0 on both sides and
      no bracket forms -- a limit of this TS coordinate, reported, never a false TS."""
    t, coords, res, _, _ = _paths(tmp_path)
    ledger = Ledger(tmp_path / "d.jsonl", arm="T", replay=False)
    ev = Batch()
    items = []
    for r in res:
        top = r.maxima[0]
        ts = FrozenWeightedDistance.build(t.numbers, (3,), b_group, r.frames[top]["positions"])
        frames = with_coordinate(r.frames, ts)
        assert all(abs(f["rho_A"] - ts.value(f["positions"])) < 1e-12 for f in frames)
        b, f = bracket_frames(frames, top)
        items.append((frames, b, f, ts))
    out = expand_brackets(t, items, ev, CONFIG, ledger, "d",
                          fmax=EF.transverse_fmax_eV_A, max_steps=EF.max_relax_steps)  # fmt: skip
    assert ledger.spent == ev.structures
    for (pb, pf, pair), (_, _, _, ts) in zip(out, items):
        if expect == "no_bracket":
            assert pair is None
            continue
        assert pair is not None
        branch = continue_branches_ef(
            t.copy(), ts, pair, lambda: TwoSaddleInvariant(), CONFIG,
            Ledger(tmp_path / f"dk{id(pair)}.jsonl", arm="T", replay=False), "k",
            branch_settings=BranchSettings(geometry_noise_A=3.3e-4), settings=EF,
            ceiling_eV=5.0,
        )  # fmt: skip
        assert branch.status == "converged_candidate", branch.status
        x, y, _ = local_coords(branch.refine["positions"])
        assert abs(x) < 1e-2 and min(abs(y - s) for s in SADDLE_Y) < 1e-2


def test_confirm_minima_batched(tmp_path):
    t, _, res, _, _ = _paths(tmp_path)
    ledger = Ledger(tmp_path / "c.jsonl", arm="T", replay=False)
    ev = Batch()
    ends = [r.frames[-1]["positions"] for r in res]
    out = confirm_minima(t, ends, ev, CONFIG, ledger, "min", rng=np.random.default_rng(0))
    assert all(o["status"] == "confirmed" for o in out), [o["status"] for o in out]
    assert sum(o["forwards"] for o in out) == ledger.spent == ev.structures
    for o in out:  # converged to quench_fmax 0.03 eV/A; the bench is soft along y
        x, y, _ = local_coords(o["positions"])
        assert abs(x - 1.0) < 1e-2 and abs(y) < 0.05


def test_narrowed_brackets_shrink_and_still_locate_the_saddle(tmp_path):
    from explore.afir_batch import narrow_brackets

    t, coords, res, _, _ = _paths(tmp_path)
    ledger = Ledger(tmp_path / "n.jsonl", arm="T", replay=False)
    ev = Batch()
    chosen = [bracket_frames(r.frames, r.maxima[0]) for r in res]
    frames_list = [[r.frames[k] for k in b + f] for r, (b, f) in zip(res, chosen)]
    pts = scan_points(t, frames_list, coords, ev, CONFIG, ledger, "s",
                      fmax=EF.transverse_fmax_eV_A, max_steps=EF.max_relax_steps)  # fmt: skip
    pairs = [pick_bracket(p[: len(b)], p[len(b) :]) for p, (b, f) in zip(pts, chosen)]
    calls_before = ev.calls
    narrowed, hist = narrow_brackets(t, coords, pairs, ev, CONFIG, ledger, "n",
                                     fmax=EF.transverse_fmax_eV_A,
                                     max_steps=EF.max_relax_steps, rounds=2)  # fmt: skip
    assert ev.calls - calls_before < 2 * 2 * EF.max_relax_steps  # batched, not per bracket
    assert ledger.spent == ev.structures
    for old, new, h in zip(pairs, narrowed, hist):
        w0 = old[1].q_target - old[0].q_target
        w1 = new[1].q_target - new[0].q_target
        assert w1 < w0 / 20 and new[0].reaction_force > 0 > new[1].reaction_force
        assert old[0].q_target <= new[0].q_target < new[1].q_target <= old[1].q_target
        assert h and all("agree" in e for e in h)
    branch = continue_branches_ef(
        t.copy(), coords[0], narrowed[0], lambda: TwoSaddleInvariant(), CONFIG,
        Ledger(tmp_path / "k.jsonl", arm="T", replay=False), "k",
        branch_settings=BranchSettings(geometry_noise_A=3.3e-4), settings=EF, ceiling_eV=5.0,
    )  # fmt: skip
    assert branch.status == "converged_candidate", branch.status
    x, y, _ = local_coords(branch.refine["positions"])
    assert abs(x) < 1e-2 and min(abs(y - s) for s in SADDLE_Y) < 1e-2
