"""Batched AFIR paths and batched (constrained) minimizations for small molecules on an MLP
(design docs/NONHESSIAN_AFIR_DESIGN_2026-11-14.md §0, "MLP-specific 2"; user 2026-11-14:
"mlp 探索的时候 是不同于 afir 我们可以大量并发 因为 分子量太小了").

GRRM pays a QM gradient per structure, so its AFIR logic is sequential. With an MLP and
small molecules one batched GPU call evaluates hundreds of structures, and the per-path
Python optimizer loop becomes the bottleneck (2026-11-08: 8 concurrent HS tasks gave only a
4.3x speed-up with the GPU not saturated). Here every path / relaxation of a wave is a row
of one coordinate array; each step is ONE ``evaluate(frames)`` call (e.g.
``parallel.MACEBatchEvaluator``) followed by a vectorised FIRE update (ASE FIRE defaults,
forces as accelerations, a per-row step cap). E/F only, no Hessian.

- ``run_paths``: AFIR paths (minimize E + sign * alpha * rho_AB per row), with the same
  per-path records and stops as ``afir_path.run_afir_path``.
- ``minimize``: unbiased minimizations, optionally each at a fixed coordinate value
  (forces projected off grad q, q restored by Newton steps along grad q, as
  ``probe.FixBondCombo``); returns energies, forces, reaction forces and bond graphs.

Every structure evaluation is charged once to ``ledger`` (one segment per call), and the
per-row forward counts are returned, so ledger total == sum of rows.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from ase.data import covalent_radii

from .afir import alpha_eV_A
from .afir_path import PathResult, _Peaks
from .ledger import BudgetExhausted

FIRE = dict(dt=0.1, dtmax=1.0, n_min=5, finc=1.1, fdec=0.5, astart=0.1, fa=0.99)


def edges_of(numbers, positions, bond_scale):
    """Bond graph by distance < bond_scale * (R_i + R_j) (the `prrs.state.encode` edge
    rule without hysteresis), for one structure."""
    x = np.asarray(positions, float)
    d = np.linalg.norm(x[:, None] - x[None, :], axis=-1)
    r = covalent_radii[numbers]
    i, j = np.nonzero(np.triu(d < bond_scale * (r[:, None] + r[None, :]), k=1))
    return sorted([int(a), int(b)] for a, b in zip(i, j))


@dataclass
class _Fire:
    n: int
    shape: tuple
    maxstep: float
    v: np.ndarray = field(init=False)
    dt: np.ndarray = field(init=False)
    a: np.ndarray = field(init=False)
    since: np.ndarray = field(init=False)

    def __post_init__(self):
        self.v = np.zeros((self.n,) + self.shape)
        self.dt = np.full(self.n, FIRE["dt"])
        self.a = np.full(self.n, FIRE["astart"])
        self.since = np.zeros(self.n, int)

    def step(self, rows, f):
        """Displacements for ``rows`` given their forces f (len(rows), N, 3)."""
        v = self.v[rows]
        p = np.einsum("kij,kij->k", f, v)
        fn = np.sqrt(np.einsum("kij,kij->k", f, f))
        vn = np.sqrt(np.einsum("kij,kij->k", v, v))
        up = p > 0
        mix = np.where(fn > 0, vn / np.where(fn > 0, fn, 1.0), 0.0)
        a = self.a[rows]
        v = np.where(
            up[:, None, None], (1 - a)[:, None, None] * v + (a * mix)[:, None, None] * f, 0.0
        )
        since = np.where(up, self.since[rows] + 1, 0)
        dt = self.dt[rows]
        grow = up & (since > FIRE["n_min"])
        dt = np.where(grow, np.minimum(dt * FIRE["finc"], FIRE["dtmax"]), dt)
        a = np.where(grow, a * FIRE["fa"], a)
        dt = np.where(up, dt, dt * FIRE["fdec"])
        a = np.where(up, a, FIRE["astart"])
        v = v + dt[:, None, None] * f
        dr = dt[:, None, None] * v
        norm = np.sqrt(np.einsum("kij,kij->k", dr, dr))
        dr = dr * np.minimum(1.0, self.maxstep / np.where(norm > 0, norm, 1.0))[:, None, None]
        self.v[rows], self.dt[rows], self.a[rows], self.since[rows] = v, dt, a, since
        return dr


def _frames(template, xs):
    out = []
    for x in xs:
        a = template.copy()
        a.set_constraint()
        a.positions = x
        out.append(a)
    return out


# Atoms per model call. POLAR-1 float32 on a 32 GB RTX 5090 ran out of memory at 128 x 19
# atoms (N1) and 256 x 16 (b05) but not at 128 x 16 (benchmarks/afir_mc/throughput.py,
# 2026-11-14); throughput was already flat from 8-32 structures, so smaller calls cost
# nothing measurable.
MAX_ATOMS_PER_CALL = int(__import__("os").environ.get("AFIR_MAX_ATOMS_PER_CALL", "1024"))


def _evaluate(evaluate, template, xs, ledger):
    """Batched calls of at most MAX_ATOMS_PER_CALL atoms; one forward charged per
    structure before the calls."""
    if ledger.budget is not None and ledger.remaining < len(xs):
        raise BudgetExhausted(f"{len(xs)} forwards requested, {ledger.remaining} left")
    for _ in xs:
        ledger.charge()
    frames = _frames(template, xs)
    per = max(1, MAX_ATOMS_PER_CALL // len(template))
    res = []
    for i in range(0, len(frames), per):
        res += evaluate(frames[i : i + per])
    e = np.array([r[0] for r in res], float)
    f = np.stack([np.asarray(r[1], float) for r in res])
    return e, f


def _ratio_min(numbers, x):
    d = np.linalg.norm(x[:, None] - x[None, :], axis=-1)
    r = covalent_radii[numbers]
    iu = np.triu_indices(len(x), 1)
    return float((d / (r[:, None] + r[None, :]))[iu].min())


def _components(n, edges):
    parent = list(range(n))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i, j in edges:
        parent[find(i)] = find(j)
    return [find(i) for i in range(n)]


def run_paths(template, starts, coordinates, evaluate, config, settings, ledger, task):
    """AFIR paths for every row (start positions, its WeightedDistance), batched.

    Returns one ``PathResult`` per row with the same frame fields and statuses as
    ``afir_path.run_afir_path`` (frames carry positions, unbiased energy, rho, bond graph,
    biased / unbiased fmax)."""
    cfg = settings
    alpha = alpha_eV_A(cfg.gamma_kJ)
    k_sign = cfg.sign * alpha
    numbers = template.numbers
    n = len(starts)
    x = np.stack([np.asarray(s, float) for s in starts])
    results = [
        PathResult(alpha_eV_A=alpha, provenance={"algorithm": "AFIR-path-EF-batch-v1",
                   "optimizer": "FIRE", "fire": FIRE, "settings": cfg.__dict__,
                   "coordinate": {"a": list(c.a), "b": list(c.b), "p": c.p}})
        for c in coordinates
    ]  # fmt: skip
    peaks = [_Peaks(cfg.noise_eV) for _ in range(n)]
    fire = _Fire(n, x.shape[1:], cfg.maxstep_A)
    active = list(range(n))
    rho0, comp0 = [None] * n, [None] * n
    spent = np.zeros(n, int)
    try:
        with ledger.segment(task=task, stage="afir_paths_batch", attempt=0) as seg:
            for step in range(cfg.max_steps + 1):
                if not active:
                    break
                try:
                    e, f = _evaluate(evaluate, template, x[active], ledger)
                except BudgetExhausted:
                    for r in active:
                        results[r].status = "budget_exhausted"
                    active = []
                    break
                spent[active] += 1
                keep, fb_all = [], []
                for idx, r in enumerate(active):
                    c = coordinates[r]
                    rho = c.value(x[r])
                    fb = f[idx] - k_sign * c.gradient(x[r])
                    edges = edges_of(numbers, x[r], config.bond_scale)
                    frame = {
                        "step": step,
                        "energy_eV": float(e[idx]),
                        "rho_A": float(rho),
                        "fmax_biased": float(np.linalg.norm(fb, axis=1).max()),
                        "fmax_unbiased": float(np.linalg.norm(f[idx], axis=1).max()),
                        "edges": edges,
                        "positions": x[r].copy(),
                    }
                    res = results[r]
                    prev = res.frames[-1] if res.frames else None
                    res.frames.append(frame)
                    if rho0[r] is None:
                        rho0[r] = rho
                        comp0[r] = len(set(_components(len(numbers), edges)))
                    peaks[r].add(step, float(e[idx]))
                    status = None
                    if _ratio_min(numbers, x[r]) < cfg.domain_min_ratio or (
                        prev is not None and prev["energy_eV"] - e[idx] > cfg.domain_max_drop_eV
                    ):
                        status = "domain_exit"
                    elif frame["fmax_biased"] <= cfg.fmax_eV_A:
                        status = "biased_minimum"
                    elif (
                        cfg.stop_after_maxima and len(peaks[r].maxima) >= cfg.stop_after_maxima
                    ):
                        status = "two_maxima"
                    elif step == cfg.early_discard_steps and cfg.sign * (rho - rho0[r]) > 0:
                        status = "not_reactive"
                    else:
                        lab = _components(len(numbers), edges)
                        if len(set(lab)) > comp0[r]:
                            d = np.linalg.norm(x[r][:, None] - x[r][None, :], axis=-1)
                            rr = covalent_radii[numbers]
                            cross = np.array(lab)[:, None] != np.array(lab)[None, :]
                            if (d / (rr[:, None] + rr[None, :]))[
                                cross
                            ].min() > cfg.dissociation_ratio:
                                status = "dissociated"
                    if status is None and step == cfg.max_steps:
                        status = "max_steps"
                    if status is None:
                        keep.append(r)
                        fb_all.append(fb)
                    else:
                        res.status = status
                if keep:
                    x[keep] += fire.step(keep, np.stack(fb_all))
                active = keep
            seg.outcome = f"{n} paths"
    finally:
        for r in range(n):
            results[r].maxima = list(peaks[r].maxima)
            results[r].spent = int(spent[r])
    return results


@dataclass
class Relaxed:
    positions: np.ndarray
    energy_eV: float
    forces: np.ndarray
    converged: bool
    steps: int
    forwards: int
    q: float | None = None
    reaction_force: float | None = None
    fmax_projected: float | None = None
    edges: list | None = None
    error: str | None = None


def minimize(template, starts, evaluate, config, ledger, task, *, fmax, max_steps,
             maxstep=0.1, constraints=None, tol=1e-8, newton_iter=50):  # fmt: skip
    """Batched FIRE minimizations of the unbiased E. ``constraints``: per row None or
    (coordinate, q): forces projected off grad q, q restored by Newton steps along grad q
    after each move (as probe.FixBondCombo). Converged: (projected) fmax <= ``fmax``."""
    numbers = template.numbers
    n = len(starts)
    x = np.stack([np.asarray(s, float) for s in starts])
    cons = constraints or [None] * n

    def restore(r):
        c = cons[r]
        if c is None:
            return
        coord, q = c
        for _ in range(newton_iter):
            err = coord.value(x[r]) - q
            if abs(err) <= tol:
                return
            g = coord.gradient(x[r])
            x[r] -= err * g / np.vdot(g, g)

    for r in range(n):
        restore(r)
    fire = _Fire(n, x.shape[1:], maxstep)
    out = [None] * n
    spent = np.zeros(n, int)
    active = list(range(n))
    with ledger.segment(task=task, stage="minimize_batch", attempt=0) as seg:
        for step in range(max_steps + 1):
            if not active:
                break
            e, f = _evaluate(evaluate, template, x[active], ledger)
            spent[active] += 1
            keep, fp_all = [], []
            for idx, r in enumerate(active):
                fr = f[idx]
                rf = fproj = None
                fp = fr
                if cons[r] is not None:
                    g = cons[r][0].gradient(x[r])
                    rf = -float(np.vdot(fr, g) / np.vdot(g, g))
                    fp = fr + rf * g
                    fproj = float(np.linalg.norm(fp, axis=1).max())
                fm = float(np.linalg.norm(fp, axis=1).max())
                done = fm <= fmax
                if done or step == max_steps:
                    out[r] = Relaxed(
                        positions=x[r].copy(),
                        energy_eV=float(e[idx]),
                        forces=fr.copy(),
                        converged=bool(done),
                        steps=step,
                        forwards=int(spent[r]),
                        q=None if cons[r] is None else cons[r][0].value(x[r]),
                        reaction_force=rf,
                        fmax_projected=fproj,
                        edges=edges_of(numbers, x[r], config.bond_scale),
                    )
                else:
                    keep.append(r)
                    fp_all.append(fp)
            if keep:
                x[keep] += fire.step(keep, np.stack(fp_all))
                for r in keep:
                    restore(r)
            active = keep
        seg.outcome = f"{n} relaxations"
    return out


def confirm_minima(template, starts, evaluate, config, ledger, task, *, rng, n_perturb=3,
                   perturb_A=0.03, max_steps=1000):  # fmt: skip
    """Batched E/F minima of ``starts`` and their perturbation checks (as
    afir_local.relax_min): every perturbed re-descent must return to the same bond graph,
    within basin_energy_eV and basin_rmsd_A (aligned). Returns one dict per row."""
    from ase import Atoms

    from prrs.state import aligned_rmsd

    fmax = config.quench_fmax_eV_A
    first = minimize(template, starts, evaluate, config, ledger, f"{task}|descent",
                     fmax=fmax, max_steps=max_steps)  # fmt: skip
    trial_rows, trial_x = [], []
    for r, m in enumerate(first):
        if m.converged:
            for _ in range(n_perturb):
                trial_rows.append(r)
                trial_x.append(m.positions + rng.normal(0.0, perturb_A, m.positions.shape))
    trials = (
        minimize(
            template,
            trial_x,
            evaluate,
            config,
            ledger,
            f"{task}|perturb",
            fmax=fmax,
            max_steps=max_steps,
        )  # fmt: skip
        if trial_x
        else []
    )
    out = []
    for r, m in enumerate(first):
        if not m.converged:
            out.append({"status": "descent_not_converged", "positions": None,
                        "forwards": m.forwards})  # fmt: skip
            continue
        a = Atoms(numbers=template.numbers, positions=m.positions)
        checks, fw = [], m.forwards
        for t, row in zip(trials, trial_rows):
            if row != r:
                continue
            fw += t.forwards
            b = Atoms(numbers=template.numbers, positions=t.positions)
            checks.append(
                bool(
                    t.converged
                    and t.edges == m.edges
                    and abs(t.energy_eV - m.energy_eV) <= config.basin_energy_eV
                    and aligned_rmsd(a, b, config.active_atoms) <= config.basin_rmsd_A
                )
            )
        out.append({"status": "confirmed" if all(checks) else "unstable",
                    "positions": m.positions, "energy_eV": m.energy_eV, "edges": m.edges,
                    "trials": checks, "forwards": fw})  # fmt: skip
    return out


TS_MIN_STEP_A = 0.05


def bracket_frames(frames, top, min_step=TS_MIN_STEP_A):
    """Path indices to relax for a bracket: from ``top`` back to the start, then from
    ``top`` forward to the end, each at least ``min_step`` in rho from the last kept one
    in its direction (the order is the walking order of afir_local.ts_from_path)."""
    rhos = [f["rho_A"] for f in frames]
    back, fwd = [top], []
    last = rhos[top]
    for k in range(top - 1, -1, -1):
        if abs(rhos[k] - last) >= min_step:
            back.append(k)
            last = rhos[k]
    last = rhos[top]
    for k in range(top + 1, len(frames)):
        if abs(rhos[k] - last) >= min_step:
            fwd.append(k)
            last = rhos[k]
    return back, fwd


def with_coordinate(frames, coordinate):
    """Path frames re-labelled with another coordinate's value (``rho_A``), for brackets
    on a TS coordinate decoupled from the push coordinate (user 2026-11-15)."""
    return [f | {"rho_A": coordinate.value(f["positions"])} for f in frames]


def pick_bracket(points_back, points_fwd):
    """First neighbouring pair (walking back from the top, then forward) that brackets a
    maximum of the rho branch: with increasing rho, dE/drho goes from + to -."""

    def ok(a, b):
        lo, hi = sorted((a, b), key=lambda p: p.q_target)
        return (
            a.converged
            and b.converged
            and hi.q_target - lo.q_target > 1e-6
            and lo.reaction_force > 0 > hi.reaction_force
        )

    for seq in (points_back, [points_back[0]] + points_fwd if points_back else points_fwd):
        for a, b in zip(seq, seq[1:]):
            if ok(a, b):
                return tuple(sorted((a, b), key=lambda p: p.q_target))
    return None


def scan_points(template, frames_list, coordinates, evaluate, config, ledger, task, *,
                fmax, max_steps):  # fmt: skip
    """Batched fixed-rho relaxations of path frames into ScanPoints. ``frames_list``: per
    path, a list of path frames; returns per path the list of ScanPoints (same order)."""
    from .probe import ScanPoint

    rows, starts, cons = [], [], []
    for p, (frames, coord) in enumerate(zip(frames_list, coordinates)):
        for j, fr in enumerate(frames):
            rows.append((p, j))
            starts.append(fr["positions"])
            cons.append((coord, coord.value(np.asarray(fr["positions"], float))))
    if not starts:
        return [[] for _ in frames_list]
    relaxed = minimize(template, starts, evaluate, config, ledger, task, fmax=fmax,
                       max_steps=max_steps, maxstep=0.05, constraints=cons)  # fmt: skip
    out = [[None] * len(fr) for fr in frames_list]
    for i, ((p, j), m) in enumerate(zip(rows, relaxed)):
        out[p][j] = ScanPoint(
            delta_A=0.0,
            q_target=float(cons[i][1]),
            q_actual=float(m.q),
            energy_eV=m.energy_eV,
            reaction_force=float(m.reaction_force),
            fmax_projected=float(m.fmax_projected),
            converged=m.converged,
            steps=m.steps,
            edges=m.edges,
            positions=m.positions,
            forces_pes=m.forces,
        )
    return out


def expand_brackets(template, items, evaluate, config, ledger, task, *, fmax, max_steps):
    """Brackets grown outward from each path's peak, stopping at the first bracket
    (user 2026-11-15: "把 TS 括区改为从峰附近向两侧逐步展开，找到括区就停止").

    ``items``: per bracket search (frames, back, fwd, coordinate), with ``back`` / ``fwd``
    from ``bracket_frames`` (back[0] is the peak). Round 0 relaxes the peak; round k the
    k-th frame on each side still available. All pending searches of a round are relaxed
    at their own fixed coordinate in ONE batched call (``scan_points``). After each round
    the pairs seen so far are checked with ``pick_bracket`` (converged neighbours, dE/dq
    from + to - with increasing q); a search stops at its first bracket or when both
    sides are exhausted. Returns per item (points_back, points_fwd, pair or None), the
    points in walking order, only the relaxed ones."""
    n = len(items)
    back_pts = [[] for _ in range(n)]
    fwd_pts = [[] for _ in range(n)]
    pairs = [None] * n
    pending = list(range(n))
    k = 0
    while pending:
        rows, frames_list, coords = [], [], []
        for i in pending:
            frames, back, fwd, coord = items[i]
            todo = []
            if k < len(back):
                todo.append(("back", back[k]))
            if k >= 1 and k - 1 < len(fwd):
                todo.append(("fwd", fwd[k - 1]))
            if todo:
                rows.append((i, [side for side, _ in todo]))
                frames_list.append([frames[f] for _, f in todo])
                coords.append(coord)
        if not rows:
            break
        pts = scan_points(template, frames_list, coords, evaluate, config, ledger, task,
                          fmax=fmax, max_steps=max_steps)  # fmt: skip
        still = []
        for (i, sides), got in zip(rows, pts):
            for side, p in zip(sides, got):
                (back_pts if side == "back" else fwd_pts)[i].append(p)
            pairs[i] = pick_bracket(back_pts[i], fwd_pts[i])
            frames, back, fwd, _ = items[i]
            more = k + 1 < len(back) or k < len(fwd)
            if pairs[i] is None and more:
                still.append(i)
        pending = still
        k += 1
    return [(back_pts[i], fwd_pts[i], pairs[i]) for i in range(n)]


NARROW_POINTS = 8
NARROW_ROUNDS = 3
NARROW_BRANCH_RMSD_A = 0.05  # same-branch test (EFSettings.branch_rmsd_A)


def narrow_brackets(template, coordinates, pairs, evaluate, config, ledger, task, *, fmax,
                    max_steps, points=NARROW_POINTS, rounds=NARROW_ROUNDS,
                    branch_rmsd=NARROW_BRANCH_RMSD_A):  # fmt: skip
    """Shrink every (lo, hi) ScanPoint bracket of a rho-branch maximum in batched rounds.

    Each round places ``points`` equally spaced q inside every bracket and relaxes each at
    fixed q from BOTH bracket ends (2 x points rows per bracket, all brackets in one batched
    call). A q counts only when both relaxations converged, agree in the sign of dE/dq and
    lie on one branch (aligned RMSD <= ``branch_rmsd``). The new bracket is the first
    neighbouring pair, in increasing q over the old ends and the agreeing interior points,
    whose dE/dq goes from + to -. A bracket whose sign change sits between disagreeing
    points (hysteresis) keeps its last width. Returns (pairs, history per bracket).
    E/F only; the HF kernel then works on a bracket where the branch is nearly linear in
    dE/dq instead of a 0.1 A wide one (b04 2026-11-14: 1 of 3 kernel runs passed its
    paired reaction-force check on the raw path bracket)."""
    from .probe import ScanPoint
    from .ts_ef import _rmsd

    current = [tuple(p) for p in pairs]
    history = [[] for _ in pairs]
    live = list(range(len(pairs)))
    for rnd in range(rounds):
        if not live:
            break
        rows, starts, cons = [], [], []
        for b in live:
            lo, hi = current[b]
            qs = np.linspace(lo.q_target, hi.q_target, points + 2)[1:-1]
            for q in qs:
                for side, src in (("lo", lo), ("hi", hi)):
                    rows.append((b, float(q), side))
                    starts.append(np.asarray(src.positions, float))
                    cons.append((coordinates[b], float(q)))
        relaxed = minimize(template, starts, evaluate, config, ledger, f"{task}|narrow{rnd}",
                           fmax=fmax, max_steps=max_steps, maxstep=0.05, constraints=cons)  # fmt: skip
        by = {}
        for (b, q, side), m in zip(rows, relaxed):
            by.setdefault((b, q), {})[side] = m
        nxt = []
        for b in live:
            lo, hi = current[b]
            seq = [lo]
            for (bb, q), sides in sorted(by.items(), key=lambda kv: kv[0][1]):
                if bb != b:
                    continue
                a, c = sides["lo"], sides["hi"]
                agree = (
                    a.converged
                    and c.converged
                    and np.sign(a.reaction_force) == np.sign(c.reaction_force)
                    and _rmsd(a.positions, c.positions) <= branch_rmsd
                )
                history[b].append({"round": rnd, "q": q, "rf_lo": a.reaction_force,
                                   "rf_hi": c.reaction_force, "agree": bool(agree)})  # fmt: skip
                if agree:
                    seq.append(ScanPoint(
                        delta_A=0.0, q_target=q, q_actual=float(a.q), energy_eV=a.energy_eV,
                        reaction_force=float(a.reaction_force),
                        fmax_projected=float(a.fmax_projected), converged=True, steps=a.steps,
                        edges=a.edges, positions=a.positions, forces_pes=a.forces,
                    ))  # fmt: skip
            seq.append(hi)
            found = None
            for x, y in zip(seq, seq[1:]):
                if x.reaction_force > 0 > y.reaction_force:
                    found = (x, y)
                    break
            if found is not None and (found[0] is not lo or found[1] is not hi):
                current[b] = found
                nxt.append(b)
        live = nxt
    return current, history
