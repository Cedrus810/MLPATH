"""Bond-change window locator for AFIR paths, E/F only (2026-11-25; user "可以" to wiring the
window rule into the runner as an option; evidence DONE_LOG 2026-11-24 / 25: on b05 the same
19 paths gave 10 / 15 / 16 recovered with windows at k = 0 / 10 / 200 against 0 / 4 / 5 with
the fixed highest-node climb, and all three references at k = 10; on b09 a first-order
saddle reactant -> (2+1)).

Per path (two E/F minima as ends, the AFIR frames in between):
  1. coarse chain: the frames from the lowest of the first third on, resampled to N nodes
     between the two minima (as benchmarks/afir_mc/e2_string_climb.py);
  2. simplified String, K iterations (afir_string, batched over all paths);
  3. one window [j-1, j+2] per interval j -> j+1 where the bond graph changes, resampled
     to M nodes with its two outer nodes fixed (not minima: a local window, as in adaptive
     NEB); S String iterations (batched over all windows);
  4. fixed-direction climb from the highest interior window node, tau from its two
     neighbours, F - 2 (F.tau) tau, FIRE, full fmax <= ci_fmax, ceiling = start + rise;
  5. converged points: E/F descent from x -/+ side_step tau on both sides.
No Hessian, no curvature estimate, no chemistry in the window choice.
"""

from __future__ import annotations

import numpy as np

from .afir_batch import _evaluate, _Fire, edges_of, minimize
from .afir_string import reparametrize, resample, tangent


def bond_windows(x, numbers, bond_scale):
    """(lo, hi, j) for every interval j -> j+1 of the chain whose bond graph changes."""
    g = [edges_of(numbers, xi, bond_scale) for xi in x]
    return [
        (max(j - 1, 0), min(j + 2, len(x) - 1), j)
        for j in range(len(x) - 1)
        if g[j] != g[j + 1]
    ]


def string_batch(template, chains, evaluate, ledger, task, iterations, *, dt=0.05,
                 max_step_A=0.05):  # fmt: skip
    """Simplified String on several equal-length chains at once (ends fixed): one batched
    call per iteration for all interior nodes. Returns (chains, interior energies)."""
    x = np.stack([reparametrize(np.asarray(c, float)) for c in chains])
    n_in = x.shape[1] - 2
    with ledger.segment(task=task, stage="string", attempt=0) as seg:
        for it in range(iterations + 1):
            inner = x[:, 1:-1].reshape((-1,) + x.shape[2:])
            e, f = _evaluate(evaluate, template, inner, ledger)
            e = np.asarray(e).reshape(len(x), n_in)
            f = np.asarray(f).reshape((len(x), n_in) + x.shape[2:])
            if it == iterations:
                break
            step = dt * f
            size = np.linalg.norm(step, axis=3).max(axis=2)
            step *= np.minimum(1.0, max_step_A / np.where(size > 0, size, 1.0))[
                :, :, None, None
            ]
            moved = x.copy()
            moved[:, 1:-1] += step
            x = np.stack([reparametrize(m) for m in moved])
        seg.outcome = f"{iterations} iterations"
    return x, e


def climb_batch(template, starts, taus, evaluate, ledger, task, *, max_step_A=0.05,
                ci_fmax=0.01, ci_steps=2000, ceiling_eV=None):  # fmt: skip
    """afir_string.climb for several (start, fixed tau) at once (copied from
    benchmarks/afir_mc/n1_tangent_variants.climb_batch, unchanged)."""
    x = np.stack([np.asarray(s, float) for s in starts])
    tau = np.stack([np.asarray(t, float) / np.linalg.norm(t) for t in taus])
    n = len(x)
    fire = _Fire(n, x.shape[1:], max_step_A)
    out = [None] * n
    spent = np.zeros(n, int)
    live = list(range(n))
    with ledger.segment(task=task, stage="climb", attempt=0) as seg:
        for k in range(ci_steps + 1):
            e, f = _evaluate(evaluate, template, x[live], ledger)
            nxt, rows, fci = [], [], []
            for j, r in enumerate(live):
                spent[r] += 1
                ec, fc = float(e[j]), f[j]
                fmax = float(np.linalg.norm(fc, axis=1).max())
                status = None
                cap = ceiling_eV[r] if np.ndim(ceiling_eV) else ceiling_eV
                if cap is not None and ec > cap:
                    status = "above_ceiling"
                elif fmax <= ci_fmax:
                    status = "converged"
                elif k == ci_steps:
                    status = "max_steps"
                if status:
                    out[r] = {"ci_status": status, "ci_steps": k, "ci_positions": x[r].copy(),
                              "ci_energy_eV": ec, "ci_fmax_eV_A": fmax,
                              "ci_forwards": int(spent[r])}  # fmt: skip
                    continue
                nxt.append(r)
                rows.append(r)
                fci.append(fc - 2.0 * np.vdot(fc, tau[r]) * tau[r])
            if not nxt:
                break
            x[rows] += fire.step(rows, np.stack(fci))
            live = nxt
        seg.outcome = "done"
    return out


def locate_windows(template, paths, evaluate, config, ledger, task, *, nodes=20, k=10, m=16,
                   s=20, ci_fmax=0.01, ci_steps=2000, side_step=0.05, ceiling_rise=3.0):  # fmt: skip
    """``paths``: dicts with frames (positions array), frame_energies, start / end
    (positions, energy_eV). Returns, per path, a list of window results: window, jump,
    ci (climb dict), sides (two minimize results or None)."""
    out = [[] for _ in paths]
    if not paths:
        return out
    chains = []
    for p in paths:
        pos, en = np.asarray(p["frames"], float), np.asarray(p["frame_energies"], float)
        low = int(np.argmin(en[: len(en) // 3 + 1]))
        mids = resample(list(pos[low:]), nodes)[1:-1]
        chains.append(
            np.stack([p["start"]["positions"]] + list(mids) + [p["end"]["positions"]])
        )
    x, _ = string_batch(template, chains, evaluate, ledger, f"{task}|string", k)
    wins, wx = [], []
    for i, xi in enumerate(x):
        for lo, hi, j in bond_windows(xi, template.numbers, config.bond_scale):
            wins.append((i, lo, hi, j))
            wx.append(resample(list(xi[lo : hi + 1]), m))
    if not wins:
        return out
    wxs, we = string_batch(template, wx, evaluate, ledger, f"{task}|window", s)
    starts, taus, ceil = [], [], []
    for (i, _, _, _), xw, ew in zip(wins, wxs, we):
        top = int(np.argmax(ew)) + 1
        starts.append(xw[top])
        taus.append(tangent(xw, top))
        ceil.append(paths[i]["start"]["energy_eV"] + ceiling_rise)
    cis = climb_batch(template, starts, taus, evaluate, ledger, f"{task}|climb",
                      ci_fmax=ci_fmax, ci_steps=ci_steps, ceiling_eV=np.array(ceil))  # fmt: skip
    conv = [r for r, c in enumerate(cis) if c["ci_status"] == "converged"]
    side_starts = []
    for r in conv:
        t = taus[r] / np.linalg.norm(taus[r])
        side_starts += [cis[r]["ci_positions"] - side_step * t,
                        cis[r]["ci_positions"] + side_step * t]  # fmt: skip
    sides = (minimize(template, side_starts, evaluate, config, ledger, f"{task}|sides",
                      fmax=config.quench_fmax_eV_A, max_steps=2000)
             if side_starts else [])  # fmt: skip
    for r, ((i, lo, hi, j), ci) in enumerate(zip(wins, cis)):
        sd = None
        if r in conv:
            q = conv.index(r)
            sd = sides[2 * q : 2 * q + 2]
        out[i].append({"window": [lo, hi], "jump": j, "ci": ci, "sides": sd})
    return out
