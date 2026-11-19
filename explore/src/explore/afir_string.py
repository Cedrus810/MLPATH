"""Simplified string method + climbing image as an E/F-only TS locator (2026-11-19, user:
"简化 String 方法交替做力驱动更新和路径重新分布，避免节点全挤到反应物盆；之后可以用爬坡像定位
鞍点。这样试试？"). Diagnostic, not wired into the runner.

Simplified string (W. E, W. Ren, E. Vanden-Eijnden, J. Chem. Phys. 126, 164103 (2007)):
nodes x_0..x_{N-1} with fixed end points (two minima). Each iteration
  1. evolve: every interior node moves along its full force, x_i += dt F(x_i), each node's
     largest atomic displacement capped at ``max_step_A``;
  2. reparametrize: nodes re-placed at equal Cartesian arc length on the piecewise-linear
     string, after each node is rigidly aligned onto its predecessor (rigid motion does
     not count as arc length).
No springs, so nodes cannot all slide into one basin (the CI-NEB failure of 2026-11-19).
Converged when the largest atomic displacement of any node in an iteration, after rigidly
aligning the node onto its previous position, is <= ``string_tol_A``.

Climbing image: the highest interior node climbs along the tangent t = (x_{i+1} - x_{i-1})
/ |...| from its fixed neighbours, with F_ci = F - 2 (F.t) t, by FIRE (afir_batch._Fire),
until the full unbiased fmax <= ``ci_fmax`` (the HF kernel's 0.01 eV/A). Stops:
``above_ceiling`` when its energy exceeds ``ceiling_eV`` (the runner's TS ceiling, start
energy + 3 eV), ``ci_max_steps``. E/F only: no Hessian, no curvature estimate.
"""

from __future__ import annotations

import numpy as np

from .afir_batch import _evaluate, _Fire
from .ts_ef import _align


def _aligned_chain(x):
    out = [x[0]]
    for k in range(1, len(x)):
        out.append(_align(out[-1], x[k]))
    return np.stack(out)


def reparametrize(x):
    """Equal arc length along the piecewise-linear string, ends kept; nodes aligned onto
    their predecessors first."""
    y = _aligned_chain(x)
    seg = np.linalg.norm((y[1:] - y[:-1]).reshape(len(y) - 1, -1), axis=1)
    s = np.concatenate([[0.0], np.cumsum(seg)])
    targets = np.linspace(0.0, s[-1], len(y))
    out = [y[0]]
    for t in targets[1:-1]:
        k = min(int(np.searchsorted(s, t, side="right")) - 1, len(y) - 2)
        w = 0.0 if s[k + 1] == s[k] else (t - s[k]) / (s[k + 1] - s[k])
        out.append((1 - w) * y[k] + w * y[k + 1])
    out.append(y[-1])
    return np.stack(out)


def resample(frames, n):
    """n structures at equal arc length along a sequence of structures (ends kept; each
    structure aligned onto its predecessor first)."""
    y = _aligned_chain(np.stack([np.asarray(f, float) for f in frames]))
    seg = np.linalg.norm((y[1:] - y[:-1]).reshape(len(y) - 1, -1), axis=1)
    s = np.concatenate([[0.0], np.cumsum(seg)])
    out = []
    for t in np.linspace(0.0, s[-1], n):
        k = min(int(np.searchsorted(s, t, side="right")) - 1, len(y) - 2)
        w = 0.0 if s[k + 1] == s[k] else (t - s[k]) / (s[k + 1] - s[k])
        out.append((1 - w) * y[k] + w * y[k + 1])
    return np.stack(out)


def tangent(x, i):
    """Unit tangent at node i from its two neighbours (one-sided at an end). The nodes must
    share one frame (a reparametrized string, or consecutive frames of one trajectory):
    re-aligning the neighbours here mixed rotation into tau and stalled the climb on the
    TwoSaddle bench (2026-11-19)."""
    lo, hi = max(i - 1, 0), min(i + 1, len(x) - 1)
    t = np.asarray(x[hi], float) - np.asarray(x[lo], float)
    return t / np.linalg.norm(t)


def climb(template, x0, tau, evaluate, ledger, task, *, max_step_A=0.05, ci_fmax=0.01,
          ci_steps=2000, ceiling_eV=None):  # fmt: skip
    """Climbing image from x0 along the FIXED unit direction tau: F_ci = F - 2 (F.tau) tau,
    FIRE, until the full unbiased fmax <= ci_fmax. Returns ci_status (converged /
    above_ceiling / max_steps), ci_steps, ci_positions, ci_energy_eV, ci_fmax_eV_A,
    ci_trace (step, E, fmax every 10 steps), ci_forwards."""
    xc = np.asarray(x0, float).copy()
    tau = np.asarray(tau, float) / np.linalg.norm(tau)
    fire = _Fire(1, xc.shape, max_step_A)
    status, trace, spent0 = "max_steps", [], ledger.spent
    with ledger.segment(task=task, stage="climb", attempt=0) as seg:
        for k in range(ci_steps + 1):
            ec, fc = _evaluate(evaluate, template, xc[None], ledger)
            ec, fc = float(ec[0]), fc[0]
            fmax = float(np.linalg.norm(fc, axis=1).max())
            if k % 10 == 0:
                trace.append((k, ec, fmax))
            if ceiling_eV is not None and ec > ceiling_eV:
                status = "above_ceiling"
                break
            if fmax <= ci_fmax:
                status = "converged"
                break
            if k == ci_steps:
                break
            fci = fc - 2.0 * np.vdot(fc, tau) * tau
            xc = xc + fire.step([0], fci[None])[0]
        seg.outcome = status
    return {"ci_status": status, "ci_steps": k, "ci_positions": xc, "ci_energy_eV": ec,
            "ci_fmax_eV_A": fmax, "ci_trace": trace, "ci_forwards": ledger.spent - spent0}  # fmt: skip


def string_climb(template, nodes, evaluate, ledger, task, *, dt=0.05, max_step_A=0.05,
                 string_tol_A=2e-3, string_steps=500, ci_fmax=0.01, ci_steps=2000,
                 ceiling_eV=None):  # fmt: skip
    """Simplified string from ``nodes`` (end points = two minima, kept fixed), then climbing
    image on the highest interior node. ``evaluate`` is a batched evaluator (one call per
    string iteration for all interior nodes); every structure charges ``ledger`` once.
    Returns a dict: string status / iterations, node energies, climbing-image status,
    steps, positions, energy, fmax, forwards."""
    x = np.stack([np.asarray(n, float) for n in nodes])
    if len(x) < 4:
        raise ValueError("need at least two interior nodes")
    x = reparametrize(x)
    spent0 = ledger.spent
    e_end, _ = _evaluate(evaluate, template, x[[0, -1]], ledger)
    out = {"string_status": "max_steps", "string_iterations": 0}
    with ledger.segment(task=task, stage="string", attempt=0) as seg:
        for it in range(1, string_steps + 1):
            e, f = _evaluate(evaluate, template, x[1:-1], ledger)
            step = dt * f
            size = np.linalg.norm(step, axis=2).max(axis=1)
            step *= np.minimum(1.0, max_step_A / np.where(size > 0, size, 1.0))[:, None, None]
            moved = x.copy()
            moved[1:-1] += step
            new = reparametrize(moved)
            # rigid motion is not a change (the end points re-align every iteration)
            change = max(
                float(np.linalg.norm(new[i] - _align(new[i], x[i]), axis=1).max())
                for i in range(len(x))
            )
            x = new
            out["string_iterations"] = it
            if change <= string_tol_A:
                out["string_status"] = "converged"
                break
        seg.outcome = out["string_status"]
    e, f = _evaluate(evaluate, template, x[1:-1], ledger)
    energies = [float(e_end[0])] + [float(v) for v in e] + [float(e_end[1])]
    out["node_energies_eV"] = energies
    top = int(np.argmax(energies[1:-1])) + 1
    out["ci_index"] = top
    ci = climb(template, x[top], tangent(x, top), evaluate, ledger, task, max_step_A=max_step_A,
               ci_fmax=ci_fmax, ci_steps=ci_steps, ceiling_eV=ceiling_eV)  # fmt: skip
    out.update(ci, string_nodes=x, forwards=ledger.spent - spent0)
    return out
