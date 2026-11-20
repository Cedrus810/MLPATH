"""Start point and fixed tangent for the climbing image: node vs two-node span vs Hermite
peak (2026-11-19, E1 of docs/TS_LOCATOR_DESIGN_2026-11-19.md; user「可以」to "把优化 1
（Hermite 峰点）和 x_{i±2} 切向加进方向研究脚本，在 N1 这 5 条路径上跑一组").

Same inputs and string as n1_direction_study.py (cached N1 DA paths, end minima, 20 nodes,
simplified string snapshots at k iterations). From every snapshot, three arms, each a
climbing image along a FIXED tangent (F - 2 (F.t) t, FIRE, fmax 0.01, ceiling +3 eV):
  node    : highest interior node i*, t = x[i*+1] - x[i*-1]           (current rule)
  span2   : highest interior node i*, t = x[i*+2] - x[i*-2]           (clamped at the ends)
  hermite : interval (i, i+1) next to i* where the tangential force f = F.t_node changes
            sign uphill-to-downhill (f_i < 0 < f_{i+1}, since dE/ds = -f); cubic
            Hermite in arc length through (E_i, -f_i),
            (E_{i+1}, -f_{i+1}); dE/ds = 0 at t* in (0, 1) by bisection; start x(t*) and
            tangent x'(t*) from the uniform Catmull-Rom spline through x[i-1..i+2]. No
            extra force evaluation. If no such sign change next to i*: falls back to
            "node" and is recorded as hermite_status = "fallback".
All arms of one path climb in one batch (one model call per step); sides of converged
images descend in one batch. Cost per arm = end minima + 18 (k + 1) string evaluations
+ its own climb (+ sides, reported apart), as in n1_direction_study.py.
Writes OUT/chains_<seed>_<label>.npz (snapshot nodes, energies, forces per k),
OUT/<arm>_k<k>/<seed>/wave_0000.json for afir_validate.py, OUT/study.json. Not resumable.

    python n1_tangent_variants.py --out OUT --case SEED:LABEL:NPZ [--case ...]
"""

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import run_mc_afir  # noqa: E402,F401  (sets sys.path for run_d2 / systems)
from run_d2 import make_atoms, make_config  # noqa: E402
from systems import SYSTEMS  # noqa: E402

from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.afir_batch import _evaluate, _Fire, minimize  # noqa: E402
from explore.afir_string import reparametrize, resample, tangent  # noqa: E402
from explore.ledger import Ledger  # noqa: E402
from explore.parallel import MACEBatchEvaluator  # noqa: E402
from explore.poses import _jsonable  # noqa: E402
from explore.ts_ef import _align  # noqa: E402

FORMING = {"C1-C6": (0, 11), "C4-C7": (3, 12), "C1-C7": (0, 12), "C4-C6": (3, 11)}
ARMS = ("node", "span2", "hermite")


def span_tangent(x, i, m):
    lo, hi = max(i - m, 0), min(i + m, len(x) - 1)
    t = x[hi] - x[lo]
    return t / np.linalg.norm(t)


def _catmull_rom(p0, p1, p2, p3, t):
    a = 2 * p1
    b = -p0 + p2
    c = 2 * p0 - 5 * p1 + 4 * p2 - p3
    d = -p0 + 3 * p1 - 3 * p2 + p3
    pos = 0.5 * (a + b * t + c * t**2 + d * t**3)
    der = 0.5 * (b + 2 * c * t + 3 * d * t**2)
    return pos, der


def hermite_start(x, energies, forces):
    """x (N, n, 3) one frame, energies (N,), forces (N, n, 3) (end forces may be zeros).
    Returns (start, tau, info)."""
    n = len(x)
    top = int(np.argmax(energies[1:-1])) + 1
    f = np.array([np.vdot(forces[i], tangent(x, i)) for i in range(n)])
    cand = []
    if f[top] < 0 < f[min(top + 1, n - 1)]:
        cand.append(top)
    if top - 1 >= 0 and f[top - 1] < 0 < f[top]:
        cand.append(top - 1)
    info = {"top": top, "f_tan": [float(v) for v in f[max(top - 2, 0) : top + 3]]}
    if not cand:
        info["hermite_status"] = "fallback"
        return x[top], tangent(x, top), info
    i = cand[0]
    h = float(np.linalg.norm(x[i + 1] - x[i]))
    m0, m1 = -f[i] * h, -f[i + 1] * h  # dE/dt at t = 0, 1 (> 0, < 0; dE/ds = -f)
    e0, e1 = energies[i], energies[i + 1]

    def dedt(t):  # derivative of the cubic Hermite E(t)
        return ((6 * t**2 - 6 * t) * e0 + (3 * t**2 - 4 * t + 1) * m0
                + (-6 * t**2 + 6 * t) * e1 + (3 * t**2 - 2 * t) * m1)  # fmt: skip

    lo, hi = 0.0, 1.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if dedt(mid) > 0 else (lo, mid)
    t = 0.5 * (lo + hi)
    p = [x[max(i - 1, 0)], x[i], x[i + 1], x[min(i + 2, n - 1)]]
    pos, der = _catmull_rom(*p, t)
    info.update(hermite_status="ok", interval=[i, i + 1], t_star=t)
    return pos, der / np.linalg.norm(der), info


def climb_batch(template, starts, taus, evaluate, ledger, task, *, max_step_A=0.05,
                ci_fmax=0.01, ci_steps=2000, ceiling_eV=None):  # fmt: skip
    """afir_string.climb for several (start, fixed tau) at once: one batched call per step
    for the still-running rows; each row stops on converged / above_ceiling / max_steps."""
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
                cap = ceiling_eV[r] if np.ndim(ceiling_eV) else ceiling_eV  # per row or scalar
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--case", action="append", required=True, help="SEEDTAG:LABEL:NPZ")
    ap.add_argument("--checkpoints", default="0,2,5,10")
    ap.add_argument("--nodes", type=int, default=20)
    ap.add_argument("--fmax", type=float, default=0.01)
    ap.add_argument("--ci-steps", type=int, default=2000)
    ap.add_argument("--side-step", type=float, default=0.05, help="A, along tau")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    ks = sorted(int(k) for k in args.checkpoints.split(","))
    t0 = time.monotonic()
    system = SYSTEMS["n1"]
    config = make_config(system)
    calc = polar1("cuda", "float32")
    atoms, _ = make_atoms(system, config)
    batch = MACEBatchEvaluator(calc)
    ledger = Ledger(out / "study.ledger.jsonl", arm="TAN", replay=False)
    report, waves = [], {}
    for case in args.case:
        tag, label, npz = case.split(":", 2)
        z = np.load(npz)
        pos, en = z["positions"], z["energies"]
        t1 = time.monotonic()
        s0 = ledger.spent
        ends = minimize(atoms, [pos[0], pos[-1]], batch, config, ledger, f"{label}|ends",
                        fmax=args.fmax, max_steps=2000)  # fmt: skip
        e0 = ends[0].energy_eV
        g_r, g_p = ends[0].edges, ends[1].edges
        ends_cost = ledger.spent - s0
        low = int(np.argmin(en[: len(en) // 3 + 1]))
        mids = resample(list(pos[low:]), args.nodes)[1:-1]
        x = reparametrize(np.stack([ends[0].positions] + list(mids) + [ends[1].positions]))
        snaps, string_cost = {}, 0
        with ledger.segment(task=f"{label}|string", stage="string", attempt=0) as seg:
            for it in range(0, ks[-1] + 1):
                e, f = _evaluate(batch, atoms, x[1:-1], ledger)
                string_cost += len(x) - 2
                if it in ks:
                    zero = np.zeros_like(f[:1])
                    snaps[it] = (x.copy(),
                                 np.array([e0] + [float(v) for v in e] + [ends[1].energy_eV]),
                                 np.concatenate([zero, f, zero]), string_cost)  # fmt: skip
                if it == ks[-1]:
                    break
                step = 0.05 * f
                size = np.linalg.norm(step, axis=2).max(axis=1)
                step *= np.minimum(1.0, 0.05 / np.where(size > 0, size, 1.0))[:, None, None]
                moved = x.copy()
                moved[1:-1] += step
                x = reparametrize(moved)
            seg.outcome = f"{ks[-1]} iterations"
        safe = label.replace("|", "_").replace("+", "_")
        np.savez(out / f"chains_{tag}_{safe}.npz", ks=np.array(ks),
                 nodes=np.stack([snaps[k][0] for k in ks]),
                 energies=np.stack([snaps[k][1] for k in ks]),
                 forces=np.stack([snaps[k][2] for k in ks]))  # fmt: skip
        arms = []
        for k in ks:
            xs, energies, forces, scost = snaps[k]
            top = int(np.argmax(energies[1:-1])) + 1
            hpos, htau, hinfo = hermite_start(xs, energies, forces)
            t_node = tangent(xs, top)
            for arm, start, tau in (("node", xs[top], t_node),
                                    ("span2", xs[top], span_tangent(xs, top, 2)),
                                    ("hermite", hpos, htau)):  # fmt: skip
                arms.append({"k": k, "arm": arm, "start": start, "tau": tau, "scost": scost,
                             "top": top, "hinfo": hinfo if arm == "hermite" else None,
                             "cos_to_node": float(abs(np.vdot(tau, t_node)))})  # fmt: skip
        cis = climb_batch(atoms, [a["start"] for a in arms], [a["tau"] for a in arms], batch,
                          ledger, f"{label}|climb", ci_fmax=args.fmax, ci_steps=args.ci_steps,
                          ceiling_eV=e0 + 3.0)  # fmt: skip
        conv = [j for j, c in enumerate(cis) if c["ci_status"] == "converged"]
        side_starts = []
        for j in conv:
            xc, tau = cis[j]["ci_positions"], arms[j]["tau"]
            side_starts += [xc - args.side_step * tau, xc + args.side_step * tau]
        c1 = ledger.spent
        sides = (minimize(atoms, side_starts, batch, config, ledger, f"{label}|sides",
                          fmax=config.quench_fmax_eV_A, max_steps=2000)
                 if side_starts else [])  # fmt: skip
        sides_cost = ledger.spent - c1
        for j, (a, ci) in enumerate(zip(arms, cis)):
            xc = ci["ci_positions"]
            rec = {"seed": tag, "label": label, "k": a["k"], "arm": a["arm"], "top": a["top"],
                   "cos_tau_to_node": a["cos_to_node"], "hermite": a["hinfo"],
                   "ci_status": ci["ci_status"], "ci_steps": ci["ci_steps"],
                   "ci_rel_eV": ci["ci_energy_eV"] - e0, "ci_fmax_eV_A": ci["ci_fmax_eV_A"],
                   "start_to_ci_rmsd_A": float(np.sqrt(np.mean(np.sum(
                       (xc - _align(xc, a["start"])) ** 2, axis=1)))),
                   "ci_forming_A": {n: round(float(np.linalg.norm(xc[jj] - xc[ii])), 3)
                                    for n, (ii, jj) in FORMING.items()},
                   "cost": {"ends": ends_cost, "string": a["scost"], "climb": ci["ci_forwards"]}}  # fmt: skip
            rec["cost"]["total_locate"] = ends_cost + a["scost"] + ci["ci_forwards"]
            if j in conv:
                q = conv.index(j)
                labels = []
                for m in sides[2 * q : 2 * q + 2]:
                    if m.edges == g_r and abs(m.energy_eV - e0) <= 0.05:
                        labels.append("reactant")
                    elif m.edges == g_p and abs(m.energy_eV - ends[1].energy_eV) <= 0.05:
                        labels.append("product")
                    elif m.edges in (g_r, g_p):
                        labels.append("reactant_graph" if m.edges == g_r else "product_graph")
                    else:
                        labels.append("other")
                rec["sides"] = labels
                key = f"{a['arm']}_k{a['k']}"
                waves.setdefault(key, {}).setdefault(tag, []).append({
                    "k": len(waves.get(key, {}).get(tag, [])), "label": label,
                    "new_events": [],
                    "ts": [{"status": "converged_candidate", "positions": xc,
                            "energy_eV": ci["ci_energy_eV"], "link": {}}],
                })  # fmt: skip
            report.append(rec)
            print(json.dumps({q: v for q, v in rec.items()}, default=str), flush=True)
        print(json.dumps({"case": label, "sides_forwards": sides_cost,
                          "wall_s": time.monotonic() - t1}), flush=True)  # fmt: skip
    ledger.close()
    for key, by_tag in waves.items():
        for tag, paths in by_tag.items():
            d = out / key / tag
            d.mkdir(parents=True, exist_ok=True)
            (d / "wave_0000.json").write_text(json.dumps(_jsonable(
                {"wave": 0, "paths": paths,
                 "spent_total": sum(r["cost"]["total_locate"] for r in report
                                    if f"{r['arm']}_k{r['k']}" == key and r["seed"] == tag)}),
                default=str) + "\n")  # fmt: skip
    (out / "study.json").write_text(json.dumps(_jsonable(
        {"args": vars(args), "wall_s": time.monotonic() - t0, "rows": report}), indent=1,
        default=str) + "\n")  # fmt: skip
    print(json.dumps({"wall_s": time.monotonic() - t0, "forwards": ledger.spent}), flush=True)


if __name__ == "__main__":
    main()
