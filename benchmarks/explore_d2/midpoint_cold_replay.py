"""Bisection-midpoint relaxations that wander at fixed q: carried BFGS curvature against a
cold start from the same first frame, on one machine (docs/PLAN_HF_MIDPOINT_COLD_2026-11-01.md;
TODO item 9; user 2026-11-01).

Inputs are frozen, read-only D2-HF runs. A midpoint relaxation is a `_relax` call of
`refine_ts_ef` with origin "from_left" / "from_right" at hf.midpoint_transverse_fmax_eV_A;
its first frame already holds the tangent predictor step, so a cold start from that frame
(`_relax(..., state=None)`, same threshold, step cap and ceiling) differs only in the
carried curvature.

  AB  A: for every intent whose recorded HS failed at a midpoint (last such run block in
      argument order), rerun the whole HS with the current code as `_hs_task` does
      (`hs_cap` ledger) and record each midpoint relaxation's first frame and outcome.
      B: cold-start every A midpoint from its first frame, REPS times, on its own ledger.
  CD  C: cold-start the first frame of every recorded failed midpoint, REPS times.
      D: the same for a SEED-sampled N_CONTROL of the recorded converged midpoints.

One JSON per (group, item, midpoint, repeat); an existing file is skipped, so rerunning
the same command resumes. LIMIT > 0 keeps only the first LIMIT items per group (smoke).

    python midpoint_cold_replay.py OUT_DIR MODE REPS N_CONTROL SEED LIMIT RUN_DIR [RUN_DIR ...]
"""

import json
import random
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT / "src"), str(ROOT / "explore" / "src")]
sys.path.insert(0, str(Path(__file__).resolve().parent))

from ase import Atoms  # noqa: E402

import explore.ts_branch as ts_branch  # noqa: E402
import explore.ts_ef as ts_ef  # noqa: E402
from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.d2 import (  # noqa: E402
    D2Settings,
    HFSettings,
    _classify,
    _ef_settings,
    _hf_bracket,
    _intent,
    _points,
    _task,
)
from explore.ledger import BudgetExhausted, Ledger, view_factory  # noqa: E402
from hs_transverse_diag import MIDPOINT_ORIGINS, segments  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402
from prrs.reliability import GateRejected  # noqa: E402

RELAX = ts_ef._relax
COLD_BUDGET = 400  # above the 200-step cap: one forward per BFGS step plus the first


def _lines(path):
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.strip():
                yield json.loads(line)


def load(run_dir):
    m = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    raw = dict(m["settings"])
    raw["magnitudes_A"] = tuple(raw["magnitudes_A"])
    poses = {}
    for r in _lines(run_dir / "poses.jsonl"):
        a = Atoms(numbers=r["numbers"], positions=r["positions"])
        a.info.update(r["info"])
        poses[r["pose_id"]] = a
    hf = HFSettings(**m["hf"])
    settled, failed, converged = {}, [], []
    for ev in _lines(run_dir / "journal.jsonl"):
        if ev.get("event") != "settle" or ev["task"][:3] not in ("L1|", "L2|", "HS|"):
            continue
        if ev["task"].startswith("HS|"):
            rsegs = segments((ev["payload"].get("refine") or {}).get("history") or [])
            mids = [s for s in rsegs if s[0]["origin"] in MIDPOINT_ORIGINS]
            for j, s in enumerate(mids):
                frame = {
                    "run": run_dir.name,
                    "task": ev["task"],
                    "midpoint": j,
                    "origin": s[0]["origin"],
                    "q": s[0]["q_target"],
                    "positions": s[0]["positions"],
                    "recorded_steps": s[-1]["relax_step"],
                    "recorded_rf": s[-1]["reaction_force"],
                }
                if s[-1]["fmax_projected"] <= hf.midpoint_transverse_fmax_eV_A:
                    converged.append(frame)
                elif s is rsegs[-1] and ev["status"] == "unresolved:transverse_not_converged":
                    failed.append(frame)
            ev["payload"] = {"status": ev["status"]}
        settled[ev["task"]] = ev
    rows = {r["intent_id"]: r for r in _lines(run_dir / "population.jsonl")}
    return SimpleNamespace(
        name=run_dir.name,
        settings=D2Settings(**raw),
        hf=hf,
        config=SearchConfig.from_dict(m["config"]),
        poses=poses,
        rows=rows,
        settled=settled,
        failed=failed,
        converged=converged,
    )


def context(run, task):
    row = run.rows[task[3:]]
    points = _points(run.settled[_task("L1", row)]["payload"])
    return row, points, points[0].energy_eV + run.settings.certify_max_rise_eV


def summary(history, n0):
    frames = history[n0:]
    pos = np.array([f["positions"] for f in frames])
    return {
        "frames": len(frames),
        "steps": frames[-1]["relax_step"] if frames else None,
        "fproj": [f["fmax_projected"] for f in frames],
        "fmax": [f["fmax_eV_A"] for f in frames],
        "energy_eV": [f["energy_eV"] for f in frames],
        "reaction_force": [f["reaction_force"] for f in frames],
        "step_max_atom_move_A": np.linalg.norm(np.diff(pos, axis=0), axis=2)
        .max(axis=1)
        .tolist()
        if len(frames) > 1
        else [],
        "last_positions": pos[-1].tolist() if frames else None,
    }


def cold(run, frame, calc, out, settings=None):
    row, _, ceiling = context(run, frame["task"])
    ledger = Ledger(out.with_suffix(".ledger.jsonl"), arm="HF", budget=COLD_BUDGET)
    history, stop, t0 = [], None, time.monotonic()
    try:
        with ledger.segment(task=frame["task"], stage="midpoint_cold") as seg:
            try:
                RELAX(
                    run.poses[row["pose_id"]].copy(),
                    _intent(row),
                    np.asarray(frame["positions"], dtype=float),
                    frame["q"],
                    frame["origin"],
                    view_factory(calc, ledger),
                    run.config,
                    settings or _ef_settings(run.hf),
                    history,
                    ceiling,
                    run.hf.midpoint_transverse_fmax_eV_A,
                    state=None,
                )
            except ts_ef._Stop as exc:
                stop = str(exc)
            except BudgetExhausted:
                stop = "budget_exhausted"
            except GateRejected as exc:
                stop = f"gate_rejected:{exc.code}"
            seg.outcome = stop or "converged"
    finally:
        ledger.close()
    record = {k: frame[k] for k in ("run", "task", "midpoint", "origin", "q", "recorded_steps")}
    record |= {
        "converged": stop is None,
        "stop": stop,
        "spent": ledger.spent,
        "wall_s": time.monotonic() - t0,
        **summary(history, 0),
    }
    out.write_text(json.dumps(record) + "\n", encoding="utf-8")
    return record


def rerun(run, task, calc, out, settings=None):
    """A: the whole HS as `_hs_task` runs it, recording every midpoint relaxation."""
    row, points, ceiling = context(run, task)
    view = SimpleNamespace(settled=run.settled, settings=run.settings)
    classified = _classify(view, row)
    ledger = Ledger(out.with_suffix(".ledger.jsonl"), arm="HF", budget=run.hf.hs_cap)
    mids = []

    def watched(*args, **kwargs):
        origin, history = args[4], args[8]
        if origin not in MIDPOINT_ORIGINS:
            return RELAX(*args, **kwargs)
        n0, stop = len(history), None
        try:
            return RELAX(*args, **kwargs)
        except ts_ef._Stop as exc:
            stop = str(exc)
            raise
        finally:
            first = history[n0] if len(history) > n0 else None
            mids.append(
                {
                    "origin": origin,
                    "q": args[3],
                    "positions": None if first is None else first["positions"].tolist(),
                    "converged": stop is None and len(history) > n0,
                    "stop": stop,
                    **summary(history, n0),
                }
            )

    ts_ef._relax = watched
    t0 = time.monotonic()
    try:
        result = ts_branch.continue_branches_ef(
            run.poses[row["pose_id"]].copy(),
            _intent(row),
            _hf_bracket(classified, points),
            view_factory(calc, ledger),
            run.config,
            ledger,
            task,
            branch_settings=ts_branch.BranchSettings(geometry_noise_A=run.hf.geometry_noise_A),
            settings=settings or _ef_settings(run.hf),
            ceiling_eV=ceiling,
        )
    finally:
        ts_ef._relax = RELAX
        ledger.close()
    record = {
        "run": run.name,
        "task": task,
        "class": classified["class"],
        "status": result.status,
        "branch_status": result.branch_status,
        "spent": result.spent,
        "wall_s": time.monotonic() - t0,
        "recorded_status": run.settled[task]["payload"]["status"],
        "midpoints": mids,
    }
    out.write_text(json.dumps(record) + "\n", encoding="utf-8")
    return record


def log(group, rec, **extra):
    keys = ("task", "status", "converged", "stop", "steps", "spent", "wall_s")
    print(
        json.dumps({"group": group} | {k: rec[k] for k in keys if k in rec} | extra),
        flush=True,
    )


def main():
    out_dir, mode = Path(sys.argv[1]), sys.argv[2]
    reps, n_control, seed, limit = map(int, sys.argv[3:7])
    runs = [load(Path(p)) for p in sys.argv[7:]]
    by_name = {r.name: r for r in runs}
    out_dir.mkdir(parents=True, exist_ok=True)
    calc = polar1("cuda", "float32")
    cut = (lambda xs: xs[:limit]) if limit > 0 else (lambda xs: xs)
    if mode == "AB":
        latest = {}
        for run in runs:
            for f in run.failed:
                latest[f["task"]] = run.name
        items = cut(sorted(latest.items(), key=lambda kv: kv[0]))
        (out_dir / "A_items.json").write_text(json.dumps(items, indent=1) + "\n")
        for i, (task, name) in enumerate(items):
            out = out_dir / f"A_{i:02d}.json"
            if out.exists():
                rec = json.loads(out.read_text(encoding="utf-8"))
            else:
                rec = rerun(by_name[name], task, calc, out)
                log("A", rec, mids=[m["converged"] for m in rec["midpoints"]])
            for j, m in enumerate(rec["midpoints"]):
                if m["positions"] is None:
                    continue
                frame = {"run": name, "task": task, "midpoint": j, "recorded_steps": m["steps"]}
                frame |= {k: m[k] for k in ("origin", "q", "positions")}
                for k in range(reps):
                    o = out_dir / f"B_{i:02d}_m{j}_r{k}.json"
                    if not o.exists():
                        log("B", cold(by_name[name], frame, calc, o), item=i, mid=j, rep=k)
    elif mode == "CD":
        failed = cut([f for r in runs for f in r.failed])
        pool = [f for r in runs for f in r.converged]
        control = cut(random.Random(seed).sample(pool, min(n_control, len(pool))))
        picks = {"C": failed, "D": control}
        (out_dir / "CD_items.json").write_text(
            json.dumps(
                {
                    g: [{k: f[k] for k in ("run", "task", "midpoint", "q")} for f in v]
                    for g, v in picks.items()
                },
                indent=1,
            )
            + "\n"
        )
        for group, frames in picks.items():
            for i, frame in enumerate(frames):
                for k in range(reps):
                    o = out_dir / f"{group}_{i:02d}_r{k}.json"
                    if not o.exists():
                        log(group, cold(by_name[frame["run"]], frame, calc, o), item=i, rep=k)
    else:
        raise SystemExit(f"unknown mode {mode}")


if __name__ == "__main__":
    main()
