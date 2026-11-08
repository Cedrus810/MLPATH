"""Near-root midpoint relaxations that BFGS cannot finish: optimizer comparison (O5;
docs/PLAN_HF_MIDPOINT_O5_2026-11-02.md; user 2026-11-02).

Every arm relaxes the same start frame at the same q with `ts_ef._relax` semantics
(constraint, guard, ceiling, fine threshold of the run manifest, 200-step cap), cold:
  A0  ASE BFGS, maxstep 0.05 A (current);
  A1  ASE BFGS, maxstep 0.01 A;
  A2  ASE FIRE, default parameters, maxstep 0.05 A (swapped in for BFGS inside ts_ef).

Groups:
  N  near-root stops of two-stage midpoints: for every failed intent (last run block) the
     whole HS is rerun with HF-EF-bracket-v4 (coarse 2e-3, bound 0.025), recording each
     fine-stage midpoint relaxation; the first frame of every fine stage that did not
     converge is a start frame. REPS_N repeats per arm.
  F  the first frame of every recorded failed midpoint (journal), 1 repeat per arm.
  S  a SEED sample of N_CONTROL recorded converged midpoints, 1 repeat per arm.

One JSON per (group, frame, arm, repeat); an existing file is skipped, so rerunning the
same command resumes. LIMIT > 0 keeps the first LIMIT items per group (smoke).

    python midpoint_optimizer_compare.py OUT_DIR GROUPS REPS_N N_CONTROL SEED LIMIT RUN_DIR ...

GROUPS is a string of group letters run in that order, e.g. NFS.
"""

import json
import random
import sys
from dataclasses import replace
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT / "src"), str(ROOT / "explore" / "src")]
sys.path.insert(0, str(Path(__file__).resolve().parent))

from ase.optimize import BFGS, FIRE  # noqa: E402

import explore.ts_ef as ts_ef  # noqa: E402
from midpoint_cold_replay import cold, load, rerun  # noqa: E402

from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.d2 import _ef_settings  # noqa: E402

COARSE, BOUND = 2e-3, 0.025
ARMS = ("A0", "A1", "A2")


def fire(atoms, logfile=None, maxstep=None):
    return FIRE(atoms, logfile=logfile, maxstep=maxstep)


def settings_for(run, arm):
    base = _ef_settings(run.hf)
    return replace(base, max_step_A=0.01) if arm == "A1" else base


def relax(run, frame, arm, calc, out):
    ts_ef.BFGS = fire if arm == "A2" else BFGS
    try:
        rec = cold(run, frame, calc, out, settings=settings_for(run, arm))
    finally:
        ts_ef.BFGS = BFGS
    rec["arm"] = arm
    out.write_text(json.dumps(rec) + "\n", encoding="utf-8")
    return rec


def log(group, rec, **extra):
    keys = ("arm", "converged", "stop", "steps", "spent", "wall_s")
    print(
        json.dumps({"group": group} | {k: rec[k] for k in keys if k in rec} | extra), flush=True
    )


def near_root(runs, by_name, calc, out_dir, reps, limit):
    latest = {}
    for run in runs:
        for f in run.failed:
            latest[f["task"]] = run.name
    items = sorted(latest.items())
    items = items[:limit] if limit > 0 else items
    (out_dir / "N_items.json").write_text(json.dumps(items, indent=1) + "\n")
    for i, (task, name) in enumerate(items):
        run = by_name[name]
        rec_path = out_dir / f"N_{i:02d}_hs.json"
        if rec_path.exists():
            hs = json.loads(rec_path.read_text(encoding="utf-8"))
        else:
            v4 = replace(
                _ef_settings(run.hf),
                midpoint_coarse_transverse_fmax_eV_A=COARSE,
                midpoint_sign_bound_eV_A=BOUND,
            )
            hs = rerun(run, task, calc, rec_path, settings=v4)
            log("N_hs", hs, mids=[m["converged"] for m in hs["midpoints"]])
        for j, m in enumerate(hs["midpoints"]):
            if m["converged"] or m["positions"] is None:
                continue
            frame = {"run": name, "task": task, "midpoint": j, "recorded_steps": m["steps"]}
            frame |= {k: m[k] for k in ("origin", "q", "positions")}
            frame["recorded_rf"] = m["reaction_force"][-1] if m["reaction_force"] else None
            for arm in ARMS:
                for k in range(reps):
                    o = out_dir / f"N_{i:02d}_m{j}_{arm}_r{k}.json"
                    if not o.exists():
                        log("N", relax(run, frame, arm, calc, o), item=i, mid=j, rep=k)


def main():
    out_dir, groups = Path(sys.argv[1]), sys.argv[2]
    reps_n, n_control, seed, limit = map(int, sys.argv[3:7])
    runs = [load(Path(p)) for p in sys.argv[7:]]
    by_name = {r.name: r for r in runs}
    out_dir.mkdir(parents=True, exist_ok=True)
    calc = polar1("cuda", "float32")
    cut = (lambda xs: xs[:limit]) if limit > 0 else (lambda xs: xs)
    picks = {
        "F": cut([f for r in runs for f in r.failed]),
        "S": cut(
            random.Random(seed).sample(
                [f for r in runs for f in r.converged],
                n_control,
            )
        ),
    }
    for group in groups:
        if group == "N":
            near_root(runs, by_name, calc, out_dir, reps_n, limit)
            continue
        frames = picks[group]
        (out_dir / f"{group}_items.json").write_text(
            json.dumps(
                [{k: f[k] for k in ("run", "task", "midpoint", "q")} for f in frames], indent=1
            )
            + "\n"
        )
        for i, frame in enumerate(frames):
            for arm in ARMS:
                o = out_dir / f"{group}_{i:02d}_{arm}.json"
                if not o.exists():
                    log(group, relax(by_name[frame["run"]], frame, arm, calc, o), item=i)


if __name__ == "__main__":
    main()
