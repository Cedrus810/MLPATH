"""HF-6 phase 1: local localisation, HF E/F search against ASE Dimer on the same frozen
brackets (docs/PLAN_HF6_DIMER_2026-11-10.md §3; user 2026-11-10 「开工」).

Items: every HS task settled in RUN_DIR (its own frozen scan, classification and bracket,
rebuilt exactly as D2 builds them). Two arms per item, each on its own ledger capped at the
run's hs_cap, through one counted view of a resident float32 POLAR-1:
  hf      the run's own frozen HS: `continue_branches_ef` with the run's HFSettings and,
          where the manifest has one, its HSBracket (as d2._hs_run);
  dimer   ASE Dimer (`ase.mep.dimer`, DimerControl defaults, fixed random seed) from the
          highest-energy scan point between the two bracket ends, initial mode the
          normalised gradient of the item's intent coordinate there (the information the
          HS also sees); MinModeTranslate until full fmax <= the run's fmax_eV_A (0.01).
No reference structure, product or negative mode is used by either arm.

One JSON per (item, arm); an existing file is skipped, so rerunning resumes.

    python hf6_local.py OUT_DIR RUN_DIR LIMIT ARMS     (ARMS: hf,dimer; LIMIT > 0: smoke)
"""

import json
import sys
import time
import traceback
from pathlib import Path
from types import SimpleNamespace

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from midpoint_cold_replay import context, load  # noqa: E402  (sets sys.path)

import explore.ts_branch as ts_branch  # noqa: E402
from ase.mep.dimer import DimerControl, MinModeAtoms, MinModeTranslate  # noqa: E402
from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.d2 import HSBracket, _classify, _ef_settings, _hf_bracket  # noqa: E402
from explore.d2 import _intent as intent_of  # noqa: E402
from explore.ledger import BudgetExhausted, Ledger, view_factory  # noqa: E402

SEED = 20261110


def items(run):
    return [t for t, ev in run.settled.items() if t.startswith("HS|")]


def frozen(run, task):
    row, points, ceiling = context(run, task)
    view = SimpleNamespace(settled=run.settled, settings=run.settings)
    classified = _classify(view, row)
    bracket = _hf_bracket(classified, points)
    ends = [next(k for k, p in enumerate(points) if p is b) for b in bracket]
    lo, hi = min(ends), max(ends)
    span = [k for k in range(lo, hi + 1) if points[k].converged] or [lo, hi]
    start = max(span, key=lambda k: points[k].energy_eV)
    return row, classified["class"], bracket, points[start], ceiling


def hf_arm(run, hs_bracket, task, row, bracket, ceiling, factory, ledger):
    result = ts_branch.continue_branches_ef(
        run.poses[row["pose_id"]].copy(),
        intent_of(row),
        bracket,
        factory,
        run.config,
        ledger,
        task,
        branch_settings=ts_branch.BranchSettings(geometry_noise_A=run.hf.geometry_noise_A),
        settings=_ef_settings(run.hf, hs_bracket),
        ceiling_eV=ceiling,
    )
    refine = result.refine or {}
    pos = refine.get("positions")
    return {
        "status": result.status,
        "algorithm": (refine.get("provenance") or {}).get("algorithm"),
        "energy_eV": refine.get("energy_eV"),
        "fmax_eV_A": refine.get("fmax_eV_A"),
        "positions": None if pos is None else np.asarray(pos).tolist(),
    }


def dimer_arm(run, row, start, factory, fmax):
    atoms = run.poses[row["pose_id"]].copy()
    atoms.positions = np.asarray(start.positions, float)
    atoms.calc = factory()
    g = np.asarray(intent_of(row).gradient(atoms.positions), float)
    mode = g / np.linalg.norm(g)
    control = DimerControl(logfile=None)
    d = MinModeAtoms(atoms, control, eigenmodes=[mode], random_seed=SEED)
    opt = MinModeTranslate(d, logfile=None)
    converged, stop = False, None
    try:
        converged = bool(opt.run(fmax=fmax, steps=100000))
    except BudgetExhausted:
        stop = "budget_exhausted"
    f = atoms.get_forces() if stop is None else None
    return {
        "status": "converged" if converged else (stop or "not_converged"),
        "energy_eV": None if stop else float(atoms.get_potential_energy()),
        "fmax_eV_A": None if f is None else float(np.linalg.norm(f, axis=1).max()),
        "curvature": None if d.curvatures is None else float(np.asarray(d.curvatures)[0]),
        "steps": int(opt.nsteps),
        "positions": atoms.positions.tolist(),
        "start_displacement_A": float(
            np.linalg.norm(atoms.positions - np.asarray(start.positions), axis=1).max()
        ),
    }


def main():
    out, run_dir = Path(sys.argv[1]), Path(sys.argv[2])
    limit, arms = int(sys.argv[3]), sys.argv[4].split(",")
    out.mkdir(parents=True, exist_ok=True)
    run = load(run_dir)
    manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    hs_bracket = HSBracket(**manifest["hs_bracket"]) if manifest.get("hs_bracket") else None
    tasks = items(run)
    tasks = tasks[:limit] if limit > 0 else tasks
    (out / "items.json").write_text(json.dumps(tasks, indent=1), encoding="utf-8")
    calc = polar1("cuda", "float32")
    for i, task in enumerate(tasks):
        row, cls, bracket, start, ceiling = frozen(run, task)
        for arm in arms:
            path = out / f"{i:03d}_{arm}.json"
            if path.exists():
                continue
            ledger = Ledger(path.with_suffix(".ledger.jsonl"), arm=arm, budget=run.hf.hs_cap)
            factory = view_factory(calc, ledger)
            t0 = time.monotonic()
            try:
                if arm == "hf":
                    res = hf_arm(run, hs_bracket, task, row, bracket, ceiling, factory, ledger)
                else:
                    res = dimer_arm(run, row, start, factory, run.settings.fmax_eV_A)
            except Exception as exc:
                res = {"status": f"error:{type(exc).__name__}", "error": str(exc),
                       "traceback": traceback.format_exc()}  # fmt: skip
            ledger.close()
            res |= {
                "run": run_dir.name,
                "task": task,
                "class": cls,
                "recorded_hs_status": run.settled[task]["status"],
                "start_q": start.q_target,
                "start_energy_eV": start.energy_eV,
                "forwards": ledger.spent,
                "wall_s": time.monotonic() - t0,
            }
            path.write_text(json.dumps(res) + "\n", encoding="utf-8")
            print(i, arm, res["status"], res["forwards"], round(res["wall_s"], 1), flush=True)


if __name__ == "__main__":
    main()
