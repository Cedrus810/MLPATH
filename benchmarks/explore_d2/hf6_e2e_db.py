"""HF-6 phase 2, DB arm: ASE Dimer from the same L0 as a D2 run, up to the same total budget
(docs/PLAN_HF6_DIMER_E2E_2026-11-11.md §2.2; user 2026-11-11 「可以」).

Starts are the run's L0 poses (poses.jsonl) as stored. One attempt: ASE Dimer with
DimerControl defaults (gauss displacement, gauss_std 0.1 A, random initial eigenmode),
`MinModeAtoms.displace()` once, then `MinModeTranslate` until full fmax <= the run's
fmax_eV_A; each attempt on its own ledger capped at min(3000, budget left). Poses are visited
in a frozen random order, round k gives every pose one attempt seeded by (seed, pose, k);
attempts continue until the cumulative forward count reaches BUDGET. No reference, product,
negative mode or D2 scan information is used.

One JSON per attempt (k_pose); existing files count towards the budget and are skipped, so
rerunning resumes.

    python hf6_e2e_db.py OUT_DIR RUN_DIR BUDGET SEED [MAX_ATTEMPTS]
"""

import json
import sys
import time
import traceback
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT / "src"), str(ROOT / "explore" / "src")]

from ase import Atoms  # noqa: E402
from ase.mep.dimer import DimerControl, MinModeAtoms, MinModeTranslate  # noqa: E402

from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.ledger import BudgetExhausted, Ledger, view_factory  # noqa: E402
from explore.streams import stream  # noqa: E402

CAP = 3000


def attempt(pose, factory, seed, fmax):
    atoms = pose.copy()
    atoms.calc = factory()
    d = MinModeAtoms(atoms, DimerControl(logfile=None), random_seed=seed)
    start = atoms.positions.copy()
    converged, stop, opt = False, None, None
    try:
        d.displace()
        opt = MinModeTranslate(d, logfile=None)
        converged = bool(opt.run(fmax=fmax, steps=100000))
    except BudgetExhausted:
        stop = "budget_exhausted"
    steps = 0 if opt is None else int(opt.nsteps)  # read also when the budget ran out
    f = None if stop else atoms.get_forces()
    return {
        "status": "converged" if converged else (stop or "not_converged"),
        "energy_eV": None if stop else float(atoms.get_potential_energy()),
        "fmax_eV_A": None if f is None else float(np.linalg.norm(f, axis=1).max()),
        "curvature": None if d.curvatures is None else float(np.asarray(d.curvatures)[0]),
        "steps": steps,
        "positions": atoms.positions.tolist(),
        "start_displacement_A": float(np.linalg.norm(atoms.positions - start, axis=1).max()),
    }


def main():
    out, run_dir = Path(sys.argv[1]), Path(sys.argv[2])
    budget, seed = int(sys.argv[3]), int(sys.argv[4])
    max_attempts = int(sys.argv[5]) if len(sys.argv) > 5 else 0
    out.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    fmax = manifest["settings"]["fmax_eV_A"]
    poses = []
    for line in (run_dir / "poses.jsonl").read_text().splitlines():
        r = json.loads(line)
        a = Atoms(numbers=r["numbers"], positions=r["positions"])
        a.info.update(r["info"])
        poses.append((r["pose_id"], a))
    order = sorted(
        range(len(poses)), key=lambda i: float(stream(seed, "db", poses[i][0]).random())
    )
    calc = polar1("cuda", "float32")
    spent, n = 0, 0
    k = 0
    while spent < budget and (max_attempts <= 0 or n < max_attempts):
        for i in order:
            if spent >= budget or (max_attempts > 0 and n >= max_attempts):
                break
            pose_id, pose = poses[i]
            path = out / f"{k:03d}_{i:03d}.json"
            n += 1
            if path.exists():
                spent += json.loads(path.read_text())["forwards"]
                continue
            cap = min(CAP, budget - spent)
            ledger = Ledger(path.with_suffix(".ledger.jsonl"), arm="DB", budget=cap)
            s = int(stream(seed, "db", pose_id, str(k)).integers(2**31))
            t0 = time.monotonic()
            try:
                res = attempt(pose, view_factory(calc, ledger), s, fmax)
            except Exception as exc:
                res = {"status": f"error:{type(exc).__name__}", "error": str(exc),
                       "traceback": traceback.format_exc()}  # fmt: skip
            ledger.close()
            spent += ledger.spent
            res |= {
                "round": k,
                "pose_index": i,
                "pose_id": pose_id,
                "seed": s,
                "cap": cap,
                "forwards": ledger.spent,
                "cumulative_forwards": spent,
                "wall_s": time.monotonic() - t0,
            }
            path.write_text(json.dumps(res) + "\n", encoding="utf-8")
            print(
                k, i, res["status"], res["forwards"], spent, round(res["wall_s"], 1), flush=True
            )
        k += 1
    (out / "done.json").write_text(
        json.dumps({"budget": budget, "spent": spent, "attempts": n})
    )


if __name__ == "__main__":
    main()
