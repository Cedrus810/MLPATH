"""How noisy is the analytic-Hessian cross-check that decides the Hessian source?

`prrs.runner._analytic_hessian` compares H v (analytic) with a central-difference
Hessian-vector product along `analytic_hessian_check_directions` seeded random directions at
step `minimum_check_step_A`, and falls back to a 6N finite-difference Hessian when the worst
relative deviation exceeds `analytic_hessian_check` (0.01). On the new b05 E2 saddle the
same geometry and the same seeded directions gave 1.086e-2 once and 6.2e-3 once. This
measures, at fixed saddle geometries of a finished run (read-only), with the run's float32
POLAR-1:
  1. force jitter: repeated force evaluations of one geometry (cache bypassed);
  2. analytic Hessian repeatability: two calls, max relative difference;
  3. the cross-check deviation (same formula as the runner) over steps h and repeats, for
     the runner's seed and a second seed.

    python hessian_check_noise.py OUT_JSON RUN_DIR REPEATS SUFFIX [SUFFIX ...]
"""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT / "src"), str(ROOT / "explore" / "src")]

from ase import Atoms  # noqa: E402
from ase.calculators.calculator import all_changes  # noqa: E402

from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402
from prrs.reliability import GuardedCalculator  # noqa: E402
from prrs.runner import _hessian_providers, hessian_vector_product  # noqa: E402

STEPS = (0.0025, 0.005, 0.01, 0.02, 0.04)


def _lines(path):
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.strip():
                yield json.loads(line)


def analytic(calc, atoms):
    mats = []
    for provider in _hessian_providers(calc):
        scale = (
            float(getattr(provider, "energy_units_to_eV", 1.0))
            / float(getattr(provider, "length_units_to_A", 1.0)) ** 2
        )
        mats.append(
            scale * np.asarray(provider.get_hessian(atoms), float).reshape(3 * len(atoms), -1)
        )
    h = np.mean(mats, axis=0)
    return 0.5 * (h + h.T)


def deviation(hessian, probe, guard, direction, step):
    measured = hessian_vector_product(probe, guard, direction.reshape(-1, 3), step).ravel()
    predicted = hessian @ direction
    return float(np.abs(predicted - measured).max() / max(np.abs(measured).max(), 1e-12))


def main():
    out, run_dir, repeats, suffixes = (
        Path(sys.argv[1]),
        Path(sys.argv[2]),
        int(sys.argv[3]),
        sys.argv[4:],
    )
    m = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    config = SearchConfig.from_dict(m["config"])
    poses = {}
    for r in _lines(run_dir / "poses.jsonl"):
        a = Atoms(numbers=r["numbers"], positions=r["positions"])
        a.info.update(r["info"])
        poses[r["pose_id"]] = a
    pose_of = {r["intent_id"]: r["pose_id"] for r in _lines(run_dir / "population.jsonl")}
    found = [
        ev
        for ev in _lines(run_dir / "journal.jsonl")
        if ev.get("event") == "settle"
        and ev["task"].startswith("HC|")
        and ev["status"] == "connected"
        and any(ev["task"].endswith(s) for s in suffixes)
    ]
    calc = polar1("cuda", "float32")
    results = []
    for ev in found:
        atoms = poses[pose_of[ev["task"][3:]]].copy()
        atoms.set_constraint()
        atoms.positions = np.asarray(ev["payload"]["saddle_positions"], float)
        size = 3 * len(atoms)
        # 1. force jitter at one geometry (bypass the calculator cache every time)
        forces = []
        for _ in range(20):
            calc.calculate(atoms.copy(), ["energy", "forces"], all_changes)
            forces.append(np.asarray(calc.results["forces"], float).copy())
        forces = np.array(forces)
        jitter = float(np.abs(forces - forces.mean(axis=0)).max())
        # 2. analytic Hessian repeatability
        atoms.calc = calc
        h1, h2 = analytic(calc, atoms), analytic(calc, atoms)
        h_rep = float(np.abs(h1 - h2).max() / np.abs(h1).max())
        # 3. cross-check deviation over steps, repeats and two seeds
        guard = GuardedCalculator(calc, config, None)
        probe = atoms.copy()
        probe.calc = guard
        grid = {}
        for seed in (config.seed, config.seed + 1):
            rng = np.random.default_rng(seed)
            dirs = []
            for _ in range(config.analytic_hessian_check_directions):
                d = rng.normal(size=size)
                dirs.append(d / np.linalg.norm(d))
            for step in STEPS:
                worst = []
                for _ in range(repeats):
                    worst.append(max(deviation(h1, probe, guard, d, step) for d in dirs))
                grid[f"seed{seed}_h{step}"] = {
                    "worst_per_repeat": worst,
                    "min": min(worst),
                    "median": float(np.median(worst)),
                    "max": max(worst),
                    "fraction_above_threshold": float(
                        np.mean(np.array(worst) > config.analytic_hessian_check)
                    ),
                }
        results.append(
            {
                "task": ev["task"],
                "force_jitter_max_eV_A": jitter,
                "force_jitter_std_eV_A": float(forces.std(axis=0).max()),
                "analytic_hessian_repeat_rel": h_rep,
                "threshold": config.analytic_hessian_check,
                "runner_step_A": config.minimum_check_step_A,
                "cross_check": grid,
            }
        )
        print(
            json.dumps({"task": ev["task"][-28:], "jitter": jitter, "h_rep": h_rep}), flush=True
        )
        for k, v in grid.items():
            print(
                f"   {k:22s} min {v['min']:.2e} med {v['median']:.2e} max {v['max']:.2e} >thr {v['fraction_above_threshold']:.2f}",
                flush=True,
            )
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(results, indent=1) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
