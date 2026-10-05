"""b05: continue the three 0.6 A quenches that ran out of steps, and nothing else.

Task 3's follow-up (user ruling, 2026-09-09). `runs/b05_probe3_budget` left three trials
at `quench_nonconverged`: stretch on each of the three symmetry-equivalent beta C-H bonds
at 0.6 A, each stopping at exactly `quench_steps = 500` with a residual force of
2.2-2.6 eV/A against a 0.01 threshold. Those three are the only unknown OUTCOMES this
grid has produced, and they sit on the elimination channel's own coordinate.

What this does, and its limits:

  * restarts from each trial's own quench STARTING point -- the last frame of its
    response trajectory -- not from its truncated endpoint, so the continuation is the
    same optimisation run longer rather than a fresh one from a different geometry
  * same potential, same fmax, same optimiser, same polish. `quench_steps` is the only
    parameter changed, to 2x the original
  * three samples kept separate. They are symmetry-equivalent and their agreement is
    evidence of reproducibility, which is worth nothing if it is assumed rather than
    measured
  * records energy, fmax and displacement per step, so "still descending" is read off a
    trajectory instead of inferred from a non-convergence code

A finding here belongs to THIS continuation, not to the original search: it is not a
channel the probe discovered, and it is not written back into probe4.
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

from ase.io import read, write  # noqa: E402

from prrs.config import SearchConfig  # noqa: E402
from prrs.runner import GuardedCalculator, confirm_minimum, quench  # noqa: E402

RUN = ROOT / "runs" / "b05_probe3_budget"
SAMPLES = {"t000002": (2, 4), "t000058": (10, 12), "t000080": (6, 8)}


def load_config(multiplier):
    text = (ROOT / "benchmarks" / "b05_tertiary_halide" / "config.json").read_text()
    raw = json.loads(
        "\n".join(line for line in text.splitlines() if not line.lstrip().startswith("//"))
    )
    original = raw["quench_steps"]
    raw["quench_steps"] = original * multiplier
    return SearchConfig(**raw), original


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--multiplier", type=int, default=2, help="quench_steps multiplier")
    parser.add_argument("--out", default=str(ROOT / "runs" / "b05_quench_continuation"))
    args = parser.parse_args(argv)

    from benchmarks.b04_sn2_chloride.calculator import polar1

    config, original_steps = load_config(args.multiplier)
    print(
        f"quench_steps {original_steps} -> {config.quench_steps} "
        f"({args.multiplier}x); fmax {config.quench_fmax_eV_A}; everything else unchanged"
    )
    factory = polar1
    physical = factory()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    summary = {}

    for trial, bond in SAMPLES.items():
        frames = read(RUN / "trials" / trial / "response.traj", index=":")
        start = frames[-1].copy()  # the quench's own starting point
        start.info["charge"], start.info["spin"] = config.total_charge, config.multiplicity
        guard = GuardedCalculator(physical, config, None)
        start.calc = guard
        history = []

        def record(step, atoms=start, origin=start.positions.copy()):
            forces = atoms.get_forces()
            history.append(
                {
                    "step": int(step),
                    "energy_eV": float(atoms.get_potential_energy()),
                    "fmax_eV_A": float(np.abs(forces).max()),
                    "displacement_A": float(
                        np.linalg.norm(atoms.positions - origin, axis=1).max()
                    ),
                    "r_bond_A": float(atoms.get_distance(*bond)),
                }
            )

        converged, diagnostics = quench(start, guard, config, observe=record)
        first, last = history[0], history[-1]
        print(f"\n{trial}  stretch{bond}@0.6")
        print(
            f"   converged={converged}  steps={diagnostics.get('steps')} "
            f"reason={diagnostics.get('reason')}"
        )
        print(
            f"   energy {first['energy_eV']:.4f} -> {last['energy_eV']:.4f} eV "
            f"(delta {last['energy_eV'] - first['energy_eV']:+.4f})"
        )
        print(
            f"   fmax   {first['fmax_eV_A']:.4f} -> {last['fmax_eV_A']:.4f} eV/A "
            f"(threshold {config.quench_fmax_eV_A})"
        )
        print(f"   r{bond} {first['r_bond_A']:.4f} -> {last['r_bond_A']:.4f} A")

        identity = None
        if converged:
            ok, checked = confirm_minimum(start, factory, config, seed=0)
            identity = {
                "confirm_minimum": bool(ok),
                "detail": {
                    k: v
                    for k, v in (checked or {}).items()
                    if k in ("saddle_order", "reason", "eigenvalues_eV_A2", "curvature_source")
                },
            }
            print(f"   confirm_minimum: {ok}  {identity['detail']}")
            write(str(out / f"{trial}_endpoint.extxyz"), start)

        (out / f"{trial}_history.json").write_text(json.dumps(history, indent=1))
        summary[trial] = {
            "bond": list(bond),
            "converged": bool(converged),
            "steps": diagnostics.get("steps"),
            "reason": diagnostics.get("reason"),
            "energy_start_eV": first["energy_eV"],
            "energy_end_eV": last["energy_eV"],
            "fmax_start_eV_A": first["fmax_eV_A"],
            "fmax_end_eV_A": last["fmax_eV_A"],
            "identity": identity,
            "quench_steps": config.quench_steps,
            "multiplier": args.multiplier,
            "source_run": RUN.name,
        }

    (out / "summary.json").write_text(json.dumps(summary, indent=1))
    print(f"\nwritten {out}/summary.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
