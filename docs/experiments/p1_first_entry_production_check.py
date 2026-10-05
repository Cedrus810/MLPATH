"""Protocol stage 2's acceptance on real data: early stop == the rule applied afterwards.

docs/P1_COMMITTOR_FIRST_ENTRY_PROTOCOL.md section 3, stage 2: "same seed, the early-stop
version and stage 1's offline application give the same verdict shot for shot".
Re-shoots every shot of runs/p1_first_entry_confirm -- same ridge frames, same seeds,
k = 160, cap 1280 -- through `shoot(first_entry=...)`, the production path, with the
run's own sample_interval (saving frames does not touch the dynamics), and compares
side and entry step against verdict.json. It also compares the confirmation quench,
which now runs from the entry geometry itself.

    python docs/experiments/p1_first_entry_production_check.py
"""

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from p1_first_entry_confirm import K, L_RUN, N, OUT as CONFIRM, SEED_OFFSET, setup  # noqa: E402
from p1_tfree_sensitivity import CALC  # noqa: E402
from dataclasses import replace  # noqa: E402
from prrs.runner import shoot  # noqa: E402
from prrs.state import encode  # noqa: E402

OUT = ROOT / "runs" / "p1_first_entry_production_check"


def main():
    if OUT.exists():
        raise SystemExit(
            f"{OUT} exists; remove it to re-check (trial directories are never reused)"
        )
    verdict = json.loads((CONFIRM / "verdict.json").read_text())
    offline = {w["scan"]: w for w in verdict["windows"]}
    _, config, registry, items = setup()
    from prrs.calculators import mace_factory

    calculator = mace_factory(**json.loads(CALC.read_text()))
    states = {s.id: s for node in registry.nodes for s in node.microstates}
    cfg = replace(config, response_steps=L_RUN)
    report, disagreements = [], 0
    for item in items:
        graphs = {
            side: frozenset(
                encode(
                    states[item[side]].structure, cfg.bond_scale, active=cfg.active_atoms
                ).edges
            )
            for side in ("A", "B")
        }
        shots = offline[item["scan"]]["shots"]
        rows = []
        for i in range(N):
            (outcome,) = shoot(
                item["frame"],
                item["probe"],
                lambda: calculator,
                cfg,
                OUT,
                1,
                item["seed0"] + SEED_OFFSET + i,
                prefix=f"{item['prefix']}{i:04d}x",
                first_entry={"graphs": graphs, "k": K},
            )
            entry = outcome.record["first_entry"]
            landed = None
            if outcome.endpoint is not None:
                _, state = registry.identify(outcome.endpoint)
                landed = state.id if state is not None else "not_in_registry"
            then = shots[i]
            same = (entry["side"], entry["step"]) == (then["side"], then["entry_step"])
            disagreements += not same
            rows.append(
                {
                    "shot": i,
                    "side": entry["side"],
                    "step": entry["step"],
                    "offline_side": then["side"],
                    "offline_step": then["entry_step"],
                    "same": same,
                    "quench_landed": landed,
                    "offline_quench_landed": then.get("quench_landed"),
                }
            )
        same_quench = sum(r["quench_landed"] == r["offline_quench_landed"] for r in rows)
        print(
            f"scan {item['scan']:2d}: verdicts same {sum(r['same'] for r in rows)}/{N}  "
            f"quench landing same {same_quench}/{N}",
            flush=True,
        )
        report.append({"scan": item["scan"], "shots": rows})
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "check.json").write_text(
        json.dumps({"disagreements": disagreements, "windows": report}, indent=1) + "\n"
    )
    print(f"disagreements: {disagreements} (must be 0)", flush=True)
    raise SystemExit(1 if disagreements else 0)


if __name__ == "__main__":
    main()
