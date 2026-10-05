"""G5 exploration (NOT acceptance): the four arms on the relay testbed.

PLAN_2026-09-28.md G5 compares discover / +I-a / +I-b / +I-a+I-b under the same seed and
budget, judged by criteria frozen in T2. T2 is not frozen, so this script only MEASURES;
nothing it prints is a criterion or a pass. It exists so the first look at whether the
interventions move anything is a recorded run rather than an impression.

Metric per run: the trial that first admitted a product-compatible node (None = not
reached inside the budget), and how many chemical nodes were found.

Seeds: at temperature_K = 0 with spectral ranking the seed never enters -- the first
batch (runs/connect_g5_explore) gave trial-for-trial identical runs for seeds 1-5 -- so
several seeds are several copies of one run, not replicates.

    python docs/experiments/connect_g5_explore.py --out runs/connect_g5_explore
    python docs/experiments/connect_g5_explore.py --variant detour --seeds 1 \
        --trials 120 --out runs/connect_g5_explore_detour
"""

import argparse
import json
import sys
from dataclasses import replace
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from ase.optimize import LBFGS  # noqa: E402

from prrs.calculators import (  # noqa: E402
    ProtonRelay,
    ProtonRelayGated,
    proton_relay_atoms,
    proton_relay_factory,
    proton_relay_gated_atoms,
    proton_relay_gated_factory,
)
from prrs.config import SearchConfig  # noqa: E402
from prrs.search import ReactionSearch  # noqa: E402

ARMS = {
    "discover": {},
    "+I-a": {"connect_queue_priority": True},
    "+I-b": {"connect_direction_priority": True},
    "+I-a+I-b": {"connect_queue_priority": True, "connect_direction_priority": True},
}


VARIANTS = {
    # name: (factory, calculator class, reactant, product)
    "chain": (
        proton_relay_factory,
        ProtonRelay,
        lambda: proton_relay_atoms((0, 0, 0)),
        lambda: proton_relay_atoms((1, 1, 1)),
    ),
    "detour": (
        proton_relay_gated_factory,
        ProtonRelayGated,
        lambda: proton_relay_gated_atoms((0, 0, 0), True),
        lambda: proton_relay_gated_atoms((1, 1, 1), True),
    ),
}


def product(variant):
    _, calculator, _, build = VARIANTS[variant]
    atoms = build()
    atoms.calc = calculator()
    LBFGS(atoms, logfile=None).run(fmax=1e-4, steps=3000)
    atoms.calc = None
    return atoms


def first_product_trial(network):
    hits = [
        int(node["discovered_by"][1:])
        for node in network["chemical_nodes"]
        if node.get("product_compatible") and node["discovered_by"]
    ]
    return min(hits) + 1 if hits else None  # t000041 is the 42nd trial


def _path_d(network):
    """d along the observed R -> P path, or None: a rise in it is the detour taken."""
    path = network["product_paths"]["observed"]
    if not path["connected"]:
        return None
    d = {
        node["id"]: min(s["d_to_product"] for s in node["microstates"])
        for node in network["chemical_nodes"]
    }
    return [d[path["edges"][0]["from"]]] + [d[edge["to"]] for edge in path["edges"]]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True)
    parser.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3, 4, 5])
    parser.add_argument("--trials", type=int, default=80)
    parser.add_argument("--depth", type=int, default=3)
    parser.add_argument("--variant", choices=sorted(VARIANTS), default="chain")
    args = parser.parse_args()
    out = Path(args.out)
    base = SearchConfig(
        families=("stretch", "compress"),
        closed_shell_only=False,
        max_depth=args.depth,
        max_trials=args.trials,
    )
    factory, _, reactant, _ = VARIANTS[args.variant]
    target = product(args.variant)
    rows = []
    for seed in args.seeds:
        for arm, switches in ARMS.items():
            config = replace(base, seed=seed, **switches)
            run_dir = out / f"seed{seed}" / arm.replace("+", "p").replace("-", "")
            status = None
            if (run_dir / "network.json").exists():
                status = json.loads((run_dir / "network.json").read_text()).get("status")
            if status == "completed":  # finished by an earlier invocation: reuse it
                network = json.loads((run_dir / "network.json").read_text())
            else:  # "running" = killed mid-run: the mainline resume picks it up
                network = ReactionSearch(factory, config).run(
                    reactant(), run_dir, product=target, resume=status == "running"
                )
            row = {
                "variant": args.variant,
                "seed": seed,
                "arm": arm,
                "first_product_trial": first_product_trial(network),
                "nodes": len(network["chemical_nodes"]),
                "min_d_reached": min(
                    s["d_to_product"] for n in network["chemical_nodes"] for s in n["microstates"]
                ),
                "trials": network["attempts"],
                "observed_path_d": _path_d(network),
            }
            rows.append(row)
            print(json.dumps(row), flush=True)
    (out / "summary.json").write_text(json.dumps({"rows": rows, "is_acceptance": False}, indent=2))


if __name__ == "__main__":
    main()
