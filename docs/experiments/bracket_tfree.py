"""Do a run's brackets move with t_free? B0 of PLAN_2026-09-25, for any finished run.

Generalises docs/experiments/p1_bracket_tfree.py (P1, fixed choices) to take the run,
the calculator factory and the lengths on the command line; first use is b04, as
evidence for decision R8 (DECISIONS_PENDING). Measured, not changed.

Every amplitude each scan sampled (coarse and refinement) is re-run from the same source
microstate with the same probe at each response_steps; the run's own length comes first
and must reproduce every published endpoint.extxyz, or nothing else runs (that pins
source, probe, charge and dynamics without replaying search.py's minimum confirmation,
whose seeds depend on registry state at the time). All lengths are judged by one rule --
where the quench lands, read with Registry.identify -- and the windows are the adjacent
(returned, departed) pairs over that fixed amplitude set.

Only windows are compared: a benchmark's check.py consumes a whole network (saddle
searches, ledgers), which a replay of samples does not produce.

    python docs/experiments/bracket_tfree.py --run runs/b04_frozen2 \\
        --calculator benchmarks.b04_sn2_chloride.calculator:polar1 \\
        --lengths 100,400 --out runs/b04_bracket_tfree
"""

import argparse
import importlib
import json
import sys
from dataclasses import replace
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
from ase.io import read  # noqa: E402

from prrs.config import SearchConfig  # noqa: E402
from prrs.network import Registry  # noqa: E402
from prrs.perturbations import Probe  # noqa: E402
from prrs.runner import run_trial  # noqa: E402
from prrs.search import _stamped_with_charge  # noqa: E402

PROBE_FIELDS = ("family", "indices", "sign", "amplitude", "seed", "priority", "priority_reason")


def windows(labelled):
    ordered = sorted(labelled)
    return [
        [low, high]
        for (low, a), (high, b) in zip(ordered, ordered[1:])
        if a == "returned" and b.startswith("departed")
    ]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--calculator", required=True, help="module:factory")
    parser.add_argument("--lengths", required=True, help="extra response_steps, e.g. 100,400")
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--scans", default="all", help="comma-separated scan indices, or all")
    args = parser.parse_args()

    run, out = args.run.resolve(), args.out.resolve()
    if out.exists():
        raise SystemExit(
            f"{out} exists; remove it to re-measure (trial directories are never reused)"
        )
    manifest = json.loads((run / "config.json").read_text(encoding="utf-8-sig"))
    config = SearchConfig.from_dict(manifest["config"])
    network = json.loads((run / "network.json").read_text(encoding="utf-8-sig"))
    registry = Registry.from_network(network, run, config)
    states = {s.id: s for node in registry.nodes for s in node.microstates}
    records = {}
    for line in (run / "attempts.jsonl").open():
        if line.strip():
            record = json.loads(line)
            records[record["attempt_id"]] = record
    module, attr = args.calculator.split(":")
    calculator = getattr(importlib.import_module(module), attr)()
    lengths = [config.response_steps] + [
        int(x) for x in args.lengths.split(",") if int(x) != config.response_steps
    ]
    scans = network["amplitude_scans"]
    indices = (
        range(len(scans)) if args.scans == "all" else [int(x) for x in args.scans.split(",")]
    )

    summary = {
        "source_run": str(run.relative_to(ROOT)),
        "calculator": args.calculator,
        "baseline_response_steps": config.response_steps,
        "temperature_K": config.temperature_K,
        "meaning": "fixed amplitude set, one judging rule, only response_steps differs",
        "scans": [],
    }
    for index in indices:
        scan = scans[index]
        source = states[scan["source_microstate"]]
        # stamps a copy and keeps the calculator run_trial reads the reference energy from
        structure = _stamped_with_charge(source.structure, config)
        entry = {
            "scan": index,
            "coordinate": f"{scan['family']} {tuple(scan['indices'])} {scan['sign']:+d}",
            "source": source.id,
            "published_brackets": [[b["low"], b["high"]] for b in scan["local_brackets"]],
            "by_length": {},
        }
        for steps in lengths:
            cfg = replace(config, response_steps=steps)
            labelled = []
            for sample in scan["samples"]:
                published = records[sample["attempt"]]
                probe = Probe(**{k: published["probe"][k] for k in PROBE_FIELDS})
                trial_id = f"s{index:02d}L{steps:04d}_{sample['attempt']}"
                outcome = run_trial(
                    structure,
                    probe,
                    lambda: calculator,
                    cfg,
                    out / "trials" / trial_id,
                    trial_id,
                )
                if steps == config.response_steps:
                    endpoint = run / "trials" / sample["attempt"] / "endpoint.extxyz"
                    if endpoint.exists() != (outcome.endpoint is not None) or (
                        outcome.endpoint is not None
                        and not np.allclose(
                            read(str(endpoint)).positions, outcome.endpoint.positions, atol=1e-5
                        )
                    ):
                        raise SystemExit(
                            f"scan {index} amplitude {sample['amplitude']}: baseline endpoint "
                            "differs from the published one -- the replay is not the run"
                        )
                if outcome.endpoint is None:
                    label = f"none:{outcome.record['status']}"
                else:
                    _, landed = registry.identify(outcome.endpoint)
                    if landed is None:
                        label = "departed:unidentified"
                    elif landed.id == source.id:
                        label = "returned"
                    else:
                        label = f"departed:{landed.id}"
                labelled.append((sample["amplitude"], label))
            entry["by_length"][str(steps)] = {
                "labels": {str(a): lab for a, lab in sorted(labelled)},
                "windows": windows(labelled),
            }
            print(
                f"scan {index:2d} {entry['coordinate']:22s} L={steps:4d} "
                f"({steps * cfg.timestep_fs:5.1f} fs): "
                + " ".join(f"{a:.4g}:{lab.split(':')[0][0]}" for a, lab in sorted(labelled))
                + f"  windows {windows(labelled)}",
                flush=True,
            )
        summary["scans"].append(entry)
        out.mkdir(parents=True, exist_ok=True)
        (out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    moved = [
        s["scan"]
        for s in summary["scans"]
        if len({json.dumps(v["windows"]) for v in s["by_length"].values()}) > 1
    ]
    print(f"done; scans whose windows move: {moved} (r = returned, d = departed, n = none)")


if __name__ == "__main__":
    main()
