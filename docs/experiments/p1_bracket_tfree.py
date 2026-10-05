"""Branch B0 of PLAN_2026-09-25: do the search's own brackets move with t_free?

The search judges "returned" / "departed" the way committor shots used to: a fixed free
segment, then a quench. runs/p1_tfree_sensitivity showed that moves P1's p_B, so the same
question stands for where the brackets land. Measured, not changed.

For scans 0 / 2 / 6 / 8 of runs/p1_committor2, every amplitude the scan sampled (coarse
and refinement) is re-run from the same source microstate with the same probe and seed at
response_steps 80 / 160 / 320. All lengths are judged by one rule -- where the quench
lands, read with Registry.identify -- and the windows are the adjacent (returned,
departed) pairs over that fixed amplitude set. No new refinement: the amplitudes are the
published ones, so a moved window shows up as a flipped amplitude.

Baseline gate: at 160 (the run's own length) every replayed endpoint must match the
published endpoint.extxyz. That pins source, probe, seed and dynamics without replaying
search.py's minimum confirmation, whose seeds depend on registry state at the time.

    python docs/experiments/p1_bracket_tfree.py
"""

import json
import sys
from dataclasses import replace
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np  # noqa: E402
from ase.io import read  # noqa: E402

from p1_tfree_sensitivity import CALC, RUN  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402
from prrs.network import Registry  # noqa: E402
from prrs.perturbations import Probe  # noqa: E402
from prrs.runner import run_trial  # noqa: E402

OUT = ROOT / "runs" / "p1_bracket_tfree"
SCANS = (0, 2, 6, 8)
LENGTHS = (160, 80, 320)  # baseline first: it gates the rest
PROBE_FIELDS = ("family", "indices", "sign", "amplitude", "seed", "priority", "priority_reason")


def windows(labelled):
    """Adjacent (returned, departed) pairs over amplitudes sorted ascending."""
    ordered = sorted(labelled)
    return [
        [low, high]
        for (low, a), (high, b) in zip(ordered, ordered[1:])
        if a == "returned" and b.startswith("departed")
    ]


def main():
    if OUT.exists():
        raise SystemExit(
            f"{OUT} exists; remove it to re-measure (trial directories are never reused)"
        )
    manifest = json.loads((RUN / "config.json").read_text(encoding="utf-8-sig"))
    config = SearchConfig.from_dict(manifest["config"])
    network = json.loads((RUN / "network.json").read_text(encoding="utf-8-sig"))
    registry = Registry.from_network(network, RUN, config)
    states = {s.id: s for node in registry.nodes for s in node.microstates}
    records = {}
    for line in (RUN / "attempts.jsonl").open():
        record = json.loads(line)
        records[record["attempt_id"]] = record

    from prrs.calculators import mace_factory

    calculator = mace_factory(**json.loads(CALC.read_text()))
    summary = {
        "source_run": str(RUN.relative_to(ROOT)),
        "meaning": "fixed amplitude set, one judging rule, only response_steps differs",
        "scans": [],
    }
    for index in SCANS:
        scan = network["amplitude_scans"][index]
        source = states[scan["source_microstate"]]
        entry = {
            "scan": index,
            "coordinate": f"{scan['family']} {tuple(scan['indices'])} {scan['sign']:+d}",
            "source": source.id,
            "published": {str(s["amplitude"]): s["target_microstate"] for s in scan["samples"]},
            "by_length": {},
        }
        for steps in LENGTHS:
            cfg = replace(config, response_steps=steps)
            labelled = []
            for sample in scan["samples"]:
                published = records[sample["attempt"]]
                probe = Probe(**{k: published["probe"][k] for k in PROBE_FIELDS})
                trial_id = f"s{index:02d}L{steps:04d}_{sample['attempt']}"
                outcome = run_trial(
                    source.structure,
                    probe,
                    lambda: calculator,
                    cfg,
                    OUT / "trials" / trial_id,
                    trial_id,
                )
                if steps == config.response_steps and outcome.endpoint is not None:
                    then = read(str(RUN / "trials" / sample["attempt"] / "endpoint.extxyz"))
                    if not np.allclose(then.positions, outcome.endpoint.positions, atol=1e-5):
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
                f"scan {index:2d} L={steps:4d} ({steps * cfg.timestep_fs:5.1f} fs): "
                + " ".join(f"{a:.4g}:{lab.split(':')[0][0]}" for a, lab in sorted(labelled))
                + f"  windows {windows(labelled)}",
                flush=True,
            )
        summary["scans"].append(entry)
        OUT.mkdir(parents=True, exist_ok=True)
        (OUT / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print("done (r = returned, d = departed, n = no endpoint)", flush=True)


if __name__ == "__main__":
    main()
