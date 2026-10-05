"""Does the P1 committor depend on how long a shot runs before it is quenched?

REVIEW_2026-09-24 §4½ / decision 7. A shot runs `response_steps` of free NVE, then is
quenched, and the quench decides the basin. A shot that has not committed by then is
judged by steepest descent from wherever it happens to be -- the geometric basin at
t_free, not where the dynamics goes. If p_B moves with t_free, the published p_B mixes
the two, and "stop on first entry into A or B" is not a matter of taste.

Paired by construction. Every shot of runs/p1_committor2 is re-run from the SAME ridge
frame with the SAME Maxwell-Boltzmann seed; only `response_steps` differs. The
baseline length (the run's own, 160) is re-shot too and must reproduce the published
landings shot for shot -- that is the check that frame, seed and judging are the ones
the run used, on today's fingerprint. A frame or seed that does not reproduce aborts
before any dynamics (the momenta are compared directly, no MD needed); a baseline whose
landings then differ is flagged, not fatal -- that would be the dynamics not replaying
bit for bit, and it sets the noise floor the other lengths are read against.

Exploratory measurement, not a protocol: nothing here feeds a criterion.

    python docs/experiments/p1_tfree_sensitivity.py            # all lengths
    python docs/experiments/p1_tfree_sensitivity.py --check    # seeds/frames only
"""

import argparse
import collections
import json
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

from ase.io import read  # noqa: E402

from prrs.calculators import mace_factory  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402
from prrs.network import Registry  # noqa: E402
from prrs.perturbations import Probe  # noqa: E402
from prrs.runner import shoot, thermal_momenta, wilson_interval  # noqa: E402
from prrs.search import _dividing_frame  # noqa: E402

RUN = ROOT / "runs" / "p1_committor2"
OUT = ROOT / "runs" / "p1_tfree_sensitivity"
CALC = ROOT / ".scratch" / "p1_committor_calc.json"
LENGTHS = (160, 40, 80, 320, 640)  # baseline first: it gates the rest


def brackets(network, config):
    """Every bracket that was shot, rebuilt the way search.py built it."""
    for scan_index, scan in enumerate(network["amplitude_scans"]):
        for bracket_index, bracket in enumerate(scan.get("local_brackets") or []):
            record = bracket.get("committor") or {}
            if not record.get("shots"):
                continue
            by_amplitude = {s["amplitude"]: s for s in scan["samples"]}
            at_high = by_amplitude[bracket["high"]]
            first = min(scan["samples"], key=lambda s: s["amplitude"])
            # probes[0] of the scan: the lowest coarse amplitude, the scan's own seed
            probe = Probe(
                family=scan["family"],
                indices=tuple(scan["indices"]),
                sign=scan["sign"],
                amplitude=first["amplitude"],
                seed=scan["seed"],
            )
            # search.py: seed0 = cfg.seed + 1009 * len(records) at the committor, and
            # the scan's trials are the latest records, so len = last trial + 1
            last = max(int(s["attempt"][1:]) for s in scan["samples"])
            yield {
                "scan": scan_index,
                "bracket": bracket_index,
                "window": [bracket["low"], bracket["high"]],
                "coordinate": f"{scan['family']} {tuple(scan['indices'])} {scan['sign']:+d}",
                "A": scan["source_microstate"],
                "B": bracket.get("target_microstate"),
                "published": record,
                "probe": probe,
                "at_high": at_high["attempt"],
                "seed0": config.seed + 1009 * (last + 1),
                "prefix": f"c{scan_index:03d}b{bracket_index:02d}s",
            }


def verify(item, frame, config):
    """Frame and seeds reproduce the published shots, compared without any dynamics."""
    problems = []
    for i in range(item["published"]["shots"]):
        directory = RUN / "trials" / f"{item['prefix']}{i:04d}"
        published = read(str(directory / "response.traj"), index=0)
        if not np.allclose(published.positions, frame.positions, atol=1e-10):
            problems.append(f"shot {i}: start geometry differs from the ridge frame")
            continue
        rows = [json.loads(line) for line in (directory / "observations.jsonl").open()]
        start = next((r for r in rows if r["phase"] == "free_start"), None)
        trial = frame.copy()
        thermal_momenta(trial, config.temperature_K, item["seed0"] + i)
        kinetic = float(trial.get_kinetic_energy())
        if start is None or not np.isclose(start["kinetic_eV"], kinetic, rtol=1e-9):
            problems.append(
                f"shot {i}: kinetic {kinetic:.9f} vs published "
                f"{None if start is None else start['kinetic_eV']}"
            )
    return problems


def judge(registry, outcome, a, b):
    record = outcome.record
    if record["status"] != "completed":
        return f"unidentified:{record['status']}"
    _, landed = registry.identify(outcome.endpoint)
    if landed is None:
        return "unidentified:not_in_registry"
    return "A" if landed.id == a else "B" if landed.id == b else f"other:{landed.id}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true", help="verify seeds/frames only")
    args = parser.parse_args()

    manifest = json.loads((RUN / "config.json").read_text(encoding="utf-8-sig"))
    config = SearchConfig.from_dict(manifest["config"])
    network = json.loads((RUN / "network.json").read_text(encoding="utf-8-sig"))
    registry = Registry.from_network(network, RUN, config)
    items = list(brackets(network, config))
    print(
        f"{len(items)} shot bracket(s) in {RUN.name}; baseline response_steps "
        f"{config.response_steps}",
        flush=True,
    )

    for item in items:
        frame, origin, _, _ = _dividing_frame(
            RUN / "trials" / item["at_high"], probe=item["probe"].to_dict()
        )
        if frame is None:
            raise SystemExit(f"scan {item['scan']}: no ridge frame ({origin})")
        item["frame"], item["frame_origin"] = frame, origin
        problems = verify(item, frame, config)
        if problems:
            raise SystemExit(f"scan {item['scan']} does not reproduce: {problems[:3]}")
        print(
            f"  scan {item['scan']:2d} {item['coordinate']}: frame + "
            f"{item['published']['shots']} seeds reproduce (seed0 {item['seed0']})",
            flush=True,
        )
    if args.check:
        return

    factory_kwargs = json.loads(CALC.read_text())
    calculator = mace_factory(**factory_kwargs)
    OUT.mkdir(parents=True, exist_ok=True)
    summary = {
        "source_run": str(RUN.relative_to(ROOT)),
        "source_fingerprint": manifest.get("implementation_sha256"),
        "meaning": (
            "paired re-shots: same ridge frame, same MB seeds, only response_steps "
            "differs; exploratory, never a criterion (REVIEW_2026-09-24 decision 7)"
        ),
        "lengths": {},
    }
    baseline = {}
    for steps in LENGTHS:
        cfg = replace(config, response_steps=steps)
        per_length = []
        for item in items:
            outcomes = shoot(
                item["frame"],
                item["probe"],
                lambda: calculator,
                cfg,
                OUT / f"L{steps:04d}",
                item["published"]["shots"],
                item["seed0"],
                prefix=item["prefix"],
            )
            verdicts = [judge(registry, o, item["A"], item["B"]) for o in outcomes]
            counts = collections.Counter(v if v in ("A", "B") else "other" for v in verdicts)
            terminated = counts["A"] + counts["B"]
            row = {
                "scan": item["scan"],
                "coordinate": item["coordinate"],
                "window": item["window"],
                "A": counts["A"],
                "B": counts["B"],
                "other": counts["other"],
                "verdicts": verdicts,
                "p_B_conditional": counts["B"] / terminated if terminated else None,
                "wilson_95": list(wilson_interval(counts["B"], terminated))
                if terminated
                else None,
                "max_nve_drift_eV_atom": max(
                    o.record.get("max_nve_drift_eV_atom") or 0.0 for o in outcomes
                ),
                "t_free_fs": steps * cfg.timestep_fs,
            }
            key = item["scan"]
            if steps == config.response_steps:
                published = item["published"]
                # Flagged, not fatal: frame and seeds were already verified exactly,
                # so a mismatch here is the dynamics itself not replaying bit for bit
                # (GPU reductions), which is a finding about the baseline's noise floor
                # -- every flip count below must then be read against it.
                row["baseline_reproduces"] = (counts["A"], counts["B"]) == (
                    published["reached_A"],
                    published["reached_B"],
                )
                if not row["baseline_reproduces"]:
                    print(
                        f"  WARNING scan {key}: baseline A {counts['A']} / B {counts['B']}"
                        f" vs published A {published['reached_A']} / B "
                        f"{published['reached_B']}",
                        flush=True,
                    )
                baseline[key] = verdicts
            else:
                row["flipped_vs_baseline"] = sum(
                    x != y for x, y in zip(verdicts, baseline[key])
                )
            per_length.append(row)
            print(
                f"  L={steps:4d} ({steps * cfg.timestep_fs:5.1f} fs) scan {key:2d}: "
                f"A {counts['A']:2d} / B {counts['B']:2d} / other {counts['other']:2d}"
                + (
                    f"  flipped {row['flipped_vs_baseline']}"
                    if "flipped_vs_baseline" in row
                    else f"  reproduces published: {row['baseline_reproduces']}"
                ),
                flush=True,
            )
        summary["lengths"][str(steps)] = per_length
        (OUT / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print("done", flush=True)


if __name__ == "__main__":
    main()
