"""Re-judge the shots of runs/p1_committor1 with the fixed attribution (2026-09-22).

That run reported A 0 / B 0 / other 16 on every bracket because judging read a key
only the scan path writes (search.py, fixed same day). The shots themselves were fine:
64 endpoints are on disk, quenched and converged. So the question the run was killed
before answering -- do 300 K shots from a ridge frame cross the barrier, or do they all
fall back to the side they started on -- is answerable from those endpoints alone, with
no new dynamics.

This rebuilds the registry from the run's own network (Registry.from_network, PLAN item
6.2) and asks `identify` where each endpoint sits. One calculator, 64 single points.

    python docs/experiments/p1_committor_rejudge.py
"""

import collections
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

from ase.io import read  # noqa: E402

from prrs.config import SearchConfig  # noqa: E402
from prrs.network import Registry  # noqa: E402
from prrs.runner import wilson_interval  # noqa: E402

RUN = ROOT / "runs" / "p1_committor1"
MODEL = "/home/ruigengji/MLP/mace/MACE-OFF24_medium.model"


def main():
    network = json.loads((RUN / "network.json").read_text())
    config = SearchConfig.from_dict(json.loads((RUN / "config.json").read_text())["config"])
    registry = Registry.from_network(network, RUN, config)

    from prrs import torch_guard

    torch_guard.install(trusted=[MODEL])
    from mace.calculators import MACECalculator

    calculator = MACECalculator(
        model_paths=MODEL, model_type="MACE", device="cuda", default_dtype="float64"
    )

    # bracket -> its shot prefix, the same one search.py builds
    brackets = [
        (scan_index, bracket_index, bracket)
        for scan_index, scan in enumerate(network["amplitude_scans"])
        for bracket_index, bracket in enumerate(scan.get("local_brackets") or [])
    ]
    results = []
    print(f"{len(brackets)} bracket(s), source microstates from the run's own scans\n")

    for scan_index, bracket_index, bracket in brackets:
        prefix = f"c{scan_index:03d}b{bracket_index:02d}s"
        directories = sorted((RUN / "trials").glob(f"{prefix}*"))
        if not directories:
            continue
        landings = collections.Counter()
        for directory in directories:
            endpoint = directory / "endpoint.extxyz"
            if not endpoint.exists():
                landings["no_endpoint"] += 1
                continue
            atoms = read(endpoint)
            atoms.calc = calculator
            _, state = registry.identify(atoms)
            landings[state.id if state is not None else "unidentified"] += 1
        target = bracket.get("target_microstate")
        shots = sum(landings.values())
        crossed = landings.get(target, 0) if target else 0
        terminated = sum(count for where, count in landings.items() if "/" in where)
        results.append(
            {
                "bracket": [bracket["low"], bracket["high"]],
                "B": target,
                "shots": shots,
                "landings": dict(landings),
                "reached_B": crossed,
                "terminated": terminated,
                "p_B": crossed / terminated if terminated else None,
                "wilson_95": list(wilson_interval(crossed, terminated)) if terminated else None,
            }
        )
        print(f"bracket [{bracket['low']:.4f}, {bracket['high']:.4f}]  B = {target}")
        print(f"  shots {shots}  landings {dict(landings)}")
        print(f"  p_B = {results[-1]['p_B']}  wilson95 {results[-1]['wilson_95']}\n")

    out = RUN / "rejudge_2026-09-22.json"
    out.write_text(
        json.dumps(
            {
                "why": "the run's own judging read a key nothing sets; these are its "
                "shots re-judged read-only with registry.identify (search.py fix, "
                "2026-09-22). The shots themselves were produced at ec32bc3df5d2bf14 "
                "and are untouched -- re-judging reads, it does not re-run.",
                "shots_implementation_sha256": "ec32bc3df5d2bf14",
                "judged_by_implementation_sha256": "63be7064f9e09f82",
                "ensemble": "maxwell-boltzmann",
                "temperature_K": config.temperature_K,
                "resolution": "N=16 answers only whether p_B is near 0 or 1; "
                "distinguishing 0.5 from 0.25 needs N>=50 (protocol, pre-registered)",
                "brackets": results,
            },
            indent=1,
        )
    )
    print(f"written {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
