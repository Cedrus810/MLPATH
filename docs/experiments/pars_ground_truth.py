"""PARS ground truth (exploration, NOT acceptance): refine EVERY candidate, not just
the representatives, so a selection rule can be judged against the full answer.

PARS V1 tests one hypothesis (pars/README.md): Pareto geometry + mechanism-aware
covering cuts refinement cost without destroying reaction diversity. Judging that
needs the refinement outcome of the paths the selection pruned, too. This script
refines all reactive candidates of a run on the run's own MLP, with the seed chain
and settings `pars refine` would use, and records one JSON per candidate; the
selection rules are scored afterwards from these records alone.

Nothing here is a criterion: the runs were already looked at while PARS was being
reviewed, so any rule compared on them is post hoc.

Seed chain: arc-length resample of the whole candidate window to refinement_images,
exactly what `pars select` writes to representatives/.

    python docs/experiments/pars_ground_truth.py --runs p1_round7 p1_committor1 \
        --out runs/pars_ground_truth

One candidate costs ~1e4 MLP calls (~15 min at 44 ms/call for a 9-atom molecule), so
the full set runs one process per candidate, e.g.

    ... --runs p1_round7 --only <candidate id> --device cpu --threads 2
"""

import argparse
import dataclasses
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "pars" / "src"))

from pars.config import ParsConfig  # noqa: E402
from pars.iox import read_json, write_json_atomic  # noqa: E402
from pars.paths.candidate import build_candidates  # noqa: E402
from pars.paths.ingest import ingest_run  # noqa: E402
from pars.paths.resample import resample_positions  # noqa: E402
from pars.refinement.ase_neb import AseNebBackend  # noqa: E402
from prrs.calculators import mace_factory  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--runs", nargs="+", required=True, help="names under runs/")
    parser.add_argument("--out", required=True)
    parser.add_argument("--limit", type=int, default=None, help="first N candidates per run")
    parser.add_argument("--only", default=None, help="refine just this candidate id")
    parser.add_argument("--device", default=None, help="override the run's device")
    parser.add_argument("--threads", type=int, default=None, help="torch CPU threads")
    args = parser.parse_args()
    if args.threads:
        import torch

        torch.set_num_threads(args.threads)
    config = ParsConfig()
    for name in args.runs:
        run_dir = ROOT / "runs" / name
        calculator_spec = read_json(run_dir / "config.json")["calculator"]
        if calculator_spec["factory"] != "prrs.calculators:mace_factory":
            raise ValueError(f"{name}: unexpected calculator {calculator_spec['factory']}")
        kwargs = dict(calculator_spec["kwargs"])
        if args.device:  # a different device changes speed, not the model or dtype
            kwargs["device"] = args.device
        calculator_spec = {**calculator_spec, "kwargs": kwargs}
        candidates = build_candidates(ingest_run(run_dir, config), config)
        if args.only:
            candidates = [c for c in candidates if c.candidate_id == args.only]
        out = Path(args.out) / name
        if all((out / f"{c.candidate_id}.json").exists() for c in candidates[: args.limit]):
            continue
        backend = AseNebBackend(mace_factory(**kwargs))
        out.mkdir(parents=True, exist_ok=True)
        for candidate in candidates[: args.limit]:
            record_path = out / f"{candidate.candidate_id}.json"
            if record_path.exists():  # resume: one finished record per candidate
                continue
            seed = dataclasses.replace(
                candidate,
                positions=resample_positions(candidate.positions, config.refinement_images),
                energies=np.zeros(config.refinement_images),
            )
            start = time.perf_counter()
            refined = backend.refine(seed, config)
            record = refined.to_summary()
            record["wall_s"] = time.perf_counter() - start
            record["energies_eV"] = (
                None if refined.energies is None else refined.energies.tolist()
            )
            record["chain_A"] = (
                None if refined.positions is None else refined.positions.tolist()
            )
            record["numbers"] = [int(z) for z in candidate.numbers]
            record["calculator"] = calculator_spec
            write_json_atomic(record_path, record)
            print(
                f"{name} {candidate.candidate_id}: {refined.status.value} "
                f"barrier={refined.barrier_eV} ({record['wall_s']:.0f}s) {refined.message}",
                flush=True,
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
