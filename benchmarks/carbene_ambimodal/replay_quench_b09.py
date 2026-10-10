"""Replay the quench of failed b09 trials with a larger step cap, through the runner's own chain.

Several b09 trials ended `quench_failed` after the 500-step FIRE cap; p001 (compress 2-15,
1.5 A) had reached the (2+1) cyclopropane when it was cut off. This restarts each listed trial
from its raw_endpoint (the frame the quench started from, momenta zeroed as run_trial does)
and runs the same steps run_trial and the search run after the free segment:
runner.quench (FIRE, then the soft-mode polish) -> state.classify against the source ->
runner.confirm_minimum. Only quench_steps changes; everything else is the run's config.
Each finished trial is written to <out>/<label>.json; rerunning the same command skips those.

Usage: python benchmarks/carbene_ambimodal/replay_quench_b09.py \
    <source extxyz> <config json> <quench_steps> <out dir> <label>=<raw_endpoint.extxyz> ...
"""

import json
import sys
import time
from dataclasses import replace
from pathlib import Path

import numpy as np
from ase.io import read, write

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path[:0] = [str(HERE), str(ROOT / "src"), str(ROOT)]
from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from classify_b09 import DFT, atoms_of_interest, heavy_rmsd, markers  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402
from prrs.io import snapshot  # noqa: E402
from prrs.reliability import GuardedCalculator  # noqa: E402
from prrs.runner import confirm_minimum, quench  # noqa: E402
from prrs.state import classify  # noqa: E402


def main():
    source_path, config_path, steps, out = sys.argv[1:5]
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    base = SearchConfig.from_dict(json.loads(Path(config_path).read_text()))
    config = replace(base, quench_steps=int(steps))
    ref = {name: read(HERE / "structures" / "ts" / f"{name}.xyz") for name in DFT}
    idx = atoms_of_interest(ref["TS-2OMe"])
    calc = polar1()
    factory = lambda: calc  # noqa: E731
    source = read(source_path)
    source.info.update(charge=config.total_charge, spin=config.multiplicity)
    source.calc = calc
    e0 = source.get_potential_energy()
    rows = []
    for arg in sys.argv[5:]:
        label, path = arg.split("=", 1)
        done = out / f"{label}.json"
        if done.exists():  # rerunning the same command skips the finished trials
            rows.append(json.loads(done.read_text()))
            continue
        atoms = read(path)
        atoms.info.update(charge=config.total_charge, spin=config.multiplicity)
        physical = GuardedCalculator(factory(), config, e0)
        atoms.calc = physical
        atoms.set_momenta(np.zeros((len(atoms), 3)))
        start = time.time()
        converged, diagnostics = quench(atoms, physical, config)
        row = dict(
            label=label,
            raw_endpoint=path,
            quench_steps_cap=config.quench_steps,
            converged=bool(converged),
            steps=diagnostics.get("steps"),
            reason=diagnostics.get("reason"),
            force_converged=diagnostics.get("force_converged"),
            polish_rounds=len(diagnostics.get("rounds", [])),
            fmax_eV_A=float(np.linalg.norm(atoms.get_forces(), axis=1).max()),
            dE_eV=float(atoms.get_potential_energy() - e0),
            markers=markers(atoms, idx),
            quench_seconds=time.time() - start,
        )
        endpoint = snapshot(atoms, physical)
        endpoint.info["quench"] = diagnostics
        write(str(out / f"{label}_endpoint.extxyz"), endpoint)
        if converged:
            classification, fingerprint = classify(source, endpoint, config)
            confirmed, check = confirm_minimum(endpoint, factory, config, config.seed)
            row.update(
                classification=classification,
                formed=fingerprint.get("formed"),
                broken=fingerprint.get("broken"),
                minimum=bool(confirmed),
                saddle_order=check.get("saddle_order"),
                rmsd_to_source=heavy_rmsd(source, endpoint),
            )
        row["total_seconds"] = time.time() - start
        rows.append(row)
        done.write_text(json.dumps(row, indent=1))
        print(json.dumps(row), flush=True)
    (out / "summary.json").write_text(json.dumps(rows, indent=1))


if __name__ == "__main__":
    main()
