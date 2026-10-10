"""Does a harder push along the carbene approach cross on b09?

The compress coverage run (runs/b09_compress_coverage) pushed carbene C -> Cb to 2.12 A at
its largest amplitude, 1.0 A; the POLAR-1 TS-2OMe has 1.98 A, and every trial came back.
This runs the two approach directions, compress and kick, beyond that grid -- one trial
each through prrs.runner.run_trial, the same call the search makes, from the run's own
source microstate and with its own config -- and reports what each reached.

Usage: python benchmarks/carbene_ambimodal/push_b09.py <run dir> <config json> <out dir>
"""

import json
import sys
from pathlib import Path

import numpy as np
from ase.io import read

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path[:0] = [str(HERE), str(ROOT / "src"), str(ROOT)]
from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from classify_b09 import DFT, atoms_of_interest, heavy_rmsd, markers  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402
from prrs.perturbations import Probe  # noqa: E402
from prrs.runner import confirm_minimum, run_trial  # noqa: E402

# (family, indices, sign, amplitude): sign -1 shortens the C...C distance for both families
PROBES = [
    ("compress", (2, 15), -1, 1.2),
    ("compress", (2, 15), -1, 1.5),
    ("kick", (2, 15), -1, 0.5),
    ("kick", (2, 15), -1, 1.0),
    ("compress", (2, 13), -1, 1.2),
    ("compress", (2, 13), -1, 1.5),
    ("kick", (2, 13), -1, 0.5),
    ("kick", (2, 13), -1, 1.0),
]


def main():
    run, out = Path(sys.argv[1]), Path(sys.argv[3])
    config = SearchConfig.from_dict(json.loads(Path(sys.argv[2]).read_text()))
    ref = {name: read(HERE / "structures" / "ts" / f"{name}.xyz") for name in DFT}
    idx = atoms_of_interest(ref["TS-2OMe"])
    pairs = {k: (idx["carbene"], idx[k]) for k in ("Ca", "Cb", "O_ald")}
    calc = polar1()
    factory = lambda: calc  # noqa: E731
    source = read(run / "chemical" / "c0000_m0000.extxyz")
    source.info.update(charge=config.total_charge, spin=config.multiplicity)
    source.calc = calc
    e0 = source.get_potential_energy()
    rows = []
    for n, (family, indices, sign, amplitude) in enumerate(PROBES):
        probe = Probe(family, indices, sign, amplitude, config.seed)
        trial = f"p{n:03d}"
        outcome = run_trial(source, probe, factory, config, out / "trials" / trial, trial)
        r = outcome.record
        frames = read(out / "trials" / trial / "response.traj", index=":")
        closest = {
            k: float(min(np.linalg.norm(f.positions[i] - f.positions[j]) for f in frames))
            for k, (i, j) in pairs.items()
        }
        row = dict(
            trial=trial,
            probe=f"{family}:{'-'.join(map(str, indices))}:{sign:+d} {amplitude}",
            status=r["status"],
            classification=r["classification"],
            formed=(r.get("fingerprint") or {}).get("formed"),
            broken=(r.get("fingerprint") or {}).get("broken"),
            failure=r["failure"],
            closest=closest,
        )
        if outcome.endpoint is not None:
            end = outcome.endpoint
            end.calc = calc
            confirmed, diagnostics = confirm_minimum(end, factory, config, config.seed)
            row.update(
                minimum=bool(confirmed),
                saddle_order=diagnostics.get("saddle_order"),
                dE_meV=1000 * (end.get_potential_energy() - e0),
                end_markers=markers(end, idx),
                rmsd_to_source=heavy_rmsd(source, end),
                heavy_rmsd_to_dft={k: heavy_rmsd(v, end) for k, v in ref.items()},
            )
        rows.append(row)
        print(json.dumps(row), flush=True)
    (out / "summary.json").write_text(json.dumps(rows, indent=1))


if __name__ == "__main__":
    main()
