"""N1 reference saddles on float32 POLAR-1 (docs/PLAN_HF5_CROSS_SYSTEM_2026-11-08.md §5 n0;
evaluation protocol §3.2 steps 2-5). Reference side only; never shown to the search.

  python reference.py BUILD_DIR

For each guess_<key>.xyz written by build.py:
  1. `prrs optimize --transition` (geomeTRIC, GAU, maxiter 300) on float32 POLAR-1, as
     carbene_ambimodal/optimize_b08_b09.sh did for b09; an existing ts_<key>/ is reused.
  2. if CONVERGED: saddle order by `prrs.runner.confirm_minimum` (central-difference Hessian,
     rigid-body modes projected out);
  3. if order 1: `descend_saddle` on both signs of the unstable mode; each end is labelled by
     the two forming distances C1-C6 / C4-C7 (both < 1.7 A: product; both > 2.8 A:
     reactant side; otherwise "other") and its energy;
  4. pairwise symmetric RMSD between the order-1 saddles (automorphisms of each saddle's own
     graph) to see whether they are distinct.
Writes BUILD_DIR/reference.json. A Hessian that cannot decide the order is recorded as
undecided, not as a pass.
"""

import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT / "src"), str(ROOT / "explore" / "src")]

from ase.io import read  # noqa: E402

from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.ledger import Ledger, view_factory  # noqa: E402
from prrs.chemistry import automorphisms, symmetric_rmsd  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402
from prrs.runner import confirm_minimum, descend_saddle, side_seed  # noqa: E402
from prrs.search import _stamped_with_charge  # noqa: E402

SEED = 20261109
PY = "/home/ruigengji/miniforge3/envs/openmm_dev/bin/python"
FORMING = ((0, 11), (3, 12))  # build.py module order: C1-C6, C4-C7


def _jsonable(x):
    if isinstance(x, dict):
        return {k: _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, np.generic):
        return x.item()
    return x


def forming(a):
    return [float(a.get_distance(i, j)) for i, j in FORMING]


def label(d):
    if all(x < 1.7 for x in d):
        return "product"
    if all(x > 2.8 for x in d):
        return "reactant_side"
    return "other"


def optimize(guess, dest, kwargs):
    if (dest / "result.json").exists():
        return json.loads((dest / "result.json").read_text())
    if dest.exists():
        return {"status": "incomplete_previous_attempt", "output": str(dest)}
    cmd = [
        PY, "-m", "prrs", "optimize", str(guess), "--output", str(dest),
        "--calculator", "benchmarks.b04_sn2_chloride.calculator:polar1",
        "--calculator-kwargs", str(kwargs), "--transition", "--convergence-set", "GAU",
        "--maxiter", "300", "--charge", "0", "--multiplicity", "1",
    ]  # fmt: skip
    with open(dest.parent / f"{dest.name}.log", "w") as log:
        env = dict(os.environ, PYTHONPATH=f"{ROOT}:{ROOT / 'src'}")  # as optimize_b08_b09.sh
        rc = subprocess.run(
            cmd, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT
        ).returncode
    if (dest / "result.json").exists():
        return json.loads((dest / "result.json").read_text())
    return {"status": f"no_result_exit_{rc}"}


def main():
    d = Path(sys.argv[1]).resolve()
    kwargs = d / "polar1_float32.json"
    kwargs.write_text(json.dumps({"device": "cuda", "default_dtype": "float32"}))
    config = SearchConfig(charge_sensitive=True, total_charge=0, multiplicity=1)
    calc = polar1("cuda", "float32")
    ledger = Ledger(d / "reference.ledger.jsonl", arm="R", replay=False)
    factory = view_factory(calc, ledger)

    def stamped(a):
        a = _stamped_with_charge(a, config)
        a.info["spin"] = 1
        a.calc = factory()
        return a

    rows, saddles = {}, {}
    for guess in sorted(d.glob("guess_*.xyz")):
        key = guess.stem[len("guess_") :]
        res = optimize(guess, d / f"ts_{key}", kwargs)
        row = {
            "guess_forming_A": forming(read(guess)),
            "optimize": {
                k: res.get(k)
                for k in ("status", "iterations", "force_evaluations", "energy_eV", "fmax_eV_A")
            },
        }
        rows[key] = row
        if res.get("status") != "CONVERGED":
            continue
        ts = stamped(read(d / f"ts_{key}" / "geometry.extxyz"))
        row["ts_forming_A"] = forming(ts)
        row["ts_energy_eV_fp32"] = float(ts.get_potential_energy())
        with ledger.segment(task=f"ref|{key}", stage="hessian", attempt=0) as seg:
            _, checked = confirm_minimum(ts, factory, config, SEED)
            seg.outcome = str(checked.get("saddle_order", checked.get("reason")))
        row["order"] = _jsonable(
            {
                k: checked.get(k)
                for k in (
                    "saddle_order",
                    "reason",
                    "imaginary_wavenumbers_icm",
                    "negative_eigenvalues",
                    "provenance",
                )
            }
        )
        mode = checked.get("unstable_mode")
        if checked.get("saddle_order") != 1 or mode is None:
            continue
        saddles[key] = ts
        ends = []
        with ledger.segment(task=f"ref|{key}", stage="descent", attempt=0) as seg:
            for sign in (1, -1):
                end, dchecked, relays = descend_saddle(
                    ts, np.asarray(mode), factory, config, sign, seed=side_seed(SEED, sign)
                )
                if end is None:
                    ends.append({"sign": sign, "reached_minimum": False, "relays": len(relays)})
                    continue
                f = forming(end)
                ends.append(
                    {
                        "sign": sign,
                        "reached_minimum": True,
                        "energy_eV": float(end.get_potential_energy()),
                        "forming_A": f,
                        "label": label(f),
                        "relays": len(relays),
                    }
                )
            seg.outcome = "done"
        row["ends"] = ends
        row["connects_reactant_and_product"] = sorted(e.get("label") for e in ends) == [
            "product",
            "reactant_side",
        ]
    keys = sorted(saddles)
    pair = {}
    for i, a in enumerate(keys):
        group, _ = automorphisms(
            saddles[a], config.bond_scale, config.active_atoms, config.automorphism_limit
        )
        for b in keys[i + 1 :]:
            pair[f"{a}|{b}"] = float(
                symmetric_rmsd(saddles[a], saddles[b], group, config.active_atoms)
            )
    ledger.close()
    out = {
        "guesses": rows,
        "pairwise_symmetric_rmsd_A": pair,
        "forwards_after_optimize": ledger.spent,
    }
    (d / "reference.json").write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8")
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
