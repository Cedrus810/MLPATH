"""HF-5 a5 / n4: reference recall and independent property check of D2-HF events
(docs/PLAN_HF5_CROSS_SYSTEM_2026-11-08.md §3, §7).

For every HC `connected` event of the given runs (journal): symmetric RMSD of its saddle to
each reference saddle under the automorphisms of that reference's own graph, and energy
difference against the reference single point on the same float32 model. A reference is
matched when RMSD <= basin_rmsd_A and |dE| <= basin_energy_eV (the runs' own HF duplicate
rule). Then the saddle order by `prrs.runner.confirm_minimum` (central-difference Hessian,
rigid-body modes projected out), as hf_ef_bounded.py did on 2026-10-25. This runs after the
search and never feeds back into it; an undecided Hessian is recorded as such.

    python hf_event_validate.py OUT_JSON --refs REF.extxyz [...] --runs RUN_DIR [...]
"""

import dataclasses
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT / "src"), str(ROOT / "explore" / "src")]

from ase import Atoms  # noqa: E402
from ase.io import read  # noqa: E402

from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.ledger import Ledger, view_factory  # noqa: E402
from prrs.chemistry import automorphisms, symmetric_rmsd  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402
from prrs.runner import confirm_minimum  # noqa: E402
from prrs.search import _stamped_with_charge  # noqa: E402

SEED = 20261109

# TODO item 12 (user 2026-11-11「3 做了」): the analytic-Hessian cross-check exists to catch
# unit / layout errors, which give O(1) relative deviations (a factor-2 error gives 0.5).
# fp32 POLAR-1 measured deviations at the default 0.01 A probe step: b04 events 0.004-0.021
# (the rejected ones 0.016-0.021, 2026-11-09), b05 median 0.0065-0.012 (2026-10-30). The
# default gate 1e-2 sits inside that noise and turned true saddles into "undecided". 0.1 is
# ~5x the largest measured fp32 deviation and ~5x below a factor-2 unit error. Used only for
# post-search validation here; search and PRRS defaults are unchanged.
VALIDATION_HESSIAN_CHECK = 0.1


def validation_config(config):
    return dataclasses.replace(config, analytic_hessian_check=VALIDATION_HESSIAN_CHECK)


_POSE_GROUPS = {}


def match_group(atoms, run_dir, config):
    """Permutations for symmetric RMSD: the automorphisms of the structure's own graph united
    with those of the run's pose (reactant) graph. The saddle graph alone misses symmetry-
    equivalent saddles whose reacting atom carries another label (b05: which methyl H
    moves; found 2026-11-12, sB events against sA events). D2's own HF duplicate rule also
    unites saddle and endpoint groups (TODO item 6); the 11-06 identity used the pose group."""
    key = str(run_dir)
    if key not in _POSE_GROUPS:
        r = json.loads((Path(run_dir) / "poses.jsonl").read_text().splitlines()[0])
        pose = Atoms(numbers=r["numbers"], positions=r["positions"])
        _POSE_GROUPS[key], _ = automorphisms(
            pose, config.bond_scale, config.active_atoms, config.automorphism_limit
        )
    own, _ = automorphisms(
        atoms, config.bond_scale, config.active_atoms, config.automorphism_limit
    )
    unique = {tuple(int(i) for i in g) for g in list(own) + list(_POSE_GROUPS[key])}
    return [np.array(g) for g in sorted(unique)]


def events(run_dir):
    out = []
    with open(run_dir / "journal.jsonl", encoding="utf-8") as fh:
        for line in fh:
            if '"HC|' not in line or '"settle"' not in line:
                continue
            ev = json.loads(line)
            if ev.get("event") != "settle" or ev["status"] != "connected":
                continue
            p = ev["payload"]
            out.append(
                {
                    "run": run_dir.name,
                    "task": ev["task"],
                    "energy_eV": p["saddle_energy_eV"],
                    "broken": p["fingerprint"]["broken"],
                    "formed": p["fingerprint"]["formed"],
                    "positions": np.asarray(p["saddle_positions"], float),
                }
            )
    return out


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


def main():
    args = sys.argv[1:]
    out = Path(args[0])
    refs = args[args.index("--refs") + 1 : args.index("--runs")]
    runs = [Path(p) for p in args[args.index("--runs") + 1 :]]
    manifest = json.loads((runs[0] / "manifest.json").read_text(encoding="utf-8"))
    config = validation_config(SearchConfig.from_dict(manifest["config"]))
    numbers = json.loads((runs[0] / "poses.jsonl").read_text().splitlines()[0])["numbers"]
    calc = polar1("cuda", "float32")
    ledger = Ledger(out.with_suffix(".ledger.jsonl"), arm="V", replay=False)
    factory = view_factory(calc, ledger)

    def atoms_at(positions):
        a = _stamped_with_charge(Atoms(numbers=numbers, positions=positions), config)
        a.info["spin"] = config.multiplicity
        a.calc = factory()
        return a

    def order(a, task):
        t0 = time.monotonic()
        with ledger.segment(task=task, stage="hessian", attempt=0) as seg:
            _, checked = confirm_minimum(a, factory, config, SEED)
            seg.outcome = str(checked.get("saddle_order", checked.get("reason")))
        rec = {k: v for k, v in checked.items() if k not in ("unstable_mode",)}
        return _jsonable(rec) | {"wall_s": time.monotonic() - t0}

    reference = []
    for path in refs:
        r = read(path)
        if list(r.numbers) != list(numbers):
            raise ValueError(f"{path}: atom order {list(r.numbers)} != run order {numbers}")
        a = atoms_at(r.positions)
        group = match_group(a, runs[0], config)
        reference.append(
            {
                "path": str(path),
                "energy_eV_fp32": float(a.get_potential_energy()),
                "group": group,
                "atoms": a,
                "order": order(a, f"ref|{Path(path).name}"),
            }
        )
    rows = []
    for run_dir in runs:
        for ev in events(run_dir):
            a = atoms_at(ev["positions"])
            match = []
            for ref in reference:
                rmsd = symmetric_rmsd(a, ref["atoms"], ref["group"], config.active_atoms)
                de = ev["energy_eV"] - ref["energy_eV_fp32"]
                match.append(
                    {
                        "ref": ref["path"],
                        "symmetric_rmsd_A": float(rmsd),
                        "dE_eV": float(de),
                        "matched": bool(
                            rmsd <= config.basin_rmsd_A and abs(de) <= config.basin_energy_eV
                        ),
                    }
                )
            rows.append(
                {k: ev[k] for k in ("run", "task", "energy_eV", "broken", "formed")}
                | {"reference": match, "order": order(a, ev["task"])}
            )
    ledger.close()
    result = {
        "basin_rmsd_A": config.basin_rmsd_A,
        "basin_energy_eV": config.basin_energy_eV,
        "references": [
            {k: r[k] for k in ("path", "energy_eV_fp32", "order")}
            | {"automorphisms": len(r["group"])}
            for r in reference
        ],
        "events": rows,
        "forwards": ledger.spent,
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=1) + "\n", encoding="utf-8")
    brief = [
        (
            e["run"],
            e["energy_eV"],
            [
                (round(m["symmetric_rmsd_A"], 4), round(m["dE_eV"], 4), m["matched"])
                for m in e["reference"]
            ],
            e["order"].get("saddle_order", e["order"].get("reason")),
        )
        for e in rows
    ]
    print(
        json.dumps({"references": result["references"], "events": brief}, indent=1, default=str)
    )


if __name__ == "__main__":
    main()
