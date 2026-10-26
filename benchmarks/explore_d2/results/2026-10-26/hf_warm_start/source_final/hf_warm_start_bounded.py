"""Bounded paired search-cost comparison on the frozen J1 brackets.

One resident fp32 POLAR-1, unchanged noise/force/branch tolerances. Each search has
its own 3000-forward Ledger. Checks run separately, only for requested warm cases.
No historical run is opened for writing; existing outputs are retained.

    python hf_warm_start_bounded.py OUT --modes cold curvature predictor warm \
        --cases b04 b05_r1_sA_6-9 --verify b05_r1_sA_6-9

The cold switches reproduce the pre-optimization relaxation strategy. This is a
small performance experiment, not a new discovery or general equivalence claim.
"""

import argparse
from contextlib import closing
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import socket
import sys
import time

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from hf_ef_bounded import CASES, ROOT, check, frozen  # noqa: E402
from hf_jump_bounded import NOISE, SIGMA_TIGHT  # noqa: E402

from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.certify import _kabsch_rmsd  # noqa: E402
from explore.d2 import _intent, _points  # noqa: E402
from explore.ledger import BudgetExhausted, Ledger, view_factory  # noqa: E402
from explore.poses import _jsonable  # noqa: E402
from explore.probe import read_scan  # noqa: E402
from explore.ts_branch import BranchSettings, continue_branches_ef  # noqa: E402
from explore.ts_ef import EFSettings  # noqa: E402
from prrs.reliability import GateRejected  # noqa: E402

MODES = {
    "cold": (False, False),
    "curvature": (True, False),
    "predictor": (False, True),
    "warm": (True, True),
}
SEARCH_CAP = 3000
CHECK_CAP = 4000


def stats(history):
    finals = []
    for p in history:
        if finals and (p["origin"], p["q_target"]) == (
            finals[-1]["origin"],
            finals[-1]["q_target"],
        ):
            finals[-1] = p
        else:
            finals.append(p)
    return {
        "points": len(finals),
        "steps": [p["relax_step"] for p in finals],
        "median_steps": float(np.median([p["relax_step"] for p in finals])) if finals else None,
        "curvature_reused_points": sum(p["curvature_reused"] for p in finals),
        "predicted_points": sum(p["predictor_step_A"] > 0 for p in finals),
    }


def run_case(case, mode, calc, out_dir, verify):
    name, run_dir, suffix = case
    output = out_dir / f"{name}-{mode}.json"
    if output.exists():
        record = json.loads(output.read_text())
        if (
            verify
            and mode == "warm"
            and record["result"]["status"] == "converged_candidate"
            and "checks" not in record
        ):
            _, config, _, pose, _, _ = frozen(run_dir, suffix)
            verify_record(record, pose, config, calc, out_dir)
            output.write_text(json.dumps(_jsonable(record), indent=2, allow_nan=False) + "\n")
        return record
    settings, config, row, pose, l1, l4 = frozen(run_dir, suffix)
    points = _points(l1["payload"])
    ridge = next(
        t for t in read_scan(points, settings.noise_eV) if t["kind"] == "ridge_suspected"
    )
    k = ridge["peak"]
    ef = EFSettings(
        reaction_force_sigma=SIGMA_TIGHT[name[:3]],
        midpoint_transverse_fmax_eV_A=1e-3,
        reuse_curvature=MODES[mode][0],
        predict_tangent=MODES[mode][1],
    )
    branch = BranchSettings(geometry_noise_A=NOISE[name[:3]][1])
    sources = [
        "explore/src/explore/ts_ef.py",
        "explore/src/explore/ts_branch.py",
        "benchmarks/explore_d2/hf_warm_start_bounded.py",
    ]
    with closing(
        Ledger(out_dir / f"{name}-{mode}-search.jsonl", "HF-perf", budget=SEARCH_CAP)
    ) as ledger:
        result = continue_branches_ef(
            pose,
            _intent(row),
            (points[k - 1], points[k + 1]),
            view_factory(calc, ledger),
            config,
            ledger,
            name,
            branch_settings=branch,
            settings=ef,
            ceiling_eV=points[0].energy_eV + settings.certify_max_rise_eV,
        )
        record = {
            "case": name,
            "mode": mode,
            "host": socket.gethostname(),
            "model": "MACE-POLAR-1-M",
            "precision": "float32",
            "source_run": str(run_dir),
            "source_sha256": {
                p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in sources
            },
            "settings": asdict(ef),
            "branch_settings": asdict(branch),
            "search_forward_cap": SEARCH_CAP,
            "search_forwards": result.spent,
            "search_wall_s": result.wall_s,
            "result": result.to_dict(),
            "continuation": stats(result.history),
            "bisection": stats((result.refine or {}).get("history", [])),
            "property_status": "not_run",
            "connection_status": "not_run",
        }
    if result.status == "converged_candidate":
        old = l4["payload"]
        pos = np.asarray(result.refine["positions"], float)
        record["vs_old_saddle"] = {
            "rmsd_kabsch_A": _kabsch_rmsd(pos, np.asarray(old["saddle_positions"])),
            "dE_eV": result.refine["energy_eV"] - old["saddle_energy_eV"],
        }
        if verify and mode == "warm":
            verify_record(record, pose, config, calc, out_dir)
    output.write_text(json.dumps(_jsonable(record), indent=2, allow_nan=False) + "\n")
    return record


def verify_record(record, pose, config, calc, out_dir):
    name, mode = record["case"], record["mode"]
    with closing(
        Ledger(out_dir / f"{name}-{mode}-checks.jsonl", "HF-check", budget=CHECK_CAP)
    ) as ledger:
        factory = view_factory(calc, ledger)
        candidate = pose.copy()
        candidate.set_constraint()
        candidate.positions = np.asarray(record["result"]["refine"]["positions"], float)
        candidate.calc = factory()
        started = time.monotonic()
        try:
            record["checks"] = check(candidate, factory, config, ledger, name)
        except BudgetExhausted:
            record["checks"] = {"status": "budget_exhausted"}
        except GateRejected as exc:
            record["checks"] = {"status": f"gate_rejected:{exc.code}"}
        record["check_forwards"] = ledger.spent
        record["check_wall_s"] = time.monotonic() - started
        order = (record["checks"].get("property") or {}).get("saddle_order")
        record["property_status"] = (
            ("index_one" if order == 1 else f"order:{order}")
            if order is not None
            else record["checks"].get("status", "unresolved")
        )
        record["connection_status"] = (record["checks"].get("connection") or {}).get(
            "status", "not_run"
        )
        record["check_driver_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--modes", nargs="+", choices=MODES, default=list(MODES))
    parser.add_argument(
        "--cases", nargs="+", choices=[c[0] for c in CASES], default=["b04", "b05_r1_sA_6-9"]
    )
    parser.add_argument("--verify", nargs="*", choices=[c[0] for c in CASES], default=[])
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    calc = polar1("cuda", "float32")
    print(f"resident fp32 model loaded in {time.monotonic() - started:.1f} s", flush=True)
    for case in CASES:
        if case[0] not in args.cases:
            continue
        for mode in args.modes:
            r = run_case(case, mode, calc, args.output, case[0] in args.verify)
            print(
                case[0],
                mode,
                r["result"]["status"],
                r["search_forwards"],
                f"{r['search_wall_s']:.1f} s",
                "steps",
                r["continuation"]["median_steps"],
                r.get("vs_old_saddle"),
                r["property_status"],
                r["connection_status"],
                flush=True,
            )


if __name__ == "__main__":
    main()
