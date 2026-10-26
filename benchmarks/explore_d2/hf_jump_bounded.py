"""J3 (J1 part) of docs/PLAN_HF_JUMP_2026-10-25.md: two-sided branch continuation on the
saved reactive brackets, real model, bounded.

Same cases, bracket choice (peak-1, peak+1), read-only journal access and independent
checks as hf_ef_bounded.py. Per-system noise inputs are the measured values recorded in
DONE_LOG 2026-10-25 (hf_ef_rf_jitter.json, hf_ef_rf_jitter_b05.json, tau = 3e-3):
reaction-force sigma and same-branch geometry noise (max pairwise RMSD).

    python hf_jump_bounded.py OUT_DIR [--mid-tau=1e-3] [CASE ...]

--mid-tau sets EFSettings.midpoint_transverse_fmax_eV_A; sigma is then taken from the
near-root repeats measured at that threshold (hf_ef_rf_jitter_root_*.json, tau 1e-3).
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from hf_ef_bounded import CASES, check, frozen  # noqa: E402  (also sets sys.path)

from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.certify import _kabsch_rmsd  # noqa: E402
from explore.d2 import _intent, _points  # noqa: E402
from explore.ledger import BudgetExhausted, Ledger, view_factory  # noqa: E402
from explore.probe import read_scan  # noqa: E402
from explore.ts_branch import BranchSettings, continue_branches_ef  # noqa: E402
from explore.ts_ef import EFSettings  # noqa: E402
from prrs.reliability import GateRejected  # noqa: E402

NOISE = {  # system: (reaction_force_sigma eV/A, geometry_noise_A), tau 3e-3
    "b04": (0.0034, 0.00062),
    "b05": (0.0025, 0.0070),
}
SIGMA_TIGHT = {"b04": 0.0042, "b05": 0.0022}  # near-root repeats at tau 1e-3
SEARCH_CAP = 8000
CHECK_CAP = 4000


def run_case(name, run_dir, suffix, calc, out_dir, mid_tau=None):
    out = out_dir / f"{name}.json"
    if out.exists():
        return json.loads(out.read_text())
    settings, config, row, pose, l1, l4 = frozen(run_dir, suffix)
    points = _points(l1["payload"])
    ridge = next(
        t for t in read_scan(points, settings.noise_eV) if t["kind"] == "ridge_suspected"
    )
    k = ridge["peak"]
    sigma, geo = NOISE[name[:3]]
    if mid_tau is not None:
        sigma = SIGMA_TIGHT[name[:3]]
    ledger = Ledger(out_dir / f"{name}_ledger.jsonl", arm="HF", budget=SEARCH_CAP + CHECK_CAP)
    factory = view_factory(calc, ledger)
    t0 = time.monotonic()
    result = continue_branches_ef(
        pose.copy(),
        _intent(row),
        (points[k - 1], points[k + 1]),
        factory,
        config,
        ledger,
        task=name,
        branch_settings=BranchSettings(geometry_noise_A=geo),
        settings=EFSettings(reaction_force_sigma=sigma, midpoint_transverse_fmax_eV_A=mid_tau),
        ceiling_eV=points[0].energy_eV + settings.certify_max_rise_eV,
    )
    record = {
        "case": name,
        "run": str(run_dir),
        "intent": row["intent_id"],
        "noise": {"reaction_force_sigma": sigma, "geometry_noise_A": geo, "mid_tau": mid_tau},
        "result": result.to_dict(),
        "forwards": result.spent,
    }
    old = l4["payload"]
    ref = np.asarray(old["saddle_positions"], float)
    record["reference"] = {
        "old_label": (old.get("event") or {}).get("label"),
        "old_saddle_energy_eV": old["saddle_energy_eV"],
    }
    if result.status == "converged_candidate":
        pos = np.asarray(result.refine["positions"], float)
        record["vs_old_saddle"] = {
            "rmsd_kabsch_A": _kabsch_rmsd(pos, ref),
            "dE_eV": result.refine["energy_eV"] - old["saddle_energy_eV"],
        }
        cand = pose.copy()
        cand.set_constraint()
        cand.positions = pos
        cand.calc = factory()
        before = ledger.spent
        try:
            record["checks"] = check(cand, factory, config, ledger, name)
        except BudgetExhausted:
            record["checks"] = {"status": "budget_exhausted"}
        except GateRejected as exc:
            record["checks"] = {"status": f"gate_rejected:{exc.code}"}
        record["check_forwards"] = ledger.spent - before
    for side in result.sides:
        if side.get("status") == "terminated":
            last = np.asarray(side["points"][-1]["positions"], float)
            side_ref = {"terminated_last_vs_old_saddle_A": _kabsch_rmsd(last, ref)}
            side_ref["terminated_last_E_minus_old_eV"] = (
                side["points"][-1]["energy_eV"] - old["saddle_energy_eV"]
            )
            record.setdefault("termination_vs_old", []).append(side_ref)
    record["wall_s"] = time.monotonic() - t0
    out.write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")
    return record


def main():
    out_dir = Path(sys.argv[1])
    out_dir.mkdir(parents=True, exist_ok=True)
    flags = [a for a in sys.argv[2:] if a.startswith("--mid-tau=")]
    mid_tau = float(flags[0].split("=")[1]) if flags else None
    only = {a for a in sys.argv[2:] if not a.startswith("--")}
    t = time.monotonic()
    calc = polar1("cuda", "float32")
    print(f"model load {time.monotonic() - t:.1f} s", flush=True)
    for name, run_dir, suffix in CASES:
        if only and name not in only:
            continue
        r = run_case(name, run_dir, suffix, calc, out_dir, mid_tau)
        res = r["result"]
        sides = [
            (s["origin"], s["status"], s.get("terminated_at_q"), len(s["points"]))
            for s in res["sides"]
        ]
        checks = r.get("checks") or {}
        print(
            name,
            res["branch_status"],
            res["status"],
            f"{r['forwards']} fw {r['wall_s']:.0f} s",
            sides,
            r.get("vs_old_saddle"),
            r.get("termination_vs_old"),
            (checks.get("property") or {}).get("saddle_order"),
            (checks.get("connection") or {}).get("status"),
            (checks.get("connection") or {}).get("label"),
            flush=True,
        )


if __name__ == "__main__":
    main()
