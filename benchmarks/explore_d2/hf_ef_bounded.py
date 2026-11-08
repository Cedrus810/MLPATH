"""Bounded real-model check of the HF-2/3 E/F locator on saved reactive brackets
(docs/HF_EF_LOCATOR_2026-10-25.md, "下一关").

For each saved reactive C2 (b04 SN2; b05 E2 from round 1 and the fp32 rerun) the L1 scan
of that intent is read from the frozen journal (read-only: no open_run, nothing appended),
the ridge bracket (peak-1, peak+1) is handed to `refine_ts_ef`, and -- only if a
`converged_candidate` comes back -- an independent full Hessian (property) and a descent
on both sides (connection) are run, each as its own ledger stage. The reference for the
comparison is the old C2 saddle of the same intent: same atom order, Kabsch RMSD and
energy difference. Old b04 / round-1 b05 saddles were refined with the fp64-upcast model.

One resident float32 POLAR-1 per system. Rerunning skips cases whose JSON exists.

    python hf_ef_bounded.py OUT_DIR [SIGMA_eV_A [CASE ...]]

SIGMA is the measured same-branch reaction-force jitter passed to EFSettings
(reaction_force_sigma); omitted, the locator's default gate is used.
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT / "src"), str(ROOT / "explore" / "src")]

from ase import Atoms  # noqa: E402

from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.certify import _kabsch_rmsd  # noqa: E402
from explore.d2 import D2Settings, _intent, _points  # noqa: E402
from explore.ledger import BudgetExhausted, Ledger, view_factory  # noqa: E402
from explore.poses import _jsonable  # noqa: E402
from explore.probe import read_scan  # noqa: E402
from explore.ts_ef import EFSettings, refine_ts_ef  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402
from prrs.reliability import GateRejected  # noqa: E402
from prrs.runner import confirm_minimum, descend_saddle, side_seed  # noqa: E402
from prrs.state import classify, same_basin  # noqa: E402

R = ROOT / "benchmarks/explore_d2/results/2026-10-23"
CASES = [
    ("b04", R / "b04_d2", "distance:0-1:+1:+:L0-v1"),
    ("b05_r1_sA_6-9", R / "b05_d2_sA", "distance:6-9:+1:+:L0-v1"),
    ("b05_r1_sA_10-13", R / "b05_d2_sA", "distance:10-13:+1:+:L0-v1"),
    ("b05_r1_sB_10-13", R / "b05_d2_sB", "distance:10-13:+1:+:L0-v1"),
    ("b05_fp32_sA_10-13", R / "b05_fp32_sA", "distance:10-13:+1:+:L0-v1"),
]
SEARCH_CAP = 6000  # forwards for refine_ts_ef per case (2 x 40 x 200 is the code's worst case)
CHECK_CAP = 4000  # forwards for the Hessian + both descents
SEED = 20261025


def frozen(run_dir, suffix):
    """Read-only view of a finished D2 run: settings, config, pose, row, L1 and L4 settles."""
    manifest = json.loads((run_dir / "manifest.json").read_text())
    raw = dict(manifest["settings"])
    raw["magnitudes_A"] = tuple(raw["magnitudes_A"])
    settings = D2Settings(**raw)
    config = SearchConfig.from_dict(manifest["config"])
    rows = [json.loads(x) for x in (run_dir / "population.jsonl").read_text().splitlines() if x]
    # The same intent suffix exists on every pose: take the one that was certified C2.
    report = json.loads((run_dir / "report.json").read_text())
    task = next(c["task"] for c in report["c2"] if c["task"].endswith("|" + suffix))
    row = next(r for r in rows if r["intent_id"] == task[3:])
    pose = None
    for line in (run_dir / "poses.jsonl").read_text().splitlines():
        r = json.loads(line) if line.strip() else None
        if r and r["pose_id"] == row["pose_id"]:
            pose = Atoms(numbers=r["numbers"], positions=r["positions"])
            pose.info.update(r["info"])
    settles = {}
    with (run_dir / "journal.jsonl").open() as fh:
        for line in fh:
            if row["intent_id"] in line and '"settle"' in line:
                ev = json.loads(line)
                if ev["event"] == "settle":
                    settles[ev["task"][:2]] = ev
    return settings, config, row, pose, settles["L1"], settles["L4"]


def check(candidate, factory, config, ledger, task):
    """Independent property (full Hessian) and connection (descent on both sides)."""
    out = {"property": None, "connection": None}
    with ledger.segment(task=task, stage="hessian", attempt=0) as s:
        _, checked = confirm_minimum(candidate, factory, config, SEED)
        s.outcome = str(checked.get("saddle_order", checked.get("reason")))
    out["property"] = _jsonable({k: v for k, v in checked.items() if k != "unstable_mode"})
    if checked.get("saddle_order") != 1 or checked.get("unstable_mode") is None:
        return out
    mode = np.asarray(checked["unstable_mode"])
    ends, sides = {}, []
    with ledger.segment(task=task, stage="descent", attempt=0) as s:
        for sign in (1, -1):
            end, dchecked, relays = descend_saddle(
                candidate, mode, factory, config, sign, seed=side_seed(SEED, sign)
            )
            ends[sign] = end
            sides.append(
                {
                    "sign": sign,
                    "reached_minimum": end is not None,
                    "energy_eV": None if end is None else float(end.get_potential_energy()),
                    "relays": len(relays),
                    "diagnostics": _jsonable(dchecked),
                }
            )
        s.outcome = "done"
    conn = {"sides": sides}
    if all(e is not None for e in ends.values()):
        if same_basin(ends[-1], ends[1], config):
            conn["status"] = "same_basin"
        else:
            label, fingerprint = classify(ends[-1], ends[1], config)
            conn["status"] = "relayed" if any(x["relays"] for x in sides) else "connected"
            conn["label"], conn["fingerprint"] = label, _jsonable(fingerprint)
    else:
        conn["status"] = "descent_failed"
    out["connection"] = conn
    return out


def run_case(name, run_dir, suffix, calc, out_dir, settings_ef=None):
    out = out_dir / f"{name}.json"
    if out.exists():
        return json.loads(out.read_text())
    settings, config, row, pose, l1, l4 = frozen(run_dir, suffix)
    points = _points(l1["payload"])
    ridge = next(t for t in read_scan(points, settings.noise_eV) if t["kind"] == "ridge_suspected")
    k = ridge["peak"]
    ledger = Ledger(out_dir / f"{name}_ledger.jsonl", arm="HF", budget=SEARCH_CAP + CHECK_CAP)
    factory = view_factory(calc, ledger)
    t0 = time.monotonic()
    search = refine_ts_ef(
        pose.copy(),
        _intent(row),
        (points[k - 1], points[k + 1]),
        factory,
        config,
        ledger,
        task=name,
        settings=settings_ef,
        ceiling_eV=points[0].energy_eV + settings.certify_max_rise_eV,
    )
    record = {
        "case": name,
        "run": str(run_dir.relative_to(ROOT)),
        "intent": row["intent_id"],
        "bracket_q": [points[k - 1].q_target, points[k + 1].q_target],
        "scan_peak_energy_eV": points[k].energy_eV,
        "search": search.to_dict(),
        "search_forwards": search.spent,
        "search_wall_s": search.wall_s,
    }
    old = l4["payload"]
    ref = np.asarray(old["saddle_positions"], float)
    record["reference"] = {
        "old_status": l4["status"],
        "old_label": (old.get("event") or {}).get("label"),
        "old_saddle_energy_eV": old["saddle_energy_eV"],
        "old_certify_forwards": l4["spent"],
    }
    if search.status == "converged_candidate":
        record["vs_old_saddle"] = {
            "rmsd_kabsch_A": _kabsch_rmsd(search.positions, ref),
            "dE_eV": search.energy_eV - old["saddle_energy_eV"],
        }
        cand = pose.copy()
        cand.set_constraint()
        cand.positions = search.positions
        cand.calc = factory()  # confirm_minimum reads the candidate's own energy first
        before = ledger.spent
        try:
            record["checks"] = check(cand, factory, config, ledger, name)
        except BudgetExhausted:
            record["checks"] = {"status": "budget_exhausted"}
        except GateRejected as exc:
            record["checks"] = {"status": f"gate_rejected:{exc.code}"}
        record["check_forwards"] = ledger.spent - before
    else:
        # Where the unresolved search stopped, against the old saddle, for orientation only.
        last = search.bracket[-1]["positions"] if search.bracket else None
        record["vs_old_saddle_last_bracket"] = (
            None if last is None else {"rmsd_kabsch_A": _kabsch_rmsd(np.asarray(last), ref)}
        )
    record["wall_s"] = time.monotonic() - t0
    out.write_text(json.dumps(_jsonable(record), indent=2, allow_nan=False) + "\n")
    return record


def main():
    out_dir = Path(sys.argv[1])
    out_dir.mkdir(parents=True, exist_ok=True)
    t = time.monotonic()
    calc = polar1("cuda", "float32")
    print(f"model load {time.monotonic() - t:.1f} s", flush=True)
    sigma = float(sys.argv[2]) if len(sys.argv) > 2 else None
    only = set(sys.argv[3:])
    settings_ef = None if sigma is None else EFSettings(reaction_force_sigma=sigma)
    for name, run_dir, suffix in CASES:
        if only and name not in only:
            continue
        r = run_case(name, run_dir, suffix, calc, out_dir, settings_ef)
        s = r["search"]
        print(
            name,
            s["status"],
            f"search {r['search_forwards']} fw {r['search_wall_s']:.0f} s",
            r.get("vs_old_saddle"),
            (r.get("checks") or {}).get("property", {}).get("saddle_order")
            if isinstance(r.get("checks"), dict) and r["checks"].get("property")
            else r.get("checks"),
            ((r.get("checks") or {}).get("connection") or {}).get("status")
            if isinstance(r.get("checks"), dict)
            else None,
            ((r.get("checks") or {}).get("connection") or {}).get("label")
            if isinstance(r.get("checks"), dict)
            else None,
            flush=True,
        )


if __name__ == "__main__":
    main()
