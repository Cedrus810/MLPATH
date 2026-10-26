"""K3 of docs/PLAN_HF_CONNECT_2026-10-26.md: E/F-only connection tracing on the five
warm-start TS candidates (results/2026-10-26/hf_warm_start/*-warm.json, read-only).

Starts are the nearest resolved frames of the candidate's own branch: the continuation
side that bracketed the maximum plus the converged bisection frames. Compared against the
old C2 of the same intent (bond changes, unordered sides) and against the old,
Hessian-dependent check cost recorded for these candidates.

    python hf_connect_bounded.py OUT_DIR [CASE ...]
"""

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from hf_ef_bounded import CASES, frozen  # noqa: E402  (also sets sys.path)

from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.ledger import Ledger, view_factory  # noqa: E402
from explore.ts_connect import ConnectSettings, branch_starts, trace_sides_ef  # noqa: E402

HERE = Path(__file__).resolve().parent
WARM = HERE / "results/2026-10-26/hf_warm_start"
OLD_CHECK = HERE / "results/2026-10-25/hf_jump_bounded_tau1e-3"
# Same-branch reaction-force sigma, larger of the measured values per system
# (DONE_LOG 2026-10-25: b04 0.0034 / 0.0042 near root; b05 0.0025 / 0.0022 near root).
SIGMA = {"b04": 0.0042, "b05": 0.0025}
CAP = 6000
TRANSVERSE = 3e-3  # endpoint / continuation threshold of the warm-start search


def edges(fp):
    return (
        frozenset(tuple(e) for e in fp.get("broken") or []),
        frozenset(tuple(e) for e in fp.get("formed") or []),
    )


def same_event(a, b):
    """Same bond changes with either side ordering."""
    if a is None or b is None:
        return None
    ea, eb = edges(a), edges(b)
    return ea == eb or ea == (eb[1], eb[0])


def run_case(name, run_dir, suffix, calc, out_dir):
    out = out_dir / f"{name}.json"
    if out.exists():
        return json.loads(out.read_text())
    warm = json.loads((WARM / f"{name}-warm.json").read_text())
    res = warm["result"]
    ref = res["refine"]
    _, config, row, pose, _, l4 = frozen(run_dir, suffix)
    frames = list(ref["history"])
    for side in res["sides"]:
        if side.get("status") == "maximum_bracketed":
            frames += side["points"]
    sigma = SIGMA[name[:3]]
    q_star = ref["metrics"]["q_actual"]
    starts = branch_starts(
        frames, q_star, sigma, transverse_fmax=TRANSVERSE, descent_fmax=config.quench_fmax_eV_A
    )
    ledger = Ledger(out_dir / f"{name}_ledger.jsonl", arm="HF-connect", budget=CAP)
    t0 = time.monotonic()
    result = trace_sides_ef(
        pose.copy(),
        ref["positions"],
        ref["energy_eV"],
        starts,
        ref["direction_cartesian"],
        view_factory(calc, ledger),
        config,
        ledger,
        name,
        settings=ConnectSettings(reaction_force_sigma=sigma),
    )
    old_fp = (l4["payload"].get("event") or {}).get("fingerprint")
    old_check = None
    if (OLD_CHECK / f"{name}.json").exists():
        old_check = json.loads((OLD_CHECK / f"{name}.json").read_text()).get("check_forwards")
    record = {
        "case": name,
        "warm_source": str((WARM / f"{name}-warm.json").relative_to(HERE)),
        "sigma": sigma,
        "connection": result.to_dict(),
        "forwards": result.spent,
        "wall_s": time.monotonic() - t0,
        "old_c2": {
            "label": (l4["payload"].get("event") or {}).get("label"),
            "broken": (old_fp or {}).get("broken"),
            "formed": (old_fp or {}).get("formed"),
        },
        "same_event_as_old_c2": same_event(result.fingerprint, old_fp),
        "old_check_forwards": {
            "warm_run": warm.get("check_forwards"),
            "tau1e-3_run": old_check,
        },
    }
    out.write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")
    return record


def main():
    out_dir = Path(sys.argv[1])
    out_dir.mkdir(parents=True, exist_ok=True)
    only = set(sys.argv[2:])
    t = time.monotonic()
    calc = polar1("cuda", "float32")
    print(f"model load {time.monotonic() - t:.1f} s", flush=True)
    for name, run_dir, suffix in CASES:
        if only and name not in only:
            continue
        r = run_case(name, run_dir, suffix, calc, out_dir)
        c = r["connection"]
        print(
            name,
            c["status"],
            c["label"],
            f"{r['forwards']} fw {r['wall_s']:.0f} s",
            "same_as_old:",
            r["same_event_as_old_c2"],
            [
                (
                    s["side"],
                    s.get("status"),
                    s.get("start", {}).get("source"),
                    (s.get("descent") or {}).get("steps"),
                )
                for s in c["sides"]
            ],
            "old check fw:",
            r["old_check_forwards"],
            flush=True,
        )


if __name__ == "__main__":
    main()
