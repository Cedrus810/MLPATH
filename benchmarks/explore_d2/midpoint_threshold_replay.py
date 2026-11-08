"""What would a looser bisection-midpoint threshold have returned? Offline, from the
recorded midpoint relaxations of finished D2-HF runs (TODO item 9; no model).

`_relax` returns the first frame whose projected (transverse) fmax meets its threshold.
For a looser threshold tau', that frame is the first recorded frame with fmax_projected <=
tau' -- the recorded trajectory up to that frame is exactly what the looser run would have
computed (same start, same BFGS, same forces), up to fp32 non-reproducibility.

Per tau' in TAUS:
  converged segments (met 1e-3): reaction force at the tau' frame minus at the 1e-3 frame,
    and the steps saved;
  converged pairs (from_left and from_right at the same q, both met 1e-3): the paired
    reaction-force difference at tau' against at 1e-3, the branch-force gate
    (3 sqrt(2) sigma, sigma from the manifest), and whether both sides keep the same sign
    class (> tol, < -tol, or unresolved) as at 1e-3;
  failed segments (never met 1e-3): whether a tau' frame exists, its step, its reaction
    force, and the spread of the reaction force over all frames at or below tau'.

    python midpoint_threshold_replay.py OUT_JSON RUN_DIR [RUN_DIR ...]
"""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / "explore" / "src"), str(ROOT / "src")]
sys.path.insert(0, str(Path(__file__).resolve().parent))

from explore.ts_ef import EFSettings  # noqa: E402
from hs_transverse_diag import MIDPOINT_ORIGINS, segments  # noqa: E402

TAUS = (1.5e-3, 2e-3, 3e-3)


def first_at(seg, tau):
    for i, f in enumerate(seg):
        if f["fmax_projected"] <= tau:
            return i
    return None


def sign(rf, tol):
    return 1 if rf > tol else -1 if rf < -tol else 0


def quant(v):
    v = np.asarray(v, dtype=float)
    if not v.size:
        return None
    return {f"p{p}": float(np.percentile(v, p)) for p in (0, 10, 50, 90, 100)} | {
        "n": int(v.size)
    }


def scan(run_dir):
    m = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    tau0 = m["hf"]["midpoint_transverse_fmax_eV_A"]
    sigma = m["hf"]["reaction_force_sigma"]
    tol = EFSettings().reaction_force_tol
    gate = 3 * np.sqrt(2) * sigma
    conv, pairs, failed = [], [], []
    with open(run_dir / "journal.jsonl", encoding="utf-8") as fh:
        for line in fh:
            if '"HS|' not in line or '"settle"' not in line:
                continue
            ev = json.loads(line)
            if ev.get("event") != "settle" or not ev["task"].startswith("HS|"):
                continue
            rsegs = segments((ev["payload"].get("refine") or {}).get("history") or [])
            mids = [s for s in rsegs if s[0]["origin"] in MIDPOINT_ORIGINS]
            ok = [s for s in mids if s[-1]["fmax_projected"] <= tau0]
            met = {id(s) for s in ok}
            conv += ok
            for a, b in zip(mids, mids[1:]):
                if (
                    a[0]["origin"] == "from_left"
                    and b[0]["origin"] == "from_right"
                    and a[0]["q_target"] == b[0]["q_target"]
                    and id(a) in met
                    and id(b) in met
                ):
                    pairs.append((a, b))
            if (
                mids
                and mids[-1] is rsegs[-1]
                and id(mids[-1]) not in met
                and ev["status"] == "unresolved:transverse_not_converged"
            ):
                failed.append(mids[-1])
    return {
        "tau0": tau0,
        "gate": gate,
        "tol": tol,
        "conv": conv,
        "pairs": pairs,
        "failed": failed,
    }


def summarize(pool):
    tau0, gate, tol = pool["tau0"], pool["gate"], pool["tol"]
    conv, pairs, failed = pool["conv"], pool["pairs"], pool["failed"]
    out = {"tau0": tau0, "gate": gate, "tol": tol}
    out["counts"] = {"converged": len(conv), "pairs": len(pairs), "failed": len(failed)}
    for tau in TAUS:
        d_rf, saved = [], []
        for s in conv:
            i = first_at(s, tau)
            d_rf.append(abs(s[i]["reaction_force"] - s[-1]["reaction_force"]))
            saved.append(s[-1]["relax_step"] - s[i]["relax_step"])
        pair_tau, pair_0, over_gate, sign_change = [], [], 0, 0
        for a, b in pairs:
            ia, ib = first_at(a, tau), first_at(b, tau)
            dt = abs(a[ia]["reaction_force"] - b[ib]["reaction_force"])
            d0 = abs(a[-1]["reaction_force"] - b[-1]["reaction_force"])
            pair_tau.append(dt)
            pair_0.append(d0)
            over_gate += dt > gate
            sign_change += (
                sign(a[ia]["reaction_force"], tol),
                sign(b[ib]["reaction_force"], tol),
            ) != (sign(a[-1]["reaction_force"], tol), sign(b[-1]["reaction_force"], tol))
        reached, steps, rf_abs, rf_spread, unresolved = 0, [], [], [], 0
        for s in failed:
            i = first_at(s, tau)
            if i is None:
                continue
            reached += 1
            steps.append(s[i]["relax_step"])
            rf_abs.append(abs(s[i]["reaction_force"]))
            unresolved += sign(s[i]["reaction_force"], tol) == 0
            rfs = [f["reaction_force"] for f in s if f["fmax_projected"] <= tau]
            rf_spread.append(max(rfs) - min(rfs))
        out[str(tau)] = {
            "converged": {"dRF_vs_tau0": quant(d_rf), "steps_saved": quant(saved)},
            "pairs": {
                "dRF_pair_at_tau": quant(pair_tau),
                "dRF_pair_at_tau0": quant(pair_0),
                "over_gate_at_tau": int(over_gate),
                "over_gate_at_tau0": int(sum(d > gate for d in pair_0)),
                "sign_class_changed": int(sign_change),
            },
            "failed": {
                "reached": reached,
                "step": quant(steps),
                "abs_RF": quant(rf_abs),
                "RF_unresolved_at_tau": int(unresolved),
                "RF_spread_over_frames_le_tau": quant(rf_spread),
            },
        }
    return out


def main():
    out, runs = Path(sys.argv[1]), [Path(p) for p in sys.argv[2:]]
    scans = [scan(r) for r in runs]
    settings = {(x["tau0"], x["gate"], x["tol"]) for x in scans}
    if len(settings) != 1:
        raise SystemExit(f"runs differ in tau0 / gate / tol: {settings}")
    pool = dict(scans[0]) | {
        k: [item for x in scans for item in x[k]] for k in ("conv", "pairs", "failed")
    }
    result = {"runs": [str(r) for r in runs], "pooled": summarize(pool)}
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=1) + "\n", encoding="utf-8")
    print(json.dumps(result["pooled"], indent=1))


if __name__ == "__main__":
    main()
