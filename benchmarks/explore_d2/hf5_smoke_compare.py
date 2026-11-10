"""HF-5 a0: run_d2.py vs run_b04_d2.py on the same b04 smoke parameters
(docs/PLAN_HF5_CROSS_SYSTEM_2026-11-08.md §3, consistency smoke). Read-only.

    python hf5_smoke_compare.py NEW_RUN_DIR OLD_RUN_DIR

Same frozen settings / HF / start-cluster blocks and the same L0 poses are expected exactly;
outcome classes are expected to match; forward counts may differ by fp32 GPU
non-determinism (10-26 vs 10-28 smokes: HS 610 / 599, HC 181 / 189)."""

import json
import sys
from pathlib import Path

import numpy as np


def load(run):
    run = Path(run)
    poses = [
        json.loads(line) for line in (run / "poses.jsonl").read_text().splitlines() if line
    ]
    return (
        json.loads((run / "manifest.json").read_text()),
        json.loads((run / "report.json").read_text()),
        poses,
    )


def main(new, old):
    mn, rn, pn = load(new)
    mo, ro, po = load(old)
    out = {
        "rule_version": [mn["rule_version"], mo["rule_version"]],
        "settings_equal": mn["settings"] == mo["settings"],
        "hf_equal": mn.get("hf") == mo.get("hf"),
        "start_cluster_equal": mn.get("start_cluster") == mo.get("start_cluster"),
        "config_equal": mn.get("config") == mo.get("config"),
        "pose_ids_equal": [p["pose_id"] for p in pn] == [p["pose_id"] for p in po],
        "pose_max_abs_diff_A": max(
            float(np.abs(np.asarray(a["positions"]) - np.asarray(b["positions"])).max())
            for a, b in zip(pn, po)
        )
        if len(pn) == len(po)
        else None,
        "hosts": [rn.get("hosts"), ro.get("hosts")],
    }
    for k in (
        "population",
        "promotion",
        "hs_status",
        "hc_status",
        "unreconciled",
        "hf_forwards",
        "ledger_total",
    ):
        out[k] = {"new": rn.get(k), "old": ro.get(k)}
    out["accounts_spent"] = {
        "new": {a: v["spent"] for a, v in rn["accounts"].items()},
        "old": {a: v["spent"] for a, v in ro["accounts"].items()},
    }
    out["hf_events"] = {
        side: [
            {
                "saddle_energy_eV": e.get("saddle_energy_eV"),
                "chemical_pair": e.get("chemical_pair"),
            }
            for e in r.get("hf_events", [])
        ]
        for side, r in (("new", rn), ("old", ro))
    }
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main(*sys.argv[1:3])
