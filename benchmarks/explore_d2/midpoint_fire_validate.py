"""Acceptance of the FIRE fine stage on frozen brackets: v4 against v5 (v4 + FIRE fine
stage), paired on one machine (docs/PLAN_HF_MIDPOINT_FINE_FIRE_2026-11-04.md; user
2026-11-04: "用原来的 22 个失败项＋15 个成功对照，同机配对跑完整 HS ... 出现候选丢失或明显改变时，
对该配对重复核查，区分算法影响与 fp32 波动").

Items are those of midpoint_coarse_validate.py (same SEED and N_CONTROL give the same 22
failed + 15 control intents). Each item's whole HS runs once per arm (`search` there).
Recheck: for every item whose two arms differ in candidate presence, or whose candidates
differ by more than the run's own duplicate criteria (basin_energy_eV, basin_rmsd_A,
Kabsch RMSD in the input atom order), both arms run REPEATS more times (files `_rK`).

One JSON per (item, arm, repeat); an existing file is skipped, so rerunning the same
command resumes. LIMIT > 0 keeps the first LIMIT items per group (smoke).

    python midpoint_fire_validate.py OUT_DIR BOUND N_CONTROL SEED LIMIT REPEATS RUN_DIR ...
"""

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from midpoint_coarse_validate import items, search  # noqa: E402  (sets sys.path)
from midpoint_cold_replay import load  # noqa: E402

from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.ts_ef import _rmsd  # noqa: E402

ARMS = ("v4", "v5")


def changed(a, b, config):
    ca, cb = a["candidate"], b["candidate"]
    if (ca is None) != (cb is None):
        return True
    if ca is None:
        return False
    de = abs(ca["energy_eV"] - cb["energy_eV"])
    rmsd = _rmsd(np.array(ca["positions"]), np.array(cb["positions"]))
    return de > config.basin_energy_eV or rmsd > config.basin_rmsd_A


def log(rec, out):
    mids = [(m["origin"][-6:], m["coarse_sign"], m["optimizer"]) for m in rec["midpoints"]]
    print(
        json.dumps(
            {"file": out.name, "status": rec["status"], "spent": rec["spent"]}
            | {"wall_s": round(rec["wall_s"], 1), "mids": mids}
        ),
        flush=True,
    )


def main():
    out_dir = Path(sys.argv[1])
    bound = float(sys.argv[2])
    n_control, seed, limit, repeats = map(int, sys.argv[3:7])
    runs = [load(Path(p)) for p in sys.argv[7:]]
    by_name = {r.name: r for r in runs}
    picks = items(runs, n_control, seed)
    if limit > 0:
        picks = [p for p in picks if p[0] == "failed"][:limit] + [
            p for p in picks if p[0] == "control"
        ][:limit]
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "items.json").write_text(json.dumps(picks, indent=1) + "\n")
    calc = polar1("cuda", "float32")

    def one(i, group, task, name, arm, k):
        out = out_dir / f"{group}_{i:02d}_{arm}{'' if k == 0 else f'_r{k}'}.json"
        if out.exists():
            return json.loads(out.read_text(encoding="utf-8"))
        rec = search(by_name[name], task, arm, bound, calc, out)
        log(rec, out)
        return rec

    first = {}
    for i, (group, task, name) in enumerate(picks):
        first[i] = {arm: one(i, group, task, name, arm, 0) for arm in ARMS}
    flagged = [
        i
        for i, (group, task, name) in enumerate(picks)
        if changed(first[i]["v4"], first[i]["v5"], by_name[name].config)
    ]
    (out_dir / "recheck.json").write_text(json.dumps(flagged) + "\n")
    print(json.dumps({"recheck": flagged}), flush=True)
    for i in flagged:
        group, task, name = picks[i]
        for k in range(1, repeats + 1):
            for arm in ARMS:
                one(i, group, task, name, arm, k)


if __name__ == "__main__":
    main()
