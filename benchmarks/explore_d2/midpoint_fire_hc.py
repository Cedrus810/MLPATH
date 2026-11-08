"""HC for the candidates the FIRE fine stage rescued (TODO item 9b-2; user 2026-11-06
"可以"; identity first: results/2026-11-06/candidate_identity/identity.json).

For every failed intent of the v4 / v5 acceptance (midpoint_coarse_validate.items, same
SEED), the whole HS is rerun with HF-EF-bracket-v5 (search records keep no branch frames,
which HC needs), and, when it gives a candidate, HC runs exactly as `d2._hc_task`:
`branch_starts` on the refine history plus bracketed side points, then `trace_sides_ef`
with ConnectSettings(sigma, task seed), each on its own ledger (hs_cap / hc_cap). Recorded:
HS status and candidate, HC status, both descent ends, chemical pair (`d2._hf_identity`),
bond changes, and whether the pair equals that of the reference new-E2 event.

One JSON per item; an existing file is skipped, so rerunning resumes. LIMIT > 0 keeps the
first LIMIT items (smoke).

    python midpoint_fire_hc.py OUT_DIR BOUND N_CONTROL SEED LIMIT RUN_DIR ...
"""

import json
import sys
import time
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent))

from midpoint_coarse_validate import items, settings_for  # noqa: E402  (sets sys.path)
from midpoint_cold_replay import context, load  # noqa: E402

import explore.ts_branch as ts_branch  # noqa: E402
from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.d2 import _classify, _hf_bracket, _hf_identity, _intent, _seed  # noqa: E402
from explore.ledger import Ledger, view_factory  # noqa: E402
from explore.ts_connect import ConnectSettings, branch_starts, trace_sides_ef  # noqa: E402

REFERENCE = ("2026-10-29", "b05_hf_v4_sA", "distance:10-12")  # the new E2 event


def reference_pair(runs_root):
    path = runs_root / REFERENCE[0] / REFERENCE[1] / "journal.jsonl"
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if '"HC|' not in line or REFERENCE[2] not in line:
                continue
            ev = json.loads(line)
            if ev.get("event") == "settle" and ev["status"] == "connected":
                return ev["payload"]["chemical_pair"], ev["payload"]["fingerprint"]
    raise SystemExit("reference new-E2 event not found")


def one(run, task, bound, calc, out):
    row, points, ceiling = context(run, task)
    view = SimpleNamespace(settled=run.settled, settings=run.settings)
    classified = _classify(view, row)
    hs_ledger = Ledger(out.with_suffix(".hs.ledger.jsonl"), arm="HF", budget=run.hf.hs_cap)
    t0 = time.monotonic()
    try:
        hs = ts_branch.continue_branches_ef(
            run.poses[row["pose_id"]].copy(),
            _intent(row),
            _hf_bracket(classified, points),
            view_factory(calc, hs_ledger),
            run.config,
            hs_ledger,
            task,
            branch_settings=ts_branch.BranchSettings(geometry_noise_A=run.hf.geometry_noise_A),
            settings=settings_for(run, "v5", bound),
            ceiling_eV=ceiling,
        )
    finally:
        hs_ledger.close()
    record = {
        "run": run.name,
        "task": task,
        "hs_status": hs.status,
        "hs_spent": hs.spent,
        "hs_wall_s": time.monotonic() - t0,
        "algorithm": (hs.refine or {}).get("provenance", {}).get("algorithm"),
    }
    found = hs.refine if hs.status == "converged_candidate" else None
    if found is None:
        out.write_text(json.dumps(record) + "\n", encoding="utf-8")
        return record
    frames = list(found["history"])
    for side in hs.sides:
        if side.get("status") == "maximum_bracketed":
            frames += side["points"]
    template = run.poses[row["pose_id"]]
    hc_task = "HC|" + task[3:]
    hc_ledger = Ledger(out.with_suffix(".hc.ledger.jsonl"), arm="HF", budget=run.hf.hc_cap)
    t1 = time.monotonic()
    try:
        starts = branch_starts(
            frames,
            found["metrics"]["q_actual"],
            run.hf.reaction_force_sigma,
            transverse_fmax=settings_for(run, "v5", bound).transverse_fmax_eV_A,
            descent_fmax=run.config.quench_fmax_eV_A,
        )
        result = trace_sides_ef(
            template.copy(),
            found["positions"],
            found["energy_eV"],
            starts,
            found["direction_cartesian"],
            view_factory(calc, hc_ledger),
            run.config,
            hc_ledger,
            hc_task,
            settings=ConnectSettings(
                reaction_force_sigma=run.hf.reaction_force_sigma,
                seed=_seed(run.settings, hc_task),
            ),
        )
    finally:
        hc_ledger.close()
    payload = result.to_dict()
    record |= {
        "saddle_energy_eV": found["energy_eV"],
        "saddle_positions": [list(map(float, r)) for r in found["positions"]],
        "hc_status": result.status,
        "hc_spent": hc_ledger.spent,
        "hc_wall_s": time.monotonic() - t1,
        "fingerprint": payload.get("fingerprint"),
        "sides": [
            {
                "side": sd["side"],
                "status": sd.get("status"),
                "end_energy_eV": (sd.get("descent") or {}).get("energy_eV"),
                "end_positions": (sd.get("descent") or {}).get("positions"),
            }
            for sd in payload.get("sides", [])
        ],
    }
    if result.status == "connected":
        ends = {sd["side"]: sd["descent"]["positions"] for sd in payload["sides"]}
        record["chemical_pair"] = sorted(
            _hf_identity(run, template, ends[k]) for k in ("low", "high")
        )
    out.write_text(json.dumps(record) + "\n", encoding="utf-8")
    return record


def main():
    out_dir = Path(sys.argv[1])
    bound = float(sys.argv[2])
    n_control, seed, limit = map(int, sys.argv[3:6])
    run_dirs = [Path(p) for p in sys.argv[6:]]
    runs = [load(p) for p in run_dirs]
    by_name = {r.name: r for r in runs}
    ref_pair, ref_fingerprint = reference_pair(run_dirs[0].parent.parent)
    picks = [p for p in items(runs, n_control, seed) if p[0] == "failed"]
    picks = picks[:limit] if limit > 0 else picks
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "items.json").write_text(
        json.dumps(
            {
                "items": picks,
                "reference_pair": ref_pair,
                "reference_fingerprint": ref_fingerprint,
            },
            indent=1,
        )
        + "\n"
    )
    calc = polar1("cuda", "float32")
    for i, (group, task, name) in enumerate(picks):
        out = out_dir / f"{group}_{i:02d}.json"
        if out.exists():
            continue
        rec = one(by_name[name], task, bound, calc, out)
        rec["pair_matches_reference"] = rec.get("chemical_pair") == sorted(ref_pair)
        out.write_text(json.dumps(rec) + "\n", encoding="utf-8")
        print(
            json.dumps(
                {"file": out.name}
                | {
                    k: rec.get(k)
                    for k in (
                        "hs_status",
                        "hs_spent",
                        "hc_status",
                        "hc_spent",
                        "saddle_energy_eV",
                        "pair_matches_reference",
                    )
                }
            ),
            flush=True,
        )


if __name__ == "__main__":
    main()
