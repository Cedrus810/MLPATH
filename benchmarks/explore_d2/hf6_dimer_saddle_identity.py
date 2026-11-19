"""Identity and connection of the order-1 Dimer saddles that match no reference (HF-6, handoff
docs/HANDOFF_HF5_HF6_2026-11-12.md §5 "Dimer 找到的参照外一阶鞍点"). Validation side only:
nothing here feeds back into any search arm, and the forward count is its own ledger.

Inputs are the candidates the HF-6 summaries already judged order 1 and unmatched:
  local  results/2026-11-10/hf6_local_summary_v3.json (rows -> hf6_local/<sys>/<i>_dimer.json)
  e2e    results/2026-11-11/hf6_e2e_summary_sA_v3.json, hf6_e2e_summary_sB_v2.json
         (candidates_detail is in the order of the converged attempt files, sorted by name)
Per system they are pooled and grouped into distinct saddles with the summaries' own rule
(|dE| <= basin_energy_eV and symmetric RMSD <= basin_rmsd_A under `match_group`). For one
representative per group:
  1. order again by `confirm_minimum` (the validation config, analytic cross-check 0.1);
  2. if order 1: an E/F-only descent on both signs of the unstable mode -- the step of
     `prrs.runner.follow_unstable_mode` (largest atom displaced `irc_step_A`), then ASE BFGS
     to fmax <= DESCENT_FMAX within DESCENT_STEPS. The ends are labelled by bond graph only:
     no Hessian minimum check, no relay. (`descend_saddle` was tried first, 2026-11-12 smoke:
     on N1 group 0 both sides stopped at its own gates, `internal_mode_recheck_disagrees` /
     `fmax_not_reached_after_polish`, so it named no end.) The ends are compared with
     `prrs.state.classify` and each with the reactant pose graph of the D2 run;
  3. nearest D2 HC-connected event of the same system (sA and sB runs), same RMSD rule.

One JSON per group, skipped if present (rerun the same command to resume).

    python hf6_dimer_saddle_identity.py OUT_DIR [--systems b05 n1] [--limit K]
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from hf_event_validate import _jsonable, match_group, validation_config  # noqa: E402  (sets sys.path)

from ase import Atoms  # noqa: E402

from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.ledger import Ledger, view_factory  # noqa: E402
from prrs.chemistry import symmetric_rmsd  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402
from prrs.runner import confirm_minimum  # noqa: E402
from ase.optimize import BFGS  # noqa: E402
from prrs.search import _stamped_with_charge  # noqa: E402
from prrs.state import classify, encode  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "benchmarks/explore_d2/results"
LOCAL = RES / "2026-11-10/hf6_local_summary_v3.json"
E2E = {
    "sA": (RES / "2026-11-11/hf6_e2e_summary_sA_v3.json", RES / "2026-11-11/hf6_e2e_db"),
    "sB": (RES / "2026-11-11/hf6_e2e_summary_sB_v2.json", RES / "2026-11-11/hf6_e2e_db_sB"),
}
D2RUNS = {
    "b04": [RES / "2026-11-09/b04_hf5_v14_sA", RES / "2026-11-09/b04_hf5_v14_sB"],
    "b05": [RES / "2026-11-08/b05_hf_v14_sA", RES / "2026-11-08/b05_hf_v14_sB"],
    "n1": [RES / "2026-11-09/n1_hf5_v12_sA_n0020", RES / "2026-11-11/n1_hf5_v12_sB_n0020"],
}
SEED = 20261112
DESCENT_FMAX = 0.01  # eV/A, the HF-6 candidate criterion
DESCENT_STEPS = 1000


def candidates(name):
    """(source label, positions, energy from the summary) of unmatched order-1 Dimer saddles."""
    out = []
    local = json.loads(LOCAL.read_text())
    for row in local[name]["rows"]:
        c = row["dimer"]
        if c.get("order") == 1 and not c.get("match"):
            rec = json.loads(
                (RES / "2026-11-10/hf6_local" / name / f"{row['item']:03d}_dimer.json").read_text()
            )
            out.append((f"local:{row['item']:03d}", rec["positions"], c["energy_eV"]))
    for tag, (summary, dbdir) in E2E.items():
        detail = json.loads(summary.read_text())[name]["arms"]["db"]["candidates_detail"]
        files = []
        for f in sorted((dbdir / name).glob("[0-9][0-9][0-9]_[0-9][0-9][0-9].json")):
            rec = json.loads(f.read_text())
            if rec["status"] == "converged":
                files.append((f.stem, rec))
        if len(files) != len(detail):
            raise RuntimeError(f"{name} {tag}: {len(files)} converged files, {len(detail)} rows")
        for (stem, rec), c in zip(files, detail):
            if abs(rec["cumulative_forwards"] - c["cost"]) > 0:
                raise RuntimeError(f"{name} {tag} {stem}: cost mismatch")
            if c.get("order") == 1 and not c.get("match"):
                out.append((f"e2e_{tag}:{stem}", rec["positions"], c["energy_eV"]))
    return out


def d2_events(run_dir):
    out = []
    with open(run_dir / "journal.jsonl", encoding="utf-8") as fh:
        for line in fh:
            if '"HC|' not in line or '"settle"' not in line:
                continue
            ev = json.loads(line)
            if ev.get("event") == "settle" and ev["status"] == "connected":
                p = ev["payload"]
                out.append(
                    {
                        "run": run_dir.name,
                        "task": ev["task"],
                        "energy_eV": p["saddle_energy_eV"],
                        "broken": p["fingerprint"]["broken"],
                        "formed": p["fingerprint"]["formed"],
                        "positions": p["saddle_positions"],
                    }
                )
    return out


def main():
    args = sys.argv[1:]
    out = Path(args[0])
    names = ["b04", "b05", "n1"]
    if "--systems" in args:
        k = args.index("--systems") + 1
        names = []
        while k < len(args) and not args[k].startswith("--"):
            names.append(args[k])
            k += 1
    limit = int(args[args.index("--limit") + 1]) if "--limit" in args else None
    out.mkdir(parents=True, exist_ok=True)
    calc = polar1("cuda", "float32")
    ledger = Ledger(out / "identity.ledger.jsonl", arm="V", replay=False)
    factory = view_factory(calc, ledger)
    t0 = time.monotonic()
    for name in names:
        run = D2RUNS[name][0]
        manifest = json.loads((run / "manifest.json").read_text())
        config = validation_config(SearchConfig.from_dict(manifest["config"]))
        pose = json.loads((run / "poses.jsonl").read_text().splitlines()[0])
        numbers = pose["numbers"]

        def atoms_at(positions):
            a = _stamped_with_charge(
                Atoms(numbers=numbers, positions=np.asarray(positions, float)), config
            )
            a.info["spin"] = config.multiplicity
            a.calc = factory()
            return a

        pose_edges = encode(
            Atoms(numbers=numbers, positions=pose["positions"]),
            config.bond_scale,
            active=config.active_atoms,
        ).edges

        def rmsd(a, b):
            return float(symmetric_rmsd(a, b, match_group(a, run, config), config.active_atoms))

        groups = []
        for label, pos, e in candidates(name):
            a = Atoms(numbers=numbers, positions=np.asarray(pos, float))
            for g in groups:
                if abs(e - g["energy_eV"]) <= config.basin_energy_eV and (
                    rmsd(a, g["atoms"]) <= config.basin_rmsd_A
                ):
                    g["members"].append(label)
                    break
            else:
                groups.append({"atoms": a, "energy_eV": e, "members": [label]})
        events = [ev for r in D2RUNS[name] for ev in d2_events(r)]
        print(name, "candidates grouped:", len(groups), "D2 events:", len(events), flush=True)
        for gi, g in enumerate(groups[:limit]):
            dest = out / f"{name}_{gi:03d}.json"
            if dest.exists():
                continue
            t1 = time.monotonic()
            spent0 = ledger.spent
            ts = atoms_at(g["atoms"].positions)
            rec = {
                "system": name,
                "group": gi,
                "members": g["members"],
                "summary_energy_eV": g["energy_eV"],
                "energy_eV": float(ts.get_potential_energy()),
            }
            near = []
            for ev in events:
                b = Atoms(numbers=numbers, positions=np.asarray(ev["positions"], float))
                near.append(
                    (rmsd(g["atoms"], b), ev["energy_eV"] - rec["energy_eV"], ev)
                )
            if near:
                d, de, ev = min(near, key=lambda t: t[0])
                rec["nearest_d2_event"] = {
                    "rmsd_A": d,
                    "dE_eV": de,
                    "same": bool(d <= config.basin_rmsd_A and abs(de) <= config.basin_energy_eV),
                } | {k: ev[k] for k in ("run", "task", "broken", "formed")}
            with ledger.segment(task=f"{name}|{gi}", stage="hessian", attempt=0) as seg:
                _, chk = confirm_minimum(ts, factory, config, SEED)
                seg.outcome = str(chk.get("saddle_order", chk.get("reason")))
            rec["order"] = chk.get("saddle_order")
            rec["order_reason"] = chk.get("reason")
            rec["imaginary_icm"] = _jsonable(chk.get("imaginary_wavenumbers_icm"))
            mode = chk.get("unstable_mode")
            if rec["order"] == 1 and mode is not None:
                ends = {}
                sides = []
                with ledger.segment(task=f"{name}|{gi}", stage="descent", attempt=0) as seg:
                    m = np.asarray(mode, float).reshape(-1, 3)
                    m = m / np.linalg.norm(m, axis=1).max()
                    for sign in (1, -1):
                        end = atoms_at(ts.positions + sign * config.irc_step_A * m)
                        opt = BFGS(end, logfile=None)
                        converged = bool(opt.run(fmax=DESCENT_FMAX, steps=DESCENT_STEPS))
                        ge = encode(end, config.bond_scale, active=config.active_atoms)
                        side = {
                            "sign": sign,
                            "converged": converged,
                            "steps": int(opt.nsteps),
                            "fmax_eV_A": float(np.linalg.norm(end.get_forces(), axis=1).max()),
                            "energy_eV": float(end.get_potential_energy()),
                            "vs_pose_broken": sorted(list(e) for e in pose_edges - ge.edges),
                            "vs_pose_formed": sorted(list(e) for e in ge.edges - pose_edges),
                            "components": int(ge.components),
                        }
                        if converged:
                            ends[sign] = end
                        sides.append(side)
                    seg.outcome = "done"
                rec["sides"] = sides
                if len(ends) == 2:
                    lab, fp = classify(ends[1], ends[-1], config)
                    rec["ends_label"] = lab
                    rec["ends_fingerprint"] = _jsonable(fp)
            rec["forwards"] = ledger.spent - spent0
            rec["wall_s"] = time.monotonic() - t1
            dest.write_text(json.dumps(rec, indent=1, default=str) + "\n", encoding="utf-8")
            print(
                name,
                gi,
                len(g["members"]),
                rec["order"],
                rec.get("ends_label"),
                rec.get("ends_fingerprint", {}).get("broken"),
                rec.get("ends_fingerprint", {}).get("formed"),
                rec["forwards"],
                round(rec["wall_s"], 1),
                flush=True,
            )
    ledger.close()
    (out / "done.json").write_text(
        json.dumps({"forwards": ledger.spent, "wall_s": time.monotonic() - t0}) + "\n"
    )


if __name__ == "__main__":
    main()
