"""MC-AFIR with the project's own E/F MIN / TS / connection (docs/PLAN_AFIR_MC_2026-11-14.md;
design docs/NONHESSIAN_AFIR_DESIGN_2026-11-14.md). Real model, float32 POLAR-1.

    python run_mc_afir.py --system NAME --out DIR --seed S --shell R --budget B
        [--max-paths N] [--n-stop 100] [--gamma 106.9] [--max-active 8] [--wave 8]

Preparation (E/F only, no Hessian): each fragment of the system's input complex is
relaxed alone and confirmed by perturbed re-descents (`afir_local.relax_min`); random
mutual orientations and approach positions come from `poses.sample_poses` (geometry only,
one keyed stream per sample, placed at contact ratio --shell).

One task per path k (`path_{k:05d}.json` + its own ledger; an existing JSON is skipped, so
rerunning resumes; an interrupted path is redone and its paid forwards stay in its ledger):
  1. AFIR path (push the two fragments together, gamma -> alpha) on E/F;
  2. MINs at the path's two ends (E/F minimization + perturbation check; GRRM's
     "CONNECTION" from the two end points);
  3. for each energy maximum on the path (at most two): TS bracket by fixed-rho_AB
     relaxations walked back from the top, HF E/F kernel, then E/F connection.
Paths run in waves through one resident batch evaluator (`parallel.run_many`).
Stops when the last --n-stop paths gave no new product (Chem. Rec. MC-AFIR), when
--budget forwards are spent, or after --max-paths paths.

noise_eV = 2 x float32 spacing at the input energy (v1.5 rule). TS / connection numerics
are the v1.5 model-level constants (systems.HF_V15, HS_BRACKET_V15); the branch geometry
noise is the v1.5 fallback (0.0104 A) -- no in-run measurement in phase 1.
Resume only on the host that wrote the directory (NFS attribute cache 6 h).
"""

import argparse
import json
import socket
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT / "src"), str(ROOT / "explore" / "src")]
sys.path.insert(0, str(ROOT / "benchmarks" / "explore_d2"))

from run_d2 import make_atoms, make_config  # noqa: E402
from systems import HF_V15, HS_BRACKET_V15, SYSTEMS  # noqa: E402

from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.afir import WeightedDistance  # noqa: E402
from explore.afir_local import link, relax_min, ts_from_path  # noqa: E402
from explore.afir_path import PathSettings, run_afir_path  # noqa: E402
from explore.d2 import _ef_settings  # noqa: E402
from explore.ledger import Ledger, tally, view_factory  # noqa: E402
from explore.parallel import MACEBatchEvaluator, TaskJob, run_many  # noqa: E402
from explore.poses import Conformer, _jsonable, fragments, sample_poses  # noqa: E402
from explore.streams import stream  # noqa: E402
from explore.ts_branch import BranchSettings  # noqa: E402

PATH_CAP = 20000  # forwards per path task (path + 2 MIN + up to 2 x (TS + connection))
CEILING_RISE_EV = 3.0  # TS search ceiling above the path start (as D2 certify_max_rise)
MAX_MAXIMA = 2


def parse():
    ap = argparse.ArgumentParser()
    ap.add_argument("--system", required=True, choices=sorted(SYSTEMS))
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--shell", type=float, required=True)
    ap.add_argument("--budget", type=int, required=True)
    ap.add_argument("--max-paths", type=int, default=100000)
    ap.add_argument("--n-stop", type=int, default=100)
    ap.add_argument("--gamma", type=float, default=106.9)
    ap.add_argument("--max-active", type=int, default=8)
    ap.add_argument("--wave", type=int, default=8)
    return ap.parse_args()


def identity(atoms, config):
    from prrs.chemistry import chemical_key

    return chemical_key(
        atoms,
        config.bond_scale,
        config.active_atoms,
        charge_sensitive=config.charge_sensitive,
        charge=config.total_charge,
        multiplicity=config.multiplicity,
        parity_tolerance=config.parity_tolerance,
    )["key"]


def prepare(out, system, config, atoms, calc):
    """E/F-relaxed, perturbation-confirmed fragments (no Hessian)."""
    path = out / "prepare.json"
    if path.exists():
        rec = json.loads(path.read_text())
        frags = fragments(atoms, config)
        confs = [
            Conformer(c["id"], f.key, f.atoms, np.asarray(c["positions"]), (), "qualified")
            for c, f in zip(rec["conformers"], frags)
        ]
        return frags, confs, rec
    frags = fragments(atoms, config)
    if len(frags) != 2:
        raise ValueError(f"MC-AFIR phase 1 needs 2 fragments, got {len(frags)}")
    ledger = Ledger(out / "prepare.ledger.jsonl", arm="P", replay=False)
    factory = view_factory(calc, ledger)
    rng = np.random.default_rng(stream(0, "prepare").integers(2**31))
    confs, rows = [], []
    for f in frags:
        sub = atoms[list(f.atoms)]
        state = system.fragment_states[f.key.split("|", 1)[0]]
        sub.info.update(charge=int(state["charge"]), spin=int(state["multiplicity"]))
        res = relax_min(sub, sub.positions, factory, config, ledger, f"prep|{f.key}", rng=rng)
        if res["positions"] is None:
            raise RuntimeError(f"fragment {f.key}: {res['status']}")
        cid = f"f{f.ordinal}|{f.key}:EF"
        confs.append(
            Conformer(cid, f.key, f.atoms, np.asarray(res["positions"]), (), "qualified")
        )
        rows.append({"id": cid, "status": res["status"], "energy_eV": res["energy_eV"],
                     "positions": np.asarray(res["positions"]).tolist()})  # fmt: skip
    ledger.close()
    rec = {"conformers": rows, "forwards": ledger.spent}
    path.write_text(json.dumps(rec, indent=1))
    return frags, confs, rec


def path_task(k, start, atoms, frags, config, settings, ef, branch, sigma, seed, gamma_noise):
    """The whole per-path pipeline as one run_many job body."""

    def go(factory):
        ledger = go.ledger
        a = atoms.copy()
        a.positions = np.asarray(start, float)
        coord = WeightedDistance.build(a.numbers, frags[0].atoms, frags[1].atoms)
        rec = {"k": k}
        res = run_afir_path(a, coord, factory, config, settings, ledger, f"path|{k}")
        frames = [{kk: v for kk, v in f.items()} for f in res.frames]
        rec["path"] = {
            "status": res.status,
            "steps": len(frames),
            "maxima": res.maxima,
            "alpha_eV_A": res.alpha_eV_A,
            "energies_eV": [f["energy_eV"] for f in frames],
            "rho_A": [f["rho_A"] for f in frames],
            "edges_first": frames[0]["edges"] if frames else None,
            "edges_last": frames[-1]["edges"] if frames else None,
            "positions_first": frames[0]["positions"] if frames else None,
            "positions_last": frames[-1]["positions"] if frames else None,
            "positions_maxima": [frames[m]["positions"] for m in res.maxima],
        }
        if res.status in ("budget_exhausted", "domain_exit", "not_reactive") or not frames:
            return rec
        rng = np.random.default_rng(seed + k)
        mins = {}
        for name, idx in (("start", 0), ("end", len(frames) - 1)):
            m = relax_min(a, frames[idx]["positions"], factory, config, ledger,
                          f"min|{k}|{name}", rng=rng)  # fmt: skip
            if m.get("positions") is not None:
                b = a.copy()
                b.positions = m["positions"]
                m["identity"] = identity(b, config)
            mins[name] = m
        rec["mins"] = mins
        ceiling = frames[0]["energy_eV"] + CEILING_RISE_EV
        rec["ts"] = []
        # A path whose bond graph never changed pushed no reaction: the rho branch around
        # its tops is the complex's own minimum, not a barrier (b04 smoke 2026-11-14: 8 of 8
        # such searches walked the whole path and found no bracket). No TS search then.
        if all(f["edges"] == frames[0]["edges"] for f in frames):
            rec["path"]["bond_change"] = False
            return rec
        rec["path"]["bond_change"] = True
        for j, top in enumerate(res.maxima[:MAX_MAXIMA]):
            status, points, br = ts_from_path(
                a, coord, frames, top, factory, config, ledger, f"ts|{k}|{j}",
                ef_settings=ef, branch_settings=branch, ceiling_eV=ceiling,
            )  # fmt: skip
            t = {
                "top": top,
                "status": status,
                "bracket_points": [
                    {"q": p.q_target, "E": p.energy_eV, "rf": p.reaction_force,
                     "converged": p.converged} for p in points
                ],
            }  # fmt: skip
            if status == "converged_candidate":
                found = br.refine
                t |= {"positions": found["positions"], "energy_eV": found["energy_eV"],
                      "fmax_eV_A": found.get("metrics", {}).get("fmax_eV_A")}  # fmt: skip
                conn = link(a, br, factory, config, ledger, f"link|{k}|{j}", ef_settings=ef,
                            sigma=sigma, seed=seed + 1000 * k + j)  # fmt: skip
                c = conn.to_dict()
                ends = {}
                for sd in c["sides"]:
                    pos = (sd.get("descent") or {}).get("positions")
                    if pos is not None:
                        b = a.copy()
                        b.positions = np.asarray(pos)
                        ends[sd["side"]] = {"identity": identity(b, config),
                                            "energy_eV": sd["descent"].get("energy_eV")}  # fmt: skip
                t["link"] = {"status": conn.status, "label": c.get("label"),
                             "fingerprint": c.get("fingerprint"), "ends": ends}  # fmt: skip
            rec["ts"].append(t)
        return rec

    return go


def main():
    args = parse()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    system = SYSTEMS[args.system]
    config = make_config(system)
    t0 = time.monotonic()
    calc = polar1("cuda", "float32")
    load_s = time.monotonic() - t0
    atoms, _ = make_atoms(system, config)
    probe = atoms.copy()
    probe.calc = calc
    e_in = float(probe.get_potential_energy())
    noise = 2 * float(np.spacing(np.float32(abs(e_in))))
    frags, confs, prep = prepare(out, system, config, atoms, calc)
    settings = PathSettings(noise_eV=noise, gamma_kJ=args.gamma)
    ef = _ef_settings(HF_V15, HS_BRACKET_V15)
    branch = BranchSettings(geometry_noise_A=HF_V15.geometry_noise_A)
    sigma = HF_V15.reaction_force_sigma
    poses = [
        p for p in sample_poses(atoms, frags, confs[0], confs[1], [args.shell], args.max_paths,
                                args.seed) if p.status == "placed"
    ]  # fmt: skip
    manifest = {
        "runner": "run_mc_afir.py",
        "system": args.system,
        "args": vars(args),
        "noise_eV": noise,
        "input_energy_eV": e_in,
        "path_settings": settings.__dict__,
        "ef_settings": ef.__dict__,
        "branch_settings": branch.__dict__,
        "path_cap": PATH_CAP,
        "host": socket.gethostname(),
        "load_s": load_s,
        "prepare": {k: v for k, v in prep.items() if k != "conformers"},
    }
    (out / "manifest.json").write_text(json.dumps(_jsonable(manifest), indent=1, default=str))
    batch = MACEBatchEvaluator(calc)
    products, last_new, done = set(), 0, 0

    def spent_total():
        total = prep.get("forwards", 0)
        for f in out.glob("path_*.ledger.jsonl"):
            total += sum(r["forwards"] for r in tally(f).values())
        return total

    def account(rec):
        nonlocal last_new
        react = (rec.get("mins") or {}).get("start", {}).get("identity")
        keys = [(rec.get("mins") or {}).get("end", {}).get("identity")]
        for t in rec.get("ts", []):
            for e in ((t.get("link") or {}).get("ends") or {}).values():
                keys.append(e["identity"])
        new = [x for x in keys if x and x != react and x not in products]
        products.update(new)
        if new:
            last_new = rec["k"]
        return new

    k = 0
    while k < len(poses):
        if spent_total() >= args.budget:
            print("stop: budget", flush=True)
            break
        if k - last_new >= args.n_stop and k > 0:
            print(f"stop: no new product in the last {args.n_stop} paths", flush=True)
            break
        wave = []
        for kk in range(k, min(k + args.wave, len(poses))):
            dest = out / f"path_{kk:05d}.json"
            if dest.exists():
                rec = json.loads(dest.read_text())
                account(rec)
                done += 1
                continue
            ledger = Ledger(dest.with_suffix(".ledger.jsonl"), arm="AFIR", budget=PATH_CAP,
                            replay=False)  # fmt: skip
            go = path_task(kk, poses[kk].positions, atoms, frags, config, settings, ef, branch,
                           sigma, args.seed, noise)  # fmt: skip
            go.ledger = ledger
            wave.append((kk, dest, ledger, TaskJob(go, ledger)))
        k += args.wave
        if not wave:
            continue
        outcomes, _ = run_many([w[3] for w in wave], batch, args.max_active)
        for (kk, dest, ledger, _), o in zip(wave, outcomes):
            ledger.close()
            rec = o.result if o.error is None else {"k": kk, "error": o.error}
            rec |= {"pose_id": poses[kk].pose_id, "forwards": ledger.spent}
            new = account(rec)
            rec["new_products"] = new
            dest.write_text(json.dumps(_jsonable(rec), default=str) + "\n")
            done += 1
            print(kk, (rec.get("path") or {}).get("status"),
                  [t["status"] for t in rec.get("ts", [])],
                  [(t.get("link") or {}).get("status") for t in rec.get("ts", [])],
                  len(new), ledger.spent, flush=True)  # fmt: skip
    summary = {
        "paths": done,
        "products": sorted(products),
        "forwards_total": spent_total(),
        "wall_s": time.monotonic() - t0,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
