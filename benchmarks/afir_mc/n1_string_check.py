"""Simplified string + climbing image as the TS locator (explore.afir_string): comparison
on N1 DA paths of PBS 29729 (2026-11-19, user: "简化 String 方法交替做力驱动更新和路径重新
分布，避免节点全挤到反应物盆；之后可以用爬坡像定位鞍点。这样试试？").

Same replay as n1_cineb_check.py (prepare.json copied, SC nodes re-minimized, only the named
two-pair items built, pair / two-pair push, no two-maxima stop); the replayed frames are
saved to OUT/frames_<label>.npz and reused on a rerun (skips the model-heavy setup). End
points: batched E/F minima of the path's first and last frame (fmax --fmax); interior nodes:
the AFIR frames from the early low-energy frame to the end, resampled at equal arc length.
Then afir_string.string_climb (all nodes batched); ceiling = start minimum + 3 eV (the
runner's TS ceiling). Writes OUT/string.json and OUT/wave_0000.json (converged climbing
images as converged_candidate) for benchmarks/afir_mc/afir_validate.py. Not resumable
beyond the cached frames.

    python n1_string_check.py --run DIR --seed S --shell R --sc-poses P --labels L1,L2 --out OUT
"""

import argparse
import json
import shutil
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from run_mc_afir import prepare  # noqa: E402  (sets sys.path)
from run_mc_afir_batch import sc_items  # noqa: E402
from run_d2 import make_atoms, make_config  # noqa: E402
from systems import SYSTEMS  # noqa: E402


from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.afir_batch import minimize, run_paths  # noqa: E402
from explore.afir_string import resample, string_climb  # noqa: E402
from explore.afir_path import PathSettings  # noqa: E402
from explore.ledger import Ledger  # noqa: E402
from explore.parallel import MACEBatchEvaluator  # noqa: E402
from explore.poses import _jsonable, sample_poses  # noqa: E402

FORMING = {"C1-C6": (0, 11), "C4-C7": (3, 12), "C1-C7": (0, 12), "C4-C6": (3, 11)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--system", default="n1")
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--shell", type=float, required=True)
    ap.add_argument("--sc-poses", type=int, required=True)
    ap.add_argument("--max-paths", type=int, default=4096)
    ap.add_argument("--gamma", type=float, default=106.9)
    ap.add_argument("--labels", required=True)
    ap.add_argument("--nodes", type=int, default=20, help="string nodes incl. end points")
    ap.add_argument(
        "--fmax", type=float, default=0.01, help="end points and climbing image, eV/A"
    )
    ap.add_argument("--string-steps", type=int, default=1000)
    ap.add_argument("--ci-steps", type=int, default=2000)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    run, out = Path(args.run), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    shutil.copy(run / "prepare.json", out / "prepare.json")
    labels = args.labels.split(",")

    t0 = time.monotonic()
    system = SYSTEMS[args.system]
    config = make_config(system)
    calc = polar1("cuda", "float32")
    atoms, _ = make_atoms(system, config)
    batch = MACEBatchEvaluator(calc)
    cache = {lb: out / f"frames_{lb.replace('|', '_').replace('+', '_')}.npz" for lb in labels}
    frames_of, replay_info = {}, {}
    if all(p.exists() for p in cache.values()):
        for lb, p in cache.items():
            z = np.load(p)
            frames_of[lb] = (z["positions"], z["energies"])
            replay_info[lb] = json.loads(str(z["info"]))
    else:
        probe = atoms.copy()
        probe.calc = calc
        noise = 2 * float(np.spacing(np.float32(abs(float(probe.get_potential_energy())))))
        frags, confs, _ = prepare(out, system, config, atoms, calc)
        settings = PathSettings(noise_eV=noise, gamma_kJ=args.gamma, early_discard_steps=1001)
        poses = [
            p
            for p in sample_poses(atoms, frags, confs[0], confs[1], [args.shell],
                                  args.max_paths, args.seed)
            if p.status == "placed"
        ]  # fmt: skip
        doubles = []
        singles, _ = sc_items(
            atoms, poses, config, batch, out,
            SimpleNamespace(seed=args.seed, sc_poses=args.sc_poses, sc_coord="pair", sc_double=1),
            doubles,
            double_filter=lambda lb: lb in labels,
        )  # fmt: skip
        by_label = {it["label"]: it for it in singles + doubles}
        chosen = [by_label[lb] for lb in labels]
        led = Ledger(out / "replay.ledger.jsonl", arm="STR", replay=False)
        paths = run_paths(atoms, [it["start"] for it in chosen], [it["coord"] for it in chosen],
                          batch, config, settings, led, "string|replay")  # fmt: skip
        led.close()
        for lb, it, p in zip(labels, chosen, paths):
            fr = p.frames
            info = {"status": p.status, "steps": len(fr), "forwards": p.spent,
                    "bond_change": bool(fr) and fr[-1]["edges"] != it["ref_edges"]}  # fmt: skip
            pos = np.stack([f["positions"] for f in fr])
            en = np.array([f["energy_eV"] for f in fr])
            np.savez(cache[lb], positions=pos, energies=en, info=json.dumps(info))
            frames_of[lb], replay_info[lb] = (pos, en), info
    t_replay = time.monotonic() - t0

    report, wave_paths = [], []
    ledger = Ledger(out / "string.ledger.jsonl", arm="STR", replay=False)
    for lb in labels:
        pos, en = frames_of[lb]
        rec = {"label": lb, "replay": replay_info[lb]}
        if not replay_info[lb]["bond_change"]:
            rec["string"] = "skipped: replay did not change the bond graph"
            report.append(rec)
            continue
        t1 = time.monotonic()
        ends = minimize(atoms, [pos[0], pos[-1]], batch, config, ledger, f"string|{lb}|ends",
                        fmax=args.fmax, max_steps=2000)  # fmt: skip
        low = int(np.argmin(en[: len(en) // 3 + 1]))
        mids = resample(list(pos[low:]), args.nodes)[1:-1]
        nodes = [ends[0].positions] + list(mids) + [ends[1].positions]
        res = string_climb(atoms, nodes, batch, ledger, f"string|{lb}",
                           string_steps=args.string_steps, ci_steps=args.ci_steps,
                           ci_fmax=args.fmax, ceiling_eV=ends[0].energy_eV + 3.0)  # fmt: skip
        ci = res["ci_positions"]
        e0 = ends[0].energy_eV
        rec.update({
            "ends_converged": [bool(m.converged) for m in ends],
            "string_status": res["string_status"], "string_iterations": res["string_iterations"],
            "node_energies_rel_eV": [round(v - e0, 3) for v in res["node_energies_eV"]],
            "ci_index": res["ci_index"], "ci_status": res["ci_status"], "ci_steps": res["ci_steps"],
            "ci_energy_rel_eV": res["ci_energy_eV"] - e0, "ci_fmax_eV_A": res["ci_fmax_eV_A"],
            "reaction_energy_eV": ends[1].energy_eV - e0,
            "ci_forming_A": {k: round(float(np.linalg.norm(ci[j] - ci[i])), 3)
                             for k, (i, j) in FORMING.items()},
            "forwards": res["forwards"] + sum(m.forwards for m in ends),
            "wall_s": time.monotonic() - t1,
        })  # fmt: skip
        report.append(rec)
        print(json.dumps(rec), flush=True)
        if res["ci_status"] == "converged":
            wave_paths.append({
                "k": len(wave_paths), "label": lb, "new_events": [],
                "ts": [{"status": "converged_candidate", "positions": ci,
                        "energy_eV": res["ci_energy_eV"], "link": {}}],
            })  # fmt: skip
    ledger.close()
    spent = sum(r.get("forwards", 0) for r in report)
    (out / "wave_0000.json").write_text(json.dumps(_jsonable(
        {"wave": 0, "paths": wave_paths, "spent_total": spent}), default=str) + "\n")  # fmt: skip
    (out / "string.json").write_text(json.dumps(_jsonable(
        {"args": vars(args), "replay_s": t_replay, "wall_s": time.monotonic() - t0,
         "paths": report}), indent=1, default=str) + "\n")  # fmt: skip
    print(json.dumps({"wall_s": time.monotonic() - t0, "replay_s": t_replay}), flush=True)


if __name__ == "__main__":
    main()
