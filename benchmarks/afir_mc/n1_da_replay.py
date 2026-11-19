"""Replay of chosen SC-AFIR paths with every frame kept (2026-11-15, user 「一起重放」).

Question (N1, DONE_LOG 2026-11-15 只读分析): paths pushing the Diels-Alder positions reached
rho below the reference DA saddles' rho (1.94-2.18 vs 2.56-2.86 A) at energies above the DA
barrier, without forming a bond. Which atom pairs were pressed together, and where did the
two forming C-C distances stop?

Re-generates the SC items of an existing run (its prepare.json copied in, SC nodes
re-minimized), replays the chosen paths with run_paths (same settings as the runner) and
writes per path: every frame (extxyz with E, rho), the forming C-C distances per frame, and
for the A x B pairs: distance and weight share omega_ij / sum omega at the early low-energy
frame and at the end, sorted by the end weight share. Not resumable (short).

    python n1_da_replay.py --run DIR --seed S --shell R --labels n28|0-12,n24|0-11 --out OUT
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

from ase.io import write  # noqa: E402

from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.afir_batch import run_paths  # noqa: E402
from explore.afir_path import PathSettings  # noqa: E402
from explore.ledger import Ledger  # noqa: E402
from explore.parallel import MACEBatchEvaluator  # noqa: E402
from explore.poses import sample_poses  # noqa: E402

FORMING = {"C1-C6": (0, 11), "C4-C7": (3, 12), "C1-C7": (0, 12), "C4-C6": (3, 11)}


def pair_table(coord, x):
    _, r, w = coord._parts(np.asarray(x, float))
    share = w / w.sum()
    return r, share


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--system", default="n1")
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--shell", type=float, required=True)
    ap.add_argument("--sc-poses", type=int, default=32)
    ap.add_argument("--max-paths", type=int, default=4096)
    ap.add_argument("--gamma", type=float, default=106.9)
    ap.add_argument("--labels", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    run, out = Path(args.run), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    shutil.copy(run / "prepare.json", out / "prepare.json")
    labels = args.labels.split(",")
    recorded = {}
    for f in sorted(run.glob("wave_*.json")):
        for r in json.loads(f.read_text())["paths"]:
            if r["label"] in labels:
                recorded[r["label"]] = r

    t0 = time.monotonic()
    system = SYSTEMS[args.system]
    config = make_config(system)
    calc = polar1("cuda", "float32")
    atoms, _ = make_atoms(system, config)
    probe = atoms.copy()
    probe.calc = calc
    noise = 2 * float(np.spacing(np.float32(abs(float(probe.get_potential_energy())))))
    frags, confs, _ = prepare(out, system, config, atoms, calc)
    base = PathSettings(noise_eV=noise, gamma_kJ=args.gamma)
    settings = PathSettings(
        noise_eV=noise,
        gamma_kJ=args.gamma,
        early_discard_steps=base.max_steps + 1,
        stop_after_maxima=2,
    )
    poses = [
        p
        for p in sample_poses(
            atoms, frags, confs[0], confs[1], [args.shell], args.max_paths, args.seed
        )
        if p.status == "placed"
    ]
    batch = MACEBatchEvaluator(calc)
    items, _ = sc_items(
        atoms,
        poses,
        config,
        batch,
        out,
        SimpleNamespace(seed=args.seed, sc_poses=args.sc_poses, sc_coord="rho"),
    )
    by_label = {it["label"]: it for it in items}
    chosen = [by_label[lb] for lb in labels]
    ledger = Ledger(out / "replay.ledger.jsonl", arm="AFIR", replay=False)
    paths = run_paths(
        atoms,
        [it["start"] for it in chosen],
        [it["coord"] for it in chosen],
        batch,
        config,
        settings,
        ledger,
        "replay|paths",
    )
    ledger.close()
    sym = atoms.get_chemical_symbols()
    report = []
    for it, p in zip(chosen, paths):
        rec = recorded[it["label"]]
        coord = it["coord"]
        fr = p.frames
        E = np.array([f["energy_eV"] for f in fr])
        third = len(fr) // 3 + 1
        low = int(np.argmin(E[:third]))
        imgs = []
        for f in fr:
            a = atoms.copy()
            a.positions = np.asarray(f["positions"], float)
            a.info.update(energy_eV=f["energy_eV"], rho_A=f["rho_A"])
            imgs.append(a)
        tag = it["label"].replace("|", "_")
        write(out / f"{tag}.extxyz", imgs)
        r_low, s_low = pair_table(coord, fr[low]["positions"])
        r_end, s_end = pair_table(coord, fr[-1]["positions"])
        rows = []
        for ia, i in enumerate(coord.a):
            for jb, j in enumerate(coord.b):
                rows.append(
                    {
                        "pair": f"{sym[i]}{i}-{sym[j]}{j}",
                        "r_low": round(float(r_low[ia, jb]), 3),
                        "r_end": round(float(r_end[ia, jb]), 3),
                        "share_low": round(float(s_low[ia, jb]), 3),
                        "share_end": round(float(s_end[ia, jb]), 3),
                    }
                )
        rows.sort(key=lambda d: -d["share_end"])
        forming = {
            k: [
                round(
                    float(np.linalg.norm(np.subtract(f["positions"][i], f["positions"][j]))), 3
                )
                for f in fr
            ]
            for k, (i, j) in FORMING.items()
        }
        report.append(
            {
                "label": it["label"],
                "a": list(coord.a),
                "b": list(coord.b),
                "status": p.status,
                "recorded_status": rec["status"],
                "steps": len(fr),
                "recorded_steps": rec["steps"],
                "rho_end": round(fr[-1]["rho_A"], 3),
                "recorded_rho_end": round(rec["rho_A"][-1], 3),
                "low_frame": low,
                "E_low_to_end_eV": round(float(E[-1] - E[low]), 3),
                "rho_low": round(fr[low]["rho_A"], 3),
                "forming_low": {k: v[low] for k, v in forming.items()},
                "forming_end": {k: v[-1] for k, v in forming.items()},
                "forming_min": {k: min(v[low:]) for k, v in forming.items()},
                "pairs_by_end_share": rows,
            }
        )
    (out / "replay.json").write_text(json.dumps(report, indent=1) + "\n")
    for r in report:
        print(json.dumps({k: v for k, v in r.items() if k != "pairs_by_end_share"}), flush=True)
        for d in r["pairs_by_end_share"][:10]:
            print("   ", d, flush=True)
    print(json.dumps({"forwards": ledger.spent, "wall_s": time.monotonic() - t0}), flush=True)


if __name__ == "__main__":
    main()
