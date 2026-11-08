"""b04 end-to-end D2 smoke (PLAN_EXPLORE_M5 Task 8).

  python run_b04_d2.py --run DIR --noise-eV X [--hf [--start-cluster]]

--hf freezes rule D2-HF-v1 (docs/PLAN_HF4_D2_2026-10-26.md H4): L4 is E/F search + E/F
connection, and ONE resident float32 model serves L0 qualification, L1/L2, release and
L4 (POLAR-1 weights are float32). Without --hf the original two-model smoke is unchanged.

One process, one GPU: a resident float32 batch model for L1/L2 and a resident float64
POLAR-1 for release and certification. Rerunning the same command resumes. Smoke
values below are NOT the frozen U1-U5 / budget choices.

Resume a run only on the host that wrote it: on NFS with attribute caching (this
project: 6 h, nocto) another host can see stale files, reuse attempt numbers, and
``repair_tail`` could truncate the real journal. The manifest and every journal event
record the host; report.json lists ``hosts``. Per-process load time goes to
``<run>/process_<pid>.json``."""

import argparse
import json
import os
import socket
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT / "src"), str(ROOT / "explore" / "src")]

from ase.io import read  # noqa: E402

from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.d2 import Backends, D2Settings, HFSettings, StartCluster, init_run, run  # noqa: E402
from explore.intents import distance_intents  # noqa: E402
from explore.ledger import Ledger, view_factory  # noqa: E402
from explore.parallel import MACEBatchEvaluator  # noqa: E402
from explore.poses import (  # noqa: E402
    conformer_grid,
    conformers_from_qualification,
    fragments,
    qualify_conformers,
    sample_poses,
)
from prrs.config import SearchConfig  # noqa: E402
from prrs.search import _stamped_with_charge  # noqa: E402

INPUT = ROOT / "benchmarks/b04_sn2_chloride/input.extxyz"
SEED = 20261023
SMOKE = {
    "shells": [2.0],
    "n_samples": 2,
    "steps_per_turn": 6,
    "conformer_cap": 1000,
    "n_far": 2,
}


def settings(noise_eV, hf=False):
    # HF: every HS / HC dispatch holds its full cap (3000), so the smoke's refine
    # account must hold several items at once; 20000 covers the b04 smoke population.
    return D2Settings(
        seed=SEED,
        search_budget=8000,
        certify_budget=20000 if hf else 3000,
        scan_cap=400,
        reverse_cap=100,
        release_cap=300,
        certify_cap=800,
        per_pose_quota=100,
        audit_fraction=0.1,
        noise_eV=noise_eV,
        certify_max_rise_eV=1.0,
        magnitudes_A=tuple(round(0.1 * k, 10) for k in range(1, 18)),
        fmax_eV_A=0.01,
        max_steps_per_point=300,
        max_rise_eV=3.0,
        max_active=4,
        wave_size=16,
    )


# b04 measured inputs (DONE_LOG 2026-10-25): same-branch reaction-force sigma (larger of
# 0.0034 far / 0.0042 near the root), max pairwise RMSD 6.2e-4 A, midpoint tau 1e-3.
HF_B04 = HFSettings(
    reaction_force_sigma=0.0042,
    geometry_noise_A=6.2e-4,
    midpoint_transverse_fmax_eV_A=1e-3,
    hs_cap=3000,
    hc_cap=3000,
)

# --start-cluster (with --hf) freezes D2-HF-v1.1 (docs/PLAN_START_CLUSTER_2026-10-28.md §3):
# radius 2 x the measured geometry noise above, one HS per start cluster in round one.
START_CLUSTER = StartCluster(noise_multiple=2.0, first_round=1)


def build_l0(atoms, config, calc64, run_dir):
    ledger = Ledger(run_dir.parent / f"{run_dir.name}_prepare.jsonl", arm="P", replay=False)
    frags = fragments(atoms, config)
    states = {
        f.key: {"charge": -1 if len(f.atoms) == 1 else 0, "multiplicity": 1} for f in frags
    }
    chosen = []
    for f in frags:
        confs, _ = conformer_grid(
            atoms, f, config, SMOKE["steps_per_turn"], SMOKE["conformer_cap"], SEED
        )
        if len(f.atoms) > 1:
            rows = qualify_conformers(
                atoms, confs, view_factory(calc64, ledger), config, ledger, SEED, states
            )
            confs = conformers_from_qualification(confs, rows)
        chosen.append(confs[0])
    ledger.close()
    poses = sample_poses(
        atoms, frags, chosen[0], chosen[1], SMOKE["shells"], SMOKE["n_samples"], SEED
    )
    out = []
    for p in poses:
        if p.status != "placed":
            continue
        a = atoms.copy()
        a.positions = p.positions
        out.append((p.pose_id, a))
    return out, ledger.spent


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--noise-eV", type=float, required=True)
    ap.add_argument("--hf", action="store_true")
    ap.add_argument(
        "--start-cluster",
        action="store_true",
        help="with --hf: freeze D2-HF-v1.1 (two-round HS dispatch, START_CLUSTER)",
    )
    args = ap.parse_args()
    run_dir = Path(args.run)
    frozen = run_dir / "manifest.json"
    if frozen.exists():
        noise = json.loads(frozen.read_text(encoding="utf-8"))["settings"]["noise_eV"]
        if noise != args.noise_eV:
            print(f"noise_eV: using the frozen {noise} from manifest.json, not {args.noise_eV}")
    config = SearchConfig(charge_sensitive=True, total_charge=-1, multiplicity=1)
    t0 = time.perf_counter()
    calc32 = polar1("cuda", "float32")
    calc64 = calc32 if args.hf else polar1("cuda", "float64")
    load_s = time.perf_counter() - t0
    if not (run_dir / "manifest.json").exists():
        atoms = _stamped_with_charge(read(INPUT), config)
        atoms.info["spin"] = 1
        poses, prepare = build_l0(atoms, config, calc64, run_dir)
        specs = {
            pid: distance_intents(a, pid, config, SMOKE["n_far"], SEED)[0] for pid, a in poses
        }
        identity = {
            "system": "b04",
            "model": "MACE-POLAR-1-M",
            "search_dtype": "float32",
            "certify_dtype": "float32" if args.hf else "float64",
            "l0": SMOKE,
            "prepare_forwards": prepare,
        }
        init_run(
            run_dir,
            poses,
            specs,
            settings(args.noise_eV, args.hf),
            config,
            identity,
            hf=HF_B04 if args.hf else None,
            start_cluster=START_CLUSTER if args.start_cluster else None,
        )
    report = run(run_dir, Backends(MACEBatchEvaluator(calc32), calc64))
    info = {"host": socket.gethostname(), "pid": os.getpid(), "load_s": load_s}
    (run_dir / f"process_{os.getpid()}.json").write_text(json.dumps(info), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
