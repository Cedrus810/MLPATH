"""b05 (tert-butyl chloride + OH-, 16 atoms, charge -1, singlet) D2 run, real model.

  python run_b05_d2.py --measure-noise OUT_DIR --seed S --shell R [R ...] --n-samples N --n-far K
  python run_b05_d2.py --run DIR --noise-eV X --seed S --shell R [R ...] --n-samples N \\
      --n-far K --per-pose-quota Q --search-budget B --certify-budget C --scan-cap ... \\
      --reverse-cap ... --release-cap ... --certify-cap ... --audit-fraction F \\
      --max-active M --wave-size W

--measure-noise: build L0, take the first placed pose and the first 4 of its single
distance intents (keyed order), scan each with the resident float32 batch model, then
repeat each intent's last converged point 5 times from coordinates perturbed by
N(0, 0.02 A) and relax back to the same q (definition of fp32_noise.py, M5-D5). Writes
OUT_DIR/noise.json (``noise_eV_rule`` = max over intents of the peak-to-peak spread) and
OUT_DIR/ledger.jsonl. Nothing is frozen by it.

--run: one process, one GPU: one resident float32 model (POLAR-1 weights are float32) for everything:
batched for L1/L2, counted views of the same model for
POLAR-1 for release and certification. init_run only if manifest.json is absent;
rerunning the same command resumes.

Resume a run only on the host that wrote it: on NFS with attribute caching (this
project: 6 h, nocto) another host can see stale files, reuse attempt numbers, and
``repair_tail`` could truncate the real journal. Per-process load time goes to
``<run>/process_<pid>.json``.

Every run parameter is an explicit CLI argument with no default. The values given are
NOT frozen U1-U5 / budget choices; the fixed constants below are the b04-smoke values."""

import argparse
import json
import os
import socket
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT / "src"), str(ROOT / "explore" / "src")]
sys.path.insert(0, str(Path(__file__).resolve().parent))

from ase.io import read  # noqa: E402

from benchmarks.b04_sn2_chloride.calculator import MODELS, polar1  # noqa: E402
from explore.d2 import (  # noqa: E402
    Backends,
    D2Settings,
    HFSettings,
    HSBracket,
    StartCluster,
    init_run,
    run,
)
from explore.intents import distance_intents  # noqa: E402
from explore.ledger import Ledger  # noqa: E402
from explore.parallel import MACEBatchEvaluator, ScanJob, scan_many  # noqa: E402
from explore.probe import ScanSettings  # noqa: E402
from explore.streams import stream  # noqa: E402
from l0 import build_l0  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402
from prrs.search import _stamped_with_charge  # noqa: E402

INPUT = ROOT / "benchmarks/b05_tertiary_halide/input.extxyz"
# Fixed constants (b04-smoke values; not frozen choices).
MAGNITUDES_A = tuple(round(0.1 * k, 10) for k in range(1, 18))
FMAX_EV_A = 0.01
MAX_STEPS_PER_POINT = 300
MAX_RISE_EV = 3.0
CERTIFY_MAX_RISE_EV = MAX_RISE_EV  # controller ruling M5-R10: not below the scan ceiling
STEPS_PER_TURN = 6
CONFORMER_CAP = 1000
MIN_QUALIFY_ATOMS = 3
# Formula strings as produced by prrs.chemistry._formula (ase Hill order).
FRAGMENT_STATES = {
    "C4H9Cl": {"charge": 0, "multiplicity": 1},
    "HO": {"charge": -1, "multiplicity": 1},
}
NOISE_INTENTS = 4
NOISE_DELTA_MAGNITUDES = MAGNITUDES_A[:5]
NOISE_REPEATS = 5
NOISE_SIGMA_A = 0.02

# --hf freezes rule D2-HF-v1 (docs/PLAN_HF4_D2_2026-10-26.md): L4 = E/F search + E/F
# connection. b05 measured inputs (DONE_LOG 2026-10-25 J0 / near-root): same-branch
# reaction-force sigma 0.0025 eV/A, max pairwise RMSD 7e-3 A, midpoint tau 1e-3.
HF_B05 = HFSettings(
    reaction_force_sigma=0.0025,
    geometry_noise_A=7e-3,
    midpoint_transverse_fmax_eV_A=1e-3,
    hs_cap=3000,
    hc_cap=3000,
)

# --start-cluster (with --hf) freezes D2-HF-v1.1 (docs/PLAN_START_CLUSTER_2026-10-28.md §3):
# radius 2 x the measured geometry noise above, one HS per start cluster in round one.
START_CLUSTER = StartCluster(noise_multiple=2.0, first_round=1)

# --same-graph-ridge (with --start-cluster) freezes D2-HF-v1.2
# (docs/PLAN_HF_SAME_GRAPH_RIDGE_2026-11-07.md); adding --hs-bracket freezes D2-HF-v1.3
# (docs/PLAN_HF_D2_V13_2026-11-07.md): HS bracket rule HF-EF-bracket-v7 with the b05
# coarse threshold and sign bound (user 2026-11-02: 2e-3, B = 0.025).
HS_BRACKET_B05 = HSBracket(
    midpoint_coarse_transverse_fmax_eV_A=2e-3, midpoint_sign_bound_eV_A=0.025
)

L0_ARGS = ("seed", "shell", "n_samples", "n_far")
RUN_ARGS = (
    "noise_eV",
    "per_pose_quota",
    "search_budget",
    "certify_budget",
    "scan_cap",
    "reverse_cap",
    "release_cap",
    "certify_cap",
    "audit_fraction",
    "max_active",
    "wave_size",
)


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--measure-noise", metavar="OUT_DIR")
    mode.add_argument("--run", metavar="DIR")
    # required per mode below, so that --measure-noise need not pass run-only values
    ap.add_argument("--seed", type=int)
    ap.add_argument("--shell", type=float, nargs="+")
    ap.add_argument("--n-samples", type=int)
    ap.add_argument("--n-far", type=int)
    ap.add_argument("--noise-eV", type=float)
    ap.add_argument("--per-pose-quota", type=int)
    ap.add_argument("--search-budget", type=int)
    ap.add_argument("--certify-budget", type=int)
    ap.add_argument("--scan-cap", type=int)
    ap.add_argument("--reverse-cap", type=int)
    ap.add_argument("--release-cap", type=int)
    ap.add_argument("--certify-cap", type=int)
    ap.add_argument("--audit-fraction", type=float)
    ap.add_argument("--max-active", type=int)
    ap.add_argument("--wave-size", type=int)
    ap.add_argument("--hf", action="store_true")
    ap.add_argument(
        "--start-cluster",
        action="store_true",
        help="with --hf: freeze D2-HF-v1.1 (two-round HS dispatch, START_CLUSTER)",
    )
    ap.add_argument(
        "--same-graph-ridge",
        action="store_true",
        help="with --start-cluster: freeze D2-HF-v1.2 (same-graph ridges dispatched last)",
    )
    ap.add_argument(
        "--hs-bracket",
        action="store_true",
        help="with --same-graph-ridge: freeze D2-HF-v1.3 (HS bracket v7, HS_BRACKET_B05)",
    )
    ap.add_argument(
        "--sample-parallel",
        action="store_true",
        help="with --hs-bracket: freeze D2-HF-v1.4 (same-graph ridge sample at "
        "--audit-fraction, HS / HC in parallel waves of --wave-size, --max-active)",
    )
    args = ap.parse_args(argv)
    needed = L0_ARGS + (RUN_ARGS if args.run else ())
    missing = [f"--{n.replace('_', '-')}" for n in needed if getattr(args, n) is None]
    if missing:
        ap.error("required: " + " ".join(missing))
    return args


def make_config():
    return SearchConfig(charge_sensitive=True, total_charge=-1, multiplicity=1)


def make_atoms(config):
    atoms = _stamped_with_charge(read(INPUT), config)
    atoms.info["spin"] = 1
    return atoms


def l0_for(args, atoms, config, calc32, prepare_path):
    return build_l0(
        atoms,
        config,
        calc32,
        prepare_path,
        FRAGMENT_STATES,
        args.shell,
        args.n_samples,
        args.seed,
        STEPS_PER_TURN,
        CONFORMER_CAP,
        MIN_QUALIFY_ATOMS,
    )


def d2_settings(args):
    return D2Settings(
        seed=args.seed,
        search_budget=args.search_budget,
        certify_budget=args.certify_budget,
        scan_cap=args.scan_cap,
        reverse_cap=args.reverse_cap,
        release_cap=args.release_cap,
        certify_cap=args.certify_cap,
        per_pose_quota=args.per_pose_quota,
        audit_fraction=args.audit_fraction,
        noise_eV=args.noise_eV,
        certify_max_rise_eV=CERTIFY_MAX_RISE_EV,
        magnitudes_A=MAGNITUDES_A,
        fmax_eV_A=FMAX_EV_A,
        max_steps_per_point=MAX_STEPS_PER_POINT,
        max_rise_eV=MAX_RISE_EV,
        max_active=args.max_active,
        wave_size=args.wave_size,
    )


def scan_settings(deltas):
    return ScanSettings(tuple(deltas), FMAX_EV_A, MAX_STEPS_PER_POINT, MAX_RISE_EV)


def measure_noise(args):
    out = Path(args.measure_noise)
    out.mkdir(parents=True, exist_ok=True)
    config = make_config()
    rec = {
        "host": socket.gethostname(),
        "system": "b05",
        "seed": args.seed,
        "sigma_A": NOISE_SIGMA_A,
        "repeats": NOISE_REPEATS,
    }
    t0 = time.perf_counter()
    calc32 = polar1("cuda", "float32")
    rec["load_s"] = time.perf_counter() - t0
    batch = MACEBatchEvaluator(calc32)
    atoms = make_atoms(config)
    poses, l0_info = l0_for(args, atoms, config, calc32, out / "prepare.jsonl")
    rec["l0"] = l0_info
    if not poses:
        raise RuntimeError("L0 produced no placed pose")
    pose_id, pose = poses[0]
    rec["pose_id"] = pose_id
    specs = distance_intents(pose, pose_id, config, args.n_far, args.seed)[0]
    specs = sorted(
        specs, key=lambda s: float(stream(args.seed, "noise", s.intent.intent_id).random())
    )[:NOISE_INTENTS]
    ledger_path = out / "ledger.jsonl"
    scan_jobs = [
        ScanJob(
            pose,
            s.intent,
            scan_settings(s.deltas(NOISE_DELTA_MAGNITUDES)),
            config,
            Ledger(ledger_path, "fp32-scan", replay=False),
            f"scan|{s.intent.intent_id}",
        )
        for s in specs
    ]
    t0 = time.perf_counter()
    scans, _ = scan_many(scan_jobs, batch, len(scan_jobs))
    rec["scan_wall_s"] = time.perf_counter() - t0
    rng = np.random.default_rng(args.seed)
    jobs, owners, rows = [], [], []
    for spec, sjob, sc in zip(specs, scan_jobs, scans):
        row = {
            "intent_id": spec.intent.intent_id,
            "sign": spec.sign,
            "scan_forwards": sjob.ledger.spent,
            "scan_status": sc.result.status if sc.result else None,
            "scan_error": sc.error,
            "energies_eV": [],
            "failed": [],
            "forwards": 0,
        }
        rows.append(row)
        converged = [p for p in (sc.result.points if sc.result else []) if p.converged]
        if not converged:
            row["failed"].append("no converged scan point")
            continue
        base = converged[-1]
        row["q_target"] = base.q_target
        row["base_eV"] = base.energy_eV
        for r in range(NOISE_REPEATS):
            a = pose.copy()
            a.positions = base.positions + rng.normal(0.0, NOISE_SIGMA_A, base.positions.shape)
            back = scan_settings((0.0, base.q_target - spec.intent.value(a.positions)))
            jobs.append(
                ScanJob(
                    a,
                    spec.intent,
                    back,
                    config,
                    Ledger(ledger_path, "fp32-noise", replay=False),
                    f"noise|{spec.intent.intent_id}|{r}",
                )
            )
            owners.append(row)
    t0 = time.perf_counter()
    # 4-way, as fp32_noise.py; --max-active belongs to --run mode and is unset here.
    outcomes, _ = scan_many(jobs, batch, 4) if jobs else ([], None)
    rec["noise_wall_s"] = time.perf_counter() - t0
    for row, job, o in zip(owners, jobs, outcomes):
        row["forwards"] += job.ledger.spent
        last = None if o.result is None or not o.result.points else o.result.points[-1]
        if last is not None and last.converged:
            row["energies_eV"].append(last.energy_eV)
        else:
            row["failed"].append(o.error or (o.result.status if o.result else None))
    for row in rows:
        e = np.asarray(row["energies_eV"])
        row["ptp_eV"] = float(np.ptp(e)) if len(e) >= 2 else None
        row["std_eV"] = float(e.std(ddof=1)) if len(e) >= 2 else None
        if "base_eV" in row:
            row["fp32_quantum_eV"] = float(np.spacing(np.float32(abs(row["base_eV"]))))
    rec["intents"] = rows
    ptps = [r["ptp_eV"] for r in rows if r["ptp_eV"] is not None]
    rec["noise_eV_rule"] = max(ptps) if ptps else None
    bases = [abs(r["base_eV"]) for r in rows if "base_eV" in r]
    rec["fp32_quantum_eV"] = float(np.spacing(np.float32(max(bases)))) if bases else None
    (out / "noise.json").write_text(json.dumps(rec, indent=2), encoding="utf-8")
    print(json.dumps({k: rec[k] for k in ("noise_eV_rule", "fp32_quantum_eV")}))


def run_mode(args):
    run_dir = Path(args.run)
    manifest = run_dir / "manifest.json"
    if manifest.exists():
        frozen = json.loads(manifest.read_text(encoding="utf-8"))["settings"]["noise_eV"]
        if frozen != args.noise_eV:
            print(
                f"noise_eV: using the frozen {frozen} from manifest.json, not {args.noise_eV}"
            )
    config = make_config()
    t0 = time.perf_counter()
    calc32 = polar1("cuda", "float32")
    load_s = time.perf_counter() - t0
    if not manifest.exists():
        atoms = make_atoms(config)
        prepare_path = run_dir.parent / f"{run_dir.name}_prepare.jsonl"
        poses, l0_info = l0_for(args, atoms, config, calc32, prepare_path)
        specs = {
            pid: distance_intents(a, pid, config, args.n_far, args.seed)[0] for pid, a in poses
        }
        path, model_type = MODELS["polar1"]
        identity = {
            "system": "b05",
            "model": path,
            "model_type": model_type,
            "search_dtype": "float32",
            "certify_dtype": "float32",
            "args": {k: v for k, v in vars(args).items() if k not in ("run", "measure_noise")},
            "l0": l0_info,
        }
        init_run(
            run_dir,
            poses,
            specs,
            d2_settings(args),
            config,
            identity,
            hf=HF_B05 if args.hf else None,
            start_cluster=START_CLUSTER if args.start_cluster else None,
            same_graph_ridge_last=args.same_graph_ridge,
            hs_bracket=HS_BRACKET_B05 if args.hs_bracket else None,
            same_graph_sample=args.sample_parallel,
            hf_parallel=args.sample_parallel,
        )
    report = run(run_dir, Backends(MACEBatchEvaluator(calc32), calc32))
    info = {"host": socket.gethostname(), "pid": os.getpid(), "load_s": load_s}
    (run_dir / f"process_{os.getpid()}.json").write_text(json.dumps(info), encoding="utf-8")
    print(json.dumps(report, indent=2))


def main():
    args = parse_args()
    if args.measure_noise:
        measure_noise(args)
    else:
        run_mode(args)


if __name__ == "__main__":
    main()
