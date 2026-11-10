"""D2 run on any system of systems.py, real model (docs/PLAN_HF5_CROSS_SYSTEM_2026-11-08.md §3).

  python run_d2.py --system NAME --measure-noise OUT_DIR --seed S --shell R [R ...] \\
      --n-samples N --n-far K
  python run_d2.py --system NAME --run DIR --noise-eV X --seed S --shell R [R ...] \\
      --n-samples N --n-far K --per-pose-quota Q --search-budget B --certify-budget C \\
      --scan-cap ... --reverse-cap ... --release-cap ... --certify-cap ... \\
      --audit-fraction F --max-active M --wave-size W \\
      [--hf [--start-cluster [--same-graph-ridge [--hs-bracket [--sample-parallel]]]]] \\
      [--rotate SEED] [--magnitude-step A] [--certify-max-rise-eV E]

Same modes, rule switches and resume behaviour as run_b05_d2.py, generalised over the
system table (input, charge, fragment states, measured HFSettings / HSBracket). Added for
the HF-5 sensitivity checks (design record §7):
  --rotate SEED        one seeded random rigid rotation of the input about its centroid
                       before L0; the matrix goes into the manifest identity.
  --magnitude-step A   scan magnitude grid A, 2A, ... up to 1.7 A (default 0.1 = the grid
                       of run_b04_d2 / run_b05_d2); frozen through D2Settings.magnitudes_A.
  --certify-max-rise-eV  default MAX_RISE_EV (controller ruling M5-R10, as run_b05_d2);
                       run_b04_d2 smokes used 1.0.

Resume a run only on the host that wrote it (NFS attribute cache 6 h, nocto)."""

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
from l0 import build_l0  # noqa: E402
from systems import SYSTEMS  # noqa: E402

from benchmarks.b04_sn2_chloride.calculator import MODELS, polar1  # noqa: E402
from explore.d2 import Backends, D2Settings, StartCluster, init_run, run  # noqa: E402
from explore.intents import distance_intents  # noqa: E402
from explore.ledger import Ledger  # noqa: E402
from explore.parallel import MACEBatchEvaluator, ScanJob, scan_many  # noqa: E402
from explore.probe import ScanSettings  # noqa: E402
from explore.streams import stream  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402
from prrs.search import _stamped_with_charge  # noqa: E402

# Fixed constants, identical to run_b05_d2.py.
MAX_MAGNITUDE_A = 1.7
FMAX_EV_A = 0.01
MAX_STEPS_PER_POINT = 300
MAX_RISE_EV = 3.0
STEPS_PER_TURN = 6
CONFORMER_CAP = 1000
MIN_QUALIFY_ATOMS = 3
NOISE_INTENTS = 4
NOISE_DELTA_COUNT = 5
NOISE_REPEATS = 5
NOISE_SIGMA_A = 0.02
# run_b04_d2 / run_b05_d2 (docs/PLAN_START_CLUSTER_2026-10-28.md §3)
START_CLUSTER = StartCluster(noise_multiple=2.0, first_round=1)

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


def magnitudes(step):
    n = int(round(MAX_MAGNITUDE_A / step))
    if n < 1 or abs(n * step - MAX_MAGNITUDE_A) > 1e-9:
        raise ValueError(f"--magnitude-step {step} must divide {MAX_MAGNITUDE_A} A")
    return tuple(round(step * k, 10) for k in range(1, n + 1))


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--system", required=True, choices=sorted(SYSTEMS))
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--measure-noise", metavar="OUT_DIR")
    mode.add_argument("--run", metavar="DIR")
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
    ap.add_argument("--rotate", type=int, metavar="SEED")
    ap.add_argument("--magnitude-step", type=float, default=0.1)
    ap.add_argument("--certify-max-rise-eV", type=float, default=MAX_RISE_EV)
    ap.add_argument("--hf", action="store_true", help="freeze D2-HF-v1")
    ap.add_argument("--start-cluster", action="store_true", help="with --hf: D2-HF-v1.1")
    ap.add_argument(
        "--same-graph-ridge", action="store_true", help="with --start-cluster: D2-HF-v1.2"
    )
    ap.add_argument(
        "--hs-bracket",
        action="store_true",
        help="with --same-graph-ridge: D2-HF-v1.3 (the system's measured HSBracket)",
    )
    ap.add_argument(
        "--sample-parallel", action="store_true", help="with --hs-bracket: D2-HF-v1.4"
    )
    args = ap.parse_args(argv)
    needed = L0_ARGS + (RUN_ARGS if args.run else ())
    missing = [f"--{n.replace('_', '-')}" for n in needed if getattr(args, n) is None]
    if missing:
        ap.error("required: " + " ".join(missing))
    system = SYSTEMS[args.system]
    if args.hf and system.hf is None:
        ap.error(f"--hf: HFSettings for {args.system} not measured yet (systems.py)")
    if args.hs_bracket and system.hs_bracket is None:
        ap.error(f"--hs-bracket: HSBracket for {args.system} not measured yet (systems.py)")
    try:
        magnitudes(args.magnitude_step)
    except ValueError as e:
        ap.error(str(e))
    return args


def make_config(system):
    return SearchConfig(
        charge_sensitive=True, total_charge=system.charge, multiplicity=system.multiplicity
    )


def rotation(seed):
    q, r = np.linalg.qr(np.random.default_rng(seed).normal(size=(3, 3)))
    q = q * np.sign(np.diag(r))
    if np.linalg.det(q) < 0:
        q[:, 0] = -q[:, 0]
    return q


def make_atoms(system, config, rotate=None):
    atoms = _stamped_with_charge(read(system.input_path), config)
    atoms.info["spin"] = system.multiplicity
    matrix = None
    if rotate is not None:
        matrix = rotation(rotate)
        c = atoms.positions.mean(axis=0)
        atoms.positions = (atoms.positions - c) @ matrix.T + c
    return atoms, matrix


def l0_for(args, system, atoms, config, calc32, prepare_path):
    return build_l0(
        atoms,
        config,
        calc32,
        prepare_path,
        system.fragment_states,
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
        certify_max_rise_eV=args.certify_max_rise_eV,
        magnitudes_A=magnitudes(args.magnitude_step),
        fmax_eV_A=FMAX_EV_A,
        max_steps_per_point=MAX_STEPS_PER_POINT,
        max_rise_eV=MAX_RISE_EV,
        max_active=args.max_active,
        wave_size=args.wave_size,
    )


def scan_settings(deltas):
    return ScanSettings(tuple(deltas), FMAX_EV_A, MAX_STEPS_PER_POINT, MAX_RISE_EV)


def measure_noise(args, system):
    """Same definition as run_b05_d2.measure_noise (fp32_noise.py, M5-D5)."""
    out = Path(args.measure_noise)
    out.mkdir(parents=True, exist_ok=True)
    config = make_config(system)
    rec = {
        "host": socket.gethostname(),
        "system": system.name,
        "seed": args.seed,
        "sigma_A": NOISE_SIGMA_A,
        "repeats": NOISE_REPEATS,
    }
    t0 = time.perf_counter()
    calc32 = polar1("cuda", "float32")
    rec["load_s"] = time.perf_counter() - t0
    batch = MACEBatchEvaluator(calc32)
    atoms, _ = make_atoms(system, config)
    poses, l0_info = l0_for(args, system, atoms, config, calc32, out / "prepare.jsonl")
    rec["l0"] = l0_info
    if not poses:
        raise RuntimeError("L0 produced no placed pose")
    pose_id, pose = poses[0]
    rec["pose_id"] = pose_id
    specs = distance_intents(pose, pose_id, config, args.n_far, args.seed)[0]
    specs = sorted(
        specs, key=lambda s: float(stream(args.seed, "noise", s.intent.intent_id).random())
    )[:NOISE_INTENTS]
    deltas = magnitudes(args.magnitude_step)[:NOISE_DELTA_COUNT]
    ledger_path = out / "ledger.jsonl"
    scan_jobs = [
        ScanJob(
            pose,
            s.intent,
            scan_settings(s.deltas(deltas)),
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


def run_mode(args, system):
    run_dir = Path(args.run)
    manifest = run_dir / "manifest.json"
    if manifest.exists():
        frozen = json.loads(manifest.read_text(encoding="utf-8"))["settings"]["noise_eV"]
        if frozen != args.noise_eV:
            print(
                f"noise_eV: using the frozen {frozen} from manifest.json, not {args.noise_eV}"
            )
    config = make_config(system)
    t0 = time.perf_counter()
    calc32 = polar1("cuda", "float32")
    load_s = time.perf_counter() - t0
    if not manifest.exists():
        atoms, matrix = make_atoms(system, config, args.rotate)
        prepare_path = run_dir.parent / f"{run_dir.name}_prepare.jsonl"
        poses, l0_info = l0_for(args, system, atoms, config, calc32, prepare_path)
        specs = {
            pid: distance_intents(a, pid, config, args.n_far, args.seed)[0] for pid, a in poses
        }
        path, model_type = MODELS["polar1"]
        identity = {
            "system": system.name,
            "runner": "run_d2.py",
            "model": path,
            "model_type": model_type,
            "search_dtype": "float32",
            "certify_dtype": "float32",
            "args": {k: v for k, v in vars(args).items() if k not in ("run", "measure_noise")},
            "input_rotation": None if matrix is None else matrix.tolist(),
            "l0": l0_info,
        }
        init_run(
            run_dir,
            poses,
            specs,
            d2_settings(args),
            config,
            identity,
            hf=system.hf if args.hf else None,
            start_cluster=START_CLUSTER if args.start_cluster else None,
            same_graph_ridge_last=args.same_graph_ridge,
            hs_bracket=system.hs_bracket if args.hs_bracket else None,
            same_graph_sample=args.sample_parallel,
            hf_parallel=args.sample_parallel,
        )
    report = run(run_dir, Backends(MACEBatchEvaluator(calc32), calc32))
    info = {"host": socket.gethostname(), "pid": os.getpid(), "load_s": load_s}
    (run_dir / f"process_{os.getpid()}.json").write_text(json.dumps(info), encoding="utf-8")
    print(json.dumps(report, indent=2))


def main():
    args = parse_args()
    system = SYSTEMS[args.system]
    if args.measure_noise:
        measure_noise(args, system)
    else:
        run_mode(args, system)


if __name__ == "__main__":
    main()
