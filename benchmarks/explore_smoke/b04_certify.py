"""Smoke: C2 certification of the 10-18 b04 scan peak (Cl- + CH3Cl, POLAR-1, charge -1).

Initial guess = the q_target == 2.28 frame of results/2026-10-18/run/fine_fwd.extxyz.
ceiling_eV / reference_eV are explicit smoke parameters (guess energy + 1.0 / guess
energy), not defaults. The reference TS is only compared after the fact (read-only).
Not an acceptance run: a result that differs from expectation is reported as is.
"""

import json
import platform
import socket
import sys
import time
from pathlib import Path

import numpy as np
from ase.io import read, write

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "explore" / "src"), str(ROOT)]

from benchmarks.b04_sn2_chloride.calculator import MODELS, polar1  # noqa: E402
from explore.certify import certify  # noqa: E402
from explore.ledger import Ledger, counted_factory  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402
from prrs.search import _stamped_with_charge  # noqa: E402

GUESS_FILE = ROOT / "benchmarks/explore_smoke/results/2026-10-18/run/fine_fwd.extxyz"
GUESS_Q = 2.28
TS = ROOT / "runs/b04_probe2_openceiling/ts_candidates/ts0000.extxyz"
C, CL_LEAVE, CL_NUC = 0, 1, 5
CEILING_RISE_EV = 1.0
BUDGET = 20000
DEVICE = "cuda"
DTYPE = "float64"
TASK = "certify"
SEED = 0


def find_guess_frame(path=GUESS_FILE, q=GUESS_Q):
    for frame in read(path, index=":"):
        if abs(float(frame.info["q_target"]) - q) < 1e-6:
            return frame
    raise LookupError(f"no frame with q_target == {q} in {path}")


def kabsch_rmsd(a, b):
    a = a - a.mean(axis=0)
    b = b - b.mean(axis=0)
    u, _, vt = np.linalg.svd(a.T @ b)
    d = np.sign(np.linalg.det(u @ vt))
    rot = u @ np.diag([1.0, 1.0, d]) @ vt
    return float(np.sqrt(np.mean(np.sum((a @ rot - b) ** 2, axis=1))))


def dist(p, i, j):
    return float(np.linalg.norm(p[i] - p[j]))


def _plain(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, np.generic):
        return obj.item()
    if isinstance(obj, dict):
        return {str(k): _plain(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_plain(v) for v in obj]
    return obj


def _vs_reference(saddle_atoms, saddle_energy):
    ts = read(TS)
    p_ts = ts.positions
    swap = list(range(len(ts)))
    swap[CL_LEAVE], swap[CL_NUC] = CL_NUC, CL_LEAVE
    p = saddle_atoms.positions
    e_ts = float(ts.get_potential_energy())  # extxyz energy= goes to the attached calc
    return {
        "reference": str(TS.relative_to(ROOT)),
        "rmsd_A": min(kabsch_rmsd(p, p_ts), kabsch_rmsd(p, p_ts[swap])),
        "rmsd_A_identity_labels": kabsch_rmsd(p, p_ts),
        "rmsd_A_swapped_cl_labels": kabsch_rmsd(p, p_ts[swap]),
        "saddle_r_C_Cl": [dist(p, C, CL_LEAVE), dist(p, C, CL_NUC)],
        "ts0000_r_C_Cl": [dist(p_ts, C, CL_LEAVE), dist(p_ts, C, CL_NUC)],
        "saddle_energy_eV": saddle_energy,
        "ts0000_energy_eV": e_ts,
        "energy_diff_saddle_minus_ts0000_eV": saddle_energy - e_ts,
    }


def _ledger_stages(ledger_path):
    stages = []
    for line in Path(ledger_path).read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        ev = json.loads(line)
        if ev["event"] == "close":
            stages.append(
                {
                    k: ev.get(k)
                    for k in (
                        "task",
                        "stage",
                        "attempt",
                        "forwards",
                        "requests",
                        "cache_hits",
                        "wall_s",
                        "outcome",
                    )
                }
            )
    return stages


def _write_json(out, record):
    (out / "certify.json").write_text(json.dumps(_plain(record), indent=2, default=str))


def main(out_dir):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=False)
    record = {"status": None, "error": None, "spent": None}
    ledger = None
    t0 = time.monotonic()
    try:
        config = SearchConfig(charge_sensitive=True, total_charge=-1, multiplicity=1)
        checkpoint, model_type = MODELS["polar1"]
        record["provenance"] = {
            "search_config": config.to_dict(),
            "seed": SEED,
            "ledger_budget": BUDGET,
            "ceiling_rise_eV": CEILING_RISE_EV,
            "model": "polar1",
            "checkpoint": checkpoint,
            "model_type": model_type,
            "device": DEVICE,
            "default_dtype": DTYPE,
            "hostname": socket.gethostname(),
            "python": sys.executable,
            "python_version": platform.python_version(),
        }
        frame = find_guess_frame()
        guess = _stamped_with_charge(frame, config)
        guess.calc = None
        ledger = Ledger(out / "ledger.jsonl", arm="smoke", budget=BUDGET)
        factory = counted_factory(lambda: polar1(DEVICE, DTYPE), ledger)
        with ledger.segment(task=TASK, stage="certify_prepare", attempt=0) as seg:
            probe = guess.copy()
            probe.calc = factory()
            guess_energy = float(probe.get_potential_energy())
            seg.outcome = "guess_energy"
        record.update(
            guess={
                "file": str(GUESS_FILE.relative_to(ROOT)),
                "q_target": GUESS_Q,
                "energy_eV": guess_energy,
                "stored_energy_eV": float(frame.info["energy_eV"])
                if "energy_eV" in frame.info
                else None,
            },
            ceiling_eV=guess_energy + CEILING_RISE_EV,
            reference_eV=guess_energy,
        )
        _write_json(out, record)

        cert = certify(
            guess,
            factory,
            config,
            ledger,
            task=TASK,
            seed=SEED,
            ceiling_eV=guess_energy + CEILING_RISE_EV,
            reference_eV=guess_energy,
        )
        record.update(
            status=cert.status,
            saddle_energy_eV=cert.saddle_energy_eV,
            unstable_mode=cert.unstable_mode,
            hessian=cert.hessian,
            sides=[{k: v for k, v in s.items() if k != "positions"} for s in cert.sides],
            event=cert.event,
            duplicate_of=cert.duplicate_of,
            chemical_pair=cert.chemical_pair,
            diagnostics=cert.diagnostics,
            certify_spent=cert.spent,
        )
        hess = cert.hessian or {}
        record["imaginary_wavenumbers"] = {
            k: v for k, v in hess.items() if "wavenumber" in k or "imag" in k or "freq" in k
        }

        template = guess.copy()
        if cert.saddle_positions is not None:
            saddle = template.copy()
            saddle.positions = cert.saddle_positions
            saddle.info["energy_eV"] = cert.saddle_energy_eV
            write(out / "saddle.extxyz", saddle)
            record["vs_reference_ts"] = _vs_reference(saddle, cert.saddle_energy_eV)
        for side in cert.sides:
            if side.get("positions") is not None:
                end = template.copy()
                end.positions = side["positions"]
                end.info["energy_eV"] = side["energy_eV"]
                name = "plus" if side["sign"] == 1 else "minus"
                write(out / f"{name}.extxyz", end)
    except BaseException as exc:
        record["error"] = repr(exc)
        raise
    finally:
        record["wall_total_s"] = time.monotonic() - t0
        if ledger is not None:
            ledger.close()
            record["spent"] = ledger.spent
            try:
                record["stages"] = _ledger_stages(out / "ledger.jsonl")
            except Exception as exc:  # keep the record even if the readback fails
                record["stages_error"] = repr(exc)
        _write_json(out, record)


if __name__ == "__main__":
    main(sys.argv[1])
