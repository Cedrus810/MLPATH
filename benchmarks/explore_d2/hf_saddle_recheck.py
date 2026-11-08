"""Diagnostic for saddles whose order `confirm_minimum` left undetermined: the full
internal eigenvalue spectrum (`curvature_spectrum`, the rigid-body block projected out)
and the internal-mode recheck (`internal_mode_recheck`, finite differences along the
deciding eigenvectors, as confirm_minimum runs it after its floor gate). The floor gate is
NOT overridden: a result here is evidence next to the gate's verdict, not a certificate.
Read-only on the run; the run's single float32 POLAR-1.

    python hf_saddle_recheck.py OUT_DIR RUN_DIR SUFFIX [SUFFIX ...]
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT / "src"), str(ROOT / "explore" / "src")]

from ase import Atoms  # noqa: E402

from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.ledger import Ledger, view_factory  # noqa: E402
from explore.poses import _jsonable  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402
from prrs.runner import curvature_spectrum, floor_residual, internal_mode_recheck  # noqa: E402

CAP = 4000


def _lines(path):
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.strip():
                yield json.loads(line)


def main():
    out_dir, run_dir, suffixes = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3:]
    out_dir.mkdir(parents=True, exist_ok=True)
    m = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    config = SearchConfig.from_dict(m["config"])
    poses = {}
    for r in _lines(run_dir / "poses.jsonl"):
        a = Atoms(numbers=r["numbers"], positions=r["positions"])
        a.info.update(r["info"])
        poses[r["pose_id"]] = a
    pose_of = {r["intent_id"]: r["pose_id"] for r in _lines(run_dir / "population.jsonl")}
    found = [
        ev
        for ev in _lines(run_dir / "journal.jsonl")
        if ev.get("event") == "settle"
        and ev["task"].startswith("HC|")
        and ev["status"] == "connected"
        and any(ev["task"].endswith(s) for s in suffixes)
    ]
    calc = polar1("cuda", "float32")
    for i, ev in enumerate(found):
        out = out_dir / f"recheck{i}.json"
        if out.exists():
            continue
        cand = poses[pose_of[ev["task"][3:]]].copy()
        cand.set_constraint()
        cand.positions = np.asarray(ev["payload"]["saddle_positions"], float)
        ledger = Ledger(out_dir / f"recheck{i}_ledger.jsonl", arm="HF", budget=CAP)
        factory = view_factory(calc, ledger)
        cand.calc = factory()
        t0 = time.monotonic()
        forces = cand.get_forces()
        eigenvalues, vectors, floor, provenance = curvature_spectrum(
            cand, factory, config, source=config.hessian_source
        )
        recheck = internal_mode_recheck(cand, factory, config, eigenvalues, vectors, provenance)
        ledger.close()
        n_rigid = 6
        near_zero = np.argsort(np.abs(eigenvalues))[:n_rigid]
        internal = np.sort(np.delete(eigenvalues, near_zero))
        record = {
            "task": ev["task"],
            "saddle_energy_eV": ev["payload"]["saddle_energy_eV"],
            "fmax_eV_A": float(np.linalg.norm(forces, axis=1).max()),
            "provenance": _jsonable(provenance),
            "eigenvalue_tol": config.minimum_check_eigenvalue_tol,
            "floor_gate_limit": config.trivial_floor_fraction
            * config.minimum_check_eigenvalue_tol,
            "trivial_mode_floor_worst_residual": floor_residual(floor),
            "internal_eigenvalues_lowest": [float(v) for v in internal[:8]],
            "negative_internal": [
                float(v) for v in internal if v < -config.minimum_check_eigenvalue_tol
            ],
            "internal_mode_recheck": _jsonable(recheck),
            "forwards": ledger.spent,
            "wall_s": time.monotonic() - t0,
        }
        out.write_text(json.dumps(record, indent=1) + "\n", encoding="utf-8")
        print(
            json.dumps(
                {
                    "task": ev["task"][-30:],
                    "fmax": record["fmax_eV_A"],
                    "floor": record["trivial_mode_floor_worst_residual"],
                    "lowest_internal": record["internal_eigenvalues_lowest"][:4],
                    "recheck_agree": (recheck or {}).get("agree"),
                    "forwards": ledger.spent,
                }
            ),
            flush=True,
        )


if __name__ == "__main__":
    main()
