"""Is a batched force noisier than a single one? (D2-HF-v1.4 measurement follow-up; user
2026-11-08 「可以继续跑 几组计算」). In the 11-08 HS measurement, max_active 16 changed 10
outcomes against serial (8: 5, 4: 4); the fine midpoint threshold 1e-3 sits inside the
force noise (DONE_LOG 11-06), so a larger batch noise would move outcomes.

Frames: b05 HS candidate geometries of the 11-07 v5 / v6 acceptance (near-saddle). For each
of TARGETS frames, REPS evaluations alone through the resident calculator (cache bypassed,
as D2's serial HS), and REPS evaluations as the first frame of a MACEBatchEvaluator batch of
size k (k in SIZES), the rest filled with other candidate frames (as D2-HF-v1.4 HS waves).
Written per target and size: per-component force std over repeats (max, median), the
per-atom max of |mean(batched) - mean(single)|, and the energy offset; no threshold.

    python batch_noise.py OUT_JSON CANDIDATE_DIR RUN_DIR
"""

import glob
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT / "src"), str(ROOT / "explore" / "src")]

from ase import Atoms  # noqa: E402
from ase.calculators.calculator import all_changes  # noqa: E402

from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.parallel import MACEBatchEvaluator  # noqa: E402

TARGETS, REPS, SIZES = 4, 8, (1, 8, 16)


def main():
    out, cand_dir, run_dir = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    pose = json.loads((run_dir / "poses.jsonl").read_text().splitlines()[0])
    frames = []
    for name in sorted({Path(f).name for f in glob.glob(str(cand_dir / "*.json"))}):
        rec = json.loads((cand_dir / name).read_text())
        if isinstance(rec, dict) and rec.get("candidate"):
            a = Atoms(numbers=pose["numbers"], positions=rec["candidate"]["positions"])
            a.info.update(pose["info"])
            frames.append((name, a))
    assert len(frames) >= TARGETS + max(SIZES) - 1, len(frames)
    calc = polar1("cuda", "float32")
    batch = MACEBatchEvaluator(calc)
    targets, fillers = frames[:TARGETS], [a for _, a in frames[TARGETS:]]
    rows = []
    for name, atoms in targets:
        single = []
        for _ in range(REPS):
            calc.calculate(atoms.copy(), ["energy", "forces"], all_changes)
            single.append((float(calc.results["energy"]), np.array(calc.results["forces"])))
        f1 = np.array([f for _, f in single])
        row = {
            "frame": name,
            "single": {
                "std_max": float(f1.std(axis=0).max()),
                "std_median": float(np.median(f1.std(axis=0))),
                "energy_ptp_eV": float(np.ptp([e for e, _ in single])),
            },
        }
        for k in SIZES:
            got = []
            for _ in range(REPS):
                (e, f), *_ = batch([atoms.copy()] + [a.copy() for a in fillers[: k - 1]])
                got.append((float(e), np.asarray(f, float)))
            fk = np.array([f for _, f in got])
            offset = fk.mean(axis=0) - f1.mean(axis=0)
            row[f"batch_{k}"] = {
                "std_max": float(fk.std(axis=0).max()),
                "std_median": float(np.median(fk.std(axis=0))),
                "mean_offset_atom_max": float(np.linalg.norm(offset, axis=1).max()),
                "energy_offset_eV": float(
                    np.mean([e for e, _ in got]) - np.mean([e for e, _ in single])
                ),
                "energy_ptp_eV": float(np.ptp([e for e, _ in got])),
            }
        rows.append(row)
        print(json.dumps(row), flush=True)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rows, indent=1) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
