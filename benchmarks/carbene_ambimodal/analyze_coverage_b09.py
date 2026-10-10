"""What did the compress-only coverage run on b09 reach, and which DFT saddle is each TS candidate?

runs/b09_compress_coverage pushed all 84 compress directions on the source microstate at four
amplitudes (336 trials). Three questions:

1. Per trial: outcome, bond changes, and how close the carbene carbon came to Ca / Cb / the
   aldehyde O anywhere along the response trajectory (not only at the quenched endpoint).
2. Per TS candidate: imaginary wavenumber, POLAR-1 energy above the source microstate, the
   classify_b09 markers, heavy-atom RMSD to the four DFT saddles.
3. The new microstates the run admitted.

Markers and RMSD are classify_b09's, so the numbers compare with the probe1 report.

Usage: python benchmarks/carbene_ambimodal/analyze_coverage_b09.py <run dir> <out json>
"""

import json
import sys
from pathlib import Path

import numpy as np
from ase.io import read

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path[:0] = [str(HERE), str(ROOT / "src"), str(ROOT)]
from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from classify_b09 import DFT, atoms_of_interest, heavy_rmsd, markers  # noqa: E402

KEYS = ("C:-Ca", "C:-Cb", "C:-O_ald")


def energy(atoms, calc):
    atoms = atoms.copy()
    atoms.info.update(charge=0, spin=1)
    atoms.calc = calc
    return float(atoms.get_potential_energy())


def closest_approach(path, idx):
    frames = read(path, index=":")
    pairs = {k: (idx["carbene"], idx[k.split(":-")[1]]) for k in KEYS}
    return {
        k: float(min(np.linalg.norm(f.positions[i] - f.positions[j]) for f in frames))
        for k, (i, j) in pairs.items()
    }


def main():
    run, out = Path(sys.argv[1]), Path(sys.argv[2])
    ref = {name: read(HERE / "structures" / "ts" / f"{name}.xyz") for name in DFT}
    idx = atoms_of_interest(ref["TS-2OMe"])
    network = json.loads((run / "network.json").read_text())
    calc = polar1()
    source = read(run / "chemical" / "c0000_m0000.extxyz")
    e0 = energy(source, calc)
    print("atoms:", idx, " source markers", {k: round(v, 3) for k, v in markers(source, idx).items()})
    for name, atoms in ref.items():
        print(f"DFT {name:14s} " + "  ".join(f"{k} {v:.3f}" for k, v in markers(atoms, idx).items()))

    trials = []
    for line in (run / "attempts.jsonl").read_text().splitlines():
        a = json.loads(line)
        p, f = a["probe"], a.get("fingerprint") or {}
        row = dict(
            attempt=a["attempt_id"],
            direction=f"{p['family']}:{'-'.join(map(str, p['indices']))}",
            amplitude=p["amplitude"],
            status=a["status"],
            classification=a["classification"],
            formed=f.get("formed", []),
            broken=f.get("broken", []),
            closest=closest_approach(run / a["trajectory"], idx),
        )
        trials.append(row)
    print(f"\n{len(trials)} trials; closest carbene approach over all trajectories:")
    for k in KEYS:
        best = min(trials, key=lambda r: r["closest"][k])
        print(f"  {k}: {best['closest'][k]:.3f} A  ({best['attempt']} {best['direction']} {best['amplitude']})")
    print("\nbond changes:")
    for r in trials:
        if r["formed"] or r["broken"] or r["classification"] == "rearrangement":
            print(f"  {r['attempt']} {r['direction']:16s} {r['amplitude']:.1f} {r['status']:13s} "
                  f"{r['classification']}  formed {r['formed']} broken {r['broken']}")
    watch = ("compress:2-13", "compress:2-15", "compress:2-1")
    print("\nwatched directions (carbene C -> Ca / Cb / aldehyde O):")
    for r in trials:
        if r["direction"] in watch:
            print(f"  {r['attempt']} {r['direction']:16s} {r['amplitude']:.1f} {r['status']:13s} "
                  + "  ".join(f"{k} {v:.2f}" for k, v in r["closest"].items()))

    by_attempt = {r["attempt"]: r for r in trials}
    candidates = []
    print("\nTS candidates:")
    for ts in network["ts_candidates"]:
        atoms = read(run / ts["structure"])
        rmsd = {k: heavy_rmsd(v, atoms) for k, v in ref.items()}
        t = by_attempt[ts["attempt"]]
        row = dict(
            id=ts["id"],
            attempt=ts["attempt"],
            direction=t["direction"],
            amplitude=t["amplitude"],
            imaginary_cm=ts["imaginary_wavenumbers_icm"],
            dE_meV=1000 * (energy(atoms, calc) - e0),
            markers=markers(atoms, idx),
            rmsd_to_source=heavy_rmsd(source, atoms),
            heavy_rmsd_to_dft=rmsd,
            closest_dft=min(rmsd, key=rmsd.get),
        )
        candidates.append(row)
        print(f"  {row['id']} {row['direction']:16s} {row['amplitude']:.1f} "
              f"imag {[round(v) for v in row['imaginary_cm']]}  dE {row['dE_meV']:+7.1f} meV  "
              + "  ".join(f"{k} {v:.2f}" for k, v in row["markers"].items())
              + f"  RMSD src {row['rmsd_to_source']:.2f}  "
              + " / ".join(f"{v:.2f}" for v in rmsd.values()) + f" -> {row['closest_dft']}")

    micro = []
    print("\nmicrostates:")
    for path in sorted((run / "chemical").glob("*.extxyz")):
        atoms = read(path)
        row = dict(name=path.stem, dE_meV=1000 * (energy(atoms, calc) - e0),
                   markers=markers(atoms, idx), rmsd_to_source=heavy_rmsd(source, atoms))
        micro.append(row)
        print(f"  {row['name']} dE {row['dE_meV']:+7.1f} meV  RMSD src {row['rmsd_to_source']:.2f}  "
              + "  ".join(f"{k} {v:.2f}" for k, v in row["markers"].items()))

    out.write_text(json.dumps(dict(atoms=idx, trials=trials, ts_candidates=candidates,
                                   microstates=micro, ledger=network["direction_ledger"]), indent=1))


if __name__ == "__main__":
    main()
