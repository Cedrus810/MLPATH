"""Where does the push go? (2026-11-15, user 「做」; main line: the frozen push stops outside
the reaction region.) For chosen SC paths of an existing run, replays the path with rho_AB
and with the frozen q (same perturbed start x0, same A / B groups, same gamma, no two-maxima
stop, so each path ends at the biased minimum / step cap) and decomposes the bias force:

    bias force on pair ij = alpha * dcoord/dr_ij     (eV/A, > 0 pushes the pair together)
    frozen q : dq/dr_ij   = c_ij  (fixed at x0)
    rho_AB   : drho/dr_ij = w_ij (6 rho / r_ij - 5) / sum w

At the biased minimum grad E = -alpha grad coord, so this is also what the unbiased forces
balance. Reported at the end frame and at the highest-energy frame after the early low:
the pushed pair (k, l of the SC label) -- its distance, the force it gets and its share of
sum |alpha dcoord/dr| -- and the pairs with the largest force, with r and r/(R_i+R_j).
Model calls: the path replays only. Not resumable (short).

    python push_force_breakdown.py --run DIR --system S --seed N --shell R --sc-poses P \
        --labels L1,L2 --out OUT
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

from ase.data import covalent_radii  # noqa: E402
from ase.io import write  # noqa: E402

from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.afir import FrozenWeightedDistance, alpha_eV_A  # noqa: E402
from explore.afir_batch import run_paths  # noqa: E402
from explore.afir_path import PathSettings  # noqa: E402
from explore.ledger import Ledger  # noqa: E402
from explore.parallel import MACEBatchEvaluator  # noqa: E402
from explore.poses import sample_poses  # noqa: E402


def pair_derivative(coord, x):
    """dcoord/dr_ij over A x B, and the distances r_ij."""
    x = np.asarray(x, float)
    a, b = list(coord.a), list(coord.b)
    r = np.linalg.norm(x[b][None, :, :] - x[a][:, None, :], axis=-1)
    if isinstance(coord, FrozenWeightedDistance):
        return coord.c.copy(), r
    w = (coord.radii_sum / r) ** coord.p
    rho = float(np.sum(w * r) / np.sum(w))
    return w * (6.0 * rho / r - 5.0) / np.sum(w), r


def breakdown(numbers, sym, coord, x, alpha, k, m, top=6):
    der, r = pair_derivative(coord, x)
    a, b = list(coord.a), list(coord.b)
    rs = covalent_radii[numbers[a]][:, None] + covalent_radii[numbers[b]][None, :]
    f = alpha * der
    tot = float(np.abs(f).sum())
    ia, jb = (a.index(k), b.index(m)) if k in a else (a.index(m), b.index(k))
    order = np.dstack(np.unravel_index(np.argsort(-np.abs(f), axis=None), f.shape))[0][:top]
    return {
        "pushed_pair": f"{sym[a[ia]]}{a[ia]}-{sym[b[jb]]}{b[jb]}",
        "pushed_r_A": round(float(r[ia, jb]), 3),
        "pushed_ratio": round(float(r[ia, jb] / rs[ia, jb]), 3),
        "pushed_force_eV_A": round(float(f[ia, jb]), 4),
        "pushed_share": round(float(abs(f[ia, jb]) / tot), 3),
        "pairs_pulled_apart": int((f < 0).sum()),
        "top_pairs": [
            {
                "pair": f"{sym[a[i]]}{a[i]}-{sym[b[j]]}{b[j]}",
                "r_A": round(float(r[i, j]), 3),
                "ratio": round(float(r[i, j] / rs[i, j]), 3),
                "force_eV_A": round(float(f[i, j]), 4),
                "share": round(float(abs(f[i, j]) / tot), 3),
            }
            for i, j in order
        ],
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--system", required=True)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--shell", type=float, required=True)
    ap.add_argument("--sc-poses", type=int, required=True)
    ap.add_argument("--max-paths", type=int, default=4096)
    ap.add_argument("--gamma", type=float, default=106.9)
    ap.add_argument("--labels", required=True)
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
    probe = atoms.copy()
    probe.calc = calc
    noise = 2 * float(np.spacing(np.float32(abs(float(probe.get_potential_energy())))))
    frags, confs, _ = prepare(out, system, config, atoms, calc)
    settings = PathSettings(noise_eV=noise, gamma_kJ=args.gamma, early_discard_steps=1001)
    alpha = alpha_eV_A(args.gamma)
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
    variants = []
    for lb in labels:
        it = by_label[lb]
        rho = it["coord"]
        frozen = FrozenWeightedDistance.build(atoms.numbers, rho.a, rho.b, it["start"])
        for cname, coord in (("rho", rho), ("frozen", frozen)):
            variants.append((lb, cname, coord, it))
    ledger = Ledger(out / "replay.ledger.jsonl", arm="AFIR", replay=False)
    paths = run_paths(
        atoms,
        [v[3]["start"] for v in variants],
        [v[2] for v in variants],
        batch,
        config,
        settings,
        ledger,
        "breakdown|paths",
    )
    ledger.close()
    sym = atoms.get_chemical_symbols()
    report = []
    for (lb, cname, coord, it), p in zip(variants, paths):
        k, m = map(int, lb.split("|")[1].split("-"))
        fr = p.frames
        E = np.array([f["energy_eV"] for f in fr])
        low = int(np.argmin(E[: len(fr) // 3 + 1]))
        hi = low + int(np.argmax(E[low:]))
        imgs = []
        for f in fr:
            a = atoms.copy()
            a.positions = np.asarray(f["positions"], float)
            a.info.update(energy_eV=f["energy_eV"], coord_A=f["rho_A"])
            imgs.append(a)
        write(out / f"{lb.replace('|', '_')}_{cname}.extxyz", imgs)
        report.append(
            {
                "label": lb,
                "coord": cname,
                "status": p.status,
                "steps": len(fr),
                "bond_change": fr[-1]["edges"] != it["ref_edges"],
                "alpha_eV_A": round(alpha, 4),
                "E_end_minus_low_eV": round(float(E[-1] - E[low]), 3),
                "E_max_minus_low_eV": round(float(E[hi] - E[low]), 3),
                "max_frame": hi,
                "at_max": breakdown(
                    atoms.numbers, sym, coord, fr[hi]["positions"], alpha, k, m
                ),
                "at_end": breakdown(
                    atoms.numbers, sym, coord, fr[-1]["positions"], alpha, k, m
                ),
            }
        )
    (out / "breakdown.json").write_text(json.dumps(report, indent=1) + "\n")
    for r in report:
        print(json.dumps(r, ensure_ascii=False), flush=True)
    print(json.dumps({"forwards": ledger.spent, "wall_s": time.monotonic() - t0}), flush=True)


if __name__ == "__main__":
    main()
