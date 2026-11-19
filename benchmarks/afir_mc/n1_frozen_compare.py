"""Frozen-weight SC coordinate vs rho_AB on the three known N1 failure paths (2026-11-15,
user: "先做冻结权重的推力坐标，保留现有 γ 和 SC 分组。用三条已知失败路径对照，直接看 H 接触是否
停止垄断、碳骨架是否真正向反应区推进").

Same items as N1 SC sA (prepare.json copied in, SC nodes re-minimized): same perturbed start
x0, same A / B groups, same gamma. Each path is run four ways in one batch:
  rho    / frozen : bias coordinate rho_AB (11-15 runs) or q = sum c_ij r_ij, c frozen at x0
  stop   / nostop : the runner's stops, or the same with the "two maxima" stop disabled
                    (noise_eV set huge so no maximum is ever counted; diagnostic only, the
                    runner is unchanged), so a path ends at the biased minimum / step cap.
Per frame: E, coordinate, rho_AB, the two forming C-C pairs (and the other regiochemistry),
the shortest intermolecular H...H, the closest pair by r/(R_i+R_j) and its force share
(share of |dcoord/dr_ij| over all A x B pairs; for frozen it is c_ij). Not resumable (short).

    python n1_frozen_compare.py --run DIR --seed S --shell R --labels L1,L2 --out OUT
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
from explore.afir import FrozenWeightedDistance  # noqa: E402
from explore.afir_batch import run_paths  # noqa: E402
from explore.afir_path import PathSettings  # noqa: E402
from explore.ledger import Ledger  # noqa: E402
from explore.parallel import MACEBatchEvaluator  # noqa: E402
from explore.poses import sample_poses  # noqa: E402

FORMING = {"C1-C6": (0, 11), "C4-C7": (3, 12), "C1-C7": (0, 12), "C4-C6": (3, 11)}
DIENE_ATOMS = set(range(11))  # cyclopentadiene 0..10, acrolein 11..18


def frame_metrics(numbers, x, coord, rho):
    x = np.asarray(x, float)
    sym_h = numbers == 1
    a, b = list(coord.a), list(coord.b)
    d = x[b][None, :, :] - x[a][:, None, :]
    r = np.linalg.norm(d, axis=-1)
    rs = covalent_radii[numbers[a]][:, None] + covalent_radii[numbers[b]][None, :]
    if isinstance(coord, FrozenWeightedDistance):
        share = coord.c / coord.c.sum()
    else:
        w = (rs / r) ** coord.p
        g = w * np.abs(6.0 * rho.value(x) / r - 5.0)
        share = g / g.sum()
    ia, jb = np.unravel_index(np.argmin(r / rs), r.shape)
    hh = [
        np.linalg.norm(x[i] - x[j])
        for i in range(len(x))
        for j in range(i + 1, len(x))
        if sym_h[i] and sym_h[j] and ((i in DIENE_ATOMS) != (j in DIENE_ATOMS))
    ]
    return {
        "rho_AB": round(rho.value(x), 3),
        "forming": {
            k: round(float(np.linalg.norm(x[i] - x[j])), 3) for k, (i, j) in FORMING.items()
        },
        "hh_min": round(float(min(hh)), 3),
        "closest_pair": f"{a[ia]}-{b[jb]}",
        "closest_ratio": round(float((r / rs)[ia, jb]), 3),
        "closest_share": round(float(share[ia, jb]), 3),
        "h_pair_share": round(
            float(
                sum(
                    share[p, q]
                    for p in range(len(a))
                    for q in range(len(b))
                    if sym_h[a[p]] or sym_h[b[q]]
                )
            ),
            3,
        ),  # fmt: skip
    }


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
    off = base.max_steps + 1
    stop = PathSettings(
        noise_eV=noise, gamma_kJ=args.gamma, early_discard_steps=off, stop_after_maxima=2
    )
    nostop = PathSettings(noise_eV=1e6, gamma_kJ=args.gamma, early_discard_steps=off)
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
            variants.append((lb, cname, coord, rho, it["start"], it["ref_edges"]))
    report = []
    ledger = Ledger(out / "compare.ledger.jsonl", arm="AFIR", replay=False)
    for sname, settings in (("stop", stop), ("nostop", nostop)):
        paths = run_paths(
            atoms,
            [v[4] for v in variants],
            [v[2] for v in variants],
            batch,
            config,
            settings,
            ledger,
            f"compare|{sname}",
        )
        for (lb, cname, coord, rho, _, ref), p in zip(variants, paths):
            fr = p.frames
            E = np.array([f["energy_eV"] for f in fr])
            low = int(np.argmin(E[: len(fr) // 3 + 1]))
            ms = [frame_metrics(atoms.numbers, f["positions"], coord, rho) for f in fr]
            tag = f"{lb.replace('|', '_')}_{cname}_{sname}"
            imgs = []
            for f in fr:
                a = atoms.copy()
                a.positions = np.asarray(f["positions"], float)
                a.info.update(energy_eV=f["energy_eV"], coord_A=f["rho_A"])
                imgs.append(a)
            write(out / f"{tag}.extxyz", imgs)
            report.append(
                {
                    "label": lb,
                    "coord": cname,
                    "stops": sname,
                    "status": p.status,
                    "steps": len(fr),
                    "bond_change": fr[-1]["edges"] != ref,
                    "low_frame": low,
                    "E_end_minus_low_eV": round(float(E[-1] - E[low]), 3),
                    "E_max_after_low_eV": round(float(E[low:].max() - E[low]), 3),
                    "coord_low_end": [round(fr[low]["rho_A"], 3), round(fr[-1]["rho_A"], 3)],
                    "low": ms[low],
                    "end": ms[-1],
                    "forming_min_after_low": {
                        k: min(m["forming"][k] for m in ms[low:]) for k in FORMING
                    },
                    "hh_min_after_low": min(m["hh_min"] for m in ms[low:]),
                    "trace_every_5": [
                        {"step": s, "E": round(float(E[s] - E[low]), 3), **ms[s]}
                        for s in range(0, len(fr), 5)
                    ],
                }
            )
    ledger.close()
    (out / "compare.json").write_text(json.dumps(report, indent=1) + "\n")
    for r in report:
        print(
            json.dumps(
                {k: v for k, v in r.items() if k != "trace_every_5"}, ensure_ascii=False
            ),
            flush=True,
        )
    print(json.dumps({"forwards": ledger.spent, "wall_s": time.monotonic() - t0}), flush=True)


if __name__ == "__main__":
    main()
