"""CI-NEB on the AFIR path as the TS locator: comparison on the N1 DA paths of PBS 29729
(2026-11-19, user 「可以试试」 after "我们是不是可以从 别的算法绕开得到突破呢").

E/F only (no Hessian, no curvature estimate): the chain is relaxed with ASE NEB
(improved tangent, climbing image, spring k) and FIRE; its endpoints are the unbiased E/F
minima of the path's start and end; the interior images are the AFIR path frames from the
early low-energy frame to the end, resampled evenly in Cartesian arc length. The climbing
image at convergence (max NEB force <= --fmax) is the TS candidate.

The chosen paths are replayed (same items as the run: prepare.json copied, SC nodes
re-minimized, two-pair items regenerated, same settings: pair / two-pair push, no
two-maxima stop) so their frames are available; a replay that does not change the bond
graph is reported and skipped.
Writes OUT/neb.json (per path: replay, NEB steps, forwards, climbing-image energy and
forming C-C distances) and OUT/wave_0000.json in the runner's wave format (every converged
climbing image as a converged_candidate), so benchmarks/afir_mc/afir_validate.py checks it
with exactly the criteria of the AFIR runs (order, reference / mirror match). Not resumable.

    python n1_cineb_check.py --run DIR --seed S --shell R --sc-poses P --labels L1,L2 --out OUT
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

from ase.calculators.calculator import Calculator, all_changes  # noqa: E402
from ase.mep import NEB  # noqa: E402
from ase.optimize import BFGS, FIRE  # noqa: E402

from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.afir_batch import run_paths  # noqa: E402
from explore.afir_path import PathSettings  # noqa: E402
from explore.ledger import Ledger  # noqa: E402
from explore.parallel import MACEBatchEvaluator  # noqa: E402
from explore.poses import _jsonable, sample_poses  # noqa: E402

FORMING = {"C1-C6": (0, 11), "C4-C7": (3, 12), "C1-C7": (0, 12), "C4-C6": (3, 11)}


class Counted(Calculator):
    """Shared unbiased calculator that counts forwards."""

    implemented_properties = ["energy", "forces"]

    def __init__(self, base):
        super().__init__()
        self.base, self.calls = base, 0

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        a = self.atoms.copy()
        a.calc = self.base
        self.calls += 1
        self.results = {"energy": float(a.get_potential_energy()), "forces": a.get_forces()}


def resample(frames, n):
    """n structures evenly spaced in cumulative Cartesian arc length along frames."""
    x = np.stack([np.asarray(f, float) for f in frames])
    s = np.concatenate(
        [[0.0], np.cumsum(np.linalg.norm((x[1:] - x[:-1]).reshape(len(x) - 1, -1), axis=1))]
    )
    targets = np.linspace(0.0, s[-1], n)
    out = []
    for t in targets:
        k = min(int(np.searchsorted(s, t, side="right")) - 1, len(x) - 2)
        w = 0.0 if s[k + 1] == s[k] else (t - s[k]) / (s[k + 1] - s[k])
        out.append((1 - w) * x[k] + w * x[k + 1])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--system", default="n1")
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--shell", type=float, required=True)
    ap.add_argument("--sc-poses", type=int, required=True)
    ap.add_argument("--max-paths", type=int, default=4096)
    ap.add_argument("--gamma", type=float, default=106.9)
    ap.add_argument("--labels", required=True)
    ap.add_argument("--images", type=int, default=12, help="interior images")
    ap.add_argument("--k", type=float, default=0.1, help="spring constant, eV/A^2")
    ap.add_argument("--fmax", type=float, default=0.01, help="NEB force criterion, eV/A")
    ap.add_argument("--steps", type=int, default=2000)
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
    poses = [
        p
        for p in sample_poses(atoms, frags, confs[0], confs[1], [args.shell], args.max_paths,
                              args.seed)
        if p.status == "placed"
    ]  # fmt: skip
    batch = MACEBatchEvaluator(calc)
    doubles = []
    singles, _ = sc_items(
        atoms, poses, config, batch, out,
        SimpleNamespace(seed=args.seed, sc_poses=args.sc_poses, sc_coord="pair", sc_double=1),
        doubles,
        double_filter=lambda lb: lb in labels,
    )  # fmt: skip
    by_label = {it["label"]: it for it in singles + doubles}
    chosen = [by_label[lb] for lb in labels]
    ledger = Ledger(out / "replay.ledger.jsonl", arm="NEB", replay=False)
    paths = run_paths(atoms, [it["start"] for it in chosen], [it["coord"] for it in chosen],
                      batch, config, settings, ledger, "neb|replay")  # fmt: skip
    ledger.close()
    t_replay = time.monotonic() - t0

    counted = Counted(calc)
    report, wave_paths = [], []
    for it, p in zip(chosen, paths):
        fr = p.frames
        rec = {"label": it["label"], "replay_status": p.status, "replay_steps": len(fr),
               "replay_forwards": p.spent,
               "bond_change": bool(fr) and fr[-1]["edges"] != it["ref_edges"]}  # fmt: skip
        if not rec["bond_change"]:
            rec["neb"] = "skipped: replay did not change the bond graph"
            report.append(rec)
            continue
        c0 = counted.calls
        ends = []
        for x in (fr[0]["positions"], fr[-1]["positions"]):
            a = atoms.copy()
            a.set_constraint()
            a.positions = np.asarray(x, float)
            a.calc = counted
            BFGS(a, logfile=None, maxstep=0.1).run(fmax=args.fmax, steps=1000)
            ends.append(a)
        e = np.array([f["energy_eV"] for f in fr])
        low = int(np.argmin(e[: len(fr) // 3 + 1]))
        mids = resample([f["positions"] for f in fr[low:]], args.images + 2)[1:-1]
        images = [ends[0]]
        for x in mids:
            a = atoms.copy()
            a.set_constraint()
            a.positions = x
            a.calc = counted
            images.append(a)
        images.append(ends[1])
        for a in images:
            a.calc = counted
        neb = NEB(images, k=args.k, climb=True, method="improvedtangent",
                  allow_shared_calculator=True)  # fmt: skip
        opt = FIRE(neb, logfile=None, maxstep=0.1)
        converged = bool(opt.run(fmax=args.fmax, steps=args.steps))
        energies = [float(a.get_potential_energy()) for a in images]
        top = int(np.argmax(energies[1:-1])) + 1
        ci = images[top]
        rec.update({
            "neb": "converged" if converged else "not_converged",
            "neb_steps": opt.nsteps, "neb_forwards": counted.calls - c0,
            "ci_index": top, "ci_energy_eV": energies[top],
            "ci_fmax_eV_A": float(np.linalg.norm(ci.get_forces(), axis=1).max()),
            "barrier_from_start_eV": energies[top] - energies[0],
            "reaction_energy_eV": energies[-1] - energies[0],
            "ci_forming_A": {k: round(float(ci.get_distance(i, j)), 3)
                             for k, (i, j) in FORMING.items()},
            "chain_energies_eV": energies,
        })  # fmt: skip
        report.append(rec)
        if converged:
            wave_paths.append({
                "k": len(wave_paths), "label": it["label"], "new_events": [],
                "ts": [{"status": "converged_candidate", "positions": ci.positions.copy(),
                        "energy_eV": energies[top], "link": {}}],
            })  # fmt: skip
    spent = sum(r.get("neb_forwards", 0) + r["replay_forwards"] for r in report)
    (out / "wave_0000.json").write_text(json.dumps(_jsonable(
        {"wave": 0, "paths": wave_paths, "spent_total": spent}), default=str) + "\n")  # fmt: skip
    (out / "neb.json").write_text(json.dumps(_jsonable(
        {"args": vars(args), "replay_s": t_replay, "wall_s": time.monotonic() - t0,
         "paths": report}), indent=1, default=str) + "\n")  # fmt: skip
    for r in report:
        print(json.dumps({k: v for k, v in r.items() if k != "chain_energies_eV"}), flush=True)
    print(json.dumps({"wall_s": time.monotonic() - t0, "replay_s": t_replay}), flush=True)


if __name__ == "__main__":
    main()
