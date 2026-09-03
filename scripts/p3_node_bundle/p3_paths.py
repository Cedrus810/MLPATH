"""HEAVY 5. Record the delivered path, so there is something to project onto the surface.

The 2D surface on its own is a picture of the potential. What makes it a picture of the
METHOD is our own points and our own path drawn on it. Until now only single frames existed
(A, B, TS), so there was nothing to draw.

PathRecorder already does this and the search already uses it; the probe scripts just never
attached one. Three stages are recorded here, each to its own file:

    climb      follow_min_mode from a seed near the dividing surface up to the saddle
    descent+1  descend_saddle off the saddle, one sign
    descent-1  the other sign

`tracked_coordinates` is set to the two projection coordinates, so every frame carries
q_PT and r_OO by name in the jsonl -- no re-derivation, and no chance of the figure using a
different definition than the record.

A FAILED descent still produces frames, and they are the point rather than a consolation.
PRRS_STATUS.md 8.19 ("the delivered path itself has to be recorded") exists because a walk
that ended in the wrong basin used to leave a verdict and no way to see why. P3's open
question is exactly that: `side +1` has never reached A, with reason
soft_mode_moved_on_the_last_round at quench 2e-3 and fmax_not_reached at 2e-4. Where that
motion goes in (q_PT, r_OO), and where it stalls, is visible in these frames and nowhere
else. So nothing here treats a failed descent as an absent one.

PathRecorder triggers no calculation of its own -- it records energy and force only when
the caller already has them or the calculator has them cached -- so attaching it does not
change the cost of the walk it observes.

Output per model under runs/p3_paths/<model>/: climb.extxyz, descent.extxyz and the
matching .jsonl. Feed the descent file to p3_pes2d.py as PATHFILE for the second panel.
"""
import sys, os, json, time
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "src"))
from _env_guard import guard, shard_arg, model_info

SHARD, NSHARD, ARGV = shard_arg(sys.argv)
MACEDIR = Path(os.environ.get("P3_MACE_DIR", "/home/ruigengji/MLP/mace"))
CANDIDATES = [MACEDIR / "MACE-OFF24_medium.model",
              MACEDIR / "mace-omol-0-extra-large-4M.model",
              MACEDIR / "MACE-POLAR-1-M.model"]
ONLY = os.environ.get("P3_ONLY_MODEL", "")
if ONLY:
    CANDIDATES = [p for p in CANDIDATES if ONLY in p.name]
OUT = Path(ARGV[1]) if len(ARGV) > 1 else HERE.parents[1] / "runs" / "p3_paths"
QUENCH = float(os.environ.get("P3_QUENCH", "0.002"))
QUENCH_STEPS = int(os.environ.get("P3_QUENCH_STEPS", "5000"))
AMP = float(os.environ.get("P3_AMP", "0.15"))

MINE = [m for k, m in enumerate(CANDIDATES) if k % NSHARD == SHARD]
env = guard(CANDIDATES[0], label=f"paths shard {SHARD}/{NSHARD}", models=CANDIDATES)

import numpy as np
from ase.io import read, write
from ase.optimize import LBFGS
from mace.calculators import MACECalculator
from prrs.chemistry import chemical_key
from prrs.config import SearchConfig
from prrs.runner import (PathRecorder, curvature_spectrum, floor_residual,
                         follow_min_mode, descend_saddle, _wavenumbers)
from prrs.state import encode

OUT.mkdir(parents=True, exist_ok=True)
runs = HERE.parents[1] / "runs"
seed_geom = read(str(runs / "p3_A_source_tight.extxyz"))
SYM = list(seed_geom.symbols)
OXY = [i for i, s in enumerate(SYM) if s == "O"]
HYD = [i for i, s in enumerate(SYM) if s == "H"]

# Located from A's own graph, then frozen into tracked_coordinates so the record and any
# figure made from it cannot disagree about what q_PT means.
edges = sorted(tuple(sorted(e)) for e in encode(seed_geom).edges)
oh = [(o, h) for o in OXY for h in HYD if tuple(sorted((o, h))) in edges]
assert len(oh) == 1, oh
DONOR, H = oh[0]
ACC = [o for o in OXY if o != DONOR][0]
TRACKED = (
    {"name": "q_PT", "kind": "bond_difference", "indices": (DONOR, H, ACC, H)},
    {"name": "r_OO", "kind": "bond", "indices": (DONOR, ACC)},
    {"name": "r_donor_H", "kind": "bond", "indices": (DONOR, H)},
    {"name": "r_acceptor_H", "kind": "bond", "indices": (ACC, H)},
)
cfg = SearchConfig(quench_fmax_eV_A=QUENCH, quench_steps=QUENCH_STEPS,
                   tracked_coordinates=TRACKED)
print(f"donor O{DONOR} acceptor O{ACC} shared H{H}", flush=True)
print(f"quench_fmax={QUENCH:.0e} quench_steps={QUENCH_STEPS} seed amplitude={AMP}", flush=True)
print(f"tracked: {[t['name'] for t in TRACKED]}", flush=True)

report = {}
for path in MINE:
    info = model_info(path)
    name = info["name"]
    mdir = OUT / name
    mdir.mkdir(parents=True, exist_ok=True)
    print(f"\n{'='*70}\n{name}\n{'='*70}", flush=True)
    calc = MACECalculator(model_paths=str(path), model_type=info["model_type"],
                          device="cuda", default_dtype="float64")
    factory = lambda c=calc: c

    def prep(atoms):
        x = atoms.copy()
        if info["charge_aware"]:
            x.info["charge"], x.info["spin"] = 0, 1
        x.calc = calc
        return x

    A = prep(seed_geom)
    LBFGS(A, logfile=None).run(fmax=QUENCH, steps=20000)
    E_A = float(A.get_potential_energy())
    kA = chemical_key(A)["key"][:16]
    print(f"  A E={E_A:.6f} key={kA}", flush=True)

    # --- climb, recorded --------------------------------------------------------------
    target = np.zeros((len(A), 3))
    u = A.positions[ACC] - A.positions[DONOR]
    target[H] = u / np.linalg.norm(u)
    seed = A.copy(); seed.positions[H] += AMP * target[H]
    seed = prep(seed)
    climb = PathRecorder(mdir / "climb", cfg, "climb")
    climb.context = {"model": name, "stage_detail": "follow_min_mode", "amplitude": AMP}
    t0 = time.time()
    ts, diag = follow_min_mode(seed, factory, cfg, direction=target,
                               ceiling=E_A + cfg.saddle_search_max_rise_eV,
                               reference=E_A, recorder=climb)
    print(f"  climb: {climb.frames} frames, reason={diag.get('reason')}, "
          f"steps={diag.get('steps')}, order={diag.get('saddle_order')} "
          f"({time.time()-t0:.0f}s)", flush=True)
    rec = {"model": name, "sha256_16": info["sha256_16"], "E_A": E_A, "key_A": kA,
           "quench_fmax_eV_A": QUENCH, "quench_steps": QUENCH_STEPS, "amplitude": AMP,
           "donor": DONOR, "acceptor": ACC, "shared_H": H,
           "tracked": [dict(t, indices=list(t["indices"])) for t in TRACKED],
           "climb_frames": climb.frames, "climb_reason": diag.get("reason"),
           "climb_steps": diag.get("steps"),
           "climb_saddle_order": diag.get("saddle_order")}
    if ts is None:
        # No saddle, but the climb frames exist and show how far the coordinate got.
        rec["status"] = "no_saddle"
        report[name] = rec
        print(f"  no saddle from amplitude {AMP}; climb frames kept anyway", flush=True)
        json.dump(dict(env=env, models=report),
                  open(OUT / f"paths_shard{SHARD:02d}.json", "w"), indent=1, default=str)
        continue

    ts.calc = calc
    lam, vec, fl, pr = curvature_spectrum(ts, factory, cfg, source="analytic")
    neg = [i for i, v in enumerate(lam) if v < -cfg.minimum_check_eigenvalue_tol]
    E_ts = float(ts.get_potential_energy())
    rec.update(status="saddle", E_TS=E_ts, saddle_order=len(neg),
               barrier_vs_A_meV=(E_ts - E_A) * 1e3,
               floor_residual=float(floor_residual(fl)))
    if len(neg) == 1:
        rec["imaginary_cm1"] = _wavenumbers([lam[0]])[0]
    write(str(mdir / "ts.extxyz"), ts)
    write(str(mdir / "A.extxyz"), A)
    print(f"  TS order={len(neg)} barrier={rec['barrier_vs_A_meV']:+.1f} meV "
          f"{rec.get('imaginary_cm1', float('nan')):.1f} cm^-1", flush=True)
    if len(neg) != 1:
        report[name] = rec
        json.dump(dict(env=env, models=report),
                  open(OUT / f"paths_shard{SHARD:02d}.json", "w"), indent=1, default=str)
        continue

    # --- both descents, recorded; a failure keeps its frames --------------------------
    mode = np.asarray(vec[:, neg[0]])
    sides = []
    for sign in (1, -1):
        descent = PathRecorder(mdir / f"descent{sign:+d}", cfg, "descent")
        descent.context = {"model": name, "sign": sign, "stage_detail": "descend_saddle"}
        t0 = time.time()
        ep, checked, relays = descend_saddle(ts, mode, factory, cfg, sign,
                                             seed=7919 * sign, recorder=descent)
        entry = {"sign": sign, "frames": descent.frames, "relays": len(relays),
                 "reached": ep is not None, "seconds": round(time.time() - t0, 1),
                 "reason": checked.get("reason"),
                 "saddle_order": checked.get("saddle_order")}
        if ep is not None:
            ep.calc = calc
            entry.update(E_eV=float(ep.get_potential_energy()),
                         key=chemical_key(ep)["key"][:16])
            write(str(mdir / f"endpoint{sign:+d}.extxyz"), ep)
        sides.append(entry)
        print(f"  descent {sign:+d}: {descent.frames} frames, reached={entry['reached']}, "
              f"reason={entry['reason']}, relays={entry['relays']} "
              f"({entry['seconds']}s)", flush=True)
    rec["sides"] = sides
    got = [x["key"] for x in sides if x.get("reached")]
    rec["both_sides_reached"] = len(got) == 2
    rec["connects_two_nodes"] = len(set(got)) == 2
    report[name] = rec
    json.dump(dict(env=env, models=report),
              open(OUT / f"paths_shard{SHARD:02d}.json", "w"), indent=1, default=str)

json.dump(dict(env=env, models=report),
          open(OUT / f"paths_shard{SHARD:02d}.json", "w"), indent=1, default=str)
print(f"\nshard {SHARD} wrote paths_shard{SHARD:02d}.json  models={list(report)}", flush=True)
print("Frames from a FAILED descent are kept on purpose: where the delivered motion goes "
      "and where it stalls is the open P3 question and is visible nowhere else.", flush=True)
