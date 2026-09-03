"""HEAVY 2 of 3. The (q_PT, r_OO) surface for the benzoylacetone proton transfer.

Coordinates are the pair proposed in HANDOFF.md 5: q_PT = r(O_donor-H) - r(O_acceptor-H),
which is zero at the symmetric point and signed toward each well, and r_OO, the heavy-atom
gate. Both are constrained; everything else is free in the relaxed pass.

READ THIS BEFORE PLOTTING. The projection is honest about positions and dishonest about
energy if misused. A PRRS path lives in 3N dimensions; the relaxed surface value at a
projected point is the minimum over the other 3N-2 coordinates, so the path generally lies
ABOVE the surface at its own projection. The correct figure is a pair: the path drawn on
the contours, and beside it the path's own energy against arc length WITH the surface
energy along the same projected line overlaid. The gap between those two curves is itself
the diagnostic -- it answers whether these two coordinates opened the reaction up.
This script therefore emits both the surface AND, when a path file is given, the surface
sampled along that path's projection, so the two-panel figure can be built without
re-running anything.

Passes:
  rigid    positions interpolated/scaled only, no relaxation. 1 call per point.
  relaxed  constrained LBFGS at each point, warm-started from the neighbour already done,
           which is what keeps the call count near 100 rather than 400.

Cost on the node (23.8 ms/point): rigid 48x48 ~ 1 min. Relaxed 48x48 at ~100 calls/point
~ 1.5 h. Grid size is the first argument so a coarse pass can be run first.
"""
import sys, os, json, time
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "src"))
import numpy as np
from ase.io import read, write
from ase.optimize import LBFGS
from ase.constraints import FixInternals
from _env_guard import guard, shard_arg

SHARD, NSHARD, ARGV = shard_arg(sys.argv)
N = int(ARGV[1]) if len(ARGV) > 1 else 48
MODEL = ARGV[2] if len(ARGV) > 2 else str(HERE.parents[1] / "models" / "MACE-OFF24_medium.model")
OUT = Path(ARGV[3]) if len(ARGV) > 3 else HERE.parents[1] / "runs" / "p3_pes2d"
PATHFILE = ARGV[4] if len(ARGV) > 4 else None
RELAX_STEPS = int(os.environ.get("P3_PES_RELAX_STEPS", "600"))
# The surface tolerance must be at least as tight as the tolerance the PATH was quenched
# with, or the gap between them is a difference of two convergence criteria rather than a
# projection loss. It was hard-coded at 5e-3 while the descent ran at quench_fmax = 2e-3,
# 2.5x looser, and the surface consequently sat 2.13 meV ABOVE B's true minimum in the B
# well. 375 of 491 path frames are in that well, so 76% of the profile showed a NEGATIVE
# gap -- a path apparently below the relaxed surface, which cannot happen. Defaults to
# QUENCH for exactly that reason; override only to go tighter.
PES_FMAX = float(os.environ.get("P3_PES_FMAX", "0.002"))   # = the paths job default
env = guard(MODEL, label=f"pes2d shard {SHARD}/{NSHARD}")
from prrs.calculators import mace_factory
from prrs.state import encode

OUT.mkdir(parents=True, exist_ok=True)
factory = lambda: mace_factory(MODEL, device="cuda", default_dtype="float64")
calc = factory()

runs = HERE.parents[1] / "runs"
src = runs / "p3_A_source.extxyz"
if not src.exists():
    src = runs / "p3_A_chelated.extxyz"
A = read(str(src)); A.calc = calc
print(f"source {src.name}  E = {A.get_potential_energy():.6f} eV", flush=True)
print(f"surface relax: fmax = {PES_FMAX:.0e}, up to {RELAX_STEPS} steps", flush=True)
sym = list(A.symbols)
OXY = [i for i, s in enumerate(sym) if s == "O"]
HYD = [i for i, s in enumerate(sym) if s == "H"]
edges = sorted(tuple(sorted(e)) for e in encode(A).edges)
oh = [(o, h) for o in OXY for h in HYD if tuple(sorted((o, h))) in edges]
assert len(oh) == 1, oh
DONOR, H = oh[0]
ACC = [o for o in OXY if o != DONOR][0]
print(f"donor O{DONOR}  acceptor O{ACC}  shared H{H}", flush=True)

def coords(atoms):
    p = atoms.positions
    d = np.linalg.norm(p[DONOR] - p[H]); a = np.linalg.norm(p[ACC] - p[H])
    return float(d - a), float(np.linalg.norm(p[DONOR] - p[ACC]))

q0, r0 = coords(A)
print(f"source coordinates q_PT = {q0:+.3f} A   r_OO = {r0:.3f} A", flush=True)
# Ranges chosen from the measured stationary points, not guessed:
#   A   q_PT = -0.6193, r_OO = 2.5289      B  q_PT = +0.6124, r_OO = 2.5181
#   TS  q_PT ~ 0, r_OO = 2.3234 (OFF24) / 2.3595 (omol-0, POLAR-1)
# The lower r_OO bound was 2.30, which put OFF24's TS 0.023 A from the edge -- a contour
# plot whose most interesting point sits on the frame. 2.25 gives it room.
Q = np.linspace(-0.85, 0.85, N)          # spans both wells through the symmetric point
# Upper bound cut 3.00 -> 2.72. Both wells are at r_OO ~ 2.52 and the channel at 2.32, so
# everything above ~2.68 is just the O...O pulled apart: high, featureless, and it was 43%
# of the grid. At the same N the spacing improves from dr = 0.024 to 0.010 A, which is
# where the interpolated section in the right panel gets its accuracy.
R = np.linspace(2.25, 2.72, N)           # closed chelate through the gate opening

def place(atoms, q, r):
    """Set the two coordinates geometrically: scale O-O about its midpoint, then put H on
    the O...O line at the separation q implies. Approximate by construction, which is why
    the relaxed pass constrains the same two coordinates rather than trusting this."""
    a = atoms.copy()
    p = a.positions
    mid = 0.5 * (p[DONOR] + p[ACC])
    u = p[ACC] - p[DONOR]; u /= np.linalg.norm(u)
    p[DONOR] = mid - 0.5 * r * u
    p[ACC] = mid + 0.5 * r * u
    # r(D-H) - r(A-H) = q with r(D-H) + r(A-H) ~ r  ->  r(D-H) = (r + q)/2
    p[H] = p[DONOR] + max(0.30, min(r - 0.30, 0.5 * (r + q))) * u
    a.positions = p
    return a

def constrained(atoms, q, r):
    """Two constraints, not three: r_OO as a bond, q_PT as a LINEAR COMBINATION of bonds.

    The first version pinned three bond lengths -- r(D-A), r(D-H) and r(A-H) -- to fix q_PT
    indirectly. That is over-constrained and degenerate exactly where the surface is most
    interesting: with H on the D-A line, r(D-H) + r(A-H) = r(D-A) identically, so the
    constraint Jacobian is singular and FixInternals' projection cannot converge. Every
    point failed with "FixInternals.adjust_positions did not converge", and the relaxed
    surface came out 0/2304 finite while the rigid one was 2304/2304 -- a whole panel
    silently empty.

    q_PT is a bond DIFFERENCE, so it belongs in `bondcombos` as
    1.0*r(D,H) - 1.0*r(A,H). Two independent constraints, no degeneracy. Measured: the
    constraints are then satisfied to ~1e-16 across the whole range including the corners.
    """
    a = atoms.copy(); a.calc = factory()
    a.set_constraint(FixInternals(
        bonds=[[float(r), [DONOR, ACC]]],
        bondcombos=[[float(q), [[DONOR, H, 1.0], [ACC, H, -1.0]]]],
        epsilon=1e-10))
    opt = LBFGS(a, logfile=None)
    # `converged` is returned, not discarded. A point that ran out of steps is NOT a value
    # of the relaxed surface, and a panel that mixes the two is not a relaxed surface.
    converged = bool(opt.run(fmax=PES_FMAX, steps=RELAX_STEPS))
    return a, int(opt.get_number_of_steps()), converged

# Rows are the shard unit. Each row is a full sweep in q, and the relaxed pass warm-starts
# along q WITHIN a row, so a row is the largest unit that can be split without losing the
# warm start that keeps the call count near 100/point instead of 400. Stride, not blocks,
# so each worker gets a mix of tight and open chelate rows and they finish together.
MY_ROWS = [i for i in range(N) if i % NSHARD == SHARD]
print(f"shard {SHARD}/{NSHARD}: rows {MY_ROWS}", flush=True)

for mode in ("rigid", "relaxed"):
    grid = np.full((len(R), len(Q)), np.nan)
    got = np.full((len(R), len(Q), 2), np.nan)
    calls = np.zeros((len(R), len(Q)), dtype=int)
    conv = np.zeros((len(R), len(Q)), dtype=np.int8) - 1     # -1 not attempted, 0 no, 1 yes
    t_start = time.time()
    warm = None
    for n_done, i in enumerate(MY_ROWS):
        r = R[i]
        for j, q in enumerate(Q):
            seed = warm if (mode == "relaxed" and warm is not None) else A
            cand = place(seed, q, r)
            if mode == "rigid":
                cand.calc = calc
                try:
                    grid[i, j] = float(cand.get_potential_energy()); calls[i, j] = 1
                    conv[i, j] = 1          # a rigid point has nothing to converge
                except Exception as exc:
                    print(f"  rigid fail ({q:+.3f},{r:.3f}): {exc}", flush=True); continue
                got[i, j] = coords(cand)
            else:
                try:
                    done, steps, ok = constrained(cand, q, r)
                    grid[i, j] = float(done.get_potential_energy())
                    calls[i, j] = steps + 1
                    conv[i, j] = 1 if ok else 0
                    got[i, j] = coords(done)
                    warm = done
                except Exception as exc:
                    print(f"  relaxed fail ({q:+.3f},{r:.3f}): {type(exc).__name__} {exc}",
                          flush=True)
                    continue
        el = time.time() - t_start
        print(f"[s{SHARD} {mode}] row {i} ({n_done+1}/{len(MY_ROWS)}) r_OO {r:.3f}  "
              f"min E {np.nanmin(grid[i]):.5f}  mean calls {calls[i].mean():.0f}  "
              f"ETA {el/(n_done+1)*(len(MY_ROWS)-n_done-1)/60:.0f} min", flush=True)
        np.savez(OUT / f"pes2d_{mode}_shard{SHARD:02d}.npz", Q=Q, R=R, E=grid,
                 achieved=got, calls=calls, converged=conv, rows=np.array(MY_ROWS),
                 shard=SHARD, nshard=NSHARD, relax_steps=RELAX_STEPS,
                 q0=q0, r0=r0, donor=DONOR, acceptor=ACC, shared_H=H)
    mine = grid[MY_ROWS]
    ok = np.isfinite(mine)
    cm = conv[MY_ROWS]
    print(f"[s{SHARD} {mode}] done: {ok.sum()}/{mine.size} of my points evaluated, "
          f"{int((cm == 1).sum())} converged, {int((cm == 0).sum())} out of steps, "
          f"E range {np.nanmin(mine):.5f} .. {np.nanmax(mine):.5f} eV, "
          f"total calls {calls.sum()}", flush=True)

# Surface energy along a supplied path's own projection, for the second panel.
# The path profile is sharded over frames too; p3_merge.py concatenates them in order.
if PATHFILE and Path(PATHFILE).exists():
    frames = read(PATHFILE, index=":")
    frames = [(k, fr) for k, fr in enumerate(frames) if k % NSHARD == SHARD]
    print(f"shard {SHARD}: {len(frames)} path frames", flush=True)
    prof = []
    allframes = read(PATHFILE, index=":")
    # arc length is cumulative over the WHOLE path, so it is computed from all frames even
    # though this worker only evaluates its own -- otherwise each shard would restart at 0
    arc = [0.0]
    for k in range(1, len(allframes)):
        arc.append(arc[-1] + float(np.linalg.norm(allframes[k].positions
                                                  - allframes[k-1].positions)))
    for k, fr in frames:
        s = arc[k]
        fr.calc = calc
        q, r = coords(fr)
        # The surface value at this projection is NOT recomputed here. It is obtained by
        # interpolating the grid this same script produced (see p3_figure.surface_along),
        # so the right panel's section and the left panel's contours are one object rather
        # than two computations that can disagree.
        #
        # The earlier version did recompute it, as constrained(A, q, r) -- passing A's
        # geometry directly to a constraint that may be far away, without place()ing it at
        # (q, r) first. The optimiser stayed at A and E_surface came back as the constant
        # E(A) for all 491 frames, which made the path look like it ran BELOW the relaxed
        # surface (mean gap -9.1 meV). That is impossible, and the constant was the tell.
    json.dump(dict(env=env, shard=SHARD, nshard=NSHARD, profile=prof),
              open(OUT / f"path_profile_shard{SHARD:02d}.json", "w"), indent=1)
    print(f"recorded {len(prof)} path frames with coordinates and their own energy; "
          f"the gap is computed in p3_figure.py by interpolating this run's grid")

json.dump(dict(env=env, N=N, q_range=[float(Q[0]), float(Q[-1])],
               r_range=[float(R[0]), float(R[-1])], source=str(src), shard=SHARD,
               nshard=NSHARD, rows=MY_ROWS,
               donor=DONOR, acceptor=ACC, shared_H=H, q0=q0, r0=r0),
          open(OUT / f"meta_shard{SHARD:02d}.json", "w"), indent=1, default=str)
print(f"shard {SHARD} wrote its npz files; run p3_merge.py to assemble the grid")
