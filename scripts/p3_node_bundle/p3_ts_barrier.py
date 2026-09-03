"""HEAVY 3 of 3. TS and barrier for the benzoylacetone chelate proton transfer.

Uses the search's own machinery -- follow_min_mode and descend_saddle out of runner.py --
so the barrier is comparable with P1's 419.14 meV and P2's +361 meV instead of being a
separate ad-hoc construction. descend_saddle is the single implementation of two-sided
descent that the search's connect_saddle and boundary continuation also use; that closure
itself lives inside search.py and is not importable, so both signs are driven here the
same way search.py drives them.

Conventions that are easy to get wrong and are gotten right here:
  curvature_spectrum returns eigenvectors as COLUMNS (`vectors[:, order]`), in the
  MASS-WEIGHTED basis. vec[i] is row i, a slice across all modes, not mode i. Cartesian
  displacement needs dx = v / sqrt(m).
  A target direction is mandatory for this system. Without one, min-mode following climbs
  the globally softest mode; on 3-oxobutanal that returned eight rotor saddles at -65..-80
  cm^-1 while the transfer sat at -3186 (PRRS_STATUS.md 8.16 item 1). Benzoylacetone is
  worse, not better: it has a phenyl twist on top of the same methyl rotor, and the source
  structure's own softest internal mode is +0.0025 eV/A^2.

Three things are reported and kept separate, per the project's evidence ladder:
  saddle order from the analytic Hessian -- a first-order saddle has exactly one negative
  two-sided descent endpoints -- these show two minima are connected
  the barrier -- only meaningful once both of the above hold
cm^-1 is printed only for a confirmed first-order saddle, per PRRS_CURVATURE_ROUTES.md E.

Cost on the node: min-mode following dominates; budget 30-90 min for the seed scan.
"""
import sys, os, json, time
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "src"))
import numpy as np
from ase.io import read, write
from _env_guard import guard, shard_arg

SHARD, NSHARD, ARGV = shard_arg(sys.argv)
MODEL = ARGV[1] if len(ARGV) > 1 else str(HERE.parents[1] / "models" / "MACE-OFF24_medium.model")
OUT = Path(ARGV[2]) if len(ARGV) > 2 else HERE.parents[1] / "runs" / "p3_ts_barrier"
env = guard(MODEL, label=f"ts shard {SHARD}/{NSHARD}")
from prrs.calculators import mace_factory
from prrs.chemistry import chemical_key, classify_transition, reaction_event_key
from prrs.config import SearchConfig
from prrs.runner import (curvature_spectrum, floor_residual, follow_min_mode,
                         descend_saddle, confirm_minimum, _wavenumbers)
from prrs.state import encode
from ase import units

OUT.mkdir(parents=True, exist_ok=True)
factory = lambda: mace_factory(MODEL, device="cuda", default_dtype="float64")
# The ONE variable this comparison changes. Everything else -- source structure, target
# direction, amplitude list, model, implementation -- is identical between the two groups,
# so a difference in outcome is attributable to the tolerance and to nothing else. P2's
# freeze lists "Q1 rising from 12/16 to 16/16 is circumstantial, not isolated attribution"
# as a non-claim precisely because three changes went in together there; this does not
# repeat that.
QUENCH = float(os.environ.get("P3_QUENCH", "0.002"))
QUENCH_STEPS = int(os.environ.get("P3_QUENCH_STEPS", "500"))
cfg = SearchConfig(quench_fmax_eV_A=QUENCH, quench_steps=QUENCH_STEPS)
print(f"quench_fmax_eV_A = {QUENCH:.0e}   quench_steps = {QUENCH_STEPS}", flush=True)
# BOTH are varied, because quench_fmax alone is not one variable. It sets the convergence
# criterion AND, at a fixed step budget, whether that criterion is reachable at all. The
# first attempt at this comparison changed only the tolerance, 2e-3 -> 2e-4, and the failure
# mode of the two-sided descent changed from "soft_mode_moved_on_the_last_round" (the quench
# converged, the soft rotor then moved) to "fmax_not_reached" (FIRE never got there in 500
# steps) -- for all three models. That measured the budget, not the tolerance. A 2x2 over
# (tolerance, steps) separates them.
calc = factory()
runs = HERE.parents[1] / "runs"

def pick(*names):
    for n in names:
        if (runs / n).exists(): return runs / n
    raise SystemExit(f"none of {names} found in {runs}")

A = read(str(pick("p3_A_source.extxyz", "p3_A_chelated.extxyz")))
B = read(str(pick("p3_B_source.extxyz", "p3_B_chelated.extxyz")))
A.calc = calc; B.calc = factory()
E_A, E_B = float(A.get_potential_energy()), float(B.get_potential_energy())
kA, kB = chemical_key(A)["key"][:16], chemical_key(B)["key"][:16]
sym = list(A.symbols)
OXY = [i for i, s in enumerate(sym) if s == "O"]
HYD = [i for i, s in enumerate(sym) if s == "H"]
edges = sorted(tuple(sorted(e)) for e in encode(A).edges)
oh = [(o, h) for o in OXY for h in HYD if tuple(sorted((o, h))) in edges]
assert len(oh) == 1, f"expected one O-H in A, got {oh}"
DONOR, H = oh[0]
ACC = [o for o in OXY if o != DONOR][0]
print(f"A {kA} E={E_A:.6f}   B {kB} E={E_B:.6f}   dE(B-A)={(E_B-E_A)*1e3:+.3f} meV")
print(f"donor O{DONOR}  acceptor O{ACC}  shared H{H}", flush=True)

def q_pt(atoms):
    p = atoms.positions
    return float(np.linalg.norm(p[DONOR]-p[H]) - np.linalg.norm(p[ACC]-p[H]))

def r_oo(atoms):
    return float(np.linalg.norm(atoms.positions[DONOR] - atoms.positions[ACC]))

# The target direction: the shared proton moving along the O...O line. Cartesian; the
# mass-weighting and rigid-body projection are done inside follow_min_mode.
target = np.zeros((len(A), 3))
u = A.positions[ACC] - A.positions[DONOR]
target[H] = u / np.linalg.norm(u)

# Eight amplitudes x two REPLICATES. The second index is not a seed and is not a knob.
#
# It started life as `climb_seed`, passed to follow_min_mode's `seed` argument, on the
# reading that PRRS_STATUS.md 8.15.4 ("the seed decides which saddle the climb reaches")
# applied to it. That was a misreading: 8.15.4 is about the seed FRAME -- probe geometry and
# amplitude -- and `seed` never appears in follow_min_mode's body at all. Confirmed by
# inspecting the source: the three occurrences are a comment, a comment, and the string
# "seed_above_ceiling". Passing it different integers cannot change anything.
#
# The pairs were then expected to come back identical, and six of eight did, energies equal
# bit for bit. Two did not:
#   amp 0.50  same verdict, energies -14637.048391676 vs -14637.048391553 (1.2e-7 eV apart)
#   amp 0.10  DIFFERENT verdicts -- lost_target_mode at step 20 vs ceiling_blocked at 45
# Identical inputs, so the difference is GPU floating-point non-determinism: MACE's scatter
# reductions accumulate with atomicAdd and the order is not fixed. At amp 0.10 that ~1e-7
# noise decides which way the walk goes, which locates that amplitude on a boundary between
# two outcomes.
#
# So the repeat stays, renamed to what it actually measures: run-to-run reproducibility of
# the climb. It is a recorded observable, not a variance knob, and it is the only thing here
# that would notice a walk sitting on a knife edge. Note what it implies for any criterion
# demanding field-for-field reproduction -- P1's round7/round8 agreement holds because those
# walks were nowhere near a boundary, not because the computation is deterministic.
SEED_GRID = [(amp, rep) for amp in (0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50)
             for rep in (0, 1)]
MINE = [(k, a, rep) for k, (a, rep) in enumerate(SEED_GRID) if k % NSHARD == SHARD]
SHARDJSON = OUT / f"ts_barrier_shard{SHARD:02d}.json"
print(f"shard {SHARD}/{NSHARD}: {len(MINE)} of {len(SEED_GRID)} seeds "
      f"{[(a, rep) for _, a, rep in MINE]}", flush=True)

records = []
for kidx, amp, replicate in MINE:
    seed = A.copy()
    seed.positions[H] += amp * target[H]
    seed.calc = factory()
    t0 = time.time()
    structure, diag = follow_min_mode(
        seed, factory, cfg, direction=target,
        ceiling=E_A + cfg.saddle_search_max_rise_eV, reference=E_A)
    dt = time.time() - t0
    # Everything the diagnostics carry is kept, including on the failure path. `wrong_index`
    # in particular means the walk CONVERGED and the index was not 1 -- so saddle_order and
    # the negative eigenvalues are the whole content of that outcome, and dropping them
    # would record "no saddle" while throwing away what was actually found.
    base = dict(seed_index=kidx, seed_amp=amp, replicate=replicate,
                seconds=round(dt, 1),
                steps=diag.get("steps"), reason=diag.get("reason"),
                target_overlap=diag.get("target_overlap"),
                final_mode_overlap=diag.get("final_mode_overlap"),
                followed_order=diag.get("followed_order"),
                selected_mode_index=diag.get("selected_mode_index"),
                initial_rise_eV=diag.get("initial_rise_eV"),
                max_rise_eV=diag.get("max_rise_eV"),
                backtracks=diag.get("backtracks"),
                ceiling_hit_step=diag.get("ceiling_hit_step"),
                fmax_eV_A=diag.get("fmax_eV_A"),
                diag_saddle_order=diag.get("saddle_order"),
                diag_negative_eigenvalues=diag.get("negative_eigenvalues"),
                diag_imaginary_icm=diag.get("imaginary_wavenumbers_icm"),
                diag_floor_residual=diag.get("trivial_mode_floor_worst_residual"))
    if structure is None:
        # Not a failure of the reaction: a failure of this climb, with its reason kept.
        base["status"] = "no_saddle"
        records.append(base)
        print(f"[s{SHARD} k{kidx:02d} amp {amp:.2f} r{replicate}] no saddle: "
              f"{diag.get('reason')}  order={diag.get('saddle_order')} "
              f"lam={[round(x, 5) for x in (diag.get('negative_eigenvalues') or [])]} "
              f"icm={[round(x, 1) for x in (diag.get('imaginary_wavenumbers_icm') or [])]} "
              f"overlap={diag.get('final_mode_overlap')} "
              f"(steps {diag.get('steps')}, {dt:.0f}s)", flush=True)
        json.dump(dict(env=env, shard=SHARD, nshard=NSHARD, quench_fmax_eV_A=QUENCH, quench_steps=QUENCH_STEPS,
                       E_A=E_A, E_B=E_B, key_A=kA, key_B=kB, records=records),
                  open(SHARDJSON, "w"), indent=1, default=str)
        continue
    structure.calc = calc
    lam, vec, floor, prov = curvature_spectrum(structure, factory, cfg, source="analytic")
    neg = [i for i, x in enumerate(lam) if x < -cfg.minimum_check_eigenvalue_tol]
    E_ts = float(structure.get_potential_energy())
    base.update(status="saddle", E_eV=E_ts,
                fmax=float(np.abs(structure.get_forces()).max()),
                saddle_order=len(neg), lam1=float(lam[0]),
                floor_residual=float(floor_residual(floor)),
                curvature_source=prov["source"],
                barrier_vs_A_meV=(E_ts - E_A) * 1e3,
                barrier_vs_B_meV=(E_ts - E_B) * 1e3,
                q_PT=q_pt(structure), r_OO=r_oo(structure))
    if len(neg) == 1:
        # runner._wavenumbers is the single implementation and carries the mass-weighting
        # factor hbar*1e10/sqrt(e*amu) ~ 0.0647. An earlier version of this script wrote
        # -sqrt(-lam)/units.invcm, dropping that factor, and reported -47736 cm^-1 for
        # lam1 = -35.03 -- 15x too large and physically impossible. The correct value is
        # -3087 cm^-1, which is the same band as P1's -3186 and P2's -3010.
        base["imaginary_cm1"] = _wavenumbers([lam[0]])[0]
        mode = np.asarray(vec[:, neg[0]])          # COLUMN, mass-weighted basis
        sides = []
        for sign in (1, -1):
            endpoint, checked, relays = descend_saddle(structure, mode, factory, cfg, sign,
                                                       seed=7919 * sign)
            if endpoint is None:
                sides.append(dict(sign=sign, reached=False,
                                  reason=checked.get("reason"),
                                  saddle_order=checked.get("saddle_order"),
                                  relays=len(relays)))
                continue
            endpoint.calc = factory()
            k = chemical_key(endpoint)
            sides.append(dict(sign=sign, reached=True, relays=len(relays),
                              E_eV=float(endpoint.get_potential_energy()),
                              key=k["key"][:16], graph_hash=k["graph_hash"],
                              fragments=k["fragments"],
                              saddle_order=checked.get("saddle_order"),
                              q_PT=q_pt(endpoint), r_OO=r_oo(endpoint)))
            write(str(OUT / f"end_k{kidx:02d}_amp{amp:.2f}_r{replicate}_sign{sign:+d}.extxyz"),
                  endpoint)
        base["sides"] = sides
        got = [s["key"] for s in sides if s.get("reached")]
        base["both_sides_reached"] = len(got) == 2
        base["connects_two_nodes"] = len(set(got)) == 2
        base["endpoints_are_A_and_B"] = sorted(set(got)) == sorted({kA, kB})
        write(str(OUT / f"ts_k{kidx:02d}_amp{amp:.2f}_r{replicate}.extxyz"), structure)
        print(f"[s{SHARD} k{kidx:02d} amp {amp:.2f} r{replicate}] order 1  E {E_ts:.6f}  barrier {base['barrier_vs_A_meV']:+.1f} "
              f"meV  {base['imaginary_cm1']:.1f} cm^-1  q {base['q_PT']:+.3f}  "
              f"r_OO {base['r_OO']:.3f}  sides {got}  ({dt:.0f}s)", flush=True)
    else:
        print(f"[s{SHARD} k{kidx:02d} amp {amp:.2f} r{replicate}] saddle order {len(neg)} (not first-order)  "
              f"E {E_ts:.6f}  ({dt:.0f}s)", flush=True)
    records.append(base)
    json.dump(dict(env=env, shard=SHARD, nshard=NSHARD, quench_fmax_eV_A=QUENCH, quench_steps=QUENCH_STEPS,
                   E_A=E_A, E_B=E_B, key_A=kA, key_B=kB, records=records),
              open(SHARDJSON, "w"), indent=1, default=str)

# Per-shard summary. The distribution over seeds is global; p3_merge.py reports that.
first = [r for r in records if r.get("saddle_order") == 1]
print(f"\n{len(records)} seeds: {len(first)} first-order saddles, "
      f"{sum(1 for r in records if r['status']=='no_saddle')} no-saddle, "
      f"{sum(1 for r in records if r.get('saddle_order',1)!=1 and r['status']=='saddle')} "
      f"higher-order")
if first:
    bs = [r["barrier_vs_A_meV"] for r in first]
    ims = [r["imaginary_cm1"] for r in first]
    print(f"barrier vs A   {min(bs):.1f} .. {max(bs):.1f} meV")
    print(f"imaginary      {min(ims):.1f} .. {max(ims):.1f} cm^-1")
    print(f"both sides reached      {sum(1 for r in first if r.get('both_sides_reached'))}/{len(first)}")
    print(f"connects two nodes      {sum(1 for r in first if r.get('connects_two_nodes'))}/{len(first)}")
    print(f"endpoints are A and B   {sum(1 for r in first if r.get('endpoints_are_A_and_B'))}/{len(first)}")
print(f"wrote {OUT}/ts_barrier.json")
