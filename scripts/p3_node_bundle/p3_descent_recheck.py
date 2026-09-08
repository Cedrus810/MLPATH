"""P3 step 3: why the sign=+1 descent fails on all six OFF24 first-order saddles.

Bounded on purpose. This re-runs ONE thing -- the +1 side of the six saddles already found
in `runs/p3_ts_barrier` -- and it exists to answer one question with three allowed answers:

  1. the soft methyl rotor takes the mode on the last polish round;
  2. the descent algorithm rejects a geometry that had in fact arrived;
  3. this saddle genuinely does not connect A on the OFF24 surface.

It is NOT a tolerance scan, a step scan or a seed scan. `quench_fmax` and `quench_steps`
keep the values the archived run used (2e-3, 500) and are read from the environment only so
the record says what they were, not so they can be swept -- the 2x2 in `run_all.sh
quench-ab` already answered the budget question and is not repeated here.

Why the archived run cannot answer it. Every one of the six sides came back
`soft_mode_moved_on_the_last_round` with `descent_failure = None`: the polish rounds are
carried by `quench`'s diagnostics, and the code that keeps them through
`follow_unstable_mode` and `descend_saddle` landed after those shards were written. So the
rounds have to be produced, not mined. They are cheap: the saddles are on disk
(`ts_k*.extxyz`), so nothing re-climbs and each side is one descent.

What is measured, beyond the reason string:

  the polish rounds themselves -- per round, whether it converged, whether it MOVED, and
  which torsion moved by how much. Alternating offsets on one coordinate is a limit cycle
  between FIRE and the polish; same-sign offsets shrinking is a budget one round short.
  Both look identical from the reason string, which is why that string could not decide
  anything.

  the geometry the descent actually ended on. `follow_unstable_mode` returns None on a
  failed quench and the structure is discarded, so this script drives `quench` directly --
  the same call, three lines of setup copied from it -- and keeps `atoms` whichever way the
  verdict goes. That geometry is what separates answer 2 from answer 3: if its
  `chemical_key` is A's, the descent reached A and the refusal is bookkeeping; if it is B's,
  the +1 side runs to B as well and this saddle does not join A and B on this surface.

  whether that geometry is a minimum, and if not, where the negative mode lives. The
  methyl fraction is the same measurement §1.2 of the probes document used to disqualify
  the first candidate A (99.83% of the negative mode on atoms 0/12/13/14).

The verdict field is a mechanical reading of those recorded numbers and claims nothing
beyond them; where they do not separate the answers it says so rather than picking one.

One implementation note for the record: `internal_mode_recheck` now fails the confirmation
CLOSED when its probe is refused, where it used to skip the gate. That path is a
`GateRejected` inside the recheck, which did not occur in the archived shards, so it cannot
be the cause of anything here -- but the implementation hash in the manifest changes because
of it and this says why.
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
SADDLES = Path(ARGV[2]) if len(ARGV) > 2 else HERE.parents[1] / "runs" / "p3_ts_barrier"
OUT = Path(ARGV[3]) if len(ARGV) > 3 else HERE.parents[1] / "runs" / "p3_descent_recheck"
env = guard(MODEL, label=f"descent-recheck shard {SHARD}/{NSHARD}")
from prrs.calculators import mace_factory
from prrs.chemistry import chemical_key
from prrs.config import SearchConfig
from prrs.reliability import GuardedCalculator
from prrs.runner import (curvature_spectrum, floor_residual, quench, confirm_minimum,
                         descend_saddle, _wavenumbers)
from prrs.state import encode

OUT.mkdir(parents=True, exist_ok=True)
factory = lambda: mace_factory(MODEL, device="cuda", default_dtype="float64")
calc = factory()

# Same as the archived shards. Read from the environment so the manifest records them,
# not so they can be swept -- see the module docstring.
QUENCH = float(os.environ.get("P3_QUENCH", "0.002"))
QUENCH_STEPS = int(os.environ.get("P3_QUENCH_STEPS", "500"))
cfg = SearchConfig(quench_fmax_eV_A=QUENCH, quench_steps=QUENCH_STEPS)

runs = HERE.parents[1] / "runs"
A = read(str(runs / "p3_A_source.extxyz")); A.calc = calc
B = read(str(runs / "p3_B_source.extxyz")); B.calc = factory()
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
# The rotor §1.2 named: the methyl carbon and its three hydrogens. Fixed here rather than
# rediscovered, so "the mode is on the methyl" is a statement about the same four atoms
# the probes document made it about.
METHYL = [0, 12, 13, 14]
assert [sym[i] for i in METHYL] == ["C", "H", "H", "H"], [sym[i] for i in METHYL]

print(f"quench_fmax_eV_A = {QUENCH:.0e}   quench_steps = {QUENCH_STEPS}", flush=True)
print(f"A {kA} E={E_A:.6f}   B {kB} E={E_B:.6f}", flush=True)
print(f"donor O{DONOR}  acceptor O{ACC}  shared H{H}  methyl {METHYL}", flush=True)


def q_pt(atoms):
    p = atoms.positions
    return float(np.linalg.norm(p[DONOR]-p[H]) - np.linalg.norm(p[ACC]-p[H]))


def r_oo(atoms):
    return float(np.linalg.norm(atoms.positions[DONOR] - atoms.positions[ACC]))


def methyl_fraction(mode, masses):
    """Share of a mass-weighted mode's squared norm carried by the methyl group."""
    weighted = np.asarray(mode, dtype=float).reshape(-1, 3)
    total = float(np.sum(weighted ** 2))
    if total <= 0:
        return None
    return float(np.sum(weighted[METHYL] ** 2) / total)


def rounds_of(diagnostics):
    """The polish rounds, with the atom symbols spelled out so a movement is readable."""
    out = []
    for report in diagnostics.get("rounds") or []:
        moves = []
        for m in report.get("movements") or []:
            idx = list(m["indices"])
            moves.append({"indices": idx, "atoms": "-".join(f"{sym[i]}{i}" for i in idx),
                          "offset_rad": float(m["offset_rad"]), "kind": m.get("kind"),
                          "is_methyl": bool(set(idx) & set(METHYL))})
        out.append({"converged": report.get("converged"), "moved": report.get("moved"),
                    "movements": moves,
                    "not_within_tolerance": [list(x["indices"]) for x in report.get("modes", [])
                                             if x.get("within_tolerance") is False]})
    return out


def cycle_or_budget(rounds):
    """Limit cycle or one round short, from the coordinate that kept moving.

    Written against what the data does rather than what was expected of it. The first
    version keyed on "exactly one coordinate moved" and gave up as soon as a second one
    appeared -- and on the first saddle a phenyl torsion moved in round 0 only while the
    methyl moved in all seven, so the one coordinate that carries the failure was thrown
    away by a coordinate that had nothing to do with it. What matters is the coordinate
    present in EVERY round that moved; the rest are transient and are reported as such.

    Same-sign offsets shrinking geometrically is a budget: the polish and FIRE are
    converging on each other and the round count ran out. Alternating signs is a limit
    cycle, which no round count fixes. The extrapolated round count is arithmetic on the
    measured decay, reported so "short by how much" is a number rather than a feeling --
    it is not a recommendation to raise the knob.
    """
    moving = [r for r in rounds if r["movements"]]
    if len(moving) < 2:
        return {"verdict": "too_few_moving_rounds", "coordinate": None, "offsets": [],
                "transient": []}
    per_round = [{tuple(m["indices"]) for m in r["movements"]} for r in moving]
    persistent = sorted(set.intersection(*per_round))
    transient = sorted(set.union(*per_round) - set(persistent))
    if len(persistent) != 1:
        return {"verdict": "no_single_persistent_coordinate", "coordinate": persistent,
                "offsets": [], "transient": [list(t) for t in transient]}
    coordinate = persistent[0]
    offsets = [m["offset_rad"] for r in moving for m in r["movements"]
               if tuple(m["indices"]) == coordinate]
    signs = {int(np.sign(o)) for o in offsets if o != 0}
    shrinking = all(abs(b) < abs(a) for a, b in zip(offsets, offsets[1:]))
    ratio, remaining = None, None
    if len(signs) == 1 and shrinking and abs(offsets[0]) > 0:
        # Geometric decay over the rounds actually observed, then how many more of them
        # the last offset needs to fall under the torsion tolerance.
        ratio = float(abs(offsets[-1] / offsets[0]) ** (1.0 / (len(offsets) - 1)))
        if 0 < ratio < 1 and abs(offsets[-1]) > cfg.torsion_tolerance_rad:
            remaining = float(np.log(cfg.torsion_tolerance_rad / abs(offsets[-1]))
                              / np.log(ratio))
    if len(signs) > 1:
        verdict = "limit_cycle"          # FIRE and the polish undoing each other
    elif shrinking:
        verdict = "budget_short"         # converging, the round count ran out
    else:
        verdict = "same_sign_not_shrinking"
    return {"verdict": verdict, "coordinate": list(coordinate), "offsets": offsets,
            "transient": [list(t) for t in transient],
            "is_methyl": bool(set(coordinate) & set(METHYL)),
            "decay_per_round": ratio,
            "polish_rounds_used": cfg.soft_polish_rounds,
            "extra_rounds_to_tolerance": remaining,
            "tolerance_rad": cfg.torsion_tolerance_rad}


saddles = sorted(SADDLES.glob("ts_k*.extxyz"))
if not saddles:
    raise SystemExit(f"no ts_k*.extxyz in {SADDLES}; run `run_all.sh ts` first")
MINE = [(k, p) for k, p in enumerate(saddles) if k % NSHARD == SHARD]
print(f"shard {SHARD}/{NSHARD}: {len(MINE)} of {len(saddles)} saddles "
      f"{[p.name for _, p in MINE]}", flush=True)
SHARDJSON = OUT / f"descent_recheck_shard{SHARD:02d}.json"

records = []
for kidx, path in MINE:
    t0 = time.time()
    ts = read(str(path))
    ts.calc = calc
    E_ts = float(ts.get_potential_energy())
    lam, vec, floor, prov = curvature_spectrum(ts, factory, cfg, source="analytic")
    neg = [i for i, x in enumerate(lam) if x < -cfg.minimum_check_eigenvalue_tol]
    rec = {"saddle_file": path.name, "index": kidx,
           "E_ts_eV": E_ts, "barrier_vs_A_meV": (E_ts - E_A) * 1e3,
           "saddle_order": len(neg), "lam1": float(lam[0]),
           "imaginary_cm1": _wavenumbers([lam[0]])[0] if len(neg) == 1 else None,
           "floor_residual": float(floor_residual(floor)),
           "curvature_source": prov["source"],
           "ts_q_PT": q_pt(ts), "ts_r_OO": r_oo(ts)}
    if len(neg) != 1:
        # The archive says these are first-order saddles. If one is not, that is the
        # finding and nothing downstream of it means anything.
        rec.update(status="not_a_first_order_saddle", seconds=round(time.time()-t0, 1))
        records.append(rec)
        print(f"[k{kidx:02d} {path.name}] order {len(neg)} -- NOT first order, skipped",
              flush=True)
        json.dump(dict(env=env, shard=SHARD, nshard=NSHARD, quench_fmax_eV_A=QUENCH,
                       quench_steps=QUENCH_STEPS, E_A=E_A, E_B=E_B, key_A=kA, key_B=kB,
                       sign=1, records=records), open(SHARDJSON, "w"), indent=1, default=str)
        continue

    # Hop 0 of descend_saddle(sign=+1), written out so the geometry survives a failure.
    # Three lines lifted from follow_unstable_mode: mass-weighted eigenvector column as the
    # displacement, scaled so the largest single-atom step is irc_step_A. The archived
    # failures all have relays=0, i.e. none of them got past this hop, so one hop is the
    # whole of what happened.
    sign = 1
    mode = np.asarray(vec[:, 0], dtype=float).reshape(-1, 3)
    extent = np.linalg.norm(mode, axis=1).max()
    atoms = ts.copy()
    guardcalc = GuardedCalculator(factory(), cfg, E_ts)
    atoms.calc = guardcalc
    atoms.set_positions(ts.positions + sign * cfg.irc_step_A / extent * mode)
    atoms.set_momenta(np.zeros((len(atoms), 3)))
    converged, diagnostics = quench(atoms, guardcalc, cfg)

    # The endpoint the archive threw away, kept whatever the verdict was.
    final = atoms.copy()
    final.calc = factory()
    E_final = float(final.get_potential_energy())
    fmax = float(np.abs(final.get_forces()).max())
    key = chemical_key(final)["key"][:16]
    settled, checked = confirm_minimum(final, factory, cfg, 7919)
    fneg = checked.get("negative_eigenvalues") or []
    fmode = checked.get("unstable_mode")
    rec.update(
        status="descended", seconds=round(time.time()-t0, 1),
        quench_converged=bool(converged), reason=diagnostics.get("reason"),
        quench_steps_used=diagnostics.get("steps"),
        rounds=rounds_of(diagnostics), motion=cycle_or_budget(rounds_of(diagnostics)),
        final_E_eV=E_final, final_vs_A_meV=(E_final - E_A) * 1e3,
        final_vs_B_meV=(E_final - E_B) * 1e3, final_vs_TS_meV=(E_final - E_ts) * 1e3,
        final_fmax_eV_A=fmax, final_key=key,
        final_geometry_note=("quench returns the geometry AFTER its last polish "
                             "displacement with no minimization following it, so fmax "
                             "above quench_fmax_eV_A is expected here and is not evidence "
                             "the descent failed to relax"),
        final_is_A=(key == kA), final_is_B=(key == kB),
        final_q_PT=q_pt(final), final_r_OO=r_oo(final),
        final_confirmed_minimum=bool(settled),
        final_confirm_reason=checked.get("reason"),
        final_saddle_order=checked.get("saddle_order"),
        final_negative_eigenvalues=fneg,
        final_recheck_performed=(checked.get("internal_mode_recheck") or {}).get("performed"),
        final_recheck_agree=(checked.get("internal_mode_recheck") or {}).get("agree"),
        final_recheck_blocked=(checked.get("internal_mode_recheck") or {}).get("blocked"),
        final_negative_mode_methyl_fraction=(
            None if fmode is None else methyl_fraction(fmode, final.get_masses())))

    # The production path, run second and recorded separately. The block above drives
    # `quench` directly because that is the only way to see the geometry a FAILED descent
    # produced; this is what `search.py` itself calls, so it is the field that says whether
    # the two-sided connection is named by the search rather than reconstructed here. Once
    # `follow_unstable_mode` admits an endpoint whose only unfinished business is a
    # contracting soft coordinate, the two should agree -- and if they do not, that
    # disagreement is the finding.
    produced, checked_p, relays = descend_saddle(ts, np.asarray(vec[:, 0]), factory, cfg,
                                                 sign, seed=7919 * sign)
    if produced is None:
        rec["production_path"] = {"reached": False, "reason": checked_p.get("reason"),
                                  "saddle_order": checked_p.get("saddle_order"),
                                  "relays": len(relays)}
    else:
        produced.calc = factory()
        info = (produced.info or {}).get("quench") or {}
        rec["production_path"] = {
            "reached": True, "relays": len(relays),
            "key": chemical_key(produced)["key"][:16],
            "is_A": chemical_key(produced)["key"][:16] == kA,
            "E_eV": float(produced.get_potential_energy()),
            "fmax_eV_A": float(np.abs(produced.get_forces()).max()),
            "quench_converged": info.get("converged"),
            "admitted_with_unsettled_soft_modes": info.get(
                "admitted_with_unsettled_soft_modes", False),
            "polish_contraction_verdict": (info.get("polish_contraction") or {}).get("verdict"),
            "confirmed_minimum": bool(checked_p.get("saddle_order") == 0),
            "saddle_order": checked_p.get("saddle_order"),
            "q_PT": q_pt(produced), "r_OO": r_oo(produced)}
        write(str(OUT / f"prod_plus_{path.stem}.extxyz"), produced)
    print(f"    [production descend_saddle] {json.dumps(rec['production_path'], default=str)}",
          flush=True)

    # Mechanical, from the fields above and nothing else.
    if key == kB:
        verdict = "plus_side_reaches_B_too"          # answer 3: this saddle does not join A
    elif key != kA:
        verdict = "plus_side_reaches_a_third_basin"  # answer 3, and the basin is named
    elif settled:
        verdict = "reached_A_and_confirms"           # answer 2 outright: the refusal was wrong
    elif (rec["final_negative_mode_methyl_fraction"] or 0) > 0.9:
        verdict = "reached_A_rotor_holds_the_mode"   # answer 1
    elif not fneg:
        verdict = "reached_A_bookkeeping_refusal"    # answer 2: a minimum by count, refused
    else:
        verdict = "reached_A_but_not_a_minimum"      # unresolved; do not force it
    rec["verdict"] = verdict
    records.append(rec)

    write(str(OUT / f"end_plus_{path.stem}.extxyz"), final)
    print(f"[k{kidx:02d} {path.name}] {verdict}", flush=True)
    print(f"    quench converged={converged} reason={diagnostics.get('reason')} "
          f"steps={diagnostics.get('steps')}  motion={rec['motion']['verdict']} "
          f"coord={rec['motion']['coordinate']} offsets="
          f"{[round(o, 5) for o in rec['motion']['offsets']]}", flush=True)
    print(f"    final key={key} ({'A' if key == kA else 'B' if key == kB else 'other'}) "
          f"E-A={rec['final_vs_A_meV']:+.1f} meV  fmax={fmax:.2e}  q_PT={rec['final_q_PT']:+.3f} "
          f"r_OO={rec['final_r_OO']:.3f}", flush=True)
    print(f"    minimum={settled} order={checked.get('saddle_order')} neg={[round(x, 6) for x in fneg]} "
          f"methyl_fraction={rec['final_negative_mode_methyl_fraction']}", flush=True)
    for n, rnd in enumerate(rec["rounds"]):
        moves = ", ".join(f"{m['atoms']}{m['offset_rad']:+.5f}({m['kind']})"
                          for m in rnd["movements"]) or "-"
        print(f"    [round {n}] converged={rnd['converged']} moved={rnd['moved']}  {moves}",
              flush=True)

    json.dump(dict(env=env, shard=SHARD, nshard=NSHARD, quench_fmax_eV_A=QUENCH,
                   quench_steps=QUENCH_STEPS, E_A=E_A, E_B=E_B, key_A=kA, key_B=kB,
                   sign=1, records=records), open(SHARDJSON, "w"), indent=1, default=str)

print(f"shard {SHARD} done: {len(records)} saddles -> {SHARDJSON}", flush=True)
