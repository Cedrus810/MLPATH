"""HEAVY 4. The same reaction under three independently trained models.

This is the one thing the OFF24-only numbers cannot supply. `+316 meV` and `-3086 cm^-1`
are MACE-OFF24's values on MACE-OFF24's surface; nothing in P0-P3 so far distinguishes
"the method is self-consistent" from "the surface is right". Three models trained on
different data, compared on RELATIVE quantities, is a real discriminant -- the same move
the SN2 probe made, where OMOL-0 and POLAR-1 agreed on the ion-dipole complex to
0.34 kcal/mol and disagreed at the five-coordinate saddle by 8.71.

What is comparable and what is not
----------------------------------
ABSOLUTE energies are NOT comparable. Different training sets, different references,
different isolated-atom baselines. Any table putting -14637 eV next to another model's
total is meaningless. Only differences within one model are reported:

    dE(B - A)          the tautomer asymmetry
    barrier            E_TS - E_A
    imaginary          from that model's own analytic Hessian
    geometry           r_OO at A and at the TS, q_PT at the TS, and the gating contraction
    saddle order       from that model's own Hessian
    chemical_key       a graph invariant, so it SHOULD agree; if it does not, the models
                       disagree about connectivity, which is itself the finding

Each model relaxes its own endpoints. Handing model B a geometry that is a minimum of
model A would measure the geometry transfer, not the reaction -- so A, B and the TS are
re-derived per model from the same starting geometry, and every endpoint gets
confirm_minimum rather than the optimizer's fmax.

Charge-aware models
-------------------
omol-0 and POLAR-1 read atoms.info["charge"] and ["spin"]. Benzoylacetone is neutral
closed-shell so charge=0, spin=1; per HANDOFF.md 5 the spin convention is NOT verified
(on omol-0, H2 at spin=3 came out 0.9 meV BELOW spin=1, when it should be ~10 eV above),
which is exactly why nothing open-shell is touched here.

POLAR-1 also emits a dipole and per-atom charges (Qs). Those are recorded along the
transfer because they are information no other model in this set can give: where the charge
sits while the proton is in flight. Recorded, not interpreted -- HANDOFF.md 5 leaves the
mechanism of POLAR-1's long-range term unexplained (it grows as ~1/r^2.6-3.3, not 1/r^2),
and for an INTRAmolecular transfer with no fragment separation that term should largely
cancel between A and B, but r_OO does move 2.53 -> 2.32 A, so "largely" is an assumption
and not a measurement.

Sharding is by model: each shard relaxes its endpoints once and then walks all amplitudes,
so the endpoint work is not repeated per amplitude.
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
OUT = Path(ARGV[1]) if len(ARGV) > 1 else HERE.parents[1] / "runs" / "p3_three_model"
QUENCH = float(os.environ.get("P3_QUENCH", "0.0002"))
AMPS = [float(x) for x in os.environ.get("P3_AMPS", "0.15,0.25,0.35,0.50").split(",")]

MINE = [m for k, m in enumerate(CANDIDATES) if k % NSHARD == SHARD]
env = guard(CANDIDATES[0], label=f"three-model shard {SHARD}/{NSHARD}", models=CANDIDATES)

import numpy as np
from ase.io import read, write
from ase.optimize import LBFGS
from mace.calculators import MACECalculator
from prrs.chemistry import chemical_key, classify_transition, reaction_event_key
from prrs.config import SearchConfig
from prrs.runner import (curvature_spectrum, floor_residual, follow_min_mode,
                         descend_saddle, _wavenumbers)
from prrs.state import encode

OUT.mkdir(parents=True, exist_ok=True)
cfg = SearchConfig(quench_fmax_eV_A=QUENCH)
print(f"quench_fmax_eV_A = {QUENCH:.0e}   amplitudes = {AMPS}", flush=True)

runs = HERE.parents[1] / "runs"
seed_geom = None
for name in ("p3_A_source_tight.extxyz", "p3_A_source.extxyz", "p3_A_chelated.extxyz"):
    if (runs / name).exists():
        seed_geom = read(str(runs / name)); print(f"starting geometry: {name}", flush=True); break
if seed_geom is None:
    raise SystemExit("no starting geometry in runs/")
SYM = list(seed_geom.symbols)
OXY = [i for i, s in enumerate(SYM) if s == "O"]
HYD = [i for i, s in enumerate(SYM) if s == "H"]

# The projection coordinate is a DEFINITION and is fixed once, here, from the shared
# starting geometry. It must not be re-derived per structure.
#
# It was, and the bug was silent and total: settle() looked for "the O-H oxygen" in each
# structure and called that the donor. In B the proton has already moved to the other
# oxygen, so B's donor was A's acceptor, and q_PT came out ~-0.6 for BOTH endpoints. On the
# figure B's marker landed exactly on A's. A transfer coordinate is only a progress
# variable if its indices are held: negative on one side, zero where the atom is shared,
# positive on the other (see internal.named_values' own docstring). Re-deriving them makes
# it measure "how far is the proton from wherever it currently sits", which is ~-0.6
# everywhere and says nothing.
_e = sorted(tuple(sorted(e)) for e in encode(seed_geom).edges)
_oh = [(o, h) for o in OXY for h in HYD if tuple(sorted((o, h))) in _e]
assert len(_oh) == 1, f"expected one O-H in the starting geometry, got {_oh}"
DONOR_REF, H_REF = _oh[0]
ACC_REF = [o for o in OXY if o != DONOR_REF][0]
print(f"projection fixed: q_PT = r(O{DONOR_REF}-H{H_REF}) - r(O{ACC_REF}-H{H_REF}), "
      f"r_OO = r(O{DONOR_REF}-O{ACC_REF})", flush=True)


def tag_charge(atoms, info):
    """Set the total charge and spin for the charge-aware models; neutral closed shell."""
    out = atoms.copy()
    if info["charge_aware"]:
        out.info["charge"], out.info["spin"] = 0, 1
    return out


def oh_bond(atoms):
    edges = sorted(tuple(sorted(e)) for e in encode(atoms).edges)
    oh = [(o, h) for o in OXY for h in HYD if tuple(sorted((o, h))) in edges]
    return oh


def coords(atoms, donor, acc, h):
    p = atoms.positions
    return (float(np.linalg.norm(p[donor] - p[h]) - np.linalg.norm(p[acc] - p[h])),
            float(np.linalg.norm(p[donor] - p[acc])))


def polar_extras(calc, atoms):
    """POLAR-1 only: dipole and per-atom charges, recorded not interpreted."""
    out = {}
    try:
        d = calc.results.get("dipole")
        if d is not None:
            out["dipole_D"] = [float(x) for x in np.asarray(d).ravel()[:3]]
    except Exception:
        pass
    try:
        qs = atoms.arrays.get("Qs")
        if qs is not None:
            q = np.asarray(qs).ravel()
            out["Qs_sum"] = float(q.sum())
            out["Qs_by_element"] = {e: float(q[[i for i, s in enumerate(SYM) if s == e]].sum())
                                    for e in sorted(set(SYM))}
            out["Q_shared_H"] = None    # filled by the caller which knows the index
            out["Qs"] = [round(float(x), 4) for x in q]
    except Exception:
        pass
    return out


results = {}
for path in MINE:
    info = model_info(path)
    name = info["name"]
    print(f"\n{'='*70}\n{name}  ({info['model_type']}, {info['sha256_16']})\n{'='*70}",
          flush=True)
    t_model = time.time()
    calc = MACECalculator(model_paths=str(path), model_type=info["model_type"],
                          device="cuda", default_dtype="float64")
    factory = lambda c=calc: c
    rec = {"model": name, "sha256_16": info["sha256_16"], "model_type": info["model_type"],
           "charge_aware": info["charge_aware"], "atomic_charges": info["atomic_charges"],
           "quench_fmax_eV_A": QUENCH, "note": info["note"]}

    def settle(atoms, label):
        x = tag_charge(atoms, info); x.calc = calc
        t0 = time.time()
        opt = LBFGS(x, logfile=None); opt.run(fmax=QUENCH, steps=20000)
        lam, vec, fl, pr = curvature_spectrum(x, factory, cfg, source="analytic")
        neg = [float(v) for v in lam if v < -cfg.minimum_check_eigenvalue_tol]
        internal = [float(v) for v in lam if abs(v) > cfg.minimum_check_eigenvalue_tol]
        k = chemical_key(x)
        d = dict(label=label, E_eV=float(x.get_potential_energy()),
                 fmax=float(np.abs(x.get_forces()).max()),
                 relax_steps=int(opt.get_number_of_steps()),
                 seconds=round(time.time() - t0, 1),
                 n_minus=len(neg), lowest_internal=internal[0] if internal else None,
                 floor_residual=float(floor_residual(fl)),
                 key=k["key"][:16], graph_hash=k["graph_hash"], fragments=k["fragments"],
                 stereo_unresolved=k["stereo_unresolved"],
                 n_locked_parity=len(k["locked_bond_parity"]))
        oh = oh_bond(x)
        # Reported on the FIXED reference indices, so A and B land on opposite sides.
        q, r = coords(x, DONOR_REF, ACC_REF, H_REF)
        d.update(donor=DONOR_REF, acceptor=ACC_REF, shared_H=H_REF, q_PT=q, r_OO=r,
                 projection="fixed reference indices from the shared starting geometry")
        # Which oxygen actually carries the proton in THIS structure is separate
        # information and is kept as such rather than redefining the coordinate.
        if len(oh) == 1:
            d["protonated_oxygen"] = oh[0][0]
        d["n_OH"] = len(oh)
        if info["atomic_charges"]:
            ex = polar_extras(calc, x)
            if "Qs" in ex and len(oh) == 1:
                ex["Q_shared_H"] = ex["Qs"][oh[0][1]]
            d["polar"] = ex
        print(f"  [{label}] E={d['E_eV']:.6f} fmax={d['fmax']:.1e} n-={d['n_minus']} "
              f"steps={d['relax_steps']} key={d['key']} "
              f"q_PT={d.get('q_PT', float('nan')):+.3f} r_OO={d.get('r_OO', float('nan')):.4f} "
              f"({d['seconds']}s)", flush=True)
        return x, d

    # --- A: this model's own minimum, from the shared starting geometry ---------------
    A, dA = settle(seed_geom, "A")
    rec["A"] = dA
    if dA["n_minus"] != 0 or dA["n_OH"] != 1:
        rec["status"] = "A_not_a_usable_minimum"
        results[name] = rec
        print(f"  {name}: A is not a confirmed single-OH minimum; stopping this model. "
              f"Recorded, not skipped.", flush=True)
        json.dump(dict(env=env, models=results),
                  open(OUT / f"three_model_shard{SHARD:02d}.json", "w"), indent=1, default=str)
        continue
    DONOR, H, ACC = DONOR_REF, H_REF, ACC_REF

    # --- B: proton moved across this model's own chelate ------------------------------
    B0 = A.copy()
    u = A.positions[DONOR] - A.positions[ACC]
    B0.positions[H] = A.positions[ACC] + 0.99 * u / np.linalg.norm(u)
    B, dB = settle(B0, "B")
    rec["B"] = dB
    rec["dE_B_minus_A_meV"] = (dB["E_eV"] - dA["E_eV"]) * 1e3
    same_key = dA["key"] == dB["key"]
    rec["A_B_keys_differ"] = not same_key
    if not same_key:
        cls = classify_transition(A, B, cfg.bond_scale, cfg.active_atoms)
        evf = reaction_event_key(A, B, cfg.bond_scale, cfg.active_atoms)
        evr = reaction_event_key(B, A, cfg.bond_scale, cfg.active_atoms)
        rec.update(classify=str(cls), event_forward=evf["key"], event_reverse=evr["key"],
                   direction_quotiented=evf["key"] == evr["key"],
                   broken=evf["broken"], formed=evf["formed"])
        print(f"  event {evf['key']}  broken={evf['broken']} formed={evf['formed']}  "
              f"dE(B-A)={rec['dE_B_minus_A_meV']:+.3f} meV", flush=True)
    else:
        print(f"  WARNING: this model puts A and B in the SAME chemical node "
              f"(key {dA['key']}). The transfer is then not a cross-node reaction under "
              f"this model, which is a disagreement about connectivity, not a bug.",
              flush=True)

    # --- TS: min-mode following on the transfer coordinate ----------------------------
    target = np.zeros((len(A), 3))
    target[H] = u / np.linalg.norm(u) * -1.0        # from donor toward acceptor
    walks = []
    for amp in AMPS:
        s0 = A.copy()
        s0.positions[H] += amp * target[H]
        s0 = tag_charge(s0, info); s0.calc = calc
        t0 = time.time()
        st, diag = follow_min_mode(s0, factory, cfg, direction=target,
                                   ceiling=dA["E_eV"] + cfg.saddle_search_max_rise_eV,
                                   reference=dA["E_eV"])
        dt = time.time() - t0
        w = dict(amplitude=amp, seconds=round(dt, 1), steps=diag.get("steps"),
                 reason=diag.get("reason"), target_overlap=diag.get("target_overlap"),
                 final_mode_overlap=diag.get("final_mode_overlap"),
                 diag_saddle_order=diag.get("saddle_order"),
                 diag_negative_eigenvalues=diag.get("negative_eigenvalues"),
                 initial_rise_eV=diag.get("initial_rise_eV"),
                 max_rise_eV=diag.get("max_rise_eV"))
        if st is None:
            w["status"] = "no_saddle"
            print(f"  [TS amp {amp:.2f}] no saddle: {diag.get('reason')} "
                  f"order={diag.get('saddle_order')} (steps {diag.get('steps')}, {dt:.0f}s)",
                  flush=True)
            walks.append(w); continue
        st.calc = calc
        lam, vec, fl, pr = curvature_spectrum(st, factory, cfg, source="analytic")
        neg = [i for i, v in enumerate(lam) if v < -cfg.minimum_check_eigenvalue_tol]
        q, r = coords(st, DONOR, ACC, H)
        w.update(status="saddle", E_eV=float(st.get_potential_energy()),
                 saddle_order=len(neg), lam1=float(lam[0]),
                 floor_residual=float(floor_residual(fl)),
                 barrier_vs_A_meV=(float(st.get_potential_energy()) - dA["E_eV"]) * 1e3,
                 q_PT=q, r_OO=r, gating_contraction_A=dA["r_OO"] - r)
        if len(neg) == 1:
            w["imaginary_cm1"] = _wavenumbers([lam[0]])[0]
            if info["atomic_charges"]:
                ex = polar_extras(calc, st)
                if "Qs" in ex:
                    ex["Q_shared_H"] = ex["Qs"][H]
                w["polar"] = ex
            sides = []
            for sign in (1, -1):
                ep, checked, relays = descend_saddle(st, np.asarray(vec[:, neg[0]]),
                                                     factory, cfg, sign, seed=7919 * sign)
                if ep is None:
                    sides.append(dict(sign=sign, reached=False,
                                      reason=checked.get("reason"),
                                      saddle_order=checked.get("saddle_order"),
                                      relays=len(relays)))
                    continue
                ep.calc = calc
                kk = chemical_key(ep)
                sides.append(dict(sign=sign, reached=True, relays=len(relays),
                                  E_eV=float(ep.get_potential_energy()),
                                  key=kk["key"][:16], graph_hash=kk["graph_hash"],
                                  q_PT=coords(ep, DONOR, ACC, H)[0],
                                  r_OO=coords(ep, DONOR, ACC, H)[1]))
            w["sides"] = sides
            got = [x["key"] for x in sides if x.get("reached")]
            w["both_sides_reached"] = len(got) == 2
            w["connects_two_nodes"] = len(set(got)) == 2
            w["endpoints_are_A_and_B"] = sorted(set(got)) == sorted({dA["key"], dB["key"]})
            write(str(OUT / f"{name}_ts_amp{amp:.2f}.extxyz"), st)
            print(f"  [TS amp {amp:.2f}] order 1  barrier {w['barrier_vs_A_meV']:+.1f} meV  "
                  f"{w['imaginary_cm1']:.1f} cm^-1  q {q:+.3f}  r_OO {r:.4f}  "
                  f"gate -{w['gating_contraction_A']:.3f} A  sides {got}  ({dt:.0f}s)",
                  flush=True)
        else:
            print(f"  [TS amp {amp:.2f}] order {len(neg)} (not first-order) "
                  f"E={w['E_eV']:.6f} ({dt:.0f}s)", flush=True)
        walks.append(w)
    rec["walks"] = walks
    first = [w for w in walks if w.get("saddle_order") == 1]
    if first:
        rec["barrier_meV"] = [w["barrier_vs_A_meV"] for w in first]
        rec["imaginary_cm1"] = [w["imaginary_cm1"] for w in first]
        rec["gating_contraction_A"] = [w["gating_contraction_A"] for w in first]
    rec["status"] = "done"
    rec["model_seconds"] = round(time.time() - t_model, 1)
    write(str(OUT / f"{name}_A.extxyz"), A)
    write(str(OUT / f"{name}_B.extxyz"), B)
    results[name] = rec
    json.dump(dict(env=env, models=results),
              open(OUT / f"three_model_shard{SHARD:02d}.json", "w"), indent=1, default=str)
    print(f"{name} finished in {rec['model_seconds']}s", flush=True)

json.dump(dict(env=env, models=results),
          open(OUT / f"three_model_shard{SHARD:02d}.json", "w"), indent=1, default=str)
print(f"\nshard {SHARD} wrote three_model_shard{SHARD:02d}.json  "
      f"(models: {list(results)})", flush=True)
print("ABSOLUTE energies across models are NOT comparable; only the within-model "
      "differences above are.", flush=True)
