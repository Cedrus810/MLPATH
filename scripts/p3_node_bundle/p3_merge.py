"""Assemble sharded P3 results and answer the questions that are global, not per-shard.

Kept separate from the workers on purpose. Each worker can only see its own stride, so
"how many distinct confirmed minima are there", "what is the seed distribution" and "what
does the surface look like" are not answerable inside a shard. A worker that printed its
own count would be printing a number that reads like the answer and is not.

Also reports what is MISSING. A merge over 7 of 8 shards is not the experiment, and a
silent merge would look identical to a complete one -- so incomplete input is stated in
the output and the exit code is nonzero.

Usage: python p3_merge.py [runs_dir]
"""
import json, sys, glob
from pathlib import Path
import numpy as np

ROOT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parents[2] / "runs"
incomplete = []


def shards(pattern):
    return sorted(glob.glob(str(pattern)))


def check(name, found, expected):
    if expected is None:
        return
    have = {int(f) for f in found}
    missing = sorted(set(range(expected)) - have)
    if missing:
        incomplete.append(f"{name}: missing shards {missing} of {expected}")
        print(f"  INCOMPLETE -- missing shards {missing} of {expected}")


# ---------------------------------------------------------------- TS / barrier
d = ROOT / "p3_ts_barrier"
files = shards(d / "ts_barrier_shard*.json")
if files:
    print(f"\n=== TS / barrier ({len(files)} shards) ===")
    recs, meta, ns, tols = [], None, None, set()
    for f in files:
        blob = json.load(open(f))
        recs += blob["records"]; meta = meta or blob; ns = blob.get("nshard")
        tols.add((blob.get("quench_fmax_eV_A"), blob.get("quench_steps")))
    # A mixed-tolerance merge would silently average two different experiments. The same
    # conflation voided P2's round 5 batch, so it is refused rather than warned about.
    if len(tols) > 1:
        print(f"  REFUSING to merge: shards carry different (quench_fmax, quench_steps) {sorted(tols, key=str)}")
        incomplete.append(f"ts_barrier: mixed quench settings {sorted(tols, key=str)}")
        files = []
    else:
        print(f"  (quench_fmax_eV_A, quench_steps) = {tols.pop()}")
    check("ts_barrier", [json.load(open(f))["shard"] for f in files], ns)
    recs.sort(key=lambda r: r.get("seed_index", -1))
    first = [r for r in recs if r.get("saddle_order") == 1]
    nosad = [r for r in recs if r.get("status") == "no_saddle"]
    higher = [r for r in recs if r.get("status") == "saddle" and r.get("saddle_order") != 1]
    print(f"{len(recs)} seeds: {len(first)} first-order, {len(nosad)} no-saddle, "
          f"{len(higher)} higher-order")
    # The three outcomes are reported as three outcomes. Per the project rule, an
    # unconverged climb is recorded as unconverged, not as a failure of the reaction.
    from collections import Counter
    if nosad:
        print("  no-saddle reasons: " + json.dumps(Counter(r.get("reason") for r in nosad)))
    if first:
        bs = np.array([r["barrier_vs_A_meV"] for r in first])
        ims = np.array([r["imaginary_cm1"] for r in first])
        qs = np.array([r["q_PT"] for r in first])
        rs = np.array([r["r_OO"] for r in first])
        print(f"  barrier vs A   {bs.min():.1f} .. {bs.max():.1f} meV   "
              f"(spread {bs.max()-bs.min():.1f}, n={len(bs)})")
        print(f"  imaginary      {ims.min():.1f} .. {ims.max():.1f} cm^-1")
        print(f"  q_PT at TS     {qs.min():+.3f} .. {qs.max():+.3f}   "
              f"r_OO {rs.min():.3f} .. {rs.max():.3f} A")
        b = sum(1 for r in first if r.get("both_sides_reached"))
        c = sum(1 for r in first if r.get("connects_two_nodes"))
        e = sum(1 for r in first if r.get("endpoints_are_A_and_B"))
        print(f"  both sides reached      {b}/{len(first)}")
        print(f"  connects two nodes      {c}/{len(first)}")
        print(f"  endpoints are A and B   {e}/{len(first)}")
        # distinct saddles by barrier, at the 0.81 meV endpoint spread P2 measured
        lvl = sorted({round(x / 0.81) for x in bs})
        print(f"  distinct barriers at P2's 0.81 meV endpoint spread: {len(lvl)}")
    json.dump(dict(source_shards=files, records=recs, E_A=meta["E_A"], E_B=meta["E_B"],
                   key_A=meta["key_A"], key_B=meta["key_B"]),
              open(d / "ts_barrier.json", "w"), indent=1, default=str)
    print(f"  merged -> {d}/ts_barrier.json")

# ---------------------------------------------------------------- landscape
d = ROOT / "p3_landscape_full"
files = shards(d / "landscape_shard*.json")
if files:
    print(f"\n=== conformer landscape ({len(files)} shards) ===")
    rows, meta, ns = [], None, None
    for f in files:
        blob = json.load(open(f))
        rows += blob["rows"]; meta = meta or blob; ns = blob.get("nshard")
    check("landscape", [json.load(open(f))["shard"] for f in files], ns)
    rows.sort(key=lambda r: r["conf"])
    conf = [r for r in rows if r["n_minus"] == 0]
    closed = [r for r in conf if r["r_OO"] < 3.0]
    openf = [r for r in conf if r["r_OO"] >= 3.0]
    print(f"{len(rows)} relaxed, {len(conf)} Hessian-confirmed minima "
          f"({len(rows)-len(conf)} had n_minus > 0 and are NOT minima)")
    print(f"  closed chelate (r_OO < 3.0)  {len(closed)}")
    print(f"  open form                    {len(openf)}")
    if closed:
        E = np.array([r["E_eV"] for r in closed])
        # 0.81 meV is the endpoint spread P2 measured for this protocol. It is NOT an
        # energy resolution -- becoming one needs every endpoint re-relaxed at one tight
        # tolerance, which was never done (P2 freeze, non-claim 4). Used here only as the
        # scale below which two energies are not being called distinct.
        print(f"  distinct at 0.1 meV   {len({round(e,4) for e in E})}")
        print(f"  distinct at 0.81 meV  {len({round(e/8.1e-4) for e in E})}"
              f"   <- the defensible count")
        print(f"  distinct chemical_key {len({r['key'] for r in closed})}")
        print(f"  E span {(E.max()-E.min())*1e3:.2f} meV, lowest {E.min():.6f} eV")
        if openf:
            Eo = np.array([r["E_eV"] for r in openf])
            print(f"  chelate stabilisation {(Eo.min()-E.min())*1e3:.1f} meV")
        n_needed = len({round(e/8.1e-4) for e in E})
        print(f"  => conformer_max_per_node for P3 must be >= {n_needed} "
              f"(P2 ran 8); reservoir would evict below that")
    json.dump(dict(source_shards=files, smiles=meta["smiles"],
                   atom_order=meta["atom_order"], rows=rows),
              open(d / "landscape.json", "w"), indent=1, default=str)
    print(f"  merged -> {d}/landscape.json")

# ---------------------------------------------------------------- 2D surface
d = ROOT / "p3_pes2d"
for mode in ("rigid", "relaxed"):
    files = shards(d / f"pes2d_{mode}_shard*.npz")
    if not files:
        continue
    print(f"\n=== 2D surface, {mode} ({len(files)} shards) ===")
    blobs = [np.load(f) for f in files]
    Q, R = blobs[0]["Q"], blobs[0]["R"]
    E = np.full((len(R), len(Q)), np.nan)
    calls = np.zeros((len(R), len(Q)), dtype=int)
    ach = np.full((len(R), len(Q), 2), np.nan)
    conv = np.full((len(R), len(Q)), -1, dtype=np.int8)
    for b in blobs:
        for i in b["rows"]:
            E[i] = b["E"][i]; calls[i] = b["calls"][i]; ach[i] = b["achieved"][i]
            if "converged" in b:
                conv[i] = b["converged"][i]
    check(f"pes2d_{mode}", [int(b["shard"]) for b in blobs],
          int(blobs[0]["nshard"]) if "nshard" in blobs[0] else None)
    ok = np.isfinite(E)
    print(f"  {ok.sum()}/{E.size} points filled, {E.size-ok.sum()} still NaN")
    # A point that ran out of optimizer steps is not a value of the relaxed surface.
    # Reported separately rather than folded into "filled", because a panel that mixes
    # converged and unconverged points is not a relaxed surface.
    if (conv >= 0).any():
        nc = int((conv == 0).sum()); yc = int((conv == 1).sum())
        print(f"  converged {yc}, ran out of steps {nc}"
              + (f"  <- {100*nc/max(1, yc+nc):.1f}% of evaluated points are NOT relaxed "
                 f"values; raise P3_PES_RELAX_STEPS" if nc else ""))
    if ok.sum() == 0:
        print(f"  ALL POINTS NaN -- this surface is empty, not flat. Check the shard logs.")
    if ok.any():
        print(f"  E range {np.nanmin(E):.6f} .. {np.nanmax(E):.6f} eV "
              f"({(np.nanmax(E)-np.nanmin(E))*1e3:.0f} meV)")
        print(f"  total force calls {calls.sum()}")
        # How well the constraint was actually met: the placement is geometric and
        # approximate, so the achieved coordinates are recorded and checked, not assumed.
        dq = np.abs(ach[..., 0] - Q[None, :]); dr = np.abs(ach[..., 1] - R[:, None])
        print(f"  constraint error: |dq| max {np.nanmax(dq):.4f} A, "
              f"|dr| max {np.nanmax(dr):.4f} A")
    np.savez(d / f"pes2d_{mode}.npz", Q=Q, R=R, E=E, achieved=ach, calls=calls, converged=conv,
             q0=blobs[0]["q0"], r0=blobs[0]["r0"], donor=blobs[0]["donor"],
             acceptor=blobs[0]["acceptor"], shared_H=blobs[0]["shared_H"])
    print(f"  merged -> {d}/pes2d_{mode}.npz")

files = shards(d / "path_profile_shard*.json")
if files:
    print(f"\n=== path profile ({len(files)} shards) ===")
    prof = []
    for f in files:
        prof += json.load(open(f))["profile"]
    prof.sort(key=lambda p: p["frame"])
    gaps = [p["E_path"] - p["E_surface"] for p in prof
            if p["E_surface"] == p["E_surface"]]          # NaN-safe
    print(f"  {len(prof)} frames, {len(gaps)} with a surface value")
    if gaps:
        g = np.array(gaps) * 1e3
        print(f"  path-above-surface gap: max {g.max():.1f} meV, mean {g.mean():.1f} meV")
        print(f"  negative gaps (path BELOW surface -- would indicate a constraint or "
              f"convergence fault, not physics): {int((g < -0.81).sum())}")
        print(f"  THIS GAP IS THE DIAGNOSTIC: it answers whether (q_PT, r_OO) opened the "
              f"reaction up. It is not an error term.")
    json.dump(dict(source_shards=files, profile=prof),
              open(d / "path_profile.json", "w"), indent=1, default=str)
    print(f"  merged -> {d}/path_profile.json")

# ---------------------------------------------------------------- three models
d = ROOT / "p3_three_model"
files = shards(d / "three_model_shard*.json")
if files:
    print(f"\n=== three models ({len(files)} shards) ===")
    models, meta = {}, None
    for f in files:
        blob = json.load(open(f)); meta = meta or blob
        models.update(blob["models"])
    print("ABSOLUTE energies are NOT comparable across models (different training sets and")
    print("references). Only the within-model differences below are.\n")
    hdr = (f"  {'model':22s} {'dE(B-A) meV':>12s} {'barrier meV':>22s} "
           f"{'imaginary cm^-1':>22s} {'r_OO A':>9s} {'gate A':>8s}")
    print(hdr); print("  " + "-" * (len(hdr) - 2))
    for name, r in models.items():
        if r.get("status") != "done":
            print(f"  {name:22s} {r.get('status')}")
            continue
        b = r.get("barrier_meV") or []
        im = r.get("imaginary_cm1") or []
        g = r.get("gating_contraction_A") or []
        bs = f"{min(b):.1f}..{max(b):.1f}" if b else "none"
        ims = f"{min(im):.1f}..{max(im):.1f}" if im else "none"
        gs = f"{min(g):.3f}" if g else "-"
        print(f"  {name:22s} {r['dE_B_minus_A_meV']:12.3f} {bs:>22s} {ims:>22s} "
              f"{r['A']['r_OO']:9.4f} {gs:>8s}")
    # A graph invariant should not depend on the model. If it does, the models disagree
    # about connectivity, and that is a finding about the models, not a bug here.
    keys = {n: (r.get("A", {}).get("key"), r.get("B", {}).get("key"))
            for n, r in models.items() if r.get("status") == "done"}
    print()
    for n, (ka, kb) in keys.items():
        print(f"  {n:22s} key(A)={ka}  key(B)={kb}  differ={ka != kb}")
    uniq = {v for v in keys.values()}
    print(f"  all models agree on the (A, B) key pair: {len(uniq) <= 1}")
    if len(uniq) > 1:
        print(f"  THEY DO NOT -- the models disagree about connectivity, which is the")
        print(f"  finding, not an error. Pairs seen: {sorted(uniq)}")
    # POLAR-1 only: where the charge sits with the proton in flight
    for n, r in models.items():
        pol = (r.get("A") or {}).get("polar") or {}
        if not pol: continue
        print(f"\n  {n} atomic charges (POLAR-1 is the only model here that emits them):")
        print(f"    at A   Q(shared H) = {pol.get('Q_shared_H')}  dipole = {pol.get('dipole_D')}")
        for w in (r.get("walks") or []):
            wp = w.get("polar") or {}
            if wp:
                print(f"    at TS  Q(shared H) = {wp.get('Q_shared_H')}  "
                      f"dipole = {wp.get('dipole_D')}  (amp {w['amplitude']})")
        print(f"    Recorded, not interpreted: HANDOFF.md 5 leaves POLAR-1's long-range")
        print(f"    term unexplained (~1/r^2.6-3.3, not 1/r^2).")
    json.dump(dict(source_shards=files, models=models),
              open(d / "three_model.json", "w"), indent=1, default=str)
    print(f"\n  merged -> {d}/three_model.json")

if incomplete:
    print("\n*** MERGE INCOMPLETE -- these are partial results, not the experiment ***")
    for line in incomplete:
        print("  " + line)
    sys.exit(1)
print("\nall shards present")
