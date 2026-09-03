"""Two-panel figure: our points and our path, on the (q_PT, r_OO) surface. CPU only.

WHY TWO PANELS, and not one contour plot with the path coloured by energy.

The projection is honest about POSITION and dishonest about ENERGY if the two are mixed.
A PRRS path lives in 3N dimensions. The relaxed surface value at a projected point is the
minimum over the other 3N-2 coordinates, so the path generally lies ABOVE the surface at
its own projection. Colouring the path by its own energy on top of a surface contoured by
the surface's energy puts two different quantities in one colour scale and reads as if the
path were following the surface.

So: left panel is positions only -- contours of the surface, our stationary points, our
path. Right panel is energy -- the path's own energy against arc length, WITH the surface
energy along the same projected line overlaid.

THE GAP BETWEEN THOSE TWO CURVES IS THE RESULT, not an error bar. It answers the question
the coordinate pair was chosen to answer: did (q_PT, r_OO) open this reaction up? A small
gap means the two coordinates carry the reaction. A large gap means the remaining 3N-2
coordinates are doing work the picture does not show, and the projection is a cartoon.

Both a rigid and a relaxed surface are drawn when both exist, because their difference is
the same statement measured on the surface itself: the rigid section lies above the relaxed
one everywhere, so a barrier read off a rigid scan is an upper bound. Measured on the
methyl rotor alone (p3_rotor_resolution) that gap was 16.8 vs 1.255 meV, a factor of 13.

Reads only what is on disk and skips what is missing, so it can be run at any stage.
"""
import json, os, sys, glob, math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

ROOT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parents[2] / "runs"
FIGDIR = Path(__file__).resolve().parents[2] / "docs" / "figures"
FIGDIR.mkdir(parents=True, exist_ok=True)
MEV = 1e3


def load_surface(mode):
    f = ROOT / "p3_pes2d" / f"pes2d_{mode}.npz"
    if not f.exists():
        merged = sorted(glob.glob(str(ROOT / "p3_pes2d" / f"pes2d_{mode}_shard*.npz")))
        if not merged:
            return None
        print(f"  {mode}: only shards present, run `run_all.sh merge` first")
        return None
    b = np.load(f)
    E = b["E"]
    # An all-NaN grid is an EMPTY surface, not a flat one. The first version reported
    # "surfaces: ['relaxed','rigid']" for a relaxed grid that was 0/2304 finite and then
    # drew an all-NaN contour -- a plot that looks like a result and contains nothing.
    finite = int(np.isfinite(E).sum())
    if finite == 0:
        print(f"  {mode}: present but 0/{E.size} finite -- treating as ABSENT")
        return None
    conv = b["converged"] if "converged" in b else None
    if conv is not None and (conv == 0).any():
        n = int((conv == 0).sum())
        print(f"  {mode}: {finite}/{E.size} finite, but {n} points ran out of optimizer "
              f"steps and are masked out (not relaxed values)")
        E = np.where(conv == 0, np.nan, E)
        finite = int(np.isfinite(E).sum())
        if finite == 0:
            print(f"  {mode}: nothing left after masking -- treating as ABSENT")
            return None
    print(f"  {mode}: {finite}/{E.size} usable points")
    return dict(Q=b["Q"], R=b["R"], E=E)


def load_paths():
    """Every recorded stage, with q_PT and r_OO read from the jsonl the recorder wrote.

    Coordinates come from the record, never re-derived here: the recorder was configured
    with tracked_coordinates so that the figure and the record cannot disagree about what
    q_PT means.
    """
    out = {}
    for mdir in sorted((ROOT / "p3_paths").glob("*")):
        if not mdir.is_dir():
            continue
        stages = {}
        for sdir in sorted(mdir.glob("*")):
            jl = list(sdir.glob("*.jsonl"))
            if not jl:
                continue
            rows = [json.loads(l) for l in open(jl[0]) if l.strip()]
            rows = [r for r in rows if r.get("tracked")]
            if not rows:
                continue
            stages[sdir.name] = dict(
                q=np.array([r["tracked"]["q_PT"] for r in rows]),
                r=np.array([r["tracked"]["r_OO"] for r in rows]),
                E=np.array([np.nan if r.get("potential_eV") is None else r["potential_eV"]
                            for r in rows]),
                xyz=sdir / f"{rows[0]['stage']}.extxyz")
        if stages:
            out[mdir.name] = stages
    return out


def load_stationary():
    """A, B and TS per model, with coordinates recomputed from the .extxyz files.

    The q_PT stored in three_model.json is NOT trusted, for a reason that already bit
    twice. The projection coordinate is only a progress variable if its atom indices are
    held fixed; an earlier p3_three_model.py re-derived "the donor" per structure, so B --
    whose proton has moved to the other oxygen -- reported q_PT ~ -0.6 like A, and B's
    marker landed on top of A's. That was patched by rewriting the json, and then
    `run_all.sh merge` rebuilt the json from the shards and silently restored the old
    values. Recomputing here from the geometries, on one fixed reference triple, makes the
    figure independent of which version of which script wrote the json.

    The reference triple comes from A's own record (A's donor IS the reference donor by
    construction); only B's was ever wrong.
    """
    pts = {}
    tm = ROOT / "p3_three_model" / "three_model.json"
    if not tm.exists():
        return pts
    from ase.io import read as _read
    blob = json.load(open(tm))["models"]
    for name, r in blob.items():
        if r.get("status") != "done":
            continue
        a = r.get("A", {})
        D, ACC, H = a.get("donor"), a.get("acceptor"), a.get("shared_H")
        if D is None:
            continue

        def coord(atoms):
            p = atoms.positions
            return (float(np.linalg.norm(p[D] - p[H]) - np.linalg.norm(p[ACC] - p[H])),
                    float(np.linalg.norm(p[D] - p[ACC])))

        d = {}
        for tag, key in (("A", "A"), ("B", "B")):
            f = ROOT / "p3_three_model" / f"{name}_{tag}.extxyz"
            if f.exists():
                q, rr = coord(_read(str(f)))
                d[key] = (q, rr, r[key]["E_eV"])
            else:                       # fall back to the record, flagged in the console
                d[key] = (r[key]["q_PT"], r[key]["r_OO"], r[key]["E_eV"])
                print(f"  {name} {tag}: no extxyz, using the recorded q_PT (may be stale)")
        ts = []
        for w in r.get("walks", []):
            if w.get("saddle_order") != 1:
                continue
            f = ROOT / "p3_three_model" / f"{name}_ts_amp{w['amplitude']:.2f}.extxyz"
            if f.exists():
                q, rr = coord(_read(str(f)))
                ts.append((q, rr, w["E_eV"]))
            else:
                ts.append((w["q_PT"], w["r_OO"], w["E_eV"]))
        if ts:
            d["TS"] = ts
        pts[name] = d
    return pts


def load_profile():
    f = ROOT / "p3_pes2d" / "path_profile.json"
    if not f.exists():
        return None
    return json.load(open(f))["profile"]


def surface_along(prof, surf):
    """Surface energy at each path frame's projection, by interpolating the GRID.

    Not read from path_profile.json's own E_surface field. That field was computed by
    re-running the constrained optimisation from the source structure A for every frame,
    and it came back as a single constant -- exactly E(A), for all 491 frames -- because
    the call passed A's geometry straight to the constraint instead of first placing it at
    (q, r). The optimiser simply stayed at A. Constant "surface" plus a descending path
    then produced a mean gap of -9.1 meV, i.e. a path apparently BELOW the relaxed surface,
    which cannot happen.

    Interpolating the grid is better than fixing that call, not just cheaper: it makes the
    right panel's blue curve and the left panel's contours the SAME object. A separately
    recomputed section could disagree with the contours and there would be no way to tell
    which was wrong.

    Bilinear on the (R, Q) grid.
    """
    # Bilinear by hand rather than scipy: this repo's declared dependencies are ase and
    # numpy, and a figure script is not a good place to add a third.
    Q, R, E = surf["Q"], surf["R"], surf["E"]
    q = np.array([p["q_PT"] for p in prof], dtype=float)
    r = np.array([p["r_OO"] for p in prof], dtype=float)
    iq = np.clip(np.searchsorted(Q, q) - 1, 0, len(Q) - 2)
    ir = np.clip(np.searchsorted(R, r) - 1, 0, len(R) - 2)
    tq = (q - Q[iq]) / (Q[iq + 1] - Q[iq])
    tr = (r - R[ir]) / (R[ir + 1] - R[ir])
    out = ((1 - tr) * ((1 - tq) * E[ir, iq] + tq * E[ir, iq + 1])
           + tr * ((1 - tq) * E[ir + 1, iq] + tq * E[ir + 1, iq + 1]))
    # Outside the grid there is no surface value; say so rather than extrapolate.
    out = np.where((q >= Q[0]) & (q <= Q[-1]) & (r >= R[0]) & (r <= R[-1]), out, np.nan)
    return out


surfaces = {m: load_surface(m) for m in ("relaxed", "rigid")}
paths = load_paths()
stat = load_stationary()
profile = load_profile()
have = [k for k, v in surfaces.items() if v]
print(f"surfaces: {have or 'none'}   paths: {list(paths) or 'none'}   "
      f"stationary: {list(stat) or 'none'}   profile: {'yes' if profile else 'no'}")
if not have and not paths and not stat:
    raise SystemExit("nothing on disk yet -- run `paths` and `pes` first")

fig = plt.figure(figsize=(13.5, 5.6))
axL = fig.add_subplot(1, 2, 1)
axR = fig.add_subplot(1, 2, 2)

# ---------------- left panel: positions only ----------------
surf = surfaces.get("relaxed") or surfaces.get("rigid")
which = "relaxed" if surfaces.get("relaxed") else ("rigid" if surfaces.get("rigid") else None)
if surf is not None:
    E = (surf["E"] - np.nanmin(surf["E"])) * MEV
    # Colour scale is capped, not linear over the full range. The surface spans ~1800 meV
    # but the reaction only uses the bottom ~320 meV (the barrier), and a linear 24-level
    # map over 1800 meV renders the entire reactive region as one flat dark patch -- both
    # wells and the saddle channel become invisible, which is the opposite of the point.
    # The cap is set from the measured barrier so it is tied to the chemistry rather than
    # picked to look good, and the colourbar is marked `extend="max"` so the truncation is
    # visible instead of implied.
    CAP = float(os.environ.get("P3_FIG_CAP_MEV", "450"))
    levels = np.linspace(0, CAP, 19)
    cs = axL.contourf(surf["Q"], surf["R"], np.clip(E, 0, CAP), levels=levels,
                      cmap="viridis", extend="max")
    axL.contour(surf["Q"], surf["R"], E, levels=np.arange(50, CAP, 50), colors="white",
                linewidths=0.45, alpha=0.55)
    # One heavier contour at the measured barrier, so the saddle channel is locatable.
    bar = [t[2] for st in stat.values() for t in st.get("TS", [])]
    if bar:
        b0 = (min(bar) - min(st["A"][2] for st in stat.values())) * MEV
        if 0 < b0 < CAP:
            axL.contour(surf["Q"], surf["R"], E, levels=[b0], colors="white",
                        linewidths=1.3, alpha=0.9)
    cb = fig.colorbar(cs, ax=axL, extend="max")
    cb.set_label(f"{which} surface, meV above its own minimum "
                 f"(capped at {CAP:.0f}; max is {np.nanmax(E):.0f})")

COL = {"MACE-OFF24_medium": "#e41a1c", "MACE-omol-0-xl-4M": "#377eb8",
       "MACE-POLAR-1-M": "#4daf4a"}
for name, st in stat.items():
    c = COL.get(name, "k")
    axL.plot(*st["A"][:2], "o", color=c, ms=9, mec="k", mew=0.8, zorder=5,
             label=f"{name}: A, B")
    axL.plot(*st["B"][:2], "s", color=c, ms=9, mec="k", mew=0.8, zorder=5)
    for q, r, _ in st.get("TS", []):
        axL.plot(q, r, "*", color=c, ms=17, mec="k", mew=0.8, zorder=6)
for name, stages in paths.items():
    c = COL.get(name, "k")
    for sname, d in stages.items():
        climbing = "climb" in sname
        axL.plot(d["q"], d["r"], "--" if climbing else "-", color=c,
                 lw=1.3 if climbing else 2.2, alpha=0.95, zorder=4,
                 path_effects=[pe.Stroke(linewidth=(2.3 if climbing else 3.4),
                                         foreground="white"), pe.Normal()])
        axL.plot(d["q"][-1], d["r"][-1], "x", color=c, ms=9, mew=2.2, zorder=6)
axL.set_xlabel(r"$q_{\rm PT}=r({\rm O_d\!-\!H})-r({\rm O_a\!-\!H})$   ($\AA$)")
axL.set_ylabel(r"$r({\rm O\cdots O})$   ($\AA$)")
axL.set_title("positions only — path, endpoints, saddles\n"
              r"$\bigstar$ TS   $\bullet$ A   $\blacksquare$ B   $\times$ where a walk stopped",
              fontsize=10)
if stat:
    # NOT upper right: B sits at (q_PT ~ +0.61, r_OO ~ 2.52), which is exactly where an
    # upper-right legend box lands, and it hid B's marker completely -- the figure looked
    # like B had not been computed. Lower left is the empty corner (the reaction runs from
    # upper-left to upper-right through a low-r_OO saddle), and the box is made transparent
    # enough to see a marker through it either way.
    axL.legend(fontsize=7, loc="upper center", framealpha=0.82, ncol=3,
               bbox_to_anchor=(0.5, -0.13))
# Label the stationary points directly, so the figure does not depend on the legend being
# readable or in a safe place.
# Labels are offset per model, because all three sit within ~0.04 A of each other and a
# common offset stacked them into unreadable overlaps ("TSTS").
for k, (name, st) in enumerate(sorted(stat.items())):
    c = COL.get(name, "k")
    dy = (10, -2, -14)[k % 3]
    for lab, key in (("A", "A"), ("B", "B")):
        q, r = st[key][:2]
        axL.annotate(lab, (q, r), textcoords="offset points", xytext=(11, dy),
                     fontsize=9, fontweight="bold", color=c, zorder=7,
                     path_effects=[pe.withStroke(linewidth=2.2, foreground="white")])
    for q, r, _ in st.get("TS", [])[:1]:
        axL.annotate("TS", (q, r), textcoords="offset points", xytext=(11, dy - 4),
                     fontsize=9, fontweight="bold", color=c, zorder=7,
                     path_effects=[pe.withStroke(linewidth=2.2, foreground="white")])
# Margin so a marker never sits on the frame.
if stat:
    qs = [v[:2][0] for st in stat.values() for k, v in st.items() if k in ("A", "B")]
    qs += [t[0] for st in stat.values() for t in st.get("TS", [])]
    rs = [v[:2][1] for st in stat.values() for k, v in st.items() if k in ("A", "B")]
    rs += [t[1] for st in stat.values() for t in st.get("TS", [])]
    if surf is None:
        axL.set_xlim(min(qs) - 0.12, max(qs) + 0.12)
        axL.set_ylim(min(rs) - 0.03, max(rs) + 0.04)

# ---------------- right panel: energy, with the gap ----------------
drew_gap = False
if profile and surf is not None:
    s = np.array([p["arc_length"] for p in profile])
    ep = np.array([p["E_path"] for p in profile])
    es = surface_along(profile, surf)
    ok = np.isfinite(es) & np.isfinite(ep)
    ref = min(np.nanmin(ep), np.nanmin(es[ok])) if ok.any() else np.nanmin(ep)
    axR.plot(s, (ep - ref) * MEV, "-", color="#e41a1c", lw=2, label="path's own energy")
    if ok.any():
        axR.plot(s[ok], (es[ok] - ref) * MEV, "-", color="#377eb8", lw=2,
                 label=f"{which} surface, interpolated from the grid")
        axR.fill_between(s[ok], (es[ok] - ref) * MEV, (ep[ok] - ref) * MEV,
                         color="0.6", alpha=0.35, label="gap")
        g = (ep[ok] - es[ok]) * MEV
        # The sign is read off the data, not asserted. A path cannot lie below a properly
        # converged relaxed surface, so a negative gap is a statement about the surface or
        # the interpolation -- never about the chemistry -- and is labelled as such.
        neg = g < -0.81          # 0.81 meV is P2's measured endpoint spread
        head = (f"gap: max {g.max():+.1f}, mean {g.mean():+.1f} meV"
                f"   ({int(neg.sum())}/{len(g)} frames below the surface)")
        axR.axhline(0, color="0.5", lw=0.5)
        # The negative side was first labelled "interpolation residual". It is not.
        # Measured: every negative frame lies in the B well (q +0.599..+0.612,
        # r 2.509..2.518) and the gap there is a flat -2.08..-0.81 meV -- a systematic
        # offset, not scatter. The surface in that well sits 2.13 meV above B's true
        # minimum because it was relaxed to fmax 5e-3 while the path was quenched to 2e-3.
        # That is a convergence-tolerance mismatch in the surface, and it is now fixed at
        # source (p3_pes2d.PES_FMAX defaults to the path's tolerance).
        axR.set_title(head + "\nabove = the other $3N-2$ coordinates do work (real);  "
                             "below = surface not relaxed as tightly as the path",
                      fontsize=9)
        drew_gap = True
        print(f"  gap: max {g.max():+.1f}  mean {g.mean():+.1f} meV, "
              f"{int(neg.sum())}/{len(g)} frames below the surface")
    axR.set_xlabel(r"arc length in full $3N$ space ($\AA$)")
    axR.legend(fontsize=8, loc="upper right")
elif paths:
    # No surface profile yet: show the recorded path energies against q_PT. This is a
    # projection of POSITION against energy, not a claim about the surface.
    for name, stages in paths.items():
        c = COL.get(name, "k")
        for sname, d in stages.items():
            if not np.isfinite(d["E"]).any():
                continue
            ref = np.nanmin(d["E"])
            axR.plot(d["q"], (d["E"] - ref) * MEV, "-" if "descent" in sname else "--",
                     color=c, lw=1.6, label=f"{name} {sname}")
    axR.set_xlabel(r"$q_{\rm PT}$ ($\AA$)")
    axR.set_title("recorded path energy (no surface profile yet)\n"
                  "run `pes` with PATHFILE to get the overlay and the gap", fontsize=10)
    axR.legend(fontsize=7)
elif stat:
    # Neither surface profile nor path frames yet, but the stationary points ARE measured,
    # so the panel shows those rather than sitting empty. Three points per model, joined
    # only to make the model comparison readable.
    #
    # Labelled for what it is: A, TS and B are stationary points on each model's own
    # surface. The joining lines are NOT a reaction path and NOT a surface section -- no
    # intermediate geometry was computed. What the panel does show honestly is the spread
    # BETWEEN models, which is the measured result: barriers differ by 3.7x.
    for name, st in sorted(stat.items()):
        c = COL.get(name, "k")
        E_A = st["A"][2]
        xs = [st["A"][0]]; ys = [0.0]
        ts = sorted(st.get("TS", []), key=lambda t: abs(t[0]))
        if ts:
            xs.append(ts[0][0]); ys.append((ts[0][2] - E_A) * MEV)
        xs.append(st["B"][0]); ys.append((st["B"][2] - E_A) * MEV)
        axR.plot(xs, ys, ":", color=c, lw=1.2, alpha=0.7)
        axR.plot(xs[0], ys[0], "o", color=c, ms=9, mec="k", mew=0.8)
        if ts:
            axR.plot(xs[1], ys[1], "*", color=c, ms=17, mec="k", mew=0.8)
            axR.annotate(f"{ys[1]:+.0f}", (xs[1], ys[1]), textcoords="offset points",
                         xytext=(11, -3), fontsize=9, fontweight="bold", color=c)
        axR.plot(xs[-1], ys[-1], "s", color=c, ms=9, mec="k", mew=0.8,
                 label=f"{name}")
        axR.annotate(f"{ys[-1]:+.1f}", (xs[-1], ys[-1]), textcoords="offset points",
                     xytext=(9, 4), fontsize=8, color=c)
    axR.axhline(0, color="0.5", lw=0.6, ls="-")
    axR.set_xlabel(r"$q_{\rm PT}$ ($\AA$)")
    axR.set_title("stationary points only, each relative to its OWN A (meV)\n"
                  "dotted lines join three points — NOT a path, NOT a surface section",
                  fontsize=10)
    axR.legend(fontsize=7, loc="center left")
axR.set_ylabel("energy above the panel minimum (meV)")

if not drew_gap:
    axR.text(0.5, 0.015,
             "GAP NOT YET MEASURED — run `paths` then `pes` to get the diagnostic",
             transform=axR.transAxes, ha="center", fontsize=8, color="#b00020")

# Crop r_OO to where the reaction actually happens. The grid runs to 3.00 A, but both
# wells sit at 2.52 and the channel at 2.32, so everything above ~2.65 is the O...O simply
# pulled apart -- high, featureless, and half the panel. Cropping is a view choice and does
# not touch the data; the full grid stays in the npz.
if surf is not None:
    axL.set_ylim(float(surf["R"].min()), float(os.environ.get("P3_FIG_RMAX", 2.68)))
for ax in (axL, axR):
    ax.grid(alpha=0.25, lw=0.4)
fig.suptitle("P3 benzoylacetone — PRRS points and path projected onto "
             r"$(q_{\rm PT},\, r_{\rm O\cdots O})$", fontsize=12)
fig.tight_layout(rect=(0, 0, 1, 0.94))
fig.savefig(FIGDIR / "p3_projection.png", dpi=170, bbox_inches="tight")
print(f"wrote {FIGDIR}/p3_projection.png")

# ---------------- rigid vs relaxed, when both exist ----------------
if surfaces.get("rigid") and surfaces.get("relaxed"):
    d = (surfaces["rigid"]["E"] - surfaces["relaxed"]["E"]) * MEV
    f2, a2 = plt.subplots(figsize=(6.4, 5.2))
    cs = a2.contourf(surfaces["relaxed"]["Q"], surfaces["relaxed"]["R"], d,
                     levels=24, cmap="magma")
    f2.colorbar(cs, ax=a2, label="rigid minus relaxed (meV)")
    a2.set_xlabel(r"$q_{\rm PT}$ ($\AA$)"); a2.set_ylabel(r"$r({\rm O\cdots O})$ ($\AA$)")
    a2.set_title("how much the other $3N-2$ coordinates absorb\n"
                 "a barrier read off the rigid scan is an upper bound", fontsize=10)
    f2.tight_layout()
    f2.savefig(FIGDIR / "p3_rigid_minus_relaxed.png", dpi=170, bbox_inches="tight")
    print(f"wrote {FIGDIR}/p3_rigid_minus_relaxed.png   "
          f"range {np.nanmin(d):.1f} .. {np.nanmax(d):.1f} meV")
