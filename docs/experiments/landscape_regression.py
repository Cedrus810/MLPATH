"""Per-conformer regression for the P3 landscape. Exit code is the conclusion.

The landscape does NOT run `confirm_minimum`. This file claimed it did, and so does
p3_landscape_full.py's own docstring; the code at p3_landscape_full.py:71 calls
`curvature_spectrum` and counts the negatives itself, which is the spectrum without any of
the gates built on it. So this regression cannot see a change to `confirm_minimum` at all --
the job that reaches it is `ts`, through `descend_saddle`.

What this file does check is real and worth keeping: a per-conformer baseline of 300 rows,
each with n_minus, lam1, E_eV, r_OO and chemical_key, which pins the spectrum and the
chemistry against arithmetic drift.

What a gate failure looks like here: `n_minus` goes from an integer to None, or from 1 to 0
(a rotor saddle silently reclassified as a minimum), or 229 confirmed minima becomes some
other number. Those are counted separately from numeric drift in lam1, because they are
different failures: the first is a gate verdict changing, the second is arithmetic moving.

Usage:
    python docs/experiments/landscape_regression.py --write   # record the baseline
    python docs/experiments/landscape_regression.py           # check against it
    python docs/experiments/landscape_regression.py --lam1-tol 1e-6
Nonzero exit on any verdict change, or on lam1 drift beyond the tolerance.
"""
import json, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
LIVE = ROOT / "runs" / "p3_landscape_full" / "landscape.json"
BASELINE = ROOT / "docs" / "experiments" / "landscape_baseline.json"
WRITE = "--write" in sys.argv
TOL = 1e-6
if "--lam1-tol" in sys.argv:
    TOL = float(sys.argv[sys.argv.index("--lam1-tol") + 1])

# Verdict fields: a change here is a gate deciding differently, not arithmetic drifting.
VERDICT = ("n_minus", "key", "graph_hash", "fragments", "stereo_unresolved",
           "n_locked_parity", "n_OH", "n_edges")
# Numeric fields: compared with a tolerance.
NUMERIC = {"lam1": None, "E_eV": None, "r_OO": None, "fmax": None, "floor_residual": None}

if not LIVE.exists():
    sys.exit(f"no landscape at {LIVE} -- run `run_all.sh landscape` first")
rows = json.load(open(LIVE))["rows"]
current = {int(r["conf"]): r for r in rows}

if WRITE:
    keep = {c: {k: r.get(k) for k in list(VERDICT) + list(NUMERIC)}
            for c, r in current.items()}
    json.dump({"n": len(keep), "rows": keep}, open(BASELINE, "w"), indent=1,
              sort_keys=True, default=str)
    import collections
    print(f"baseline written: {len(keep)} conformers, "
          f"n_minus distribution {dict(collections.Counter(r['n_minus'] for r in current.values()))}")
    print(f"  -> {BASELINE}")
    sys.exit(0)

if not BASELINE.exists():
    sys.exit(f"no baseline at {BASELINE}; run with --write against the OLD run first")
base = json.load(open(BASELINE))["rows"]
base = {int(k): v for k, v in base.items()}

missing = sorted(set(base) - set(current))
extra = sorted(set(current) - set(base))
verdict_changes, numeric = [], {k: 0.0 for k in NUMERIC}
numeric_worst = {k: None for k in NUMERIC}
for c, want in base.items():
    got = current.get(c)
    if got is None:
        continue
    for f in VERDICT:
        if want.get(f) != got.get(f):
            verdict_changes.append(f"conf{c:04d}: {f} {want.get(f)!r} -> {got.get(f)!r}")
    for f in NUMERIC:
        a, b = want.get(f), got.get(f)
        if a is None or b is None:
            if a != b:
                verdict_changes.append(f"conf{c:04d}: {f} {a!r} -> {b!r} (one is None)")
            continue
        d = abs(float(b) - float(a))
        if d > numeric[f]:
            numeric[f] = d
            numeric_worst[f] = c

n = len(base)
kept = sum(1 for c in base if c in current)
print(f"{kept}/{n} conformers present" + (f"  MISSING {missing[:8]}" if missing else ""))
if extra:
    print(f"  {len(extra)} conformers not in the baseline: {extra[:8]}")
print(f"verdict fields unchanged: {n - len({l.split(':')[0] for l in verdict_changes})}/{n}"
      if verdict_changes else f"verdict fields unchanged: {kept}/{kept}")
for f in NUMERIC:
    w = numeric_worst[f]
    print(f"  max |Δ{f}| = {numeric[f]:.3e}" + (f"  (conf{w:04d})" if w is not None else ""))

# The recheck gate's own output, when the run carries it.
have_recheck = [r for r in current.values() if "internal_mode_recheck" in r]
if have_recheck:
    diag = [r["internal_mode_recheck"] for r in have_recheck]
    agree = sum(1 for d in diag if d.get("agree"))
    stable = sum(1 for d in diag if d.get("resolution_stable"))
    print(f"internal_mode_recheck present on {len(have_recheck)}/{kept}: "
          f"agree {agree}, resolution_stable {stable}")
    bad = [r for r, d in zip(have_recheck, diag) if not d.get("agree")]
    for r in bad[:6]:
        print(f"  DISAGREE conf{int(r['conf']):04d}: "
              f"{r['internal_mode_recheck'].get('resolution_diagnosis')}")
else:
    print("internal_mode_recheck: absent (this run predates the gate)")

fail = []
if missing:
    fail.append(f"{len(missing)} conformers missing")
if verdict_changes:
    fail.append(f"{len(verdict_changes)} verdict change(s)")
over = [f for f in NUMERIC if numeric[f] > TOL]
if over:
    fail.append(f"numeric drift beyond {TOL:g} in {over}")
if not fail:
    print(f"\nlandscape unchanged (lam1 tolerance {TOL:g})")
    sys.exit(0)
print(f"\nLANDSCAPE REGRESSION: {'; '.join(fail)}")
for line in verdict_changes[:30]:
    print("  " + line)
if len(verdict_changes) > 30:
    print(f"  ... and {len(verdict_changes)-30} more")
sys.exit(1)
