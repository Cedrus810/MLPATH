"""Golden-sample regression for chemical_key / reaction_event_key digests.

Purpose: let a change to src/prrs/ prove it did not move any persisted digest. The keys in
runs/ are load-bearing -- P1 and P2 are frozen on them and P3's probes are recorded against
them -- so a refactor that shifts a digest by one byte invalidates finished experiments
silently.

The P3 corpus is the useful part of this set and is deliberately awkward:
  22 atoms with an aromatic ring, so admissible_bond_orders really enumerates a Kekule
  pair (2 solutions) and locked_edges has to reject the ring bonds while keeping the enol
  C=C. Anything that perturbs index handling shows up here first.
  300 relaxed conformers of the same substance, so a change that is conformation-dependent
  (it must not be) shows up as a split in what is currently a single key.
  Three independently trained models, which agree on the (A, B) key pair -- a graph
  invariant should not depend on the potential, so a disagreement is a real defect.

A mismatch is not automatically a defect. Some changes to chemistry.py are SUPPOSED to
move digests -- the stereochemistry parity tolerance is the known example: a near-planar
centre should not split into a new macrostate because of floating-point jitter, so fixing
it moves the parity of exactly those frames and nothing else. "All digests unchanged" is
therefore the wrong acceptance criterion for that kind of change.

So a mismatch is classified rather than just counted:

    graph_invariant  graph_hash or fragments moved. The molecular graph itself changed.
                     This is severe under any intent: those are what "same substance"
                     means, and P1/P2 are frozen on them.
    parity_only      only tetrahedral/locked parity counts or stereo_unresolved moved,
                     with graph_hash and fragments intact. Consistent with a deliberate
                     stereochemistry change -- still requires chemical judgement, but it
                     is the signature a parity fix should have.
    key_only         the key moved while every recorded component is identical. The key is
                     a hash OF those components, so this means the component set or its
                     serialisation changed, not the chemistry.

Two baselines, on purpose (2026-09-04):

    digest_baseline.json          v1. FROZEN. P1/P2 are recorded against it and it is not
                                  rewritten. Comparison projects onto the fields it holds,
                                  so adding a field cannot turn it red.
    digest_baseline_stereo.json   v2. Full parity CONTENT, not counts. v1 records only
                                  len(locked_bond_parity) and len(tetrahedral_parity), so a
                                  parity flipping SIGN leaves the counts intact and is
                                  invisible to it. An identity migration cannot be reviewed
                                  with an instrument that cannot see the thing being migrated.

`--allow parity_only` therefore covers count changes ONLY. A sign flip is classified
`parity_content` and has to be allowed explicitly, because it means a configuration
inverted -- a different substance, not a tolerance nudge.

Usage:
    python docs/experiments/digest_regression.py --write         # record v1 (do not)
    python docs/experiments/digest_regression.py --write-stereo  # record v2
    python docs/experiments/digest_regression.py                # check against it
    python docs/experiments/digest_regression.py --allow parity_only
                                                               # accept a parity-only diff
Exit code is nonzero unless every mismatch falls in an --allow'ed class.
"""

import collections
import glob
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from ase.io import read
from prrs.chemistry import chemical_key, reaction_event_key, classify_transition
from prrs.config import SearchConfig

ROOT = Path(__file__).resolve().parents[2]
BASELINE = ROOT / "docs" / "experiments" / "digest_baseline.json"
# v1 is frozen: P1/P2 are recorded against it and it is never rewritten. Comparison
# against it projects onto the fields it actually holds, so adding a field cannot turn
# it red. v2 carries the full stereo content and is what an identity migration is
# reviewed against.
BASELINE_V2 = ROOT / "docs" / "experiments" / "digest_baseline_stereo.json"
cfg = SearchConfig()
WRITE = "--write" in sys.argv
ALLOW = set()
if "--allow" in sys.argv:
    ALLOW = set(sys.argv[sys.argv.index("--allow") + 1].split(","))
GRAPH_FIELDS = ("graph_hash", "fragments")
PARITY_FIELDS = ("n_locked_parity", "n_tetrahedral", "stereo_unresolved")
# Content, not counts. A sign flip lands here and nowhere else -- it is a configuration
# inverting, which is a different substance, and it must never be swept in with a count
# change under `--allow parity_only`.
PARITY_CONTENT_FIELDS = ("locked_bond_parity", "tetrahedral_parity")
# Fields whose recorded values came from a defective computation. Kept in the baseline
# file as the record of the defect, skipped on comparison so a wrong number cannot go on
# being the standard. See the note at the `transition_class` assignment.
RETIRED_EVENT_FIELDS = {"classify"}


def structures():
    """Every persisted geometry whose key a frozen or recorded result depends on."""
    seen = set()
    for pat in (
        "runs/p3_*.extxyz",
        "runs/p3_three_model/*.extxyz",
        "runs/p3_landscape_full/conf*.extxyz",
        "runs/p3_paths/*/*/*.extxyz",
        "runs/p1_round8/chemical/*.extxyz",
        "runs/p2_round7/chemical/*.extxyz",
    ):
        for f in sorted(glob.glob(str(ROOT / pat))):
            rel = str(Path(f).relative_to(ROOT))
            if rel in seen:
                continue
            seen.add(rel)
            yield rel, f


def digest_of(path):
    frames = read(path, index=":")
    out = []
    for i, a in enumerate(frames):
        k = chemical_key(a, cfg.bond_scale, cfg.active_atoms)
        out.append(
            {
                "frame": i,
                "key": k["key"],
                "graph_hash": k["graph_hash"],
                "fragments": k["fragments"],
                "stereo_unresolved": k["stereo_unresolved"],
                # Counts are the v1 fields. None when the annotator could not decide
                # which bonds carry configuration: v1 recorded 0 there, which is the
                # same value a molecule with genuinely no locked bonds gives, and that
                # conflation is what the v2 identity change removes. Reporting 0 here
                # would carry the old lie forward past the fix.
                "n_locked_parity": (
                    None if k["locked_bond_parity"] is None else len(k["locked_bond_parity"])
                ),
                "n_tetrahedral": len(k["tetrahedral_parity"]),
                # v2, added 2026-09-04. Counts are blind to a parity flipping SIGN:
                # the number of recorded stereocentres does not change, so a
                # configuration inverting is invisible to v1. An identity migration
                # cannot be reviewed with an instrument that cannot see it.
                "locked_bond_parity": k["locked_bond_parity"],
                "tetrahedral_parity": k["tetrahedral_parity"],
            }
        )
    return out


# Event keys: the pairs whose reaction_event_key is quoted in the frozen/recorded results.
PAIRS = [
    ("runs/p3_A_source.extxyz", "runs/p3_B_source.extxyz", "P3 A->B, event 6e1ee45c"),
    (
        "runs/p3_three_model/MACE-OFF24_medium_A.extxyz",
        "runs/p3_three_model/MACE-OFF24_medium_B.extxyz",
        "P3 OFF24",
    ),
    (
        "runs/p3_three_model/MACE-omol-0-xl-4M_A.extxyz",
        "runs/p3_three_model/MACE-omol-0-xl-4M_B.extxyz",
        "P3 omol-0",
    ),
    (
        "runs/p3_three_model/MACE-POLAR-1-M_A.extxyz",
        "runs/p3_three_model/MACE-POLAR-1-M_B.extxyz",
        "P3 POLAR-1",
    ),
]

# --------------------------------------------------------------------------------
# Identity relations. These are assertions about what the identity layer MUST hold,
# stated independently of any baseline -- which is the point: a baseline can be
# rewritten, and then it no longer constrains anything. These survive that.
#
# Added 2026-09-04 alongside the v2 stereo baseline, for the same reason: an identity
# migration has to be reviewable, and "the digests all moved" is not a review.
RELATIONS = [
    # (label, kind, argument, expected)
    (
        "300+ relaxed conformers of one substance collapse to one key",
        "one_key",
        "runs/p3_landscape_full/conf*.extxyz",
        1,
    ),
    (
        "P1's two chemical structures are one substance (degenerate transfer)",
        "one_key",
        "runs/p1_round8/chemical/*.extxyz",
        1,
    ),
    (
        "P2's two chemical nodes really are two substances",
        "n_keys",
        "runs/p2_round7/chemical/*.extxyz",
        2,
    ),
    (
        "P3 A and B are different substances",
        "differ",
        ("runs/p3_A_source.extxyz", "runs/p3_B_source.extxyz"),
        True,
    ),
]


# Reviewed identity migrations. NOT `--allow`: each is an enumerated, permanent record of
# a difference that was examined frame by frame and accepted, so the frozen baselines keep
# guarding everything else. A frame whose diff matches a migration exactly is counted as
# `migrated`; any OTHER movement on the same frames still fails.
MIGRATIONS = [
    {
        "id": "2026-09-04-unknown-is-not-empty",
        "why": (
            "Two changes to how an unresolved configuration is represented, landing on the "
            "same frames.\n"
            "(1) `locked_bond_parity` was computed with an empty locked set when the "
            "annotator could not decide which bonds carry configuration, producing [] -- "
            "the same value a molecule with genuinely no locked bonds gives. Unknown is "
            "now None.\n"
            "(2) The hashed `stereo_unresolved` carried the REASON. Hashing a diagnostic "
            "makes it part of the substance's name: improve the resolver's failure "
            "reporting and every affected structure silently becomes a different species. "
            "The hash now carries only whether it resolved; the reason is recorded outside "
            "it. Resolved frames hash None exactly as before, so all 2009 of them keep "
            "their keys and every key quoted in a frozen acceptance document still "
            "recomputes -- checked: P3 A 8710064d12ce6a5a, P3 B b1a4338a45f1d300."
        ),
        "old_key": "65673127b9550ac7023a643a71c78956c9a8b2d964b7058d187cf224035be4f7",
        "new_key": "4b2be76140f37b5cd9c3a1a05f0b4c7657d6db958a6a3747959811bd38aa7dd4",
        # Every affected frame had exactly this state and exactly these field moves.
        "requires": {"stereo_unresolved": "no_assignment_satisfies_valences"},
        "moves": {"locked_bond_parity": ([], None), "n_locked_parity": (0, None)},
        "frames": 127,
        "review": (
            "127 frames, 120 in p3_paths and 7 in p3_three_model. ONE (old key -> new key) "
            "pair, a bijection: nothing merged, nothing split, and the new key collides "
            "with none of the 6 unaffected keys. Graph hash, fragments, tetrahedral parity "
            "and every event field are unchanged, and all four identity relations still "
            "hold. See docs/DECISIONS_PENDING.md."
        ),
    }
]


def migrated(want, got):
    """Does this frame's difference match a reviewed migration exactly?

    Exactly: same old and new key, same prior state, and the moved fields are precisely
    the ones the migration declares, with precisely the declared values. A frame that
    moved further than the migration says is not migrated -- it is a mismatch.
    """
    for m in MIGRATIONS:
        if want.get("key") != m["old_key"] or got.get("key") != m["new_key"]:
            continue
        if any(want.get(f) != v for f, v in m["requires"].items()):
            continue
        moved = {f for f in want if f != "frame" and want[f] != got.get(f)}
        declared = set(m["moves"]) & set(want)
        if moved != declared | {"key"} - {"key"} | ({"key"} if "key" in want else set()):
            continue
        if all(
            want.get(f) == a and got.get(f) == b
            for f, (a, b) in m["moves"].items()
            if f in want
        ):
            return m["id"]
    return None


def relation_keys(pattern):
    found = collections.Counter()
    for f in sorted(glob.glob(str(ROOT / pattern))):
        for a in read(f, index=":"):
            found[chemical_key(a, cfg.bond_scale, cfg.active_atoms)["key"]] += 1
    return found


def check_relations():
    """Returns a list of (ok, label, detail). Never consults a baseline."""
    out = []
    for label, kind, arg, expect in RELATIONS:
        if kind in ("one_key", "n_keys"):
            found = relation_keys(arg)
            if not found:
                out.append((None, label, "corpus absent, skipped"))
                continue
            ok = len(found) == expect
            out.append(
                (
                    ok,
                    label,
                    f"{len(found)} distinct key(s), expected {expect}: "
                    + ", ".join(f"{k[:16]}x{n}" for k, n in found.most_common(4)),
                )
            )
        elif kind == "differ":
            a, b = (ROOT / arg[0]), (ROOT / arg[1])
            if not (a.exists() and b.exists()):
                out.append((None, label, "structures absent, skipped"))
                continue
            ka = chemical_key(read(str(a)), cfg.bond_scale, cfg.active_atoms)["key"]
            kb = chemical_key(read(str(b)), cfg.bond_scale, cfg.active_atoms)["key"]
            out.append(((ka != kb) == expect, label, f"{ka[:16]} vs {kb[:16]}"))
    return out


current = {"structures": {}, "events": {}}
n = 0
for rel, f in structures():
    try:
        current["structures"][rel] = digest_of(f)
        n += 1
    except Exception as exc:
        current["structures"][rel] = {"error": f"{type(exc).__name__}: {exc}"}
for a, b, label in PAIRS:
    pa, pb = ROOT / a, ROOT / b
    if not (pa.exists() and pb.exists()):
        continue
    A, B = read(str(pa)), read(str(pb))
    ef = reaction_event_key(A, B, cfg.bond_scale, cfg.active_atoms)
    er = reaction_event_key(B, A, cfg.bond_scale, cfg.active_atoms)
    current["events"][label] = {
        "forward": ef["key"],
        "reverse": er["key"],
        "direction_quotiented": ef["key"] == er["key"],
        "broken": ef["broken"],
        "formed": ef["formed"],
        # `classify` is RETIRED, not recomputed. It was called as
        #     classify_transition(A, B, cfg.bond_scale, cfg.active_atoms)
        # against a signature of (source_key, endpoint_key, broken, formed, ...), so the
        # bond_scale 1.2 landed in the `broken` slot as a truthy value. Verified: that
        # call returns "degenerate_reaction" for a pair where NOTHING changed, which the
        # correct call classifies as "conformational_transition". The recorded values in
        # the v1 baseline are therefore not a classification regression and must not be
        # treated as one. They stay in the file as the record of the defect; the field is
        # skipped on comparison (RETIRED_EVENT_FIELDS) rather than deleted.
        "transition_class": str(
            classify_transition(
                chemical_key(A, cfg.bond_scale, cfg.active_atoms)["key"],
                chemical_key(B, cfg.bond_scale, cfg.active_atoms)["key"],
                ef["broken"],
                ef["formed"],
            )
        ),
    }

if "--write-stereo" in sys.argv:
    # v2 only. v1 is never rewritten from here -- that is the whole point of keeping two.
    json.dump(current, open(BASELINE_V2, "w"), indent=1, sort_keys=True, default=str)
    print(
        f"stereo baseline written: {n} files, "
        f"{sum(len(v) for v in current['structures'].values() if isinstance(v, list))} frames"
    )
    print(f"  -> {BASELINE_V2}")
    sys.exit(0)

if WRITE:
    if BASELINE.exists():
        sys.exit(
            f"refusing to overwrite the frozen v1 baseline at {BASELINE}.\n"
            "P1 and P2 are recorded against it. If a new definition is intended, write\n"
            "a new baseline with --write-stereo and keep this one."
        )
    json.dump(current, open(BASELINE, "w"), indent=1, sort_keys=True, default=str)
    keys = {
        d["key"][:16] for v in current["structures"].values() if isinstance(v, list) for d in v
    }
    print(
        f"baseline written: {n} files, {sum(len(v) for v in current['structures'].values() if isinstance(v, list))} frames, "
        f"{len(keys)} distinct chemical_key, {len(current['events'])} event pairs"
    )
    print(f"  -> {BASELINE}")
    sys.exit(0)

if not BASELINE.exists():
    sys.exit(f"no baseline at {BASELINE}; run with --write first")
base = json.load(open(BASELINE))
bad, classified = [], collections.Counter()


def classify(want, got):
    """Which fields moved decides how serious a frame-level mismatch is."""
    if any(want.get(f) != got.get(f) for f in GRAPH_FIELDS):
        return "graph_invariant"
    if any(f in want and want[f] != got.get(f) for f in PARITY_CONTENT_FIELDS):
        return "parity_content"
    if any(want.get(f) != got.get(f) for f in PARITY_FIELDS):
        return "parity_only"
    return "key_only"


def compare_frames(want, got, rel, bad, classified):
    """Compare one file's frames, projected onto the fields the baseline actually has.

    Projection is what lets v1 stay frozen while v2 records more. A baseline written
    before a field existed must not be turned red by that field appearing.
    """
    if not (isinstance(want, list) and isinstance(got, list) and len(want) == len(got)):
        bad.append(("shape", f"{rel}: shape changed"))
        classified["shape"] += 1
        return
    for w, g in zip(want, got):
        moved = [k for k in w if k != "frame" and w[k] != g.get(k)]
        if not moved:
            continue
        known = migrated(w, g)
        if known:
            classified[f"migrated:{known}"] += 1
            continue
        cls = classify(w, g)
        classified[cls] += 1
        bad.append(
            (
                cls,
                f"[{cls}] {rel} frame {w['frame']}: "
                + "; ".join(f"{k}: {w[k]!r} -> {g.get(k)!r}" for k in moved),
            )
        )


for rel, want in base["structures"].items():
    got = current["structures"].get(rel)
    if got is None:
        bad.append(("missing", f"{rel}: MISSING now"))
        classified["missing"] += 1
    else:
        compare_frames(want, got, rel, bad, classified)
if not BASELINE_V2.exists():
    sys.exit(
        f"no stereo baseline at {BASELINE_V2}.\n"
        "v1 records only the COUNT of parity entries, so a configuration inverting is\n"
        "invisible to it. Running without v2 is not a weaker check, it is a check with a\n"
        "hole in exactly the place an identity migration lives. Write it with\n"
        "--write-stereo before relying on this script."
    )
base2 = json.load(open(BASELINE_V2))
if True:
    for rel, want in base2["structures"].items():
        got = current["structures"].get(rel)
        if got is None:
            bad.append(("missing", f"{rel}: MISSING now (stereo baseline)"))
            classified["missing"] += 1
        else:
            compare_frames(want, got, rel, bad, classified)

for label, want in base["events"].items():
    got = current["events"].get(label) or {}
    moved = {
        k: (v, got.get(k))
        for k, v in want.items()
        if k not in RETIRED_EVENT_FIELDS and got.get(k) != v
    }
    if moved:
        bad.append(
            (
                "event",
                f"[event] {label}: "
                + "; ".join(f"{k}: {a!r} -> {b!r}" for k, (a, b) in moved.items()),
            )
        )
        classified["event"] += 1
if BASELINE_V2.exists():
    for label, want in base2.get("events", {}).items():
        got = current["events"].get(label) or {}
        moved = {k: (v, got.get(k)) for k, v in want.items() if got.get(k) != v}
        if moved:
            bad.append(
                (
                    "event",
                    f"[event v2] {label}: "
                    + "; ".join(f"{k}: {a!r} -> {b!r}" for k, (a, b) in moved.items()),
                )
            )
            classified["event"] += 1

relations = check_relations()
broken = [r for r in relations if r[0] is False]
for ok, label, detail in relations:
    mark = "SKIP" if ok is None else ("ok  " if ok else "FAIL")
    print(f"  [{mark}] {label}\n         {detail}")
print()
if broken:
    bad.extend(("relation", f"[relation] {label}: {detail}") for _, label, detail in broken)
    classified["relation"] += len(broken)

nf = sum(len(v) for v in current["structures"].values() if isinstance(v, list))
migrations = {c: k for c, k in classified.items() if c.startswith("migrated:")}
if not bad:
    print(f"digests unchanged: {n} files, {nf} frames, {len(current['events'])} event pairs")
    for cls, count in sorted(migrations.items()):
        # Never report a migrated frame as "unchanged". It changed; the change was
        # reviewed. Saying otherwise would hide the one class of movement that a reader
        # most needs to see when deciding whether a result still means what it meant.
        #
        # The count is per COMPARISON, so a frame present in both baselines is counted
        # twice. Said plainly rather than halved, because which baseline saw it is part
        # of what a reader is checking.
        print(f"  plus {count} comparison(s) on reviewed migration {cls.split(':', 1)[1]}")
        print("        (v1 and v2 both hold these frames, so each is counted twice)")
    sys.exit(0)

print(f"DIGEST DIFF: {len(bad)} frame/event mismatch(es) over {n} files / {nf} frames")
for cls, cnt in classified.most_common():
    mark = (
        "REVIEWED MIGRATION"
        if cls.startswith("migrated:")
        else "ALLOWED"
        if cls in ALLOW
        else "NOT ALLOWED"
    )
    print(f"  {cls:16s} {cnt:5d}   {mark}")
print()
for cls, line in bad[:40]:
    print("  " + line)
if len(bad) > 40:
    print(f"  ... and {len(bad) - 40} more")

# A reviewed migration is neither a mismatch nor an --allow: it was examined frame by
# frame and written down, so it is expected by name and cannot be widened by a flag.
unexpected = {c for c in classified if c not in ALLOW and not c.startswith("migrated:")}
if not unexpected:
    print(f"\nall mismatches are in --allow {sorted(ALLOW)}; treating as intended")
    print(
        "NOTE: an allowed class still needs chemical judgement -- this only confirms the "
        "SHAPE of the diff, not that the new values are right."
    )
    sys.exit(0)
print(f"\nunexpected classes: {sorted(unexpected)}")
if "graph_invariant" in unexpected:
    print(
        "  graph_hash / fragments moved: the molecular graph changed. P1 and P2 are "
        "frozen on these."
    )
sys.exit(1)
