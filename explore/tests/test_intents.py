import numpy as np
from molecules import chloride_methyl_chloride, ethanol, relabel
from prrs.config import SearchConfig

from explore.intents import IntentSpec, distance_intents, pair_intents
from explore.probe import Intent, _feasible

CONFIG = SearchConfig(closed_shell_only=False)


def _mapped(specs, perm):
    """Intent identity in ORIGINAL atom indices: perm[k] is the original index of new atom k.

    Each intent is normalised: (pair, coef) tuples sorted, and if the first coef is
    negative every coef and the sign are flipped ((1,-1,s) == (-1,1,-s))."""
    out = set()
    for s in specs:
        pairs = tuple(tuple(sorted((int(perm[i]), int(perm[j])))) for i, j in s.intent.pairs)
        items = sorted(zip(pairs, s.intent.coefs))
        sign = s.sign
        if items[0][1] < 0:
            items = [(p, -c) for p, c in items]
            sign = -sign
        out.add((tuple(items), sign, s.layer))
    return out


def _plain(specs):
    return _mapped(specs, np.arange(100))


def test_bonds_get_both_signs_and_close_intramolecular_pairs_are_excluded():
    atoms = ethanol()
    singles, report = distance_intents(atoms, "p", CONFIG, n_far=0, seed=0)
    bonded = [s for s in singles if s.layer == "bonded"]
    assert len(bonded) == 2 * 8  # ethanol has 8 bonds
    assert report["excluded"]["bend_in_disguise"] > 0
    assert all(s.family == "distance" and len(s.intent.pairs) == 1 for s in singles)


def test_cross_fragment_pairs_are_scoped_inter():
    atoms = chloride_methyl_chloride()
    singles, _ = distance_intents(atoms, "p", CONFIG, n_far=0, seed=0)
    inter = [s for s in singles if s.scope == "inter"]
    assert inter and all(5 in s.intent.pairs[0] for s in inter)


def test_far_layer_is_sampled_by_colour_class_and_counted():
    atoms = chloride_methyl_chloride()
    atoms.positions[5] = [0, 0, -9.0]  # push Cl- beyond pair_cutoff_A
    none, _ = distance_intents(atoms, "p", CONFIG, n_far=0, seed=0)
    some, rep1 = distance_intents(atoms, "p", CONFIG, n_far=1, seed=0)
    assert not [s for s in none if s.layer == "far"]
    far = [s for s in some if s.layer == "far"]
    assert far and rep1["far_classes_total"] >= rep1["far_classes_sampled"] == 1
    # a sampled colour class comes in whole: the three H...Cl- pairs together
    hcl = [s for s in far if {int(atoms.numbers[i]) for i in s.intent.pairs[0]} == {1, 17}]
    assert len(hcl) in (0, 6)  # whole-class property is pinned exactly in the relabel test


def test_relabelling_maps_the_intent_set():
    atoms = ethanol()
    singles, _ = distance_intents(atoms, "p", CONFIG, n_far=2, seed=0)
    moved, perm = relabel(atoms, seed=11)
    msingles, _ = distance_intents(moved, "p", CONFIG, n_far=2, seed=0)
    assert _mapped(msingles, perm) == _plain(singles)
    combos, _ = pair_intents(singles, atoms, CONFIG, "adjacent", "p")
    mcombos, _ = pair_intents(msingles, moved, CONFIG, "adjacent", "p")
    assert _mapped(mcombos, perm) == _plain(combos)


def test_localities_are_nested_and_ids_are_unique():
    atoms = ethanol()
    singles, _ = distance_intents(atoms, "p", CONFIG, n_far=0, seed=0)
    counts = {}
    for loc in ("shared_atom", "adjacent", "all"):
        combos, rep = pair_intents(singles, atoms, CONFIG, loc, "p")
        ids = [c.intent.intent_id for c in combos]
        assert len(ids) == len(set(ids))
        assert all(len(c.intent.pairs) == 2 and c.family == "distance_pair" for c in combos)
        counts[loc] = len(combos)
    assert counts["shared_atom"] <= counts["adjacent"] <= counts["all"]


def test_same_sign_pair_below_overlap_floor_is_infeasible():
    numbers = np.array([6, 6, 6])
    intent = Intent("i", ((0, 1), (1, 2)), (1.0, 1.0))
    floor = 2 * CONFIG.min_distance_ratio * (2 * 0.76)
    assert not _feasible(intent, floor + 0.5, (0.0, -1.0), numbers, CONFIG)
    assert _feasible(intent, floor + 0.5, (0.0, -0.4), numbers, CONFIG)


def test_mixed_sign_pair_is_not_bounded_by_the_floor():
    numbers = np.array([6, 6, 6, 6])
    intent = Intent("i", ((0, 1), (2, 3)), (1.0, -1.0))
    assert _feasible(intent, 0.0, (0.0, -5.0, 5.0), numbers, CONFIG)


def test_pair_intent_ids_keep_the_full_pose_id_when_it_contains_bars():
    atoms = ethanol()
    ids = {}
    for k in (0, 1):
        pose = f"f0|C2H6O|x:t0:L0-v1+f1|C2H6O|y:k{k}:f2:L0-v1"
        singles, _ = distance_intents(atoms, pose, CONFIG, n_far=0, seed=0)
        combos, _ = pair_intents(singles, atoms, CONFIG, "adjacent", pose)
        assert combos and all(c.intent.intent_id.startswith(pose + "|") for c in combos)
        ids[k] = {c.intent.intent_id for c in combos}
    assert not ids[0] & ids[1]


def test_far_sampling_is_relabel_invariant_and_takes_whole_classes():
    atoms = chloride_methyl_chloride()
    atoms.positions[5] = [0, 0, -9.0]
    moved, perm = relabel(atoms, seed=5)
    sizes = set()
    for seed in range(8):
        singles, rep = distance_intents(atoms, "p", CONFIG, n_far=1, seed=seed)
        assert rep["far_classes_total"] == 3 > rep["far_classes_sampled"] == 1
        msingles, _ = distance_intents(moved, "p", CONFIG, n_far=1, seed=seed)
        far = [s for s in singles if s.layer == "far"]
        assert _mapped([s for s in msingles if s.layer == "far"], perm) == _plain(far)
        # classes: C...Cl- (1 pair), Cl...Cl- (1 pair), H...Cl- (3 pairs); both signs each
        assert len(far) in (2, 6)
        sizes.add(len(far))
    assert sizes == {2, 6}  # the three-pair class does get drawn, and then comes whole


def test_deltas_apply_the_sign_and_start_at_zero():
    base = Intent("i", ((0, 1),), (1.0,))
    minus = IntentSpec(base, -1, "distance", "bonded", "inter")
    plus = IntentSpec(base, 1, "distance", "bonded", "inter")
    assert minus.deltas((0.5, 1.0, 1.5)) == (0.0, -0.5, -1.0, -1.5)
    assert plus.deltas((0.5, 1.0, 1.5)) == (0.0, 0.5, 1.0, 1.5)
    assert minus.deltas((-0.5,)) == (0.0, -0.5)
