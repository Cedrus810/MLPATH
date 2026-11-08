"""L0 (poses) builder shared by the D2 real-model scripts; generalises run_b04_d2.build_l0.

No model is imported here: the caller passes the (float32) calculator. Fragment electronic
states come from ``fragment_states`` keyed by the formula part of ``Fragment.key``
(``prrs.chemistry._formula``, e.g. ``"C4H9Cl"``); a missing formula raises before any
forward. Fragments with fewer than ``min_qualify_atoms`` atoms are not qualified (the
generated conformer is used, as b04 did for Cl-); that is recorded in the returned info.
"""

from collections import Counter

from explore.ledger import Ledger, view_factory
from explore.poses import (
    conformer_grid,
    conformers_from_qualification,
    fragments,
    qualify_conformers,
    sample_poses,
)


def build_l0(
    atoms,
    config,
    calc64,
    prepare_ledger_path,
    fragment_states,
    shells,
    n_samples,
    seed,
    steps_per_turn,
    conformer_cap,
    min_qualify_atoms,
):
    frags = fragments(atoms, config)
    if len(frags) != 2:
        raise ValueError(f"build_l0 needs exactly 2 fragments, got {len(frags)}")
    by_key = {}
    for f in frags:
        formula = f.key.split("|", 1)[0]
        if formula not in fragment_states:
            raise ValueError(f"fragment_states has no entry for fragment formula {formula!r}")
        by_key[f.key] = dict(fragment_states[formula])
    ledger = Ledger(prepare_ledger_path, arm="P", replay=False)
    chosen, unqualified, conformer_counts = [], [], {}
    try:
        for f in frags:
            confs, _ = conformer_grid(atoms, f, config, steps_per_turn, conformer_cap, seed)
            conformer_counts[f.key] = len(confs)
            if len(f.atoms) < min_qualify_atoms:
                unqualified.append(f.key)
            else:
                rows = qualify_conformers(
                    atoms, confs, view_factory(calc64, ledger), config, ledger, seed, by_key
                )
                confs = conformers_from_qualification(confs, rows)
            if not confs:
                raise ValueError(f"fragment {f.key}: no usable conformer")
            chosen.append(confs[0])
    finally:
        ledger.close()
    poses = sample_poses(atoms, frags, chosen[0], chosen[1], list(shells), n_samples, seed)
    out = []
    for p in poses:
        if p.status != "placed":
            continue
        a = atoms.copy()
        a.positions = p.positions
        out.append((p.pose_id, a))
    info = {
        "prepare_forwards": ledger.spent,
        "unqualified_fragments": unqualified,
        "conformer_counts": conformer_counts,
        "pose_status_counts": dict(Counter(p.status for p in poses)),
    }
    return out, info
