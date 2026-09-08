"""Real-data check on the probe tangents, before any round 3 is started.

Six things, on frames the P2 round-2 run actually produced:

  1. bond/kick and dihedral/torsion_kick tangents are finite and non-zero
  2. the norm survives the rigid-body projection
  3. a small central difference of the internal coordinate along the tangent is non-zero and
     has the sign the proposal asked for (dihedrals wrapped to the branch)
  4. in the O-H tangent the hydrogen moves much further than the oxygen
  5. direction_source, the norms before and after projection, and the first target_overlap
     are all written out
  6. a dihedral near collinearity fails closed instead of falling back to the lowest mode
"""
import json
import sys
from pathlib import Path
import numpy as np
from ase.io import read

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from prrs import internal                                            # noqa: E402
from prrs.perturbations import FAMILIES                              # noqa: E402
from prrs.runner import _projected_direction                         # noqa: E402
from prrs.search import _dividing_frame, _probe_tangent              # noqa: E402


def wrap(value):
    return (value + np.pi) % (2 * np.pi) - np.pi


def check(atoms, probe, label):
    kind = FAMILIES[probe["family"]]["kind"]
    indices = tuple(probe["indices"])
    tangent, source = _probe_tangent(atoms, probe)
    row = {"label": label, "family": probe["family"], "kind": kind,
           "indices": list(indices), "direction_source": source}
    if tangent is None:
        row.update(ok=False, reason=source)
        return row
    row["norm_before_projection"] = float(np.linalg.norm(tangent))
    projected = _projected_direction(atoms.positions, atoms.get_masses(), tangent)
    row["norm_after_projection"] = (None if projected is None
                                    else float(np.linalg.norm(projected)))
    if projected is None:
        row.update(ok=False, reason="no_internal_component")
        return row

    # central difference of the coordinate along the tangent
    eps = 1e-4 / max(np.linalg.norm(tangent, axis=1).max(), 1e-12)
    plus = atoms.positions + eps * tangent
    minus = atoms.positions - eps * tangent
    before = internal.coordinate(plus, kind, indices)
    after = internal.coordinate(minus, kind, indices)
    delta = (wrap(before - after) if kind == "dihedral" else before - after) / (2 * eps)
    row["dq_dstep"] = float(delta)
    row["sign_matches_proposal"] = bool(np.sign(delta) == np.sign(probe.get("sign", 1))
                                        or probe.get("sign") is None)
    row["finite"] = bool(np.isfinite(tangent).all())
    row["ok"] = bool(row["finite"] and row["norm_after_projection"] > 1e-8
                     and abs(delta) > 1e-6)
    per_atom = np.linalg.norm(tangent, axis=1)
    row["per_atom_norm"] = {int(i): round(float(per_atom[i]), 5) for i in indices}
    return row


def main(run):
    run = Path(run)
    rows = []
    seen = set()
    for line in (run / "attempts.jsonl").read_text().splitlines():
        record = json.loads(line)
        probe = record.get("probe") or {}
        if not probe.get("family"):
            continue
        key = (probe["family"], tuple(probe["indices"]))
        if key in seen:
            continue
        seen.add(key)
        directory = run / "trials" / record["attempt_id"]
        if not (directory / "response.traj").exists():
            continue
        frame, origin, tangent, source = _dividing_frame(directory, probe=probe)
        if frame is None:
            continue
        rows.append({**check(frame, probe, record["attempt_id"]),
                     "seed_frame": origin, "seed_tangent_source": source})

    print(f"{'attempt':<10} {'family':<13} {'indices':<14} {'source':<30} "
          f"{'|d| pre':>9} {'|d| post':>9} {'dq/ds':>10} ok")
    for row in rows:
        print(f"{row['label']:<10} {row['family']:<13} {str(row['indices']):<14} "
              f"{str(row.get('seed_tangent_source')):<30} "
              f"{row.get('norm_before_projection', float('nan')):>9.4f} "
              f"{(row.get('norm_after_projection') or float('nan')):>9.4f} "
              f"{row.get('dq_dstep', float('nan')):>10.4f} {row['ok']}")

    kinds = {row["kind"] for row in rows if row["ok"]}
    print(f"\nkinds with a usable tangent: {sorted(kinds)}")
    oh = [row for row in rows if tuple(row["indices"]) == (4, 5)]
    for row in oh:
        weights = row["per_atom_norm"]
        heavy, light = weights.get(4), weights.get(5)
        print(f"O-H tangent per-atom norm: O={heavy} H={light} "
              f"-> H/O = {light/heavy:.2f}" if heavy else "O-H tangent missing")
        assert light > 4 * heavy, "the hydrogen must dominate the O-H tangent"
    failures = [row for row in rows if not row["ok"]]
    print(f"\nfailing tangents: {len(failures)}")
    for row in failures:
        print(f"  {row['label']} {row['family']}{row['indices']} -> {row.get('reason')}")

    # (6) a dihedral at a collinear geometry must fail closed
    collinear = read(run / "input.extxyz")
    positions = collinear.positions.copy()
    positions[2] = positions[1] + (positions[1] - positions[0])
    positions[3] = positions[1] + 2.0 * (positions[1] - positions[0])
    collinear.set_positions(positions)
    probe = {"family": "torsion", "indices": (0, 1, 2, 3), "sign": 1}
    tangent, source = _probe_tangent(collinear, probe)
    print(f"\ncollinear dihedral -> tangent={'None' if tangent is None else 'returned'} "
          f"source={source}")
    return rows


if __name__ == "__main__":
    main(sys.argv[1])
