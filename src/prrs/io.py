"""Atomic manifests, append-only attempt records, physical snapshots."""

import hashlib
import json
import math
import os
from pathlib import Path
import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator


def nonfinite_paths(data, path=""):
    """Every location in a nested record holding a non-finite float.

    `allow_nan=False` is the right guard -- a NaN or an infinity in a record is a defect,
    not a value -- but on its own it raises "Out of range float values are not JSON
    compliant: -inf" and names nothing. Three seeds of the T=300 scan died on exactly that
    message, and because the write that failed was publish(), the run could not even record
    that it had failed: status stayed "running", which reads like being killed by a signal.
    That cost a wrong diagnosis before the log was read. So the guard now says where.
    """
    found = []
    if isinstance(data, dict):
        for key, value in data.items():
            found += nonfinite_paths(value, f"{path}.{key}" if path else str(key))
    elif isinstance(data, (list, tuple)):
        for index, value in enumerate(data):
            found += nonfinite_paths(value, f"{path}[{index}]")
    elif isinstance(data, float) and not math.isfinite(data):
        found.append((path or "<root>", data))
    return found


def dumps_strict(data):
    """json.dumps with allow_nan=False, and an error that names the offending fields."""
    try:
        return json.dumps(data, allow_nan=False)
    except ValueError as exc:
        offenders = nonfinite_paths(data)
        if not offenders:
            raise
        detail = ", ".join(f"{where} = {value}" for where, value in offenders[:8])
        more = "" if len(offenders) <= 8 else f" (and {len(offenders) - 8} more)"
        raise ValueError(
            f"{exc}. Non-finite values at: {detail}{more}. A record may not carry NaN or "
            "infinity; fix the source rather than relaxing this check."
        ) from exc


def atomic_json(path, data):
    path = Path(path)
    temp = path.with_name(path.name + ".tmp")
    # Same strict guard as append_json, and for the same reason: atomic_json is what
    # publish() calls, so an unlocated failure here is the one that takes a whole run down.
    try:
        text = json.dumps(data, indent=2, allow_nan=False)
    except ValueError:
        dumps_strict(data)  # raises with the field names
        raise  # unreachable; kept so the intent is plain
    temp.write_text(text + "\n", encoding="utf-8")
    os.replace(temp, path)


def append_json(path, data):
    with Path(path).open("a", encoding="utf-8") as handle:
        handle.write(dumps_strict(data) + "\n")
        handle.flush()


def structure_hash(atoms):
    payload = {
        "numbers": atoms.numbers.tolist(),
        "positions": atoms.positions.tolist(),
        "masses": atoms.get_masses().tolist(),
        "pbc": atoms.pbc.tolist(),
        "cell": atoms.cell.tolist(),
        "initial_charges": atoms.get_initial_charges().tolist(),
        "initial_magmoms": atoms.get_initial_magnetic_moments().tolist(),
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


def snapshot(atoms, physical):
    energy = physical.get_potential_energy(atoms)
    forces = physical.get_forces(atoms)
    result = atoms.copy()
    result.calc = SinglePointCalculator(result, energy=energy, forces=np.array(forces))
    return result
