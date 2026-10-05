"""Standalone I/O helpers.

PARS cannot use ``prrs.io``, so this mirrors its guarantees: JSON is written
atomically (temporary file in the target directory, then ``os.replace``) and
NaN/Inf are rejected at serialization time so a corrupt number can never enter
a run record silently.
"""

from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path


def read_json(path) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def read_jsonl(path) -> list[dict]:
    rows = []
    for line in Path(path).read_text(encoding="utf-8-sig").splitlines():
        line = line.strip()
        if line:
            rows.append(json.loads(line))
    return rows


def _atomic_write(path, text: str) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", dir=destination.parent, delete=False, suffix=".tmp"
    )
    try:
        handle.write(text)
        handle.flush()
        os.fsync(handle.fileno())
        handle.close()
        os.replace(handle.name, destination)
    except BaseException:
        handle.close()
        Path(handle.name).unlink(missing_ok=True)
        raise


def write_json_atomic(path, payload) -> None:
    _atomic_write(path, json.dumps(payload, indent=2, allow_nan=False) + "\n")


def write_jsonl_atomic(path, rows) -> None:
    text = "".join(json.dumps(row, allow_nan=False) + "\n" for row in rows)
    _atomic_write(path, text)
