"""Calculator identity for resume; no calculator is executed by this module.

Explicit metadata is the contract for external resources. Local files named
in metadata, defaults or closures are hashed by content, not just by name.
Opaque captured objects make the identity unverifiable and prevent resume.
"""

import functools
import hashlib
import inspect
import json
import marshal
from pathlib import Path
from types import ModuleType

import numpy as np


def calculator_identity(factory, metadata):
    files, opaque, seen = {}, [], set()

    def freeze(value):
        if value is None or isinstance(value, (bool, int, float)):
            return value
        if isinstance(value, (str, Path)):
            text = str(value)
            try:
                path = Path(text)
                if path.is_file():
                    files[str(path.resolve())] = hashlib.sha256(path.read_bytes()).hexdigest()
            except (OSError, ValueError):
                pass  # strings need not be filesystem paths
            return text
        if isinstance(value, np.ndarray):
            return {"array": value.tolist(), "dtype": str(value.dtype)}
        if isinstance(value, np.generic):
            return value.item()
        if isinstance(value, (tuple, list)):
            return [freeze(v) for v in value]
        if isinstance(value, dict):
            return {
                str(k): freeze(v) for k, v in sorted(value.items(), key=lambda kv: str(kv[0]))
            }
        if isinstance(value, ModuleType):
            return {"module": value.__name__, "source": source(value)}
        if callable(value):
            return callable_record(value)
        opaque.append(f"{type(value).__module__}.{type(value).__qualname__}")
        return {"opaque": opaque[-1]}

    def source(value):
        try:
            path = inspect.getsourcefile(value)
        except TypeError:
            path = None
        return (
            hashlib.sha256(Path(path).read_bytes()).hexdigest()
            if path and Path(path).is_file()
            else None
        )

    def callable_record(value):
        if isinstance(value, functools.partial):
            return {
                "partial": callable_record(value.func),
                "args": freeze(value.args),
                "kwargs": freeze(value.keywords),
            }
        record = {
            "module": getattr(value, "__module__", type(value).__module__),
            "qualname": getattr(value, "__qualname__", type(value).__qualname__),
            "source_sha256": source(value),
        }
        if id(value) in seen:
            return record
        seen.add(id(value))
        if inspect.isfunction(value):
            record["code_sha256"] = hashlib.sha256(marshal.dumps(value.__code__)).hexdigest()
            record["defaults"] = freeze(value.__defaults__)
            record["kwdefaults"] = freeze(value.__kwdefaults__)
            record["closure"] = {
                name: freeze(cell.cell_contents)
                for name, cell in zip(value.__code__.co_freevars, value.__closure__ or ())
            }
        elif not inspect.isclass(value) and not inspect.isbuiltin(value):
            record["state"] = freeze(vars(value))
        return record

    result = {"factory": callable_record(factory), "metadata": freeze(metadata)}
    result["local_file_sha256"] = files
    result["verifiable"] = not opaque
    result["opaque_captures"] = sorted(set(opaque))
    # Check declared hashes as well as re-hashing all local inputs. A caller
    # cannot reuse stale metadata after a checkpoint was replaced.
    for name, expected in metadata.get("local_file_sha256", {}).items():
        path = Path(name)
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise ValueError(f"calculator source file missing or changed: {name}")
    return json.loads(json.dumps(result, allow_nan=False))
