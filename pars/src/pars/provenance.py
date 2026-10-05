"""Independent refinement provenance, including actual calculator/model state.

Factories may export ``calculator_metadata()`` or the caller may supply JSON
metadata. File paths in that metadata are hashed every invocation. Loaded
model state_dicts are hashed too, so replacing weights behind the same
factory string does not reuse a result from the old PES.
"""

import hashlib
import inspect
import json
import platform
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import numpy as np


def file_sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def implementation_sha256():
    digest = hashlib.sha256()
    root = Path(__file__).parent
    for path in sorted(root.rglob("*.py")):
        digest.update(str(path.relative_to(root)).encode())
        digest.update(path.read_bytes())
    return digest.hexdigest()


def numerical_versions():
    result = {"python": platform.python_version()}
    for name in ("numpy", "scipy", "ase", "torch", "mace-torch", "openmm"):
        try:
            result[name] = version(name)
        except PackageNotFoundError:
            result[name] = None
    return result


def calculator_identity(factory, calculator, metadata=None):
    files = {}

    def source(obj):
        try:
            path = inspect.getsourcefile(obj)
        except TypeError:
            path = None
        return file_sha256(path) if path and Path(path).is_file() else None

    def freeze(value):
        if value is None or isinstance(value, (bool, int, float)):
            return value
        if isinstance(value, (str, Path)):
            text = str(value)
            try:
                path = Path(text)
                if path.is_file():
                    files[str(path.resolve())] = file_sha256(path)
            except (OSError, ValueError):
                pass
            return text
        if isinstance(value, np.generic):
            return value.item()
        if isinstance(value, np.ndarray) or hasattr(value, "detach"):
            array = value if isinstance(value, np.ndarray) else value.detach().cpu().numpy()
            return {
                "shape": list(array.shape),
                "dtype": str(array.dtype),
                "sha256": hashlib.sha256(array.tobytes()).hexdigest(),
            }
        if isinstance(value, (tuple, list)):
            return [freeze(v) for v in value]
        if isinstance(value, dict):
            return {
                str(k): freeze(v) for k, v in sorted(value.items(), key=lambda kv: str(kv[0]))
            }
        if hasattr(value, "state_dict"):
            return {"type": type(value).__qualname__, "weights": freeze(value.state_dict())}
        if type(value).__module__ == "torch" and type(value).__name__ in ("device", "dtype"):
            return str(value)
        if hasattr(value, "__dict__"):
            # ASE caches change while computing; neither coordinates nor
            # last energy/forces define a calculator's frozen configuration.
            state = {
                k: v
                for k, v in vars(value).items()
                if k not in {"atoms", "results", "_directory", "prefix"}
            }
            return {
                "type": f"{type(value).__module__}.{type(value).__qualname__}",
                "source_sha256": source(type(value)),
                "state": freeze(state),
            }
        raise ValueError(
            f"cannot fingerprint calculator state {type(value).__module__}.{type(value).__qualname__}; supply an inspectable calculator"
        )

    if metadata is None:
        module = inspect.getmodule(factory)
        hook = getattr(module, "calculator_metadata", None)
        metadata = hook() if callable(hook) else {}
    if not isinstance(metadata, dict):
        raise ValueError("calculator metadata must be a JSON object")
    result = {
        "factory_source_sha256": source(factory),
        "metadata": freeze(metadata),
        "calculator": freeze(calculator),
        "local_file_sha256": files,
    }
    for name, expected in metadata.get("local_file_sha256", {}).items():
        if not Path(name).is_file() or file_sha256(name) != expected:
            raise ValueError(f"calculator checkpoint missing or changed: {name}")
    return json.loads(json.dumps(result, allow_nan=False))


def condition_groups(records):
    """Seed hashes identify candidates; all other conditions identify a cohort."""
    groups = {}
    for record in records:
        conditions = record.get("conditions")
        if conditions is None:
            raise ValueError(
                "existing refined record has no conditions; redo it before aggregation"
            )
        common = {k: v for k, v in conditions.items() if k != "seed_sha256"}
        key = hashlib.sha256(
            json.dumps(common, sort_keys=True, allow_nan=False).encode()
        ).hexdigest()
        groups.setdefault(
            key, {"conditions_sha256": key, "conditions": common, "candidate_ids": []}
        )["candidate_ids"].append(record["candidate_id"])
        record["conditions_sha256"] = key
    return {"homogeneous": len(groups) <= 1, "groups": [groups[k] for k in sorted(groups)]}
