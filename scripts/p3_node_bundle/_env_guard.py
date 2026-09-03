"""Precision, noise and environment guard. Every heavy P3 script imports this FIRST.

Import order matters: this module changes how torch.load behaves, which must happen before
e3nn or mace are imported. So `from _env_guard import guard` goes above
`from prrs.calculators import ...` in every script here, and that is not a style preference.

Three things it handles.

1. The torch.load noise, fixed at the cause rather than filtered. torch 2.6 flipped the
   default of `weights_only` in torch.load to True. Two things then break:

     e3nn/o3/_wigner.py loads constants.pt, which contains a pickled `slice`, not an
     allowed global -- so `import e3nn.o3` raises UnpicklingError outright.
     mace/calculators/mace.py:226 calls `torch.load(f=model_path, map_location=device)`
     for a full pickled model object, which weights-only loading cannot reconstruct.

   The workaround is TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD=1, and mace sets it ITSELF -- see
   mace/__init__.py:5 and mace/calculators/mace.py:15. That is where the noise comes from:
   torch/serialization.py:1499 warns on EVERY torch.load where the variable is set and the
   callsite did not pass `weights_only` itself. At 8 workers that is twice per worker start.

   Because mace re-asserts the variable on import, unsetting it does not work -- an earlier
   version of this guard popped it and the warning came back the moment mace was imported.
   So the value of the variable is not depended on at all. Instead:

     `slice` is registered as a safe global, so e3nn's constants.pt loads under
     weights_only=True for real -- it genuinely is a weights-only load.
     torch.load is wrapped to pass `weights_only` EXPLICITLY on every call: False for the
     one file whose sha256 this module verified, torch's safe default otherwise. Since the
     callsite now always sets the argument, torch's warning condition
     (`weights_only_not_set`) is permanently false and nothing is emitted.

   The result is narrower than the environment variable, not just quieter: arbitrary
   pickles keep torch's safe default, the single exemption is bound to the checkpoint hash,
   and a future dependency that really needs weights_only=False will raise an
   UnpicklingError naming itself instead of being silently permitted.

2. TF32 is forced OFF. TF32 accelerates fp32 matmuls; the whole PRRS stack runs
   default_dtype="float64", which TF32 does not touch, so enabling it buys nothing.
   Dropping to fp32 to reach it would be destructive: frozen results carry rigid-mode
   floor residuals near 1e-16, the analytic Hessian is relied on to machine precision, and
   the truncation analysis in PRRS_STATUS.md 2.6f is written for float64. An fp32 run would
   look completely normal and be worthless. torch ships backends.cudnn.allow_tf32 = True by
   DEFAULT, so this is a real change, not a no-op.

3. Identity of the run is recorded, including hostname -- the field PRRS_STATUS.md 9 lists
   as missing from the search manifest, which is why cross-host liveness had to be guessed.
"""
import hashlib
import json
import os
import platform
import socket
import sys
from pathlib import Path

# Every checkpoint this bundle is allowed to load, by sha256 prefix, with the model_type
# MACECalculator needs. Unknown checkpoints are refused rather than run: a number produced
# by an unrecorded model cannot be compared with anything.
#
# model_type matters and is not cosmetic. POLAR-1 must be constructed with
# model_type="PolarMACE" or inference never reaches the polarisation-field path, and the
# model then silently behaves like a plain MACE -- see HANDOFF.md 5.
KNOWN_MODELS = {
    "e5ccf5837f685899": {"name": "MACE-OFF24_medium", "model_type": "MACE",
                         "charge_aware": False, "atomic_charges": False,
                         "note": "P0-P2 all used this one; the frozen reference"},
    "1876b02448337595": {"name": "MACE-omol-0-xl-4M", "model_type": "MACE",
                         "charge_aware": True, "atomic_charges": False,
                         "note": "responds to total charge, emits none per atom"},
    "fab8b8713c832f31": {"name": "MACE-POLAR-1-M", "model_type": "PolarMACE",
                         "charge_aware": True, "atomic_charges": True,
                         "note": "needs graph_longrange 0.4.0; emits dipole and Qs"},
}
FROZEN_MODEL_PREFIX = "e5ccf5837f685899"


def model_info(model_path):
    """Identify a checkpoint by hash. Refuses anything not in KNOWN_MODELS."""
    model = Path(model_path)
    if not model.exists():
        raise SystemExit(f"model not found: {model}")
    digest = hashlib.sha256(model.read_bytes()).hexdigest()
    for prefix, info in KNOWN_MODELS.items():
        if digest.startswith(prefix):
            return dict(info, sha256_16=digest[:16], path=str(model))
    raise SystemExit(
        f"checkpoint {model.name} has sha256 {digest[:16]}, which is not in KNOWN_MODELS. "
        f"Refusing: a result from an unrecorded model cannot be compared with anything. "
        f"Add it to KNOWN_MODELS with its model_type if it is intended.")


def shard_arg(argv, flag="--shard"):
    """Parse `--shard i/n` out of argv, returning (i, n) and the argv with it removed.

    Sharding rather than a thread pool because the workers are separate processes: the
    measured contention (PRRS_STATUS.md 8.17, scripts/_saturation_probe.py) is between
    CUDA contexts, and one process per shard is what the P2 scan already does.
    """
    rest, i, n = [], 0, 1
    k = 0
    while k < len(argv):
        if argv[k] == flag and k + 1 < len(argv):
            i, n = (int(x) for x in argv[k + 1].split("/"))
            k += 2
            continue
        rest.append(argv[k])
        k += 1
    if not (0 <= i < n):
        raise SystemExit(f"bad shard {i}/{n}")
    return i, n, rest


def _quiet_torch_load(trusted_sha_prefix, trusted_paths):
    """Make torch.load silent and safe without an environment variable.

    Registers `slice` so e3nn's constants.pt is a genuine weights-only load, then wraps
    torch.load so that ONLY the hash-verified checkpoint gets an explicit
    weights_only=False. Explicit is the whole point: torch warns only when the variable is
    set and the callsite left the argument unset, so setting it at the callsite removes the
    warning without suppressing anything.
    """
    import torch

    # `slice` first: with it registered, e3nn's constants.pt is a genuine weights-only load,
    # so the wrapper below can hand it torch's safe default instead of an exemption.
    torch.serialization.add_safe_globals([slice])
    if getattr(torch.load, "_prrs_wrapped", False):
        return
    original = torch.load
    trusted = {str(Path(p).resolve()) for p in trusted_paths}

    def load(*args, **kwargs):
        """Always pass weights_only EXPLICITLY, which is what silences torch.

        torch/serialization.py warns only when TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD is set AND
        the callsite left `weights_only` unset. Setting it here makes that condition false
        permanently, so the warning cannot come back no matter who sets the variable.

        And someone does keep setting it: mace sets it itself, in mace/__init__.py:5 and
        mace/calculators/mace.py:15. An earlier version of this guard popped the variable
        instead, on the theory that it was inherited from the launching shell -- that was
        wrong. Popping worked right up until `import mace`, which put it straight back, so
        the warning returned from inside this very wrapper. Removing a variable that a
        dependency re-asserts on import is not a fix; not depending on its value is.

        The trusted checkpoint is loaded as a full pickle because its sha256 was verified
        against the frozen model a few lines above. Everything else gets torch's safe
        default, so a future dependency that genuinely needs weights_only=False will raise
        an UnpicklingError naming itself rather than being silently permitted.
        """
        target = kwargs.get("f", args[0] if args else None)
        if "weights_only" not in kwargs:
            resolved = None
            if isinstance(target, (str, Path)):
                try:
                    resolved = str(Path(target).resolve())
                except OSError:
                    resolved = None
            kwargs["weights_only"] = resolved not in trusted if resolved else True
        return original(*args, **kwargs)

    load._prrs_wrapped = True
    torch.load = load


def shard_hint():
    return os.environ.get("CUDA_VISIBLE_DEVICES", "<all>")


def guard(model_path, label="", models=None):
    """Verify the environment and the checkpoint(s), and report what was actually used.

    `models` accepts several paths for a cross-model comparison; `model_path` stays for the
    single-model scripts. All of them are identified by hash and all are trusted for the
    full-pickle load, since each was checked against KNOWN_MODELS.
    """
    paths = [Path(p) for p in (models if models is not None else [model_path])]
    infos = [model_info(p) for p in paths]

    import numpy as np
    import torch

    # Only after the hashes are verified: the exemption below is bound to that check.
    _quiet_torch_load(None, paths)

    was_on = (torch.backends.cuda.matmul.allow_tf32, torch.backends.cudnn.allow_tf32)
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False
    torch.set_float32_matmul_precision("highest")
    if any(was_on):
        # Stated once, not silently corrected: whether TF32 was on decides whether a run is
        # comparable with the frozen ones. cudnn.allow_tf32 is True by default in torch.
        print(f"note: TF32 was on at import (matmul={was_on[0]}, cudnn={was_on[1]}); "
              f"forced off. float64 is unaffected either way.", flush=True)
    if not torch.cuda.is_available():
        raise SystemExit("no CUDA device visible")
    probe = torch.zeros(1, dtype=torch.float64, device="cuda")
    assert probe.dtype == torch.float64

    src = Path(__file__).resolve().parents[2] / "src" / "prrs"
    impl = hashlib.sha256()
    for path in sorted(src.glob("*.py")):
        impl.update(path.name.encode())
        impl.update(path.read_bytes())

    env = {
        "label": label,
        "hostname": socket.gethostname(),
        "platform": platform.platform(),
        "python": sys.version.split()[0],
        "torch": torch.__version__,
        "numpy": np.__version__,
        "cuda_visible": shard_hint(),
        "cuda_device": torch.cuda.get_device_name(0),
        "omp_threads": os.environ.get("OMP_NUM_THREADS", "<unset>"),
        "tf32_matmul": torch.backends.cuda.matmul.allow_tf32,
        "tf32_cudnn": torch.backends.cudnn.allow_tf32,
        "tf32_was_on_at_import": list(was_on),
        "default_dtype": "float64",
        "torch_force_no_weights_only_env": os.environ.get(
            "TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD", "<removed by guard; not needed>"),
        "models": [{k: v for k, v in i.items() if k != "path"} for i in infos],
        "model_sha256_16": infos[0]["sha256_16"],
        "implementation_sha256": impl.hexdigest(),
    }
    print("env " + json.dumps(env), flush=True)
    return env
