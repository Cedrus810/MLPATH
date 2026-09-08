"""Charge-aware calculator factory for b04, resolved through prrs.cli's --calculator hook.

`prrs.calculators.mace_factory` cannot load these checkpoints: it calls MACECalculator
without `model_type`, so POLAR-1 (whose checkpoint class is PolarMACE) would be loaded
as a plain MACE and would compute the wrong thing WITHOUT RAISING. That is gap 2 in
README.md. Rather than widen mace_factory before anything has exercised the path, the
factory lives here and reaches the search through the documented
`--calculator module:function` hook.

Nothing here sets the charge. The search stamps it from the configuration
(prrs.search._stamped_with_charge) and the calculator reads atoms.info["charge"], so
there is exactly one place that decides it.

    prrs search input.extxyz --output runs/... \
        --calculator benchmarks.b04_sn2_chloride.calculator:polar1 \
        --config benchmarks/b04_sn2_chloride/config.json
"""

from pathlib import Path

MODELS = {
    "polar1": ("/home/ruigengji/MLP/mace/MACE-POLAR-1-M.model", "PolarMACE"),
    "omol0": ("/home/ruigengji/MLP/mace/mace-omol-0-extra-large-4M.model", "MACE"),
}


def _build(name, device="cuda", default_dtype="float64"):
    path, model_type = MODELS[name]
    if not Path(path).is_file():
        raise ValueError(f"checkpoint missing: {path}")
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
    from prrs import torch_guard

    torch_guard.install(trusted=[path])
    from mace.calculators import MACECalculator

    return MACECalculator(
        model_paths=path, model_type=model_type, device=device, default_dtype=default_dtype
    )


def polar1(device="cuda", default_dtype="float64"):
    """MACE-POLAR-1-M. Charge-aware; measured 16.67 eV spread across charge 0/-1/+1."""
    return _build("polar1", device, default_dtype)


def omol0(device="cuda", default_dtype="float64"):
    """MACE-omol-0-4M. Second reading only -- it is off by 10.5 kcal on this barrier."""
    return _build("omol0", device, default_dtype)
