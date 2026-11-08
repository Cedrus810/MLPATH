"""Measure: does float32 change Hessian signs / negative-mode counts on soft-mode structures?

Structures (POLAR-1 stationary points re-optimized 2026-10-08 in float64, plus the b09
complex minimum): TS-2OMe, TS-2OMe-cis, Dimethoxycarbene, 1-Cyclohexene-1-carboxaldehyde,
probe_b09 input complex, and b08's TS-CCl2-C4 (-4 cm^-1, nearly flat) and TS-CCl2.
`run DTYPE OUT` evaluates, per structure, prrs `hessian_spectrum` (12 lowest) and
`confirm_minimum` (saddle order, internal recheck) in one dtype per process.
`compare OUT` writes compare.json. Measurement only; no threshold is introduced.
"""

import json
import sys
from pathlib import Path

import numpy as np
from ase.io import read

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "explore" / "src"), str(ROOT)]

from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402
from prrs.io import snapshot  # noqa: E402
from prrs.reliability import GuardedCalculator  # noqa: E402
from prrs.runner import confirm_minimum, hessian_spectrum  # noqa: E402
from prrs.search import _stamped_with_charge  # noqa: E402

OPT = ROOT / "benchmarks/carbene_ambimodal/results/2026-10-08/opt"
STRUCTURES = {
    "TS-2OMe": OPT / "TS-2OMe/geometry.extxyz",
    "TS-2OMe-cis": OPT / "TS-2OMe-cis/geometry.extxyz",
    "Dimethoxycarbene": OPT / "Dimethoxycarbene/geometry.extxyz",
    "Cyclohexene-carboxaldehyde": OPT / "1-Cyclohexene-1-carboxaldehyde_6O/geometry.extxyz",
    "b09-complex": ROOT / "benchmarks/carbene_ambimodal/probe_b09/input.extxyz",
    "TS-CCl2-C4": OPT / "TS-CCl2-C4/geometry.extxyz",
    "TS-CCl2": OPT / "TS-CCl2/geometry.extxyz",
}
CONFIG = SearchConfig(charge_sensitive=True, total_charge=0, multiplicity=1)


def plain(x):
    if isinstance(x, dict):
        return {k: plain(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [plain(v) for v in x]
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, np.generic):
        return x.item()
    return x


def run(dtype, out):
    out.mkdir(parents=True, exist_ok=True)
    calc = polar1("cuda", dtype)

    def factory():
        return calc

    rec = {"dtype": dtype, "structures": {}}
    for name, path in STRUCTURES.items():
        atoms = _stamped_with_charge(read(path), CONFIG)
        entry = {}
        try:
            snap = snapshot(atoms, GuardedCalculator(factory(), CONFIG))
            entry["fmax_eV_A"] = float(np.max(np.linalg.norm(snap.get_forces(), axis=1)))
            eig = np.asarray(hessian_spectrum(snap, factory, CONFIG))
            entry["lowest12"] = eig[:12].tolist()
            entry["n_below_tol"] = int(np.sum(eig < -CONFIG.minimum_check_eigenvalue_tol))
            ok, checked = confirm_minimum(snap, factory, CONFIG, seed=0)
            entry["confirm"] = {
                "is_minimum": bool(ok),
                **{
                    k: checked.get(k)
                    for k in (
                        "saddle_order",
                        "reason",
                        "imaginary_wavenumbers_icm",
                        "lowest_eigenvalues",
                        "internal_mode_recheck",
                    )
                },
            }
        except Exception as exc:  # recorded, not hidden
            entry["error"] = repr(exc)
        rec["structures"][name] = plain(entry)
        (out / f"run_{dtype}.json").write_text(json.dumps(rec, indent=2, default=str))


def compare(out):
    a, b = (json.loads((out / f"run_{d}.json").read_text()) for d in ("float64", "float32"))
    res = {}
    for name in STRUCTURES:
        x, y = a["structures"][name], b["structures"][name]
        if "error" in x or "error" in y:
            res[name] = {"error64": x.get("error"), "error32": y.get("error")}
            continue
        e64, e32 = np.array(x["lowest12"]), np.array(y["lowest12"])
        res[name] = {
            "fmax_eV_A": [x["fmax_eV_A"], y["fmax_eV_A"]],
            "n_below_tol": [x["n_below_tol"], y["n_below_tol"]],
            "saddle_order": [x["confirm"]["saddle_order"], y["confirm"]["saddle_order"]],
            "is_minimum": [x["confirm"]["is_minimum"], y["confirm"]["is_minimum"]],
            "reason": [x["confirm"]["reason"], y["confirm"]["reason"]],
            "imag_icm": [
                x["confirm"]["imaginary_wavenumbers_icm"],
                y["confirm"]["imaginary_wavenumbers_icm"],
            ],
            # rigid-body modes are ~1e-15 with random sign; count flips only where the
            # fp64 eigenvalue is outside the existing confirm_minimum tolerance
            "sign_flips_beyond_tol": int(
                np.sum(
                    (np.abs(e64) > CONFIG.minimum_check_eigenvalue_tol)
                    & (np.sign(e64) != np.sign(e32))
                )
            ),
            "n_soft_within_10tol_fp64": int(
                np.sum(np.abs(e64) < 10 * CONFIG.minimum_check_eigenvalue_tol)
            ),
            "max_abs_eig_diff_lowest12": float(np.max(np.abs(e64 - e32))),
            "lowest12_fp64": e64.tolist(),
            "lowest12_fp32": e32.tolist(),
        }
    res["tol"] = CONFIG.minimum_check_eigenvalue_tol
    (out / "compare.json").write_text(json.dumps(res, indent=2, default=str))


if __name__ == "__main__":
    if sys.argv[1] == "run":
        run(sys.argv[2], Path(sys.argv[3]))
    else:
        compare(Path(sys.argv[2]))
