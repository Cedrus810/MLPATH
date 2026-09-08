"""Step 0 for any ionic chemistry: does the model respond to total charge at all?

MACE-OFF24 does not. Measured on the documented channel (`atoms.info["charge"]`, which the
calculator maps to the model's total_charge input): the energy is bit-identical at charge
0, -1 and +1, and the loaded model has no module whose name mentions charge or spin. That
is not an approximation, it is an absence -- the model cannot tell Cl- from Cl radical, so
every SN1/SN2/E2 path is one object to it.

Before designing a single acceptance criterion for ionic chemistry, the replacement model
has to be shown to move when the knob moves. This is the same discipline the temperature
scan taught the hard way: a variance experiment designed before checking that the varied
knob perturbs the system buys an answer that was already determined.

Reports, per model and per structure:
  * whether E changes with total charge, and by how much
  * whether E changes with total spin
  * whether the model declares charge/spin machinery at all
  * whether an analytic Hessian is available and agrees with finite differences,
    because every energy-resolution argument in this project rests on that
"""
import json
import sys
import time
from pathlib import Path
import numpy as np
from ase import Atoms
from ase.io import read

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from prrs.calculators import load_factory                            # noqa: E402
from prrs.config import SearchConfig                                 # noqa: E402
from prrs.runner import hessian_vector_product                       # noqa: E402
from prrs.reliability import GuardedCalculator                       # noqa: E402

MACE = Path("/home/ruigengji/MLP/mace")
MODELS = [
    ("MACE-OFF24_medium", MACE / "MACE-OFF24_medium.model"),
    ("MACE-omol-0-xl-4M", MACE / "mace-omol-0-extra-large-4M.model"),
    ("MACE-POLAR-1-M", MACE / "MACE-POLAR-1-M.model"),
]

def structures():
    """One neutral closed-shell molecule and one genuinely ionic species."""
    yield "3-oxobutanal enol (neutral, the P2 reactant)", read(
        str(Path("runs/p2_oxobutanal/p2_input.extxyz")))
    # Chloromethane at a sane geometry, and the chloride ion on its own. The second is the
    # species OFF24 cannot represent and the whole of SN2 depends on.
    yield "CH3Cl (neutral)", Atoms("CClH3", positions=[
        [0.000, 0.000, 0.000], [0.000, 0.000, 1.785],
        [1.028, 0.000, -0.373], [-0.514, 0.890, -0.373], [-0.514, -0.890, -0.373]])
    yield "Cl- (a bare chloride)", Atoms("Cl", positions=[[0.0, 0.0, 0.0]])


def charge_response(calc, atoms, charges=(0, -1, 1), spins=(1, 2)):
    out = {}
    for q in charges:
        probe = atoms.copy()
        probe.info["charge"] = q
        probe.info["spin"] = 1
        probe.calc = calc
        calc.results.clear()
        try:
            out[f"charge {q:+d}"] = float(probe.get_potential_energy())
        except Exception as exc:                                     # noqa: BLE001
            out[f"charge {q:+d}"] = f"error: {type(exc).__name__}: {exc}"[:90]
    for s in spins:
        probe = atoms.copy()
        probe.info["charge"] = 0
        probe.info["spin"] = s
        probe.calc = calc
        calc.results.clear()
        try:
            out[f"spin {s}"] = float(probe.get_potential_energy())
        except Exception as exc:                                     # noqa: BLE001
            out[f"spin {s}"] = f"error: {type(exc).__name__}: {exc}"[:90]
    return out


def hessian_check(calc, atoms, config):
    """Is there an analytic Hessian, and does one column agree with finite differences?"""
    provider = getattr(calc, "get_hessian", None)
    if provider is None:
        return {"analytic": False, "reason": "the calculator exposes no get_hessian"}
    probe = atoms.copy()
    probe.info.setdefault("charge", 0)
    probe.info.setdefault("spin", 1)
    probe.calc = calc
    try:
        started = time.perf_counter()
        raw = np.asarray(provider(probe), dtype=float)
        seconds = time.perf_counter() - started
    except Exception as exc:                                         # noqa: BLE001
        return {"analytic": False, "reason": f"{type(exc).__name__}: {exc}"[:110]}
    size = 3 * len(probe)
    matrix = raw.reshape(size, size)
    direction = np.zeros((len(probe), 3))
    direction[0, 0] = 1.0
    guard = GuardedCalculator(calc, config)
    numeric = hessian_vector_product(probe, guard, direction,
                                     config.minimum_check_step_A).ravel()
    worst = float(np.abs(matrix[:, 0] - numeric).max())
    return {"analytic": True, "seconds": round(seconds, 3),
            "shape": list(matrix.shape),
            "symmetry": float(np.abs(matrix - matrix.T).max()),
            "worst_column_deviation_vs_fd": worst,
            "fd_step_A": config.minimum_check_step_A}


def main():
    config = SearchConfig(closed_shell_only=False)
    report = {}
    for name, path in MODELS:
        if not path.exists():
            print(f"== {name}: file missing at {path}"); continue
        print(f"\n{'='*78}\n{name}  ({path.stat().st_size/1e6:.0f} MB)")
        try:
            calc = load_factory("prrs.calculators:mace_factory",
                                {"model_paths": str(path), "device": "cuda",
                                 "default_dtype": "float64"})()
        except Exception as exc:                                     # noqa: BLE001
            print(f"   could not load: {type(exc).__name__}: {exc}")
            report[name] = {"load_error": str(exc)[:200]}
            continue
        model = calc.models[0]
        machinery = [n for n, _ in model.named_modules()
                     if "charge" in n.lower() or "spin" in n.lower()]
        heads = getattr(calc, "available_heads", None)
        zs = getattr(model, "atomic_numbers", None)
        print(f"   heads={heads}  r_max={float(getattr(model,'r_max',float('nan'))):.1f} A")
        print(f"   elements={None if zs is None else sorted(int(z) for z in zs)}")
        print(f"   modules naming charge/spin: {len(machinery)}"
              + (f"  e.g. {machinery[:3]}" if machinery else "  <-- none: cannot consume them"))
        entry = {"heads": heads, "charge_spin_modules": len(machinery),
                 "r_max_A": float(getattr(model, "r_max", float("nan"))),
                 "elements": None if zs is None else sorted(int(z) for z in zs),
                 "structures": {}}
        for label, atoms in structures():
            values = charge_response(calc, atoms)
            numbers = [v for v in values.values() if isinstance(v, float)]
            spread = (max(numbers) - min(numbers)) if len(numbers) > 1 else 0.0
            responds = spread > 1e-9
            print(f"   {label}")
            for key, value in values.items():
                print(f"      {key:<11} "
                      + (f"{value:.9f} eV" if isinstance(value, float) else value))
            print(f"      -> spread {spread:.6f} eV   "
                  f"{'RESPONDS to charge/spin' if responds else 'NO RESPONSE'}")
            entry["structures"][label] = {"values": values, "spread_eV": spread,
                                          "responds": bool(responds)}
        first = next(structures())[1]
        entry["hessian"] = hessian_check(calc, first, config)
        print(f"   analytic Hessian: {entry['hessian']}")
        report[name] = entry
        del calc
    Path("runs/charge_knob.json").write_text(json.dumps(report, indent=1), encoding="utf-8")
    print("\nwrote runs/charge_knob.json")


if __name__ == "__main__":
    main()
