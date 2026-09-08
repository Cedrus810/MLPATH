"""Gas-phase identity SN2, Cl- + CH3Cl, with charge-aware models.

The first ionic chemistry this project can touch. MACE-OFF24 returns bit-identical energies
for charge 0, -1 and +1 and has no charge machinery at all, so every point on this path was
one object to it.

Two earlier attempts are retracted and the reasons are worth keeping:

  1. A hand-interpolated geometry generator whose umbrella angle flipped mid-scan. One of
     its points came out at +90.9 kcal/mol -- a broken geometry, not a barrier. It measured
     the generator.
  2. Scanning the symmetric structure and taking the MAXIMUM along r(C-Cl). Wrong direction:
     the D3h structure is a MINIMUM along the symmetric stretch and a maximum only along the
     asymmetric one, which is the reaction coordinate. The scan's "top" landed on the edge
     of the grid and its Hessian had two negative eigenvalues -- both signs that it was not
     a stationary point at all.

The fix uses symmetry properly. At a D3h geometry the gradient along every non-symmetric
direction vanishes by symmetry, so minimising the energy inside the two-parameter symmetric
subspace (r_CCl, r_CH) lands on a genuine stationary point of the full surface. Its Hessian
then has to show exactly one negative eigenvalue, and that is the check.

The reference is ONE calculation at total charge -1 with the ion 20 A away, so the constant
size-consistency offset of the global charge embedding cancels out of every reported number.
"""
import json
import sys
from pathlib import Path
import numpy as np
from ase import Atoms
from ase.optimize import BFGS
from scipy.optimize import minimize

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from prrs.config import SearchConfig                                 # noqa: E402
from prrs.runner import confirm_minimum, curvature_spectrum, _wavenumbers  # noqa: E402
from mace.calculators import MACECalculator                          # noqa: E402

MACE = Path("/home/ruigengji/MLP/mace")
MODELS = [("MACE-omol-0-4M", MACE / "mace-omol-0-extra-large-4M.model", "MACE"),
          ("MACE-POLAR-1-M", MACE / "MACE-POLAR-1-M.model", "PolarMACE")]
KCAL = 23.060548
CONFIG = SearchConfig(closed_shell_only=False, quench_fmax_eV_A=0.01,
                      minimum_check_step_A=0.01)
# Literature, gas phase, high level. NOT verified inside this session; quoted so the
# comparison is falsifiable and flagged so nobody reads it as measured here.
LIT = {"complex_kcal": -10.5, "barrier_kcal": 3.0}


def charged(atoms, q=-1, s=1):
    out = atoms.copy()
    out.info["charge"], out.info["spin"] = q, s
    return out


def d3h(r_ccl, r_ch):
    """Planar CH3 midway between two equal C-Cl bonds: the symmetric subspace."""
    positions = [[0, 0, 0], [0, 0, r_ccl], [0, 0, -r_ccl]]
    for k in range(3):
        phi = 2 * np.pi * k / 3
        positions.append([r_ch * np.cos(phi), r_ch * np.sin(phi), 0.0])
    return Atoms("CClClH3", positions=positions)


def run(name, path, model_type):
    calc = MACECalculator(model_paths=str(path), model_type=model_type,
                          device="cuda", default_dtype="float64")

    def factory():
        return calc

    def energy(atoms, q=-1, s=1):
        probe = charged(atoms, q, s)
        probe.calc = calc
        calc.results.clear()
        return float(probe.get_potential_energy())

    def relax(atoms, q=-1, fmax=0.01):
        probe = charged(atoms, q)
        probe.calc = calc
        BFGS(probe, logfile=None).run(fmax=fmax, steps=400)
        return probe

    print("=" * 78)
    print(f"{name}   model_type={model_type}")

    ch3cl = relax(Atoms("CClH3", positions=[
        [0, 0, 0], [0, 0, 1.785], [1.028, 0, -0.373],
        [-0.514, 0.890, -0.373], [-0.514, -0.890, -0.373]]), q=0)
    apart = ch3cl.copy()
    apart += Atoms("Cl", positions=[[0, 0, -20.0]])
    reference = energy(apart)
    print(f"   relaxed CH3Cl r(C-Cl) = {ch3cl.get_distance(0, 1):.4f} A")
    print(f"   reference (ion 20 A away, one calculation at charge -1) = {reference:.6f} eV")

    # ---- ion-dipole complex: full relaxation, backside attack ----
    start = ch3cl.copy()
    start += Atoms("Cl", positions=[[0, 0, -3.2]])
    complexed = relax(start)
    e_complex = complexed.get_potential_energy()
    settled, checked = confirm_minimum(complexed, factory, CONFIG, 0)
    print(f"   complex: r(C-Cl) = {complexed.get_distance(0, 1):.4f}  "
          f"r(C...Cl-) = {complexed.get_distance(0, 5):.4f} A")
    print(f"            {(e_complex - reference) * KCAL:+7.2f} kcal/mol vs reference   "
          f"(literature about {LIT['complex_kcal']:+.1f}, unverified here)")
    print(f"            confirm_minimum = {settled}  order = {checked.get('saddle_order')}")

    # ---- D3h stationary point: MINIMISE inside the symmetric subspace ----
    calls = {"n": 0}

    def symmetric_energy(parameters):
        calls["n"] += 1
        r_ccl, r_ch = float(parameters[0]), float(parameters[1])
        if not (1.9 < r_ccl < 3.2 and 0.95 < r_ch < 1.25):
            return 1e3
        return energy(d3h(r_ccl, r_ch))

    result = minimize(symmetric_energy, x0=[2.32, 1.070], method="Nelder-Mead",
                      options={"xatol": 1e-4, "fatol": 1e-7, "maxiter": 400})
    r_ccl, r_ch = float(result.x[0]), float(result.x[1])
    saddle = charged(d3h(r_ccl, r_ch))
    saddle.calc = calc
    e_saddle = saddle.get_potential_energy()
    print(f"   D3h stationary point ({calls['n']} energy calls):")
    print(f"            r(C-Cl) x2 = {r_ccl:.4f} A   r(C-H) = {r_ch:.4f} A")
    print(f"            {(e_saddle - reference) * KCAL:+7.2f} kcal/mol vs reference   "
          f"(literature about {LIT['barrier_kcal']:+.1f}, unverified here)")
    print(f"            {(e_saddle - e_complex) * KCAL:+7.2f} kcal/mol vs the complex")

    values, _, floor, provenance = curvature_spectrum(saddle, factory, CONFIG,
                                                      source=CONFIG.hessian_source)
    negative = [float(v) for v in values if v < -CONFIG.minimum_check_eigenvalue_tol]
    icm = [round(v, 1) for v in _wavenumbers(negative)]
    verdict = ("first-order saddle" if len(negative) == 1
               else f"NOT a first-order saddle: {len(negative)} negative eigenvalues")
    print(f"            Hessian ({provenance.get('source')}): {verdict}")
    print(f"            imaginary wavenumbers {icm} cm^-1")
    del calc
    return {"reference_eV": reference,
            "ch3cl_r_CCl_A": round(float(ch3cl.get_distance(0, 1)), 4),
            "complex": {"r_CCl_A": round(float(complexed.get_distance(0, 1)), 4),
                        "r_C_Clminus_A": round(float(complexed.get_distance(0, 5)), 4),
                        "relative_kcal": round((e_complex - reference) * KCAL, 2),
                        "confirmed_minimum": bool(settled),
                        "saddle_order": checked.get("saddle_order")},
            "saddle": {"r_CCl_A": round(r_ccl, 4), "r_CH_A": round(r_ch, 4),
                       "energy_calls": calls["n"],
                       "vs_reference_kcal": round((e_saddle - reference) * KCAL, 2),
                       "vs_complex_kcal": round((e_saddle - e_complex) * KCAL, 2),
                       "negative_eigenvalues": negative,
                       "imaginary_cm": icm, "verdict": verdict,
                       "curvature_source": provenance.get("source")}}


def main():
    report = {"literature_unverified_here": LIT,
              "reference": "one calculation at total charge -1, ion 20 A away",
              "models": {}}
    for name, path, model_type in MODELS:
        if not path.exists():
            print(f"{name}: missing at {path}")
            continue
        try:
            report["models"][name] = run(name, path, model_type)
        except Exception as exc:                                     # noqa: BLE001
            print(f"{name}: {type(exc).__name__}: {exc}"[:200])
            report["models"][name] = {"error": f"{type(exc).__name__}: {exc}"[:200]}
    Path("runs/sn2_probe.json").write_text(json.dumps(report, indent=1), encoding="utf-8")
    print("\nwrote runs/sn2_probe.json")


if __name__ == "__main__":
    main()
