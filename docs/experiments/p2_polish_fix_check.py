"""Does the free/flat_biased fix actually recover the 9.97 meV, on the real geometry?

The audit located the discrepancy in one methyl torsion that polish_soft_modes had tiered
"free" on a curvature of 0.0041 eV/rad^2 while its gradient was -0.014 eV/rad. This runs the
same quench on the same structure with the fix in place and reports where it lands. It also
checks the source structure A, whose energy must not move: if it does, the P2 comparison is
against a different reactant than rounds one to three used.
"""
import json
import sys
from pathlib import Path
import numpy as np
from ase.io import read

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from prrs.calculators import load_factory                             # noqa: E402
from prrs.config import SearchConfig                                  # noqa: E402
from prrs.reliability import GuardedCalculator                        # noqa: E402
from prrs.runner import confirm_minimum, quench, relax_source         # noqa: E402
from prrs.network import free_bonds_of                                # noqa: E402
from prrs import internal                                             # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
manifest = json.loads((ROOT / "runs/p2_round3/config.json").read_text())
cfg = SearchConfig.from_dict({**manifest["config"],
                              **{k: v for k, v in SearchConfig().to_dict().items()
                                 if k not in manifest["config"]}})
factory = load_factory(manifest["calculator"]["factory"], manifest["calculator"]["kwargs"])
calculator = factory()
pre = json.loads((ROOT / "runs/p2_oxobutanal/p2_preflight2.json").read_text())["report"]

print(f"soft_mode_scan_points = {cfg.soft_mode_scan_points}, "
      f"floor = {cfg.soft_mode_curvature_floor_eV_rad2} eV/rad^2, "
      f"fmax = {cfg.quench_fmax_eV_A} eV/A")

def polish_and_report(label, path, expected):
    atoms = read(str(path))
    guard = GuardedCalculator(calculator, cfg)
    atoms.calc = guard
    start = atoms.get_potential_energy()
    ok, diag = quench(atoms, guard, cfg)
    energy = atoms.get_potential_energy()
    fmax = float(np.linalg.norm(atoms.get_forces(), axis=1).max())
    print(f"\n{label}")
    print(f"  start   E = {start:.9f}")
    print(f"  quench  converged={ok} steps={diag.get('steps')} reason={diag.get('reason')}")
    print(f"  end     E = {energy:.9f}   fmax = {fmax:.3e}")
    print(f"  vs the tight reference {expected:.9f}: {(energy - expected)*1e3:+.4f} meV")
    print(f"  methyl dihedral(2,3,8,9) = "
          f"{internal.coordinate(atoms.positions, 'dihedral', (2, 3, 8, 9)):+.4f} rad")
    for round_ in diag.get("rounds", []):
        for mode in round_["modes"]:
            if mode["tier"] in ("flat_biased", "free"):
                print(f"    {mode['indices']} tier={mode['tier']} "
                      f"k={mode['curvature_eV_rad2']:.6f} "
                      f"g={mode['initial_gradient_eV_rad']:+.6f} "
                      f"bound={mode.get('free_gradient_bound_eV_rad')} "
                      f"scan_offset={mode.get('scan_offset_rad')} "
                      f"gain_meV={None if mode.get('scan_energy_gain_eV') is None else round(mode['scan_energy_gain_eV']*1e3, 3)}")
    atoms.info["quench"] = diag
    print(f"  free_bonds now: {free_bonds_of(atoms)}")
    return atoms, energy

runB, E_B = polish_and_report("round 3's product microstate c0001/m0000",
                              ROOT / "runs/p2_round3/chemical/c0001_m0000.extxyz",
                              pre["B"]["E"])
ok, checked = confirm_minimum(runB, lambda: calculator, cfg, 0)
print(f"  confirm_minimum: {ok} order={checked.get('saddle_order')}")

srcA, E_A = polish_and_report("the source structure A (must not move)",
                              ROOT / "runs/p2_oxobutanal/p2_input.extxyz",
                              pre["A"]["E"])

print("\n=== the number the benchmark reports ===")
print(f"  dE(B-A) with the fix   = {(E_B - E_A)*1e3:+.3f} meV")
print(f"  dE(B-A) preflight      = {(pre['B']['E'] - pre['A']['E'])*1e3:+.3f} meV")
print(f"  dE(B-A) round 3        = -1.750 meV")
