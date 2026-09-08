"""Where does P2's B-minus-A energy come from?

Two numbers are in the record and they disagree by a factor of six:

    preflight (p2_preflight2.json)   dE(B-A) = -11.72 meV
    round 3   (network.json)         dE(c0001-c0000) = -1.75 meV

Both nodes carry the graph hashes of the preflight A and B, so the disagreement is not a
mislabelled state. It is either a different conformer of B or a looser stationary point.
The torsion fingerprint decides which, and it decides it without the calculator.
"""
import json
import sys
from pathlib import Path
import numpy as np
from ase.io import read

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from prrs import internal
from prrs.config import SearchConfig                                 # noqa: E402
from prrs.network import torsion_profile                             # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
cfg = SearchConfig.from_dict(json.loads((ROOT / "runs/p2_round3/config.json").read_text())["config"])
net = json.loads((ROOT / "runs/p2_round3/network.json").read_text())
pre = json.loads((ROOT / "runs/p2_oxobutanal/p2_preflight2.json").read_text())["report"]

nodes = {n["id"]: n for n in net["chemical_nodes"]}
run_A = nodes["c0000"]["microstates"][0]
run_B = nodes["c0001"]["microstates"][0]

print("=== 1. the arithmetic ===")
dE_run = run_B["energy_eV"] - run_A["energy_eV"]
dE_pre = pre["B"]["E"] - pre["A"]["E"]
print(f"  run   A {run_A['energy_eV']:.9f}   B {run_B['energy_eV']:.9f}   dE {dE_run*1e3:+.3f} meV")
print(f"  pre   A {pre['A']['E']:.9f}   B {pre['B']['E']:.9f}   dE {dE_pre*1e3:+.3f} meV")
print(f"  A side  run - pre = {(run_A['energy_eV']-pre['A']['E'])*1e3:+.6f} meV")
print(f"  B side  run - pre = {(run_B['energy_eV']-pre['B']['E'])*1e3:+.6f} meV")
print(f"  the whole discrepancy sits on B: {(dE_run-dE_pre)*1e3:+.3f} meV")

print("\n=== 2. the same torsions, measured on the preflight structures ===")
paths = {"A": ROOT / "runs/p2_oxobutanal/p2_input.extxyz",
         "B": ROOT / "runs/p2_oxobutanal/p2_product_reference.extxyz"}
measured = {}
for label, path in paths.items():
    atoms = read(str(path))
    measured[label] = torsion_profile(atoms, cfg.bond_scale, cfg.active_atoms)
    print(f"  {label} ({path.name})")
    for indices, value in measured[label].items():
        print(f"      {list(indices)!s:>14}  {value:+.6f} rad")

def compare(recorded, reference, label):
    print(f"\n  {label}: recorded (run) vs measured (preflight)")
    worst = 0.0
    for key, value in recorded.items():
        indices = tuple(json.loads(key))
        if indices not in reference:
            print(f"      {key:>14}  {value:+.4f}   NOT A GENUINE TORSION of the preflight structure")
            continue
        delta = abs((value - reference[indices] + np.pi) % (2 * np.pi) - np.pi)
        worst = max(worst, delta)
        verdict = "same" if delta < 0.05 else "DIFFERENT"
        print(f"      {key:>14}  {value:+.4f} vs {reference[indices]:+.4f}   |d|={delta:.4f} rad  {verdict}")
    missing = [list(k) for k in reference if k not in {tuple(json.loads(s)) for s in recorded}]
    if missing:
        print(f"      recorded set is missing {missing}")
    print(f"      worst circular difference: {worst:.4f} rad")
    return worst

worst_A = compare(run_A["torsions_rad"], measured["A"], "c0000/m0000 vs preflight A")
worst_B = compare(run_B["torsions_rad"], measured["B"], "c0001/m0000 vs preflight B")

print("\n=== 3. verdict ===")
if worst_B < 0.05:
    print("  same conformer -> the 9.97 meV is a stationary-point tightness difference,")
    print("  not a different structure. The descent quench stopped at a looser point than")
    print("  the preflight polish did.")
else:
    print("  DIFFERENT conformer -> the run's B is another microstate of the same chemical")
    print("  state and the two numbers are not comparable as 'the' reaction energy.")
print(f"  free_bonds  run B {run_B['free_bonds']}   run A {run_A['free_bonds']}")
print(f"  quench fmax in force during the run: {cfg.quench_fmax_eV_A} eV/A")
print(f"  preflight fmax actually reached:  A {pre['A']['fmax']:.2e}  B {pre['B']['fmax']:.2e} eV/A")

print("\n=== 4. how much energy can fmax = 0.02 eV/A hide? ===")
print("  Displace the preflight B minimum along each projected Hessian mode until the")
print("  largest atomic force reaches the run's tolerance, and read off the energy rise.")
print("  This is the run's own energy resolution, measured rather than assumed.")
from prrs.calculators import load_factory                            # noqa: E402
from prrs.reliability import GuardedCalculator                       # noqa: E402
from prrs.runner import hessian_spectrum                             # noqa: E402

meta = json.loads((ROOT / "runs/p2_round3/config.json").read_text())["calculator"]
factory = load_factory(meta["factory"], meta["kwargs"])
calculator = factory()
B = read(str(paths["B"]))
guard = GuardedCalculator(calculator, cfg)
B.calc = guard
E0 = B.get_potential_energy()
f0 = float(np.linalg.norm(B.get_forces(), axis=1).max())
print(f"  B single point: E = {E0:.9f} eV   fmax = {f0:.3e} eV/A")
print(f"  matches the preflight record to {abs(E0 - pre['B']['E'])*1e3:.2e} meV")

values, vectors = hessian_spectrum(B, lambda: calculator, cfg, return_vectors=True)
masses = B.get_masses()
positions = B.positions.copy()

def rise_at_tolerance(mode_index):
    """Energy above the minimum at the displacement where fmax hits the tolerance."""
    step = vectors[:, mode_index].reshape(-1, 3) / np.sqrt(masses)[:, None]
    step = step / np.linalg.norm(step, axis=1).max()
    def probe(scale):
        B.set_positions(positions + scale * step)
        energy = B.get_potential_energy()
        return energy - E0, float(np.linalg.norm(B.get_forces(), axis=1).max())
    low, high = 0.0, 0.01
    for _ in range(30):
        _, fmax = probe(high)
        if fmax >= cfg.quench_fmax_eV_A or high > 2.0:
            break
        low, high = high, high * 2
    for _ in range(30):
        middle = (low + high) / 2
        _, fmax = probe(middle)
        low, high = (low, middle) if fmax > cfg.quench_fmax_eV_A else (middle, high)
    rise, fmax = probe((low + high) / 2)
    return (low + high) / 2, rise, fmax

nontrivial = [i for i, v in enumerate(values) if abs(v) > 1e-8]
picks = nontrivial[:4] + nontrivial[-1:]
print(f"  {'mode':>5} {'eigenvalue':>13} {'displacement':>13} {'energy rise':>13}")
for index in picks:
    scale, rise, fmax = rise_at_tolerance(index)
    print(f"  {index:>5} {values[index]:>13.5f} {scale:>10.4f} A {rise*1e3:>10.3f} meV"
          f"   (fmax {fmax:.4f})")
B.set_positions(positions)
print(f"\n  reference: the discrepancy to explain is {(dE_run-dE_pre)*1e3:+.3f} meV,")
print(f"  and the reaction energy itself is {dE_pre*1e3:+.3f} meV.")

print("\n=== 5. the run's own B, measured directly ===")
print("  The round-3 run directory was recovered, so the 9.97 meV can be measured on the")
print("  geometry that carries it instead of being attributed by elimination.")
from prrs.runner import confirm_minimum, quench                      # noqa: E402
from prrs.chemistry import free_aligned_symmetric_rmsd               # noqa: E402
from prrs.network import free_bonds_of                               # noqa: E402

runB = read(str(ROOT / "runs/p2_round3/chemical/c0001_m0000.extxyz"))
runB.calc = GuardedCalculator(calculator, cfg)
E_run = runB.get_potential_energy()
f_run = float(np.linalg.norm(runB.get_forces(), axis=1).max())
print(f"  c0001/m0000 single point: E = {E_run:.9f} eV   fmax = {f_run:.3e} eV/A")
print(f"  network.json says {run_B['energy_eV']:.9f}  ->  {abs(E_run-run_B['energy_eV'])*1e3:.2e} meV apart")
print(f"  fmax is {'inside' if f_run <= cfg.quench_fmax_eV_A else 'OUTSIDE'} the run tolerance"
      f" {cfg.quench_fmax_eV_A}")
print(f"  above the preflight B minimum by {(E_run - E0)*1e3:+.3f} meV")

reference = read(str(paths["B"]))
plain = float(np.sqrt(np.mean(np.sum((runB.positions - reference.positions) ** 2, axis=1))))
print(f"  raw coordinate RMSD to preflight B: {plain:.4f} A  (no alignment; includes"
      " translation/rotation)")

print("\n  re-quench the run's B with the run's own settings, then tighten:")
tight = runB.copy()
guard2 = GuardedCalculator(calculator, cfg)
tight.calc = guard2
ok, diag = quench(tight, guard2, cfg)
print(f"    run settings (fmax {cfg.quench_fmax_eV_A}):  converged={ok}  "
      f"E={tight.get_potential_energy():.9f}  "
      f"fmax={float(np.linalg.norm(tight.get_forces(), axis=1).max()):.3e}")
from dataclasses import replace as dc_replace                        # noqa: E402
strict = dc_replace(cfg, quench_fmax_eV_A=1e-3)
tight2 = runB.copy()
guard3 = GuardedCalculator(calculator, strict)
tight2.calc = guard3
ok2, diag2 = quench(tight2, guard3, strict)
E_tight = tight2.get_potential_energy()
print(f"    fmax 1e-3:                    converged={ok2}  E={E_tight:.9f}  "
      f"fmax={float(np.linalg.norm(tight2.get_forces(), axis=1).max()):.3e}")
print(f"    tightening recovered {(E_run - E_tight)*1e3:+.3f} meV")
print(f"    and lands {(E_tight - E0)*1e3:+.3f} meV from the preflight B minimum")
print(f"    torsions after tightening: "
      f"{ {str(list(k)): round(v,4) for k,v in torsion_profile(tight2, cfg.bond_scale, cfg.active_atoms).items()} }")

print("\n=== 6. is the run's B a different minimum, or an unfinished one? ===")
print("  Decisive: give it enough steps and a tight tolerance. Converging onto the")
print("  preflight energy means it was unfinished; converging elsewhere means it is a")
print("  different stationary point and the two energies were never comparable.")
print(f"  elements: {list(zip(range(len(runB)), runB.get_chemical_symbols()))}")
long_run = dc_replace(cfg, quench_fmax_eV_A=1e-4, quench_steps=4000, soft_polish_rounds=6)
far = runB.copy()
guard4 = GuardedCalculator(calculator, long_run)
far.calc = guard4
ok3, diag3 = quench(far, guard4, long_run)
E_far = far.get_potential_energy()
print(f"  fmax 1e-4, 4000 steps: converged={ok3} steps={diag3.get('steps')} "
      f"reason={diag3.get('reason')}")
print(f"  E = {E_far:.9f}   fmax = {float(np.linalg.norm(far.get_forces(), axis=1).max()):.3e}")
print(f"  vs preflight B: {(E_far - E0)*1e3:+.4f} meV")
print(f"  torsions: { {str(list(k)): round(v,4) for k,v in torsion_profile(far, cfg.bond_scale, cfg.active_atoms).items()} }")
ok4, checked4 = confirm_minimum(far, lambda: calculator, long_run, 0)
print(f"  confirm_minimum: {ok4}  order={checked4.get('saddle_order')} "
      f"icm={checked4.get('imaginary_wavenumbers_icm')}")
print("\n  the same question asked of the methyl rotor, which the torsion fingerprint")
print("  excludes by design (a 120 degree turn is a symmetry copy):")
for tag, atoms_ in (("run B", runB), ("tightened", far), ("preflight B", reference)):
    ang = internal.coordinate(atoms_.positions, "dihedral", (2, 3, 8, 9)) if len(atoms_) > 9 else None
    print(f"    {tag:>12}: dihedral(2,3,8,9) = {ang}")
