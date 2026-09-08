"""Does the Fukui-equilibrated charge localise correctly when fragments separate?

MACE-POLAR-1 satisfies the global total-charge constraint through a learnable Fukui
equilibration. That is the same family as charge equilibration, and the family's classic
pathology is charge smearing at dissociation: infinitely separated fragments come out with
fractional charges, because the equilibration has no way to know the charge belongs on one
side. SN1 is literally heterolysis into R+ and X-, so if this fails, SN1 is not describable
however good the energies look elsewhere.

Four tests, in increasing severity:
  1. does the energy respond to total charge at all
  2. atomic electron affinities, where the answer is known -- weak evidence either way,
     since single atoms are far outside a molecular training set
  3. size consistency: one calculation with the fragments far apart, against the sum of
     two separate calculations
  4. CHARGE LOCALISATION: with the fragments 40 A apart and a total charge of -1, does the
     predicted atomic charge put -1 on the chloride and 0 on the neutral molecule
"""
import json
import sys
from pathlib import Path
import numpy as np
from ase import Atoms
from ase.io import read
from mace.calculators import MACECalculator

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

MACE = Path("/home/ruigengji/MLP/mace")
KCAL = 23.060548
MODELS = [("MACE-POLAR-1-S", MACE / "MACE-POLAR-1-S.model", "PolarMACE"),
          ("MACE-POLAR-1-M", MACE / "MACE-POLAR-1-M.model", "PolarMACE"),
          ("MACE-omol-0-4M", MACE / "mace-omol-0-extra-large-4M.model", "MACE")]
EA = {"Cl": (3.613, 2), "Br": (3.364, 2), "F": (3.401, 2)}


def put(calc, atoms, charge=0, spin=1):
    probe = atoms.copy()
    probe.info["charge"], probe.info["spin"] = charge, spin
    probe.calc = calc
    calc.results.clear()
    return probe


def charges_of(calc, probe):
    """Whatever the model reports as atomic charges, by whichever route it exposes them."""
    probe.get_potential_energy()
    for source, key in ((probe.arrays, "Qs"), (probe.arrays, "charges"),
                        (calc.results, "charges"), (calc.results, "Qs")):
        if key in source:
            return np.asarray(source[key], dtype=float)
    return None


def main():
    report = {}
    for name, path, model_type in MODELS:
        if not path.exists():
            continue
        print("=" * 78)
        print(f"{name}   ({path.stat().st_size/1e6:.0f} MB, model_type={model_type})")
        try:
            calc = MACECalculator(model_paths=str(path), model_type=model_type,
                                  device="cuda", default_dtype="float64")
        except Exception as exc:                                     # noqa: BLE001
            print(f"   load failed: {type(exc).__name__}: {exc}"[:160])
            continue
        entry = {"r_max_A": float(getattr(calc.models[0], "r_max", float("nan")))}
        print(f"   r_max = {entry['r_max_A']:.1f} A")

        # 1 + 2: charge response and electron affinities
        entry["electron_affinity_eV"] = {}
        for symbol, (experiment, spin_neutral) in EA.items():
            atom = Atoms(symbol, positions=[[0, 0, 0]])
            neutral = put(calc, atom, 0, spin_neutral).get_potential_energy()
            anion = put(calc, atom, -1, 1).get_potential_energy()
            value = neutral - anion
            entry["electron_affinity_eV"][symbol] = {
                "model": round(value, 3), "experiment": experiment,
                "deviation": round(value - experiment, 3)}
            print(f"   EA({symbol}) = {value:6.3f} eV   experiment {experiment:.3f}   "
                  f"deviation {value - experiment:+.3f}")

        # 3: size consistency
        ch3cl = Atoms("CClH3", positions=[[0, 0, 0], [0, 0, 1.785],
                                          [1.028, 0, -0.373], [-0.514, 0.890, -0.373],
                                          [-0.514, -0.890, -0.373]])
        e_neutral = put(calc, ch3cl, 0).get_potential_energy()
        e_anion = put(calc, Atoms("Cl", positions=[[0, 0, 0]]), -1).get_potential_energy()
        separate = e_neutral + e_anion
        entry["size_consistency_kcal"] = {}
        # POLAR-1 sums the long-range part in reciprocal space (graph_longrange.kspace),
        # so the k-vector count grows with the bounding box, not the atom count: a 40 A
        # separation asks for a 4.24 GiB allocation and dies on an 11 GB card. Each
        # distance is therefore attempted independently and a failure is recorded, not
        # fatal -- the box dependence is itself part of what this measures.
        for distance in (8.0, 12.0, 20.0, 30.0, 40.0):
            joined = ch3cl.copy()
            joined += Atoms("Cl", positions=[[0, 0, -distance]])
            try:
                together = put(calc, joined, -1).get_potential_energy()
            except Exception as exc:                                 # noqa: BLE001
                entry["size_consistency_kcal"][f"{distance:.0f} A"] = \
                    f"{type(exc).__name__}"
                print(f"   size consistency at {distance:5.1f} A: "
                      f"{type(exc).__name__} (reciprocal-space memory)")
                continue
            gap = (together - separate) * KCAL
            entry["size_consistency_kcal"][f"{distance:.0f} A"] = round(gap, 3)
            print(f"   size consistency at {distance:5.1f} A: {gap:+8.3f} kcal/mol")

        # 4: charge localisation -- the test that decides whether SN1 is describable
        # 20 A, not 40: far outside the 6 A cutoff so the fragments cannot see each
        # other, while the reciprocal-space box still fits in memory.
        joined = ch3cl.copy()
        joined += Atoms("Cl", positions=[[0, 0, -20.0]])
        probe = put(calc, joined, -1)
        try:
            q = charges_of(calc, probe)
        except Exception as exc:                                     # noqa: BLE001
            q = None
            print(f"   charge localisation: {type(exc).__name__}")
        if q is None:
            entry["charge_localisation"] = {"available": False,
                                            "note": "the model reports no atomic charges"}
            print("   charge localisation: the model reports no atomic charges")
        else:
            molecule = float(q[:5].sum())
            ion = float(q[5])
            entry["charge_localisation"] = {
                "available": True, "on_CH3Cl": round(molecule, 4),
                "on_chloride": round(ion, 4), "total": round(float(q.sum()), 4),
                "per_atom": [round(float(v), 4) for v in q]}
            print(f"   charge localisation at 20 A, total -1:")
            print(f"      on CH3Cl    {molecule:+.4f}")
            print(f"      on chloride {ion:+.4f}")
            print(f"      sum         {q.sum():+.4f}")
            verdict = ("localised" if abs(ion + 1.0) < 0.05 and abs(molecule) < 0.05
                       else "SMEARED -- heterolysis (SN1) is not describable this way")
            entry["charge_localisation"]["verdict"] = verdict
            print(f"      -> {verdict}")
        report[name] = entry
        del calc
    Path("runs/polar_charge_battery.json").write_text(json.dumps(report, indent=1),
                                                      encoding="utf-8")
    print("\nwrote runs/polar_charge_battery.json")


if __name__ == "__main__":
    main()
