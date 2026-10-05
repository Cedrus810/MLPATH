"""Does the fragment-rotation coordinate move b05 off the -76 cm^-1 stationary point?

Idea B of the 2026-09-13 brainstorm, measured on the case that produced it. b05's three
elimination endpoints ended the continuation run (`runs/b05_quench_continuation/`) on
order-1 stationary points at -75 to -77 cm^-1 whose unstable mode is, measured rather
than argued, 99.4% on the water and 99.92% a RIGID ROTATION of it
(`runs/b05_saddle_descent/summary.json`, `within_water`).

The reading recorded at the time was that this is a domain gap, not a budget: the soft
mode polish iterates `rotatable_torsions`, which enumerates bonds of the graph, and the
nascent water shares no bond with the isobutene, so its libration was a coordinate nothing
in the polish could reach. `fragment_rotations` supplies it. This checks the claim the way
the project requires -- by measuring, on the same three structures, with everything else
held fixed:

    same potential, same config, same fmax, same round budget. The only difference is
    whether `polish_soft_modes` is given the fragment rotations.

What would falsify it: the saddle order does not change, or it changes but the polish
moved the structure without lowering the energy. Both are printed.

IT WAS FALSIFIED, and the second half of this script says why. The coordinate is now in
the domain -- six are offered and one Newton step is taken, worth -0.43 meV on all three
-- but the saddle order stays 1 and the imaginary frequency stays near -76 cm^-1, because
`polish_soft_modes` deliberately leaves NEGATIVE curvature alone: not being a minimum
along a coordinate is the minimum check's verdict to reach, not the polish's to erase.
So the second half asks a different question, the one the answer turns on: a fragment
rotation is exact AND periodic, so a full-period scan is well defined wherever the local
curvature is unusable -- the same argument the `flat_biased` tier already makes. Scanning
finds a water orientation 28.5-29.1 meV lower, at a quarter turn, on all three samples.
And the structure is STILL order 1, now at -54 cm^-1. So b05's blocker is not a polish
domain problem after all; the residual mode is a coupled inter-fragment motion that one
coordinate at a time cannot remove, and fragment TRANSLATIONS are the half deliberately
left out here. Whether the polish should scan off negative curvature for an exact periodic
coordinate is a semantic change to what `confirm_minimum` is handed, so it is measured
here and decided elsewhere.

    python docs/experiments/b05_fragment_rotation_polish.py
"""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

from prrs.config import SearchConfig  # noqa: E402
from prrs.perturbations import fragment_rotations, step_fragment_rotation  # noqa: E402
from prrs.reliability import GuardedCalculator  # noqa: E402
from prrs.runner import confirm_minimum, polish_soft_modes  # noqa: E402
from prrs.state import encode  # noqa: E402

ENDPOINTS = ROOT / "runs" / "b05_quench_continuation"


def main(argv):
    from ase.io import read

    from benchmarks.b04_sn2_chloride.calculator import polar1

    # The grid's config carries // comments, as b05_quench_continuation.py already has to
    # strip. Same stripping, same file, so this runs against the same knobs that run did.
    text = (ROOT / "benchmarks/b05_tertiary_halide/config.json").read_text()
    raw = json.loads(
        "\n".join(line for line in text.splitlines() if not line.lstrip().startswith("//"))
    )
    config = SearchConfig(**raw)
    physical = polar1()

    results = {}
    for path in sorted(ENDPOINTS.glob("*_endpoint.extxyz")):
        name = path.stem.replace("_endpoint", "")
        atoms = read(str(path))
        atoms.info["charge"], atoms.info["spin"] = config.total_charge, config.multiplicity
        graph = encode(atoms, config.bond_scale, active=config.active_atoms)
        modes = fragment_rotations(atoms, graph)

        def checked(structure):
            probe = structure.copy()
            probe.calc = GuardedCalculator(physical, config, None)
            probe.get_potential_energy()  # curvature_spectrum reads the gradient for its floor
            return confirm_minimum(probe, polar1, config, seed=0)

        before_ok, before = checked(atoms)
        guard = GuardedCalculator(physical, config, None)
        working = atoms.copy()
        working.calc = guard
        start_energy = float(working.get_potential_energy())

        reports = []
        for _ in range(config.soft_polish_rounds):
            report = polish_soft_modes(working, guard, config)
            reports.append(report)
            if report["converged"] and not report["moved"]:
                break
        end_energy = float(working.get_potential_energy())
        after_ok, after = checked(working)

        rotated = [
            movement
            for report in reports
            for movement in report["movements"]
            if len(movement["indices"]) != 4  # a fragment, not a torsion
        ]
        results[name] = {
            "fragment_rotation_coordinates": len(modes),
            "rounds": len(reports),
            "saddle_order_before": before.get("saddle_order"),
            "saddle_order_after": after.get("saddle_order"),
            "icm_before": before.get("imaginary_wavenumbers_icm"),
            "icm_after": after.get("imaginary_wavenumbers_icm"),
            "energy_change_eV": end_energy - start_energy,
            "fragment_rotations_taken": len(rotated),
            "largest_rotation_rad": max((abs(m["offset_rad"]) for m in rotated), default=0.0),
            "confirmed_minimum_before": bool(before_ok),
            "confirmed_minimum_after": bool(after_ok),
        }
        entry = results[name]
        print(f"{name}")
        print(f"  fragment rotation coordinates offered : {entry['fragment_rotation_coordinates']}")
        print(f"  polish rounds used                    : {entry['rounds']}")
        print(
            f"  saddle order   before -> after        : "
            f"{entry['saddle_order_before']} -> {entry['saddle_order_after']}"
        )
        print(
            f"  imaginary cm^-1                       : "
            f"{np.round(entry['icm_before'] or [], 1).tolist()} -> "
            f"{np.round(entry['icm_after'] or [], 1).tolist()}"
        )
        print(f"  energy change from the polish         : {entry['energy_change_eV'] * 1000:+.3f} meV")
        print(
            f"  fragment rotations actually taken     : {entry['fragment_rotations_taken']}"
            f"  (largest {entry['largest_rotation_rad']:.4f} rad)"
        )
        print(
            f"  confirmed minimum                     : "
            f"{entry['confirmed_minimum_before']} -> {entry['confirmed_minimum_after']}"
        )

        # Second question: would a full-period SCAN do what the Newton step cannot?
        scanned = atoms.copy()
        scanned.calc = GuardedCalculator(physical, config, None)
        scan_start = float(scanned.get_potential_energy())
        turned = 0.0
        for mode in modes:
            offsets = np.linspace(-np.pi, np.pi, config.soft_mode_scan_points * 8 + 1)[:-1]
            energies = [
                float(
                    scanned.calc.get_potential_energy(
                        step_fragment_rotation(scanned, mode, float(d), in_place=False)[0]
                    )
                )
                for d in offsets
            ]
            best = int(np.argmin(energies))
            if energies[best] < float(scanned.calc.get_potential_energy(scanned)) - 1e-6:
                step_fragment_rotation(scanned, mode, float(offsets[best]))
                turned += abs(float(offsets[best]))
        scan_end = float(scanned.calc.get_potential_energy(scanned))
        scan_ok, scan_check = checked(scanned)
        entry.update(
            scan_gain_meV=1000.0 * (scan_start - scan_end),
            scan_rotation_rad=turned,
            scan_saddle_order=scan_check.get("saddle_order"),
            scan_icm=scan_check.get("imaginary_wavenumbers_icm"),
            scan_confirmed_minimum=bool(scan_ok),
        )
        print(
            f"  full-period scan instead              : {entry['scan_gain_meV']:+.1f} meV "
            f"over {entry['scan_rotation_rad']:.3f} rad -> order "
            f"{entry['scan_saddle_order']}, "
            f"{np.round(entry['scan_icm'] or [], 1).tolist()} cm^-1, "
            f"minimum={entry['scan_confirmed_minimum']}"
        )
        print()

    out = ENDPOINTS / "fragment_rotation_polish.json"
    out.write_text(json.dumps(results, indent=2) + "\n")
    print(f"written {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
