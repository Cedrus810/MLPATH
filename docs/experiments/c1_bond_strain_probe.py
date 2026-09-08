#!/usr/bin/env python3
"""What the chord probe costs in bond strain, on real saddles.

`response_curvature` measures curvature along a direction by displacing to R +/- a v in
straight Cartesian lines. That is the right map for a curvature -- a chord has no
second-order term, so its a -> 0 limit is v^T H v exactly -- but it is NOT the arc the
search itself walks, and at the top of the default amplitude band a probe on hydrogen
moves about 0.1 Angstrom, twice `internal`'s own single-step cap.

`kappa` does not care. `anharmonicity` and `sign_stable` are read off the finite-amplitude
samples, and a sample whose bonds have been stretched is not reporting on the motion it is
named after. `bond_strain_by_amplitude` is that bill. It costs no evaluation to measure
and until now nobody had read it on a real surface.

Run it on the saddles the ts job left on disk:

    PYTHONPATH=src python docs/experiments/c1_bond_strain_probe.py

Reports, does not gate.
"""
import glob
import json
import os
import sys
from pathlib import Path

import numpy as np
from ase.io import read

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "src"))

MODEL = os.environ.get("MODEL", str(ROOT / "models" / "MACE-OFF24_medium.model"))
DEVICE = os.environ.get("DEVICE", "cuda")
PATTERN = os.environ.get("SADDLES", str(ROOT / "runs" / "p3_ts_barrier" / "ts_k*.extxyz"))

from prrs.config import SearchConfig                                  # noqa: E402
from prrs.runner import confirm_minimum                               # noqa: E402


def main():
    paths = sorted(glob.glob(PATTERN))
    if not paths:
        print(f"no saddle geometries matched {PATTERN}")
        return 2

    from prrs.calculators import mace_factory
    factory = lambda: mace_factory(MODEL, device=DEVICE, default_dtype="float64")
    config = SearchConfig(closed_shell_only=False, quench_fmax_eV_A=0.002)

    print(f"{len(paths)} saddle geometries, band "
          f"{config.response_band_A[0]}..{config.response_band_A[1]} A*sqrt(amu)\n")
    rows = []
    for path in paths:
        atoms = read(path)
        atoms.calc = factory()
        confirmed, checked = confirm_minimum(atoms, factory, config, 0)
        response = checked.get("response")
        name = Path(path).name
        if not isinstance(response, dict) or "bond_strain_by_amplitude" not in response:
            reason = (response or {}).get("reason") if isinstance(response, dict) else None
            print(f"{name}: no response report "
                  f"(order={checked.get('saddle_order')}, reason={reason})")
            continue
        strain = response["bond_strain_by_amplitude"]
        amplitudes = [a for a, _ in response["samples"]]
        rows.append({"file": name, "saddle_order": checked.get("saddle_order"),
                     "kappa": response["kappa"], "sigma": response["sigma"],
                     "anharmonicity": response["anharmonicity"],
                     "sign_stable": response["sign_stable"],
                     "retraction": response["retraction"],
                     "amplitudes": amplitudes,
                     "bond_strain_by_amplitude": strain,
                     "worst_bond_strain": response["worst_bond_strain"]})
        print(f"{name}  order={checked.get('saddle_order')}  "
              f"kappa={response['kappa']:+.5g}  sign_stable={response['sign_stable']}")
        for amplitude, value in zip(amplitudes, strain):
            bar = "#" * int(round(value * 200))
            print(f"    a={amplitude:.4f}  bond strain {value:7.2%}  {bar}")

    if not rows:
        print("\nno response reports collected")
        return 1

    worst = max(r["worst_bond_strain"] for r in rows)
    top = [r["bond_strain_by_amplitude"][-1] for r in rows]
    print(f"\nworst strain anywhere: {worst:.2%}")
    print(f"at the top of the band: {min(top):.2%} to {max(top):.2%} "
          f"(median {float(np.median(top)):.2%})")
    print("\nThe chord is the correct map for kappa_0 regardless of this number. What it "
          "bounds is\nhow far the finite-amplitude samples have travelled from the motion "
          "they are named\nafter, which is what anharmonicity and sign_stable are read "
          "off.")
    out = ROOT / "runs" / "c1_bond_strain.json"
    out.write_text(json.dumps({"model": MODEL, "rows": rows}, indent=1))
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
