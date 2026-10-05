"""A narrow intermediate well skipped by the default finite negative-mode kick."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.calculators.calculator import Calculator, all_changes

from pars.config import ParsConfig
from pars.iox import write_json_atomic
from pars.provenance import implementation_sha256, numerical_versions
from pars.refinement.ase_neb import AseNebBackend
from pars.schema import PathCandidate, Provenance


def energy(x):
    return (x * x - 1) ** 2 - 0.02 * np.exp(-(((x - 0.015) / 0.002) ** 2))


def force(x):
    u = (x - 0.015) / 0.002
    return -4 * x * (x * x - 1) - 20 * u * np.exp(-u * u)


class NarrowIntermediate(Calculator):
    implemented_properties = ["energy", "forces"]

    def calculate(
        self, atoms=None, properties=("energy", "forces"), system_changes=all_changes
    ):
        super().calculate(atoms, properties, system_changes)
        x, y, z = self.atoms.positions[0]
        self.results = {
            "energy": float(energy(x) + 0.5 * (y * y + z * z)),
            "forces": np.array([[force(x), -y, -z]]),
        }


def root(left, right):
    assert force(left) * force(right) < 0
    for _ in range(60):
        middle = (left + right) / 2
        if force(left) * force(middle) > 0:
            left = middle
        else:
            right = middle
    return (left + right) / 2


def run(output):
    output = Path(output)
    if output.exists():
        raise FileExistsError(output)
    positions = np.zeros((9, 1, 3))
    positions[:, 0, 0] = np.linspace(-1, 1, 9)
    provenance = Provenance(
        "analytic-narrow-well",
        "straight-chain",
        "outer-to-outer",
        None,
        "x=-1",
        "x=+1",
        0,
        9,
        extra={
            "endpoint_basin_references": {
                "reactant": positions[0].tolist(),
                "product": positions[-1].tolist(),
            }
        },
    )
    candidate = PathCandidate(
        "narrow",
        "outer-to-outer",
        "straight-chain",
        np.array([1]),
        positions,
        np.array([energy(x) for x in positions[:, 0, 0]]),
        None,
        None,
        provenance,
    )
    config = ParsConfig(refinement_remove_rigid_motion=False)
    result = AseNebBackend(NarrowIntermediate()).refine(candidate, config)
    minimum = root(0.015, 0.016)
    step = 1e-6
    curvature = (
        energy(minimum + step) - 2 * energy(minimum) + energy(minimum - step)
    ) / step**2
    assert curvature > 0
    assert all(force(x) > 0 for x in np.linspace(0.0001, minimum - 1e-6, 10000))
    report = {
        "evidence_scope": "analytic_external_field_resolution_counterexample",
        "potential": "E=(x^2-1)^2+(y^2+z^2)/2-0.02*exp(-((x-0.015)/0.002)^2)",
        "reference": {
            "first_right_minimum_x": minimum,
            "next_saddle_x": root(0.017, 0.03),
            "minimum_energy_difference_curvature_eV_A2": curvature,
            "main_saddle": [0, 0, 0],
            "coarse_kick_A": 0.05,
        },
        "config": config.to_dict(),
        "versions": numerical_versions(),
        "runner_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "pars_implementation_sha256": implementation_sha256(),
        "result": result.to_summary(),
        "band_positions_A": None if result.positions is None else result.positions.tolist(),
        "band_energies_eV": None if result.energies is None else result.energies.tolist(),
    }
    write_json_atomic(output, report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True)
    report = run(parser.parse_args().output)
    print(json.dumps({k: report["result"][k] for k in ("status", "barrier_eV", "message")}))
