"""Three-minimum counterexample: a local index-one saddle is not an A->B witness."""

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


class ThreeMinima(Calculator):
    implemented_properties = ["energy", "forces"]

    def calculate(
        self, atoms=None, properties=("energy", "forces"), system_changes=all_changes
    ):
        super().calculate(atoms, properties, system_changes)
        x, y, z = self.atoms.positions[0]
        self.results = {
            "energy": float(x * x * (x * x - 1) ** 2 + 0.5 * (y * y + z * z)),
            "forces": np.array([[-2 * x * (x * x - 1) * (3 * x * x - 1), -y, -z]]),
        }


def run(output):
    output = Path(output)
    if output.exists():
        raise FileExistsError(output)
    config = ParsConfig(
        refinement_remove_rigid_motion=False, refinement_images=11, refinement_fmax_neb=0.001
    )
    x = np.linspace(-1, 1, 11)
    positions = np.zeros((11, 1, 3))
    positions[:, 0, 0] = x
    provenance = Provenance(
        run_dir="analytic-three-minima",
        attempt_id="straight-chain",
        channel_id="outer-to-outer",
        event_key=None,
        reactant_state_id="x=-1",
        product_state_id="x=+1",
        frame_start=0,
        frame_stop=11,
        extra={
            "endpoint_basin_references": {
                "reactant": positions[0].tolist(),
                "product": positions[-1].tolist(),
            }
        },
    )
    candidate = PathCandidate(
        "straight-chain",
        "outer-to-outer",
        "straight-chain",
        np.array([1]),
        positions,
        x * x * (x * x - 1) ** 2,
        None,
        None,
        provenance,
    )
    result = AseNebBackend(ThreeMinima()).refine(candidate, config)
    report = {
        "evidence_scope": "analytic_external_field_counterexample",
        "potential": "E=x^2*(x^2-1)^2+(y^2+z^2)/2",
        "reference": {
            "minima_x": [-1, 0, 1],
            "saddles_x": [-1 / np.sqrt(3), 1 / np.sqrt(3)],
            "saddle_eigenvalues_eV_A2": [-8 / 3, 1, 1],
            "local_barrier_eV": 4 / 27,
            "connectivity": [[-1, 0], [0, 1]],
        },
        "config": config.to_dict(),
        "pars_implementation_sha256": implementation_sha256(),
        "runner_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "versions": numerical_versions(),
        "result": result.to_summary(),
        "band_positions_A": None if result.positions is None else result.positions.tolist(),
        "band_energies_eV": None if result.energies is None else result.energies.tolist(),
    }
    write_json_atomic(output, report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    report = run(args.output)
    print(
        json.dumps({key: report["result"][key] for key in ("status", "barrier_eV", "message")})
    )
