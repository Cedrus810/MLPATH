"""Frozen two-saddle reference and same-budget selection/refinement experiment.

E=(x²+y²-1)²+a*y²+b*y³+z_k*z²/2 for a mobile H, with two
harmonically anchored O observers. Contacts are perceived from distances,
never supplied as saddle labels. Two index-one saddles connect identical
mapped minima. This is an external-field analytic stress test, not chemistry.
"""

import argparse
import copy
import dataclasses
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.data import covalent_radii

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "pars" / "src"))
from pars.config import ParsConfig  # noqa: E402
from pars.frontier.archive import ParetoArchive  # noqa: E402
from pars.iox import write_json_atomic  # noqa: E402
from pars.provenance import implementation_sha256, numerical_versions  # noqa: E402
from pars.refinement.ase_neb import AseNebBackend  # noqa: E402
from pars.schema import PathCandidate, Provenance  # noqa: E402

PROTOCOL = Path(__file__).with_name("protocol.json")
ANCHORS = np.array([[0.0, 1.8, 0.0], [0.0, -1.8, 0.0]])


class ForceBudgetExceeded(RuntimeError):
    pass


class TwoSaddles(Calculator):
    implemented_properties = ["energy", "forces"]

    def __init__(self, potential, limit=None):
        super().__init__()
        self.potential = potential
        self.limit = limit
        self.force_evaluations = 0

    def calculate(
        self, atoms=None, properties=("energy", "forces"), system_changes=all_changes
    ):
        if self.limit is not None and self.force_evaluations >= self.limit:
            raise ForceBudgetExceeded("common force-evaluation budget exhausted")
        super().calculate(atoms, properties, system_changes)
        self.force_evaluations += 1  # every calculate evaluates both E and F
        a, b = self.potential["a"], self.potential["b"]
        x, y, z = self.atoms.positions[0]
        r = x * x + y * y - 1
        offsets = self.atoms.positions[1:] - ANCHORS
        energy = (
            r * r
            + a * y * y
            + b * y**3
            + 0.5 * self.potential["z_k"] * z * z
            + 0.5 * self.potential["anchor_k"] * np.sum(offsets**2)
        )
        gradient = np.array(
            [4 * x * r, 4 * y * r + 2 * a * y + 3 * b * y * y, self.potential["z_k"] * z]
        )
        forces = np.vstack([-gradient, -self.potential["anchor_k"] * offsets])
        self.results = {"energy": float(energy), "forces": forces}


def geometry(x, y=0.0):
    return np.vstack([[x, y, 0.0], ANCHORS])


def atoms_at(positions, calculator):
    atoms = Atoms("HOO", positions=positions)
    atoms.calc = calculator
    return atoms


def hessian(positions, potential, step=1e-4):
    base = np.asarray(positions)
    columns = []
    for j in range(base.size):
        plus, minus = base.copy(), base.copy()
        plus.flat[j] += step
        minus.flat[j] -= step
        columns.append(
            -(
                atoms_at(plus, TwoSaddles(potential)).get_forces().ravel()
                - atoms_at(minus, TwoSaddles(potential)).get_forces().ravel()
            )
            / (2 * step)
        )
    matrix = np.asarray(columns).T
    return np.linalg.eigh((matrix + matrix.T) * 0.5)


def references(protocol):
    p, limits = protocol["potential"], protocol["reference"]
    roots = np.roots([4.0, 3 * p["b"], 2 * p["a"] - 4.0])
    result = []
    for y in sorted(float(r) for r in roots):
        positions = geometry(0.0, y)
        atoms = atoms_at(positions, TwoSaddles(p))
        eig, modes = hessian(positions, p)
        fmax = float(np.linalg.norm(atoms.get_forces(), axis=1).max())
        negatives = int(np.sum(eig < limits["negative_eigenvalue_cutoff"]))
        if fmax > limits["stationary_fmax"] or negatives != 1:
            raise ValueError("reference is not an independently checked first-order saddle")
        endpoints = []
        descent_cost = 0
        for sign in (-1, 1):
            start = positions + sign * limits["descent_displacement_A"] * modes[:, 0].reshape(
                3, 3
            )
            landed = atoms_at(start, TwoSaddles(p))
            # Positive-step gradient descent preserves the downhill direction;
            # quasi-Newton minimizers can instead return to the saddle from a
            # small negative-mode displacement. This is an independent descent,
            # not a force-only stationary-point acceptance.
            for _ in range(10000):
                forces = landed.get_forces()
                if np.linalg.norm(forces, axis=1).max() < 1e-7:
                    break
                energy_before = landed.get_potential_energy()
                landed.positions += 0.02 * forces
                if landed.get_potential_energy() > energy_before + 1e-12:
                    raise ValueError("reference descent increased energy")
            else:
                raise ValueError("reference descent failed")
            nearest = min((-1, 1), key=lambda x: np.linalg.norm(landed.positions - geometry(x)))
            rmsd = float(
                np.sqrt(np.mean(np.sum((landed.positions - geometry(nearest)) ** 2, axis=1)))
            )
            endpoint_eigenvalues, _ = hessian(landed.positions, p)
            if rmsd > limits["endpoint_match_rmsd_A"] or min(endpoint_eigenvalues) <= 0:
                raise ValueError("reference does not connect the expected confirmed minimum")
            endpoints.append(
                {
                    "minimum": nearest,
                    "rmsd_A": rmsd,
                    "eigenvalues": endpoint_eigenvalues.tolist(),
                }
            )
            descent_cost += landed.calc.force_evaluations
        if {e["minimum"] for e in endpoints} != {-1, 1}:
            raise ValueError("reference saddle does not connect both minima")
        result.append(
            {
                "id": f"saddle-{len(result)}",
                "positions_A": positions.tolist(),
                "energy_eV": atoms.get_potential_energy(),
                "fmax": fmax,
                "eigenvalues": eig.tolist(),
                "descents": endpoints,
                "force_evaluations": atoms.calc.force_evaluations
                + 6 * positions.size
                + descent_cost,
            }
        )
    return result


def contacts(positions, numbers):
    return frozenset(
        (i, j)
        for i in range(len(numbers))
        for j in range(i + 1, len(numbers))
        if np.linalg.norm(positions[i] - positions[j])
        <= 1.2 * (covalent_radii[numbers[i]] + covalent_radii[numbers[j]])
    )


def population(protocol):
    candidates = []
    numbers = np.array([1, 8, 8])
    for sign in (-1, 1):
        for amplitude in protocol["pool"]["amplitudes"]:
            theta = np.linspace(np.pi, 0.0, protocol["pool"]["frames"])
            positions = np.asarray(
                [geometry(np.cos(t), sign * amplitude * np.sin(t)) for t in theta]
            )
            cid = "pc-" + hashlib.sha256(positions.tobytes()).hexdigest()[:12]
            energies = np.array(
                [
                    atoms_at(frame, TwoSaddles(protocol["potential"])).get_potential_energy()
                    for frame in positions
                ]
            )
            provenance = Provenance(
                run_dir="analytic_seed_pool",
                channel_id="same-mapped-minima",
                event_key="analytic_external_field",
                reactant_state_id="x=-1",
                product_state_id="x=+1",
                attempt_id=cid,
                frame_start=0,
                frame_stop=len(positions),
                perturbation_amplitude=amplitude,
                extra={
                    "endpoint_basin_references": {
                        "reactant": geometry(-1).tolist(),
                        "product": geometry(1).tolist(),
                    }
                },
            )
            candidates.append(
                PathCandidate(
                    cid,
                    "same-mapped-minima",
                    cid,
                    numbers,
                    positions,
                    energies,
                    None,
                    tuple(contacts(frame, numbers) for frame in positions),
                    provenance,
                )
            )
    return candidates


def saddle_label(refined, reference, protocol):
    if refined.positions is None or refined.status.value != "ci_neb_converged":
        return None
    top = refined.positions[int(np.argmax(refined.energies))]
    matched = [
        r["id"]
        for r in reference
        if np.sqrt(np.mean(np.sum((top - np.array(r["positions_A"])) ** 2, axis=1)))
        <= protocol["reference"]["saddle_match_rmsd_A"]
    ]
    return matched[0] if len(matched) == 1 else None


def run(out, protocol_path=PROTOCOL):
    out = Path(out)
    if out.exists():
        raise FileExistsError(f"will not reuse {out}")
    protocol = json.loads(Path(protocol_path).read_text())
    out.mkdir(parents=True)
    config = ParsConfig(**protocol["refinement"])
    started = time.perf_counter()
    selection_start = time.perf_counter()
    candidates = population(protocol)
    archive = ParetoArchive("same-mapped-minima", candidates, config)
    legacy = ParetoArchive(
        "same-mapped-minima",
        copy.deepcopy(candidates),
        dataclasses.replace(config, preserve_mechanism_diversity=False),
    )
    selection_wall = time.perf_counter() - selection_start
    k = protocol["selection"]["max_candidates"]
    ids = sorted(c.candidate_id for c in candidates)
    plans = {
        "pars": archive.refinement_order()[:k],
        "channel_pruning": legacy.refinement_order()[:k],
        "cheap_energy": sorted(
            ids, key=lambda cid: (archive.objectives_raw[cid]["energy_excursion"], cid)
        )[:k],
    }
    for seed in protocol["selection"]["random_seeds"]:
        plans[f"random-{seed:02d}"] = (
            np.random.default_rng(seed).choice(ids, size=k, replace=False).tolist()
        )
    # Freeze plans before measuring ANY refined outcome or reference labels.
    manifest = {
        "protocol": protocol,
        "protocol_sha256": hashlib.sha256(Path(protocol_path).read_bytes()).hexdigest(),
        "runner_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "pars_implementation_sha256": implementation_sha256(),
        "effective_refinement_config": config.to_dict(),
        "versions": numerical_versions(),
        "selection_wall_s": selection_wall,
        "pool_generation_force_evaluations": sum(len(c.positions) for c in candidates),
        "candidate_inputs_sha256": hashlib.sha256(
            b"".join(c.positions.tobytes() for c in candidates)
        ).hexdigest(),
        "plans": plans,
        "selection": archive.selection_dict(),
    }
    write_json_atomic(out / "manifest.json", manifest)
    reference = references(protocol)
    write_json_atomic(out / "reference.json", reference)
    oracle, oracle_cost = {}, 0
    for candidate in candidates:
        calc = TwoSaddles(protocol["potential"])
        refined = AseNebBackend(calc).refine(candidate, config)
        oracle[candidate.candidate_id] = {
            "status": refined.status.value,
            "saddle": saddle_label(refined, reference, protocol),
            "barrier_eV": refined.barrier_eV,
            "force_evaluations": calc.force_evaluations,
            "saddle_check": refined.provenance.get("saddle_check"),
            "connectivity_check": refined.provenance.get("connectivity_check"),
        }
        oracle_cost += calc.force_evaluations
    write_json_atomic(out / "oracle.json", oracle)
    by_id = {c.candidate_id: c for c in candidates}
    rows = []
    for name, plan in plans.items():
        calc = TwoSaddles(protocol["potential"], protocol["cost"]["force_evaluations_per_arm"])
        arm_start = time.perf_counter()
        outcomes = []
        for cid in plan:
            refined = AseNebBackend(calc).refine(by_id[cid], config)
            outcomes.append(
                {
                    "candidate": cid,
                    "status": refined.status.value,
                    "saddle": saddle_label(refined, reference, protocol),
                    "barrier_eV": refined.barrier_eV,
                    "connectivity_status": refined.provenance.get("connectivity_check", {}).get(
                        "status"
                    ),
                }
            )
        found = {r["saddle"] for r in outcomes if r["saddle"] is not None}
        rows.append(
            {
                "arm": name,
                "outcomes": outcomes,
                "coverage": len(found) / len(reference),
                "force_evaluations": calc.force_evaluations,
                "refinement_wall_s": time.perf_counter() - arm_start,
                "failed": sum(r["saddle"] is None for r in outcomes),
            }
        )
        write_json_atomic(out / "arms.json", rows)
        print(
            f"{name}: coverage={len(found)}/{len(reference)}, forces={calc.force_evaluations}",
            flush=True,
        )
    random_rows = [r for r in rows if r["arm"].startswith("random-")]
    pars_row = next(r for r in rows if r["arm"] == "pars")
    pool_labels = {v["saddle"] for v in oracle.values() if v["saddle"] is not None}
    accepted = (
        len(reference) == protocol["acceptance"]["reference_first_order_saddles"]
        and len(pool_labels) == len(reference)
        and pars_row["coverage"] == 1.0
    )
    report = {
        "status": "passed" if accepted else "failed",
        "evidence_scope": protocol["evidence_scope"],
        "reference_saddles": len(reference),
        "pool_saddles": len(pool_labels),
        "candidate_count": len(candidates),
        "reference_oracle_force_evaluations": oracle_cost,
        "reference_certification_force_evaluations": sum(
            r["force_evaluations"] for r in reference
        ),
        "shared_pool_generation_force_evaluations": sum(len(c.positions) for c in candidates),
        "cost_accounting": "per-arm cold refinement calls include endpoint checks and failures; shared pool generation, reference certification and all-pool oracle are reported separately; total wall includes preparation, all arms and oracle",
        "rows": rows,
        "random_mean_coverage": float(np.mean([r["coverage"] for r in random_rows])),
        "random_full_coverage_fraction": float(
            np.mean([r["coverage"] == 1 for r in random_rows])
        ),
        "total_wall_s": time.perf_counter() - started,
        "no_general_superiority_claim": True,
    }
    write_json_atomic(out / "report.json", report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = run(args.output)
    print(
        json.dumps({k: result[k] for k in ("status", "random_mean_coverage", "total_wall_s")})
    )
    raise SystemExit(0 if result["status"] == "passed" else 2)
