"""Shared fixtures: a synthetic PRRS run directory and analytic potentials.

Tests never read the repository's runs/ and never import prrs: the synthetic
run builder writes the exact file contract documented in
pars/adapters/prrs_run.py, and the potentials are small analytic ASE
calculators defined here.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import write as ase_write

NUMBERS = [1, 6, 8]
PHASES = ("source", "perturbed", "free_start", "free", "free", "free", "quench", "quenched")


def _trajectory(base: np.ndarray, step: np.ndarray, frames: int = 8) -> np.ndarray:
    return base + step * np.arange(frames)[:, None, None]


def _trial_record(
    attempt_id: str, amplitude: float, classification: str, uncertainty_status: str
) -> dict:
    return {
        "attempt_id": attempt_id,
        "probe": {
            "family": "stretch",
            "indices": [0, 1],
            "sign": 1,
            "amplitude": amplitude,
            "seed": 1234 + int(attempt_id[1:]),
            "unit": "Angstrom",
            "kind": "bond",
            "mode": "displace",
        },
        "status": "completed",
        "classification": classification,
        "fingerprint": {"broken": [], "formed": [], "active_atoms": [0, 1]},
        "failure": None,
        "trajectory": f"trials/{attempt_id}/response.traj",
        "observations": f"trials/{attempt_id}/observations.jsonl",
        "endpoint": f"trials/{attempt_id}/endpoint.extxyz",
        "uncertainty_status": uncertainty_status,
        "max_uncertainty_eV_A": 0.05 if uncertainty_status == "active" else None,
    }


def _write_trial(
    run: Path,
    attempt_id: str,
    positions: np.ndarray,
    energies: np.ndarray,
    edges_per_frame: list,
    uncertainties,
    endpoint: np.ndarray,
    endpoint_energy: float,
    classification: str = "reactive",
) -> None:
    # PRRS stamps the declared electronic state on the structure it evolves,
    # so every response frame and the endpoint carry it (search.py)
    state = {"charge": 0, "spin": 1}
    trial = run / "trials" / attempt_id
    trial.mkdir(parents=True)
    numbers = np.asarray(NUMBERS)
    frames = []
    for index, frame_positions in enumerate(positions):
        atoms = Atoms(numbers=numbers, positions=frame_positions)
        atoms.calc = SinglePointCalculator(
            atoms, energy=float(energies[index]), forces=np.zeros((len(numbers), 3))
        )
        atoms.info = {"phase": PHASES[index], "time_fs": float(index), **state}
        frames.append(atoms)
    ase_write(str(trial / "response.traj"), frames)

    rows = []
    for index in range(len(positions)):
        rows.append(
            {
                "frame": index,
                "phase": PHASES[index],
                "step": index * 10,
                "time_fs": float(index),
                "potential_eV": float(energies[index]),
                "kinetic_eV": 0.0,
                "physical_total_eV": float(energies[index]),
                "max_force_eV_A": 0.1,
                "uncertainty_eV_A": (
                    None if uncertainties is None else float(uncertainties[index])
                ),
                "coordination": [1.0] * len(numbers),
                "edges": edges_per_frame[index],
                "events": {"formed": [], "broken": []},
            }
        )
    (trial / "observations.jsonl").write_text(
        "".join(json.dumps(row) + "\n" for row in rows), encoding="utf-8"
    )

    endpoint_atoms = Atoms(numbers=numbers, positions=endpoint, info=dict(state))
    endpoint_atoms.calc = SinglePointCalculator(
        endpoint_atoms, energy=endpoint_energy, forces=np.zeros((len(numbers), 3))
    )
    ase_write(str(trial / "endpoint.extxyz"), endpoint_atoms)

    record = _trial_record(
        attempt_id,
        0.1 if classification == "elastic" else 0.6,
        classification,
        "active" if uncertainties is not None else "unavailable",
    )
    (trial / "result.json").write_text(json.dumps(record, indent=2), encoding="utf-8")


def build_run(root: Path, with_uncertainty: bool = False, subdir: str = "seed1/pIa") -> Path:
    """Write a minimal but contract-faithful PRRS run directory.

    Two reactive trials in two channels (rc0000, rc0001) and one elastic trial
    that ingestion must skip. Trial energies differ so the channels have
    distinct cheap-objective landscapes.
    """
    run = root / subdir
    run.mkdir(parents=True, exist_ok=True)
    numbers = np.asarray(NUMBERS)
    base = np.zeros((len(numbers), 3))
    base[:, 0] = [0.0, 1.1, 2.2]

    records = [
        _trial_record("t000000", 0.1, "elastic", "unavailable"),
        _trial_record("t000001", 0.6, "reactive", "unavailable"),
        _trial_record(
            "t000002", 0.6, "reactive", "active" if with_uncertainty else "unavailable"
        ),
    ]
    (run / "attempts.jsonl").write_text(
        "".join(json.dumps(record) + "\n" for record in records), encoding="utf-8"
    )

    bonded = [[0, 1], [1, 2]]
    reacted = [[0, 2]]
    full_early = [bonded, bonded, bonded]
    full_late = [reacted, reacted, reacted, reacted, reacted]

    energies_low = np.asarray([-1.0, 0.2, 0.9, 0.4, -0.2, -0.8, -0.8, -0.8])
    positions_low = _trajectory(base, np.asarray([[[0.3, 0.0, 0.0]]]))
    _write_trial(
        run,
        "t000001",
        positions_low,
        energies_low,
        full_early + full_late,
        None,
        base + np.asarray([[2.5, 0.0, 0.0]]),
        -0.9,
    )

    energies_high = np.asarray([-1.0, 0.4, 1.7, 1.2, 0.6, 0.1, 0.1, 0.1])
    positions_high = _trajectory(base, np.asarray([[[0.0, 0.3, 0.0]]]))
    uncertainties = np.linspace(0.01, 0.06, 8) if with_uncertainty else None
    _write_trial(
        run,
        "t000002",
        positions_high,
        energies_high,
        full_early + full_late,
        uncertainties,
        base + np.asarray([[0.0, 2.5, 0.0]]),
        0.0,
    )

    # elastic trial: complete files so a bug cannot hide behind missing paths
    _write_trial(
        run,
        "t000000",
        positions_low * 0.01,
        energies_low * 0.01,
        [bonded] * 8,
        None,
        base,
        -1.0,
        classification="elastic",
    )

    network = {
        "schema_version": 2,
        "reaction_channels": [
            {
                "id": "rc0000",
                "event_key": "aaaabbbbccccdddd",
                "class": "reaction",
                "ends": ["c0000", "c0001"],
                "self_loop": False,
                "broken": [[0, 1]],
                "formed": [[0, 2]],
                "symmetry_enumerated": True,
                "undirected": True,
                "observed_by": ["t000001"],
                "ts_support": [],
            },
            {
                "id": "rc0001",
                "event_key": "eeeeffffgggghhhh",
                "class": "reaction",
                "ends": ["c0000", "c0002"],
                "self_loop": False,
                "broken": [[1, 2]],
                "formed": [[0, 1]],
                "symmetry_enumerated": True,
                "undirected": True,
                "observed_by": ["t000002"],
                "ts_support": [],
            },
        ],
        "reactions": [
            {
                "id": "r0000",
                "source": "c0000",
                "target": "c0001",
                "observed_only": True,
                "attempts": ["t000001"],
                "classifications": ["reactive"],
                "best_by_family": {
                    "stretch": {"amplitude": 0.6, "unit": "Angstrom", "attempt": "t000001"}
                },
            },
            {
                "id": "r0001",
                "source": "c0000",
                "target": "c0002",
                "observed_only": True,
                "attempts": ["t000002"],
                "classifications": ["reactive"],
                "best_by_family": {
                    "stretch": {"amplitude": 0.6, "unit": "Angstrom", "attempt": "t000002"}
                },
            },
        ],
    }
    (run / "network.json").write_text(json.dumps(network, indent=2), encoding="utf-8")

    manifest = {
        "schema_version": 2,
        "units": {
            "energy": "eV",
            "length": "Angstrom",
            "force": "eV/Angstrom",
            "time": "fs",
        },
        "config": {"charge_sensitive": False, "total_charge": None, "multiplicity": None},
        "calculator": {"factory": "synthetic"},
        "implementation_sha256": "0" * 64,
        "versions": {},
    }
    (run / "config.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return run


# -- analytic potentials --------------------------------------------------------


class TripleWellCalculator(Calculator):
    """Single-atom 2D triple well with analytic gradients.

    E(r) = k/2 |r|^2 - sum_i D_i exp(-|r_xy - c_i|^2 / (2 s^2))

    The wells are narrow enough (s=0.4 vs d=1.5) that the three minima near
    c0=(0,0), c1=(d,0), c2=(0,d) are genuinely separated by real saddles. The
    c0->c1 saddle sits at about (0.774, 0.0005): c2 tilts it slightly off the
    x axis, so y=0 is close to, but not exactly, the MEP.
    """

    implemented_properties = ["energy", "forces"]

    def __init__(self, d=1.5, depth=3.0, width=0.4, stiffness=0.5, **kwargs):
        super().__init__(**kwargs)
        self.centers = np.asarray([[0.0, 0.0], [d, 0.0], [0.0, d]])
        self.depth = depth
        self.width = width
        self.stiffness = stiffness

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        positions = self.atoms.positions
        energy = 0.0
        forces = np.zeros_like(positions)
        for center in self.centers:
            delta = positions[:, :2] - center
            sq = np.sum(delta**2, axis=1)
            gaussian = self.depth * np.exp(-sq / (2.0 * self.width**2))
            energy -= float(gaussian.sum())
            forces[:, :2] -= (gaussian / self.width**2)[:, None] * delta
        energy += 0.5 * self.stiffness * float(np.sum(positions**2))
        forces -= self.stiffness * positions
        self.results["energy"] = energy
        self.results["forces"] = forces


class NarrowBarrierCalculator(Calculator):
    """Single atom, smooth double well along x with a narrow 0.4 eV barrier.

    E = h exp(-x^2 / 2w^2) + c (x^2 - a^2)^2 + (y^2 + z^2) / 2

    Minima at x = +-1.0342, saddle at x = 0, barrier 0.41947 eV (dense 1D
    scan); barrier curvature h/w^2 ~ 5 eV/A^2, about a proton transfer's. The
    basins are flat enough that weak NEB springs leave images where the seed
    put them.
    """

    implemented_properties = ["energy", "forces"]

    def __init__(self, height=0.4, width=0.28, wall=0.02, a=1.0, **kwargs):
        super().__init__(**kwargs)
        self.height, self.width, self.wall, self.a = height, width, wall, a

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        positions = self.atoms.positions
        x = positions[:, 0]
        gaussian = self.height * np.exp(-(x**2) / (2.0 * self.width**2))
        energy = np.sum(gaussian + self.wall * (x**2 - self.a**2) ** 2)
        energy += 0.5 * np.sum(positions[:, 1:] ** 2)
        forces = np.zeros_like(positions)
        forces[:, 0] = gaussian * x / self.width**2 - 4.0 * self.wall * x * (x**2 - self.a**2)
        forces[:, 1:] = -positions[:, 1:]
        self.results["energy"] = float(energy)
        self.results["forces"] = forces


class ProtonHopCalculator(Calculator):
    """Three atoms A(0) B(1) H(2), a function of interatomic distances only, so
    the energy is invariant under rigid motion like a gas-phase molecule.

    E = V(q) + kab/2 (r_AB - d0)^2 + kb (r_AH + r_BH - r_AB)^2,  q = r_AH - r_BH
    V(q) = h exp(-q^2 / 2w^2) + c (q^2 - a^2)^2   (NarrowBarrierCalculator's x profile)

    With A-B at d0 and H on the A-B line the extra terms vanish, so minima are at
    q = +-1.0342 and the saddle at q = 0, 0.41947 eV above them. Forces are
    central finite differences of the energy (step 1e-6 A).
    """

    implemented_properties = ["energy", "forces"]
    height, width, wall, a, d0, k_ab, k_bend = 0.4, 0.28, 0.02, 1.0, 2.5, 5.0, 2.0

    def _energy(self, p: np.ndarray) -> float:
        r_ab = np.linalg.norm(p[0] - p[1])
        r_ah = np.linalg.norm(p[0] - p[2])
        r_bh = np.linalg.norm(p[1] - p[2])
        q = r_ah - r_bh
        return float(
            self.height * np.exp(-(q**2) / (2.0 * self.width**2))
            + self.wall * (q**2 - self.a**2) ** 2
            + 0.5 * self.k_ab * (r_ab - self.d0) ** 2
            + self.k_bend * (r_ah + r_bh - r_ab) ** 2
        )

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        p = self.atoms.positions.copy()
        forces = np.zeros_like(p)
        step = 1e-6
        for i in range(len(p)):
            for j in range(3):
                p[i, j] += step
                up = self._energy(p)
                p[i, j] -= 2.0 * step
                down = self._energy(p)
                p[i, j] += step
                forces[i, j] = -(up - down) / (2.0 * step)
        self.results["energy"] = self._energy(p)
        self.results["forces"] = forces


class UncertaintyCalculator(Calculator):
    """Wrapper that adds a constant force_uncertainty_eV_A to a base calculator,
    mirroring the results convention of PRRS's GuardedCalculator."""

    implemented_properties = ["energy", "forces"]

    def __init__(self, base: Calculator, uncertainty: float = 0.5, **kwargs):
        super().__init__(**kwargs)
        self.base = base
        self.uncertainty_value = float(uncertainty)

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        self.base.calculate(atoms, properties, system_changes)
        self.results["energy"] = self.base.results["energy"]
        self.results["forces"] = self.base.results["forces"]
        self.results["force_uncertainty_eV_A"] = self.uncertainty_value


class ExplodingCalculator(Calculator):
    """Raises on every call, to exercise backend failure reporting."""

    implemented_properties = ["energy", "forces"]

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        raise RuntimeError("synthetic calculator failure")


# -- fixtures --------------------------------------------------------------------


@pytest.fixture
def run_dir(tmp_path) -> Path:
    return build_run(tmp_path, with_uncertainty=False)


@pytest.fixture
def run_dir_with_uncertainty(tmp_path) -> Path:
    return build_run(tmp_path, with_uncertainty=True, subdir="seed1/pIa-unc")


def make_candidate(
    candidate_id: str,
    channel_id: str,
    energies,
    signatures,
    positions=None,
    uncertainty=None,
    amplitude=0.5,
):
    """Test helper: a PathCandidate without touching a run directory."""
    from pars.schema import PathCandidate, Provenance

    n_atoms = 3
    if positions is not None:
        n_atoms = int(np.asarray(positions).shape[1])
    if positions is None:
        positions = np.linspace(0.0, 1.0, len(energies))[:, None, None] * np.ones(
            (1, n_atoms, 3)
        )
    positions = np.asarray(positions, dtype=float)
    return PathCandidate(
        candidate_id=candidate_id,
        channel_id=channel_id,
        trace_id=f"trace-{candidate_id}",
        numbers=np.asarray(NUMBERS[:n_atoms]),
        positions=positions,
        energies=np.asarray(energies, dtype=float),
        uncertainty=uncertainty,
        signatures=tuple(frozenset(map(tuple, sig)) for sig in signatures),
        provenance=Provenance(
            run_dir="synthetic",
            attempt_id=f"t-{candidate_id}",
            channel_id=channel_id,
            event_key=None,
            reactant_state_id="cA",
            product_state_id="cB",
            frame_start=0,
            frame_stop=len(energies),
            perturbation_family="stretch",
            perturbation_amplitude=amplitude,
            perturbation_unit="Angstrom",
        ),
    )


TRIO_EARLY = frozenset({(0, 1), (1, 2)})
TRIO_LATE = frozenset({(0, 2)})
