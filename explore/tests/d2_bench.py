"""TwoSaddle bench for the D2 tests and the kill->resume driver.

Pose: H at the left minimum (-1, 0). H-O1 compression and H-O0 stretch both drive H
over to the right minimum; O0-O1 stretch only loads the anchor frame (no trigger)."""

import os

from ase.calculators.calculator import all_changes

from explore.d2 import Backends, D2Settings, HFSettings, init_run
from explore.intents import IntentSpec
from explore.probe import Intent
from prrs.config import SearchConfig
from twosaddle import TwoSaddleInvariant, two_saddle_atoms

CONFIG = SearchConfig(closed_shell_only=False)
POSE = "p0"


def spec(pair, sign, pose=POSE, layer="near"):
    iid = f"{pose}|distance:{pair[0]}-{pair[1]}:{'+' if sign > 0 else '-'}"
    return IntentSpec(Intent(iid, (pair,), (1.0,)), sign, "distance", layer, "inter")


SPECS = (spec((1, 3), -1), spec((0, 3), +1), spec((0, 1), +1))


def settings(**over):
    base = dict(
        seed=7,
        search_budget=20000,
        certify_budget=20000,
        scan_cap=600,
        reverse_cap=200,
        release_cap=600,
        certify_cap=3000,
        per_pose_quota=10,
        audit_fraction=0.0,
        noise_eV=1e-6,
        certify_max_rise_eV=5.0,
        magnitudes_A=tuple(round(0.1 * k, 10) for k in range(1, 18)),
        fmax_eV_A=1e-3,
        max_steps_per_point=300,
        max_rise_eV=3.0,
        max_active=2,
        wave_size=4,
    )
    base.update(over)
    return D2Settings(**base)


def pose():
    atoms = two_saddle_atoms(-1.0, 0.0)
    atoms.info.update(charge=0, spin=1)
    return atoms


# HF-4 bench inputs, measured on this bench like the molecular ones: 8 re-relaxations
# (0.01 A perturbations) at a near-root point of p0|distance:0-3:+ (q 2.7594) to the
# midpoint threshold 1e-3 gave reaction-force sigma 6.7e-4 eV/A and max pairwise RMSD
# 3.3e-4 A (2026-10-26). The jitter is the start-dependent transverse residual, not the
# (float64 central-difference) force error.
HF = HFSettings(
    reaction_force_sigma=7e-4,
    geometry_noise_A=3.3e-4,
    midpoint_transverse_fmax_eV_A=1e-3,
    hs_cap=3000,
    hc_cap=3000,
)


def init(
    run_dir,
    specs=SPECS,
    hf=None,
    start_cluster=None,
    same_graph_ridge_last=False,
    hs_bracket=None,
    same_graph_sample=False,
    hf_parallel=False,
    **over,
):
    init_run(
        run_dir,
        [(POSE, pose())],
        {POSE: list(specs)},
        settings(**over),
        CONFIG,
        identity={"bench": "twosaddle"},
        hf=hf,
        start_cluster=start_cluster,
        same_graph_ridge_last=same_graph_ridge_last,
        hs_bracket=hs_bracket,
        same_graph_sample=same_graph_sample,
        hf_parallel=hf_parallel,
    )


class Tick:
    """External evaluation counter: one byte per model evaluation, written before it
    runs; evaluation number ``crash_after`` (in this process) kills it like SIGKILL."""

    def __init__(self, path=None, crash_after=-1):
        self.path, self.crash_after, self.count = path, crash_after, 0

    def __call__(self):
        self.count += 1
        if self.path is not None:
            with open(self.path, "ab") as fh:
                fh.write(b".")
        if self.count == self.crash_after:
            os._exit(137)


class CountedTwoSaddle(TwoSaddleInvariant):
    def __init__(self, tick):
        super().__init__()
        self.tick = tick

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        self.tick()
        super().calculate(atoms, properties, system_changes)


def backends(tick):
    def batch(frames):
        out = []
        for frame in frames:
            tick()
            atoms = frame.copy()
            atoms.calc = TwoSaddleInvariant()
            out.append((atoms.get_potential_energy(), atoms.get_forces()))
        return out

    return Backends(search_batch=batch, certify_calculator=CountedTwoSaddle(tick))
