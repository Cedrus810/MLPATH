"""One AFIR path, E/F only (docs/NONHESSIAN_AFIR_DESIGN_2026-11-14.md §3; plan
docs/PLAN_AFIR_MC_2026-11-14.md §2-3).

The AFIR path is the trajectory of a minimization of F = E + sign * alpha * rho_AB from a
start structure (Chem. Rec. 2016: "Such a path obtained by minimization of the AFIR
function is called AFIR path"). GRRM minimizes with RFO and an exact Hessian every 50
steps; here ASE BFGS on E/F only. Every step records the unbiased energy, rho_AB, the bond
graph and the geometry. Stops (GRRM17): a minimum of F, two passed energy maxima (off by
default since 2026-11-15, see PathSettings.stop_after_maxima), or dissociation; early discard after 10 steps if the force term has not decreased; plus a
model-domain guard (the force can push atoms where the MLP has no data) and step / budget
caps. Every evaluation is one unbiased forward on the caller's counted factory.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field

import numpy as np
from ase.data import covalent_radii
from ase.optimize import BFGS

from .afir import AFIRCalculator, alpha_eV_A
from .ledger import BudgetExhausted
from .poses import _jsonable


@dataclass(frozen=True)
class PathSettings:
    noise_eV: float  # energy resolution for maxima (v1.5 rule: 2 x float32 spacing)
    gamma_kJ: float = 106.9  # Chem. Rec.: T 298.15 K, t 10 days (user 2026-11-14)
    sign: int = +1
    maxstep_A: float = 0.1
    fmax_eV_A: float = 0.03
    max_steps: int = 1000
    early_discard_steps: int = 10
    dissociation_ratio: float = 3.0
    domain_min_ratio: float = 0.5
    domain_max_drop_eV: float = 5.0
    # GRRM17 stops a path after two passed energy maxima. With FIRE (inertial, no Hessian)
    # the "maxima" include overshoots: on N1 the frozen-coordinate paths were cut while
    # still climbing (a 0.013 eV overshoot, then a 0.0029 eV wiggle; DONE_LOG 2026-11-15).
    # 0 = off (user 2026-11-15 「先去掉“两个极大”停止」); 2 reproduces runs before that.
    stop_after_maxima: int = 0

    def __post_init__(self):
        if self.sign not in (+1, -1):
            raise ValueError("sign must be +1 or -1")
        if (
            isinstance(self.stop_after_maxima, bool)
            or not isinstance(self.stop_after_maxima, (int, np.integer))
            or self.stop_after_maxima < 0
        ):
            raise ValueError("stop_after_maxima must be a non-negative integer")
        for name, value in asdict(self).items():
            if name in ("sign", "stop_after_maxima"):
                continue
            if isinstance(value, bool) or not np.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be finite and positive")
        for name in ("max_steps", "early_discard_steps"):
            if not isinstance(getattr(self, name), (int, np.integer)):
                raise ValueError(f"{name} must be an integer")


@dataclass
class PathResult:
    status: str = "unresolved:not_started"
    frames: list = field(default_factory=list)
    maxima: list = field(default_factory=list)
    alpha_eV_A: float = 0.0
    provenance: dict = field(default_factory=dict)
    spent: int = 0

    def to_dict(self):
        return _jsonable(asdict(self))


class _Peaks:
    """Energy maxima along a trajectory, each confirmed by a rise and a fall > noise."""

    def __init__(self, noise):
        self.noise, self.maxima = noise, []
        self.trough = None
        self.peak = None  # (energy, index)

    def add(self, k, e):
        if self.trough is None:
            self.trough = e
            return
        if self.peak is None:
            if e > self.trough + self.noise:
                self.peak = (e, k)
            else:
                self.trough = min(self.trough, e)
            return
        if e >= self.peak[0]:
            self.peak = (e, k)
        elif e < self.peak[0] - self.noise:
            self.maxima.append(self.peak[1])
            self.peak, self.trough = None, e


def _ratios(numbers, positions):
    x = np.asarray(positions, float)
    d = np.linalg.norm(x[:, None] - x[None, :], axis=-1)
    rs = covalent_radii[numbers][:, None] + covalent_radii[numbers][None, :]
    iu = np.triu_indices(len(x), 1)
    return d[iu] / rs[iu], d, rs


def _dissociated(numbers, positions, edges, n_start, ratio):
    """More bond-graph components than at the start, and the closest atoms of two
    components farther apart than ``ratio`` x (R_i + R_j)."""
    n = len(numbers)
    parent = list(range(n))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i, j in edges:
        parent[find(i)] = find(j)
    label = [find(i) for i in range(n)]
    if len(set(label)) <= n_start:
        return False
    _, d, rs = _ratios(numbers, positions)
    cross = np.array(label)[:, None] != np.array(label)[None, :]
    return bool(np.min((d / rs)[cross]) > ratio)


def run_afir_path(atoms, coordinate, factory, config, settings, ledger, task, attempt=0):
    """Minimize F from ``atoms`` (positions = start structure); see module docstring."""
    from prrs.state import encode

    cfg = settings
    alpha = alpha_eV_A(cfg.gamma_kJ)
    result = PathResult(
        alpha_eV_A=alpha,
        provenance={
            "algorithm": "AFIR-path-EF-v1",
            "settings": asdict(cfg),
            "coordinate": {"a": list(coordinate.a), "b": list(coordinate.b), "p": coordinate.p},
        },
    )
    work = atoms.copy()
    work.set_constraint()
    numbers = work.numbers
    calc = AFIRCalculator(factory(), coordinate, alpha, cfg.sign)
    work.calc = calc
    peaks = _Peaks(cfg.noise_eV)
    start = ledger.spent
    n_start = None
    rho0 = None
    try:
        with ledger.segment(task=task, stage="afir_path", attempt=attempt) as seg:
            opt = BFGS(work, logfile=None, maxstep=cfg.maxstep_A)
            status = "max_steps"
            for k, converged in enumerate(opt.irun(fmax=cfg.fmax_eV_A, steps=cfg.max_steps)):
                ub = calc.unbiased
                graph = encode(work, config.bond_scale, active=config.active_atoms)
                edges = sorted([int(i), int(j)] for i, j in graph.edges)
                biased_f = float(np.linalg.norm(work.get_forces(), axis=1).max())
                frame = {
                    "step": k,
                    "energy_eV": ub["energy"],
                    "rho_A": ub["rho"],
                    "fmax_biased": biased_f,
                    "fmax_unbiased": float(np.linalg.norm(ub["forces"], axis=1).max()),
                    "edges": edges,
                    "positions": work.positions.copy(),
                }
                prev = result.frames[-1] if result.frames else None
                result.frames.append(frame)
                if n_start is None:
                    n_start, rho0 = graph.components, ub["rho"]
                peaks.add(k, ub["energy"])
                ratios, _, _ = _ratios(numbers, work.positions)
                if ratios.min() < cfg.domain_min_ratio or (
                    prev is not None
                    and prev["energy_eV"] - ub["energy"] > cfg.domain_max_drop_eV
                ):
                    status = "domain_exit"
                    break
                if converged:
                    status = "biased_minimum"
                    break
                if cfg.stop_after_maxima and len(peaks.maxima) >= cfg.stop_after_maxima:
                    status = "two_maxima"
                    break
                if _dissociated(
                    numbers, work.positions, edges, n_start, cfg.dissociation_ratio
                ):
                    status = "dissociated"
                    break
                if k == cfg.early_discard_steps and cfg.sign * (ub["rho"] - rho0) > 0:
                    status = "not_reactive"
                    break
            seg.outcome = status
            result.status = status
    except BudgetExhausted:
        result.status = "budget_exhausted"
    finally:
        result.maxima = list(peaks.maxima)
        result.spent = ledger.spent - start
    return result
