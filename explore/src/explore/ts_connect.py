"""E/F-only two-sided connection tracing from a located TS candidate (HF record F7, §6).

docs/PLAN_HF_CONNECT_2026-10-26.md K1. Replaces, for HF results only, the
Hessian-dependent chain confirm_minimum -> descend_saddle (unstable-mode start,
soft-mode polish, Hessian endpoint check, relays). Nothing here reads a model Hessian.

1. Starts come from the search itself: on each side of the candidate, the nearest
   converged constrained-branch frame whose reaction force is resolved,
   |dE/dq| >= resolve_factor * sigma. Only if a side has none, a displacement along the
   search tangent ``direction_cartesian`` is used, with the smallest listed amplitude
   whose energy drop exceeds three float32 energy quanta.
2. Each start is released (no constraint) and minimized with BFGS on E/F to the force
   tolerance; no torsion polish, since the outcome read off is the basin / bond graph.
3. Endpoint qualification is a finite-perturbation re-descent: every perturbed copy must
   come back with the same bond graph and without a resolved energy drop. In a convex
   basin E(x + d) >= E(x) - |F(x)|_2 |d|_2, so a drop larger than three float32 quanta
   plus |F_end|_2 times the aligned displacement norm means the re-descent left the
   endpoint's basin (the endpoint itself is converged only to the descent tolerance). Leaving a saddle drops the energy (b04 / b05
   2026-10-26: 0.4-2.4 eV); a floppy minimum may return to a neighbouring conformer at
   the same energy, which is the same chemical endpoint (b05 product complex: 0.21 A,
   +0.002 eV). ``same_basin`` is recorded alongside but does not decide. This is
   E/F-level stability evidence only, never a certified minimum. A perturbation whose
   energy rise is not resolved (<= 3 quanta) does not count as a test.
4. Both endpoints stable and in different basins -> ``connected`` (with classify's label
   and fingerprint); same basin -> ``same_basin``; anything else stays explicit.

A full Hessian of candidate or endpoints is a separate, optional property check and is
never consulted here.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
import time

import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator
from ase.optimize import BFGS

from .ledger import BudgetExhausted
from .poses import _jsonable
from .ts_ef import _rmsd


@dataclass(frozen=True)
class ConnectSettings:
    # Measured same-branch reaction-force jitter (1 sigma) of the system and model.
    reaction_force_sigma: float
    resolve_factor: float = 3.0
    # None: SearchConfig.quench_fmax_eV_A, the force tolerance of the old descent.
    descent_fmax_eV_A: float | None = None
    max_descent_steps: int = 1000
    max_step_A: float = 0.1
    tangent_amplitudes_A: tuple = (0.02, 0.05, 0.1, 0.2)
    perturb_A: float = 0.03
    n_perturb: int = 3
    seed: int = 0

    def __post_init__(self):
        for name in ("reaction_force_sigma", "resolve_factor", "max_step_A", "perturb_A"):
            value = getattr(self, name)
            if not np.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be finite and positive")
        if self.descent_fmax_eV_A is not None and not self.descent_fmax_eV_A > 0:
            raise ValueError("descent_fmax_eV_A must be positive")
        if not self.tangent_amplitudes_A or any(a <= 0 for a in self.tangent_amplitudes_A):
            raise ValueError("tangent_amplitudes_A must be positive")
        for name in ("max_descent_steps", "n_perturb", "seed"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (int, np.integer)) or value < 0:
                raise ValueError(f"{name} must be a non-negative integer")
        if self.max_descent_steps == 0:
            raise ValueError("max_descent_steps must be positive")


@dataclass
class ConnectionResult:
    status: str = "unresolved:not_started"
    sides: list = field(default_factory=list)
    label: str | None = None
    fingerprint: dict | None = None
    provenance: dict = field(default_factory=dict)
    spent: int = 0
    wall_s: float = 0.0

    def to_dict(self):
        return _jsonable(asdict(self))


class _Fail(RuntimeError):
    pass


def _quanta(energy):
    """Three float32 energy quanta at this magnitude: the smallest resolvable drop."""
    return 3 * float(np.spacing(np.float32(abs(energy))))


def _frozen(template, positions, energy, forces):
    atoms = template.copy()
    atoms.set_constraint()
    atoms.positions = np.array(positions, float)
    atoms.calc = SinglePointCalculator(atoms, energy=energy, forces=forces)
    return atoms


def branch_starts(
    frames, q_star, sigma, resolve_factor=3.0, transverse_fmax=None, descent_fmax=None
):
    """Nearest resolved converged branch frame on each side of q_star.

    ``frames`` are search-history dicts (q_actual, reaction_force, positions,
    fmax_projected, fmax_eV_A, ...). Low side: q < q_star with dE/dq > +limit (force
    points to lower q); high side: q > q_star with dE/dq < -limit. With ``descent_fmax``
    a start must not already meet the descent's own force tolerance -- such a start is a
    zero-step "descent" that stays next to the TS (b04 / b05 2026-10-26)."""
    limit = resolve_factor * sigma
    ok = [
        f
        for f in frames
        if (transverse_fmax is None or f.get("fmax_projected", np.inf) <= transverse_fmax)
        and (descent_fmax is None or f.get("fmax_eV_A", np.inf) > descent_fmax)
    ]
    low = [f for f in ok if f["q_actual"] < q_star and f["reaction_force"] > limit]
    high = [f for f in ok if f["q_actual"] > q_star and f["reaction_force"] < -limit]
    return (
        max(low, key=lambda f: f["q_actual"]) if low else None,
        min(high, key=lambda f: f["q_actual"]) if high else None,
    )


def trace_sides_ef(
    atoms,
    candidate_positions,
    candidate_energy,
    starts,
    direction,
    factory,
    config,
    ledger,
    task,
    *,
    settings,
    attempt=0,
):
    """Release both sides of a candidate, minimize on E/F, qualify endpoints, compare.

    ``starts`` is (low, high) from :func:`branch_starts`, either may be None; ``direction``
    is the search tangent used only for a missing side (+ for high, - for low)."""
    from prrs.reliability import GateRejected, GuardedCalculator
    from prrs.state import classify, encode, same_basin

    cfg = settings
    fmax = cfg.descent_fmax_eV_A or config.quench_fmax_eV_A
    candidate_positions = np.asarray(candidate_positions, float)
    rng = np.random.default_rng(cfg.seed)
    result = ConnectionResult(
        provenance={
            "algorithm": "HF-EF-connect-v1",
            "settings": asdict(cfg),
            "descent_fmax_eV_A": fmax,
            "numbers": atoms.numbers.tolist(),
            "candidate_energy_eV": float(candidate_energy),
        }
    )
    start_spent, started = ledger.spent, time.monotonic()

    def evaluate(positions):
        work = atoms.copy()
        work.set_constraint()
        work.positions = np.array(positions, float)
        work.calc = GuardedCalculator(factory(), config)
        return float(work.get_potential_energy()), work.get_forces()

    def descend(positions, hessian=None):
        work = atoms.copy()
        work.set_constraint()
        work.positions = np.array(positions, float)
        work.calc = GuardedCalculator(factory(), config)
        opt = BFGS(work, logfile=None, maxstep=cfg.max_step_A)
        if hessian is not None:
            opt.H0 = hessian.copy()
        energies = []
        converged = False
        for _ in opt.irun(fmax=fmax, steps=cfg.max_descent_steps):
            energies.append(float(work.get_potential_energy()))
            if float(np.linalg.norm(work.get_forces(), axis=1).max()) <= fmax:
                converged = True
                break
        forces = work.get_forces()
        rises = np.diff(energies) if len(energies) > 1 else np.zeros(1)
        learned = getattr(opt, "H", None)
        return {
            "converged": converged,
            "steps": int(opt.nsteps),
            "energy_eV": energies[-1],
            "fmax_eV_A": float(np.linalg.norm(forces, axis=1).max()),
            "max_energy_rise_eV": float(max(0.0, rises.max())),
            "energy_drop_eV": energies[0] - energies[-1],
            "positions": work.positions.copy(),
            "forces": forces.copy(),
            "hessian": None if learned is None else learned.copy(),
        }

    def tangent_start(sign):
        if direction is None:
            raise _Fail("no_start")
        d = np.asarray(direction, float).reshape(-1, 3)
        d = d / np.linalg.norm(d, axis=1).max()
        tried = []
        for a in cfg.tangent_amplitudes_A:
            x = candidate_positions + sign * a * d
            e, _ = evaluate(x)
            tried.append({"amplitude_A": a, "energy_eV": e})
            if candidate_energy - e > _quanta(candidate_energy):
                return x, {"source": "tangent", "sign": sign, "tried": tried}
        raise _Fail("tangent_no_resolved_drop")

    def stability(end):
        """Finite-perturbation re-descents from the endpoint, reusing its BFGS curvature."""
        ref = _frozen(atoms, end["positions"], end["energy_eV"], end["forces"])
        graph = encode(ref, config.bond_scale, active=config.active_atoms).edges
        floor = _quanta(end["energy_eV"])
        trials, tested, back = [], 0, 0
        for _ in range(cfg.n_perturb):
            x = end["positions"] + rng.normal(0.0, cfg.perturb_A, end["positions"].shape)
            e0, _ = evaluate(x)
            rise = e0 - end["energy_eV"]
            resolved = rise > _quanta(end["energy_eV"])
            again = descend(x, end["hessian"])
            other = _frozen(atoms, again["positions"], again["energy_eV"], again["forces"])
            same_graph = encode(other, config.bond_scale, active=config.active_atoms).edges
            moved = _rmsd(end["positions"], again["positions"]) * np.sqrt(len(atoms))
            allowed = floor + float(np.linalg.norm(end["forces"])) * moved
            dropped = again["energy_eV"] - end["energy_eV"] < -allowed
            same = bool(again["converged"] and same_graph == graph and not dropped)
            trials.append(
                {
                    "energy_rise_eV": rise,
                    "resolved": bool(resolved),
                    "converged": again["converged"],
                    "steps": again["steps"],
                    "returned": same,
                    "same_graph": bool(same_graph == graph),
                    "resolved_drop": bool(dropped),
                    "allowed_drop_eV": allowed,
                    "same_basin": bool(again["converged"] and same_basin(ref, other, config)),
                    "rmsd_A": _rmsd(end["positions"], again["positions"]),
                    "dE_eV": again["energy_eV"] - end["energy_eV"],
                }
            )
            if resolved:
                tested += 1
                back += same
        if tested == 0:
            verdict = "perturbation_unresolved"
        elif back == tested:
            verdict = "ef_stable"
        else:
            verdict = "ef_unstable"
        return verdict, trials

    ends = {}
    try:
        with ledger.segment(task=task, stage="connect_ef", attempt=attempt) as seg:
            try:
                for name, sign, frame in (("low", -1, starts[0]), ("high", 1, starts[1])):
                    side = {"side": name}
                    result.sides.append(side)
                    if frame is not None:
                        x0 = np.asarray(frame["positions"], float)
                        side["start"] = {
                            "source": "branch_frame",
                            "q_actual": frame["q_actual"],
                            "reaction_force": frame["reaction_force"],
                            "energy_eV": frame["energy_eV"],
                        }
                    else:
                        x0, side["start"] = tangent_start(sign)
                    side["start"]["rmsd_to_candidate_A"] = _rmsd(candidate_positions, x0)
                    end = descend(x0)
                    side["descent"] = {
                        k: v for k, v in end.items() if k not in ("forces", "hessian")
                    }
                    if not end["converged"]:
                        side["status"] = "descent_not_converged"
                        continue
                    verdict, trials = stability(end)
                    side["stability"] = {"verdict": verdict, "trials": trials}
                    side["status"] = verdict
                    ends[name] = _frozen(
                        atoms, end["positions"], end["energy_eV"], end["forces"]
                    )
                statuses = [s["status"] for s in result.sides]
                if any(s == "descent_not_converged" for s in statuses):
                    result.status = "failed:descent_not_converged"
                elif any(s == "ef_unstable" for s in statuses):
                    result.status = "unresolved:endpoint_unstable"
                elif any(s == "perturbation_unresolved" for s in statuses):
                    result.status = "unresolved:endpoint_untested"
                elif same_basin(ends["low"], ends["high"], config):
                    result.status = "same_basin"
                else:
                    label, fingerprint = classify(ends["low"], ends["high"], config)
                    result.status = "connected"
                    result.label, result.fingerprint = label, _jsonable(fingerprint)
            except _Fail as exc:
                result.status = f"failed:{exc}"
            except BudgetExhausted:
                result.status = "budget_exhausted"
            except GateRejected as exc:
                result.status = f"gate_rejected:{exc.code}"
            finally:
                seg.outcome = result.status
    finally:
        result.spent = ledger.spent - start_spent
        result.wall_s = time.monotonic() - started
    return result
