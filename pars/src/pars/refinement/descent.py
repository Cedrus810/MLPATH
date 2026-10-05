"""Energy-monotone Cartesian descent for a local MLP connectivity witness.

This is neither mass-weighted IRC nor QM validation. The caller supplies a
signed negative-mode start and verifies curvature and mapped basin identity.
The report is populated in place so refused/incomplete branches remain auditable.
All trials use the same calculator and its existing force-evaluation budget.
"""

from __future__ import annotations

import numpy as np


def descend(atoms, saddle_energy, config, report):
    report.update(
        status="running",
        start_positions_A=atoms.positions.tolist(),
        energies_eV=[],
        accepted_steps=0,
        backtracks=0,
        uncertainty_max_eV_A=None,
        uncertainty_observations=0,
        uncertainty_missing=0,
    )
    fmax = min(config.refinement_descent_fmax, config.refinement_fmax_endpoints)
    report["fmax_target_eV_A"] = fmax

    def observe():
        energy = float(atoms.get_potential_energy())
        forces = np.asarray(atoms.get_forces(), dtype=float)
        if not np.isfinite(energy) or not np.isfinite(forces).all():
            raise ValueError("non-finite descent energy or forces")
        uncertainty = atoms.calc.results.get("force_uncertainty_eV_A")
        report["uncertainty_observations"] += 1
        report["uncertainty_missing"] += int(uncertainty is None)
        missing = report["uncertainty_missing"]
        report["uncertainty_coverage"] = (
            "unavailable"
            if missing == report["uncertainty_observations"]
            else "partial"
            if missing
            else "complete"
        )
        if uncertainty is not None:
            uncertainty = float(uncertainty)
            if not np.isfinite(uncertainty) or uncertainty < 0:
                raise ValueError("invalid descent force uncertainty")
            previous = report["uncertainty_max_eV_A"]
            report["uncertainty_max_eV_A"] = max(previous or 0.0, uncertainty)
        refusal = config.refinement_uncertainty_refusal
        if refusal is not None and (uncertainty is None or uncertainty > refusal):
            report.update(
                status="unresolved_mlp",
                message=f"descent force uncertainty {report['uncertainty_coverage']} or above refusal threshold",
            )
            return None
        return energy, forces

    try:
        observed = observe()
        if observed is None:
            return False
        energy, forces = observed
        report["energies_eV"].append(energy)
        if energy >= saddle_energy:
            report.update(
                status="failed", message="signed negative-mode start did not lower energy"
            )
            return False
        for step in range(config.refinement_descent_max_steps + 1):
            report["positions_A"] = atoms.positions.tolist()
            report["max_force_eV_A"] = float(np.linalg.norm(forces, axis=1).max())
            if report["max_force_eV_A"] <= fmax:
                report["status"] = "force_converged"
                return True
            if step == config.refinement_descent_max_steps:
                report.update(status="failed", message="descent accepted-step budget exhausted")
                return False
            original = atoms.positions.copy()
            displacement = config.refinement_descent_step * forces
            largest = float(np.linalg.norm(displacement, axis=1).max())
            displacement *= min(1.0, config.refinement_descent_max_displacement / largest)
            accepted = False
            for backtrack in range(config.refinement_descent_max_backtracks + 1):
                atoms.positions = original + 2.0 ** (-backtrack) * displacement
                observed = observe()
                if observed is None:
                    return False
                next_energy, next_forces = observed
                if next_energy < energy:
                    energy, forces = next_energy, next_forces
                    report["energies_eV"].append(energy)
                    report["accepted_steps"] += 1
                    accepted = True
                    break
                report["backtracks"] += 1
            if not accepted:
                atoms.positions = original
                report.update(
                    status="failed",
                    message="no energy-decreasing descent step within backtracking budget",
                )
                return False
    except Exception as error:  # noqa: BLE001 - persist failures including calculator budgets
        report.update(status="failed", message=f"{type(error).__name__}: {error}")
        return False
