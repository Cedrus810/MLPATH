"""ASE NEB adapter: the first RefinementBackend (design section 18).

Staged logic: checked endpoint minima -> ordinary NEB -> climbing-image
continuation. Guards keep the result honest:

- Pathology guard: a climbing image started on an unstable band can run away
  onto a repulsive wall. The CI stage is only accepted when it converges, its
  highest-energy image is interior, and that image stayed bracketed by the
  stable band (see ``_climbing_stage_accepted``); otherwise the plain-NEB band
  is reported with an explicit note (design 18: monitor pathological image
  crossing, never silently declare success). The barrier itself is NOT a
  guard: plain-NEB images straddle the saddle, so a correct climbing image
  raises the barrier by construction.
- Resolution guard: a band that falls back to plain NEB (no climbing image, or
  the climbing stage was rejected) reports a barrier only when the gaps next
  to its highest image are comparable to the band's typical spacing
  (``_plain_band_resolved``). Otherwise the saddle can lie between two images,
  the "barrier" is an artefact, and the result is FAILED without a barrier.
- MLP reliability gate: an explicit threshold requires complete image coverage
  of ``force_uncertainty_eV_A``; absent/partial coverage or excessive uncertainty
  is UNRESOLVED_MLP. Nonfinite or negative values are FAILED even without a threshold.
- Saddle gate: an accepted climbing image must have exactly one negative
  projected Hessian eigenvalue beyond the configured noise floor. Otherwise
  refinement fails without a barrier.
- Hessian sampling gate: every displaced force evaluation of every curvature
  check (endpoints, basin witnesses, climbing image, descent minima) passes the
  same uncertainty checks as band images; its maximum enters the reported one.
- Electronic state: every structure handed to the calculator carries the
  candidate's charge/spin in ``atoms.info``; a hand-built candidate without one
  is recorded as ``unspecified`` (calculator defaults), never guessed.
- Connectivity gate: independent signed negative-mode kicks must descend with
  strictly decreasing energy to checked minima uniquely matching the two
  refined endpoints. This Cartesian MLP witness is not mass-weighted IRC or
  validation of the PES against QM.

Endpoints must relax (to refinement_fmax_endpoints within the step budget),
pass a projected force-difference Hessian check, and match independently relaxed
source/quench basin witnesses when present. They must not fall into the same minimum:
a last frame that slides back into
the reactant basin would otherwise yield a converged A->A band reported as a
reaction path. Both cases are FAILED with the reason.

The initial chain is the arc-length resampled candidate trajectory with its
endpoints replaced by freshly minimized copies; interpolation between raw
endpoints is NOT used because the candidate already contains the geometric
memory of the reaction (design 17.4/17.5).
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass

import numpy as np
from ase import Atoms
from ase.calculators.calculator import Calculator
from ase.optimize import FIRE, LBFGS

try:  # ase >= 3.29
    from ase.mep.neb import NEB
except ImportError:  # pragma: no cover - older ase
    from ase.neb import NEB

from .. import electronic
from ..config import ParsConfig
from ..paths.resample import align_frames, kabsch_align, resample_positions
from ..schema import PathCandidate, RefinedPath, RefinementStatus
from .descent import descend

_UNCERTAINTY_KEY = "force_uncertainty_eV_A"


def _mapped_geometry_errors(positions, reference, remove_rigid_motion):
    """Mapped distances; the maximum cannot be diluted by spectator atoms."""
    compared = kabsch_align(positions, reference) if remove_rigid_motion else positions
    distances = np.linalg.norm(compared - reference, axis=1)
    return float(np.sqrt(np.mean(distances**2))), float(distances.max())


# Henkelman-Jonsson improved tangent: springs act along the tangent only, so a stiff
# spring constant equalizes spacing without cutting corners
_NEB_METHOD = "improvedtangent"


def _connection_branches(report):
    """Uncertainty evidence includes partially completed finer-scale branches."""
    return report["branches"] + report.get("resolution_check", {}).get("fine", {}).get(
        "branches", []
    )


def _uncertainty_max(provenance, band):
    """Largest uncertainty behind the result: band images, every Hessian
    sample and every descent evaluation recorded so far."""
    reports = []
    for check in provenance.get("endpoint_checks", {}).values():
        reports += [check, check.get("basin_reference_check", {})]
    reports.append(provenance.get("saddle_check", {}))
    connection = provenance.get("connectivity_check", {})
    if "branches" in connection:
        for branch in _connection_branches(connection):
            reports += [branch, branch.get("minimum_check", {})]
    values = [
        r["uncertainty_max_eV_A"] for r in reports if r.get("uncertainty_max_eV_A") is not None
    ]
    if band is not None:
        values.append(band)
    return max(values, default=None)


def _plain_band_resolved(positions: np.ndarray, energies: np.ndarray, ratio: float) -> bool:
    """Whether a plain-NEB band resolves its top: both gaps next to the
    highest-energy image are at most ``ratio`` times the median image gap.

    Plain NEB images straddle the saddle; their highest energy is a barrier
    estimate only when images near the top are as close as elsewhere. Big
    gaps at the ends (relaxed endpoints far from their raw frames) do not
    matter, a big gap at the top does.
    """
    flat = positions.reshape(len(positions), -1)
    gaps = np.linalg.norm(np.diff(flat, axis=0), axis=1)
    top = int(np.argmax(energies))
    if not 0 < top < len(energies) - 1:
        return False
    adjacent = [gaps[i] for i in (top - 1, top) if 0 <= i < len(gaps)]
    return max(adjacent) <= ratio * float(np.median(gaps))


def _climbing_stage_accepted(
    plain_positions: np.ndarray,
    ci_positions: np.ndarray,
    ci_energies: np.ndarray,
    converged: bool,
) -> bool:
    """Accept a climbing-image band only if it converged, its top image is
    interior, and that image stayed bracketed by the stable plain-NEB band.

    The saddle lies between the plain band's images around the top, so a
    healthy climbing image moves at most about one image spacing away from
    its plain-band counterpart; a runaway up a repulsive wall moves further.
    """
    if not converged:
        return False
    top = int(np.argmax(ci_energies))
    if not 0 < top < len(ci_energies) - 1:
        return False
    flat_plain = plain_positions.reshape(len(plain_positions), -1)
    moved = float(np.linalg.norm(ci_positions[top].ravel() - flat_plain[top]))
    bracket = max(
        float(np.linalg.norm(flat_plain[top] - flat_plain[top - 1])),
        float(np.linalg.norm(flat_plain[top + 1] - flat_plain[top])),
    )
    return moved <= bracket


@dataclass
class _BandMetrics:
    energies: np.ndarray
    barrier_eV: float | None
    max_image_force_eV_A: float
    neb_residual_eV_A: float


def _curvature(atoms, config, step, tolerance, report, include_mode=False):
    """Independent force-difference Hessian, with optional SE(3) projection.

    Each displaced force evaluation is checked where it happens: a nonfinite or
    negative uncertainty raises; with a refusal threshold, a missing or excessive
    one stops sampling with status ``unresolved_mlp`` and no eigenvalues.
    ``report`` is filled in place so partial evidence survives either outcome;
    the original coordinates are always restored.
    """
    report.update(
        status="sampling",
        step_A=step,
        negative_tolerance_eV_A2=tolerance,
        uncertainty_max_eV_A=None,
        uncertainty_observations=0,
        uncertainty_missing=0,
    )
    refusal = config.refinement_uncertainty_refusal
    original = atoms.positions.copy()
    size = original.size
    hessian = np.empty((size, size))

    def sampled_forces(offset, column, sign):
        atoms.positions = original + sign * offset.reshape(original.shape)
        forces = atoms.get_forces().ravel().copy()
        value = getattr(atoms.calc, "results", {}).get(_UNCERTAINTY_KEY)
        report["uncertainty_observations"] += 1
        report["uncertainty_missing"] += int(value is None)
        if value is not None:
            value = float(value)
            if not np.isfinite(value) or value < 0:
                report.update(
                    status="failed", message=f"invalid uncertainty at coordinate {column}"
                )
                raise ValueError(
                    "invalid Hessian sample force uncertainty: must be finite and nonnegative"
                )
            report["uncertainty_max_eV_A"] = max(report["uncertainty_max_eV_A"] or 0.0, value)
        if refusal is not None and (value is None or value > refusal):
            measured = "missing" if value is None else f"{value:.3g} eV/A"
            report.update(
                status="unresolved_mlp",
                message=(
                    f"Hessian sample force uncertainty {measured} at coordinate {column} "
                    f"(sign {sign:+d}) fails refusal threshold {refusal:.3g}"
                ),
            )
            return None
        return forces

    try:
        for column in range(size):
            offset = np.zeros(size)
            offset[column] = step
            plus = sampled_forces(offset, column, 1)
            if plus is None:
                return report
            minus = sampled_forces(offset, column, -1)
            if minus is None:
                return report
            hessian[:, column] = -(plus - minus) / (2 * step)
    finally:
        atoms.positions = original
        missing = report["uncertainty_missing"]
        report["uncertainty_coverage"] = (
            "unavailable"
            if missing == report["uncertainty_observations"]
            else "partial"
            if missing
            else "complete"
        )
    if not np.isfinite(hessian).all():
        report.update(status="failed", message="non-finite forces")
        raise ValueError("Hessian contains non-finite forces")
    hessian = (hessian + hessian.T) / 2
    removed = 0
    basis = np.eye(size)
    if config.refinement_remove_rigid_motion:
        centered = original - original.mean(axis=0)
        translations = [np.tile(axis, (len(atoms), 1)).ravel() for axis in np.eye(3)]
        rotations = [np.cross(axis, centered).ravel() for axis in np.eye(3)]
        u, singular, _ = np.linalg.svd(np.array(translations + rotations).T, full_matrices=True)
        removed = int(np.sum(singular > 1e-10 * singular.max()))
        basis = u[:, removed:]
        hessian = basis.T @ hessian @ basis
    if include_mode:
        eigenvalues, modes = np.linalg.eigh(hessian)
    else:
        eigenvalues = np.linalg.eigvalsh(hessian)
    lowest = float(eigenvalues.min()) if eigenvalues.size else 0.0
    report.update(
        status="checked",
        lowest_eigenvalue_eV_A2=lowest,
        eigenvalues_eV_A2=eigenvalues.tolist(),
        negative_modes=int(np.sum(eigenvalues < -tolerance)),
        removed_rigid_modes=removed,
    )
    if include_mode and report["negative_modes"] == 1:
        mode = basis @ modes[:, 0]
        report["negative_mode_cartesian"] = mode.reshape(original.shape).tolist()
    return report


def _endpoint_curvature(atoms, config, report):
    return _curvature(
        atoms,
        config,
        config.refinement_endpoint_hessian_step,
        config.refinement_endpoint_curvature_tolerance,
        report,
    )


class AseNebBackend:
    backend = "ase-neb"

    def __init__(self, calculator: Calculator):
        self.calculator = calculator

    def refine(self, candidate: PathCandidate, config: ParsConfig) -> RefinedPath:
        climbing = config.refinement_climbing_image
        provenance = {
            "backend": self.backend,
            "images": config.refinement_images,
            "fmax_endpoints": config.refinement_fmax_endpoints,
            "fmax_neb": config.refinement_fmax_neb,
            "climbing_image": climbing,
            "neb_method": _NEB_METHOD,
            "spring_k": config.refinement_spring_k,
            "remove_rigid_motion": config.refinement_remove_rigid_motion,
            "candidate_id": candidate.candidate_id,
            "source": candidate.provenance.to_dict(),
        }
        try:
            state = candidate.provenance.electronic_state
            if state is not None:
                state = electronic.checked_state(state, candidate.candidate_id)
            provenance["electronic_state"] = "unspecified" if state is None else state
            checks = provenance["endpoint_checks"] = {}
            chain, endpoint_problem, unresolved = self._build_chain(
                candidate, config, state, checks
            )
            if endpoint_problem:
                return self._result(
                    candidate,
                    chain,
                    RefinementStatus.UNRESOLVED_MLP if unresolved else RefinementStatus.FAILED,
                    climbing,
                    provenance,
                    uncertainty=_uncertainty_max(provenance, None),
                    message=(
                        f"{endpoint_problem}; not reporting success; no barrier reported"
                        if unresolved
                        else endpoint_problem
                    ),
                )
            neb = NEB(
                chain,
                k=config.refinement_spring_k,
                climb=False,
                allow_shared_calculator=True,
                method=_NEB_METHOD,
                remove_rotation_and_translation=config.refinement_remove_rigid_motion,
            )
            # FIRE: first-order and robust to the force discontinuities that
            # tangent reconstruction and climbing-image switching introduce
            optimizer = FIRE(neb, logfile=None)
            optimizer.run(
                fmax=config.refinement_fmax_neb,
                steps=config.refinement_max_steps_neb,
            )
            if not optimizer.converged():
                return self._result(
                    candidate,
                    chain,
                    RefinementStatus.FAILED,
                    climbing,
                    provenance,
                    message=(
                        f"NEB did not reach fmax={config.refinement_fmax_neb} "
                        f"within {config.refinement_max_steps_neb} steps"
                    ),
                )
            plain = self._metrics(neb, chain)
            snapshot = np.array([atoms.positions.copy() for atoms in chain])

            status = RefinementStatus.NEB_CONVERGED
            message = ""
            metrics = plain
            if climbing:
                neb.climb = True
                ci_optimizer = FIRE(neb, logfile=None)
                ci_optimizer.run(
                    fmax=config.refinement_fmax_neb,
                    steps=config.refinement_max_steps_neb,
                )
                ci = self._metrics(neb, chain)
                ci_positions = np.array([atoms.positions.copy() for atoms in chain])
                if _climbing_stage_accepted(
                    snapshot, ci_positions, ci.energies, ci_optimizer.converged()
                ):
                    status = RefinementStatus.CI_NEB_CONVERGED
                    metrics = ci
                else:
                    message = (
                        "climbing-image stage rejected (not converged, top image at "
                        "an endpoint, or the climbing image left the stable band; "
                        f"barrier {ci.barrier_eV:.3g} vs {plain.barrier_eV:.3g} eV); "
                        "reporting the plain NEB band"
                    )
                    neb.climb = False
                    for atoms, positions in zip(chain, snapshot):
                        atoms.positions = positions
                    metrics = self._metrics(neb, chain)

            # Reject unreliable bands before spending on Hessian/descent witnesses.
            uncertainty, coverage = self._path_uncertainty(chain)
            provenance["uncertainty_coverage"] = coverage
            refusal = config.refinement_uncertainty_refusal
            if refusal is not None and (coverage != "complete" or uncertainty > refusal):
                reason = (
                    f"force uncertainty {coverage}; every image is required by the configured refusal threshold"
                    if coverage != "complete"
                    else f"transition-region force uncertainty {uncertainty:.3g} eV/A exceeds refusal threshold {refusal:.3g}"
                )
                return self._result(
                    candidate,
                    chain,
                    RefinementStatus.UNRESOLVED_MLP,
                    climbing,
                    provenance,
                    metrics=dataclasses.replace(metrics, barrier_eV=None),
                    uncertainty=_uncertainty_max(provenance, uncertainty),
                    message=f"{reason}; not reporting success; no barrier reported",
                )

            if status is RefinementStatus.CI_NEB_CONVERGED:
                top = int(np.argmax(metrics.energies))
                saddle_check = {"image_index": top, "connectivity": "not_checked"}
                provenance["saddle_check"] = saddle_check
                _curvature(
                    chain[top],
                    config,
                    config.refinement_saddle_hessian_step,
                    config.refinement_saddle_curvature_tolerance,
                    saddle_check,
                    include_mode=True,
                )
                if saddle_check["status"] == "unresolved_mlp":
                    return self._result(
                        candidate,
                        chain,
                        RefinementStatus.UNRESOLVED_MLP,
                        climbing,
                        provenance,
                        metrics=dataclasses.replace(metrics, barrier_eV=None),
                        uncertainty=_uncertainty_max(provenance, uncertainty),
                        message=(
                            f"climbing image {saddle_check['message']}; not reporting "
                            "success; no barrier reported"
                        ),
                    )
                if saddle_check["negative_modes"] != 1:
                    return self._result(
                        candidate,
                        chain,
                        RefinementStatus.FAILED,
                        climbing,
                        provenance,
                        metrics=dataclasses.replace(metrics, barrier_eV=None),
                        message=(
                            "climbing image is not a resolved index-one saddle: "
                            f"{saddle_check['negative_modes']} negative modes beyond "
                            f"{config.refinement_saddle_curvature_tolerance:g} eV/A^2; "
                            "no barrier reported"
                        ),
                    )

                connection = {"status": "checking", "branches": []}
                provenance["connectivity_check"] = connection
                try:
                    self._check_connectivity(
                        chain[top],
                        saddle_check["negative_mode_cartesian"],
                        [chain[0], chain[-1]],
                        config,
                        connection,
                    )
                except Exception as error:  # noqa: BLE001 - preserve partial descent/budget evidence
                    connection.update(
                        status="failed", message=f"{type(error).__name__}: {error}"
                    )
                saddle_check["connectivity"] = connection["status"]
                if connection["status"] != "matched":
                    return self._result(
                        candidate,
                        chain,
                        RefinementStatus.UNRESOLVED_MLP
                        if connection["status"] == "unresolved_mlp"
                        else RefinementStatus.FAILED,
                        climbing,
                        provenance,
                        metrics=dataclasses.replace(metrics, barrier_eV=None),
                        uncertainty=_uncertainty_max(provenance, uncertainty),
                        message=(
                            f"climbing image connectivity check: {connection['message']}; "
                            "not reporting success; no barrier reported"
                        ),
                    )
                connection["endpoint_identity_scope"] = (
                    "source_basin_witnesses"
                    if all(
                        checks[name].get("basin_identity") == "matched"
                        for name in ("reactant", "product")
                    )
                    else "refined_endpoints_only"
                )

            uncertainty = _uncertainty_max(provenance, uncertainty)
            if status is RefinementStatus.NEB_CONVERGED:
                provenance["connectivity_check"] = {
                    "status": "not_checked",
                    "reason": "plain_neb_band_estimate",
                }
                band = np.array([atoms.positions for atoms in chain])
                if not _plain_band_resolved(
                    band, metrics.energies, config.refinement_max_top_gap_ratio
                ):
                    note = f"{message}; " if message else ""
                    return self._result(
                        candidate,
                        chain,
                        RefinementStatus.FAILED,
                        climbing,
                        provenance,
                        metrics=dataclasses.replace(metrics, barrier_eV=None),
                        message=(
                            f"{note}plain NEB band unresolved at its top: the highest image "
                            "is at an endpoint or an adjacent gap exceeds "
                            f"{config.refinement_max_top_gap_ratio}x the median image "
                            "gap, so the saddle can lie between images; no barrier "
                            "reported"
                        ),
                    )

            return self._result(
                candidate,
                chain,
                status,
                climbing,
                provenance,
                metrics=metrics,
                uncertainty=uncertainty,
                message=message,
            )
        except Exception as error:  # noqa: BLE001 - the backend reports, not raises
            return RefinedPath(
                candidate_id=candidate.candidate_id,
                backend=self.backend,
                status=RefinementStatus.FAILED,
                numbers=candidate.numbers,
                message=f"{type(error).__name__}: {error}",
                climbing_image=climbing,
                provenance=provenance,
            )

    # -- internals ------------------------------------------------------------

    def _check_connectivity(self, saddle, mode, endpoints, config, report):
        """Require matching endpoint witnesses at two frozen descent resolutions."""
        self._check_connectivity_at_scale(saddle, mode, endpoints, config, report)
        report["continuous_connection_certified"] = False
        if report["status"] != "matched":
            return
        factor = config.refinement_descent_resolution_factor
        finer = dataclasses.replace(
            config,
            refinement_descent_displacement=config.refinement_descent_displacement * factor,
            refinement_descent_step=config.refinement_descent_step * factor,
            refinement_descent_max_displacement=config.refinement_descent_max_displacement
            * factor,
        )
        fine = {"status": "checking", "branches": []}
        resolution = {"factor": factor, "status": "checking", "fine": fine}
        report["resolution_check"] = resolution
        try:
            self._check_connectivity_at_scale(saddle, mode, endpoints, finer, fine)
        except Exception as error:  # noqa: BLE001 - retain resolution/budget diagnostics
            fine.update(status="failed", message=f"{type(error).__name__}: {error}")
        labels = [b.get("matched_endpoint") for b in report["branches"]]
        fine_labels = [b.get("matched_endpoint") for b in fine["branches"]]
        if fine["status"] != "matched" or labels != fine_labels:
            resolution["status"] = "sensitive_or_unresolved"
            report.update(
                status="unresolved_mlp" if fine["status"] == "unresolved_mlp" else "failed",
                message=f"connectivity changes or is unresolved at finer descent resolution: {fine.get('message', 'endpoint labels differ')}",
            )
            return
        resolution["status"] = "stable_on_tested_scales"

    def _check_connectivity_at_scale(self, saddle, mode, endpoints, config, report):
        """Independently descend both signs, check minima and a unique pairing."""
        report["method"] = "energy_monotone_cartesian_gradient_descent"
        report["settings"] = {
            name: value
            for name, value in config.to_dict().items()
            if name.startswith("refinement_descent_")
        }
        report["endpoint_match_rmsd_A"] = config.refinement_endpoint_identity_rmsd
        saddle_energy = float(saddle.get_potential_energy())
        report["saddle_energy_eV"] = saddle_energy
        mode = np.asarray(mode, dtype=float)
        mode = mode / np.linalg.norm(mode)
        for sign in (-1, 1):
            branch = {"sign": sign}
            report["branches"].append(branch)
            atoms = saddle.copy()
            atoms.calc = self.calculator
            atoms.positions = (
                saddle.positions + sign * config.refinement_descent_displacement * mode
            )
            if not descend(atoms, saddle_energy, config, branch):
                report.update(status=branch["status"], message=branch["message"])
                return
            minimum = branch["minimum_check"] = {}
            _endpoint_curvature(atoms, config, minimum)
            if minimum["status"] == "unresolved_mlp":
                branch.update(
                    status="unresolved_mlp", message=f"descent minimum {minimum['message']}"
                )
                report.update(status="unresolved_mlp", message=branch["message"])
                return
            if minimum["negative_modes"]:
                branch.update(status="failed", message="descent stopped at negative curvature")
                report.update(status="failed", message=branch["message"])
                return
            branch["status"] = "minimum_checked"
            branch["endpoint_rmsd_A"] = {}
            branch["endpoint_max_displacement_A"] = {}
            for name, endpoint in zip(("reactant", "product"), endpoints):
                rmsd, maximum = _mapped_geometry_errors(
                    atoms.positions, endpoint.positions, config.refinement_remove_rigid_motion
                )
                branch["endpoint_rmsd_A"][name] = rmsd
                branch["endpoint_max_displacement_A"][name] = maximum
        branches = report["branches"]
        pairings = [
            names
            for names in (("reactant", "product"), ("product", "reactant"))
            if all(
                b["endpoint_rmsd_A"][name] <= config.refinement_endpoint_identity_rmsd
                and b["endpoint_max_displacement_A"][name]
                <= config.refinement_endpoint_identity_rmsd
                for b, name in zip(branches, names)
            )
        ]
        if len(pairings) != 1:
            report.update(
                status="failed",
                message="signed descents do not uniquely match both endpoint basins",
            )
            return
        for branch, name in zip(branches, pairings[0]):
            branch["matched_endpoint"] = name
        report.update(
            status="matched", message="two checked minima match distinct endpoint basins"
        )

    def _build_chain(
        self, candidate: PathCandidate, config: ParsConfig, state: dict | None, checks: dict
    ) -> tuple[list[Atoms], str, bool]:
        """The seed chain, a non-empty reason when its endpoints are unusable, and
        whether every such reason is an MLP-reliability refusal. ``checks`` is
        filled in place with the endpoint curvature evidence."""
        template = Atoms(numbers=candidate.numbers, positions=candidate.positions[0])
        if state is not None:
            # Atoms.copy() carries info: chain images, witnesses and descents inherit it
            electronic.stamp(template.info, state)

        def atoms_for(positions: np.ndarray) -> Atoms:
            atoms = template.copy()
            atoms.calc = self.calculator
            atoms.positions = positions
            return atoms

        rigid = config.refinement_remove_rigid_motion
        frames = align_frames(candidate.positions) if rigid else candidate.positions
        endpoints = []
        problems = []
        refusals = 0
        references = candidate.provenance.extra.get("endpoint_basin_references", {})
        for name, frame in (("reactant", frames[0]), ("product", frames[-1])):
            atoms = atoms_for(frame)
            converged = LBFGS(atoms, logfile=None).run(
                fmax=config.refinement_fmax_endpoints,
                steps=config.refinement_max_steps_endpoints,
            )
            if not converged:
                problems.append(
                    f"{name} endpoint relaxation did not reach "
                    f"fmax={config.refinement_fmax_endpoints} within "
                    f"{config.refinement_max_steps_endpoints} steps"
                )
            else:
                checked = checks[name] = {}
                _endpoint_curvature(atoms, config, checked)
                if checked["status"] == "unresolved_mlp":
                    refusals += 1
                    problems.append(f"{name} endpoint {checked['message']}")
                elif (
                    checked["lowest_eigenvalue_eV_A2"]
                    < -config.refinement_endpoint_curvature_tolerance
                ):
                    problems.append(
                        f"{name} endpoint has negative curvature ({checked['lowest_eigenvalue_eV_A2']:.6g} eV/A^2); not a minimum"
                    )
                checked["basin_identity"] = "unavailable"
                if name in references and not problems:
                    reference = np.asarray(references[name], dtype=float)
                    if reference.shape != frame.shape or not np.isfinite(reference).all():
                        raise ValueError(f"invalid {name} endpoint basin reference")
                    witness = atoms_for(reference)
                    if not LBFGS(witness, logfile=None).run(
                        fmax=config.refinement_fmax_endpoints,
                        steps=config.refinement_max_steps_endpoints,
                    ):
                        problems.append(f"{name} basin reference did not converge")
                    else:
                        witness_check = checked["basin_reference_check"] = {}
                        _endpoint_curvature(witness, config, witness_check)
                        if witness_check["status"] == "unresolved_mlp":
                            refusals += 1
                            problems.append(
                                f"{name} basin reference {witness_check['message']}"
                            )
                        elif (
                            witness_check["lowest_eigenvalue_eV_A2"]
                            < -config.refinement_endpoint_curvature_tolerance
                        ):
                            problems.append(f"{name} basin reference has negative curvature")
                        rmsd, maximum = _mapped_geometry_errors(
                            witness.positions, atoms.positions, rigid
                        )
                        checked["basin_reference_rmsd_A"] = rmsd
                        checked["basin_reference_max_displacement_A"] = maximum
                        checked["basin_identity"] = (
                            "matched"
                            if rmsd <= config.refinement_endpoint_identity_rmsd
                            and maximum <= config.refinement_endpoint_identity_rmsd
                            else "mismatch"
                        )
                        if checked["basin_identity"] == "mismatch":
                            problems.append(
                                f"{name} endpoint basin identity mismatch (mapped RMSD {rmsd:.6g} A, maximum displacement {maximum:.6g} A)"
                            )
            # relaxation may drift the molecule; keep it in the seed's frame
            relaxed = atoms.positions.copy()
            endpoints.append(kabsch_align(relaxed, frame) if rigid else relaxed)
        interior = resample_positions(frames, config.refinement_images)[1:-1]
        chain = [atoms_for(endpoints[0])]
        chain.extend(atoms_for(frame) for frame in interior)
        chain.append(atoms_for(endpoints[1]))
        if not problems:
            # rigid motion only separates minima when it is not an invariance
            rmsd, maximum = _mapped_geometry_errors(endpoints[1], endpoints[0], rigid)
            if (
                rmsd < config.refinement_endpoint_min_rmsd
                and maximum < config.refinement_endpoint_min_rmsd
            ):
                problems.append(
                    f"reactant and product endpoints relaxed into the same minimum "
                    f"(aligned RMSD {rmsd:.3g} A < {config.refinement_endpoint_min_rmsd} A)"
                )
        return chain, "; ".join(problems), bool(problems) and refusals == len(problems)

    def _metrics(self, neb: NEB, chain: list[Atoms]) -> _BandMetrics:
        energies = np.asarray([atoms.get_potential_energy() for atoms in chain])
        max_force = float(max(np.abs(atoms.get_forces()).max() for atoms in chain))
        neb_forces = np.asarray(neb.get_forces())
        residual = max(
            (
                float(np.linalg.norm(image_forces, axis=1).max())
                for image_forces in neb_forces.reshape(len(neb_forces), -1, 3)
            ),
            default=0.0,
        )
        return _BandMetrics(
            energies=energies,
            barrier_eV=float(energies.max() - energies[0]),
            max_image_force_eV_A=max_force,
            neb_residual_eV_A=residual,
        )

    def _path_uncertainty(self, chain: list[Atoms]) -> tuple[float | None, str]:
        """Maximum and complete/partial/unavailable coverage; invalid values raise."""
        values = []
        for atoms in chain:
            atoms.get_forces()
            value = getattr(self.calculator, "results", {}).get(_UNCERTAINTY_KEY)
            if value is not None and (not np.isfinite(float(value)) or float(value) < 0):
                raise ValueError("invalid force uncertainty: must be finite and nonnegative")
            values.append(None if value is None else float(value))
        if all(value is None for value in values):
            return None, "unavailable"
        coverage = "partial" if any(v is None for v in values) else "complete"
        return float(max(v for v in values if v is not None)), coverage

    def _result(
        self,
        candidate: PathCandidate,
        chain: list[Atoms],
        status: RefinementStatus,
        climbing: bool,
        provenance: dict,
        metrics: _BandMetrics | None = None,
        uncertainty: float | None = None,
        message: str = "",
    ) -> RefinedPath:
        positions = np.asarray([atoms.positions for atoms in chain])
        return RefinedPath(
            candidate_id=candidate.candidate_id,
            backend=self.backend,
            status=status,
            numbers=candidate.numbers,
            positions=positions,
            energies=None if metrics is None else metrics.energies,
            barrier_eV=None if metrics is None else metrics.barrier_eV,
            max_image_force_eV_A=(None if metrics is None else metrics.max_image_force_eV_A),
            neb_residual_eV_A=(None if metrics is None else metrics.neb_residual_eV_A),
            uncertainty_max_eV_A=uncertainty,
            climbing_image=climbing,
            message=message,
            provenance=provenance,
        )
