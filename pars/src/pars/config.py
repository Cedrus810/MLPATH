"""PARS configuration: one frozen dataclass, JSON in and out, unknown keys rejected.

Defaults are prospective (design sections 30/36): they are frozen before any
benchmark outcome is examined, and every non-obvious value carries a RATIONALE
entry. The configuration is echoed into every output record; changing it
silently would change the frontier.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, fields
import json
import math

RATIONALE: dict[str, str] = {
    "margin_gamma": (
        "design section 30 example value; few QM labels => strong frontier prior, "
        "more evidence => widen the rescue region. To be frozen prospectively "
        "before benchmark evaluation, not tuned afterwards"
    ),
    "cover_objective_tolerance": (
        "normalized-objective tolerance for the greedy cover; 0.10 keeps "
        "representatives within a tenth of the robust scale of the frontier"
    ),
    "tube_distance_threshold": (
        "Jaccard/Frechet threshold under which two path signature sequences count "
        "as the same reaction tube; mechanism diversity is preserved per tube, so "
        "this value trades compression against duplicate refinement risk"
    ),
    "dominance_tolerance": (
        "normalized-objective improvement below which two candidates count as "
        "indistinguishable; 0.0 is strict Pareto. Larger values only widen the "
        "frontier, and hard deletion then also needs at least this separation on "
        "every objective"
    ),
    "uncertainty_quantile": (
        "1.0 uses the worst frame (max) for the mlp_uncertainty objective; the "
        "design section 9 warns against the mean, which can hide a bad transition "
        "region inside a quiet path"
    ),
    "refinement_images": (
        "NEB image count, and the resolution of the path_curvature objective: "
        "curvature is measured on the aligned free segment resampled to this "
        "many uniform arc-length points, so vibrations finer than one image "
        "spacing do not count. "
        "Changing it changes curvature values and therefore the frontier"
    ),
    "refinement_remove_rigid_motion": (
        "remove rigid translation/rotation from the seed (frames aligned onto the "
        "first before resampling) and during NEB. Valid only when the energy is "
        "invariant under rigid motion, i.e. an isolated molecule, which is the "
        "scope of the PRRS runs ('finite nonperiodic'). Must be False for an "
        "embedded core (QM/MM with a fixed environment) or any external field. "
        "Without it, MD drift inside the seed made a p1 band spend its arc on "
        "rigid motion and 'converge' at 0.0002 eV for a 0.42 eV saddle"
    ),
    "refinement_spring_k": (
        "NEB spring constant in eV/A^2, the VTST Tools default (SPRING = -5). The "
        "springs are projected on the tangent (improved-tangent NEB), so k sets "
        "only the image spacing, not the converged path. The ASE default 0.1 let "
        "images settle into the basins with the saddle between two of them: 22 of "
        "27 p1 bands 'converged' at ~0.01 eV for a 0.42 eV proton transfer"
    ),
    "refinement_max_top_gap_ratio": (
        "a plain-NEB barrier is reported only if both gaps next to the highest "
        "image are at most this multiple of the band's median image gap; a "
        "coarser top means the saddle can lie between images and the barrier is "
        "an artefact. Not applied to an accepted climbing image, which sits on "
        "a stationary point; a separate Hessian check requires index one"
    ),
    "refinement_endpoint_min_rmsd": (
        "RMSD in Angstrom (Kabsch-aligned when refinement_remove_rigid_motion) "
        "below which two relaxed endpoints count as the same minimum only if "
        "the maximum mapped per-atom displacement is also below this tolerance; "
        "spectator atoms cannot dilute a local difference. This is a numerical "
        "geometry guard, not a chemical basin identity certificate"
    ),
    "refinement_uncertainty_refusal": (
        "when set, every refined image must provide finite nonnegative force "
        "uncertainty; missing/partial coverage or exceeding the threshold is "
        "UNRESOLVED_MLP. None disables the threshold, not invalid-value checks"
    ),
    "refinement_endpoint_hessian_step": "central force-difference step in Angstrom for independent endpoint curvature checks; 0.001 resolves local curvature without using PRRS identity code",
    "refinement_endpoint_curvature_tolerance": "negative-curvature noise floor in eV/A^2 after removing rigid modes when appropriate; 1e-4 is an explicit numerical acceptance floor, not QM certification",
    "refinement_endpoint_identity_rmsd": "both mapped RMSD and maximum per-atom displacement must be within this Angstrom tolerance for source/quench witnesses and signed-descent endpoint matching; spectators cannot dilute local mismatches",
    "refinement_saddle_hessian_step": "central force-difference step in Angstrom at the climbing image; adds 6N force evaluations for a local MLP saddle-index check, without certifying connectivity or QM accuracy",
    "refinement_saddle_curvature_tolerance": "negative-curvature noise floor in eV/A^2 for accepting exactly one unstable mode after optional rigid-mode projection; weak or additional unstable modes require further investigation",
    "refinement_descent_displacement": "Cartesian norm in Angstrom of each signed negative-mode kick from the climbing image; both starts must lower its energy",
    "refinement_descent_step": "initial steepest-descent step in A^2/eV, halved by energy backtracking; descent is independent of the NEB chain",
    "refinement_descent_max_displacement": "maximum per-atom displacement in Angstrom for a trial descent step, limiting excursions before backtracking",
    "refinement_descent_fmax": "descent force tolerance in eV/A, capped by the endpoint force tolerance; landed geometries must also pass curvature and mapped endpoint matching",
    "refinement_descent_max_steps": "hard accepted-step limit per signed descent; unfinished descent refuses CI success",
    "refinement_descent_max_backtracks": "maximum halvings per descent step after the initial trial; no decreasing step means refusal, not a stationary-point success",
    "refinement_descent_resolution_factor": "mandatory second signed-descent replay scales kick, gradient step and displacement cap by 0.25; matched endpoint labels must agree. Stability on two tested scales does not prove continuous IRC connectivity",
}


@dataclass(frozen=True)
class ParsConfig:
    """All PARS knobs. JSON is the configuration surface (no YAML)."""

    # --- candidate extraction (design 6.2/17.1) -------------------------------
    include_classifications: tuple[str, ...] = ("reactive",)
    include_phases: tuple[str, ...] = ("source", "free_start", "free")
    min_path_frames: int = 3

    # --- pre-refinement objectives (design 8) ---------------------------------
    uncertainty_quantile: float = 1.0

    # --- normalization (design 10) --------------------------------------------
    normalization: str = "robust"
    objective_scales: tuple[tuple[str, float], ...] = ()

    # --- dominance, margin, cover (design 11-16) ------------------------------
    dominance_mode: str = "uncertainty_aware"
    dominance_tolerance: float = 0.0
    margin_gamma: float = 0.15
    cover_objective_tolerance: float = 0.10
    preserve_mechanism_diversity: bool = True
    tube_distance_threshold: float = 0.30

    # --- audit sample for the later QM feedback experiment (design 20) --------
    audit_size: int = 5
    audit_seed: int = 0

    # --- refinement defaults (design 18) ---------------------------------------
    refinement_images: int = 9
    refinement_fmax_endpoints: float = 0.01
    refinement_fmax_neb: float = 0.05
    refinement_climbing_image: bool = True
    refinement_max_steps_endpoints: int = 500
    refinement_max_steps_neb: int = 1000
    refinement_endpoint_min_rmsd: float = 0.05
    refinement_spring_k: float = 5.0
    refinement_remove_rigid_motion: bool = True
    refinement_max_top_gap_ratio: float = 2.0
    refinement_uncertainty_refusal: float | None = None
    refinement_endpoint_hessian_step: float = 0.001
    refinement_endpoint_curvature_tolerance: float = 0.0001
    refinement_endpoint_identity_rmsd: float = 0.05
    refinement_saddle_hessian_step: float = 0.001
    refinement_saddle_curvature_tolerance: float = 0.0001
    refinement_descent_displacement: float = 0.05
    refinement_descent_step: float = 0.05
    refinement_descent_max_displacement: float = 0.05
    refinement_descent_fmax: float = 0.001
    refinement_descent_max_steps: int = 2000
    refinement_descent_max_backtracks: int = 20
    refinement_descent_resolution_factor: float = 0.25

    def __post_init__(self) -> None:
        for name, value in self.to_dict().items():
            if (
                name.startswith("refinement_")
                and isinstance(value, float)
                and not math.isfinite(value)
            ):
                raise ValueError(f"{name} must be finite")
        if (
            self.refinement_endpoint_hessian_step <= 0
            or self.refinement_endpoint_curvature_tolerance < 0
            or self.refinement_endpoint_identity_rmsd <= 0
        ):
            raise ValueError(
                "endpoint Hessian step/identity tolerance must be positive and curvature tolerance nonnegative"
            )
        if (
            self.refinement_saddle_hessian_step <= 0
            or self.refinement_saddle_curvature_tolerance < 0
        ):
            raise ValueError("saddle Hessian step must be positive and tolerance nonnegative")
        for name in (
            "refinement_descent_displacement",
            "refinement_descent_step",
            "refinement_descent_max_displacement",
            "refinement_descent_fmax",
        ):
            if getattr(self, name) <= 0:
                raise ValueError(f"{name} must be positive")
        for name in ("refinement_descent_max_steps", "refinement_descent_max_backtracks"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                raise ValueError(f"{name} must be a nonnegative integer")
        if not 0 < self.refinement_descent_resolution_factor < 1:
            raise ValueError("refinement_descent_resolution_factor must be in (0, 1)")
        if self.margin_gamma < 0:
            raise ValueError("margin_gamma must be >= 0")
        if self.dominance_tolerance < 0:
            raise ValueError("dominance_tolerance must be >= 0")
        if self.cover_objective_tolerance < 0:
            raise ValueError("cover_objective_tolerance must be >= 0")
        if self.tube_distance_threshold <= 0 or self.tube_distance_threshold > 1:
            raise ValueError("tube_distance_threshold must be in (0, 1]")
        if not 0.0 <= self.uncertainty_quantile <= 1.0:
            raise ValueError("uncertainty_quantile must be in [0, 1]")
        if self.refinement_images < 3:
            raise ValueError("refinement_images must be >= 3")
        if self.refinement_fmax_endpoints <= 0 or self.refinement_fmax_neb <= 0:
            raise ValueError("refinement fmax values must be > 0")
        if self.refinement_spring_k <= 0:
            raise ValueError("refinement_spring_k must be > 0")
        if self.refinement_max_top_gap_ratio < 1:
            raise ValueError("refinement_max_top_gap_ratio must be >= 1")
        if self.refinement_endpoint_min_rmsd <= 0:
            raise ValueError("refinement_endpoint_min_rmsd must be > 0")
        if self.min_path_frames < 2:
            raise ValueError("min_path_frames must be >= 2")
        if self.dominance_mode not in ("nominal", "uncertainty_aware"):
            raise ValueError("dominance_mode must be 'nominal' or 'uncertainty_aware'")
        if self.normalization != "robust":
            raise ValueError("normalization must be 'robust' (the only V1 method)")
        if self.refinement_uncertainty_refusal is not None and (
            self.refinement_uncertainty_refusal <= 0
        ):
            raise ValueError("refinement_uncertainty_refusal must be > 0 when set")

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, payload: dict) -> ParsConfig:
        if not isinstance(payload, dict):
            raise ValueError("ParsConfig payload must be a JSON object")
        known = {f.name for f in fields(cls)}
        unknown = set(payload) - known
        if unknown:
            raise ValueError(f"unknown ParsConfig keys: {sorted(unknown)}")
        values = dict(payload)
        for name in ("include_classifications", "include_phases"):
            if name in values:
                values[name] = tuple(values[name])
        if "objective_scales" in values:
            values["objective_scales"] = tuple(
                (str(k), float(v)) for k, v in values["objective_scales"]
            )
        return cls(**values)


def load_config(path) -> ParsConfig:
    from .iox import read_json

    return ParsConfig.from_dict(read_json(path))


def config_to_json_text(config: ParsConfig) -> str:
    return json.dumps(config.to_dict(), indent=2, allow_nan=False)
