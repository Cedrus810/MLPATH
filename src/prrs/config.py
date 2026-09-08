"""Validated, JSON-serializable search parameters; eV, Angstrom, fs, kelvin."""

from dataclasses import asdict, dataclass, fields
import math

# Which numeric fields get which bound. Kept at module level rather than inside
# __post_init__ so a test can read them: every int and float field of SearchConfig must
# appear in exactly one of these, and test_validation pins that. The lists had drifted
# behind the dataclass once already -- ten fields added for the min-mode walk, the
# response descent and the saddle search were never added to any of them, and
# min_mode_refresh = 0 reached `step % config.min_mode_refresh` in follow_min_mode as a
# plain ZeroDivisionError deep inside a run. Naming the buckets here is what makes the
# omission a failing test instead of a crash an hour in.
INTEGER_POSITIVE = (
    "max_trials",
    "max_directions",
    "response_steps",
    "irc_max_per_run",
    "irc_relay_hops",
    "boundary_continuation_max_steps",
    "boundary_continuation_seeds",
    "boundary_continuation_per_run",
    "sample_interval",
    "quench_steps",
    "chemical_max_nodes",
    "chemical_trials_per_node",
    "conformer_max_per_node",
    "conformer_trials_per_microstate",
    "automorphism_limit",
    "response_samples",
    "response_descent_max_steps",
    "analytic_hessian_check_directions",
    "min_mode_max_steps",
    "min_mode_refresh",
    "saddle_search_per_run",
    "saddle_search_seeds",
    "free_alignment_samples",
    "min_pair_hops",
    "min_directions_per_family",
)
INTEGER_NONNEGATIVE = (
    "seed",
    "max_depth",
    "refinement_steps",
    "pulse_steps",
    "minimum_check_probes",
    "soft_polish_rounds",
    "soft_mode_scan_points",
    "soft_polish_steps",
    # Zero backtracks is a legitimate configuration: take the step or give up.
    "min_mode_backtracks",
)
FLOAT_POSITIVE = (
    "irc_step_A",
    "torsion_tolerance_rad",
    "soft_mode_curvature_max_eV_rad2",
    "soft_mode_curvature_floor_eV_rad2",
    "soft_mode_step_rad",
    "torsion_symmetry_margin_rad",
    "collateral_bond_A",
    "timestep_fs",
    "quench_fmax_eV_A",
    "minimum_check_displacement_A",
    "minimum_check_step_A",
    "minimum_check_eigenvalue_tol",
    "pair_cutoff_A",
    "hbond_acceptor_A",
    "hbond_angle_deg",
    "direction_quota_fraction",
    "bond_scale",
    "parity_tolerance",
    "basin_rmsd_A",
    "basin_energy_eV",
    "min_distance_ratio",
    "analytic_hessian_check",
    "trivial_floor_fraction",
    "internal_mode_recheck_step_scale",
    "response_significance",
    "response_descent_step",
    "response_descent_tol",
    "min_mode_step_A",
    "min_mode_fmax_eV_A",
    "min_mode_overlap",
    "saddle_search_max_rise_eV",
    "boundary_continuation_step_A",
    "max_force_eV_A",
    "max_energy_rise_eV_atom",
    "max_nve_drift_eV_atom",
)
FLOAT_NONNEGATIVE = (
    "temperature_K",
    "conformer_temperature_K",
    "bond_hysteresis",
    "reservoir_weight_energy",
    "reservoir_weight_diversity",
    "reservoir_weight_response",
)
# Fractions that are meaningless above one: a quota that exceeds the whole, an overlap
# beyond a unit dot product, a relative tolerance that accepts anything. Positivity is
# already required above; this is only the upper end.
FRACTION_AT_MOST_ONE = ("direction_quota_fraction", "min_mode_overlap")
FRACTION_BELOW_ONE = ("response_descent_tol", "min_distance_ratio", "parity_tolerance")
# Numeric fields the bucket lists deliberately skip, because a dedicated branch below
# validates them against something the buckets cannot express.
_UNBUCKETED_NUMERIC = frozenset(
    {
        "total_charge",  # int | None, tied to charge_sensitive
        "multiplicity",  # int | None, tied to charge_sensitive
        "max_uncertainty_eV_A",  # float | None, positive only when set
    }
)


@dataclass(frozen=True)
class SearchConfig:
    seed: int = 42
    max_trials: int = 96
    max_depth: int = 1
    chemical_max_nodes: int = 16
    chemical_trials_per_node: int = 96
    conformer_max_per_node: int = 8
    conformer_trials_per_microstate: int = 32
    conformer_temperature_K: float = 300.0
    reservoir_weight_energy: float = 1.0
    reservoir_weight_diversity: float = 1.0
    reservoir_weight_response: float = 1.0
    automorphism_limit: int = 20000
    charge_sensitive: bool = False
    total_charge: int | None = None
    multiplicity: int | None = None
    max_directions: int = 24
    families: tuple[str, ...] = (
        "stretch",
        "compress",
        "kick",
        "bend",
        "torsion",
        "torsion_kick",
    )
    geometry_amplitudes_A: tuple[float, ...] = (0.1, 0.3, 0.6)
    kick_energies_eV: tuple[float, ...] = (0.1, 0.5, 1.0)
    pulse_forces_eV_A: tuple[float, ...] = (0.5, 1.0, 2.0)
    bend_amplitudes_rad: tuple[float, ...] = (0.09, 0.17, 0.35)
    torsion_amplitudes_rad: tuple[float, ...] = (0.52, 1.05, 2.09)
    # How far a rotor amplitude must sit from any symmetry copy of zero to be informative.
    torsion_symmetry_margin_rad: float = 0.2
    refinement_steps: int = 3
    timestep_fs: float = 0.25
    response_steps: int = 200
    pulse_steps: int = 20
    sample_interval: int = 5
    temperature_K: float = 0.0
    quench_steps: int = 500
    quench_fmax_eV_A: float = 0.03
    torsion_tolerance_rad: float = 0.0175  # one degree
    soft_mode_curvature_max_eV_rad2: float = 50.0  # above this, fmax already suffices
    soft_mode_curvature_floor_eV_rad2: float = 0.01  # below this, a free coordinate
    soft_mode_step_rad: float = 0.02
    soft_mode_scan_points: int = 12
    soft_polish_rounds: int = 6
    soft_polish_steps: int = 6
    minimum_check: str = "hessian"
    hessian_source: str = "auto"  # auto | analytic | fd
    analytic_hessian_check_directions: int = 3
    analytic_hessian_check: float = 1e-2  # relative, vs one finite difference
    minimum_check_step_A: float = 0.01
    minimum_check_eigenvalue_tol: float = 1e-3
    trivial_floor_fraction: float = 0.5
    internal_mode_recheck: bool = True
    internal_mode_recheck_step_scale: float = 2.0
    free_alignment_samples: int = 24
    response_band_enabled: bool = True
    response_band_A: tuple = (0.01, 0.1)  # A*sqrt(amu), mass-weighted
    response_samples: int = 5
    response_significance: float = 3.0  # |kappa| must exceed this many sigma
    response_descent_enabled: bool = False
    response_descent_max_steps: int = 40
    response_descent_step: float = 0.2
    response_descent_tol: float = 1e-3
    min_mode_max_steps: int = 200
    min_mode_refresh: int = 5
    min_mode_step_A: float = 0.05
    min_mode_fmax_eV_A: float = 0.01
    min_mode_backtracks: int = 8
    min_mode_overlap: float = 0.5
    saddle_search_enabled: bool = True
    saddle_search_per_run: int = 8
    saddle_search_seeds: int = 2
    saddle_search_max_rise_eV: float = 1.5
    closed_shell_only: bool = True
    irc_enabled: bool = True
    irc_step_A: float = 0.15
    irc_max_per_run: int = 16
    irc_relay_hops: int = 3
    boundary_continuation_enabled: bool = True
    boundary_continuation_step_A: float = 0.05
    boundary_continuation_max_steps: int = 8
    boundary_continuation_descend_saddles: bool = True
    boundary_continuation_seeds: int = 2
    boundary_continuation_per_run: int = 4
    tracked_coordinates: tuple = ()
    minimum_check_probes: int = 2
    minimum_check_displacement_A: float = 0.05
    pair_cutoff_A: float = 4.0
    min_pair_hops: int = 4
    hbond_acceptor_A: float = 2.5
    hbond_angle_deg: float = 120.0
    direction_quota_fraction: float = 0.5
    min_directions_per_family: int = 2
    parity_tolerance: float = 0.01
    bond_scale: float = 1.2
    bond_hysteresis: float = 0.08
    basin_rmsd_A: float = 0.15
    basin_energy_eV: float = 0.05
    min_distance_ratio: float = 0.45
    collateral_gate: bool = True
    collateral_bond_A: float = 0.05
    max_force_eV_A: float = 100.0
    max_energy_rise_eV_atom: float = 20.0
    max_nve_drift_eV_atom: float = 0.1
    max_uncertainty_eV_A: float | None = None
    active_atoms: tuple[int, ...] | None = None

    def __post_init__(self):
        if not isinstance(self.irc_enabled, bool):
            raise ValueError("irc_enabled must be a bool")
        integer_positive = INTEGER_POSITIVE
        integer_nonnegative = INTEGER_NONNEGATIVE
        for name in integer_positive + integer_nonnegative:
            value = getattr(self, name)
            lower = 1 if name in integer_positive else 0
            if type(value) is not int or value < lower:
                raise ValueError(f"{name} must be an integer >= {lower}")
        for name in (
            "collateral_gate",
            "response_band_enabled",
            "internal_mode_recheck",
            "response_descent_enabled",
            "saddle_search_enabled",
            "closed_shell_only",
            "boundary_continuation_enabled",
            "boundary_continuation_descend_saddles",
        ):
            if not isinstance(getattr(self, name), bool):
                raise ValueError(f"{name} must be a bool")
        if (
            not 0
            < self.soft_mode_curvature_floor_eV_rad2
            < self.soft_mode_curvature_max_eV_rad2
        ):
            raise ValueError("soft-mode curvature floor must be positive and below the maximum")
        for name in FLOAT_POSITIVE:
            value = getattr(self, name)
            if isinstance(value, bool) or not math.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be finite and positive")
        for name in FRACTION_AT_MOST_ONE:
            if getattr(self, name) > 1:
                raise ValueError(f"{name} must be at most 1")
        for name in FRACTION_BELOW_ONE:
            if getattr(self, name) >= 1:
                raise ValueError(f"{name} must be below 1")
        for name in FLOAT_NONNEGATIVE:
            value = getattr(self, name)
            if isinstance(value, bool) or not math.isfinite(value) or value < 0:
                raise ValueError(f"{name} must be finite and nonnegative")
        if (
            sum(
                (
                    self.reservoir_weight_energy,
                    self.reservoir_weight_diversity,
                    self.reservoir_weight_response,
                )
            )
            <= 0
        ):
            raise ValueError("at least one reservoir weight must be positive")
        if not isinstance(self.charge_sensitive, bool):
            raise ValueError("charge_sensitive must be a bool")
        if not self.charge_sensitive and self.total_charge not in (None, 0):
            raise ValueError(
                "charge_sensitive is false, so total_charge cannot be set; "
                "the potential does not respond to it"
            )
        if self.minimum_check not in ("hessian", "probe", "none"):
            raise ValueError("minimum_check must be 'hessian', 'probe' or 'none'")
        if self.hessian_source not in ("auto", "analytic", "fd"):
            raise ValueError("hessian_source must be 'auto', 'analytic' or 'fd'")
        if self.response_samples < 3:
            raise ValueError("response_samples must be at least 3 to fit two coefficients")
        if len(tuple(self.response_band_A)) != 2:
            raise ValueError("response_band_A must be a (low, high) pair")
        low, high = self.response_band_A
        if not 0 < low < high:
            raise ValueError("response_band_A must satisfy 0 < low < high")
        if self.bond_hysteresis >= 0.5:
            raise ValueError("bond_hysteresis must be below 0.5")
        if self.max_uncertainty_eV_A is not None:
            if not math.isfinite(self.max_uncertainty_eV_A) or self.max_uncertainty_eV_A <= 0:
                raise ValueError("max_uncertainty_eV_A must be finite and positive")
        if not self.families or len(set(self.families)) != len(self.families):
            raise ValueError("families must be nonempty and unique")
        from .internal import atom_index
        from .perturbations import FAMILIES

        unknown = set(self.families) - set(FAMILIES)
        if unknown:
            raise ValueError(f"unknown perturbation family: {sorted(unknown)}")
        if "pulse" in self.families and self.pulse_steps < 1:
            raise ValueError("pulse requires pulse_steps >= 1")
        for name in (
            "geometry_amplitudes_A",
            "kick_energies_eV",
            "pulse_forces_eV_A",
            "bend_amplitudes_rad",
            "torsion_amplitudes_rad",
        ):
            values = tuple(getattr(self, name))
            if not values or any(
                isinstance(v, bool) or not math.isfinite(v) or v <= 0 for v in values
            ):
                raise ValueError(f"{name} must contain finite positive values")
            if list(values) != sorted(set(values)):
                raise ValueError(f"{name} must be strictly increasing")
            object.__setattr__(self, name, values)
        if self.active_atoms is not None:
            index = tuple(self.active_atoms)
            if len(index) < 2 or len(set(index)) != len(index):
                raise ValueError("active_atoms must list at least two distinct indices")
            if any(type(i) is not int or i < 0 for i in index):
                raise ValueError("active_atoms must be nonnegative integers")
            object.__setattr__(self, "active_atoms", tuple(sorted(index)))
        object.__setattr__(self, "families", tuple(self.families))
        tracked = tuple(dict(spec) for spec in self.tracked_coordinates)
        seen = set()
        for spec in tracked:
            if set(spec) != {"name", "kind", "indices"}:
                raise ValueError(
                    "tracked_coordinates entries need exactly name, kind and indices"
                )
            if not isinstance(spec["name"], str) or not spec["name"]:
                raise ValueError("tracked coordinate name must be a nonempty string")
            if spec["name"] in seen:
                raise ValueError(f"duplicate tracked coordinate name: {spec['name']}")
            seen.add(spec["name"])
            widths = {"bond": 2, "angle": 3, "dihedral": 4, "bond_difference": 4}
            if spec["kind"] not in widths:
                raise ValueError(f"unknown tracked coordinate kind: {spec['kind']}")
            # Same validator the coordinate functions use, so a spec that is accepted
            # here cannot be rejected per frame. No atom count is available yet -- the
            # config is built before any molecule -- so only sign and integrality are
            # settled now; the upper bound is checked against the geometry.
            indices = tuple(atom_index(i) for i in spec["indices"])
            if len(indices) != widths[spec["kind"]]:
                raise ValueError(
                    f"{spec['kind']} needs {widths[spec['kind']]} nonnegative atom indices"
                )
            # A bond difference shares its middle atom by construction, so distinctness is
            # required per bond rather than across all four.
            if spec["kind"] == "bond_difference":
                if indices[0] == indices[1] or indices[2] == indices[3]:
                    raise ValueError("each bond of a bond_difference needs two distinct atoms")
            elif len(set(indices)) != len(indices):
                raise ValueError(f"{spec['kind']} needs distinct atom indices")
            spec["indices"] = indices
        object.__setattr__(self, "tracked_coordinates", tracked)

    def to_dict(self):
        return asdict(self)

    @classmethod
    def from_dict(cls, data):
        unknown = set(data) - {f.name for f in fields(cls)}
        if unknown:
            raise ValueError(f"Unknown config keys: {sorted(unknown)}")
        return cls(**data)


RATIONALE = {
    # Why each non-obvious default is the value it is. Every entry is the record of a
    # measurement or a failure, not a preference, so retiring a default means retiring the
    # entry that justified it. Kept beside the dataclass rather than inside it so the field
    # table stays readable; test_validation pins every key to a real field.
    "chemical_max_nodes": (
        "Two budgets, deliberately not one. A ninth conformer of a known substance is not a "
        "chemical discovery and must not consume a chemical node slot; but conformers stay "
        "searchable, because reaction accessibility is conditional on conformer."
    ),
    "charge_sensitive": (
        "Charge enters chemical identity only when the potential responds to it. MACE-OFF "
        "returns identical energies for total charge 0, -1 and +1."
    ),
    "families": (
        "The default basis spans stretch, bend and torsion. A basis of pair directions alone "
        "cannot reach a torsional basin at any amplitude."
    ),
    "torsion_tolerance_rad": (
        "A force tolerance says how small the force is, not how far the geometry is from "
        "stationary: |dq| ~ |g_q| / k_q. Measured on ethanol, the hydroxyl torsion is 0.21 "
        "eV/rad^2 against 33-53 eV/A^2 for the bond stretches, so one tolerance cannot serve "
        "both. Soft coordinates get their own criterion instead of dragging the global one "
        "down to a precision the stiff coordinates do not need."
    ),
    "soft_mode_scan_points": (
        "A torsion whose local curvature is below the floor but whose gradient is not zero is "
        "not flat -- the Newton scale is simply unusable there. Its own path is an exact "
        "rigid rotation, so the symmetry-reduced period is scanned instead of trusting a "
        "local model that has just been shown not to hold. Below three points no scan is done "
        "and the coordinate is polished by capped steps alone, which is cheaper and less "
        "reliable."
    ),
    "soft_polish_rounds": (
        "Five rounds, measured. A polish round that moves the geometry invalidates the force "
        "tolerance it was standing on, so the quench re-minimizes and polishes again until "
        "neither moves; on the 3-oxobutanal product the methyl torsion needed five of those "
        "alternations to settle -- its curvature rose from 0.0041 eV/rad^2 at the inflection "
        "to 0.107 at the minimum, so the first steps are taken on a scale 25 times too small. "
        "Rounds cost nothing where nothing moves: the source structure exits at the first."
    ),
    "minimum_check": (
        '"hessian" is decisive; "probe" is the cheap fallback. The Hessian itself has two '
        "sources: a model whose forces are an autograd gradient already carries the second "
        "derivative, so rebuilding it from 6N force differences is numerical differentiation "
        'of an analytic derivative. "auto" prefers analytic, verifies it against one finite '
        "difference, and records which path was taken."
    ),
    "analytic_hessian_check_directions": (
        "Several seeded random probe directions, not one fixed pattern: a structured probe "
        "can be blind to a structured layout error, a random one cannot except on a "
        "measure-zero set. Two force evaluations each."
    ),
    "trivial_floor_fraction": (
        "The six rigid-body modes have exactly zero curvature on any invariant potential, so "
        "their measured size is the error floor of the whole calculation and needs no "
        "reference data. Measured at 3e-6 to 3e-5 on MACE-OFF24, against a tolerance of 1e-3."
    ),
    "internal_mode_recheck": (
        "The rigid-body floor is a necessary condition only -- an error confined to the "
        "internal subspace leaves it at exactly zero -- so the eigenvalues that decide the "
        "count are re-measured by a direct HVP at a different step. Different on purpose: at "
        "the same step the central difference is the same computation that built the matrix."
    ),
    "free_alignment_samples": (
        "A torsion the polish labelled free connects the structures either side of it by a "
        "barrierless path, so they are one microstate and the comparison has to be quotiented "
        "by that coordinate; otherwise a free rotor mints unbounded spurious conformers. Pure "
        "geometry, no force evaluations: this is only how finely the rotation is scanned."
    ),
    "response_band_enabled": (
        "A symmetric energy difference at amplitude a is the curvature of the potential "
        "low-passed at wavenumber ~2/a, not an approximation that improves as a shrinks. So "
        "the band is a resolution choice and belongs in the manifest. The thermal amplitude "
        "sqrt(kT/kappa) is the physical anchor: about 0.3 A*sqrt(amu) for a soft mode at 300 "
        "K."
    ),
    "response_descent_enabled": (
        "Sphere descent on the response curvature: two force evaluations per step, the price "
        "of one Hessian-vector product, and no spectrum required. Reduces to dimer rotation "
        "as the amplitude goes to zero. The a > 0 branch as an independent reading rather "
        "than a refinement. It seeds the descent from a random direction, so where it "
        "converges on the Hessian's eigenvector that agreement is corroboration -- seeding it "
        "with that eigenvector, which is what response_curvature does, cannot corroborate "
        "anything. Off by default: it costs up to 2 * response_descent_max_steps evaluations "
        "per confirmed first-order saddle, and every result recorded so far was produced "
        "without it."
    ),
    "min_mode_max_steps": (
        "Min-mode following: invert the force along the lowest mode and minimize the rest. "
        "Finding the unstable direction is not finding the saddle; the centre has to move "
        "too."
    ),
    "min_mode_backtracks": (
        "A step that overshoots the energy ceiling is halved, not fatal: one overshoot is a "
        "step-length problem and must not be reported as a domain problem."
    ),
    "min_mode_overlap": (
        "Following a target coordinate is only meaningful while some eigenvector still "
        "resembles it. Below this mass-weighted overlap the walk has changed subject."
    ),
    "saddle_search_enabled": (
        "Climbing to a saddle from a probe frame, once a scan has already shown the reaction "
        "happens. Separate from the trial budget: it spends force evaluations, not trials, "
        "and the seed comes from work already done."
    ),
    "closed_shell_only": (
        "Refuse endpoints with no neutral closed-shell valence assignment. MACE-OFF is "
        "trained on neutral closed-shell molecules and ignores total charge outright (section "
        "2.5), so a radical pair reached by an over-driven probe is not a product it can "
        "describe."
    ),
    "irc_enabled": (
        "Following a first-order saddle's unstable mode downhill on both sides names the two "
        "minima it connects, which is what turns an observed response into a channel with "
        "evidence at both ends. Costs two quenches and two minimum checks per candidate."
    ),
    "irc_relay_hops": (
        "Where a descent stops on a further saddle, relay from that saddle instead of pushing "
        "harder from the first one: a larger displacement can cross several dividing surfaces "
        "and name a basin the saddle does not connect."
    ),
    "boundary_continuation_enabled": (
        "Continuing the delivered probe past the frame where the topology switched. A saddle "
        "and its two-sided descent establish that two minima are joined; they do not "
        "establish that the probe that was delivered reaches the product, and those are "
        "different claims. The continuation answers the second one by staying on the probe's "
        "own tangent and quenching each rung independently: a rung that returns to the source "
        "says the crossing had not yet entered the product's basin, a rung that reaches "
        "another chemical state is a directed observation, and a rung that stops on a saddle "
        "again says the delivered path runs along the dividing surface."
    ),
    "boundary_continuation_step_A": (
        "A cap on the largest single-atom displacement per rung, not an RMS. Small on "
        "purpose: a long step can cross several dividing surfaces and the rung's verdict "
        "would then name a basin the probe never passed through."
    ),
    "boundary_continuation_descend_saddles": (
        "A rung whose quench stops on a first-order saddle has not necessarily failed to "
        "arrive: a saddle lies on the boundary of two basins, so if both of its descents name "
        "the same chemical state then every neighbourhood of that rung except the ridge "
        "belongs to that state. Measured on 3-oxobutanal round four, every rung of every "
        "continuation stopped on a rotor saddle 11 meV above the product minimum at -65 to "
        "-79 cm^-1 while the proton had crossed. Costs two quenches per stalled rung and is "
        "reported as a chain."
    ),
    "tracked_coordinates": (
        "Named internal coordinates written into every path record. Which coordinate measures "
        "reaction progress is chemistry, so it is declared rather than guessed. Each entry is "
        '{"name", "kind", "indices"} with kind in bond / angle / dihedral / bond_difference.'
    ),
    "min_pair_hops": (
        "A pair-approach probe between atoms two or three bonds apart is a bend or a torsion "
        "in disguise: the covalent skeleton already holds them, so the continuation must drag "
        "the connecting bonds and the collateral gate refuses the delivery. Measured on "
        "3-oxobutanal, every refusal was three hops or fewer and they were 45% of the budget."
    ),
    "hbond_acceptor_A": (
        "An X-H bond whose hydrogen already has an acceptor within reach is marked unusual by "
        "the reactant's own geometry, so it goes in the front band. This is an ordering, "
        "never a filter, and it reads nothing about the product."
    ),
    "direction_quota_fraction": (
        "The per-family direction quota grows with the candidate count. A constant budget "
        "cannot hold a constant chance of sampling any given direction while candidates grow "
        "with the molecule; where the global cap still binds, the manifest says so."
    ),
    "parity_tolerance": (
        "Below this the configuration at a stereocentre or a locked bond is recorded as "
        '"planar" instead of a sign. Both discriminants are dimensionless in [-1, 1], so one '
        "number serves both. Calibrated from the slope near the crossing, d|cos|/dangle is "
        "about 0.017 per degree, so 0.01 is roughly 0.6 degrees of angular slack -- the scale "
        "a quench on a soft torsion actually leaves. A converged minimum never sits there: 90 "
        "degrees about a locked double bond is the torsional barrier top, not a basin."
    ),
    "collateral_gate": (
        "Collateral gate. A probe may bend and twist what it couples to, but it must not "
        "stretch a covalent bond, drive atoms into each other, or tear the molecule apart."
    ),
    "active_atoms": (
        "None means every atom. A subset scopes graphs, RMSD and basin identity to a reactive "
        "region so an environment cannot masquerade as a new basin."
    ),
}
