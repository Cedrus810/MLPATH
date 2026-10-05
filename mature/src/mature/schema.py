"""Strict data objects for the maturation layer (PLAN_GEN_COMPAS_MATURATION
2026-10-03, sections 3.3 and M1).

Every object round-trips through JSON and validates on construction: unknown
keys, missing identity, NaN/inf and malformed hashes are refused, never
repaired. Two rules come straight from the plan's fail-closed stance:

- A region whose A and B sides share one graph is refused at definition time
  (``unsupported_same_graph``): first-entry shooting cannot tell the sides
  apart, and the M3 rule says report it, not quietly fall back to a fixed
  t_free verdict.
- ``p_B_given_AB`` is ``None`` while ``A + B == 0``: a conditional ratio
  without a denominator is not zero, it is unasked.

Units are pinned per protocol (eV / Angstrom / fs / K) and validated to be
exactly that set, so a record can never be read in another unit system.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass, field
from datetime import datetime, timezone

SCHEMA_VERSION = 1

#: ``z`` for a two-sided 95% Wilson interval, the same constant
#: ``prrs.runner.wilson_interval`` uses, so both layers quote identical bounds.
WILSON_Z = 1.959963985

PROPOSAL_KINDS = (
    "observed_frame",  # a real frame of an already-recorded path
    "path_interpolation",  # linear interpolation between parent frames
    "endpoint_idpp",  # IDPP between endpoints; never a real MD frame
    "generated",  # only legal once the B stage exists (plan section 6)
)
BATCH_ROLES = ("zero_measure", "explore", "confirm", "sidecheck")
SHOT_OUTCOMES = (
    "A",
    "B",
    "other",  # entered neither side's graph (third basin) or unidentified
    "unabsorbed",  # ran to the cap without k consecutive steps on one graph
    "gate_rejected",  # refused by a pre-run gate (drift, budget, domain)
    "failed",  # the trial itself failed (quench, calculator, ...)
    "crashed",  # the process died before the trial recorded a verdict
)
DIAGNOSIS_LABELS = (
    "transition_region_compatible",
    "A_side",
    "B_side",
    "unresolved",
    "gate_rejected",
)
# Deliberately no "verified"/"confirmed_ts" value: q ~ 0.5 is not a saddle
# proof (feasibility assessment, corrected judgement 3). Verification lives
# in the refinement/validation handoff, not in this layer's vocabulary.
RECOMMENDED_USES = ("neb_refinement", "independent_confirmation", "no_action")
EXPECTED_UNITS = {
    "energy": "eV",
    "length": "Angstrom",
    "time": "fs",
    "temperature": "K",
}


class SchemaError(ValueError):
    """A record violated the maturation schema; the message says which field."""


# ---- validation helpers -----------------------------------------------------


def _keys(data, required, optional, what):
    if not isinstance(data, dict):
        raise SchemaError(f"{what}: expected an object, got {type(data).__name__}")
    missing = [k for k in required if k not in data]
    if missing:
        raise SchemaError(f"{what}: missing required fields {sorted(missing)}")
    unknown = set(data) - set(required) - set(optional)
    if unknown:
        raise SchemaError(f"{what}: unknown fields {sorted(unknown)}")


def _finite(value, what):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise SchemaError(f"{what}: expected a number, got {value!r}")
    if not math.isfinite(value):
        raise SchemaError(f"{what}: not finite ({value!r}); NaN/inf are never recorded")
    return float(value)


def _positive(value, what):
    value = _finite(value, what)
    if value <= 0:
        raise SchemaError(f"{what}: must be > 0, got {value}")
    return value


def _ratio(value, what):
    """A number in [0, 1], or None where the quantity is legitimately absent."""
    if value is None:
        return None
    value = _finite(value, what)
    if not 0.0 <= value <= 1.0:
        raise SchemaError(f"{what}: must be within [0, 1], got {value}")
    return value


def _int_at_least(value, what, minimum):
    if isinstance(value, bool) or not isinstance(value, int):
        raise SchemaError(f"{what}: expected an integer, got {value!r}")
    if value < minimum:
        raise SchemaError(f"{what}: must be >= {minimum}, got {value}")
    return value


def _sha256_hex(value, what):
    if not isinstance(value, str) or len(value) != 64:
        raise SchemaError(f"{what}: expected a 64-hex sha256, got {value!r}")
    try:
        int(value, 16)
    except ValueError:
        raise SchemaError(f"{what}: expected a 64-hex sha256, got {value!r}") from None
    return value


def _identifier(value, what):
    if not isinstance(value, str) or not value.strip():
        raise SchemaError(f"{what}: expected a non-empty identifier, got {value!r}")
    return value


def _enum(value, choices, what):
    if value not in choices:
        raise SchemaError(f"{what}: must be one of {list(choices)}, got {value!r}")
    return value


def _bond_pairs(value, what):
    """A sorted list of unique [i, j] atom-index pairs with i < j."""
    if not isinstance(value, list):
        raise SchemaError(f"{what}: expected a list of [i, j] pairs, got {value!r}")
    pairs = []
    for pair in value:
        if (
            not isinstance(pair, (list, tuple))
            or len(pair) != 2
            or any(isinstance(x, bool) or not isinstance(x, int) for x in pair)
        ):
            raise SchemaError(f"{what}: each entry must be [i, j] integers, got {pair!r}")
        i, j = pair
        if i < 0 or j < 0 or i == j:
            raise SchemaError(f"{what}: pair {pair!r} is not a valid bond")
        pairs.append((min(i, j), max(i, j)))
    if pairs != sorted(set(pairs)):
        raise SchemaError(f"{what}: pairs must be sorted and unique, got {value!r}")
    return pairs


def canonical_json(obj) -> str:
    """The one serialisation hashes are computed over (no NaN, sorted keys)."""
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), allow_nan=False)


def sha256_object(obj) -> str:
    return hashlib.sha256(canonical_json(obj).encode("utf-8")).hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def wilson_interval(successes, total, z=WILSON_Z):
    """Two-sided 95% Wilson score interval; [] when the denominator is 0.

    The interval travels with the point estimate everywhere in this layer
    (plan M3: 16/16 is not "close to 1" -- the lower bound is the claim).
    """
    if total <= 0:
        return []
    p = successes / total
    denom = 1 + z * z / total
    centre = (p + z * z / (2 * total)) / denom
    half = z * math.sqrt(p * (1 - p) / total + z * z / (4 * total * total)) / denom
    return [max(0.0, centre - half), min(1.0, centre + half)]


# ---- protocol ---------------------------------------------------------------


@dataclass(frozen=True)
class MicrostateRef:
    """One side of an A/B region, identified by its mapped microstate.

    Identity is the run's own microstate id (``c0000/m0001`` style) plus the
    full bond graph of its structure -- NOT the symmetry-quotiented
    ``chemical_key``, which for P1 is identical on both sides and therefore
    cannot define a region (plan M0: "识别 A/B 用映射后的微态与键集合").
    """

    node_id: str
    microstate_id: str
    chemical_key: str
    structure_path: str  # relative to the source run directory
    structure_sha256: str
    edges: tuple[tuple[int, int], ...]
    mapped_broken: tuple[tuple[int, int], ...] | None = None
    mapped_formed: tuple[tuple[int, int], ...] | None = None

    def __post_init__(self):
        _identifier(self.node_id, "MicrostateRef.node_id")
        _identifier(self.microstate_id, "MicrostateRef.microstate_id")
        _identifier(self.chemical_key, "MicrostateRef.chemical_key")
        _identifier(self.structure_path, "MicrostateRef.structure_path")
        _sha256_hex(self.structure_sha256, "MicrostateRef.structure_sha256")
        object.__setattr__(
            self, "edges", tuple(_bond_pairs(list(self.edges), "MicrostateRef.edges"))
        )
        if not self.edges:
            raise SchemaError("MicrostateRef.edges: a graph with no bonds is not a region side")
        for name in ("mapped_broken", "mapped_formed"):
            value = getattr(self, name)
            if value is not None:
                object.__setattr__(
                    self, name, tuple(_bond_pairs(list(value), f"MicrostateRef.{name}"))
                )

    @property
    def graph_sha256(self) -> str:
        return sha256_object({"edges": [list(pair) for pair in self.edges]})

    def to_dict(self) -> dict:
        return {
            "node_id": self.node_id,
            "microstate_id": self.microstate_id,
            "chemical_key": self.chemical_key,
            "structure_path": self.structure_path,
            "structure_sha256": self.structure_sha256,
            "edges": [list(pair) for pair in self.edges],
            "mapped_broken": (
                [list(pair) for pair in self.mapped_broken] if self.mapped_broken else None
            ),
            "mapped_formed": (
                [list(pair) for pair in self.mapped_formed] if self.mapped_formed else None
            ),
        }

    @classmethod
    def from_dict(cls, data) -> "MicrostateRef":
        _keys(
            data,
            [
                "node_id",
                "microstate_id",
                "chemical_key",
                "structure_path",
                "structure_sha256",
                "edges",
            ],
            ["mapped_broken", "mapped_formed"],
            "MicrostateRef",
        )
        return cls(
            node_id=data["node_id"],
            microstate_id=data["microstate_id"],
            chemical_key=data["chemical_key"],
            structure_path=data["structure_path"],
            structure_sha256=data["structure_sha256"],
            edges=[tuple(pair) for pair in data["edges"]],
            mapped_broken=(
                [tuple(pair) for pair in data["mapped_broken"]]
                if data.get("mapped_broken")
                else None
            ),
            mapped_formed=(
                [tuple(pair) for pair in data["mapped_formed"]]
                if data.get("mapped_formed")
                else None
            ),
        )


@dataclass(frozen=True)
class RegionDefinition:
    """Exactly two sides, A and B, whose graphs differ.

    ``same graph`` is refused here rather than discovered at M3: the region
    is what makes first-entry shooting meaningful, so a region that cannot
    support it is not a region (plan M3: "图相同明确报 unsupported_same_graph").
    """

    sides: dict  # {"A": MicrostateRef, "B": MicrostateRef}
    method: str = "mapped_microstate"
    version: str = "1"

    def __post_init__(self):
        if set(self.sides) != {"A", "B"}:
            raise SchemaError(
                f"RegionDefinition.sides: must have exactly the keys ['A', 'B'], "
                f"got {sorted(self.sides)}"
            )
        object.__setattr__(self, "sides", {k: _as_ref(v) for k, v in self.sides.items()})
        if self.sides["A"].graph_sha256 == self.sides["B"].graph_sha256:
            raise SchemaError(
                "RegionDefinition: A and B share one graph (unsupported_same_graph); "
                "first-entry shooting cannot tell the sides apart, so this region "
                "must be refused at definition time, not at shooting time"
            )
        _enum(self.method, ("mapped_microstate",), "RegionDefinition.method")
        _identifier(self.version, "RegionDefinition.version")

    @property
    def graph_sha256s(self) -> dict:
        return {side: ref.graph_sha256 for side, ref in sorted(self.sides.items())}

    def to_dict(self) -> dict:
        return {
            "sides": {side: ref.to_dict() for side, ref in sorted(self.sides.items())},
            "method": self.method,
            "version": self.version,
        }

    @classmethod
    def from_dict(cls, data) -> "RegionDefinition":
        _keys(data, ["sides"], ["method", "version"], "RegionDefinition")
        return cls(
            sides={side: MicrostateRef.from_dict(ref) for side, ref in data["sides"].items()},
            method=data.get("method", "mapped_microstate"),
            version=data.get("version", "1"),
        )


def _as_ref(value) -> MicrostateRef:
    return value if isinstance(value, MicrostateRef) else MicrostateRef.from_dict(value)


@dataclass(frozen=True)
class SourceRef:
    """Where the candidates come from, with the hashes that pin them.

    ``implementation_sha256`` is the SOURCE RUN's fingerprint -- provenance,
    not a claim that this layer runs under it. The fingerprint this batch
    itself runs under is recorded in the manifest by the io layer.
    """

    run_dir: str
    implementation_sha256: str
    model_factory: str
    model_kwargs: dict
    input_structure_sha256: str | None = None
    model_checkpoint_path: str | None = None
    model_checkpoint_sha256: str | None = None
    pars_selection_dir: str | None = None
    pars_selection_sha256: str | None = None
    files: dict = field(default_factory=dict)  # relative path -> sha256

    def __post_init__(self):
        _identifier(self.run_dir, "SourceRef.run_dir")
        _sha256_hex(self.implementation_sha256, "SourceRef.implementation_sha256")
        _identifier(self.model_factory, "SourceRef.model_factory")
        if not isinstance(self.model_kwargs, dict):
            raise SchemaError("SourceRef.model_kwargs: expected an object")
        if self.input_structure_sha256 is not None:
            _sha256_hex(self.input_structure_sha256, "SourceRef.input_structure_sha256")
        if self.model_checkpoint_path is not None:
            _identifier(self.model_checkpoint_path, "SourceRef.model_checkpoint_path")
        if self.model_checkpoint_sha256 is not None:
            _sha256_hex(self.model_checkpoint_sha256, "SourceRef.model_checkpoint_sha256")
        if self.pars_selection_dir is not None:
            _identifier(self.pars_selection_dir, "SourceRef.pars_selection_dir")
        if self.pars_selection_sha256 is not None:
            _sha256_hex(self.pars_selection_sha256, "SourceRef.pars_selection_sha256")
        if not isinstance(self.files, dict) or not self.files:
            raise SchemaError(
                "SourceRef.files: at least the consumed source files must be "
                "recorded as {relative path: sha256}; nothing else pins the inputs"
            )
        object.__setattr__(
            self,
            "files",
            {str(k): _sha256_hex(v, f"SourceRef.files[{k}]") for k, v in self.files.items()},
        )

    def to_dict(self) -> dict:
        return {
            "run_dir": self.run_dir,
            "implementation_sha256": self.implementation_sha256,
            "model_factory": self.model_factory,
            "model_kwargs": dict(self.model_kwargs),
            "input_structure_sha256": self.input_structure_sha256,
            "model_checkpoint_path": self.model_checkpoint_path,
            "model_checkpoint_sha256": self.model_checkpoint_sha256,
            "pars_selection_dir": self.pars_selection_dir,
            "pars_selection_sha256": self.pars_selection_sha256,
            "files": dict(self.files),
        }

    @classmethod
    def from_dict(cls, data) -> "SourceRef":
        _keys(
            data,
            ["run_dir", "implementation_sha256", "model_factory", "model_kwargs", "files"],
            [
                "input_structure_sha256",
                "model_checkpoint_path",
                "model_checkpoint_sha256",
                "pars_selection_dir",
                "pars_selection_sha256",
            ],
            "SourceRef",
        )
        return cls(
            run_dir=data["run_dir"],
            implementation_sha256=data["implementation_sha256"],
            model_factory=data["model_factory"],
            model_kwargs=dict(data["model_kwargs"]),
            input_structure_sha256=data.get("input_structure_sha256"),
            model_checkpoint_path=data.get("model_checkpoint_path"),
            model_checkpoint_sha256=data.get("model_checkpoint_sha256"),
            pars_selection_dir=data.get("pars_selection_dir"),
            pars_selection_sha256=data.get("pars_selection_sha256"),
            files=dict(data["files"]),
        )


@dataclass(frozen=True)
class DynamicsSpec:
    """The dynamics the shots will run; units are the protocol's.

    v1 is NVE with Maxwell-Boltzmann momenta and no forcing: the pulse stage
    of ``run_trial`` is never entered for shots (``deliver=False``), and a
    protocol that claims pulse steps is recording an intent the code must
    not have -- so it is refused here (plan section 1: "不把位置 q 模型当成
    相空间 Doob 桥的 h"; feasibility: "施力段不作为无偏动力学").
    """

    temperature_K: float
    timestep_fs: float
    first_entry_k: int
    free_step_cap: int
    ensemble: str = "NVE"
    momenta: str = "maxwell-boltzmann"
    integrator: str = "velocity-verlet"
    pulse_steps: int = 0

    def __post_init__(self):
        _positive(self.temperature_K, "DynamicsSpec.temperature_K")
        _positive(self.timestep_fs, "DynamicsSpec.timestep_fs")
        _int_at_least(self.first_entry_k, "DynamicsSpec.first_entry_k", 1)
        _int_at_least(self.free_step_cap, "DynamicsSpec.free_step_cap", 1)
        _enum(self.ensemble, ("NVE",), "DynamicsSpec.ensemble")
        _enum(self.momenta, ("maxwell-boltzmann",), "DynamicsSpec.momenta")
        _enum(self.integrator, ("velocity-verlet",), "DynamicsSpec.integrator")
        if self.pulse_steps != 0:
            raise SchemaError(
                "DynamicsSpec.pulse_steps: v1 shots are delivered without forcing "
                "(deliver=False); a protocol with pulse steps is refused"
            )

    def to_dict(self) -> dict:
        return {
            "temperature_K": self.temperature_K,
            "timestep_fs": self.timestep_fs,
            "first_entry_k": self.first_entry_k,
            "free_step_cap": self.free_step_cap,
            "ensemble": self.ensemble,
            "momenta": self.momenta,
            "integrator": self.integrator,
            "pulse_steps": self.pulse_steps,
        }

    @classmethod
    def from_dict(cls, data) -> "DynamicsSpec":
        _keys(
            data,
            ["temperature_K", "timestep_fs", "first_entry_k", "free_step_cap"],
            ["ensemble", "momenta", "integrator", "pulse_steps"],
            "DynamicsSpec",
        )
        return cls(
            temperature_K=data["temperature_K"],
            timestep_fs=data["timestep_fs"],
            first_entry_k=data["first_entry_k"],
            free_step_cap=data["free_step_cap"],
            ensemble=data.get("ensemble", "NVE"),
            momenta=data.get("momenta", "maxwell-boltzmann"),
            integrator=data.get("integrator", "velocity-verlet"),
            pulse_steps=data.get("pulse_steps", 0),
        )


@dataclass(frozen=True)
class CandidateRule:
    """Which candidates a channel may contribute (plan M2 候选规则)."""

    max_paths_per_channel: int = 2
    inner_candidates_per_path: int = 5
    arc_fractions: tuple = (1 / 6, 2 / 6, 3 / 6, 4 / 6, 5 / 6)
    allow_interpolation: bool = True
    preserve_all_hydrogens: bool = True
    rigid_alignment: str = "isolated_rigid_invariant_only"

    def __post_init__(self):
        _int_at_least(self.max_paths_per_channel, "CandidateRule.max_paths_per_channel", 1)
        _int_at_least(
            self.inner_candidates_per_path, "CandidateRule.inner_candidates_per_path", 1
        )
        if not isinstance(self.arc_fractions, (list, tuple)) or not self.arc_fractions:
            raise SchemaError("CandidateRule.arc_fractions: expected a non-empty list")
        fractions = [_ratio(f, "CandidateRule.arc_fractions[i]") for f in self.arc_fractions]
        if any(f <= 0.0 or f >= 1.0 for f in fractions):
            raise SchemaError(
                "CandidateRule.arc_fractions: interior positions only -- endpoints are "
                "not internal transition frames and must not be proposed as such"
            )
        if fractions != sorted(fractions) or len(set(fractions)) != len(fractions):
            raise SchemaError("CandidateRule.arc_fractions: must be strictly increasing")
        object.__setattr__(self, "arc_fractions", tuple(fractions))
        if not self.preserve_all_hydrogens:
            raise SchemaError(
                "CandidateRule.preserve_all_hydrogens: v1 keeps every hydrogen "
                "(plan M2: 全 H 保留); refusing anything else"
            )
        _enum(
            self.rigid_alignment,
            ("isolated_rigid_invariant_only",),
            "CandidateRule.rigid_alignment",
        )

    def to_dict(self) -> dict:
        return {
            "max_paths_per_channel": self.max_paths_per_channel,
            "inner_candidates_per_path": self.inner_candidates_per_path,
            "arc_fractions": list(self.arc_fractions),
            "allow_interpolation": self.allow_interpolation,
            "preserve_all_hydrogens": self.preserve_all_hydrogens,
            "rigid_alignment": self.rigid_alignment,
        }

    @classmethod
    def from_dict(cls, data) -> "CandidateRule":
        _keys(
            data,
            [],
            [
                "max_paths_per_channel",
                "inner_candidates_per_path",
                "arc_fractions",
                "allow_interpolation",
                "preserve_all_hydrogens",
                "rigid_alignment",
            ],
            "CandidateRule",
        )
        return cls(
            max_paths_per_channel=data.get("max_paths_per_channel", 2),
            inner_candidates_per_path=data.get("inner_candidates_per_path", 5),
            arc_fractions=data.get("arc_fractions", [1 / 6, 2 / 6, 3 / 6, 4 / 6, 5 / 6]),
            allow_interpolation=data.get("allow_interpolation", True),
            preserve_all_hydrogens=data.get("preserve_all_hydrogens", True),
            rigid_alignment=data.get("rigid_alignment", "isolated_rigid_invariant_only"),
        )


@dataclass(frozen=True)
class ShootingRule:
    """Shot counts, gates and seeds, per batch role (plan M3/M5)."""

    batch_role: str
    seed_base: int
    shots_per_candidate: int | None = None  # zero_measure / explore
    confirm_shots: int | None = None  # confirm
    sidecheck_shots: int | None = None  # sidecheck
    selection_rule: str | None = None  # required for confirm/sidecheck
    diagnostic_band: tuple = (0.4, 0.6)
    other_fraction_max: float = 0.10
    entry_quench_mismatch_max: int = 0
    drift_tolerance_eV_atom: float | None = None

    def __post_init__(self):
        _enum(self.batch_role, BATCH_ROLES, "ShootingRule.batch_role")
        _int_at_least(self.seed_base, "ShootingRule.seed_base", 0)
        for name in ("shots_per_candidate", "confirm_shots", "sidecheck_shots"):
            value = getattr(self, name)
            if value is not None:
                object.__setattr__(self, name, _int_at_least(value, f"ShootingRule.{name}", 1))
        if self.batch_role in ("zero_measure", "explore"):
            if self.shots_per_candidate is None:
                raise SchemaError(
                    f"ShootingRule: batch_role {self.batch_role!r} needs shots_per_candidate"
                )
        elif self.batch_role == "confirm":
            if self.confirm_shots is None:
                raise SchemaError("ShootingRule: batch_role 'confirm' needs confirm_shots")
        elif self.batch_role == "sidecheck":
            if self.sidecheck_shots is None:
                raise SchemaError("ShootingRule: batch_role 'sidecheck' needs sidecheck_shots")
        if self.batch_role in ("confirm", "sidecheck"):
            if not self.selection_rule or not str(self.selection_rule).strip():
                raise SchemaError(
                    "ShootingRule: confirm/sidecheck candidates are chosen by a rule "
                    "frozen before the run; selection_rule is required and must be "
                    "non-empty"
                )
        if self.batch_role != "zero_measure" and self.drift_tolerance_eV_atom is None:
            raise SchemaError(
                "ShootingRule.drift_tolerance_eV_atom: only the zero-measurement batch "
                "may leave the drift tolerance open -- its purpose is to measure it; "
                "every other batch freezes a value before running"
            )
        if self.drift_tolerance_eV_atom is not None:
            object.__setattr__(
                self,
                "drift_tolerance_eV_atom",
                _positive(self.drift_tolerance_eV_atom, "ShootingRule.drift_tolerance_eV_atom"),
            )
        band = self.diagnostic_band
        if (
            not isinstance(band, (list, tuple))
            or len(band) != 2
            or not (0.0 <= band[0] < band[1] <= 1.0)
        ):
            raise SchemaError(
                f"ShootingRule.diagnostic_band: expected [low, high] with "
                f"0 <= low < high <= 1, got {band!r}"
            )
        object.__setattr__(self, "diagnostic_band", tuple(band))
        object.__setattr__(
            self,
            "other_fraction_max",
            _ratio(self.other_fraction_max, "ShootingRule.other_fraction_max"),
        )
        _int_at_least(
            self.entry_quench_mismatch_max, "ShootingRule.entry_quench_mismatch_max", 0
        )

    def shots_for_role(self) -> int:
        """The one shot count this batch role uses; None means 'not this role'."""
        return {
            "zero_measure": self.shots_per_candidate,
            "explore": self.shots_per_candidate,
            "confirm": self.confirm_shots,
            "sidecheck": self.sidecheck_shots,
        }[self.batch_role]

    def to_dict(self) -> dict:
        return {
            "batch_role": self.batch_role,
            "seed_base": self.seed_base,
            "shots_per_candidate": self.shots_per_candidate,
            "confirm_shots": self.confirm_shots,
            "sidecheck_shots": self.sidecheck_shots,
            "selection_rule": self.selection_rule,
            "diagnostic_band": list(self.diagnostic_band),
            "other_fraction_max": self.other_fraction_max,
            "entry_quench_mismatch_max": self.entry_quench_mismatch_max,
            "drift_tolerance_eV_atom": self.drift_tolerance_eV_atom,
        }

    @classmethod
    def from_dict(cls, data) -> "ShootingRule":
        _keys(
            data,
            ["batch_role", "seed_base"],
            [
                "shots_per_candidate",
                "confirm_shots",
                "sidecheck_shots",
                "selection_rule",
                "diagnostic_band",
                "other_fraction_max",
                "entry_quench_mismatch_max",
                "drift_tolerance_eV_atom",
            ],
            "ShootingRule",
        )
        return cls(
            batch_role=data["batch_role"],
            seed_base=data["seed_base"],
            shots_per_candidate=data.get("shots_per_candidate"),
            confirm_shots=data.get("confirm_shots"),
            sidecheck_shots=data.get("sidecheck_shots"),
            selection_rule=data.get("selection_rule"),
            diagnostic_band=data.get("diagnostic_band", [0.4, 0.6]),
            other_fraction_max=data.get("other_fraction_max", 0.10),
            entry_quench_mismatch_max=data.get("entry_quench_mismatch_max", 0),
            drift_tolerance_eV_atom=data.get("drift_tolerance_eV_atom"),
        )


@dataclass(frozen=True)
class Budget:
    """The accounting limits; at least one must be set (plan M4 续跑与成本)."""

    max_shots: int | None = None
    max_free_steps: int | None = None
    max_force_evaluations: int | None = None
    max_wall_minutes: float | None = None

    def __post_init__(self):
        for name in ("max_shots", "max_free_steps", "max_force_evaluations"):
            value = getattr(self, name)
            if value is not None:
                object.__setattr__(self, name, _int_at_least(value, f"Budget.{name}", 1))
        if self.max_wall_minutes is not None:
            object.__setattr__(
                self,
                "max_wall_minutes",
                _positive(self.max_wall_minutes, "Budget.max_wall_minutes"),
            )
        if all(
            getattr(self, name) is None
            for name in (
                "max_shots",
                "max_free_steps",
                "max_force_evaluations",
                "max_wall_minutes",
            )
        ):
            raise SchemaError(
                "Budget: at least one limit must be set; a batch with no accounting "
                "is not a batch this layer runs"
            )

    def to_dict(self) -> dict:
        return {
            "max_shots": self.max_shots,
            "max_free_steps": self.max_free_steps,
            "max_force_evaluations": self.max_force_evaluations,
            "max_wall_minutes": self.max_wall_minutes,
        }

    @classmethod
    def from_dict(cls, data) -> "Budget":
        _keys(
            data,
            [],
            ["max_shots", "max_free_steps", "max_force_evaluations", "max_wall_minutes"],
            "Budget",
        )
        return cls(
            max_shots=data.get("max_shots"),
            max_free_steps=data.get("max_free_steps"),
            max_force_evaluations=data.get("max_force_evaluations"),
            max_wall_minutes=data.get("max_wall_minutes"),
        )


@dataclass(frozen=True)
class MaturationProtocol:
    """The frozen contract of one maturation batch (plan sections 3.2-3.3).

    Frozen means frozen: the protocol's canonical-JSON hash is written into
    the manifest, the resume check compares it, and any field change -- model,
    region definition, budget, seed table -- produces a different hash and is
    refused as a resume (M1 acceptance: "换模型/区域定义不能续在原目录").
    """

    batch_role: str
    source: SourceRef
    regions: RegionDefinition
    dynamics: DynamicsSpec
    candidates: CandidateRule
    shooting: ShootingRule
    budget: Budget
    title: str = ""
    notes: str = ""
    created_utc: str = field(default_factory=utc_now)
    schema_version: int = SCHEMA_VERSION
    units: dict = field(default_factory=lambda: dict(EXPECTED_UNITS))

    def __post_init__(self):
        _enum(self.batch_role, BATCH_ROLES, "MaturationProtocol.batch_role")
        if not isinstance(self.source, SourceRef):
            raise SchemaError("MaturationProtocol.source: expected a SourceRef")
        if not isinstance(self.regions, RegionDefinition):
            raise SchemaError("MaturationProtocol.regions: expected a RegionDefinition")
        for name, kind in (
            ("dynamics", DynamicsSpec),
            ("candidates", CandidateRule),
            ("shooting", ShootingRule),
            ("budget", Budget),
        ):
            if not isinstance(getattr(self, name), kind):
                raise SchemaError(f"MaturationProtocol.{name}: expected a {kind.__name__}")
        if self.shooting.batch_role != self.batch_role:
            raise SchemaError(
                f"MaturationProtocol: batch_role {self.batch_role!r} does not match "
                f"ShootingRule.batch_role {self.shooting.batch_role!r}"
            )
        if self.schema_version != SCHEMA_VERSION:
            raise SchemaError(
                f"MaturationProtocol.schema_version: expected {SCHEMA_VERSION}, "
                f"got {self.schema_version}"
            )
        if self.units != EXPECTED_UNITS:
            raise SchemaError(
                f"MaturationProtocol.units: v1 records are eV / Angstrom / fs / K; "
                f"got {self.units}"
            )
        _identifier(self.created_utc, "MaturationProtocol.created_utc")

    @property
    def sha256(self) -> str:
        return sha256_object(self.to_dict())

    def to_dict(self) -> dict:
        return {
            "schema_version": self.schema_version,
            "batch_role": self.batch_role,
            "title": self.title,
            "created_utc": self.created_utc,
            "notes": self.notes,
            "source": self.source.to_dict(),
            "regions": self.regions.to_dict(),
            "dynamics": self.dynamics.to_dict(),
            "candidates": self.candidates.to_dict(),
            "shooting": self.shooting.to_dict(),
            "budget": self.budget.to_dict(),
            "units": dict(self.units),
        }

    @classmethod
    def from_dict(cls, data) -> "MaturationProtocol":
        _keys(
            data,
            [
                "batch_role",
                "source",
                "regions",
                "dynamics",
                "candidates",
                "shooting",
                "budget",
            ],
            ["title", "notes", "schema_version", "units", "created_utc"],
            "MaturationProtocol",
        )
        return cls(
            batch_role=data["batch_role"],
            source=SourceRef.from_dict(data["source"]),
            regions=RegionDefinition.from_dict(data["regions"]),
            dynamics=DynamicsSpec.from_dict(data["dynamics"]),
            candidates=CandidateRule.from_dict(data["candidates"]),
            shooting=ShootingRule.from_dict(data["shooting"]),
            budget=Budget.from_dict(data["budget"]),
            title=data.get("title", ""),
            notes=data.get("notes", ""),
            created_utc=data.get("created_utc") or utc_now(),
            schema_version=data.get("schema_version", SCHEMA_VERSION),
            units=dict(data.get("units", EXPECTED_UNITS)),
        )


# ---- records ----------------------------------------------------------------


@dataclass(frozen=True)
class PathProposal:
    """One candidate geometry taken from an already-discovered path (plan 3.3).

    ``proposal_kind`` is the honest provenance label: an interpolated or IDPP
    geometry is a legal PROPOSAL but is never a real unbiased MD frame, and
    the record must not blur that line (plan 3.3: "插值候选是合法提议几何,
    不能标成真实无偏 MD 帧").
    """

    proposal_id: str
    batch_role: str
    channel_id: str
    geometry_file: str  # relative to the maturation output directory
    geometry_sha256: str
    n_atoms: int
    event_key: str | None = None
    tube_id: str | None = None
    source_candidate_id: str | None = None
    source_run_dir: str | None = None
    source_attempt_id: str | None = None
    parent_frame_start: int | None = None
    parent_frame_stop: int | None = None
    parent_path_sha256: str | None = None
    proposal_kind: str = "observed_frame"
    arc_fraction: float | None = None
    frame_index: int | None = None
    interpolation: dict | None = None  # {"method": ..., "weights": [...]}
    region_graph_sha256s: dict = field(default_factory=dict)  # {"A": ..., "B": ...}

    def __post_init__(self):
        _identifier(self.proposal_id, "PathProposal.proposal_id")
        _enum(self.batch_role, BATCH_ROLES, "PathProposal.batch_role")
        _identifier(self.channel_id, "PathProposal.channel_id")
        _identifier(self.geometry_file, "PathProposal.geometry_file")
        _sha256_hex(self.geometry_sha256, "PathProposal.geometry_sha256")
        _int_at_least(self.n_atoms, "PathProposal.n_atoms", 1)
        _enum(self.proposal_kind, PROPOSAL_KINDS, "PathProposal.proposal_kind")
        for name in ("event_key", "tube_id", "source_candidate_id", "source_run_dir"):
            value = getattr(self, name)
            if value is not None:
                _identifier(value, f"PathProposal.{name}")
        for name in ("source_attempt_id",):
            value = getattr(self, name)
            if value is not None:
                _identifier(value, f"PathProposal.{name}")
        for name in (
            "parent_frame_start",
            "parent_frame_stop",
            "frame_index",
        ):
            value = getattr(self, name)
            if value is not None:
                _int_at_least(value, f"PathProposal.{name}", 0)
        if self.parent_path_sha256 is not None:
            _sha256_hex(self.parent_path_sha256, "PathProposal.parent_path_sha256")
        object.__setattr__(
            self, "arc_fraction", _ratio(self.arc_fraction, "PathProposal.arc_fraction")
        )
        if self.proposal_kind == "observed_frame":
            if self.frame_index is None:
                raise SchemaError(
                    "PathProposal: an observed_frame must name the frame it came from"
                )
            if self.interpolation is not None:
                raise SchemaError(
                    "PathProposal: an observed_frame cannot carry interpolation weights"
                )
        elif self.proposal_kind == "path_interpolation":
            if not self.allowable_interpolation():
                raise SchemaError(
                    "PathProposal: path_interpolation needs {'method', 'weights'} with "
                    "finite positive weights"
                )
        elif self.proposal_kind == "endpoint_idpp":
            method = (self.interpolation or {}).get("method")
            if method != "idpp":
                raise SchemaError(
                    "PathProposal: endpoint_idpp needs interpolation={'method': 'idpp', ...}"
                )
        if not self.region_graph_sha256s:
            raise SchemaError(
                "PathProposal.region_graph_sha256s: a proposal must name the region it "
                "was proposed against; identity comes from the region definition, not "
                "from geometric similarity guesses"
            )
        if set(self.region_graph_sha256s) != {"A", "B"}:
            raise SchemaError(
                f"PathProposal.region_graph_sha256s: expected keys ['A', 'B'], "
                f"got {sorted(self.region_graph_sha256s)}"
            )
        for side, digest in self.region_graph_sha256s.items():
            _sha256_hex(digest, f"PathProposal.region_graph_sha256s[{side}]")

    def allowable_interpolation(self) -> bool:
        data = self.interpolation
        if not isinstance(data, dict) or set(data) - {"method", "weights"}:
            return False
        method = data.get("method")
        if not isinstance(method, str) or not method.strip():
            return False
        weights = data.get("weights")
        if not isinstance(weights, list) or not weights:
            return False
        try:
            values = [_finite(w, "PathProposal.interpolation.weights[i]") for w in weights]
        except SchemaError:
            return False
        return all(w >= 0 for w in values) and sum(values) > 0

    def to_dict(self) -> dict:
        return {
            "proposal_id": self.proposal_id,
            "batch_role": self.batch_role,
            "channel_id": self.channel_id,
            "event_key": self.event_key,
            "tube_id": self.tube_id,
            "source_candidate_id": self.source_candidate_id,
            "source_run_dir": self.source_run_dir,
            "source_attempt_id": self.source_attempt_id,
            "parent_frame_start": self.parent_frame_start,
            "parent_frame_stop": self.parent_frame_stop,
            "parent_path_sha256": self.parent_path_sha256,
            "proposal_kind": self.proposal_kind,
            "arc_fraction": self.arc_fraction,
            "frame_index": self.frame_index,
            "interpolation": (
                {**self.interpolation, "weights": list(self.interpolation["weights"])}
                if self.interpolation
                else None
            ),
            "geometry_file": self.geometry_file,
            "geometry_sha256": self.geometry_sha256,
            "n_atoms": self.n_atoms,
            "region_graph_sha256s": dict(self.region_graph_sha256s),
        }

    @classmethod
    def from_dict(cls, data) -> "PathProposal":
        _keys(
            data,
            [
                "proposal_id",
                "batch_role",
                "channel_id",
                "geometry_file",
                "geometry_sha256",
                "n_atoms",
                "proposal_kind",
                "region_graph_sha256s",
            ],
            [
                "event_key",
                "tube_id",
                "source_candidate_id",
                "source_run_dir",
                "source_attempt_id",
                "parent_frame_start",
                "parent_frame_stop",
                "parent_path_sha256",
                "arc_fraction",
                "frame_index",
                "interpolation",
            ],
            "PathProposal",
        )
        return cls(
            proposal_id=data["proposal_id"],
            batch_role=data["batch_role"],
            channel_id=data["channel_id"],
            event_key=data.get("event_key"),
            tube_id=data.get("tube_id"),
            source_candidate_id=data.get("source_candidate_id"),
            source_run_dir=data.get("source_run_dir"),
            source_attempt_id=data.get("source_attempt_id"),
            parent_frame_start=data.get("parent_frame_start"),
            parent_frame_stop=data.get("parent_frame_stop"),
            parent_path_sha256=data.get("parent_path_sha256"),
            proposal_kind=data["proposal_kind"],
            arc_fraction=data.get("arc_fraction"),
            frame_index=data.get("frame_index"),
            interpolation=dict(data["interpolation"]) if data.get("interpolation") else None,
            geometry_file=data["geometry_file"],
            geometry_sha256=data["geometry_sha256"],
            n_atoms=data["n_atoms"],
            region_graph_sha256s=dict(data["region_graph_sha256s"]),
        )


@dataclass(frozen=True)
class ShotRecord:
    """One shot's full verdict trail (plan 3.3; the E3 statistics gate).

    Outcomes ``A``/``B`` name a side; everything else -- third basin, running
    to the cap, a gate rejection, a failed quench, a dead process -- keeps its
    own label and is excluded from the conditional denominator. Nothing here
    is ever folded into A or B to make a ratio look cleaner.
    """

    proposal_id: str
    shot_id: str
    batch_role: str
    seed: int
    protocol_sha256: str
    implementation_sha256: str  # the prrs fingerprint the shot actually ran under
    outcome: str
    trial_dir: str
    outcome_reason: str | None = None
    entry_side: str | None = None
    entry_step: int | None = None
    entry_time_fs: float | None = None
    landed_microstate: str | None = None
    entry_quench_consistent: bool | None = None
    max_nve_drift_eV_atom: float | None = None
    quench_converged: bool | None = None
    cost_force_evaluations: int | None = None
    cost_wall_s: float | None = None
    endpoint_sha256: str | None = None
    model_checkpoint_sha256: str | None = None

    def __post_init__(self):
        _identifier(self.proposal_id, "ShotRecord.proposal_id")
        _identifier(self.shot_id, "ShotRecord.shot_id")
        _enum(self.batch_role, BATCH_ROLES, "ShotRecord.batch_role")
        if isinstance(self.seed, bool) or not isinstance(self.seed, int):
            raise SchemaError(f"ShotRecord.seed: expected an integer, got {self.seed!r}")
        _sha256_hex(self.protocol_sha256, "ShotRecord.protocol_sha256")
        _sha256_hex(self.implementation_sha256, "ShotRecord.implementation_sha256")
        _enum(self.outcome, SHOT_OUTCOMES, "ShotRecord.outcome")
        _identifier(self.trial_dir, "ShotRecord.trial_dir")
        if self.outcome in ("other", "unabsorbed", "gate_rejected", "failed", "crashed"):
            if not self.outcome_reason:
                raise SchemaError(
                    f"ShotRecord: outcome {self.outcome!r} without a reason is a dead "
                    "end; outcome_reason is required"
                )
        if (
            self.outcome in ("A", "B")
            and self.entry_side is None
            and (self.landed_microstate is None)
        ):
            raise SchemaError(
                "ShotRecord: an A/B verdict needs its evidence -- the first-entry "
                "side, or (fixed-rule fallback) the landed microstate"
            )
        if self.entry_side is not None:
            _enum(self.entry_side, ("A", "B"), "ShotRecord.entry_side")
            if self.entry_step is None:
                raise SchemaError(
                    "ShotRecord: an entry side without entry_step cannot be checked "
                    "against the protocol's k"
                )
            _int_at_least(self.entry_step, "ShotRecord.entry_step", 1)
        for name in ("entry_time_fs", "max_nve_drift_eV_atom", "cost_wall_s"):
            value = getattr(self, name)
            if value is not None:
                object.__setattr__(self, name, _finite(value, f"ShotRecord.{name}"))
        if self.cost_force_evaluations is not None:
            _int_at_least(self.cost_force_evaluations, "ShotRecord.cost_force_evaluations", 0)
        if self.model_checkpoint_sha256 is not None:
            _sha256_hex(self.model_checkpoint_sha256, "ShotRecord.model_checkpoint_sha256")
        if self.endpoint_sha256 is not None:
            _sha256_hex(self.endpoint_sha256, "ShotRecord.endpoint_sha256")

    def to_dict(self) -> dict:
        return {
            "proposal_id": self.proposal_id,
            "shot_id": self.shot_id,
            "batch_role": self.batch_role,
            "seed": self.seed,
            "protocol_sha256": self.protocol_sha256,
            "implementation_sha256": self.implementation_sha256,
            "outcome": self.outcome,
            "outcome_reason": self.outcome_reason,
            "trial_dir": self.trial_dir,
            "entry_side": self.entry_side,
            "entry_step": self.entry_step,
            "entry_time_fs": self.entry_time_fs,
            "landed_microstate": self.landed_microstate,
            "entry_quench_consistent": self.entry_quench_consistent,
            "max_nve_drift_eV_atom": self.max_nve_drift_eV_atom,
            "quench_converged": self.quench_converged,
            "cost_force_evaluations": self.cost_force_evaluations,
            "cost_wall_s": self.cost_wall_s,
            "endpoint_sha256": self.endpoint_sha256,
            "model_checkpoint_sha256": self.model_checkpoint_sha256,
        }

    @classmethod
    def from_dict(cls, data) -> "ShotRecord":
        _keys(
            data,
            [
                "proposal_id",
                "shot_id",
                "batch_role",
                "seed",
                "protocol_sha256",
                "implementation_sha256",
                "outcome",
                "trial_dir",
            ],
            [
                "outcome_reason",
                "entry_side",
                "entry_step",
                "entry_time_fs",
                "landed_microstate",
                "entry_quench_consistent",
                "max_nve_drift_eV_atom",
                "quench_converged",
                "cost_force_evaluations",
                "cost_wall_s",
                "endpoint_sha256",
                "model_checkpoint_sha256",
            ],
            "ShotRecord",
        )
        return cls(
            proposal_id=data["proposal_id"],
            shot_id=data["shot_id"],
            batch_role=data["batch_role"],
            seed=data["seed"],
            protocol_sha256=data["protocol_sha256"],
            implementation_sha256=data["implementation_sha256"],
            outcome=data["outcome"],
            outcome_reason=data.get("outcome_reason"),
            trial_dir=data["trial_dir"],
            entry_side=data.get("entry_side"),
            entry_step=data.get("entry_step"),
            entry_time_fs=data.get("entry_time_fs"),
            landed_microstate=data.get("landed_microstate"),
            entry_quench_consistent=data.get("entry_quench_consistent"),
            max_nve_drift_eV_atom=data.get("max_nve_drift_eV_atom"),
            quench_converged=data.get("quench_converged"),
            cost_force_evaluations=data.get("cost_force_evaluations"),
            cost_wall_s=data.get("cost_wall_s"),
            endpoint_sha256=data.get("endpoint_sha256"),
            model_checkpoint_sha256=data.get("model_checkpoint_sha256"),
        )


@dataclass(frozen=True)
class ShootingEstimate:
    """The conditional statistics of one proposal's shots (plan M3 统计与分类).

    ``p_B_given_AB`` is None while A + B == 0 -- it is conditional on reaching
    a side, and without conditions met it has no value (not 0, not 0.5). The
    band classification never upgrades an interval to a certainty: overlap
    with the band is "compatible", entirely below/above is a side, and
    anything inconclusive stays ``unresolved``.
    """

    proposal_id: str
    shots: int
    counts: dict  # outcome -> count, including every excluded label
    other_fraction: float | None
    conditional_denominator: int
    p_B_given_AB: float | None
    wilson_95: list | None
    entry_quench_mismatches: int
    band: tuple
    meaning: str

    @classmethod
    def from_records(cls, records, band=(0.4, 0.6)) -> "ShootingEstimate":
        if not records:
            raise SchemaError("ShootingEstimate.from_records: no shots to summarise")
        ids = {r.proposal_id for r in records}
        if len(ids) != 1:
            raise SchemaError(
                f"ShootingEstimate.from_records: records span {len(ids)} proposals; "
                "one estimate summarises one proposal"
            )
        counts = {outcome: 0 for outcome in SHOT_OUTCOMES}
        for record in records:
            counts[record.outcome] += 1
        denominator = counts["A"] + counts["B"]
        p_b = counts["B"] / denominator if denominator else None
        interval = wilson_interval(counts["B"], denominator) or None
        return cls(
            proposal_id=records[0].proposal_id,
            shots=len(records),
            counts=counts,
            other_fraction=(counts["other"] / len(records)) if records else None,
            conditional_denominator=denominator,
            p_B_given_AB=p_b,
            wilson_95=interval,
            entry_quench_mismatches=sum(
                1 for r in records if r.entry_quench_consistent is False
            ),
            band=tuple(band),
            meaning=(
                "p_B_given_AB is conditional on reaching A or B; unabsorbed, third-basin "
                "and failed shots keep their own labels and stay out of the denominator"
            ),
        )

    @property
    def band_class(self) -> str:
        """Plan M3's vocabulary; never a saddle claim."""
        low, high = self.band
        if self.p_B_given_AB is None or self.wilson_95 is None:
            return "unresolved"
        w_low, w_high = self.wilson_95
        if w_high < low:
            return "A_side"
        if w_high < high and w_low <= high:
            return "transition_region_compatible"
        if w_low > high:
            return "B_side"
        if w_low <= low and w_high >= high:
            return "transition_region_compatible"  # interval spans the whole band
        return "unresolved"

    def to_dict(self) -> dict:
        return {
            "proposal_id": self.proposal_id,
            "shots": self.shots,
            "counts": dict(self.counts),
            "other_fraction": self.other_fraction,
            "conditional_denominator": self.conditional_denominator,
            "p_B_given_AB": self.p_B_given_AB,
            "wilson_95": self.wilson_95,
            "entry_quench_mismatches": self.entry_quench_mismatches,
            "band": list(self.band),
            "band_class": self.band_class,
            "meaning": self.meaning,
        }


@dataclass(frozen=True)
class TransitionDiagnosis:
    """What the shooting says about one candidate -- and what it does not (E6).

    The label vocabulary contains no "verified" and no "TS": q near 0.5 is
    compatible with a transition region, not proof of a saddle. Saddle
    evidence stays where it always was -- refinement, curvature, two-sided
    connectivity, QM (plan section 7, gate E6).
    """

    proposal_id: str
    label: str
    basis: str
    checks: dict  # {"entry_quench_consistent": ..., "drift_within_tolerance": ..., ...}
    recommended_use: str = "no_action"
    parent_path_sha256: str | None = None

    def __post_init__(self):
        _identifier(self.proposal_id, "TransitionDiagnosis.proposal_id")
        _enum(self.label, DIAGNOSIS_LABELS, "TransitionDiagnosis.label")
        _identifier(self.basis, "TransitionDiagnosis.basis")
        if not isinstance(self.checks, dict):
            raise SchemaError("TransitionDiagnosis.checks: expected an object")
        _enum(self.recommended_use, RECOMMENDED_USES, "TransitionDiagnosis.recommended_use")
        if self.parent_path_sha256 is not None:
            _sha256_hex(self.parent_path_sha256, "TransitionDiagnosis.parent_path_sha256")

    def to_dict(self) -> dict:
        return {
            "proposal_id": self.proposal_id,
            "label": self.label,
            "basis": self.basis,
            "checks": dict(self.checks),
            "recommended_use": self.recommended_use,
            "parent_path_sha256": self.parent_path_sha256,
        }

    @classmethod
    def from_dict(cls, data) -> "TransitionDiagnosis":
        _keys(
            data,
            ["proposal_id", "label", "basis", "checks"],
            ["recommended_use", "parent_path_sha256"],
            "TransitionDiagnosis",
        )
        return cls(
            proposal_id=data["proposal_id"],
            label=data["label"],
            basis=data["basis"],
            checks=dict(data["checks"]),
            recommended_use=data.get("recommended_use", "no_action"),
            parent_path_sha256=data.get("parent_path_sha256"),
        )
