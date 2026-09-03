"""Engine-neutral execution contract shared by the ASE and OpenMM backends.

The search core speaks only this module: eV, Angstrom, eV/Angstrom, fs, amu, K.
Any other unit system (OpenMM's nm/kJ/mol/ps) must be converted inside a backend
adapter and never leak across this boundary.

The contract exists because the two backends differ in ways that cannot be hidden:
per-step Python gating is free under ASE and expensive under a native OpenMM
Context, committee force disagreement is reachable through an ASE calculator but
not through a single TorchForce, and a bias term can be wrapped around an ASE
calculator at will but must be present in an OpenMM System before its Context is
created. Each backend therefore declares what it can actually deliver, and the
declaration is persisted in the run manifest instead of being assumed.
"""
from dataclasses import dataclass, field
from typing import Protocol, runtime_checkable
import math
import numpy as np

UNITS = {"energy": "eV", "length": "Angstrom", "force": "eV/Angstrom",
         "time": "fs", "mass": "amu", "temperature": "K"}

BIAS_KINDS = ("pair_linear",)
UNCERTAINTY_MODES = ("available", "sampled", "unavailable")


class GateRejected(RuntimeError):
    """A reliability gate refused a configuration; the trial must not yield a basin."""

    def __init__(self, code, message, positions=None):
        super().__init__(message)
        self.code = code
        self.positions = None if positions is None else np.array(positions, dtype=float)


class EngineError(RuntimeError):
    """The backend itself failed (model, platform, memory); distinct from a gate refusal."""


@dataclass(frozen=True)
class EngineCapabilities:
    """What a backend promises. Recorded in the manifest; never silently widened."""

    name: str
    gate_delay_steps: int
    quench_gate_delay_steps: int
    uncertainty: str
    bias_kinds: tuple[str, ...] = ()
    bitwise_replay: bool = False
    platform: str | None = None
    precision: str | None = None
    notes: dict = field(default_factory=dict)

    def __post_init__(self):
        for name in ("gate_delay_steps", "quench_gate_delay_steps"):
            value = getattr(self, name)
            if type(value) is not int or value < 1:
                raise ValueError(f"{name} must be an integer >= 1")
        if self.uncertainty not in UNCERTAINTY_MODES:
            raise ValueError(f"uncertainty must be one of {UNCERTAINTY_MODES}")
        unknown = set(self.bias_kinds) - set(BIAS_KINDS)
        if unknown:
            raise ValueError(f"unknown bias kinds: {sorted(unknown)}")
        object.__setattr__(self, "bias_kinds", tuple(self.bias_kinds))

    def to_dict(self):
        return {"name": self.name, "gate_delay_steps": self.gate_delay_steps,
                "quench_gate_delay_steps": self.quench_gate_delay_steps,
                "uncertainty": self.uncertainty, "bias_kinds": list(self.bias_kinds),
                "bitwise_replay": self.bitwise_replay, "platform": self.platform,
                "precision": self.precision, "notes": dict(self.notes),
                "units": dict(UNITS)}


@dataclass(frozen=True)
class Eval:
    """One physical evaluation. uncertainty_eV_A is None when unavailable, never 0.0."""

    energy_eV: float
    forces_eV_A: np.ndarray
    uncertainty_eV_A: float | None = None

    @property
    def fmax_eV_A(self):
        return float(np.linalg.norm(self.forces_eV_A, axis=1).max())


@dataclass(frozen=True)
class BiasSpec:
    """Conservative pair bias U = -sign * amplitude * (r - r0), constant during a pulse.

    A backend that reports pair_linear must expose bias_energy_eV in every StateView
    and accumulate external work; a backend that cannot must omit the kind entirely
    rather than approximating it.
    """

    kind: str
    pair: tuple[int, int]
    sign: int
    amplitude_eV_A: float
    reference_distance_A: float

    def __post_init__(self):
        if self.kind not in BIAS_KINDS:
            raise ValueError(f"bias kind must be one of {BIAS_KINDS}")
        i, j = self.pair
        if type(i) is not int or type(j) is not int or i == j or min(i, j) < 0:
            raise ValueError("bias pair must be two distinct nonnegative atom indices")
        if self.sign not in (-1, 1):
            raise ValueError("bias sign must be +1 or -1")
        for name in ("amplitude_eV_A", "reference_distance_A"):
            value = getattr(self, name)
            if isinstance(value, bool) or not math.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be finite and positive")
        object.__setattr__(self, "pair", (i, j))

    def to_dict(self):
        return {"kind": self.kind, "pair": list(self.pair), "sign": self.sign,
                "amplitude_eV_A": self.amplitude_eV_A,
                "reference_distance_A": self.reference_distance_A,
                "potential": "U = -sign * amplitude * (r - r0)"}


@dataclass(frozen=True)
class StateView:
    """Everything an observer may inspect, so observers stay backend-agnostic."""

    phase: str
    step: int
    time_fs: float
    numbers: np.ndarray
    masses_amu: np.ndarray
    positions_A: np.ndarray
    momenta: np.ndarray
    evaluation: Eval
    bias_energy_eV: float = 0.0

    @property
    def kinetic_eV(self):
        return float(0.5 * np.sum(self.momenta ** 2 / self.masses_amu[:, None]))

    @property
    def physical_total_eV(self):
        return self.evaluation.energy_eV + self.kinetic_eV


@dataclass(frozen=True)
class SegmentResult:
    steps_run: int
    observed_steps: int
    stop_reason: str


@dataclass(frozen=True)
class MinimizeResult:
    """Convergence is always judged on fmax in eV/Angstrom by the controller's criterion.

    A backend whose native minimizer uses a different criterion must run in chunks
    and report the externally measured fmax here.
    """

    converged: bool
    steps: int
    fmax_eV_A: float


@runtime_checkable
class Engine(Protocol):
    """Owns positions, momenta and the potential for one trial; hands back plain arrays.

    Observers are called at least every ``gate_delay_steps`` steps and may raise
    GateRejected to abort a segment. A backend must abort promptly and leave its
    state inspectable; it must not swallow the exception.
    """

    def capabilities(self) -> EngineCapabilities: ...

    def model_identity(self) -> dict: ...

    def get_positions(self) -> np.ndarray: ...

    def set_positions(self, positions_A) -> None: ...

    def get_momenta(self) -> np.ndarray: ...

    def set_momenta(self, momenta) -> None: ...

    def evaluate(self) -> Eval: ...

    def view(self, phase: str, step: int, time_fs: float) -> StateView: ...

    def run_nve(self, steps: int, timestep_fs: float, observer, observe_every: int = 1,
                bias: BiasSpec | None = None) -> SegmentResult: ...

    def minimize(self, fmax_eV_A: float, max_steps: int, observer=None) -> MinimizeResult: ...

    def close(self) -> None: ...


def check_capabilities(capabilities, config, required_bias_kinds=()):
    """Fail closed before a run starts rather than mid-trial.

    Uncertainty demanded by the config but absent from the backend is a startup
    failure, per the unavailable-is-not-measured-zero rule.
    """
    if getattr(config, "max_uncertainty_eV_A", None) is not None:
        if capabilities.uncertainty == "unavailable":
            raise EngineError(
                f"config requires an uncertainty gate but backend '{capabilities.name}' "
                "supplies none; uncertainty is unavailable, not zero")
    missing = set(required_bias_kinds) - set(capabilities.bias_kinds)
    if missing:
        raise EngineError(
            f"backend '{capabilities.name}' cannot realize bias kinds {sorted(missing)}")
    return dict(capabilities.to_dict())
