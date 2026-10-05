"""Reader for a PRRS run directory. This is the entire integration surface.

It parses files that PRRS wrote and nothing else. It must stay free of
``import prrs`` so the ``pars`` package can be lifted out of this repository
intact; a boundary test pins that rule.

Expected layout (the directory that contains ``network.json``, i.e. the
``<seed>/<arm>`` directory of a connect run or the run root of a discover run):

    network.json      reaction channels (undirected identities) and reactions
    config.json       run manifest: units, calculator, implementation_sha256
    attempts.jsonl    one record per trial, append-only chronological order
    trials/<id>/      response.traj, observations.jsonl, result.json, endpoint.extxyz

Field notes verified against real runs:

- ``observations.jsonl`` rows correspond 1:1 with the frames of
  ``response.traj`` (row ``frame`` == trajectory index). Each row carries
  ``phase`` (source/perturbed/free_start/free/quench/quenched), ``potential_eV``,
  ``uncertainty_eV_A`` (null when the run has no committee) and the per-frame
  bond list ``edges``.
- ``network.json`` ``reaction_channels`` entries carry ``id``, ``event_key`` and
  ``observed_by`` (attempt ids). Directed ``reactions`` entries carry
  ``attempts`` and the ``source``/``target`` chemical node ids.

Electronic state: ``config.json`` ``config.total_charge`` / ``multiplicity``
(None meaning PRRS's neutral-singlet default) is the declared state. PRRS stamps
it as ``info["charge"]`` / ``["spin"]`` on every response frame and endpoint;
each frame must agree with the declaration, and unstamped frames are accepted
only when the declared state is the default the calculator would assume anyway.

Fail-closed: old runs without ``reaction_channels`` or a run config, frame/row
counts or indices that disagree, phase/species-sequence and electronic-state
mismatches, and conflicting ledger/result classifications are rejected. A missing ``edges`` field disables signatures
rather than guessing bonds. Species-sequence checks cannot detect permutations
of identical atoms; the producer must preserve its mapped trajectory order.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from ase.io import read as ase_read

from .. import electronic
from ..iox import read_json, read_jsonl
from ..schema import DiscoveryTrace, Provenance

_RECORD_DROP_KEYS = frozenset({"probe", "trajectory", "observations", "endpoint", "attempt_id"})


class PrrsRun:
    """Read-only view of one PRRS run directory."""

    def __init__(self, run_dir):
        self.run_dir = Path(run_dir)
        if not self.run_dir.is_dir():
            raise FileNotFoundError(f"run directory does not exist: {self.run_dir}")
        self.network = read_json(self.run_dir / "network.json")
        if "reaction_channels" not in self.network:
            raise ValueError(
                f"{self.run_dir}: network.json has no reaction_channels; PARS "
                "requires a new-format run whose channel identities are already "
                "materialized. Re-run or upgrade the run rather than re-deriving "
                "symmetry logic here."
            )
        self.manifest = read_json(self.run_dir / "config.json")
        self.electronic_state = self._declared_electronic_state()
        self._channels = {c["id"]: c for c in self.network["reaction_channels"]}
        self._channels_by_attempt: dict[str, list[str]] = {}
        for channel_id, channel in self._channels.items():
            for attempt_id in channel.get("observed_by", ()):
                self._channels_by_attempt.setdefault(attempt_id, []).append(channel_id)
        self._reactions_by_attempt: dict[str, list[dict]] = {}
        for reaction in self.network.get("reactions", ()):
            for attempt_id in set(reaction.get("attempts", ())):
                self._reactions_by_attempt.setdefault(attempt_id, []).append(reaction)
        self._attempts = read_jsonl(self.run_dir / "attempts.jsonl")
        self._attempts_by_id = {r["attempt_id"]: r for r in self._attempts}

    # -- lookups ------------------------------------------------------------

    @property
    def units(self) -> dict:
        return dict(self.manifest.get("units", {}))

    @property
    def implementation_sha256(self) -> str | None:
        value = self.manifest.get("implementation_sha256")
        return value if isinstance(value, str) else None

    def attempts(self) -> list[dict]:
        return list(self._attempts)

    def checked_attempts(self) -> list[dict]:
        """Authoritative per-trial records in ledger order, each cross-checked
        against its ledger row. Reads result.json only, never trajectories, so
        filtering on it cannot skip a trial whose two records disagree."""
        return [self._trial_record(record["attempt_id"]) for record in self._attempts]

    def attempt(self, attempt_id: str) -> dict:
        try:
            return self._attempts_by_id[attempt_id]
        except KeyError:
            raise KeyError(f"unknown attempt id: {attempt_id}") from None

    def channel_ids_for(self, attempt_id: str) -> list[str]:
        return list(self._channels_by_attempt.get(attempt_id, ()))

    def reaction_for(self, attempt_id: str) -> dict | None:
        reactions = self._reactions_by_attempt.get(attempt_id, ())
        if len(reactions) > 1:
            raise ValueError(
                f"{attempt_id}: listed by multiple reactions "
                f"{[r.get('id') for r in reactions]}; refusing to pick one silently"
            )
        return reactions[0] if reactions else None

    # -- trace loading --------------------------------------------------------

    def load_trace(self, attempt_id: str) -> DiscoveryTrace:
        """Load the full raw trajectory of one attempt as a DiscoveryTrace.

        frame_start/frame_stop in the provenance cover the whole trajectory;
        trimming to the reaction segment happens in paths.candidate.
        """
        record = self._trial_record(attempt_id)
        trajectory_path = self.run_dir / record["trajectory"]
        observations_path = self.run_dir / record["observations"]
        frames = ase_read(str(trajectory_path), index=":")
        rows = read_jsonl(observations_path)
        if len(frames) != len(rows):
            raise ValueError(
                f"{attempt_id}: response trajectory has {len(frames)} frames but "
                f"{len(rows)} observation rows; refusing to guess the pairing"
            )
        if not frames:
            raise ValueError(f"{attempt_id}: empty response trajectory")
        numbers = np.asarray(frames[0].get_atomic_numbers(), dtype=int)
        for index, (frame, row) in enumerate(zip(frames, rows)):
            if type(row.get("frame")) is not int or row["frame"] != index:
                raise ValueError(
                    f"{attempt_id}: observation frame index mismatch at row {index}; refusing to guess the pairing"
                )
            if "phase" in frame.info and frame.info["phase"] != row.get("phase"):
                raise ValueError(f"{attempt_id}: observation phase mismatch at frame {index}")
            if not np.array_equal(frame.get_atomic_numbers(), numbers):
                raise ValueError(f"{attempt_id}: atomic mapping mismatch at frame {index}")
        positions = np.asarray([f.positions for f in frames], dtype=float)
        states = [
            electronic.info_state(f.info, f"{attempt_id} frame {i}")
            for i, f in enumerate(frames)
        ]
        energies = np.asarray([float(row["potential_eV"]) for row in rows], dtype=float)
        phases = tuple(str(row["phase"]) for row in rows)

        signatures = None
        if all(row.get("edges") is not None for row in rows):
            signatures = tuple(
                frozenset((int(i), int(j)) for i, j in row["edges"]) for row in rows
            )

        uncertainty_values = [row.get("uncertainty_eV_A") for row in rows]
        uncertainty = None
        if any(value is not None for value in uncertainty_values):
            uncertainty = np.asarray(
                [np.nan if v is None else float(v) for v in uncertainty_values],
                dtype=float,
            )

        endpoint_positions, endpoint_energy = None, None
        if record.get("endpoint"):
            endpoint = ase_read(str(self.run_dir / record["endpoint"]))
            if not np.array_equal(endpoint.get_atomic_numbers(), numbers):
                raise ValueError(f"{attempt_id}: endpoint atomic mapping mismatch")
            states.append(electronic.info_state(endpoint.info, f"{attempt_id} endpoint"))
            endpoint_positions = np.asarray(endpoint.positions, dtype=float)
            endpoint_energy = (
                float(endpoint.get_potential_energy()) if endpoint.calc is not None else None
            )

        self._check_recorded_states(attempt_id, states)
        provenance = self._provenance(attempt_id, record, 0, len(frames))
        return DiscoveryTrace(
            trace_id=attempt_id,
            numbers=numbers,
            positions=positions,
            energies=energies,
            phases=phases,
            signatures=signatures,
            uncertainty=uncertainty,
            endpoint_positions=endpoint_positions,
            endpoint_energy=endpoint_energy,
            provenance=provenance,
        )

    # -- internals ------------------------------------------------------------

    def _declared_electronic_state(self) -> dict:
        config = self.manifest.get("config")
        if not isinstance(config, dict):
            raise ValueError(
                f"{self.run_dir}: config.json has no run config, so the electronic "
                "state (total_charge/multiplicity) the potential saw is unknown"
            )
        declared = {
            "charge": config.get("total_charge"),
            "spin": config.get("multiplicity"),
        }
        return electronic.checked_state(
            {k: electronic.DEFAULTS[k] if v is None else v for k, v in declared.items()},
            f"{self.run_dir}/config.json",
        )

    def _check_recorded_states(self, attempt_id: str, states: list[dict]) -> None:
        """Every recorded structure carries the declared state, or (default state
        only) none of them carries the key at all."""
        for key in electronic.KEYS:
            declared = self.electronic_state[key]
            values = [state[key] for state in states]
            if all(value is None for value in values):
                if declared != electronic.DEFAULTS[key]:
                    raise ValueError(
                        f"{attempt_id}: frames record no electronic state {key} but the "
                        f"run declares {declared}; cannot certify which PES produced the path"
                    )
                continue
            if any(value != declared for value in values):
                raise ValueError(
                    f"{attempt_id}: recorded electronic state {key} "
                    f"{sorted(set(values), key=str)} disagrees with the declared {declared}"
                )

    def _trial_record(self, attempt_id: str) -> dict:
        """result.json is the authoritative per-trial record when present."""
        result_path = self.run_dir / "trials" / attempt_id / "result.json"
        if result_path.is_file():
            result = read_json(result_path)
            ledger = self.attempt(attempt_id)
            for key in ("attempt_id", "status", "classification"):
                if result.get(key) != ledger.get(key):
                    raise ValueError(f"{attempt_id}: result/attempt ledger mismatch in {key}")
            return result
        return self.attempt(attempt_id)

    def _provenance(
        self, attempt_id: str, record: dict, frame_start: int, frame_stop: int
    ) -> Provenance:
        probe = record.get("probe") or {}
        reaction = self.reaction_for(attempt_id)
        channel_ids = self.channel_ids_for(attempt_id)
        if len(channel_ids) > 1:
            raise ValueError(
                f"{attempt_id}: observed by multiple reaction channels "
                f"{channel_ids}; refusing to pick one silently"
            )
        return Provenance(
            run_dir=str(self.run_dir),
            attempt_id=attempt_id,
            channel_id=channel_ids[0] if channel_ids else None,
            event_key=(
                self._channels[channel_ids[0]].get("event_key") if channel_ids else None
            ),
            reactant_state_id=reaction.get("source") if reaction else None,
            product_state_id=reaction.get("target") if reaction else None,
            frame_start=frame_start,
            frame_stop=frame_stop,
            perturbation_family=probe.get("family"),
            perturbation_amplitude=probe.get("amplitude"),
            perturbation_unit=probe.get("unit"),
            seed=probe.get("seed"),
            implementation_sha256=self.implementation_sha256,
            calculator=dict(self.manifest.get("calculator", {})),
            electronic_state=dict(self.electronic_state),
            extra={k: v for k, v in record.items() if k not in _RECORD_DROP_KEYS},
        )
