"""DiscoveryTrace -> trimmed PathCandidate (design sections 6.2 and 17.1).

The event window keeps the frames whose phase the run itself recorded as the
free response, drops the artificial forcing stage (its parameters stay in
provenance) and the quench/relaxation tail. The product endpoint snapshot is
NOT appended to the candidate: endpoints are re-minimized on the refinement
PES anyway (design 17.2), and the endpoint snapshot carries no bond signature,
which would corrupt path distances.

One candidate is emitted per (trial, channel) pair; a trial observed by
several channels is rejected upstream in the adapter.
"""

from __future__ import annotations

import numpy as np

from ..config import ParsConfig
from ..schema import DiscoveryTrace, Provenance, PathCandidate


def build_candidates(traces: list[DiscoveryTrace], config: ParsConfig) -> list[PathCandidate]:
    candidates = []
    for trace in traces:
        window = [i for i, phase in enumerate(trace.phases) if phase in config.include_phases]
        if trace.provenance.channel_id is None:
            raise ValueError(
                f"trace {trace.trace_id}: no reaction channel; PARS organizes "
                "paths per channel and cannot place this trace"
            )
        if len(window) < config.min_path_frames:
            raise ValueError(
                f"trace {trace.trace_id}: only {len(window)} frames match "
                f"phases {config.include_phases}; need >= {config.min_path_frames}"
            )
        first, last = window[0], window[-1] + 1
        kept = np.asarray(window, dtype=int)
        candidate_id = f"pc-{trace.provenance.channel_id}-{trace.trace_id}"
        provenance_dict = trace.provenance.to_dict()
        extra = dict(provenance_dict.get("extra") or {})
        # exact original frame indices of every path image: the window need not
        # be contiguous (the forcing stage sits between source and free_start)
        extra["kept_frames"] = [int(i) for i in window]
        # Basin witnesses come from the source run, not the trimmed free
        # segment. Preserve atom mapping without re-deriving chemical keys.
        if trace.endpoint_positions is not None:
            extra["endpoint_basin_references"] = {
                "reactant": trace.positions[0].tolist(),
                "product": trace.endpoint_positions.tolist(),
            }
        provenance_dict["extra"] = extra
        provenance_dict["frame_start"] = first
        provenance_dict["frame_stop"] = last
        candidates.append(
            PathCandidate(
                candidate_id=candidate_id,
                channel_id=trace.provenance.channel_id,
                trace_id=trace.trace_id,
                numbers=trace.numbers,
                positions=np.ascontiguousarray(trace.positions[kept]),
                energies=np.ascontiguousarray(trace.energies[kept]),
                uncertainty=(None if trace.uncertainty is None else trace.uncertainty[kept]),
                signatures=(
                    None
                    if trace.signatures is None
                    else tuple(trace.signatures[i] for i in window)
                ),
                provenance=Provenance(**provenance_dict),
            )
        )
    return candidates
