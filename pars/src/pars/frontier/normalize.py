"""Robust normalization (design section 10).

Raw objectives have incompatible units, so the frontier operates on
median/IQR-normalized values. The fitted parameters are part of the result
record: changing them silently would change the frontier. A zero IQR (all
candidates identical on an objective) carries no discrimination; the scale
falls back to a configured physical scale if one was given, else to 1.0.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..config import ParsConfig


@dataclass(frozen=True)
class Normalization:
    method: str
    center: dict[str, float]
    scale: dict[str, float]

    def normalize(self, raw: dict[str, float | None]) -> dict[str, float | None]:
        out: dict[str, float | None] = {}
        for name, value in raw.items():
            if value is None:
                out[name] = None
            else:
                out[name] = float((value - self.center[name]) / self.scale[name])
        return out

    def to_dict(self) -> dict:
        return {
            "method": self.method,
            "center": dict(self.center),
            "scale": dict(self.scale),
        }


def fit_normalization(
    raw_by_candidate: dict[str, dict[str, float | None]], config: ParsConfig
) -> Normalization:
    configured_scales = dict(config.objective_scales)
    names: list[str] = []
    for raw in raw_by_candidate.values():
        for name in raw:
            if name not in names:
                names.append(name)
    center: dict[str, float] = {}
    scale: dict[str, float] = {}
    for name in names:
        values = np.asarray(
            [raw[name] for raw in raw_by_candidate.values() if raw.get(name) is not None],
            dtype=float,
        )
        if values.size == 0:
            # The objective is unavailable for the whole channel (e.g. no
            # committee ran): it can never discriminate, so it is left out of
            # the normalized space entirely instead of being fabricated.
            continue
        center[name] = float(np.median(values))
        iqr = float(np.quantile(values, 0.75) - np.quantile(values, 0.25))
        if iqr < 1e-12:
            iqr = float(configured_scales.get(name, 1.0))
        scale[name] = iqr
    return Normalization(method=config.normalization, center=center, scale=scale)
