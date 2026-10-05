"""Backend-independent refinement contract.

A RefinementBackend turns a PathCandidate into a RefinedPath. It receives the
calculator from the caller, so PARS never imports a physics engine: the CLI
resolves a user-named factory (``module:attribute``) at runtime and a test can
pass an analytic potential directly.

Status semantics (design section 18): a refinement that converges but leaves
the transition region more uncertain than the configured refusal threshold is
reported UNRESOLVED_MLP, never silently as success. Failures keep the chain
data they produced.
"""

from __future__ import annotations

from typing import Protocol

from ..config import ParsConfig
from ..schema import PathCandidate, RefinedPath, RefinementStatus


class RefinementBackend(Protocol):
    def refine(self, candidate: PathCandidate, config: ParsConfig) -> RefinedPath:
        """Refine one candidate into a RefinedPath."""
        ...


__all__ = ["RefinementBackend", "RefinedPath", "RefinementStatus"]
