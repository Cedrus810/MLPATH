"""Maturation layer for PRRS reaction channels.

An independent module next to PRRS and PARS. PRRS discovers reaction
channels; PARS selects and refines representative paths; the maturation
layer adds the stage between selection and refinement: candidate
geometries taken from already-discovered paths, diagnostic shooting with
independent thermal momenta, first-entry verdicts, transition-region
diagnosis, and a traceable refinement handoff (construction plan
``docs/PLAN_GEN_COMPAS_MATURATION_2026-10-03.md``, sections 1 and 3).

Module contract: ``mature`` imports ``prrs`` one-way -- dynamics and
verdicts are reused (``runner.shoot``, read-only ``Registry.identify``),
never reimplemented. Nothing imports ``mature`` back: neither ``src/prrs``
nor ``pars/src/pars`` may depend on this layer, enforced statically by
``tests/test_boundary.py``. PARS's own file-only boundary is unaffected.
"""

__all__ = ["__version__"]

__version__ = "0.1.0"
