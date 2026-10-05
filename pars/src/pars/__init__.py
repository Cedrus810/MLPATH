"""PARS: Pareto-Adaptive Reaction Path Search.

An independent module next to PRRS. PRRS discovers reactions; PARS organizes
and selects reaction paths, then hands representative reaction tubes to NEB /
string refinement and QM validation (design document sections 1-2).

Independence contract: PARS consumes PRRS run directories as files and never
imports ``prrs``. The whole ``pars`` folder can be lifted out of the repository
intact; see README.md.
"""

__all__ = ["__version__"]

__version__ = "0.1.0"
