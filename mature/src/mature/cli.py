"""Command line shape for the maturation layer: ``mature prepare|shoot|report``.

The three verbs and the order they land in are fixed by the construction
plan (``docs/PLAN_GEN_COMPAS_MATURATION_2026-10-03.md`` sections 1, 3.1
and 8): ``prepare`` extracts candidate geometries from PARS-selected
paths (M2), ``shoot`` runs diagnostic shooting with independent thermal
momenta (M3), ``report`` writes per-channel diagnosis, budgets and the
refinement handoff (M4).

Until then every subcommand refuses with exit code 2 instead of
pretending to run: a stage that does not exist must never look like a
stage that succeeded (fail closed, HANDOFF section 2).
"""

from __future__ import annotations

import argparse
import sys

_STAGES = {
    "prepare": "extract candidate geometries from PARS-selected paths (M2)",
    "shoot": "diagnostic shooting with independent thermal momenta (M3)",
    "report": "per-channel diagnosis, budgets and refinement handoff (M4)",
}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="mature",
        description=(
            "Path maturation: shooting-based transition-region diagnosis and "
            "refinement handoff for PRRS-discovered reaction channels."
        ),
    )
    sub = parser.add_subparsers(dest="command", required=True)
    for name, text in _STAGES.items():
        sub.add_parser(name, help=text)

    args = parser.parse_args(argv)
    print(
        f"mature {args.command}: not implemented yet "
        f"(lands with M1-M4 of docs/PLAN_GEN_COMPAS_MATURATION_2026-10-03.md)",
        file=sys.stderr,
    )
    return 2
