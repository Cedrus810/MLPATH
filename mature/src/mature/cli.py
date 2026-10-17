"""Command line for the maturation layer: ``mature prepare|shoot|report``.

The three verbs are fixed by the construction plan
(``docs/PLAN_GEN_COMPAS_MATURATION_2026-10-03.md`` sections 1, 3.1 and 8):

    prepare  --protocol P.json --out DIR   candidates from PARS-selected paths (M2)
    shoot    --out DIR                     diagnostic shooting, resumable (M3/M4)
    report   --out DIR                     diagnosis, budgets, refinement handoff (M4)

``--dry-run`` on prepare and shoot says what would happen and writes nothing.
A batch directory is owned by its protocol: prepare refuses a non-empty
directory, shoot and report refuse a changed protocol, source, geometry or
prrs implementation (M1).
"""

from __future__ import annotations

import argparse
import json
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
    prepare = sub.add_parser("prepare", help=_STAGES["prepare"])
    prepare.add_argument("--protocol", required=True, help="frozen MaturationProtocol JSON")
    prepare.add_argument("--out", required=True, help="new batch directory")
    prepare.add_argument("--dry-run", action="store_true")
    shoot = sub.add_parser("shoot", help=_STAGES["shoot"])
    shoot.add_argument("--out", required=True, help="batch directory from prepare")
    shoot.add_argument("--max-new-shots", type=int, default=None)
    shoot.add_argument("--dry-run", action="store_true")
    report = sub.add_parser("report", help=_STAGES["report"])
    report.add_argument("--out", required=True, help="batch directory")

    args = parser.parse_args(argv)
    from . import maturation
    from .io import read_json_strict
    from .schema import MaturationProtocol

    if args.command == "prepare":
        protocol = MaturationProtocol.from_dict(read_json_strict(args.protocol))
        result = maturation.prepare(protocol, args.out, dry_run=args.dry_run)
        print(f"mature prepare: {len(result['proposals'])} candidates, {result['tasks']} tasks")
    elif args.command == "shoot":
        result = maturation.shoot_batch(
            args.out, max_new_shots=args.max_new_shots, dry_run=args.dry_run
        )
        print(f"mature shoot: {json.dumps(result, sort_keys=True)}")
    else:
        result = maturation.report(args.out)
        labels = {}
        for row in result["candidates"]:
            label = row.get("diagnosis", {}).get("label", row.get("status"))
            labels[label] = labels.get(label, 0) + 1
        print(f"mature report: {json.dumps(labels, sort_keys=True)}", file=sys.stdout)
    return 0
