"""Command line entry points; use python -m prrs or the prrs executable."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
from ase.io import read, write
from .calculators import demo_atoms, double_well_factory, load_factory
from .config import SearchConfig
from .search import ReactionSearch


def load_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def factory_metadata(spec, kwargs):
    metadata = {"factory": spec, "kwargs": kwargs, "local_file_sha256": {}}
    for value in kwargs.values():
        values = value if isinstance(value, list) else [value]
        for item in values:
            if isinstance(item, str) and Path(item).is_file():
                with Path(item).open("rb") as handle:
                    digest = hashlib.file_digest(handle, "sha256").hexdigest() if hasattr(hashlib, "file_digest") else hashlib.sha256(handle.read()).hexdigest()
                metadata["local_file_sha256"][str(Path(item).resolve())] = digest
    return metadata


def export_path(run, attempt, output, images):
    if images < 3:
        raise ValueError("At least three images required")
    destination = Path(output)
    if destination.exists():
        raise FileExistsError(f"Will not overwrite {destination}")
    run = Path(run)
    network = load_json(run / "network.json")
    if not any(attempt in edge["attempts"] for edge in network["reactions"]):
        raise ValueError("Path export requires an accepted channel attempt")
    record = load_json(run / "trials" / attempt / "result.json")
    frames = read(str(run / record["trajectory"]), index=":")
    # Quench iterates are not MD; use only source + physical trajectory + endpoint.
    frames = [f for f in frames if f.info.get("phase") not in ("quench", "quenched")]
    endpoint = read(str(run / record["endpoint"]))
    frames.append(endpoint)
    coordinates = np.asarray([f.positions for f in frames])
    arc = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(coordinates, axis=0).reshape(len(frames)-1, -1), axis=1))])
    keep = np.concatenate([[True], np.diff(arc) > 1e-12])
    arc, coordinates = arc[keep], coordinates[keep]
    if len(arc) < 2:
        raise ValueError("Trajectory has zero path length")
    targets = np.linspace(0, arc[-1], images)
    flat = coordinates.reshape(len(arc), -1)
    interpolated = np.column_stack([np.interp(targets, arc, flat[:, i]) for i in range(flat.shape[1])])
    result = []
    for position in interpolated.reshape(images, len(endpoint), 3):
        frame = endpoint.copy()
        frame.calc = None
        frame.set_positions(position)
        frame.set_momenta(np.zeros((len(frame), 3)))
        frame.info = {"purpose": "unrefined_path_seed", "parent_attempt": attempt}
        result.append(frame)
    destination.parent.mkdir(parents=True, exist_ok=True)
    write(str(destination), result)
    return {"output": str(destination), "images": images,
            "warning": "Unrefined interpolated path seed; not a TS, MEP or barrier"}


def main(argv=None):
    parser = argparse.ArgumentParser(description="PRRS: perturb, respond, quench, discover")
    sub = parser.add_subparsers(dest="command", required=True)
    demo = sub.add_parser("demo", help="Offline analytic double-well benchmark (not chemical validation)")
    demo.add_argument("--output", required=True)
    demo.add_argument("--seed", type=int, default=42)
    search = sub.add_parser("search", help="Search with an explicitly trusted ASE calculator factory")
    search.add_argument("input")
    search.add_argument("--output", required=True)
    search.add_argument("--config")
    search.add_argument("--calculator", required=True, help="module:function returning an ASE Calculator")
    search.add_argument("--calculator-kwargs", help="JSON file with factory keyword arguments")
    path = sub.add_parser("export-path", help="Export an observed channel as a local refinement seed")
    path.add_argument("--run", required=True)
    path.add_argument("--attempt", required=True)
    path.add_argument("--output", required=True)
    path.add_argument("--images", type=int, default=9)
    args = parser.parse_args(argv)
    try:
        if args.command == "demo":
            config = SearchConfig(seed=args.seed, max_trials=48, max_depth=2,
                                  geometry_amplitudes_A=(0.1, 0.3, 0.7, 1.0),
                                  kick_energies_eV=(0.1, 0.5, 1.0), timestep_fs=0.1,
                                  # The demo's potential is a synthetic radial double well,
                                  # explicitly not carbon chemistry, so a valence-based
                                  # domain gate has nothing to say about it.
                                  closed_shell_only=False,
                                  response_steps=150, quench_fmax_eV_A=0.001)
            result = ReactionSearch(double_well_factory, config,
                                    {"factory": "prrs.calculators:double_well_factory",
                                     "purpose": "synthetic software benchmark"}).run(demo_atoms(), args.output)
        elif args.command == "search":
            config = SearchConfig.from_dict(load_json(args.config)) if args.config else SearchConfig()
            kwargs = load_json(args.calculator_kwargs) if args.calculator_kwargs else {}
            factory = load_factory(args.calculator, kwargs)
            result = ReactionSearch(factory, config, factory_metadata(args.calculator, kwargs)).run(
                read(args.input), args.output)
        else:
            print(json.dumps(export_path(args.run, args.attempt, args.output, args.images), indent=2))
            return 0
        print(json.dumps({"output": str(Path(args.output).resolve()), "status": result["status"],
                          "chemical_nodes": len(result["chemical_nodes"]),
                          "microstates": result["budgets"]["microstates_total"],
                          "reactions": len(result["reactions"]),
                          "ts_candidates": len(result["ts_candidates"]),
                          "attempts": result["attempts"], "outcome_counts": result["outcome_counts"]}, indent=2))
        return 0
    except Exception as exc:
        print(f"prrs: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 2

