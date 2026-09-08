"""Read a potential's declared domain out of its own checkpoint.

Two different questions get confused because both produce a list of elements:

  1. "Which elements can this potential evaluate at all?"  -- a DOMAIN question. The
     checkpoint answers it: `model.atomic_numbers` is what the model was built with.
     Nothing in prrs currently checks a structure against it, so a run on an element
     the model never saw produces numbers with no warning.

  2. "What is this element's neutral closed-shell valence?" -- a LEWIS STRUCTURE
     question, answered by prrs.chemistry.VALENCE. The model has NOTHING to say about
     it. Palladium has no single-valued valence; that is chemistry, not metadata.
     VALENCE is deliberately small and reports anything outside it as unresolved
     rather than guessing, and reading it off a checkpoint would be exactly the guess
     it refuses to make.

This script answers (1) only. Its output is what a benchmark's model.json should record
and what a preflight should gate on. It never touches (2).

    python scripts/model_registry.py <checkpoint> [<checkpoint> ...]
    python scripts/model_registry.py --json <checkpoint>        # machine-readable
"""

import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def describe(path):
    import torch
    from ase.data import chemical_symbols

    from prrs import torch_guard

    torch_guard.install()
    path = Path(path)
    model = torch.load(str(path), map_location="cpu", weights_only=False)
    numbers = [int(z) for z in getattr(model, "atomic_numbers", [])]
    r_max = getattr(model, "r_max", None)
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return {
        "path": str(path),
        "sha256": digest.hexdigest(),
        "checkpoint_class": type(model).__name__,
        "heads": list(getattr(model, "heads", []) or []),
        "r_max_A": None if r_max is None else float(r_max),
        "element_count": len(numbers),
        "atomic_numbers": numbers,
        "symbols": [chemical_symbols[z] for z in numbers],
    }


def coverage(spec):
    """What prrs can and cannot say about this model's elements."""
    from prrs.chemistry import VALENCE

    inside = sorted(z for z in spec["atomic_numbers"] if z in VALENCE)
    outside = sorted(z for z in spec["atomic_numbers"] if z not in VALENCE)
    return {
        "valence_table_covers": len(inside),
        "valence_table_missing": len(outside),
        "meaning": (
            "Elements the model can evaluate but prrs.chemistry.VALENCE has no valence "
            "for are NOT refused: chemical_key records stereo_unresolved = "
            "'unsupported_elements' and falls back to graph plus fragments. So they are "
            "searchable without stereochemistry. Only a graph that HAS valences and "
            "cannot satisfy them is refused as out_of_domain."
        ),
    }


def main(argv):
    as_json = "--json" in argv
    paths = [a for a in argv if not a.startswith("--")]
    if not paths:
        print(__doc__.strip())
        return 2
    out = []
    for path in paths:
        spec = describe(path)
        spec["prrs_identity_coverage"] = coverage(spec)
        out.append(spec)
        if as_json:
            continue
        from ase.data import chemical_symbols

        from prrs.chemistry import VALENCE

        print(f"{Path(path).name}")
        print(f"  class        {spec['checkpoint_class']}   heads={spec['heads']}")
        print(f"  sha256       {spec['sha256']}")
        print(f"  r_max        {spec['r_max_A']} A")
        print(f"  elements     {spec['element_count']}")
        print(f"               {' '.join(spec['symbols'][:20])}"
              f"{' ...' if spec['element_count'] > 20 else ''}")
        cov = spec["prrs_identity_coverage"]
        print(f"  prrs VALENCE covers {cov['valence_table_covers']} of them; "
              f"{cov['valence_table_missing']} fall back to graph+fragments only")
        missing_here = sorted(set(VALENCE) - set(spec["atomic_numbers"]))
        if missing_here:
            print(f"  NOTE: prrs has valences for {[chemical_symbols[z] for z in missing_here]} "
                  f"which this model cannot evaluate")
        print()
    if as_json:
        print(json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
