"""Print the implementation_sha256 the way a run manifest computes it.

The hash in docs/HANDOFF.md used to be copied by hand and went stale: it read
9c13fdb5cec8e735 while the tree hashed to dba26be31fd909f5, two rounds of fixes later.
A wrong hash in the handoff is worse than none -- it names a code state that produced
none of the runs on disk.

The algorithm is duplicated from prrs/search.py rather than imported, and that is the
point: importing it would make this script agree with a broken implementation. If the
two ever disagree, the run manifests are the truth and this script is the bug.

    python scripts/implementation_hash.py            # the current tree
    python scripts/implementation_hash.py --check    # compare against every run manifest
"""
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def implementation_sha256():
    digest = hashlib.sha256()
    for path in sorted((ROOT / "src" / "prrs").glob("*.py")):
        digest.update(path.name.encode())
        digest.update(path.read_bytes())
    return digest.hexdigest()


def main(argv):
    current = implementation_sha256()
    print(f"implementation_sha256 = {current}")
    print(f"            short     = {current[:16]}")
    if "--check" not in argv:
        return 0
    seen = {}
    for manifest in sorted((ROOT / "runs").glob("*/config.json")):
        try:
            recorded = json.loads(manifest.read_text())["implementation_sha256"]
        except (ValueError, KeyError):
            continue
        seen.setdefault(recorded[:16], []).append(manifest.parent.name)
    print(f"\n{len(seen)} implementation(s) recorded under runs/:")
    for short, runs in sorted(seen.items(), key=lambda kv: -len(kv[1])):
        mark = "  <- current tree" if short == current[:16] else ""
        print(f"  {short}  {len(runs):3d} run(s){mark}")
        print(f"      {', '.join(sorted(runs)[:6])}{' ...' if len(runs) > 6 else ''}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
