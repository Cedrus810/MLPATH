"""Do all the seeds of a scan share one implementation?

A scan whose seeds ran under different code is a mixed-implementation dataset and its
ratios are not one measurement. The first T=300 batch had 13 seeds on one implementation
hash and 3 on another -- the three that had crashed on a non-finite diagnostic value and
were re-run after the fix. The fix looked confined to the record-writing layer, but this
directory has no git history, so there was no way to prove that independently. "It looks
confined" is a judgement, not a checkable fact, so the batch had to be called mixed.
"""
import json
import sys
from collections import Counter
from pathlib import Path

prefix = sys.argv[1]
seeds = sys.argv[2:]
seen = {}
for seed in seeds:
    manifest = Path(f"{prefix}{seed}/config.json")
    if not manifest.exists():
        seen[seed] = "missing"
        continue
    seen[seed] = json.loads(manifest.read_text())["implementation_sha256"][:8]

counts = Counter(seen.values())
for value, number in counts.most_common():
    print(f"  {value}  {number} seed(s): {[s for s in seen if seen[s] == value]}")
real = [value for value in counts if value != "missing"]
if len(real) > 1:
    print("\nMIXED IMPLEMENTATION: these ratios are not one measurement. Re-run the whole")
    print("batch under a single implementation before quoting them.")
    sys.exit(1)
if not real:
    print("\nno runs found")
    sys.exit(2)
print(f"\nsingle implementation across {len(seeds)} seeds: {real[0]}")
