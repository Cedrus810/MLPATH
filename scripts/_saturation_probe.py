"""How many concurrent processes does this GPU actually reward?

Aggregate throughput stops rising once the per-call time that cannot be overlapped is
filled. Measured on an RTX 2080 Ti the ceiling was 3.2x (reached by 7 workers, not improved
by 20); on an RTX 5080 it was 1.70x while per-call time grew 13.9x from 1 to 24 workers.
A faster card makes the GPU part shorter, so a LOWER ceiling on the faster card says the
contended resource is not GPU compute -- the prime suspect is CPU thread oversubscription,
since a 12-atom model gets nothing from intra-op parallelism yet torch opens a pool per
process anyway.

So by default each worker gets one thread and one dedicated core. Compare against the
library defaults with THREADS=0, and drop the affinity with PIN=0.
"""
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

kwargs_path = sys.argv[1]
calls = int(sys.argv[2]) if len(sys.argv) > 2 else 80
THREADS = os.environ.get("THREADS", "1")
PIN = os.environ.get("PIN", "1") != "0"
CORES = list(range(os.cpu_count() or 1))
HAVE_TASKSET = shutil.which("taskset") is not None
pin_env = {} if THREADS == "0" else {
    "OMP_NUM_THREADS": THREADS, "MKL_NUM_THREADS": THREADS,
    "OPENBLAS_NUM_THREADS": THREADS, "NUMEXPR_NUM_THREADS": THREADS}

worker = Path(os.environ.get("TMPDIR", "/tmp")) / f"p2_probe_{os.getpid()}.py"
lines = ["import json, os, sys, time"]
if THREADS != "0":
    lines.append(f"import torch; torch.set_num_threads({THREADS})")
lines += [
    "sys.path.insert(0, 'src')",
    "from ase.io import read",
    "from prrs.calculators import load_factory",
    "atoms = read('runs/p2_oxobutanal/p2_input.extxyz')",
    f"atoms.calc = load_factory('prrs.calculators:mace_factory', json.load(open({kwargs_path!r})))()",
    "atoms.get_potential_energy()",
    "n = int(sys.argv[1])",
    "t = time.perf_counter()",
    "for _ in range(n):",
    "    atoms.positions[0, 0] += 1e-6",
    "    atoms.get_potential_energy(); atoms.get_forces()",
    "print(json.dumps({'calls': n, 'seconds': time.perf_counter() - t}))",
]
worker.write_text("\n".join(lines) + "\n", encoding="utf-8")

print(f"cores {len(CORES)}   threads/process "
      + ("library default" if THREADS == "0" else THREADS)
      + f"   affinity {'on' if (PIN and HAVE_TASKSET) else 'off'}"
      + ("" if HAVE_TASKSET or not PIN else "  (taskset not found)"))
print(f"{'workers':>7} {'per-call ms':>12} {'calls/s':>10} {'speedup':>8}")

baseline = None
for workers in (1, 2, 4, 8, 16, 24, 32):
    if workers > len(CORES) and PIN and HAVE_TASKSET:
        print(f"{workers:>7}  skipped: more workers than cores ({len(CORES)})")
        continue
    env = dict(os.environ, **pin_env)
    handles = []
    for index in range(workers):
        cmd = [sys.executable, str(worker), str(calls)]
        if PIN and HAVE_TASKSET:
            cmd = ["taskset", "-c", str(CORES[index % len(CORES)])] + cmd
        handles.append(subprocess.Popen(cmd, env=env, stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE, text=True))
    results, errors = [], []
    for handle in handles:
        out, err = handle.communicate()
        found = [l for l in out.splitlines() if l.startswith("{")]
        if found:
            results.append(json.loads(found[-1]))
        else:
            errors.append((err or "")[-300:])
    if not results:
        print(f"{workers:>7}  all workers failed:\n    {errors[0] if errors else '?'}")
        break
    if errors:
        print(f"{workers:>7}  note: {len(errors)}/{workers} failed (likely memory); "
              f"row uses the {len(results)} that returned")
    slowest = max(r["seconds"] for r in results)
    rate = sum(r["calls"] for r in results) / slowest
    baseline = baseline or rate
    per_call = sum(r["seconds"] / r["calls"] for r in results) / len(results) * 1e3
    print(f"{workers:>7} {per_call:>12.1f} {rate:>10.1f} {rate / baseline:>7.2f}x")

worker.unlink(missing_ok=True)
print("\nPER_GPU = the largest concurrency still gaining, divided by the number of GPUs.")
print("A row whose speedup barely moves from the one above it is the ceiling.")
print("Compare THREADS=0 PIN=0 against the default to see how much of the ceiling was")
print("thread contention rather than the GPU.")
