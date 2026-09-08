"""Does CUDA MPS help or hurt many small MACE processes sharing one GPU?

Two earlier attempts at this were wrong and are worth naming, because both failures are
about measurement rather than about MPS:

1. Comparing PRRS seeds launched as MPS clients against seeds already running without it is
   confounded by run phase. A seed's trial rate depends strongly on how far in it is -- the
   same non-MPS seeds went from 4.2 min/trial in their first hour to 21.6 min/trial in their
   second with no configuration change -- so seeds of different ages are incomparable
   whatever their MPS status. That A/B suggested MPS was 2.7x faster. It measured age.
2. Running the two arms back to back on a shared GPU lets background load drift between
   them, and the first version also sent worker stderr to /dev/null, so when an arm returned
   nothing there was no way to see why.

So: identical fixed workload, arms INTERLEAVED (A B A B ...) so a drifting background hits
both equally, stderr captured, and the background load recorded per repeat. Two workloads,
because they answer different questions:

  kernel   2000 trivial `x = x + 1` launches -- pure per-launch cost
  mace     N single-point energy+force evaluations on the P2 structure -- the real thing

The only difference between arms is whether the worker can find the MPS server. Pointing
CUDA_MPS_PIPE_DIRECTORY at a path with no control pipe makes the client build its own
context; verified separately that this does not break CUDA.
"""
import json
import os
import statistics
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PYTHON = "/home/ruigengji/miniforge3/envs/openmm_dev/bin/python"
NO_MPS = {"CUDA_MPS_PIPE_DIRECTORY": "/tmp/claude-1000/mps-absent"}

KERNEL = '''
import json, time, torch
x = torch.zeros(8, device="cuda")
torch.cuda.synchronize()
start = time.perf_counter()
for _ in range(2000):
    x = x + 1.0
torch.cuda.synchronize()
print(json.dumps({"calls": 2000, "seconds": time.perf_counter() - start}))
'''

MACE = '''
import json, sys, time
from pathlib import Path
sys.path.insert(0, "%s/src")
from ase.io import read
from prrs.calculators import load_factory
meta = json.loads(Path("%s/runs/p2_round6/config.json").read_text())["calculator"]
atoms = read("%s/runs/p2_oxobutanal/p2_input.extxyz")
atoms.calc = load_factory(meta["factory"], meta["kwargs"])()
atoms.get_potential_energy()
n = int(sys.argv[1])
start = time.perf_counter()
for _ in range(n):
    atoms.positions[0, 0] += 1e-6
    atoms.get_potential_energy()
    atoms.get_forces()
print(json.dumps({"calls": n, "seconds": time.perf_counter() - start}))
''' % (ROOT, ROOT, ROOT)


def background():
    gpu = subprocess.run(["nvidia-smi", "--query-gpu=memory.used,utilization.gpu",
                          "--format=csv,noheader,nounits"],
                         capture_output=True, text=True).stdout.strip().replace("\n", " ")
    procs = subprocess.run(["bash", "-c", "ps aux | grep -c '[p]rrs search'"],
                           capture_output=True, text=True).stdout.strip()
    return f"gpu {gpu} MiB/%, searches {procs}"


def run_arm(script, args, workers, env_extra):
    env = dict(os.environ, **env_extra)
    handles = [subprocess.Popen([PYTHON, str(script), *args], env=env, text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
               for _ in range(workers)]
    started = time.perf_counter()
    seconds, failures = [], []
    for handle in handles:
        out, err = handle.communicate()
        lines = [l for l in out.splitlines() if l.startswith("{")]
        if lines:
            seconds.append(json.loads(lines[-1]))
        else:
            failures.append((err or "")[-400:])
    return {"wall_s": time.perf_counter() - started, "results": seconds,
            "failures": failures, "background": background()}


def report(name, workers, repeats, source, args):
    script = ROOT / f"docs/experiments/_mps_{name}.py"
    script.write_text(source, encoding="utf-8")
    arms = {"no_mps": [], "mps": []}
    for repeat in range(repeats):
        for label, env in (("no_mps", NO_MPS), ("mps", {})):
            row = run_arm(script, args, workers, env)
            arms[label].append(row)
            if row["failures"]:
                print(f"  [{name} {label} repeat {repeat}] "
                      f"{len(row['failures'])}/{workers} workers failed:")
                print("    " + row["failures"][0].replace("\n", "\n    "))
    print(f"\n=== {name}: {workers} workers, {repeats} interleaved repeats ===")
    summary = {}
    for label, rows in arms.items():
        per_call = [r["seconds"] / r["calls"] * 1e3
                    for row in rows for r in row["results"]]
        aggregate = [sum(r["calls"] for r in row["results"])
                     / max((r["seconds"] for r in row["results"]), default=float("nan"))
                     for row in rows if row["results"]]
        if not per_call:
            print(f"  {label}: no successful workers")
            continue
        summary[label] = {"per_call_ms": statistics.median(per_call),
                          "aggregate_calls_s": statistics.median(aggregate),
                          "n": len(per_call)}
        print(f"  {label:>7}  median per-call {statistics.median(per_call):8.3f} ms   "
              f"median aggregate {statistics.median(aggregate):8.1f} calls/s   "
              f"n={len(per_call)}")
        print(f"           backgrounds: " + " | ".join(r["background"] for r in rows))
    if len(summary) == 2:
        ratio = summary["mps"]["aggregate_calls_s"] / summary["no_mps"]["aggregate_calls_s"]
        verdict = "MPS faster" if ratio > 1.05 else (
            "MPS slower" if ratio < 0.95 else "no difference beyond noise")
        print(f"  --> MPS / no-MPS aggregate = {ratio:.2f}x   ({verdict})")
    return summary


if __name__ == "__main__":
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    calls = int(sys.argv[2]) if len(sys.argv) > 2 else 60
    repeats = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    out = {"workers": workers, "mace_calls": calls, "repeats": repeats,
           "kernel": report("kernel", workers, repeats, KERNEL, []),
           "mace": report("mace", workers, repeats, MACE, [str(calls)])}
    Path("runs/mps_ab.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    print("\nwrote runs/mps_ab.json")
