
import json, sys, time
from pathlib import Path
sys.path.insert(0, "/home/ruigengji/MLPATH/src")
from ase.io import read
from prrs.calculators import load_factory
meta = json.loads(Path("/home/ruigengji/MLPATH/runs/p2_round6/config.json").read_text())["calculator"]
atoms = read("/home/ruigengji/MLPATH/runs/p2_oxobutanal/p2_input.extxyz")
atoms.calc = load_factory(meta["factory"], meta["kwargs"])()
atoms.get_potential_energy()
n = int(sys.argv[1])
start = time.perf_counter()
for _ in range(n):
    atoms.positions[0, 0] += 1e-6
    atoms.get_potential_energy()
    atoms.get_forces()
print(json.dumps({"calls": n, "seconds": time.perf_counter() - start}))
