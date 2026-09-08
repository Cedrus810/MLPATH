"""Where does the time go, and can batching absorb it?

The search is a long chain of single-point evaluations on a 12-atom molecule. That is the
worst possible shape for a GPU: almost none of the 35 ms per call is arithmetic. Before
proposing anything, measure the per-call floor and what a batched forward pass costs.
"""
import sys, time; sys.path.insert(0, "/home/ruigengji/MLPATH/src")
import numpy as np, torch
from ase.io import read
from mace.calculators import MACECalculator
from prrs.calculators import mace_factory

MODEL = "/home/kasuga/.cache/mace/MACE-OFF24_medium.model"
atoms = read("p2_input.extxyz")
calc = MACECalculator(model_paths=MODEL, device="cuda", default_dtype="float64")
atoms.calc = calc
atoms.get_potential_energy()

def timed(fn, reps):
    fn(); torch.cuda.synchronize(); t = time.perf_counter()
    for _ in range(reps): fn()
    torch.cuda.synchronize(); return (time.perf_counter() - t) / reps

rng = np.random.default_rng(0)
frames = []
for _ in range(64):
    a = atoms.copy(); a.set_positions(atoms.positions + rng.normal(scale=0.02, size=atoms.positions.shape))
    frames.append(a)

def single():
    a = frames[0]; a.calc = calc; calc.results.clear(); a.get_forces()
per_call = timed(single, 30)
print(f"single point, {len(atoms)} atoms      {per_call*1e3:7.2f} ms")

# batched forward: one graph containing N independent configurations
from mace import data
from mace.tools import torch_geometric, utils
model = calc.models[0]
z_table = utils.AtomicNumberTable([int(z) for z in model.atomic_numbers])
cutoff = float(model.r_max)

def batched(n):
    configs = [data.AtomicData.from_config(data.config_from_atoms(f), z_table=z_table,
                                          cutoff=cutoff) for f in frames[:n]]
    loader = torch_geometric.dataloader.DataLoader(configs, batch_size=n, shuffle=False)
    batch = next(iter(loader)).to("cuda")
    dtype = next(model.parameters()).dtype
    for key in batch.keys:
        value = batch[key]
        if torch.is_tensor(value) and torch.is_floating_point(value):
            batch[key] = value.to(dtype=dtype)
    bd = batch.to_dict()
    bd["positions"].requires_grad_(True)
    def run():
        out = model(bd, compute_force=True, compute_stress=False, training=False)
        return out["forces"]
    return timed(run, 20)

print(f"{'batch':>6} {'total ms':>10} {'ms/config':>11} {'speedup':>9}")
for n in (1, 4, 8, 16, 32, 64):
    total = batched(n)
    print(f"{n:>6} {total*1e3:>10.2f} {total/n*1e3:>11.3f} {per_call/(total/n):>8.1f}x")
print(f"\npeak CUDA mem {torch.cuda.max_memory_allocated()/2**20:.0f} MiB / 11264 MiB")
