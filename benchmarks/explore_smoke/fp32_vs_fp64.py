"""Measure: what does running POLAR-1 in float32 instead of float64 cost and buy?

`run DTYPE OUT` (one dtype per process: MACECalculator sets torch's global default dtype)
  1. timing   median wall per forward on b04 (6 atoms) and b09 (29 atoms), cache busted
  2. values   energy and forces on fixed geometry sets (b04 10-18 coarse + 0.01 A fine
              scans; b09 input + 10 seeded perturbations), saved for comparison
  3. repeat   the same b04 geometry evaluated by a fresh calculator: determinism
  4. certify  the b04 C2 certification from the 10-18 q = 2.28 A point
`compare OUT` reads both dtypes and writes compare.json.
Measurement only: nothing in explore/ or prrs/ is changed; no threshold is applied.
"""

import json
import socket
import statistics
import sys
import time
from pathlib import Path

import numpy as np
from ase.io import read

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "explore" / "src"), str(ROOT)]

from b04_certify import find_guess_frame, kabsch_rmsd  # noqa: E402

from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.certify import certify  # noqa: E402
from explore.ledger import Ledger, reused_counted_factory  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402
from prrs.search import _stamped_with_charge  # noqa: E402

B04_INPUT = ROOT / "benchmarks/b04_sn2_chloride/input.extxyz"
B09_INPUT = ROOT / "benchmarks/carbene_ambimodal/probe_b09/input.extxyz"
SCAN = ROOT / "benchmarks/explore_smoke/results/2026-10-18/run"
TS = ROOT / "runs/b04_probe2_openceiling/ts_candidates/ts0000.extxyz"
SYSTEMS = {"b04": -1, "b09": 0}
N_TIMED = 20


def config_for(charge):
    return SearchConfig(charge_sensitive=True, total_charge=charge, multiplicity=1)


def geometry_sets():
    b04 = [_stamped_with_charge(a, config_for(-1)) for a in read(B04_INPUT, index=":")]
    for name in ("forward", "fine_fwd"):
        for f in read(SCAN / f"{name}.extxyz", index=":"):
            b04.append(_stamped_with_charge(f, config_for(-1)))
    base = _stamped_with_charge(read(B09_INPUT), config_for(0))
    rng = np.random.default_rng(0)
    b09 = [base]
    for _ in range(10):
        a = base.copy()
        a.positions = a.positions + rng.normal(scale=0.05, size=a.positions.shape)
        b09.append(a)
    return {"b04": b04, "b09": b09}


def evaluate(calc, atoms):
    a = atoms.copy()
    a.calc = calc
    return float(a.get_potential_energy()), np.asarray(a.get_forces(), float)


def run(dtype, out):
    out.mkdir(parents=True, exist_ok=True)
    rec = {"dtype": dtype, "host": socket.gethostname(), "timing": {}, "repeat": None}
    t = time.monotonic()
    calc = polar1("cuda", dtype)
    rec["load_s"] = time.monotonic() - t
    sets = geometry_sets()

    for name, geoms in sets.items():
        for _ in range(3):
            evaluate(calc, geoms[0])
        walls = []
        for k in range(N_TIMED):
            a = geoms[0].copy()
            a.positions[0, 0] += 1e-4 * (k + 1)  # distinct positions: no cache hits
            t = time.monotonic()
            evaluate(calc, a)
            walls.append(time.monotonic() - t)
        rec["timing"][name] = {"n_atoms": len(geoms[0]), "median_s": statistics.median(walls)}

    values = {}
    for name, geoms in sets.items():
        es, fs = zip(*(evaluate(calc, g) for g in geoms))
        values[f"{name}_E"] = np.array(es)
        values[f"{name}_F"] = np.array(fs)
    np.savez(out / f"values_{dtype}.npz", **values)

    fresh = polar1("cuda", dtype)
    e1, f1 = evaluate(calc, sets["b04"][5])
    e2, f2 = evaluate(fresh, sets["b04"][5])
    rec["repeat"] = {"dE_eV": e2 - e1, "max_dF_eV_A": float(np.max(np.abs(f2 - f1)))}

    config = config_for(-1)
    guess = _stamped_with_charge(find_guess_frame(), config)
    ledger = Ledger(out / f"certify_{dtype}.jsonl", arm=f"fp-{dtype}", budget=20000)
    factory = reused_counted_factory(lambda: polar1("cuda", dtype), ledger)
    atoms = guess.copy()
    atoms.calc = factory()
    with ledger.segment(task="certify", stage="certify_prepare"):
        e0 = float(atoms.get_potential_energy())
    t = time.monotonic()
    cert = certify(
        atoms,
        factory,
        config,
        ledger,
        task="certify",
        seed=0,
        ceiling_eV=e0 + 1.0,
        reference_eV=e0,
    )
    wall = time.monotonic() - t
    ledger.close()
    ts = read(TS)
    rec["certify"] = {
        "status": cert.status,
        "wall_s": wall,
        "spent": cert.spent,
        "saddle_energy_eV": cert.saddle_energy_eV,
        "imaginary_wavenumbers_icm": (cert.hessian or {}).get("imaginary_wavenumbers_icm"),
        "lowest_eigenvalues": (cert.hessian or {}).get("lowest_eigenvalues"),
        "side_energies_eV": [s.get("energy_eV") for s in cert.sides],
        "rmsd_to_ts0000_A": None
        if cert.saddle_positions is None
        else kabsch_rmsd(np.asarray(cert.saddle_positions), ts.positions),
    }
    (out / f"run_{dtype}.json").write_text(json.dumps(rec, indent=2, default=str))


def compare(out):
    a, b = (json.loads((out / f"run_{d}.json").read_text()) for d in ("float64", "float32"))
    va, vb = (np.load(out / f"values_{d}.npz") for d in ("float64", "float32"))
    res = {
        "timing": {},
        "values": {},
        "certify": {"float64": a["certify"], "float32": b["certify"]},
    }
    for name in SYSTEMS:
        t64, t32 = a["timing"][name]["median_s"], b["timing"][name]["median_s"]
        res["timing"][name] = {
            "n_atoms": a["timing"][name]["n_atoms"],
            "fp64_s": t64,
            "fp32_s": t32,
            "speedup": t64 / t32,
        }
        dE = vb[f"{name}_E"] - va[f"{name}_E"]
        rel = (vb[f"{name}_E"] - vb[f"{name}_E"][0]) - (va[f"{name}_E"] - va[f"{name}_E"][0])
        dF = np.abs(vb[f"{name}_F"] - va[f"{name}_F"])
        res["values"][name] = {
            "n_geometries": int(len(dE)),
            "abs_energy_diff_eV": {
                "max": float(np.max(np.abs(dE))),
                "mean": float(np.mean(dE)),
            },
            "relative_energy_diff_eV_max": float(np.max(np.abs(rel))),
            "force_diff_eV_A": {"max": float(np.max(dF)), "median": float(np.median(dF))},
        }
    res["repeat"] = {"float64": a["repeat"], "float32": b["repeat"]}
    res["load_s"] = {"float64": a["load_s"], "float32": b["load_s"]}
    (out / "compare.json").write_text(json.dumps(res, indent=2, default=str))


if __name__ == "__main__":
    if sys.argv[1] == "run":
        run(sys.argv[2], Path(sys.argv[3]))
    else:
        compare(Path(sys.argv[2]))
