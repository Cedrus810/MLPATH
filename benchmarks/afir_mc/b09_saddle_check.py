"""Post-hoc check of the b09 saddle found by bond-window refinement (DONE_LOG 2026-11-25):
order and the two sides. Validation only (design §4.4: Hessian after the search, reference
only, never fed back). Real model, float32 POLAR-1.

  1. mass-weighted Hessian spectrum (prrs.runner.hessian_spectrum), internal modes only;
  2. displace +/- STEP A along the lowest mode (Cartesian), E/F minimize both sides
     (explore.afir_batch.minimize, the runner's quench fmax);
  3. name each side by the carbene-C bonds to Ca, Cb and the aldehyde O (bond_scale 1.2):
     none = reactant side, Cb only = open C:-Cb intermediate, Ca+Cb = (2+1), Cb+O = (4+1);
     any other bond change from the input complex is reported as "other".

  python b09_saddle_check.py REFINE_JSON LABEL K LO HI OUT_JSON
"""

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
from run_mc_afir_batch import identity  # noqa: E402,F401  (sets sys.path)
from run_d2 import make_atoms, make_config  # noqa: E402
from systems import SYSTEMS  # noqa: E402

from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.afir_batch import edges_of, minimize  # noqa: E402
from explore.ledger import Ledger  # noqa: E402
from explore.parallel import MACEBatchEvaluator  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402
from prrs.runner import _trivial_modes, _wavenumbers, hessian_spectrum  # noqa: E402

CARBENE, CA, CB, O_ALD = 2, 13, 15, 1  # TS-2OMe atom order (classify_b09 / build inputs)
STEP_A = 0.05
KEY = {(CARBENE, CA), (CARBENE, CB), (O_ALD, CARBENE)}


def name(edges, ref_edges):
    e = {tuple(sorted(x)) for x in edges}
    has = {k: tuple(sorted(k)) in e for k in KEY}
    a, b, o = has[(CARBENE, CA)], has[(CARBENE, CB)], has[(O_ALD, CARBENE)]
    base = {"none": not (a or b or o), "Cb": b and not a and not o,
            "(2+1)": a and b and not o, "(4+1)": b and o and not a}  # fmt: skip
    label = next((k for k, v in base.items() if v), "mixed")
    rest = (e ^ {tuple(sorted(x)) for x in ref_edges}) - {tuple(sorted(k)) for k in KEY}
    return label if not rest else f"{label}+other{sorted(rest)}"


def main():
    rj, label, k, lo, hi, out = sys.argv[1:7]
    rows = json.loads(Path(rj).read_text())
    row = next(r for r in rows if r["label"] == label and r["k"] == int(k)
               and r["window"] == [int(lo), int(hi)])  # fmt: skip
    system = SYSTEMS["b09"]
    config = make_config(system)
    atoms, _ = make_atoms(system, config)
    calc = polar1("cuda", "float32")
    ts = atoms.copy()
    ts.positions = np.asarray(row["ci_positions"])
    ts.info.update({"charge": 0, "spin": 1})
    ev, vec = hessian_spectrum(ts, lambda: calc, SearchConfig(), return_vectors=True)
    trivial = len(_trivial_modes(ts.positions, ts.get_masses()))
    keep = np.argsort(np.abs(ev))[trivial:]
    ev_int, vec_int = ev[keep], vec[:, keep]
    o = np.argsort(ev_int)
    ev_int, vec_int = ev_int[o], vec_int[:, o]
    n_neg = int((ev_int < 0).sum())
    mode = (vec_int[:, 0] / np.sqrt(np.repeat(ts.get_masses(), 3))).reshape(-1, 3)
    mode /= np.abs(np.linalg.norm(mode, axis=1)).max()
    starts = [ts.positions + STEP_A * mode, ts.positions - STEP_A * mode]
    led = Ledger(Path(out).with_suffix(".ledger.jsonl"), arm="check", replay=False)
    mins = minimize(atoms, starts, MACEBatchEvaluator(calc), config, led, "check|sides",
                    fmax=config.quench_fmax_eV_A, max_steps=3000)  # fmt: skip
    led.close()
    ref_edges = edges_of(atoms.numbers, atoms.positions, config.bond_scale)
    res = {"label": label, "k": int(k), "window": [int(lo), int(hi)],
           "ts_energy_eV": row["dE_vs_ref_eV"], "n_negative": n_neg,
           "lowest_wavenumbers_cm": [float(w) for w in _wavenumbers(ev_int[:4])],
           "sides": [{"converged": m.converged, "energy_eV": m.energy_eV,
                      "rel_to_ts_eV": None, "name": name(m.edges, ref_edges)} for m in mins]}  # fmt: skip
    print(json.dumps(res), flush=True)
    Path(out).write_text(json.dumps(res, indent=1) + "\n")


if __name__ == "__main__":
    main()
