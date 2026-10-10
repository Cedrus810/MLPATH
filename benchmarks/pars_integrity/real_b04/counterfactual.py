"""Refine the b04 representative seeds as the pre-fix code did: no charge reaches POLAR-1."""
import dataclasses, json, sys
from pathlib import Path
from pars.cli import _candidate_from_seed
from pars.config import ParsConfig
from pars.refinement.ase_neb import AseNebBackend
from benchmarks.b04_sn2_chloride.calculator import polar1

selection = Path(sys.argv[1])
config = ParsConfig.from_dict(json.loads((selection / "pars_selection.json").read_text())["config"])
calc = polar1()
seen = set()
original = calc.calculate
def calculate(atoms=None, *a, **k):
    seen.add((atoms.info.get("charge"), atoms.info.get("spin")))
    return original(atoms, *a, **k)
calc.calculate = calculate
out = {}
for seed in sorted((selection / "representatives").glob("rep-*.extxyz")):
    candidate = _candidate_from_seed(seed)
    candidate.provenance = dataclasses.replace(candidate.provenance, electronic_state=None)
    result = AseNebBackend(calc).refine(candidate, config)
    p = result.provenance
    out[candidate.candidate_id] = {
        "status": result.status.value, "barrier_eV": result.barrier_eV,
        "message": result.message[:300],
        "negative_modes": p.get("saddle_check", {}).get("negative_modes"),
        "connectivity": p.get("connectivity_check", {}).get("status"),
    }
    print(candidate.candidate_id, out[candidate.candidate_id], flush=True)
out["states_seen"] = sorted(map(str, seen))
Path(sys.argv[2]).write_text(json.dumps(out, indent=1))
print(out["states_seen"])
