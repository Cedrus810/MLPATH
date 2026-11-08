"""Measure: does rebuilding the POLAR-1 calculator on every factory() call dominate wall time?

Same b04 C2 certification as b04_certify.py, run twice in one process:
  fresh   counted_factory(polar1): every factory() call builds a new MACECalculator
  cached  one MACECalculator built once, instrumented once, returned by every call
Before both, three back-to-back polar1() builds are timed (load cost alone; the first
one includes the cold start). Results must agree (status, saddle energy) or the
comparison is void. Measurement only; nothing in explore/ or prrs/ is changed.
"""

import json
import socket
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "explore" / "src"), str(ROOT)]

from b04_certify import find_guess_frame  # noqa: E402

from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.certify import certify  # noqa: E402
from explore.ledger import Ledger, counted_factory, instrument  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402
from prrs.search import _stamped_with_charge  # noqa: E402

CEILING_RISE_EV = 1.0


def counting(factory):
    calls = {"n": 0}

    def make():
        calls["n"] += 1
        return factory()

    return make, calls


def run_arm(name, out, config, guess, make_factory):
    ledger = Ledger(out / f"{name}.jsonl", arm=name, budget=20000)
    factory, calls = counting(make_factory(ledger))
    atoms = guess.copy()
    atoms.calc = factory()
    with ledger.segment(task=name, stage="certify_prepare"):
        e0 = float(atoms.get_potential_energy())
    t0 = time.monotonic()
    cert = certify(
        atoms,
        factory,
        config,
        ledger,
        task=name,
        seed=0,
        ceiling_eV=e0 + CEILING_RISE_EV,
        reference_eV=e0,
    )
    wall = time.monotonic() - t0
    ledger.close()
    stages = [
        json.loads(x)
        for x in (out / f"{name}.jsonl").read_text().splitlines()
        if json.loads(x)["event"] == "close"
    ]
    return {
        "status": cert.status,
        "saddle_energy_eV": cert.saddle_energy_eV,
        "certify_wall_s": wall,
        "certify_spent": cert.spent,
        "factory_calls_total": calls["n"],
        "stages": [
            {k: s[k] for k in ("stage", "forwards", "wall_s", "outcome")} for s in stages
        ],
    }


def main(out_dir):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=False)
    record = {"host": socket.gethostname(), "loads_s": [], "arms": {}}
    try:
        for _ in range(3):
            t = time.monotonic()
            polar1()
            record["loads_s"].append(time.monotonic() - t)
        config = SearchConfig(charge_sensitive=True, total_charge=-1, multiplicity=1)
        guess = _stamped_with_charge(find_guess_frame(), config)
        record["arms"]["fresh"] = run_arm(
            "fresh", out, config, guess, lambda ledger: counted_factory(polar1, ledger)
        )

        def cached(ledger):
            calc = instrument(polar1(), ledger)
            return lambda: calc

        record["arms"]["cached"] = run_arm("cached", out, config, guess, cached)
        a, b = record["arms"]["fresh"], record["arms"]["cached"]
        record["same_result"] = (
            a["status"] == b["status"]
            and abs(a["saddle_energy_eV"] - b["saddle_energy_eV"]) < 1e-6
        )
    finally:
        (out / "factory_reuse.json").write_text(json.dumps(record, indent=2, default=str))


if __name__ == "__main__":
    main(sys.argv[1])
