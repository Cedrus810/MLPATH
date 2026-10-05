"""Stage 1 of the first-entry committor protocol: the confirmation batch.

PLAN_2026-09-25 steps A2-A3; docs/P1_COMMITTOR_FIRST_ENTRY_PROTOCOL.md (frozen
2026-09-25). Everything a verdict depends on is fixed there and copied here as
constants, not chosen here:

- the ridge frames of runs/p1_committor2, scans 0 / 2 / 6 / 8 (proton transfer);
- N = 50 fresh Maxwell-Boltzmann seeds per window, seed0_published + 100000 + i,
  disjoint from the 16 stage 0 used (seed0_published + i);
- every shot integrated for L_max_run = 1280 free steps and never stopped early, with
  sample_interval = 1 so every step's hysteretic graph and frame is on disk -- that is
  how the entry frame exists for the quench check, and why src/prrs/ is untouched;
- the frozen rule (k = 160, graph on A or B) applied afterwards, with k = 80 and 320
  for V2 and L_max = 640 for V1.

Shooting and judging are separate passes, so a kill costs only the shot in flight:
a shot directory without result.json is removed and redone.

    python docs/experiments/p1_first_entry_confirm.py --check     # frames only, no MD
    python docs/experiments/p1_first_entry_confirm.py             # shoot (hours)
    python docs/experiments/p1_first_entry_confirm.py --analyse   # V1-V5 per window
"""

import argparse
import hashlib
import json
import shutil
import sys
from dataclasses import replace
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np  # noqa: E402
from ase.io import read  # noqa: E402

from p1_first_entry_stage0 import entry, labels_of, reexits, rows_of  # noqa: E402
from p1_tfree_sensitivity import CALC, RUN, brackets  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402
from prrs.io import atomic_json, snapshot  # noqa: E402
from prrs.network import Registry  # noqa: E402
from prrs.reliability import GuardedCalculator  # noqa: E402
from prrs.runner import quench, shoot, wilson_interval  # noqa: E402
from prrs.search import _dividing_frame  # noqa: E402
from prrs.state import encode  # noqa: E402

PROTOCOL = ROOT / "docs" / "P1_COMMITTOR_FIRST_ENTRY_PROTOCOL.md"
PROTOCOL_SHA256_16 = "b685120fe25665a4"  # the frozen text; any edit refuses the run
OUT = ROOT / "runs" / "p1_first_entry_confirm"
SCANS = (0, 2, 6, 8)
N = 50
SEED_OFFSET = 100000
L_RUN = 1280
K = 160
K_V2 = (80, 160, 320)
L_V1 = 640
V2_MAX_FLIP, V3_MAX, V4_MAX, V5_MAX = 0.10, 0.10, 0.05, 0.20


def frozen_protocol():
    digest = hashlib.sha256(PROTOCOL.read_bytes()).hexdigest()
    if not digest.startswith(PROTOCOL_SHA256_16):
        raise SystemExit(f"protocol text changed since freezing: sha256 {digest[:16]}")
    return digest


def implementation_sha256():
    """Same algorithm as a run manifest (search.py): name + bytes of src/prrs/*.py."""
    digest = hashlib.sha256()
    for path in sorted((ROOT / "src" / "prrs").glob("*.py")):
        digest.update(path.name.encode())
        digest.update(path.read_bytes())
    return digest.hexdigest()


def setup():
    manifest = json.loads((RUN / "config.json").read_text(encoding="utf-8-sig"))
    config = SearchConfig.from_dict(manifest["config"])
    network = json.loads((RUN / "network.json").read_text(encoding="utf-8-sig"))
    registry = Registry.from_network(network, RUN, config)
    items = []
    for item in brackets(network, config):
        if item["scan"] not in SCANS:
            continue
        frame, origin, _, _ = _dividing_frame(
            RUN / "trials" / item["at_high"], probe=item["probe"].to_dict()
        )
        if frame is None:
            raise SystemExit(f"scan {item['scan']}: no ridge frame ({origin})")
        # the frame the published shots started from, compared without any dynamics
        published = read(str(RUN / "trials" / f"{item['prefix']}0000" / "response.traj"), 0)
        if not np.allclose(published.positions, frame.positions, atol=1e-10):
            raise SystemExit(f"scan {item['scan']}: ridge frame differs from the published one")
        item["frame"] = frame
        items.append(item)
    if [i["scan"] for i in items] != list(SCANS):
        raise SystemExit(f"windows found: {[i['scan'] for i in items]}")
    return manifest, config, registry, items


def shot_dir(item, i):
    return OUT / "trials" / f"{item['prefix']}{i:04d}x0000"


def run_shots(config, items):
    from prrs.calculators import mace_factory

    calculator = mace_factory(**json.loads(CALC.read_text()))
    cfg = replace(config, response_steps=L_RUN, sample_interval=1)
    for item in items:
        for i in range(N):
            directory = shot_dir(item, i)
            if (directory / "result.json").exists():
                continue
            if directory.exists():
                shutil.rmtree(directory)  # killed mid-shot: redo it whole
            # one shot per call so a kill loses one shot; the trial is named
            # <prefix><i>x0000 because shoot numbers from 0 within each call
            (outcome,) = shoot(
                item["frame"],
                item["probe"],
                lambda: calculator,
                cfg,
                OUT,
                1,
                item["seed0"] + SEED_OFFSET + i,
                prefix=f"{item['prefix']}{i:04d}x",
            )
            # run_trial does not write result.json (search.py's execute does); it is
            # this script's "shot finished" marker, written last
            atomic_json(directory / "result.json", outcome.record)
            print(f"scan {item['scan']:2d} shot {i + 1:2d}/{N} done", flush=True)


def analyse(config, registry, items):
    from prrs.calculators import mace_factory

    calculator = mace_factory(**json.loads(CALC.read_text()))
    states = {s.id: s for node in registry.nodes for s in node.microstates}

    def graph_of(state_id):
        return frozenset(
            encode(states[state_id].structure, config.bond_scale, active=config.active_atoms).edges
        )

    def quench_at(directory, free, step):
        row = next(r for r in free if r["step"] == step)
        atoms = read(str(directory / "response.traj"), index=row["frame"])
        guard = GuardedCalculator(calculator, config)
        atoms.calc = guard
        atoms.set_momenta(np.zeros(atoms.positions.shape))
        try:
            converged, _ = quench(atoms, guard, config)
        except Exception as exc:  # noqa: BLE001 - a gate rejection is a result here
            return f"gate_{type(exc).__name__}"
        if not converged:
            return "quench_unconverged"
        _, landed = registry.identify(snapshot(atoms, calculator))
        return landed.id if landed is not None else "not_in_registry"

    windows, v1_ok = [], True
    for item in items:
        graph_a, graph_b = graph_of(item["A"]), graph_of(item["B"])
        shots = []
        for i in range(N):
            directory = shot_dir(item, i)
            if not (directory / "result.json").exists():
                raise SystemExit(f"{directory.name}: not shot yet")
            free = rows_of(directory, L_RUN)
            if len({r["step"] for r in free if r["phase"] == "free"}) != L_RUN:
                raise SystemExit(f"{directory.name}: not every step was saved")
            labels = labels_of(free, graph_a, graph_b, L_RUN)
            shots.append((directory, free, labels))

        sides = {k: [entry(labels, k, L_RUN) for _, _, labels in shots] for k in K_V2}
        main = sides[K]
        short = [entry(labels, K, L_V1) for _, _, labels in shots]
        # V1: an entry by L_V1 is the same entry at L_RUN (a code check)
        v1 = all(s == m for s, m in zip(short, main) if s[0] is not None)
        v1 = v1 and sum(s[0] is None for s in short) >= sum(m[0] is None for m in main)
        v1_ok = v1_ok and v1

        per_shot, entered = [], []
        for (directory, free, labels), (side, at) in zip(shots, main):
            row = {"shot": directory.name, "side": side, "entry_step": at}
            if side is not None:
                row["reexits"] = reexits(labels, K, side, at)
                row["quench_landed"] = quench_at(directory, free, at)
                row["quench_agrees"] = row["quench_landed"] == (
                    item["A"] if side == "A" else item["B"]
                )
                at_row = next(r for r in free if r["step"] == at)
                row["q_PT"] = (at_row.get("tracked") or {}).get("q_PT")
                entered.append(row)
            per_shot.append(row)

        a = sum(r["side"] == "A" for r in per_shot)
        b = sum(r["side"] == "B" for r in per_shot)
        uncommitted = N - a - b
        flips = {
            f"{k1}-{k2}": sum(x[0] != y[0] for x, y in zip(sides[k1], sides[k2]))
            for k1, k2 in zip(K_V2, K_V2[1:])
        }
        intervals = {}
        for k, v in sides.items():
            a_k = sum(side == "A" for side, _ in v)
            b_k = sum(side == "B" for side, _ in v)
            intervals[str(k)] = list(wilson_interval(b_k, a_k + b_k)) if a_k + b_k else None
        reexit = sum(r["reexits"] for r in entered) / len(entered) if entered else None
        mismatch = sum(not r["quench_agrees"] for r in entered) / len(entered) if entered else None
        checks = {
            "V2": all(f <= V2_MAX_FLIP * N for f in flips.values()),
            "V3": reexit is not None and reexit <= V3_MAX,
            "V4": mismatch is not None and mismatch <= V4_MAX,
            "V5": uncommitted <= V5_MAX * N,
        }
        overlap = all(
            x is not None and y is not None and max(x[0], y[0]) <= min(x[1], y[1])
            for x, y in zip(list(intervals.values()), list(intervals.values())[1:])
        )
        windows.append(
            {
                "scan": item["scan"],
                "coordinate": item["coordinate"],
                "A": item["A"],
                "B": item["B"],
                "counts": {"A": a, "B": b, "uncommitted": uncommitted},
                "p_B_given_AB": b / (a + b) if a + b else None,
                "wilson_95": list(wilson_interval(b, a + b)) if a + b else None,
                "V1_window": v1,
                "V2_flips": flips,
                "V2_interval_overlap_reported": overlap,
                "V2_wilson_by_k": intervals,
                "V3_reexit_rate": reexit,
                "V4_quench_mismatch_rate": mismatch,
                "V4_landings": sorted({r["quench_landed"] for r in entered}),
                "V5_uncommitted_rate": uncommitted / N,
                "checks": checks,
                "q_PT_at_entry": {
                    s: [r["q_PT"] for r in entered if r["side"] == s] for s in ("A", "B")
                },
                "shots": per_shot,
            }
        )
    for window in windows:
        window["citable"] = v1_ok and all(window["checks"].values())
        print(
            f"scan {window['scan']:2d}: A {window['counts']['A']} / B {window['counts']['B']} / "
            f"uncommitted {window['counts']['uncommitted']}  p_B_given_AB "
            f"{window['p_B_given_AB']}  {window['wilson_95']}  checks {window['checks']}  "
            f"flips {window['V2_flips']}  citable {window['citable']}",
            flush=True,
        )
    print(f"V1 (whole batch): {v1_ok}", flush=True)
    return {"V1": v1_ok, "windows": windows}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true", help="frames and protocol only")
    parser.add_argument("--analyse", action="store_true", help="apply the frozen rule")
    args = parser.parse_args()

    protocol = frozen_protocol()
    manifest, config, registry, items = setup()
    for item in items:
        print(
            f"scan {item['scan']:2d} {item['coordinate']}: A {item['A']} B {item['B']}  "
            f"seeds {item['seed0'] + SEED_OFFSET}..{item['seed0'] + SEED_OFFSET + N - 1}",
            flush=True,
        )
    if args.check:
        return
    OUT.mkdir(parents=True, exist_ok=True)
    record = {
        "protocol": str(PROTOCOL.relative_to(ROOT)),
        "protocol_sha256": protocol,
        "source_run": str(RUN.relative_to(ROOT)),
        "source_fingerprint": manifest.get("implementation_sha256"),
        "implementation_sha256": implementation_sha256(),
        "N": N,
        "seed_offset": SEED_OFFSET,
        "L_max_run": L_RUN,
        "k": K,
        "k_V2": list(K_V2),
        "L_max_V1": L_V1,
        "thresholds": {"V2_flip": V2_MAX_FLIP, "V3": V3_MAX, "V4": V4_MAX, "V5": V5_MAX},
    }
    if not args.analyse:
        (OUT / "batch.json").write_text(json.dumps(record, indent=1) + "\n")
        run_shots(config, items)
        print("shots done", flush=True)
        return
    record.update(analyse(config, registry, items))
    (OUT / "verdict.json").write_text(json.dumps(record, indent=1) + "\n")
    print(f"wrote {OUT / 'verdict.json'}", flush=True)


if __name__ == "__main__":
    main()
