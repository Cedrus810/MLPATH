"""Stage 0 of the first-entry committor protocol: choose k and L_max, report no p_B.

PLAN_2026-09-25 step A0; docs/P1_COMMITTOR_FIRST_ENTRY_PROTOCOL.md section 3. Zero MD.
`runs/p1_tfree_sensitivity/L0640` already holds 16 shots per window integrated for 640
free steps (160 fs), and `observe()` wrote a row at every step whose hysteretic bond
graph changed -- so the graph at EVERY step is recoverable exactly, and this script
refuses to go on if the rows say otherwise (each row's edges must equal the previous
row's edges plus its own events; a change nobody saved would break that).

For each window, each k (consecutive steps on the A or B graph that count as entry)
and each L_max, it reports: A / B / uncommitted, the conditional p_B, how many shots
re-exit (spend k consecutive steps on the OTHER graph after entering), q_PT at entry,
and whether one quench from the first saved frame at or after entry lands on the
microstate the graph named. The parameters are chosen on these 64 shots, so nothing
here is citable as a p_B -- that is stage 1, on fresh seeds.

    python docs/experiments/p1_first_entry_stage0.py              # with quench checks
    python docs/experiments/p1_first_entry_stage0.py --no-quench  # graph analysis only
"""

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np  # noqa: E402
from ase.io import read  # noqa: E402

from p1_tfree_sensitivity import CALC, RUN, brackets  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402
from prrs.io import snapshot  # noqa: E402
from prrs.network import Registry  # noqa: E402
from prrs.reliability import GuardedCalculator  # noqa: E402
from prrs.runner import quench, wilson_interval  # noqa: E402
from prrs.state import encode  # noqa: E402

SOURCE = ROOT / "runs" / "p1_tfree_sensitivity"
LENGTH = 640
OUT = ROOT / "runs" / "p1_first_entry_stage0"
KS = (10, 20, 40, 80, 120, 160)
LMAXES = (160, 320, 640)
SCANS = (0, 2, 6, 8)  # proton transfer only; scan 15 (torsion) keeps one graph


def rows_of(directory, length=LENGTH):
    """The free-segment rows, with the exactness gate applied."""
    rows = [json.loads(line) for line in (directory / "observations.jsonl").open()]
    free = [r for r in rows if r["phase"] in ("free_start", "free")]
    previous = None
    for row in free:
        edges = frozenset(tuple(e) for e in row["edges"])
        if previous is not None:
            expected = (previous - {tuple(e) for e in row["events"]["broken"]}) | {
                tuple(e) for e in row["events"]["formed"]
            }
            if edges != expected:
                raise SystemExit(
                    f"{directory.name} step {row['step']}: edges do not follow from the "
                    "previous row and this row's events -- a graph change was not saved"
                )
        previous = edges
        row["_edges"] = edges
    if free[-1]["step"] != length:
        raise SystemExit(f"{directory.name}: free segment ends at {free[-1]['step']}")
    return free


def labels_of(free, graph_a, graph_b, length=LENGTH):
    """Label of every step 1..length: 'A', 'B' or None, from the last saved row."""
    by_step = {row["step"]: row["_edges"] for row in free}
    labels, current = [None], None
    for step in range(0, length + 1):
        if step in by_step:
            current = by_step[step]
        if step >= 1:
            labels.append("A" if current == graph_a else "B" if current == graph_b else None)
    return labels  # labels[s] for s = 1..length; labels[0] unused


def entry(labels, k, limit):
    run, last = 0, None
    for step in range(1, limit + 1):
        label = labels[step]
        run = run + 1 if (label is not None and label == last) else (1 if label else 0)
        last = label
        if run >= k:
            return label, step
    return None, None


def reexits(labels, k, side, at):
    other = "B" if side == "A" else "A"
    run = 0
    for step in range(at + 1, len(labels)):
        run = run + 1 if labels[step] == other else 0
        if run >= k:
            return True
    return False


def q_pt_at(free, step):
    row = max((r for r in free if r["step"] <= step), key=lambda r: r["step"])
    return (row.get("tracked") or {}).get("q_PT"), step - row["step"]


def summarise(values):
    values = [v for v in values if v is not None]
    if not values:
        return None
    return {"n": len(values), "min": min(values), "max": max(values)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--no-quench", action="store_true")
    args = parser.parse_args()

    manifest = json.loads((RUN / "config.json").read_text(encoding="utf-8-sig"))
    config = SearchConfig.from_dict(manifest["config"])
    network = json.loads((RUN / "network.json").read_text(encoding="utf-8-sig"))
    registry = Registry.from_network(network, RUN, config)
    states = {s.id: s for node in registry.nodes for s in node.microstates}
    fixed = json.loads((SOURCE / "summary.json").read_text())["lengths"]
    fixed = {str(L): {r["scan"]: r["verdicts"] for r in rows} for L, rows in fixed.items()}

    calculator = None
    if not args.no_quench:
        from prrs.calculators import mace_factory

        calculator = mace_factory(**json.loads(CALC.read_text()))

    def graph_of(state_id):
        return frozenset(
            encode(states[state_id].structure, config.bond_scale, active=config.active_atoms).edges
        )

    quenched = {}

    def confirm(directory, free, step, side, item):
        """One quench from the first saved frame at or after entry."""
        row = min((r for r in free if r["step"] >= step), key=lambda r: r["step"])
        key = (directory.name, row["frame"])
        if key not in quenched:
            atoms = read(str(directory / "response.traj"), index=row["frame"])
            guard = GuardedCalculator(calculator, config)
            atoms.calc = guard
            atoms.set_momenta(np.zeros(atoms.positions.shape))
            try:
                converged, _ = quench(atoms, guard, config)
            except Exception as exc:  # noqa: BLE001 - a gate rejection is a result here
                quenched[key] = f"gate_{type(exc).__name__}"
            else:
                if not converged:
                    quenched[key] = "quench_unconverged"
                else:
                    _, landed = registry.identify(snapshot(atoms, calculator))
                    quenched[key] = landed.id if landed is not None else "not_in_registry"
        landed = quenched[key]
        expected = item["A"] if side == "A" else item["B"]
        return {"frame_step": row["step"], "landed": landed, "agrees": landed == expected}

    summary = {
        "source": str((SOURCE / f"L{LENGTH:04d}").relative_to(ROOT)),
        "meaning": (
            "stage 0 of P1_COMMITTOR_FIRST_ENTRY_PROTOCOL: parameter choice only, no "
            "citable p_B -- chosen and read on the same 64 shots"
        ),
        "bond_scale": config.bond_scale,
        "bond_hysteresis": config.bond_hysteresis,
        "timestep_fs": config.timestep_fs,
        "windows": [],
    }
    for item in brackets(network, config):
        if item["scan"] not in SCANS:
            continue
        graph_a, graph_b = graph_of(item["A"]), graph_of(item["B"])
        if graph_a == graph_b:
            raise SystemExit(f"scan {item['scan']}: A and B share one graph")
        shots = []
        for i in range(item["published"]["shots"]):
            directory = SOURCE / f"L{LENGTH:04d}" / "trials" / f"{item['prefix']}{i:04d}"
            free = rows_of(directory)
            shots.append((directory, free, labels_of(free, graph_a, graph_b)))
        window = {
            "scan": item["scan"],
            "coordinate": item["coordinate"],
            "A": item["A"],
            "B": item["B"],
            "fixed_tfree_p_B": {},
            "by_k": {},
        }
        for L, verdicts in fixed.items():
            v = verdicts[item["scan"]]
            a, b = v.count("A"), v.count("B")
            window["fixed_tfree_p_B"][L] = b / (a + b) if a + b else None
        for k in KS:
            per_k = {"by_L_max": {}}
            entries = [entry(labels, k, LENGTH) for _, _, labels in shots]
            per_shot = []
            for (directory, free, labels), (side, at) in zip(shots, entries):
                row = {"shot": directory.name, "side": side, "entry_step": at}
                if side is not None:
                    row["q_PT"], row["q_PT_row_gap"] = q_pt_at(free, at)
                    row["reexits"] = reexits(labels, k, side, at)
                    row["horizon_steps"] = LENGTH - at
                    if calculator is not None:
                        row["quench"] = confirm(directory, free, at, side, item)
                per_shot.append(row)
            per_k["shots"] = per_shot
            entered = [r for r in per_shot if r["side"]]
            per_k["reexit"] = sum(r["reexits"] for r in entered)
            per_k["entered"] = len(entered)
            if calculator is not None:
                per_k["quench_disagrees"] = sum(not r["quench"]["agrees"] for r in entered)
            per_k["q_PT_at_entry"] = {
                side: summarise([r["q_PT"] for r in entered if r["side"] == side])
                for side in ("A", "B")
            }
            for limit in LMAXES:
                sides = [r["side"] if r["side"] and r["entry_step"] <= limit else None for r in per_shot]
                a, b = sides.count("A"), sides.count("B")
                per_k["by_L_max"][str(limit)] = {
                    "A": a,
                    "B": b,
                    "uncommitted": len(sides) - a - b,
                    "p_B_given_AB": b / (a + b) if a + b else None,
                    "wilson_95": list(wilson_interval(b, a + b)) if a + b else None,
                }
            window["by_k"][str(k)] = per_k
            # per-shot flips against the previous k, at the longest L_max
            if len(window["by_k"]) > 1:
                earlier = list(window["by_k"].values())[-2]["shots"]
                per_k["flips_vs_previous_k"] = sum(
                    a["side"] != b["side"] for a, b in zip(earlier, per_shot)
                )
            last = per_k["by_L_max"][str(LENGTH)]
            print(
                f"scan {item['scan']:2d} k={k:3d} ({k * config.timestep_fs:5.1f} fs): "
                + "  ".join(
                    f"L{lim}: {per_k['by_L_max'][str(lim)]['A']}/"
                    f"{per_k['by_L_max'][str(lim)]['B']}/"
                    f"{per_k['by_L_max'][str(lim)]['uncommitted']}"
                    for lim in LMAXES
                )
                + f"  p_B@640 {last['p_B_given_AB']}"
                + f"  reexit {per_k['reexit']}/{per_k['entered']}"
                + f"  flips {per_k.get('flips_vs_previous_k', '-')}"
                + f"  entry-step median {sorted(r['entry_step'] for r in entered)[len(entered) // 2] if entered else '-'}"
                + f"  q_PT A {per_k['q_PT_at_entry']['A']} B {per_k['q_PT_at_entry']['B']}"
                + (
                    f"  quench-disagree {per_k['quench_disagrees']}"
                    if calculator is not None
                    else ""
                ),
                flush=True,
            )
        summary["windows"].append(window)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(f"wrote {OUT / 'summary.json'}", flush=True)


if __name__ == "__main__":
    main()
