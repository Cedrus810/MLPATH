"""Branch C3 of PLAN_2026-09-25, exploration batch: p_B along each P1 saddle's path.

REVIEW_2026-09-24 section 4 / decisions R3-R5: shoot from several frames of ONE actual
path through the saddle -- the two-sided descent of each ts candidate, not frames pooled
across perturbation directions -- with the frozen first-entry rule (k = 160, cap 1280,
docs/P1_COMMITTOR_FIRST_ENTRY_PROTOCOL.md) through the production path
(`shoot(first_entry=...)`, implementation b601b124).

Exploration only (R5: N = 16 per frame). It chooses where a fixed N = 50 batch should
shoot; its numbers are not citable, and the fixed batch needs its own protocol frozen
first. "Other" is reported, never in the denominator (R4).

Frames per saddle, ordered by the proton-transfer coordinate q_PT:
  the saddle itself (ts_candidates/<id>.extxyz);
  saddle + f * delta for f in +-1/3, +-2/3, +-1, where delta is the first descent step
  (|q_PT| ~ 0.30). The descent does not record the stretch between the saddle and its
  first step, and that stretch is where p_B moves (entries in the confirmation batch all
  had |q_PT| >= 0.24); a straight segment of one mode step is the unstable mode itself,
  to first order.

Why mirrored: in runs/p1_committor2 (written before 99edddde) descent.extxyz holds only
the minus branch. search.py recorded the two
branches with two PathRecorders on one directory, each counting frames from 0, so the
second branch's first write (append=False) overwrote the first's geometries; the jsonl
kept both branches' rows. The plus side is therefore saddle - delta, and its q_PT is
printed next to the logged plus-branch q_PT so the mirror can be checked.

    python docs/experiments/p1_path_committor_explore.py
"""

import json
import sys
from dataclasses import replace
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from ase.io import read  # noqa: E402

from p1_tfree_sensitivity import CALC, RUN  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402
from prrs.internal import named_values  # noqa: E402
from prrs.network import Registry  # noqa: E402
from prrs.perturbations import Probe  # noqa: E402
from prrs.runner import shoot, wilson_interval  # noqa: E402
from prrs.state import encode  # noqa: E402

OUT = ROOT / "runs" / "p1_path_committor_explore"
N = 16
K = 160
L_MAX = 1280
SEED0 = 200000  # disjoint from every batch before (committor2, stage 0, confirmation)


def frames_of(ts, config):
    """(label, atoms) along the unstable mode, saddle in the middle, ordered by q_PT."""
    saddle = read(str(RUN / ts["structure"]))
    rows = [json.loads(line) for line in (RUN / "paths" / ts["id"] / "descent.jsonl").open()]
    path = read(str(RUN / "paths" / ts["id"] / "descent.extxyz"), index=":")
    if {f.info.get("phase") for f in path} != {"descent_minus"}:
        raise SystemExit(f"{ts['id']}: descent.extxyz is not the minus branch alone")
    minus = min((r for r in rows if r["sign"] == -1), key=lambda r: (r["step"], r["frame"]))
    plus = min((r for r in rows if r["sign"] == 1), key=lambda r: (r["step"], r["frame"]))
    delta = path[minus["frame"]].positions - saddle.positions
    frames = []
    for f in (-1.0, -2 / 3, -1 / 3, 0.0, 1 / 3, 2 / 3, 1.0):
        atoms = saddle.copy()
        atoms.positions = saddle.positions + f * delta
        frames.append((f"saddle{f:+.2f}delta", atoms))
    q = {
        label: named_values(atoms.positions, config.tracked_coordinates).get("q_PT")
        for label, atoms in frames
    }
    # saddle + delta IS the minus branch's first frame; saddle - delta stands in for the
    # lost plus branch, so it should land near the plus branch's logged q_PT
    print(
        f"{ts['id']}: q_PT minus logged {minus['tracked'].get('q_PT'):+.3f} vs saddle+delta "
        f"{q['saddle+1.00delta']:+.3f}; plus logged {plus['tracked'].get('q_PT'):+.3f} vs "
        f"mirror saddle-delta {q['saddle-1.00delta']:+.3f}",
        flush=True,
    )
    return sorted(frames, key=lambda f: q[f[0]]), q


def main():
    if OUT.exists():
        raise SystemExit(
            f"{OUT} exists; remove it to re-shoot (trial directories are never reused)"
        )
    manifest = json.loads((RUN / "config.json").read_text(encoding="utf-8-sig"))
    config = SearchConfig.from_dict(manifest["config"])
    network = json.loads((RUN / "network.json").read_text(encoding="utf-8-sig"))
    registry = Registry.from_network(network, RUN, config)
    states = {s.id: s for node in registry.nodes for s in node.microstates}
    from prrs.calculators import mace_factory

    calculator = mace_factory(**json.loads(CALC.read_text()))
    cfg = replace(config, response_steps=L_MAX)
    probe = Probe(family="stretch", indices=(4, 5), sign=1, amplitude=0.1, seed=0)  # label only
    summary = {
        "source_run": str(RUN.relative_to(ROOT)),
        "meaning": "exploration batch (N = 16), first-entry rule k = 160 / cap 1280; not citable",
        "saddles": [],
    }
    for t, ts in enumerate(network["ts_candidates"]):
        ends = {end["sign"]: end.get("microstate") for end in ts.get("connects") or []}
        if None in ends.values() or len(set(ends.values())) != 2:
            print(
                f"{ts['id']}: descent ends {ends} -- skipped, needs two distinct ends",
                flush=True,
            )
            continue
        a_id, b_id = sorted(ends.values())
        graphs = {
            side: frozenset(
                encode(states[i].structure, cfg.bond_scale, active=cfg.active_atoms).edges
            )
            for side, i in (("A", a_id), ("B", b_id))
        }
        if graphs["A"] == graphs["B"]:
            print(f"{ts['id']}: both ends share one graph -- skipped", flush=True)
            continue
        frames, q = frames_of(ts, cfg)
        rows = []
        for f, (label, atoms) in enumerate(frames):
            atoms.calc = calculator
            atoms.get_potential_energy()
            outcomes = shoot(
                atoms,
                probe,
                lambda: calculator,
                cfg,
                OUT / ts["id"],
                N,
                SEED0 + 100000 * t + 1000 * f,
                prefix=f"f{f:02d}s",
                first_entry={"graphs": graphs, "k": K},
            )
            sides = [(o.record.get("first_entry") or {}).get("side") for o in outcomes]
            a, b = sides.count("A"), sides.count("B")
            row = {
                "frame": label,
                "q_PT": q[label],
                "energy_eV": float(atoms.get_potential_energy()),
                "A": a,
                "B": b,
                "other": N - a - b,
                "p_B_given_AB": b / (a + b) if a + b else None,
                "wilson_95": list(wilson_interval(b, a + b)) if a + b else None,
            }
            rows.append(row)
            print(
                f"{ts['id']} {label:16s} q_PT {q[label]:+.3f}  A {a:2d} / B {b:2d} / other "
                f"{row['other']:2d}  p_B {row['p_B_given_AB']}",
                flush=True,
            )
        summary["saddles"].append({"ts": ts["id"], "A": a_id, "B": b_id, "frames": rows})
        OUT.mkdir(parents=True, exist_ok=True)
        (OUT / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print("done", flush=True)


if __name__ == "__main__":
    main()
