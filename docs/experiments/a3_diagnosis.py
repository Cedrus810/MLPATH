"""Read the delivered path out of a run directory and say where it stopped.

Round three of the 3-oxobutanal benchmark had three transition states supporting the proton
transfer channel and zero trials in which A reached B. That is not a contradiction -- the
channel was named by descending from the saddles, which is undirected structural evidence --
but the record could not say why the delivered probe stopped, because it kept only the
endpoint classification. This reads the per-frame path data instead.

Usage: python a3_diagnosis.py runs/p2_round4
"""
import json
import sys
from pathlib import Path


def rows(path):
    return [json.loads(line) for line in path.read_text().splitlines()] if path.exists() else []


def summarise_track(frames, name):
    """Range and turning point of one tracked coordinate over a stage."""
    values = [(row["frame"], row["phase"], row["tracked"].get(name)) for row in frames
              if row.get("tracked", {}).get(name) is not None]
    if not values:
        return None
    numbers = [v for _, _, v in values]
    lowest = min(values, key=lambda item: item[2])
    highest = max(values, key=lambda item: item[2])
    return {"first": numbers[0], "last": numbers[-1],
            "min": lowest[2], "min_phase": lowest[1],
            "max": highest[2], "max_phase": highest[1], "samples": len(numbers)}


def main(run):
    run = Path(run)
    network = json.loads((run / "network.json").read_text())
    config = json.loads((run / "config.json").read_text())["config"]
    tracked = [spec["name"] for spec in config.get("tracked_coordinates") or ()]
    print(f"=== {run} ===")
    print(f"status {network['status']}  attempts {network['attempts']}  "
          f"nodes {len(network['chemical_nodes'])}  "
          f"reactions {len(network['reactions'])}  "
          f"channels {len(network['reaction_channels'])}  "
          f"ts {len(network['ts_candidates'])}")
    print(f"tracked coordinates: {tracked or 'NONE DECLARED'}")

    print("\n--- directed edges, by mechanism ---")
    if not network["reactions"]:
        print("  none. A3 as written (len(reactions) >= 1) fails.")
    for edge in network["reactions"]:
        print(f"  {edge['id']}  {edge['source']} -> {edge['target']}  "
              f"mechanisms={edge.get('mechanisms')}")
        print(f"      attempts (trial observed it itself): {edge['attempts']}")
        print(f"      continuations (probe advanced past the switch): "
              f"{edge.get('continuations')}")

    print("\n--- boundary continuation, rung by rung ---")
    any_continuation = False
    for scan in network.get("amplitude_scans", []):
        report = scan.get("boundary_continuation")
        if not report:
            continue
        if not report.get("attempted"):
            print(f"  {scan['family']} {scan['indices']}: not attempted "
                  f"({report.get('reason')})")
            continue
        any_continuation = True
        print(f"  {scan['family']} {scan['indices']} from "
              f"{scan['source_microstate']}  reached={report.get('reached')}")
        for seed in report["seeds"]:
            print(f"    seed {seed['seed_attempt']}  frame={seed.get('seed_frame')}  "
                  f"tangent={seed.get('direction_source')}  "
                  f"amplitude={seed.get('amplitude')}  "
                  f"transient={seed.get('transient')}")
            if seed.get("reason"):
                print(f"      refused: {seed['reason']}")
                continue
            walk = seed.get("continuation") or {}
            print(f"      outcome={walk.get('outcome')}  step={walk.get('step_A')} A")
            for rung in walk.get("rungs", []):
                extra = ""
                if rung.get("saddle_order") is not None:
                    extra += f" order={rung['saddle_order']}"
                if rung.get("gate"):
                    extra += f" gate={rung['gate']}"
                print(f"        rung {rung['rung']} at {rung['displacement_A']:.3f} A: "
                      f"{rung['verdict']}{extra}"
                      + (f"  E={rung['energy_eV']:.4f}" if "energy_eV" in rung else ""))
            if seed.get("target"):
                print(f"      admitted as {seed['target']} / "
                      f"{seed.get('target_microstate')}  "
                      f"class={seed.get('transition_class')}  "
                      f"channel={seed.get('reaction_channel')}")
    if not any_continuation:
        print("  no scan produced a continuation attempt")

    print("\n--- what the tracked coordinate did on each recorded stage ---")
    stages = sorted((run / "paths").glob("*/*.jsonl")) if (run / "paths").exists() else []
    if not stages:
        print("  no path records (paths/ is empty)")
    for path in stages:
        frames = rows(path)
        label = f"{path.parent.name}/{path.stem}"
        for name in tracked[:1]:
            summary = summarise_track(frames, name)
            if summary is None:
                continue
            print(f"  {label:<28} {name}: {summary['first']:+.3f} -> "
                  f"{summary['last']:+.3f}   "
                  f"[{summary['min']:+.3f} .. {summary['max']:+.3f}]  "
                  f"{summary['samples']} frames")

    print("\n--- trials whose delivered probe crossed the topology ---")
    for line in (run / "attempts.jsonl").read_text().splitlines():
        record = json.loads(line)
        crossing = (record.get("crossing") or {}).get("observed")
        if not (crossing or record["status"] == "ts_candidate"):
            continue
        frames = rows(run / "trials" / record["attempt_id"] / "observations.jsonl")
        summary = summarise_track(frames, tracked[0]) if tracked else None
        print(f"  {record['attempt_id']}  {record['status']:<14} "
              f"amp={record['probe']['amplitude']}  "
              f"target={record.get('target_microstate')}"
              + (f"   {tracked[0]}: {summary['first']:+.3f} -> {summary['last']:+.3f} "
                 f"[{summary['min']:+.3f} .. {summary['max']:+.3f}]" if summary else ""))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "runs/p2_round4")
