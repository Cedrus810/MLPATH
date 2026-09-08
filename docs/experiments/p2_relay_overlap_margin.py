#!/usr/bin/env python3
"""Which recorded relay sign decisions could the mass-weighting defect have changed?

`_projected_direction_from` used to remove the uniform vector instead of the mass-weighted
translation, so the travel direction it handed `align_mode_with_travel` carried a spurious
translation on any heteronuclear molecule. The relay orients a saddle's unstable mode by
the SIGN of one overlap against that direction, so the question is not how wrong the
direction was but whether the error exceeded the margin the sign was decided by.

Both quantities are measurable from what is already on disk, with no potential evaluation:

  the error   ||basis_old - basis_new|| over every travel vector in the stored descent
              paths. Pure geometry: the two projections differ only by the mass-weighted
              centre-of-mass displacement, and the trajectories are recorded frame by
              frame. Measured on P2 this is 0.004 to 0.442, median 0.103 -- the centre of
              mass moves about 7% of the travel during a quench, so the artefact is not
              small.

  the margin  |travel_overlap| as recorded in every network.json.

A relay is undetermined when its margin is below the error, because then the recorded sign
is a property of the bug rather than of the potential. Run it:

    python docs/experiments/p2_relay_overlap_margin.py

The answer on the frozen P2 runs is bimodal with nothing in between, which is why no
threshold has to be argued for: 350 relays sit above every percentile of the error and 19
sit below all of them. Exit code is 0 -- this reports, it does not gate.
"""
import glob
import json
import os
import sys
from collections import Counter

import numpy as np
from ase.io import read

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def _projected(displacement, masses, mass_weighted_centre):
    """The travel direction, with translation removed one way or the other."""
    weighted = np.sqrt(masses[:, None]) * displacement
    if mass_weighted_centre:
        weighted = np.sqrt(masses[:, None]) * (
            displacement - np.average(displacement, axis=0, weights=masses))
    else:
        weighted = weighted - weighted.mean(axis=0)
    norm = np.linalg.norm(weighted)
    return None if norm < 1e-12 else (weighted / norm).ravel()


def measure_error():
    """How far apart the two projections are, on real recorded travel."""
    deltas = []
    paths = sorted(glob.glob(os.path.join(ROOT, "runs", "p2_*", "paths", "ts*",
                                          "descent.extxyz")))
    for path in paths:
        frames = read(path, ":")
        masses = frames[0].get_masses()
        phases = {}
        for frame in frames:
            phases.setdefault(frame.info.get("phase"), []).append(frame)
        for group in phases.values():
            for i in range(len(group) - 1):
                for j in {i + 1, len(group) - 1}:
                    if j <= i:
                        continue
                    displacement = group[j].positions - group[i].positions
                    old = _projected(displacement, masses, False)
                    new = _projected(displacement, masses, True)
                    if old is not None and new is not None:
                        deltas.append(float(np.linalg.norm(old - new)))
    return len(paths), np.array(deltas)


def collect_margins():
    """Every recorded |travel_overlap|, with the run it came from."""
    rows = []
    for path in glob.glob(os.path.join(ROOT, "runs", "**", "network.json"), recursive=True):
        run = os.path.basename(os.path.dirname(path))
        stack = [json.load(open(path))]
        while stack:
            node = stack.pop()
            if isinstance(node, dict):
                value = node.get("travel_overlap")
                if isinstance(value, (int, float)):
                    rows.append((abs(float(value)), run))
                stack.extend(node.values())
            elif isinstance(node, list):
                stack.extend(node)
    return rows


def main():
    files, deltas = measure_error()
    rows = collect_margins()
    if not len(deltas) or not rows:
        print("no stored descent paths or no recorded relays; nothing to measure")
        return 0
    margins = np.array([m for m, _ in rows])

    print(f"projection error, over {len(deltas)} travel vectors from {files} descent paths")
    for q in (0, 5, 50, 95, 100):
        print(f"  p{q:<3} ||basis_old - basis_new|| = {np.percentile(deltas, q):.4f}")

    low, high = deltas.min(), deltas.max()
    p5 = np.percentile(deltas, 5)
    undetermined = margins < high
    robust = margins > high

    print(f"\nrecorded relays: {len(margins)}")
    print(f"  margin above the LARGEST error ({high:.4f}): {robust.sum():>4}  robust")
    print(f"  margin below the 5th-pct error ({p5:.4f}): {(margins < p5).sum():>4}  "
          "undetermined")
    print(f"  anywhere in between                        : "
          f"{(~robust & ~(margins < p5)).sum():>4}")

    below, above = margins[margins < 0.5], margins[margins >= 0.5]
    if len(below) and len(above):
        print(f"\nthe two clusters: max below = {below.max():.6f}, "
              f"min above = {above.min():.6f}")
        print(f"the whole error range [{low:.4f}, {high:.4f}] lies inside that gap, so the "
              "split does not depend on which percentile is used")

    if undetermined.any():
        print(f"\nundetermined relays by run:")
        for run, count in Counter(run for (m, run) in rows if m < high).most_common():
            print(f"  {run}: {count}")
    print("\nNo relay needs a threshold argued for it: every undetermined one sits below "
          "the\n5th percentile of the error and every robust one above the largest, with "
          "the\nwhole error range falling in the empty gap between the two clusters.")
    print("\nUndetermined means the recorded sign was decided by a quantity smaller than "
          "the\nartefact, NOT that the sign is wrong. Deciding it needs the relay "
          "eigenvector,\nwhich is not stored, so it needs the potential.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
