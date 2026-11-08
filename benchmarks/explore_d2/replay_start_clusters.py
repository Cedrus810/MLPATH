"""Replay D2-HF-v1.1 start clustering (explore.d2._start_clusters, the shipped code) on
finished D2-HF-v1 runs: cluster sizes, the refine forwards round one alone would have
spent, and whether every distinct event (symmetry-merged, hf_event_symmetry.json) keeps
a round-one member. Also times the clustering itself. Reads journals only; no model.

    python replay_start_clusters.py OUT_JSON CENSUS_DIR NOISE_MULTIPLE FIRST_ROUND RUN_DIR [...]
"""

import json
import sys
import time
from collections import Counter
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / "explore" / "src"), str(ROOT / "src")]

from ase import Atoms  # noqa: E402
from explore.d2 import (  # noqa: E402
    _CLASS_RANK,
    D2Settings,
    HFSettings,
    StartCluster,
    _classify,
    _key,
    _start_clusters,
    _task,
)
from prrs.config import SearchConfig  # noqa: E402


def _lines(path):
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.strip():
                yield json.loads(line)


def replay(run_dir, census_dir, multiple, first_round):
    run_dir = Path(run_dir)
    m = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    raw = dict(m["settings"])
    raw["magnitudes_A"] = tuple(raw["magnitudes_A"])
    poses = {}
    for r in _lines(run_dir / "poses.jsonl"):
        a = Atoms(numbers=r["numbers"], positions=r["positions"])
        a.info.update(r["info"])
        poses[r["pose_id"]] = a
    settled = {}
    for ev in _lines(run_dir / "journal.jsonl"):
        if ev.get("event") != "settle":
            continue
        if ev["task"][:3] in ("HS|", "HC|"):
            p = ev["payload"]
            ev["payload"] = {"duplicate_of": p.get("duplicate_of")}
        settled[ev["task"]] = ev
    run = SimpleNamespace(
        settled=settled,
        settings=D2Settings(**raw),
        config=SearchConfig.from_dict(m["config"]),
        hf=HFSettings(**m["hf"]),
        start_cluster=StartCluster(multiple, first_round),
        poses=poses,
        population=list(_lines(run_dir / "population.jsonl")),
    )
    chosen = []
    for row in run.population:
        if _task("L1", row) in settled:
            c = _classify(run, row)
            if c["class"] is not None:
                chosen.append((row, c))
    s = run.settings
    chosen.sort(key=lambda rc: (_CLASS_RANK[rc[1]["class"]], _key(s.seed, "promote", rc[0]["intent_id"])))
    t0 = time.perf_counter()
    clusters = _start_clusters(run, chosen)
    wall = time.perf_counter() - t0

    rep = {}
    for r in json.loads((Path(census_dir) / "hf_event_symmetry.json").read_text()):
        if Path(r["run"]).name == run_dir.name:
            ev, d = r["events"], r["pairwise_symmetric_rmsd_A"]
            for i, e in enumerate(ev):
                rep[e["task"]] = next(
                    (rep[ev[j]["task"]] for j in range(i) if d[i][j] <= r["basin_rmsd_A"]),
                    e["task"],
                )

    def cost(iid):
        return sum(settled[k]["spent"] for k in (f"HS|{iid}", f"HC|{iid}") if k in settled)

    def event(iid):
        hc = settled.get(f"HC|{iid}")
        if hc is None:
            return None
        if hc["status"] == "connected":
            return rep[hc["task"]]
        if hc["status"] == "hf_duplicate":
            return rep[hc["payload"]["duplicate_of"]]
        return None

    ran = [iid for iid in clusters if f"HS|{iid}" in settled]
    first = [iid for iid in ran if clusters[iid]["rank"] < first_round]
    total = sum(cost(i) for i in ran)
    events = {event(i) for i in ran} - {None}
    sizes = Counter(Counter(c["leader"] for c in clusters.values()).values())
    return {
        "run": str(run_dir),
        "radius_A": multiple * run.hf.geometry_noise_A,
        "first_round": first_round,
        "promoted": len(chosen),
        "clusters": len({c["leader"] for c in clusters.values()}),
        "cluster_sizes": dict(sorted(sizes.items())),
        "clustering_wall_s": wall,
        "hs_ran_in_v1": len(ran),
        "of_which_round_one": len(first),
        "refine_forwards_v1": total,
        "refine_forwards_round_one": sum(cost(i) for i in first),
        "round_two_fraction": 1 - sum(cost(i) for i in first) / total,
        "events_distinct": len(events),
        "events_with_round_one_member": len({event(i) for i in first} - {None}),
    }


def main():
    out, census, mult, k, runs = Path(sys.argv[1]), sys.argv[2], float(sys.argv[3]), int(sys.argv[4]), sys.argv[5:]
    results = [replay(r, census, mult, k) for r in runs]
    out.write_text(json.dumps(results, indent=1), encoding="utf-8")
    print(json.dumps(results, indent=1))


if __name__ == "__main__":
    main()
