"""Aggregate the P2 multi-seed scan against the questions frozen in
docs/P2_MULTISEED_PROTOCOL.md.

The denominator is always the number of seeds LAUNCHED, never the number that finished. A
seed that crashed or is still running is reported as such and stays in the denominator: the
whole point of the scan is that one hit is not a hit rate.

Usage: python p2_multiseed.py 17:runs/p2_round6 19:runs/p2_seed19 ...
"""
import json
import sys
from pathlib import Path

TRUTH = ([[0, 5]], [[4, 5]])       # the pre-registered proton transfer, either orientation


def load(path):
    """A missing network.json is not the same thing as a broken run.

    "unreadable: FileNotFoundError" conflated three states that need different handling: a
    seed still queued, a seed that died before writing anything, and a seed whose file is
    corrupt. Only the last is a failure; the first is simply not an answer yet.
    """
    path = Path(path)
    network = path / "network.json"
    if not path.exists():
        return {"status": "not_started"}
    if not network.exists():
        log = Path("docs/experiments") / f"{path.name}.out"
        tail = ""
        if log.exists():
            lines = [l for l in log.read_text(errors="replace").splitlines() if l.strip()]
            tail = lines[-1][:120] if lines else ""
        return {"status": "started_but_no_network_json", "note": tail}
    try:
        return json.loads(network.read_text())
    except Exception as exc:                                  # noqa: BLE001
        return {"status": f"corrupt: {type(exc).__name__}"}


def truth_channel(network):
    for channel in network.get("reaction_channels", []):
        pair = ([sorted(b) for b in channel["broken"]],
                [sorted(f) for f in channel["formed"]])
        if {tuple(map(tuple, pair[0])), tuple(map(tuple, pair[1]))} == {
                tuple(map(tuple, TRUTH[0])), tuple(map(tuple, TRUTH[1]))}:
            return channel
    return None


def fmt_dE(row):
    return f"{row['dE_meV']:+.3f}" if row.get("dE_meV") is not None else "n/a"


def summarise(seed, path):
    network = load(path)
    row = {"seed": seed, "path": str(path), "status": network.get("status", "missing")}
    if row["status"] != "completed":
        return row
    nodes = network["chemical_nodes"]
    directed = [e for e in network["reactions"] if e["source"] != e["target"]]
    channel = truth_channel(network)
    ts = network["ts_candidates"]
    supporting = [] if channel is None else [t for t in ts
                                             if t["id"] in channel["ts_support"]]
    row.update(
        nodes=len(nodes),
        # Q1 needs the directed edge to be ON the transfer channel, not merely
        # some directed edge somewhere in the network.
        q1=bool(channel is not None
                and any({e["source"], e["target"]} == set(channel["ends"])
                        for e in directed)),
        q2=channel is not None,
        q4=bool([t for t in supporting if t["found_by"] == "min_mode_following"]),
        attempts_direct=sum(len(e["attempts"]) for e in directed),
        continuation_edges=sum(len(e.get("continuations") or []) for e in directed),
        mechanisms=sorted({m for e in directed for m in e.get("mechanisms", [])}),
        transfer_icm=sorted(round(v, 1) for t in supporting
                            if t["found_by"] == "min_mode_following"
                            for v in (t["imaginary_wavenumbers_icm"] or [])),
        trials=network["attempts"],
        outcomes=network["outcome_counts"],
    )
    energies = {n["id"]: n["microstates"][0]["energy_eV"] for n in nodes}
    if channel is not None and len(set(channel["ends"])) == 2:
        low, high = (energies.get(e) for e in channel["ends"])
        if low is not None and high is not None:
            row["dE_meV"] = round((high - low) * 1e3, 3)
    rungs, arrivals, stalled = 0, 0, 0
    attempted_walks = 0
    for scan in network.get("amplitude_scans", []):
        for seed_entry in (scan.get("boundary_continuation") or {}).get("seeds", []):
            walk = seed_entry.get("continuation")
            if walk is None:
                continue
            attempted_walks += 1
            rungs += len(walk["rungs"])
            arrivals += sum(1 for r in walk["rungs"]
                            if r["verdict"].startswith("reached"))
            stalled += sum(1 for r in walk["rungs"]
                           if r["verdict"] == "stopped_on_saddle")
            row.setdefault("via_saddle", 0)
            row["via_saddle"] += sum(1 for r in walk["rungs"]
                                     if r.get("arrival") == "via_saddle")
    row.update(walks=attempted_walks, rungs=rungs, rung_arrivals=arrivals,
               rung_stalled=stalled)
    return row


def main(argv):
    rows = []
    for item in argv:
        seed, _, path = item.partition(":")
        rows.append(summarise(int(seed), path))
    launched = len(rows)
    done = [r for r in rows if r["status"] == "completed"]
    unfinished = [r for r in rows if r["status"] not in ("completed",)]

    print(f"launched {launched} seeds, completed {len(done)}, unfinished {len(unfinished)}")
    if unfinished:
        print("  unfinished: " + ", ".join(f"{r['seed']}({r['status']})"
                                           for r in unfinished))
    print("\nseed  status     nodes  Q1 A3  Q2 truth  Q4 TS  dE(meV)   "
          "direct  cont  transfer icm")
    for row in rows:
        if row["status"] != "completed":
            print(f"{row['seed']:>4}  {row['status']}")
            continue
        print(f"{row['seed']:>4}  completed  {row['nodes']:>5}  "
              f"{'PASS' if row['q1'] else 'fail':>5}  "
              f"{'PASS' if row['q2'] else 'fail':>8}  "
              f"{'PASS' if row['q4'] else 'fail':>5}  "
              # "nan" reads as a computation that went wrong. The truth is that this seed
              # has no transfer channel, so there is no pair of endpoints to difference.
              f"{fmt_dE(row):>8}  "
              f"{row['attempts_direct']:>6}  {row['continuation_edges']:>4}  "
              f"{row['transfer_icm']}")

    print("\n--- the frozen questions ---")
    print("  The protocol denominator is seeds LAUNCHED. While seeds are still running that")
    print("  ratio is a LOWER BOUND, not a rate: an unfinished seed scores 0 because it has")
    print("  not answered, not because it failed. The second column is the interim rate over")
    print("  finished seeds only -- useful now, NOT the pre-registered number.")
    for tag, key, text in (("Q1", "q1", "A3: directed edge on the transfer channel"),
                           ("Q2", "q2", "truth-set recall: the transfer channel exists"),
                           ("Q4", "q4", "a min-mode TS supports that channel")):
        hits = sum(1 for r in rows if r.get(key))
        among = sum(1 for r in done if r.get(key))
        frozen = f"{hits}/{launched}"
        interim = f"{among}/{len(done)}" if done else "n/a"
        flag = "  <-- final" if not unfinished else "  (lower bound)"
        print(f"  {tag}  frozen {frozen:>7}{flag:<15} interim {interim:>7}   {text}")
    print("\n--- recorded, not criteria ---")
    if done:
        icm = sorted(v for r in done for v in r["transfer_icm"])
        des = [r.get("dE_meV") for r in done if r.get("dE_meV") is not None]
        print(f"  transfer icm across seeds: n={len(icm)} "
              + (f"min {icm[0]} max {icm[-1]}" if icm else "none"))
        print(f"  dE(B-A) across seeds: "
              + (f"min {min(des)} max {max(des)} meV" if des else "none"))
        print(f"  direct observations (attempts) total: "
              f"{sum(r['attempts_direct'] for r in done)}")
        print(f"  continuation walks {sum(r['walks'] for r in done)}, "
              f"rungs {sum(r['rungs'] for r in done)}, "
              f"arrivals {sum(r['rung_arrivals'] for r in done)}, "
              f"stalled {sum(r['rung_stalled'] for r in done)}, "
              f"via_saddle {sum(r.get('via_saddle', 0) for r in done)}")
        print(f"  trials per seed: {[r['trials'] for r in done]}")
    Path("runs/p2_multiseed.json").write_text(
        json.dumps({"launched": launched, "rows": rows}, indent=1), encoding="utf-8")
    print("\nwrote runs/p2_multiseed.json")


if __name__ == "__main__":
    main(sys.argv[1:] or ["17:runs/p2_round6"] +
         [f"{s}:runs/p2_seed{s}" for s in
          (19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79,
           83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151,
           157, 163)])
