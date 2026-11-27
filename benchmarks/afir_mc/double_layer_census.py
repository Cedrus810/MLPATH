"""What does the two-pair layer pay for and what does it find? (2026-11-26, plan for freezing
the rules; read-only, no model.) Each two-pair item n|k-l+m-q is classified by its two pairs:
II both inter-fragment, IX one inter + one intra, XX both intra (fragments = bond-graph
components of the system input). Counted: items run, forwards (wave ledgers), new events and
order-1 reference matches credited to the class (afir_validate ts_rows, layer double).

  python double_layer_census.py OUT_JSON RUN:SYSTEM[:VALIDATE_JSON:KEY] ...
"""

import json
import os
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / "explore/src"), str(ROOT / "src"), str(ROOT), str(Path(__file__).parent)]
from ase.io import read  # noqa: E402

from explore.afir_batch import edges_of  # noqa: E402

INPUTS = {
    "b04": "benchmarks/b04_sn2_chloride/input.extxyz",
    "b05": "benchmarks/b05_tertiary_halide/input.extxyz",
    "n1": "benchmarks/n1_diels_alder/input.extxyz",
    "b09": "benchmarks/carbene_ambimodal/results/2026-11-20/b09_afir_input/input.extxyz",
}


def comps(system):
    a = read(ROOT / INPUTS[system])
    parent = list(range(len(a)))

    def f(u):
        while parent[u] != u:
            u = parent[u]
        return u

    for i, j in edges_of(a.numbers, a.positions, 1.2):
        parent[f(i)] = f(j)
    return [f(u) for u in range(len(a))]


def cls(label, comp):
    pairs = label.split("|")[1].split("+")
    inter = sum(comp[int(p.split("-")[0])] != comp[int(p.split("-")[1])] for p in pairs)
    return {2: "II", 1: "IX", 0: "XX"}[inter]


def main():
    out = {}
    for spec in sys.argv[2:]:
        parts = spec.split(":")
        run, system = Path(parts[0]), parts[1]
        comp = comps(system)
        n, ev, fw = Counter(), Counter(), Counter()
        names = sorted({x for x in os.listdir(run)})
        for fn in names:
            if fn.startswith("dwave_") and fn.endswith(".json"):
                w = json.loads((run / fn).read_text())
                for p in w["paths"]:
                    c = cls(p["label"], comp)
                    n[c] += 1
                    ev[c] += len(p.get("new_events") or [])
                    fw[c] += p.get("forwards_path") or 0
        refs = Counter()
        if len(parts) >= 4:
            v = json.loads(Path(parts[2]).read_text())[parts[3]]
            for r in v.get("ts_rows") or []:
                if r.get("layer") == "double" and r.get("order") == 1 and (r.get("match") or r.get("mirror_match")):
                    refs[cls(r["label"].split("#")[0], comp)] += 1
        out[str(run)] = {"items": dict(n), "new_events": dict(ev), "path_forwards": dict(fw),
                         "order1_ref_matches": dict(refs)}  # fmt: skip
        print(run, out[str(run)], flush=True)
    Path(sys.argv[1]).write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()
