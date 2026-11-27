"""MC-AFIR / SC-AFIR in batched waves (MLP-specific: many concurrent paths; design §0 and
explore.afir_batch). Real model, float32 POLAR-1.

    python run_mc_afir_batch.py --system NAME --out DIR --seed S --shell R --budget B
        [--mode mc|sc] [--sc-poses N] [--wave 32] [--max-paths N] [--n-stop-waves 2]
        [--gamma 106.9] [--max-active 16]

Preparation as run_mc_afir.py (fragments E/F-relaxed and perturbation-confirmed; random
orientations from poses.sample_poses).
  mc: one path per random orientation, force between the two whole fragments.
  sc: reactant-complex MINs = the input complex and the first --sc-poses random
      orientations, each E/F-minimized, duplicates merged (bond graph + symmetric RMSD +
      energy); on every MIN, one path per candidate atom pair (r/(R_i+R_j) <= 3.0, one per
      symmetry orbit), started from the SC-AFIR perturbed structure with the SC-AFIR
      fragments A / B around the pair (explore.afir_fragments). Push only (rho = +1).
      Push coordinate (--sc-coord): pair (default, user 2026-11-15) = the selected pair
      only, V_bias = alpha r_kl, alpha(gamma) as AFIR; frozen = q = sum c_ij r_ij over A x B
      with weights frozen at the perturbed start (design §2.1, now a control); rho = rho_AB
      (runs before 2026-11-15). TS coordinate (--ts-coord), decoupled from the push:
      frozen-peak (default) = q over the same A x B with weights frozen at the bracket's
      peak frame; rho = rho_AB over A x B; push = the push coordinate itself (coupled).
      Two-pair layer (--sc-double 1, default; user 2026-11-15), run after the single-pair
      base layer with its own waves / ledgers (dwave_*, dts_*) and event set, under the
      same --budget: on every MIN all two-contact combinations (two non-bonded contacts,
      r/(R_i+R_j) <= 3.0, four different atoms, one per bond-graph symmetry orbit, no
      ring-size filter; afir_fragments.two_contact_combos), shuffled with --seed; both
      distances shortened 20 % in one least-squares fit; push V_bias = alpha (r_kl +
      r_mn)/sqrt 2 (gamma unchanged); TS coordinate q over the pair set
      P = canon(A_kl x B_kl) union canon(A_mn x B_mn), weights frozen at the peak frame.
Paths are processed in waves of W:
  1. all AFIR paths of the wave in one batched FIRE run (one model call per step);
  2. both end points of every path: batched E/F minimization; only an end whose chemical
     identity is new (not the path's reactant, not seen before) gets the perturbation
     check (a confirmed new product);
  3. paths whose bond graph changed: for each of their two highest energy maxima, path
     frames (>= 0.05 A apart in the coordinate) relaxed at their own fixed coordinate,
     growing outward from the peak one frame per side per round (all searches of a round
     in one batch) until the first + -> - pair of dE/dq (expand_brackets, 2026-11-15);
  4. each bracket: HF E/F kernel (continue_branches_ef) and E/F connection, as parallel
     tasks (parallel.run_many) with one ledger each.
Batched stages charge the wave ledger; per-bracket tasks their own ledgers. A wave writes
`wave_{w:04d}.json` when complete; rerunning skips complete waves (an interrupted wave is
redone, its paid forwards stay in its ledgers). Stops when the last --n-stop-waves waves
gave no new product (the wave form of "no new product in the last n paths"), when --budget
forwards are spent, or after --max-paths paths.
Resume only on the host that wrote the directory (NFS attribute cache 6 h).
"""

import argparse
import json
import socket
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from run_mc_afir import CEILING_RISE_EV, MAX_MAXIMA, identity, prepare  # noqa: E402  (sets sys.path)
from run_d2 import make_atoms, make_config  # noqa: E402
from systems import HF_V15, HS_BRACKET_V15, SYSTEMS  # noqa: E402

from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.afir import (  # noqa: E402
    FrozenPairSet,
    FrozenWeightedDistance,
    PairSumDistance,
    WeightedDistance,
)
from explore.afir_batch import (  # noqa: E402
    bracket_frames,
    confirm_minima,
    edges_of,
    expand_brackets,
    minimize,
    narrow_brackets,
    run_paths,
    with_coordinate,
)
from explore.afir_fragments import (  # noqa: E402
    PAIR_RATIO_MAX,
    candidate_pairs,
    perturbed_structure,
    perturbed_structure_pairs,
    ratios,
    sc_fragments,
    two_contact_combos,
)
from explore.afir_local import link  # noqa: E402
from explore.afir_window import locate_windows  # noqa: E402
from explore.afir_path import PathSettings  # noqa: E402
from explore.d2 import _ef_settings  # noqa: E402
from explore.ledger import BudgetExhausted, Ledger, tally  # noqa: E402
from explore.parallel import MACEBatchEvaluator, TaskJob, run_many  # noqa: E402
from explore.poses import _jsonable, sample_poses  # noqa: E402
from explore.ts_branch import BranchSettings, continue_branches_ef  # noqa: E402


def event_key(atoms, config, source_positions, end_positions):
    """Reaction-event identity (bonds broken / formed, canonical under the two ends'
    symmetry, direction-free): counts degenerate reactions such as b04's Cl exchange,
    whose product has the reactant's chemical identity (prrs.chemistry.reaction_event_key)."""
    from prrs.chemistry import IdentityUnavailable, reaction_event_key

    a, b = atoms.copy(), atoms.copy()
    a.positions, b.positions = np.asarray(source_positions), np.asarray(end_positions)
    try:
        return reaction_event_key(a, b, config.bond_scale, config.active_atoms,
                                  config.automorphism_limit)["key"]  # fmt: skip
    except IdentityUnavailable:
        # b05 (2026-11-14): the group generated by the two ends' automorphisms exceeds the
        # limit. Fall back to the unordered pair of chemical identities plus the element
        # pairs of the broken / formed bonds: it can merge channels, never split one.
        import hashlib

        ea = set(map(tuple, edges_of(atoms.numbers, a.positions, config.bond_scale)))
        eb = set(map(tuple, edges_of(atoms.numbers, b.positions, config.bond_scale)))
        z = atoms.numbers

        def kinds(es):
            return sorted(tuple(sorted((int(z[i]), int(z[j])))) for i, j in es)

        ends = sorted((identity(a, config), identity(b, config)))
        bonds = sorted([str(kinds(ea - eb)), str(kinds(eb - ea))])
        text = "|".join(ends + bonds)
        return "fallback:" + hashlib.sha256(text.encode()).hexdigest()[:16]


TS_CAP = 12000  # forwards per bracket task (up to 3 HF kernel attempts + connection)
# Step cap of the fixed-rho transverse relaxations of the bracket scan and narrowing (not of
# the HF kernel, which starts from converged points). The kernel's 200 (ef.max_relax_steps)
# left 91% of b05 scan points unconverged and every bracket no_bracket (PBS 28756); with
# 2000, 93/95 converged (steps p10 / median / p90 201 / 539 / 954), and the reaction-force
# sign changed after step 200 on 16/95 points, so unconverged signs are not usable
# (benchmarks/afir_mc/bracket_convergence.py, 2026-11-15; user 「可以」).
SCAN_MAX_STEPS = 2000
KERNEL_ATTEMPTS = 3  # retries only on the kernel's paired reaction-force check


def parse():
    ap = argparse.ArgumentParser()
    ap.add_argument("--system", required=True, choices=sorted(SYSTEMS))
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--shell", type=float, required=True)
    ap.add_argument("--budget", type=int, required=True)
    ap.add_argument("--wave", type=int, default=32)
    ap.add_argument("--max-paths", type=int, default=4096)
    ap.add_argument("--n-stop-waves", type=int, default=2)
    ap.add_argument("--gamma", type=float, default=106.9)
    ap.add_argument("--max-active", type=int, default=16)
    ap.add_argument("--mode", choices=("mc", "sc"), default="mc")
    ap.add_argument("--sc-poses", type=int, default=8)
    # GRRM17 discards a path whose force term has not decreased after 10 steps; with FIRE
    # starting from rest that discarded ~30 of 32 SC paths per wave on N1 (2026-11-14
    # overnight). 0 = off (default, user 2026-11-15 "可以").
    ap.add_argument("--early-discard-steps", type=int, default=0)
    # SC push coordinate (user 2026-11-15 「SC 每条路径只推被选中的原子对……冻结权重退回对照
    # 方案」): "pair" = r_kl only (V_bias = alpha r_kl); "frozen" = q with weights frozen at
    # the perturbed start (design §2.1, control); "rho" = rho_AB over the SC groups.
    ap.add_argument("--sc-coord", choices=("pair", "frozen", "rho"), default="pair")
    # TS bracket / kernel coordinate, decoupled from the push (user 2026-11-15 「推力坐标和
    # TS 定位坐标同时解耦」): "frozen-peak" = q over the SC groups A x B with weights frozen
    # at the bracket's peak frame (Claude's choice: multi-pair, no self-reinforcement,
    # weights from the contacts near the TS); "rho" = rho_AB over A x B; "push" = coupled.
    ap.add_argument("--ts-coord", choices=("frozen-peak", "rho", "push"), default="frozen-peak")
    # two-pair layer after the single-pair base (user 2026-11-15); 0 = single pairs only
    ap.add_argument("--sc-double", type=int, choices=(0, 1), default=1)
    # GRRM17's "two passed energy maxima" stop; 0 = off (user 2026-11-15), 2 = runs before
    ap.add_argument("--stop-after-maxima", type=int, default=0)
    # 2026-11-20 trial rules (docs/PLAN_AFIR_PAIR_SPLIT_2026-11-20.md §7.2 A, C, D; user
    # "可以动手试试"); all default 0 = the 11-15 rules, so earlier runs reproduce.
    # A: single-pair items on complex nodes: a pair inside one bonded fragment is pushed only
    #    on the first node (the input complex); every other node pushes inter-fragment pairs.
    # 2 (2026-11-25, after the regression lost b05 event2 found by an intra pair on node n3,
    #    an OH- assisted elimination): an intra-fragment pair is also pushed on a later node
    #    when another fragment has an atom within the candidate cut-off (r/(R+R) <= 3.0,
    #    PAIR_RATIO_MAX) of either atom of the pair; no new threshold.
    ap.add_argument("--sc-split", type=int, choices=(0, 1, 2), default=0)
    # C: no E/F minimization of a path's start when the perturbed start keeps the node's bond
    #    graph; the node MIN stands in as the start MIN.
    ap.add_argument("--skip-start-min", type=int, choices=(0, 1), default=0)
    # D: confirm_minima only for end points whose reaction event is not yet in the layer's
    #    event set, one representative per new event first, the rest only if it fails.
    ap.add_argument("--confirm-dedup", type=int, choices=(0, 1), default=0)
    # peak de-duplication before brackets (docs/PLAN_AFIR_PEAK_DEDUP_2026-11-20.md; user
    # 2026-11-22 "可以 做" for d* = 0.25 A): d* in A, 0 = off. 0.25 = about half the nearest
    # different-saddle distance of a peak frame measured on b04 / b05 / N1 (0.52 A), above
    # the largest peak -> own-saddle distance on b04 / N1 (0.22 A); peak_ts_distance.py.
    ap.add_argument("--peak-dedup", type=float, default=0.0)
    # TS locator for bond-changing paths: window = explore.afir_window (default since
    # 2026-11-26, user "可以 可以" after the window regression, DONE_LOG 11-26); hf = bracket
    # + HF kernel (the 11-15 rules, kept as the control arm; runs before 11-26 used it, so
    # reproduce or resume them with --locator hf). String iterations: --window-k.
    ap.add_argument("--locator", choices=("hf", "window"), default="window")
    ap.add_argument("--window-k", type=int, default=10)
    # 2026-11-26 (docs/PLAN_FREEZE_AND_HELDOUT_2026-11-26.md §2; user "做起来"); defaults =
    # the rules before. D1: two intra-fragment pairs (XX) only on the first complex node.
    ap.add_argument("--double-xx", choices=("all", "first"), default="all")
    # D2: two-pair items in class order II -> IX -> XX (shuffled within a class)
    ap.add_argument("--double-order", choices=("shuffle", "class"), default="shuffle")
    # B: conformer reflow, one layer. New fragment conformers found at the ends of
    # single-pair paths that kept the bond graph are placed against the other fragment's
    # input conformer in --reflow random orientations each (0 = off); their nodes push
    # inter-fragment single pairs (layer "r") and II / IX two-pair items.
    ap.add_argument("--reflow", type=int, default=0)
    return ap.parse_args()


def peak_rmsd(numbers, x, y, config):
    """Symmetric RMSD of two geometries over the union of their automorphism groups, mirror
    image included (as peak_ts_distance.py, the measurement d* comes from)."""
    from ase import Atoms

    from prrs.chemistry import automorphisms, symmetric_rmsd

    a, b = Atoms(numbers=numbers, positions=x), Atoms(numbers=numbers, positions=y)
    perms = []
    for m in (a, b):
        g, _ = automorphisms(
            m, config.bond_scale, config.active_atoms, config.automorphism_limit
        )
        perms += list(g)
    best = float("inf")
    for sign in (1.0, -1.0):
        b.positions = sign * np.asarray(y)
        best = min(best, symmetric_rmsd(a, b, perms, config.active_atoms))
    return best


def reflow_poses(atoms, frags, confs, config, out, args):
    """Conformer reflow, one layer (2026-11-26, PLAN_FREEZE_AND_HELDOUT §2.3). Reads the
    finished single-pair waves; an end MIN that converged with the node's bond graph gives
    each fragment a candidate conformer; one is new when its symmetric RMSD (fragment
    automorphisms, mirror included) to every conformer already in the library exceeds
    basin_rmsd_A. Each new conformer is placed against the other fragment's input
    conformer in --reflow orientations (poses.sample_poses, the run's shell). Returns the
    placed poses and the library record."""
    import os

    from ase import Atoms

    from explore.poses import Conformer
    from prrs.chemistry import automorphisms, symmetric_rmsd

    lib = [[(c.conformer_id, np.asarray(c.positions), None)] for c in confs]
    groups = []
    for f, c in zip(frags, confs):
        sub = Atoms(numbers=atoms.numbers[list(f.atoms)], positions=np.asarray(c.positions))
        g, _ = automorphisms(sub, config.bond_scale)
        groups.append((sub.numbers, [np.asarray(x) for x in g] or [np.arange(len(sub))]))
    new = [[] for _ in frags]
    names = sorted(
        n for n in set(os.listdir(out)) if n.startswith("wave_") and n.endswith(".json")
    )
    for name in names:
        for r in json.loads((out / name).read_text())["paths"]:
            en = (r.get("mins") or {}).get("end") or {}
            if en.get("status") != "converged" or en.get("edges") != r.get("ref_edges"):
                continue
            x = np.asarray(en["positions"], float)
            for fi, f in enumerate(frags):
                nums, g = groups[fi]
                a = Atoms(numbers=nums, positions=x[list(f.atoms)])
                best = float("inf")
                for _, pos, _ in lib[fi]:
                    for sign in (1.0, -1.0):
                        b = Atoms(numbers=nums, positions=sign * pos)
                        best = min(best, symmetric_rmsd(a, b, g))
                if best > config.basin_rmsd_A:
                    cid = f"f{f.ordinal}|{f.key}:R{len(new[fi])}"
                    lib[fi].append((cid, a.positions.copy(), r["label"]))
                    new[fi].append(
                        Conformer(cid, f.key, f.atoms, a.positions.copy(), (), "qualified")
                    )
    poses = []
    for fi, conf_list in enumerate(new):
        for c in conf_list:
            pair = (c, confs[1]) if fi == 0 else (confs[0], c)
            poses += [p for p in sample_poses(atoms, frags, pair[0], pair[1], [args.shell],
                                              args.reflow, args.seed + 101) if p.status == "placed"]  # fmt: skip
    record = [{"fragment": fi, "id": cid, "source": src} for fi, L in enumerate(lib)
              for cid, _, src in L[1:]]  # fmt: skip
    return poses, record


def components(n_atoms, edges):
    """Connected-component label per atom of a bond graph (list of [i, j])."""
    parent = list(range(n_atoms))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for i, j in edges:
        parent[find(i)] = find(j)
    return [find(a) for a in range(n_atoms)]


def spent_total(out):
    return sum(
        sum(r["forwards"] for r in tally(f).values())
        for f in sorted(set(out.glob("*.ledger.jsonl")))  # NFS readdir can repeat names
    )


def bracket_task(atoms, coord, pair, config, ef, branch, sigma, seed, ceiling, task):
    def go(factory):
        ledger = go.ledger
        tries = []
        for attempt in range(KERNEL_ATTEMPTS):
            br = continue_branches_ef(
                atoms.copy(), coord, pair, factory, config, ledger, task,
                branch_settings=branch, settings=ef, ceiling_eV=ceiling, attempt=attempt,
            )  # fmt: skip
            tries.append(br.status)
            if br.status != "unresolved:branch_force_disagreement":
                break
        t = {"status": br.status, "kernel_attempts": tries}
        if br.status != "converged_candidate":
            return t
        found = br.refine
        t |= {"positions": found["positions"], "energy_eV": found["energy_eV"],
              "metrics": found.get("metrics")}  # fmt: skip
        conn = link(atoms, br, factory, config, ledger, task + "|link", ef_settings=ef,
                    sigma=sigma, seed=seed)  # fmt: skip
        c = conn.to_dict()
        ends = {}
        for sd in c["sides"]:
            pos = (sd.get("descent") or {}).get("positions")
            if pos is not None:
                b = atoms.copy()
                b.positions = np.asarray(pos)
                ends[sd["side"]] = {"identity": identity(b, config),
                                    "energy_eV": sd["descent"].get("energy_eV"),
                                    "positions": pos}  # fmt: skip
        t["link"] = {"status": conn.status, "label": c.get("label"),
                     "fingerprint": c.get("fingerprint"), "ends": ends}  # fmt: skip
        return t

    return go


def run_wave(w, items, atoms, config, settings, ef, branch, sigma, args, batch, out,
             budget_left, known, prefix="", registry=None):  # fmt: skip
    """``items``: dicts with k, label, start (positions), coord (WeightedDistance).
    ``registry``: saddles found so far (--peak-dedup), appended in place."""
    registry = registry if registry is not None else []
    t0 = time.monotonic()
    ledger = Ledger(out / f"{prefix}wave_{w:04d}.ledger.jsonl", arm="AFIR", budget=budget_left,
                    replay=False)  # fmt: skip
    coords = [it["coord"] for it in items]
    starts = [it["start"] for it in items]
    paths = run_paths(atoms, starts, coords, batch, config, settings, ledger, f"wave|{w}|paths")
    t_paths = time.monotonic() - t0
    recs = []
    for it, p in zip(items, paths):
        fr = p.frames
        recs.append({
            "k": it["k"], "label": it["label"], "a": list(it["coord"].a),
            "b": list(it["coord"].b), "groups": it.get("groups"),
            "ts_pairs": it.get("ts_pairs"), "double": bool(it.get("double")),
            "status": p.status, "steps": len(fr),
            "maxima": p.maxima, "forwards_path": p.spent,
            "ref_edges": it.get("ref_edges") or (fr[0]["edges"] if fr else None),
            "bond_change": bool(fr)
            and fr[-1]["edges"] != (it.get("ref_edges") or fr[0]["edges"]),
            "energies_eV": [f["energy_eV"] for f in fr], "rho_A": [f["rho_A"] for f in fr],
            "edges_first": fr[0]["edges"] if fr else None,
            "edges_last": fr[-1]["edges"] if fr else None,
        })  # fmt: skip
    # 2. end-point MINs: minimize all; confirm only ends with a new identity
    rows, ends = [], []
    skip_start = bool(getattr(args, "skip_start_min", 0))
    for it, r, p in zip(items, recs, paths):
        if p.frames and p.status not in ("budget_exhausted", "domain_exit", "not_reactive"):
            if (skip_start and it.get("node_positions") is not None
                    and p.frames[0]["edges"] == it["ref_edges"]):  # fmt: skip
                # --skip-start-min: the perturbed start kept the node's bond graph; the
                # node MIN (already E/F-minimized in sc_items) is the start MIN
                nb = atoms.copy()
                nb.positions = np.asarray(it["node_positions"])
                r.setdefault("mins", {})["start"] = {
                    "status": "converged", "energy_eV": it["node_energy_eV"],
                    "edges": it["ref_edges"], "forwards": 0, "identity": identity(nb, config),
                    "positions": nb.positions, "from_node": True,
                }  # fmt: skip
                todo_ends = (("end", len(p.frames) - 1),)
            else:
                todo_ends = (("start", 0), ("end", len(p.frames) - 1))
            for name, j in todo_ends:
                rows.append((r, name))
                ends.append(p.frames[j]["positions"])
    rng = np.random.default_rng(args.seed + w)
    stopped = None
    if ends:
        try:
            relaxed = minimize(atoms, ends, batch, config, ledger, f"wave|{w}|mins",
                               fmax=config.quench_fmax_eV_A, max_steps=1000)  # fmt: skip
        except BudgetExhausted:
            relaxed, stopped = [], "budget_exhausted:mins"
        for (r, name), m in zip(rows, relaxed):
            b = atoms.copy()
            b.positions = m.positions
            r.setdefault("mins", {})[name] = {
                "status": "converged" if m.converged else "descent_not_converged",
                "energy_eV": m.energy_eV, "edges": m.edges, "forwards": m.forwards,
                "identity": identity(b, config), "positions": m.positions,
            }  # fmt: skip
        fresh = []
        for r in recs:
            mins = r.get("mins") or {}
            st, en = mins.get("start"), mins.get("end")
            if st and en and en["status"] == "converged" and en["edges"] != st["edges"]:
                fresh.append(r)
        groups = [[r] for r in fresh]
        if fresh and getattr(args, "confirm_dedup", 0):
            # --confirm-dedup: the event key decides what a confirmation can add (account()
            # only counts keys not yet in the layer's set), so known events are not
            # re-confirmed and a new event is confirmed once, the others only as fallback
            by_key = {}
            for r in fresh:
                try:
                    key = event_key(atoms, config, r["mins"]["start"]["positions"],
                                    r["mins"]["end"]["positions"])  # fmt: skip
                except Exception:  # noqa: BLE001  (unkeyable: confirm on its own)
                    key = ("unkeyed", r["k"])
                r["mins"]["end"]["event_key_precheck"] = str(key)
                if key in known:
                    r["mins"]["end"]["confirmation"] = {"status": "skipped_known_event"}
                    continue
                by_key.setdefault(key, []).append(r)
            groups = list(by_key.values())
        if groups and stopped is None:
            try:
                rounds = [[g[0] for g in groups]]
                while rounds[-1]:
                    batch_rows = rounds[-1]
                    conf = confirm_minima(atoms, [r["mins"]["end"]["positions"] for r in batch_rows],
                                          batch, config, ledger, f"wave|{w}|confirm", rng=rng)  # fmt: skip
                    for r, c in zip(batch_rows, conf):
                        r["mins"]["end"]["confirmation"] = {
                            k: c.get(k) for k in ("status", "trials", "forwards")
                        }
                    # next representative of every group whose current one failed
                    nxt = []
                    for g in groups:
                        if any((x["mins"]["end"].get("confirmation") or {}).get("status")
                               == "confirmed" for x in g):  # fmt: skip
                            continue
                        rest = [x for x in g if "confirmation" not in x["mins"]["end"]]
                        if rest:
                            nxt.append(rest[0])
                    rounds.append(nxt)
                for g in groups:  # members left unconfirmed because a sibling confirmed
                    for x in g:
                        x["mins"]["end"].setdefault(
                            "confirmation", {"status": "skipped_same_event_in_wave"}
                        )
            except BudgetExhausted:
                stopped = "budget_exhausted:confirm"
    t_mins = time.monotonic() - t0 - t_paths
    if getattr(args, "locator", "hf") == "window":
        # --locator window (2026-11-25, user "可以"): bond-change windows on a short String,
        # fixed-direction climb, E/F descent on both sides (explore.afir_window). The HF
        # bracket / kernel path below is not run; --peak-dedup does not apply.
        t1 = time.monotonic()
        sel = []
        for r, p in zip(recs, paths) if stopped is None else ():  # noqa: B007
            if not r["bond_change"]:
                continue
            mins = r.get("mins") or {}
            st, en = mins.get("start") or {}, mins.get("end") or {}
            if st.get("status") != "converged" or en.get("status") != "converged":
                r.setdefault("ts", []).append({"locator": "window", "status": "window_no_ends"})
                continue
            sel.append((r, p, st, en))
        found = []
        if sel:
            try:
                found = locate_windows(
                    atoms,
                    [{"frames": np.stack([f["positions"] for f in p.frames]),
                      "frame_energies": [f["energy_eV"] for f in p.frames],
                      "start": {"positions": np.asarray(st["positions"]), "energy_eV": st["energy_eV"]},
                      "end": {"positions": np.asarray(en["positions"]), "energy_eV": en["energy_eV"]}}
                     for r, p, st, en in sel],
                    batch, config, ledger, f"wave|{w}|window", k=args.window_k,
                    ceiling_rise=CEILING_RISE_EV,
                )  # fmt: skip
            except BudgetExhausted:
                stopped, found = "budget_exhausted:window", []
        n_conv = 0
        for (r, p, st, en), res in zip(sel, found):
            for wr in res:
                ci = wr["ci"]
                ok = ci["ci_status"] == "converged"
                t = {"locator": "window", "window": wr["window"], "jump": wr["jump"],
                     "status": "converged_candidate" if ok else f"window:{ci['ci_status']}",
                     "ci_steps": ci["ci_steps"], "energy_eV": ci["ci_energy_eV"],
                     "ci_fmax_eV_A": ci["ci_fmax_eV_A"], "forwards": ci["ci_forwards"]}  # fmt: skip
                if ok:
                    n_conv += 1
                    t["positions"] = ci["ci_positions"]
                    sd = wr["sides"]
                    ends = {}
                    if all(m.converged for m in sd):
                        lo_m, hi_m = sorted(sd, key=lambda m: m.energy_eV)
                        for side, m in (("low", lo_m), ("high", hi_m)):
                            b = atoms.copy()
                            b.positions = m.positions
                            ends[side] = {"identity": identity(b, config),
                                          "energy_eV": m.energy_eV, "positions": m.positions}  # fmt: skip
                    status = ("side_unconverged" if len(ends) < 2
                              else "connected" if sd[0].edges != sd[1].edges
                              else "same_basin")  # fmt: skip
                    t["link"] = {"status": status, "ends": ends,
                                 "sides_match_path": [m.edges == st["edges"] or m.edges == en["edges"]
                                                      for m in sd]}  # fmt: skip
                r.setdefault("ts", []).append(t)
        ledger.close()
        return recs, {"paths_s": t_paths, "mins_s": t_mins, "bracket_s": 0.0,
                      "ts_s": time.monotonic() - t1, "wave_ledger": ledger.spent,
                      "brackets": 0, "ts_tasks": 0, "windows": sum(len(x) for x in found),
                      "windows_converged": n_conv, "peaks_skipped": 0, "stopped": stopped}  # fmt: skip
    # 3-4. candidate peaks of bond-changing paths -> brackets -> HF kernel + connection
    cands = []
    for r, p, c in zip(recs, paths, coords) if stopped is None else ():  # noqa: B007
        if not r["bond_change"]:
            continue
        # the highest maxima, not the first: without the two-maxima stop a path also has
        # small FIRE overshoot maxima (2026-11-15)
        tops = sorted(p.maxima, key=lambda k: -p.frames[k]["energy_eV"])[:MAX_MAXIMA]
        for j, top in enumerate(tops):
            cands.append((r, p, c, j, top))
    d_star = float(getattr(args, "peak_dedup", 0.0) or 0.0)
    acc = {"items": 0, "jobs": 0, "bracket_s": 0.0, "ts_s": 0.0, "skipped": 0}

    def stage(sel):
        """Steps 3-4 for the candidate peaks ``sel``; returns the stop reason or None."""
        nonlocal stopped
        t1 = time.monotonic()
        items = []
        for r, p, c, j, top in sel:
            ts = ts_coordinate(atoms.numbers, c, r, p.frames[top], args)
            frames = p.frames if ts is c else with_coordinate(p.frames, ts)
            back, fwd = bracket_frames(frames, top)
            items.append((r, p, ts, j, top, back, fwd, frames))
        expanded = []
        if items:
            try:
                expanded = expand_brackets(
                    atoms, [(it[7], it[5], it[6], it[2]) for it in items], batch, config,
                    ledger, f"wave|{w}|bracket", fmax=ef.transverse_fmax_eV_A,
                    max_steps=SCAN_MAX_STEPS,
                )  # fmt: skip
            except BudgetExhausted:
                stopped, items = "budget_exhausted:bracket", []
        picked = []
        for it, (pb, pf, pair) in zip(items, expanded):
            picked.append((it, pb + pf, pair))
        todo = [x for x in picked if x[2] is not None]
        narrowed = {}
        if todo and stopped is None:
            try:
                new_pairs, hist = narrow_brackets(
                    atoms, [x[0][2] for x in todo], [x[2] for x in todo], batch, config,
                    ledger, f"wave|{w}|narrow", fmax=ef.transverse_fmax_eV_A,
                    max_steps=SCAN_MAX_STEPS,
                )  # fmt: skip
                for x, pr, h in zip(todo, new_pairs, hist):
                    narrowed[id(x[0])] = (pr, h)
            except BudgetExhausted:
                stopped = "budget_exhausted:narrow"
        t2 = time.monotonic()
        jobs, meta, done_pairs = [], [], set()
        for it, p, pair in picked:
            r, path, coord, j, top, back, fwd, _ = it
            raw = pair
            if id(it) in narrowed:
                pair, hist = narrowed[id(it)]
            else:
                hist = None
            key = (
                None
                if raw is None
                else (r["k"], round(raw[0].q_target, 6), round(raw[1].q_target, 6))
            )
            if key is not None and key in done_pairs:
                r.setdefault("ts", []).append({"top": top, "status": "same_bracket"})
                continue
            if key is not None:
                done_pairs.add(key)
            t = {"top": top, "bracket_raw": None if raw is None else [q.q_target for q in raw],
                 "bracket": None if pair is None else [q.q_target for q in pair],
                 "narrow": hist,
                 "scan": [{"q": q.q_target, "E": q.energy_eV, "rf": q.reaction_force,
                           "converged": q.converged} for q in p]}  # fmt: skip
            r.setdefault("ts", []).append(t)
            if pair is None:
                t["status"] = "no_bracket"
                continue
            task = f"ts|{r['k']}|{j}"
            led = Ledger(out / f"{prefix}ts_{r['k']:05d}_{j}.ledger.jsonl", arm="AFIR",
                         budget=TS_CAP, replay=False)  # fmt: skip
            go = bracket_task(atoms, coord, pair, config, ef, branch, sigma,
                              args.seed + 1000 * r["k"] + j, path.frames[0]["energy_eV"]
                              + CEILING_RISE_EV, task)  # fmt: skip
            go.ledger = led
            jobs.append(TaskJob(go, led))
            meta.append((t, led, r, path, top))
        if jobs:
            outcomes, _ = run_many(jobs, batch, args.max_active)
            for (t, led, r, path, top), o in zip(meta, outcomes):
                led.close()
                t |= o.result if o.error is None else {"status": f"error:{o.error}"}
                t["forwards"] = led.spent
                if d_star > 0 and t.get("status") == "converged_candidate" \
                        and (t.get("link") or {}).get("status") == "connected":  # fmt: skip
                    entry = {"event": r.get("peak_event"), "label": r["label"], "top": top,
                             "wave": f"{prefix}{w}", "energy_eV": t.get("energy_eV"),
                             "positions": np.asarray(t["positions"]).tolist()}  # fmt: skip
                    registry.append(entry)
                    with open(out / "ts_registry.jsonl", "a", encoding="utf-8") as fh:
                        fh.write(json.dumps(entry) + "\n")
        acc["items"] += len(items)
        acc["jobs"] += len(jobs)
        acc["bracket_s"] += t2 - t1
        acc["ts_s"] += time.monotonic() - t2

    if d_star > 0 and cands:
        # --peak-dedup (docs/PLAN_AFIR_PEAK_DEDUP_2026-11-20.md): a peak within d* (symmetric
        # RMSD, mirror included) of a registered saddle of the same path event is not
        # bracketed; one representative per event first, the rest after it is registered
        for r, p, c, j, top in cands:
            if "peak_event" in r:
                continue
            mins = r.get("mins") or {}
            st, en = mins.get("start") or {}, mins.get("end") or {}
            key = ("path", r["k"])
            if st.get("positions") is not None and en.get("positions") is not None:
                try:
                    key = event_key(atoms, config, st["positions"], en["positions"])
                except Exception:  # noqa: BLE001
                    pass
            r["peak_event"] = key if isinstance(key, str) else None
            r["_peak_group"] = key

        def near_known(r, p, top):
            if r["peak_event"] is None:
                return None
            peak = np.asarray(p.frames[top]["positions"])
            best = None
            for e in registry:
                if e["event"] != r["peak_event"]:
                    continue
                d = peak_rmsd(atoms.numbers, peak, np.asarray(e["positions"]), config)
                if best is None or d < best[0]:
                    best = (d, e["label"], e["wave"])
            return best if best is not None and best[0] <= d_star else None

        def split_known(sel):
            keep = []
            for r, p, c, j, top in sel:
                hit = near_known(r, p, top)
                if hit is None:
                    keep.append((r, p, c, j, top))
                else:
                    r.setdefault("ts", []).append(
                        {"top": top, "status": "skipped_near_known_ts",
                         "near_ts": {"rmsd_A": hit[0], "label": hit[1], "wave": hit[2]}}
                    )  # fmt: skip
                    acc["skipped"] += 1
            return keep

        rest = split_known(cands)
        first, seen = [], set()
        for cand in sorted(rest, key=lambda x: -x[1].frames[x[4]]["energy_eV"]):
            g = cand[0]["_peak_group"]
            if g not in seen:
                seen.add(g)
                first.append(cand)
        stage(first)
        if stopped is None:
            chosen = {(id(x[0]), x[4]) for x in first}
            stage(split_known([x for x in rest if (id(x[0]), x[4]) not in chosen]))
        for r in recs:
            r.pop("_peak_group", None)
    elif cands:
        stage(cands)
    ledger.close()
    return recs, {"paths_s": t_paths, "mins_s": t_mins, "bracket_s": acc["bracket_s"],
                  "ts_s": acc["ts_s"], "wave_ledger": ledger.spent, "brackets": acc["items"],
                  "ts_tasks": acc["jobs"], "peaks_skipped": acc["skipped"],
                  "stopped": stopped}  # fmt: skip


def push_coordinate(numbers, i, j, a_grp, b_grp, x0, kind):
    """SC push coordinate (--sc-coord): pair = r_ij alone (WeightedDistance on A = {i},
    B = {j} is exactly r_ij, so V_bias = alpha r_ij); frozen = weights frozen at x0; rho."""
    if kind == "pair":
        return WeightedDistance.build(numbers, (i,), (j,), intent_id="afir-pair")
    if kind == "frozen":
        return FrozenWeightedDistance.build(numbers, a_grp, b_grp, x0)
    return WeightedDistance.build(numbers, a_grp, b_grp)


def ts_coordinate(numbers, push, rec, peak_frame, args):
    """TS bracket / kernel coordinate (--ts-coord), decoupled from the push. frozen-peak:
    FrozenPairSet over the record's pair set P (single pair: canon(A x B); two pairs: the
    union of both canon(A x B)), weights at the peak frame. rho: rho_AB over A x B, single
    pairs only (two-pair paths use frozen-peak). MC paths (no pair set) keep the push."""
    kind = getattr(args, "ts_coord", "push")
    pairs = rec.get("ts_pairs")
    if kind == "push" or not pairs:
        return push
    if kind == "rho" and not rec.get("double"):
        a_grp, b_grp = rec["groups"]
        return WeightedDistance.build(numbers, a_grp, b_grp)
    return FrozenPairSet.build(numbers, pairs, peak_frame["positions"],
                               intent_id="afir-ts-frozen-peak")  # fmt: skip


def _canon_product(a_grp, b_grp):
    return [list(p) for p in FrozenPairSet.canon((i, j) for i in a_grp for j in b_grp)]


def sc_items(atoms, poses, config, batch, out, args, doubles=None, double_filter=None, *,
             node_prefix="n", include_input=True, intra="rule",
             ledger_name="sc_nodes.ledger.jsonl", k_offset=0, n_poses=None):  # fmt: skip
    """SC-AFIR paths on the reactant-complex MINs (module doc, mode sc). Two-pair items
    are appended to ``doubles`` (a list) when given and --sc-double is on;
    ``double_filter(label) -> bool`` (diagnostic replays only) skips building the others.
    Keys k and the shuffled order are then those of the filtered pool, not of a run.
    Conformer reflow (2026-11-26) calls it with node_prefix "r", include_input False and
    intra "none" (only inter-fragment single pairs, no XX two-pair combinations)."""
    from ase import Atoms

    from prrs.chemistry import automorphisms, symmetric_rmsd

    n_poses = args.sc_poses if n_poses is None else n_poses
    starts = ([atoms.positions] if include_input else []) + [
        p.positions for p in poses[:n_poses]
    ]
    if not starts:
        return [], {"nodes": [], "forwards": 0, "paths": 0, "double_paths": 0}
    led = Ledger(out / ledger_name, arm="AFIR", replay=False)
    relaxed = minimize(atoms, starts, batch, config, led, "sc|nodes",
                       fmax=config.quench_fmax_eV_A, max_steps=1000)  # fmt: skip
    led.close()
    nodes = []
    for j, m in enumerate(relaxed):
        if not m.converged:
            continue
        a = Atoms(numbers=atoms.numbers, positions=m.positions)
        group, _ = automorphisms(a, config.bond_scale, config.active_atoms,
                                 config.automorphism_limit)  # fmt: skip
        dup = any(
            abs(m.energy_eV - n["energy_eV"]) <= config.basin_energy_eV
            and m.edges == n["edges"]
            and symmetric_rmsd(a, n["atoms"], group, config.active_atoms) <= config.basin_rmsd_A
            for n in nodes
        )
        if not dup:
            nodes.append({"source": j, "energy_eV": m.energy_eV, "edges": m.edges, "atoms": a})
    items = []
    split = int(getattr(args, "sc_split", 0))
    for ni, n in enumerate(nodes):
        x = n["atoms"].positions
        comp = components(len(atoms), n["edges"])
        rr = ratios(atoms.numbers, x)[0] if split == 2 else None
        for i, j in candidate_pairs(n["atoms"], config):
            if intra == "none" and comp[i] == comp[j]:
                continue
            if split and ni > 0 and comp[i] == comp[j]:
                # intra-fragment pair: first node only (--sc-split 1); with 2 also where a
                # partner fragment is in contact with i or j
                if split == 1 or not any(
                    comp[m] != comp[i] and min(rr[i, m], rr[j, m]) <= PAIR_RATIO_MAX
                    for m in range(len(atoms))
                ):
                    continue
            px = perturbed_structure(
                atoms.numbers, x, i, j, sign=+1, seed=args.seed + i * 97 + j
            )
            try:
                a_grp, b_grp = sc_fragments(atoms.numbers, px, i, j)
            except ValueError:
                continue
            items.append({"k": k_offset + len(items), "label": f"{node_prefix}{ni}|{i}-{j}",
                          "start": px,
                          "coord": push_coordinate(atoms.numbers, i, j, a_grp, b_grp, px,
                                                   args.sc_coord),
                          "groups": [list(a_grp), list(b_grp)],
                          "ts_pairs": _canon_product(a_grp, b_grp),
                          "ref_edges": n["edges"], "node_positions": x,
                          "node_energy_eV": n["energy_eV"]})  # fmt: skip
    info = [{"node": ni, "source": n["source"], "energy_eV": n["energy_eV"]} for ni, n in
            enumerate(nodes)]  # fmt: skip
    n_double = 0
    if doubles is not None and getattr(args, "sc_double", 0):
        pool, rank = [], []
        xx_first = getattr(args, "double_xx", "all") == "first"
        for ni, n in enumerate(nodes):
            x = n["atoms"].positions
            comp = components(len(atoms), n["edges"])
            multi = len(set(comp)) > 1
            for (k, l), (m, q) in two_contact_combos(n["atoms"], config):  # noqa: E741
                if double_filter is not None and not double_filter(
                    f"{node_prefix}{ni}|{k}-{l}+{m}-{q}"
                ):
                    continue
                n_inter = (comp[k] != comp[l]) + (comp[m] != comp[q])
                # --double-xx first (2026-11-26, D1): on a complex, two intra-fragment pairs
                # only on the first node; reflow nodes (intra "none") never
                if n_inter == 0 and multi and (intra == "none" or (xx_first and ni > 0)):
                    continue
                px = perturbed_structure_pairs(
                    atoms.numbers, x, [(k, l), (m, q)], sign=+1,
                    seed=args.seed + 7919 * (k * 97 + l) + m * 97 + q,
                )  # fmt: skip
                try:
                    g1 = sc_fragments(atoms.numbers, px, k, l)
                    g2 = sc_fragments(atoms.numbers, px, m, q)
                except ValueError:
                    continue
                ts_pairs = FrozenPairSet.canon(
                    [(i, j) for i in g1[0] for j in g1[1]]
                    + [(i, j) for i in g2[0] for j in g2[1]]
                )
                rank.append(2 - n_inter)  # II 0, IX 1, XX 2
                pool.append({"label": f"{node_prefix}{ni}|{k}-{l}+{m}-{q}", "start": px,
                             "coord": PairSumDistance.build([(k, l), (m, q)]),
                             "groups": [[list(g1[0]), list(g1[1])], [list(g2[0]), list(g2[1])]],
                             "ts_pairs": [list(p) for p in ts_pairs], "double": True,
                             "ref_edges": n["edges"], "node_positions": x,
                             "node_energy_eV": n["energy_eV"]})  # fmt: skip
        order = np.random.default_rng(args.seed).permutation(len(pool))
        if getattr(args, "double_order", "shuffle") == "class":
            # --double-order class (2026-11-26, D2): II, then IX, then XX; shuffled within
            order = sorted(order, key=lambda idx: rank[idx])
        for pos, idx in enumerate(order):
            doubles.append(
                pool[idx] | {"k": k_offset + len(items) + pos, "pair_class": rank[idx]}
            )
        n_double = len(pool)
    return items, {"nodes": info, "forwards": led.spent, "paths": len(items),
                   "double_paths": n_double}  # fmt: skip


def main():
    args = parse()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    system = SYSTEMS[args.system]
    config = make_config(system)
    t0 = time.monotonic()
    calc = polar1("cuda", "float32")
    load_s = time.monotonic() - t0
    atoms, _ = make_atoms(system, config)
    probe = atoms.copy()
    probe.calc = calc
    e_in = float(probe.get_potential_energy())
    noise = 2 * float(np.spacing(np.float32(abs(e_in))))
    frags, confs, prep = prepare(out, system, config, atoms, calc)
    base = PathSettings(noise_eV=noise, gamma_kJ=args.gamma)
    settings = PathSettings(
        noise_eV=noise,
        gamma_kJ=args.gamma,
        early_discard_steps=args.early_discard_steps or base.max_steps + 1,
        stop_after_maxima=args.stop_after_maxima,
    )
    ef = _ef_settings(HF_V15, HS_BRACKET_V15)
    branch = BranchSettings(geometry_noise_A=HF_V15.geometry_noise_A)
    sigma = HF_V15.reaction_force_sigma
    if len(frags) == 1:
        # unimolecular input (2026-11-26, PLAN_FREEZE_AND_HELDOUT §3.0): no orientations, the
        # input structure is the only SC node; every candidate pair is intramolecular
        if args.mode == "mc":
            raise ValueError("--mode mc needs two fragments; use --mode sc")
        poses = []
    else:
        poses = [p for p in sample_poses(atoms, frags, confs[0], confs[1], [args.shell],
                                         args.max_paths, args.seed) if p.status == "placed"]  # fmt: skip
    batch = MACEBatchEvaluator(calc)
    if args.mode == "mc":
        items = [
            {"k": k, "label": p.pose_id, "start": p.positions,
             "coord": WeightedDistance.build(atoms.numbers, frags[0].atoms, frags[1].atoms)}
            for k, p in enumerate(poses)
        ]  # fmt: skip
        nodes_info, doubles = None, []
    else:
        doubles = []
        items, nodes_info = sc_items(atoms, poses, config, batch, out, args, doubles)
    manifest = {"runner": "run_mc_afir_batch.py", "system": args.system, "args": vars(args),
                "noise_eV": noise, "input_energy_eV": e_in, "path_settings": settings.__dict__,
                "ef_settings": ef.__dict__, "branch_settings": branch.__dict__,
                "ts_cap": TS_CAP, "scan_max_steps": SCAN_MAX_STEPS, "host": socket.gethostname(), "load_s": load_s,
                "placed_poses": len(poses), "paths_planned": len(items),
                "double_paths_planned": len(doubles),
                "sc_nodes": nodes_info}  # fmt: skip
    (out / "manifest.json").write_text(json.dumps(_jsonable(manifest), indent=1, default=str))

    def account(recs, products):
        """New reaction events: a connected TS's two ends, or a path whose confirmed end
        differs in bonds from its start MIN; keyed by reaction_event_key (not the
        product's chemical identity: degenerate reactions keep the identity). ``products``
        is the layer's own event set (single / two-pair layers are accounted apart)."""
        new = []
        for r in recs:
            keys = []
            mins = r.get("mins") or {}
            st, en = mins.get("start") or {}, mins.get("end") or {}
            if (en.get("confirmation") or {}).get("status") == "confirmed" and st.get("positions") \
                    is not None and st.get("edges") != en.get("edges"):  # fmt: skip
                keys.append(event_key(atoms, config, st["positions"], en["positions"]))
            for t in r.get("ts", []):
                link = t.get("link") or {}
                ends = link.get("ends") or {}
                if link.get("status") == "connected" and {"low", "high"} <= set(ends):
                    keys.append(event_key(atoms, config, ends["low"]["positions"],
                                          ends["high"]["positions"]))  # fmt: skip
            fresh = list(dict.fromkeys(x for x in keys if x not in products))
            products.update(fresh)
            r["new_events"] = fresh
            new += fresh
        return new

    def run_layer(layer_items, prefix, products):
        quiet, w = 0, 0
        while w * args.wave < len(layer_items):
            dest = out / f"{prefix}wave_{w:04d}.json"
            if dest.exists():
                rec = json.loads(dest.read_text())
                new = account(rec["paths"], products)
            else:
                left = args.budget - spent_total(out)
                if left <= 0:
                    print(f"stop: budget ({prefix or 'single'} layer)", flush=True)
                    break
                chunk = layer_items[w * args.wave : (w + 1) * args.wave]
                recs, timing = run_wave(w, chunk, atoms, config, settings, ef, branch, sigma,
                                        args, batch, out, left, set(products), prefix,
                                        registry)  # fmt: skip
                # write the wave first: an accounting failure must not lose paid results
                rec = {"wave": w, "paths": recs, "timing": timing,
                       "spent_total": spent_total(out)}  # fmt: skip
                dest.write_text(json.dumps(_jsonable(rec), default=str) + "\n")
                try:
                    new = account(recs, products)
                except Exception as exc:  # recorded, the run goes on
                    new, rec["accounting_error"] = [], f"{type(exc).__name__}: {exc}"
                rec["new_products"] = new
                dest.write_text(json.dumps(_jsonable(rec), default=str) + "\n")
                st = {}
                for r in recs:
                    st[r["status"]] = st.get(r["status"], 0) + 1
                ts = [t.get("status") for r in recs for t in r.get("ts", [])]
                print(prefix or "s", w, st, "bond_change", sum(r["bond_change"] for r in recs),
                      "ts", ts, "new", len(new),
                      {k: round(v, 1) if isinstance(v, float) else v for k, v in timing.items()},
                      rec["spent_total"], flush=True)  # fmt: skip
            quiet = 0 if new else quiet + 1
            w += 1
            if args.mode == "mc" and quiet >= args.n_stop_waves:
                print(f"stop: no new product in the last {args.n_stop_waves} waves", flush=True)
                break
        return w

    def layer_forwards(prefix):
        pats = (f"{prefix}wave_*.ledger.jsonl", f"{prefix}ts_*.ledger.jsonl")
        return sum(sum(r["forwards"] for r in tally(f).values()) for pat in pats
                   for f in sorted(set(out.glob(pat))))  # fmt: skip

    # saddles registered by --peak-dedup, shared by both layers; read back on resume
    registry = []
    reg_file = out / "ts_registry.jsonl"
    if reg_file.exists():
        registry = [json.loads(x) for x in reg_file.read_text().splitlines() if x.strip()]
    single, double = set(), set()
    waves = run_layer(items, "", single)
    r_waves = 0
    if args.mode == "sc" and args.reflow > 0 and len(frags) == 2:
        r_poses, lib = reflow_poses(atoms, frags, confs, config, out, args)
        r_doubles = [] if doubles is not None else None
        r_items, r_info = sc_items(
            atoms, r_poses, config, batch, out, args, r_doubles, node_prefix="r",
            include_input=False, intra="none", ledger_name="reflow_nodes.ledger.jsonl",
            k_offset=len(items) + len(doubles), n_poses=len(r_poses),
        )  # fmt: skip
        (out / "reflow.json").write_text(json.dumps(_jsonable(
            {"conformers": lib, "poses": len(r_poses), "nodes": r_info}), default=str))  # fmt: skip
        r_waves = run_layer(r_items, "r", single)
        if r_doubles:
            doubles = doubles + r_doubles
            if args.double_order == "class":
                doubles = sorted(doubles, key=lambda it: it.get("pair_class", 1))
    d_waves = run_layer(doubles, "d", double) if doubles else 0
    summary = {"waves": waves, "reflow_waves": r_waves, "events": sorted(single),
               "double_waves": d_waves,
               "double_events": sorted(double), "double_new_events": sorted(double - single),
               "forwards_single_layer": layer_forwards(""), "forwards_reflow_layer": layer_forwards("r"),
               "forwards_double_layer": layer_forwards("d"),
               "forwards_total": spent_total(out), "wall_s": time.monotonic() - t0}  # fmt: skip
    (out / "summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
