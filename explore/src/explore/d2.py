"""D2 orchestration, discovery v1 (DISCOVERY_LAYER_DESIGN_2026-10-13 §1.1, §4.2-4.3,
§5, §6.1-6.3 v1 rows, §9; EVALUATION_PROTOCOL_DESIGN_2026-10-13 §2.2-2.3, §4.1).

A run directory holds manifest.json (frozen settings, config, identity), poses.jsonl,
population.jsonl (the intent table: order keys, audit draws, scan domains),
journal.jsonl (dispatch / settle / reconcile / skip in coordinator order) and
ledger.jsonl (every forward, one Ledger per task attempt). ``run`` is a function of
these files and the two backends: rerunning it after a crash skips settled tasks,
charges interrupted attempts (their cap when a segment was left open) and redoes them
under a new attempt number. Nothing is reused by hash.

Resume a run only on the host that wrote it. On NFS with attribute caching (this
project: 6 h, nocto) another host can see stale journal/ledger files, reuse attempt
numbers, and ``repair_tail`` could then truncate the real file. Every journal event
and the manifest record the writing host, and the report lists the hosts seen, so a
cross-host resume is visible."""

from __future__ import annotations

import json
import socket
import traceback
from collections import Counter
from dataclasses import asdict, dataclass, replace
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator

from .certify import Certification, certify, pair_symmetries
from .journal import Account, Journal, repair_tail
from .ledger import Ledger, tally, view_factory
from .parallel import ScanJob, TaskJob, run_many, scan_many
from .poses import _jsonable
from .probe import (
    Intent,
    ScanPoint,
    ScanSettings,
    bracket_jump,
    feasible_prefix,
    read_scan,
    release,
    reverse_plan,
)
from .streams import stream
from .ts_branch import BranchSettings, continue_branches_ef
from .ts_connect import ConnectSettings, branch_starts, trace_sides_ef
from .ts_ef import EFSettings

ARM = "D2"
RULE_VERSION = "D2-v1.1"  # v1.1: F3 read_scan fix + F6 scan tangent to the refiner (2026-10-24)
# HF-4 (docs/PLAN_HF4_D2_2026-10-26.md): L4 replaced by E/F search (HS) + E/F connection
# (HC); L0-L3, promotion and strict phases unchanged. Chosen at freeze by passing ``hf``.
HF_RULE_VERSION = "D2-HF-v1"
# D2-HF-v1 + start-cluster dispatch order (docs/PLAN_START_CLUSTER_2026-10-28.md): HS of
# promoted items whose starting guesses coincide is dispatched in two rounds. Chosen at
# freeze by passing ``start_cluster`` together with ``hf``.
HF_SC_RULE_VERSION = "D2-HF-v1.1"
# D2-HF-v1.1 + same-graph ridges last (docs/PLAN_HF_SAME_GRAPH_RIDGE_2026-11-07.md; user
# 2026-11-07, R3a-2 of docs/HF_HS_NO_MAXIMUM_OPTIONS_2026-11-07.md): a ridge whose scan
# peak and next point share their bond graph gets its HS only after both start-cluster
# rounds, on what they left; never dropped, skipped as ``HS_deferred:same_graph_ridge``.
# On eight b05 HF runs such ridges took 92% of ridge HS forwards and gave no candidate.
HF_SG_RULE_VERSION = "D2-HF-v1.2"
# D2-HF-v1.2 + HS bracket rule HF-EF-bracket-v7 (TODO item 9d; user 2026-11-07 「准备」):
# two-stage midpoints with a FIRE fine stage (v5, user 2026-11-04), coarse-stage candidate
# (v6, user 2026-11-06) and endpoint candidate (v7, user 2026-11-07). Chosen at freeze by
# passing ``hs_bracket`` together with ``same_graph_ridge_last``.
HF_BR_RULE_VERSION = "D2-HF-v1.3"
# D2-HF-v1.3 + same-graph ridge sample (S3) + HS / HC in parallel waves
# (docs/PLAN_HF_D2_V14_2026-11-07.md; user 2026-11-07). Same-graph ridges get their HS
# only when a seeded draw falls below ``audit_fraction``; the rest are skipped as
# ``HS_skipped:same_graph_ridge``. HS and HC of each dispatch group run in waves of
# ``wave_size`` through one resident batch evaluator (``max_active`` concurrent tasks)
# and settle in input order; HC duplicate checks happen at settlement, in order.
# Chosen at freeze by passing ``same_graph_sample`` and ``hf_parallel`` on top of v1.3.
HF_V14_RULE_VERSION = "D2-HF-v1.4"


@dataclass(frozen=True)
class D2Settings:
    seed: int
    search_budget: int
    certify_budget: int
    scan_cap: int
    reverse_cap: int
    release_cap: int
    certify_cap: int
    per_pose_quota: int
    audit_fraction: float
    noise_eV: float
    certify_max_rise_eV: float
    magnitudes_A: tuple[float, ...]
    fmax_eV_A: float
    max_steps_per_point: int
    max_rise_eV: float
    max_active: int
    wave_size: int

    def __post_init__(self):
        for name in (
            "search_budget",
            "certify_budget",
            "scan_cap",
            "reverse_cap",
            "release_cap",
            "certify_cap",
            "per_pose_quota",
            "max_active",
            "wave_size",
        ):
            if getattr(self, name) < 1:
                raise ValueError(f"{name} must be positive")
        if not 0.0 <= self.audit_fraction <= 1.0:
            raise ValueError("audit_fraction must be in [0, 1]")
        mags = list(self.magnitudes_A)
        if not mags or any(m <= 0 for m in mags) or mags != sorted(mags):
            raise ValueError("magnitudes_A must be positive and ascending")


@dataclass(frozen=True)
class HFSettings:
    """Per-system measured noise inputs and per-task caps of the HF L4 (no defaults on
    purpose: a missing value must not silently borrow another system's measurement).
    sigma / geometry noise: same-branch reaction-force jitter and max pairwise RMSD of
    independent re-relaxations (benchmarks/explore_d2/hf_ef_rf_jitter.py)."""

    reaction_force_sigma: float
    geometry_noise_A: float
    midpoint_transverse_fmax_eV_A: float
    hs_cap: int
    hc_cap: int

    def __post_init__(self):
        for name, value in asdict(self).items():
            if isinstance(value, bool) or not np.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be finite and positive")
        for name in ("hs_cap", "hc_cap"):
            if not isinstance(getattr(self, name), (int, np.integer)):
                raise ValueError(f"{name} must be an integer")


@dataclass(frozen=True)
class StartCluster:
    """Two-round HS dispatch of D2-HF-v1.1 (docs/PLAN_START_CLUSTER_2026-10-28.md §3).
    Promoted items join the first earlier cluster whose leader's starting guess is within
    ``noise_multiple`` x HFSettings.geometry_noise_A (symmetric RMSD under the pose-graph
    automorphisms); the first ``first_round`` members of every cluster get their HS in
    round one, the rest only with the refine budget left after it. Nothing is dropped:
    a round-two item the budget cannot reach is skipped as ``HS_deferred:same_start``."""

    noise_multiple: float
    first_round: int

    def __post_init__(self):
        m = self.noise_multiple
        if isinstance(m, bool) or not np.isfinite(m) or m <= 0:
            raise ValueError("noise_multiple must be finite and positive")
        k = self.first_round
        if isinstance(k, bool) or not isinstance(k, (int, np.integer)) or k < 1:
            raise ValueError("first_round must be a positive integer")


@dataclass(frozen=True)
class HSBracket:
    """Per-system inputs of the D2-HF-v1.3 HS bracket rule (HF-EF-bracket-v7); no defaults
    on purpose. Coarse midpoint threshold and the sign bound measured for that threshold
    (b05: 2e-3 and 0.025 eV/A, docs/PLAN_HF_MIDPOINT_COARSE_2026-11-02.md §2). The FIRE fine
    stage, the coarse-stage candidate and the endpoint candidate are fixed by the rule."""

    midpoint_coarse_transverse_fmax_eV_A: float
    midpoint_sign_bound_eV_A: float

    def __post_init__(self):
        for name, value in asdict(self).items():
            if isinstance(value, bool) or not np.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be finite and positive")


@dataclass
class Backends:
    """``search_batch``: ``evaluate_batch`` for ``scan_many`` (resident float32 model in
    production). ``certify_calculator``: one resident physical calculator for release
    and L4 (float64 in production)."""

    search_batch: object
    certify_calculator: object


def _key(seed, *parts):
    return float(stream(seed, *parts).random())


def init_run(
    run_dir,
    poses,
    specs_by_pose,
    settings,
    config,
    identity,
    hf=None,
    start_cluster=None,
    same_graph_ridge_last=False,
    hs_bracket=None,
    same_graph_sample=False,
    hf_parallel=False,
):
    """Freeze a D2 run: write poses, the intent table and the manifest (last, so its
    presence marks a complete freeze). Refuses before any evaluation when a pose lacks
    its electronic state. ``hf`` (HFSettings) freezes rule D2-HF-v1; ``hf`` together with
    ``start_cluster`` (StartCluster) freezes D2-HF-v1.1; adding ``same_graph_ridge_last``
    freezes D2-HF-v1.2; adding ``hs_bracket`` (HSBracket) freezes D2-HF-v1.3; adding both
    ``same_graph_sample`` and ``hf_parallel`` freezes D2-HF-v1.4."""
    run_dir = Path(run_dir)
    if start_cluster is not None and hf is None:
        raise ValueError("start_cluster orders HS dispatch and needs hf")
    if not isinstance(same_graph_ridge_last, bool):
        raise ValueError("same_graph_ridge_last must be a bool")
    if same_graph_ridge_last and start_cluster is None:
        raise ValueError("same_graph_ridge_last extends D2-HF-v1.1 and needs start_cluster")
    if hs_bracket is not None:
        if not same_graph_ridge_last:
            raise ValueError("hs_bracket extends D2-HF-v1.2 and needs same_graph_ridge_last")
        _ef_settings(hf, hs_bracket)  # EFSettings checks the thresholds against hf
    for name, value in (("same_graph_sample", same_graph_sample), ("hf_parallel", hf_parallel)):
        if not isinstance(value, bool):
            raise ValueError(f"{name} must be a bool")
    if same_graph_sample != hf_parallel:
        raise ValueError("same_graph_sample and hf_parallel freeze D2-HF-v1.4 together")
    if same_graph_sample and hs_bracket is None:
        raise ValueError("D2-HF-v1.4 extends D2-HF-v1.3 and needs hs_bracket")
    if (run_dir / "manifest.json").exists():
        raise FileExistsError(f"{run_dir} already holds a D2 run; resume it with run()")
    for pose_id, atoms in poses:
        missing = [k for k in ("charge", "spin") if k not in atoms.info]
        if missing:
            raise ValueError(f"pose {pose_id} lacks atoms.info {missing}")
    rows = []
    for pose_id, atoms in poses:
        pose_key = _key(settings.seed, "pose_order", pose_id)
        specs = sorted(
            specs_by_pose.get(pose_id, ()),
            key=lambda sp: _key(settings.seed, "population", sp.intent.intent_id),
        )
        for position, spec in enumerate(specs):
            intent = spec.intent
            wanted = spec.deltas(settings.magnitudes_A)
            q0 = intent.value(atoms.positions)
            deltas = feasible_prefix(intent, q0, wanted, atoms.numbers, config)
            sampled = position < settings.per_pose_quota
            if not sampled:
                status = "not_sampled"
            else:
                status = "ready" if len(deltas) >= 2 else "infeasible"
            rows.append(
                {
                    "intent_id": intent.intent_id,
                    "pose_id": pose_id,
                    "pairs": [list(p) for p in intent.pairs],
                    "coefs": list(intent.coefs),
                    "sign": spec.sign,
                    "family": spec.family,
                    "layer": spec.layer,
                    "scope": spec.scope,
                    "deltas": [float(d) for d in deltas],
                    "domain_truncated": len(deltas) < len(wanted),
                    "sampled": sampled,
                    "status": status,
                    "order": [position, pose_key],
                    "audit": _key(settings.seed, "audit", intent.intent_id)
                    < settings.audit_fraction,
                }
            )
    rows.sort(key=lambda r: tuple(r["order"]))
    run_dir.mkdir(parents=True, exist_ok=True)
    with (run_dir / "poses.jsonl").open("w", encoding="utf-8") as fh:
        for pose_id, atoms in poses:
            record = {
                "pose_id": pose_id,
                "numbers": atoms.numbers.tolist(),
                "positions": atoms.positions.tolist(),
                "info": {"charge": int(atoms.info["charge"]), "spin": int(atoms.info["spin"])},
            }
            fh.write(json.dumps(record) + "\n")
    with (run_dir / "population.jsonl").open("w", encoding="utf-8") as fh:
        for row in rows:
            fh.write(json.dumps(row) + "\n")
    if hf is None:
        rule = RULE_VERSION
    elif start_cluster is None:
        rule = HF_RULE_VERSION
    else:
        rule = HF_SC_RULE_VERSION
        if same_graph_ridge_last:
            rule = HF_SG_RULE_VERSION if hs_bracket is None else HF_BR_RULE_VERSION
            if same_graph_sample:
                rule = HF_V14_RULE_VERSION
    manifest = {
        "rule_version": rule,
        "arm": ARM,
        "settings": asdict(settings),
        "hf": None if hf is None else asdict(hf),
        "start_cluster": None if start_cluster is None else asdict(start_cluster),
        "config": config.to_dict(),
        "identity": identity,
        "host": socket.gethostname(),
    }
    if same_graph_ridge_last:  # earlier rules keep their manifest unchanged
        manifest["same_graph_ridge_last"] = True
    if hs_bracket is not None:
        manifest["hs_bracket"] = asdict(hs_bracket)
    if same_graph_sample:
        manifest["same_graph_sample"] = manifest["hf_parallel"] = True
    (run_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")


def _task(kind, row):
    return f"{kind}|{row['intent_id']}"


def _seed(settings, task):
    return int(stream(settings.seed, "task", task).integers(2**31 - 1))


def _intent(row):
    return Intent(row["intent_id"], tuple(tuple(p) for p in row["pairs"]), tuple(row["coefs"]))


def _scan_settings(settings, deltas):
    return ScanSettings(
        tuple(float(d) for d in deltas),
        settings.fmax_eV_A,
        settings.max_steps_per_point,
        settings.max_rise_eV,
        reuse_curvature=True,
    )


def _scan_payload(outcome):
    if outcome.error is not None:
        return f"error:{outcome.error.split(':')[0]}", {"error": outcome.error}
    result = outcome.result
    rejected = result.rejected_positions
    return result.status, {
        "q0": result.q0,
        "points": [_jsonable(vars(p)) for p in result.points],
        "rejected_positions": None if rejected is None else rejected.tolist(),
    }


def _points(payload):
    out = []
    for d in payload.get("points", []):
        forces = d["forces_pes"]
        out.append(
            ScanPoint(
                **{
                    **d,
                    "positions": np.asarray(d["positions"], float),
                    "forces_pes": None if forces is None else np.asarray(forces, float),
                }
            )
        )
    return out


class _Run:
    """Run state rebuilt from the run directory. Construction reconciles every
    dispatched-but-unsettled attempt and journals it (see ``open_run``)."""

    def __init__(self, run_dir):
        from prrs.config import SearchConfig

        self.dir = Path(run_dir)
        manifest = json.loads((self.dir / "manifest.json").read_text(encoding="utf-8"))
        raw = dict(manifest["settings"])
        raw["magnitudes_A"] = tuple(raw["magnitudes_A"])
        self.settings = D2Settings(**raw)
        self.config = SearchConfig.from_dict(manifest["config"])
        self.rule_version = manifest["rule_version"]
        self.hf = HFSettings(**manifest["hf"]) if manifest.get("hf") else None
        sc = manifest.get("start_cluster")
        self.start_cluster = StartCluster(**sc) if sc else None
        self.start_clusters = None  # filled by _promote under D2-HF-v1.1
        self.same_graph_ridge_last = bool(manifest.get("same_graph_ridge_last", False))
        hb = manifest.get("hs_bracket")
        self.hs_bracket = HSBracket(**hb) if hb else None
        self.same_graph_sample = bool(manifest.get("same_graph_sample", False))
        self.hf_parallel = bool(manifest.get("hf_parallel", False))
        # In HF mode the second account pays E/F search + connection ("refine"); its
        # budget is the frozen certify_budget. Stages are told apart in the ledger.
        self.l4_account = "certify" if self.hf is None else "refine"
        self.poses = {}
        for line in (self.dir / "poses.jsonl").read_text(encoding="utf-8").splitlines():
            if line.strip():
                r = json.loads(line)
                atoms = Atoms(numbers=r["numbers"], positions=r["positions"])
                atoms.info.update(r["info"])
                self.poses[r["pose_id"]] = atoms
        self.population = [
            json.loads(line)
            for line in (self.dir / "population.jsonl").read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
        self.ledger_path = self.dir / "ledger.jsonl"
        repair_tail(self.ledger_path)
        self.journal = Journal(self.dir / "journal.jsonl")
        self.accounts = {
            "search": Account("search", self.settings.search_budget),
            self.l4_account: Account(self.l4_account, self.settings.certify_budget),
        }
        self.settled, self.skipped, self.attempts = {}, {}, {}
        unsettled = {}
        for ev in self.journal.events:
            kind = ev["event"]
            if kind == "dispatch":
                self.attempts[ev["task"]] = max(
                    self.attempts.get(ev["task"], 0), ev["attempt"] + 1
                )
                unsettled[(ev["task"], ev["attempt"])] = ev
            elif kind in ("settle", "reconcile"):
                unsettled.pop((ev["task"], ev["attempt"]), None)
                self.accounts[ev["account"]].charge(ev["spent"])
                if kind == "settle":
                    self.settled[ev["task"]] = ev
            elif kind == "skip":
                self.skipped[ev["task"]] = ev["reason"]
        self._reconcile(unsettled)

    def _reconcile(self, unsettled):
        counts = tally(self.ledger_path)
        for (task, attempt), ev in unsettled.items():
            row = counts.get((task, attempt), {"forwards": 0, "open": 0})
            exact = row["open"] == 0
            spent = row["forwards"] if exact else ev["cap"]
            self.journal.append(
                {
                    "event": "reconcile",
                    "task": task,
                    "attempt": attempt,
                    "account": ev["account"],
                    "cap": ev["cap"],
                    "spent": spent,
                    "reconciled": exact,
                    "known_forwards": row["forwards"],
                }
            )
            self.accounts[ev["account"]].charge(spent)

    def total_spent(self):
        return sum(a.spent for a in self.accounts.values())

    def task_ledger(self, cap):
        return Ledger(self.ledger_path, arm=ARM, budget=cap, replay=False)

    def dispatch(self, task, account, cap, kind):
        if not self.accounts[account].hold(cap):
            return None
        attempt = self.attempts.get(task, 0)
        self.attempts[task] = attempt + 1
        self.journal.append(
            {
                "event": "dispatch",
                "task": task,
                "attempt": attempt,
                "account": account,
                "cap": cap,
                "kind": kind,
            }
        )
        return attempt

    def settle(self, task, attempt, account, cap, ledger, status, payload):
        ledger.close()
        self.accounts[account].settle(cap, ledger.spent)
        event = {
            "event": "settle",
            "task": task,
            "attempt": attempt,
            "account": account,
            "cap": cap,
            "spent": ledger.spent,
            "status": status,
            "payload": payload,
            "cumulative_spent": self.total_spent(),
        }
        self.journal.append(event)
        self.settled[task] = event

    def skip(self, task, reason):
        if task in self.skipped:
            return
        self.journal.append(
            {
                "event": "skip",
                "task": task,
                "reason": reason,
                "cumulative_spent": self.total_spent(),
            }
        )
        self.skipped[task] = reason

    def pending(self, task):
        return task not in self.settled and task not in self.skipped

    def known_c2(self):
        """Certified saddles so far, in journal order, for certify's duplicate rule."""
        out = []
        for task, ev in self.settled.items():
            if task.startswith("L4|") and ev["status"] == "c2":
                p = ev["payload"]
                cert = Certification(
                    status="c2",
                    saddle_positions=np.asarray(p["saddle_positions"], float),
                    saddle_energy_eV=p["saddle_energy_eV"],
                    chemical_pair=tuple(p["chemical_pair"]),
                )
                out.append((task, cert))
        return out


def open_run(run_dir):
    """Load a frozen run and reconcile interrupted attempts (journals reconcile events)."""
    return _Run(run_dir)


def _run_scans(run, backends, items, cap, kind):
    """Batched scans in waves of ``wave_size``; each job holds its cap from the search
    account before dispatch and settles its real spend after the wave."""
    s = run.settings
    todo = list(items)
    while todo:
        wave = []
        for task, atoms, intent, scan_settings in todo[: s.wave_size]:
            attempt = run.dispatch(task, "search", cap, kind)
            if attempt is None:
                break
            job = ScanJob(
                atoms,
                intent,
                scan_settings,
                run.config,
                run.task_ledger(cap),
                task,
                attempt=attempt,
            )
            wave.append((task, attempt, job))
        if not wave:
            for task, *_ in todo:
                run.skip(task, f"{kind}_skipped:search_budget")
            return
        outcomes, _ = scan_many(
            [job for _, _, job in wave], backends.search_batch, s.max_active
        )
        for (task, attempt, job), outcome in zip(wave, outcomes):
            status, payload = _scan_payload(outcome)
            run.settle(task, attempt, "search", cap, job.ledger, status, payload)
        todo = todo[len(wave) :]


def _l1(run, backends):
    s = run.settings
    items = []
    for row in run.population:
        task = _task("L1", row)
        if row["status"] == "ready" and run.pending(task):
            atoms = run.poses[row["pose_id"]].copy()
            items.append((task, atoms, _intent(row), _scan_settings(s, row["deltas"])))
    _run_scans(run, backends, items, s.scan_cap, "L1")


def _first(triggers, kind):
    return next((t for t in triggers if t["kind"] == kind), None)


def _l2(run, backends):
    s = run.settings
    items = []
    for row in run.population:
        ev = run.settled.get(_task("L1", row))
        task = _task("L2", row)
        if ev is None or not run.pending(task):
            continue
        points = _points(ev["payload"])
        triggers = read_scan(points, s.noise_eV)
        jump = _first(triggers, "jump")
        if _first(triggers, "ridge_suspected") is not None or jump is None:
            continue
        start, back = reverse_plan(points, jump, _scan_settings(s, row["deltas"]))
        atoms = run.poses[row["pose_id"]].copy()
        atoms.positions = start
        items.append((task, atoms, _intent(row), back))
    _run_scans(run, backends, items, s.reverse_cap, "L2")


_CLASS_RANK = {"ridge_suspected": 0, "branch_bracketed": 1, "audit": 2}


def _classify(run, row):
    """Promotion class or deferral reason of one scanned intent (pure, from the journal).
    Audit eligibility = deferred with an L1 response (at least two points), drawn
    once at freeze time (``row["audit"]``)."""
    ev = run.settled[_task("L1", row)]
    points = _points(ev["payload"])
    if len(points) < 2:
        return {"class": None, "info": None, "reason": f"no_response:{ev['status']}"}
    triggers = read_scan(points, run.settings.noise_eV)
    ridge = _first(triggers, "ridge_suspected")
    if ridge is not None:
        return {"class": "ridge_suspected", "info": ridge, "reason": None}
    jump = _first(triggers, "jump")
    reason = "no_trigger"
    if jump is not None:
        rev = run.settled.get(_task("L2", row))
        if rev is None:
            reason = "jump_unchecked"
        else:
            back = _points(rev["payload"]) if rev["status"] == "completed" else []
            bracket = bracket_jump(points, jump, back[-1]) if len(back) >= 2 else None
            if bracket is not None:
                # the scan's own point after the jump: the product branch end of the
                # forward jump interval, which the HF search continues (_hf_bracket)
                info = {**bracket, "product_side_scan": jump["after"]}
                return {"class": "branch_bracketed", "info": info, "reason": None}
            reason = "jump_not_bracketed"
    if row["audit"]:
        return {"class": "audit", "info": None, "reason": reason}
    return {"class": None, "info": None, "reason": reason}


def _promotion_points(classified, points):
    """(release positions, certification guess, scan tangent at the guess) by the frozen
    rule M5-D1. The tangent is the Cartesian displacement of the scan across the guess
    (F6, TS_HESSIAN_FREE_DESIGN_RECORD_2026-10-24 §4): it is handed to the refiner so it
    follows the scanned coordinate instead of the softest mode."""
    info = classified["info"]
    if classified["class"] == "ridge_suspected":
        k = info["peak"]
        tangent = points[k + 1].positions - points[k - 1].positions
        return points[info["release"]].positions, points[k].positions, tangent
    if classified["class"] == "branch_bracketed":
        product = np.asarray(info["product_side_positions"], float)
        reactant = points[info["reactant_side"]].positions
        return product, reactant, product - reactant
    good = [k for k, p in enumerate(points) if p.converged] or list(range(len(points)))
    interior = [k for k in good if k >= 1] or [len(points) - 1]
    peak = max(interior, key=lambda k: points[k].energy_eV)
    after = points[min(peak + 1, len(points) - 1)].positions
    before = points[max(peak - 1, 0)].positions
    return points[good[-1]].positions, points[peak].positions, after - before


def _source(template, point):
    atoms = template.copy()
    atoms.positions = point.positions
    atoms.calc = SinglePointCalculator(atoms, energy=point.energy_eV, forces=point.forces_pes)
    return atoms


def _release_json(out):
    end = out.get("endpoint")
    data = {k: v for k, v in out.items() if k != "endpoint"}
    data["endpoint_positions"] = None if end is None else end.positions.tolist()
    data["endpoint_energy_eV"] = None if end is None else float(end.get_potential_energy())
    return _jsonable(data)


def _release_task(run, backends, row, start, positions, tags):
    s = run.settings
    task = _task("L3", row)
    if not run.pending(task):
        return
    attempt = run.dispatch(task, "search", s.release_cap, "L3")
    if attempt is None:
        run.skip(task, "L3_skipped:search_budget")
        return
    ledger = run.task_ledger(s.release_cap)
    template = run.poses[row["pose_id"]]
    try:
        out = release(
            positions,
            template,
            _source(template, start),
            view_factory(backends.certify_calculator, ledger),
            run.config,
            ledger,
            task,
            seed=_seed(s, task),
            attempt=attempt,
        )
        status, payload = out["status"], _release_json(out)
    except Exception as exc:
        status, payload = (
            f"error:{type(exc).__name__}",
            {"error": str(exc), "traceback": traceback.format_exc()},
        )
    run.settle(task, attempt, "search", s.release_cap, ledger, status, {**payload, **tags})


def _certify_task(run, backends, row, start, positions, tags, direction=None):
    s = run.settings
    task = _task("L4", row)
    if not run.pending(task):
        return
    attempt = run.dispatch(task, "certify", s.certify_cap, "L4")
    if attempt is None:
        run.skip(task, "L4_skipped:certify_budget")
        return
    ledger = run.task_ledger(s.certify_cap)
    guess = run.poses[row["pose_id"]].copy()
    guess.positions = np.asarray(positions, float)
    try:
        cert = certify(
            guess,
            view_factory(backends.certify_calculator, ledger),
            run.config,
            ledger,
            task,
            seed=_seed(s, task),
            ceiling_eV=start.energy_eV + s.certify_max_rise_eV,
            reference_eV=start.energy_eV,
            direction=direction,
            known=run.known_c2(),
            attempt=attempt,
        )
        status, payload = cert.status, _jsonable(vars(cert))
    except Exception as exc:
        status, payload = (
            f"error:{type(exc).__name__}",
            {"error": str(exc), "traceback": traceback.format_exc()},
        )
    run.settle(task, attempt, "certify", s.certify_cap, ledger, status, {**payload, **tags})


def _hf_bracket(classified, points):
    """Two bracket ends for the HF search (PLAN_HF4 §2.1); neither end's saved flags are
    trusted, continue_branches_ef re-relaxes both.

    A bracketed jump is searched on the forward scan's own jump interval: the last point
    on the reactant branch and the first on the product branch, at two different q, each
    continued on its own branch (HF record §5.2). The L2 reverse point is not an end: it
    sits at the reactant point's q and only shows that both branches coexist there --
    using it gave a zero-width interval (b05 HF v1.1, 11 HS errors, 2026-10-28)."""
    info = classified["info"]
    if classified["class"] == "ridge_suspected":
        k = info["peak"]
        return points[k - 1], points[k + 1]
    if classified["class"] == "branch_bracketed":
        return points[info["reactant_side"]], points[info["product_side_scan"]]
    good = [k for k, p in enumerate(points) if p.converged] or list(range(len(points)))
    interior = [k for k in good if k >= 1] or [len(points) - 1]
    peak = max(interior, key=lambda k: points[k].energy_eV)
    return points[max(peak - 1, 0)], points[min(peak + 1, len(points) - 1)]


def _lean(frame):
    """Search frame without per-frame forces: positions stay (connection starts use
    them), forces would multiply the journal size by about two."""
    return {k: v for k, v in frame.items() if k != "forces_pes"}


def _hs_payload(result):
    data = result.to_dict()
    data["history"] = [_lean(f) for f in data["history"]]
    for side in data["sides"]:
        side["points"] = [_lean(f) for f in side["points"]]
    if data.get("refine"):
        data["refine"]["history"] = [_lean(f) for f in data["refine"]["history"]]
    return data


def _ef_settings(hf, hs_bracket=None):
    settings = EFSettings(
        reaction_force_sigma=hf.reaction_force_sigma,
        midpoint_transverse_fmax_eV_A=hf.midpoint_transverse_fmax_eV_A,
    )
    if hs_bracket is None:
        return settings
    return replace(
        settings,
        midpoint_coarse_transverse_fmax_eV_A=hs_bracket.midpoint_coarse_transverse_fmax_eV_A,
        midpoint_sign_bound_eV_A=hs_bracket.midpoint_sign_bound_eV_A,
        midpoint_fine_optimizer="fire",
        midpoint_coarse_candidate=True,
        endpoint_candidate=True,
    )


def _hs_run(run, row, start, bracket, factory, ledger, task, attempt):
    """One HS search with any counted factory; (status, payload), errors kept."""
    s, hf = run.settings, run.hf
    try:
        result = continue_branches_ef(
            run.poses[row["pose_id"]].copy(),
            _intent(row),
            bracket,
            factory,
            run.config,
            ledger,
            task,
            branch_settings=BranchSettings(geometry_noise_A=hf.geometry_noise_A),
            settings=_ef_settings(hf, run.hs_bracket),
            ceiling_eV=start.energy_eV + s.certify_max_rise_eV,
            attempt=attempt,
        )
        return result.status, _hs_payload(result)
    except Exception as exc:
        return (
            f"error:{type(exc).__name__}",
            {"error": str(exc), "traceback": traceback.format_exc()},
        )


def _hs_task(run, backends, row, start, bracket, tags, skip_reason="HS_skipped:refine_budget"):
    hf = run.hf
    task = _task("HS", row)
    if not run.pending(task):
        return
    attempt = run.dispatch(task, "refine", hf.hs_cap, "HS")
    if attempt is None:
        run.skip(task, skip_reason)
        return
    ledger = run.task_ledger(hf.hs_cap)
    factory = view_factory(backends.certify_calculator, ledger)
    status, payload = _hs_run(run, row, start, bracket, factory, ledger, task, attempt)
    run.settle(task, attempt, "refine", hf.hs_cap, ledger, status, {**payload, **tags})


def _hf_identity(run, template, positions):
    from prrs.chemistry import chemical_key

    atoms = template.copy()
    atoms.positions = np.asarray(positions, float)
    c = run.config
    return chemical_key(
        atoms,
        c.bond_scale,
        c.active_atoms,
        charge_sensitive=c.charge_sensitive,
        charge=c.total_charge,
        multiplicity=c.multiplicity,
        parity_tolerance=c.parity_tolerance,
    )["key"]


def _hf_duplicate(run, template, saddle, energy, pair, ends):
    """First earlier connected HF event with the same unordered chemical pair, saddle
    energy within basin_energy_eV and symmetric RMSD within basin_rmsd_A under the saddle
    and endpoint symmetries (the duplicate rule of certify step 5, Hessian-free)."""
    from prrs.chemistry import symmetric_rmsd

    c, group = run.config, None
    atoms = template.copy()
    atoms.positions = np.asarray(saddle, float)
    for task, ev in run.settled.items():
        if not task.startswith("HC|") or ev["status"] != "connected":
            continue
        p = ev["payload"]
        if tuple(p["chemical_pair"]) != pair:
            continue
        if abs(p["saddle_energy_eV"] - energy) > c.basin_energy_eV:
            continue
        if group is None:
            minima = []
            for positions in ends:
                end = template.copy()
                end.positions = np.asarray(positions, float)
                minima.append(end)
            group, _ = pair_symmetries(atoms, minima, c)
        other = atoms.copy()
        other.positions = np.asarray(p["saddle_positions"], float)
        if symmetric_rmsd(atoms, other, group, c.active_atoms) <= c.basin_rmsd_A:
            return task
    return None


def _hc_ready(run, row):
    hs = run.settled.get(_task("HS", row))
    return (
        hs is not None
        and hs["status"] == "converged_candidate"
        and run.pending(_task("HC", row))
    )


def _hc_run(run, row, factory, ledger, task, attempt):
    """Connection trace of a settled HS candidate; (status, payload) before the
    duplicate check, which needs the settled events and runs in ``_hc_finish``."""
    hf = run.hf
    hs = run.settled[_task("HS", row)]
    template = run.poses[row["pose_id"]]
    found = hs["payload"]["refine"]
    try:
        frames = list(found["history"])
        for side in hs["payload"]["sides"]:
            if side.get("status") == "maximum_bracketed":
                frames += side["points"]
        starts = branch_starts(
            frames,
            found["metrics"]["q_actual"],
            hf.reaction_force_sigma,
            transverse_fmax=_ef_settings(hf).transverse_fmax_eV_A,
            descent_fmax=run.config.quench_fmax_eV_A,
        )
        result = trace_sides_ef(
            template.copy(),
            found["positions"],
            found["energy_eV"],
            starts,
            found["direction_cartesian"],
            factory,
            run.config,
            ledger,
            task,
            settings=ConnectSettings(
                reaction_force_sigma=hf.reaction_force_sigma, seed=_seed(run.settings, task)
            ),
            attempt=attempt,
        )
        status, payload = result.status, result.to_dict()
        payload["saddle_positions"] = found["positions"]
        payload["saddle_energy_eV"] = found["energy_eV"]
        return status, payload
    except Exception as exc:
        return (
            f"error:{type(exc).__name__}",
            {"error": str(exc), "traceback": traceback.format_exc()},
        )


def _hc_finish(run, row, status, payload):
    """Chemical pair and the duplicate rule against events settled so far."""
    if status != "connected":
        return status, payload
    template = run.poses[row["pose_id"]]
    found = run.settled[_task("HS", row)]["payload"]["refine"]
    try:
        ends = {sd["side"]: sd["descent"]["positions"] for sd in payload["sides"]}
        pair = tuple(sorted(_hf_identity(run, template, ends[k]) for k in ("low", "high")))
        payload["chemical_pair"] = list(pair)
        dup = _hf_duplicate(
            run, template, found["positions"], found["energy_eV"], pair, ends.values()
        )
        if dup is not None:
            status, payload["duplicate_of"] = "hf_duplicate", dup
    except Exception as exc:
        status, payload = (
            f"error:{type(exc).__name__}",
            {"error": str(exc), "traceback": traceback.format_exc()},
        )
    return status, payload


def _hc_task(run, backends, row, tags):
    hf = run.hf
    task = _task("HC", row)
    if not _hc_ready(run, row):
        return
    attempt = run.dispatch(task, "refine", hf.hc_cap, "HC")
    if attempt is None:
        run.skip(task, "HC_skipped:refine_budget")
        return
    ledger = run.task_ledger(hf.hc_cap)
    factory = view_factory(backends.certify_calculator, ledger)
    status, payload = _hc_run(run, row, factory, ledger, task, attempt)
    status, payload = _hc_finish(run, row, status, payload)
    run.settle(task, attempt, "refine", hf.hc_cap, ledger, status, {**payload, **tags})


def _job(fn, *args):
    """TaskJob runner: ``fn(..., factory, ledger, task, attempt)`` with run_many's factory."""
    run, row, *rest, ledger, task, attempt = args

    def go(factory):
        return fn(run, row, *rest, factory, ledger, task, attempt)

    return TaskJob(go, ledger)


def _outcome(outcome):
    if outcome.error is not None:  # failures outside the task's own try
        return f"error:{outcome.error.split(':', 1)[0]}", {"error": outcome.error}
    return outcome.result


def _hf_parallel(run, backends, items, skip_reason):
    """HS / HC of one dispatch group in waves (D2-HF-v1.4). ``items`` are
    (row, points, bracket, tags) in dispatch order. Each wave holds its caps before
    dispatch, runs through the resident batch evaluator and settles in input order;
    an item the budget cannot reach waits for the next wave, and only a wave that can
    dispatch nothing skips what is left (with the serial rules' reasons)."""
    s, hf = run.settings, run.hf

    def hc_waves():
        while True:
            ready = [it for it in items if _hc_ready(run, it[0])]
            if not ready:
                return
            wave = []
            for row, _, _, tags in ready[: s.wave_size]:
                task = _task("HC", row)
                attempt = run.dispatch(task, "refine", hf.hc_cap, "HC")
                if attempt is None:
                    break
                wave.append((row, tags, task, attempt, run.task_ledger(hf.hc_cap)))
            if not wave:
                for row, *_ in ready:
                    run.skip(_task("HC", row), "HC_skipped:refine_budget")
                return
            jobs = [_job(_hc_run, run, row, led, task, att) for row, _, task, att, led in wave]
            outcomes, _ = run_many(jobs, backends.search_batch, s.max_active)
            for (row, tags, task, attempt, ledger), out in zip(wave, outcomes):
                status, payload = _hc_finish(run, row, *_outcome(out))
                run.settle(
                    task, attempt, "refine", hf.hc_cap, ledger, status, {**payload, **tags}
                )

    while True:
        pending = [it for it in items if run.pending(_task("HS", it[0]))]
        wave = []
        for row, points, bracket, tags in pending[: s.wave_size]:
            task = _task("HS", row)
            attempt = run.dispatch(task, "refine", hf.hs_cap, "HS")
            if attempt is None:
                break
            ledger = run.task_ledger(hf.hs_cap)
            wave.append((row, points, bracket, tags, task, attempt, ledger))
        if wave:
            jobs = [
                _job(_hs_run, run, row, points[0], bracket, ledger, task, attempt)
                for row, points, bracket, _, task, attempt, ledger in wave
            ]
            outcomes, _ = run_many(jobs, backends.search_batch, s.max_active)
            for (row, _, _, tags, task, attempt, ledger), out in zip(wave, outcomes):
                status, payload = _outcome(out)
                run.settle(
                    task, attempt, "refine", hf.hs_cap, ledger, status, {**payload, **tags}
                )
        else:
            for row, *_ in pending:
                run.skip(_task("HS", row), skip_reason)
        hc_waves()
        if not wave:
            return


def _same_graph_ridge(classified, points):
    """A ridge whose scan peak and the next point share their bond graph (D2-HF-v1.2):
    the rise-fall-flip of ``read_scan`` met on a conformational jump (F4 of
    TS_HESSIAN_FREE_DESIGN_RECORD_2026-10-24). Same ``edges`` as the jump trigger."""
    if classified["class"] != "ridge_suspected":
        return False
    k = classified["info"]["peak"]
    return points[k].edges == points[k + 1].edges


def _promote(run, backends):
    s = run.settings
    chosen = []
    for row in run.population:
        if _task("L1", row) in run.settled:
            classified = _classify(run, row)
            if classified["class"] is not None:
                chosen.append((row, classified))
    chosen.sort(
        key=lambda rc: (
            _CLASS_RANK[rc[1]["class"]],
            _key(s.seed, "promote", rc[0]["intent_id"]),
        )
    )
    late = set()
    if run.same_graph_ridge_last:
        late = {
            row["intent_id"]
            for row, classified in chosen
            if _same_graph_ridge(classified, _points(run.settled[_task("L1", row)]["payload"]))
        }
    clusters = first_round = None
    if run.start_cluster is not None:
        # Same-graph ridges join clusters last, so none of them leads one and pushes a
        # bond-changing ridge into round two.
        ordered = [rc for rc in chosen if rc[0]["intent_id"] not in late]
        ordered += [rc for rc in chosen if rc[0]["intent_id"] in late]
        clusters = run.start_clusters = _start_clusters(run, ordered)
        first_round = run.start_cluster.first_round
    second, last, first = [], [], []
    for row, classified in chosen:
        points = _points(run.settled[_task("L1", row)]["payload"])
        release_at, guess_at, tangent = _promotion_points(classified, points)
        tags = _tags(classified)
        _release_task(run, backends, row, points[0], release_at, tags)
        if run.hf is None:
            _certify_task(run, backends, row, points[0], guess_at, tags, direction=tangent)
            continue
        if run.same_graph_ridge_last:
            tags = {**tags, "same_graph_ridge": row["intent_id"] in late}
            if row["intent_id"] in late:
                if run.same_graph_sample:
                    drawn = _key(s.seed, "same_graph_sample", row["intent_id"])
                    tags["same_graph_sample"] = drawn < s.audit_fraction
                    if not tags["same_graph_sample"]:
                        if run.pending(_task("HS", row)):
                            run.skip(_task("HS", row), "HS_skipped:same_graph_ridge")
                        continue
                last.append((row, classified, points, tags))
                continue
        if clusters is not None and clusters[row["intent_id"]]["rank"] >= first_round:
            second.append((row, classified, points) + ((tags,) if run.hf_parallel else ()))
            continue
        if run.hf_parallel:
            first.append((row, classified, points, tags))
            continue
        _hf_tasks(run, backends, row, classified, points, tags, clusters)
    if run.hf is not None and run.hf_parallel:
        # D2-HF-v1.4: the same three groups, each in parallel waves
        for group, reason in (
            (first, "HS_skipped:refine_budget"),
            (second, "HS_deferred:same_start"),
            (last, "HS_deferred:same_graph_ridge"),
        ):
            items = []
            for row, classified, points, tags in group:
                hf_tags = {k: v for k, v in tags.items() if k != "direction"}
                if clusters is not None:
                    hf_tags["start_cluster"] = clusters[row["intent_id"]]
                items.append((row, points, _hf_bracket(classified, points), hf_tags))
            _hf_parallel(run, backends, items, reason)
        return
    # round two (D2-HF-v1.1): the rest of every start cluster, on what round one left
    for row, classified, points in second:
        _hf_tasks(
            run,
            backends,
            row,
            classified,
            points,
            _tags(classified),
            clusters,
            skip_reason="HS_deferred:same_start",
        )
    # D2-HF-v1.2: same-graph ridges, on what both rounds left
    for row, classified, points, tags in last:
        _hf_tasks(
            run,
            backends,
            row,
            classified,
            points,
            tags,
            clusters,
            skip_reason="HS_deferred:same_graph_ridge",
        )


def _tags(classified):
    return {
        "class": classified["class"],
        "audit_release": classified["class"] == "audit",
        "deferred_reason": classified["reason"],
        "direction": "scan_tangent",
    }


def _hf_tasks(run, backends, row, classified, points, tags, clusters, skip_reason=None):
    hf_tags = {k: v for k, v in tags.items() if k != "direction"}
    if clusters is not None:
        hf_tags["start_cluster"] = clusters[row["intent_id"]]
    bracket = _hf_bracket(classified, points)
    extra = {} if skip_reason is None else {"skip_reason": skip_reason}
    _hs_task(run, backends, row, points[0], bracket, hf_tags, **extra)
    _hc_task(run, backends, row, hf_tags)


def _rmsd_lower_bound(a, b, index):
    """Lower bound of the Kabsch RMSD between the ``index`` blocks of a and b over every
    like-atom relabelling: about the centroids |x_i - R y_j| >= | |x_i| - |y_j| | for any
    rotation R, and matching sorted radii per element minimises the sum of squares."""
    x, y = a.positions[index], b.positions[index]
    rx = np.linalg.norm(x - x.mean(axis=0), axis=1)
    ry = np.linalg.norm(y - y.mean(axis=0), axis=1)
    z = a.numbers[index]
    total = sum(np.sum((np.sort(rx[z == e]) - np.sort(ry[z == e])) ** 2) for e in np.unique(z))
    return float(np.sqrt(total / len(index)))


def _start_clusters(run, chosen):
    """Start cluster of every promoted item, in promotion order (D2-HF-v1.1). The
    starting guess is the M5-D1 guess point (``_promotion_points``); an item joins the
    first earlier leader within noise_multiple x geometry_noise_A, compared by symmetric
    RMSD under the union of both poses' graph automorphisms. Pure function of the
    journal, so a resumed run rebuilds the same clusters."""
    from prrs.chemistry import symmetric_rmsd
    from prrs.state import resolve_active

    c = run.config
    radius = run.start_cluster.noise_multiple * run.hf.geometry_noise_A
    groups, leaders, out = {}, [], {}
    for row, classified in chosen:
        points = _points(run.settled[_task("L1", row)]["payload"])
        atoms = run.poses[row["pose_id"]].copy()
        atoms.positions = np.asarray(_promotion_points(classified, points)[1], float)
        index = resolve_active(atoms, c.active_atoms)
        hit = None
        for leader in leaders:
            other = leader["atoms"]
            if not np.array_equal(atoms.numbers, other.numbers):
                continue
            if _rmsd_lower_bound(atoms, other, index) > radius:
                continue
            key = tuple(sorted((row["pose_id"], leader["pose"])))
            if key not in groups:
                groups[key], _ = pair_symmetries(run.poses[key[0]], (run.poses[key[1]],), c)
            if symmetric_rmsd(other, atoms, groups[key], c.active_atoms) <= radius:
                hit = leader
                break
        if hit is None:
            hit = {"atoms": atoms, "pose": row["pose_id"], "id": row["intent_id"], "size": 0}
            leaders.append(hit)
        out[row["intent_id"]] = {"leader": hit["id"], "rank": hit["size"]}
        hit["size"] += 1
    return out


def _counts(run, prefix):
    return dict(Counter(ev["status"] for t, ev in run.settled.items() if t.startswith(prefix)))


def _report(run):
    rows = run.population
    promotion, deferred = Counter(), Counter()
    for row in rows:
        if _task("L1", row) in run.settled:
            classified = _classify(run, row)
            if classified["class"] is not None:
                promotion[classified["class"]] += 1
            else:
                deferred[classified["reason"]] += 1
    population = dict(Counter(r["status"] for r in rows))
    population["domain_truncated"] = sum(
        1 for r in rows if r["sampled"] and r["domain_truncated"]
    )
    events = run.journal.events
    report = {
        "rule_version": run.rule_version,
        "population": population,
        "l1_status": _counts(run, "L1|"),
        "l2_status": _counts(run, "L2|"),
        "promotion": dict(promotion),
        "deferred": dict(deferred),
        "release_status": _counts(run, "L3|"),
        "certify_status": _counts(run, "L4|"),
        "c1": [
            {
                "task": t,
                "label": ev["payload"].get("label"),
                "class": ev["payload"]["class"],
                "audit_release": ev["payload"]["audit_release"],
            }
            for t, ev in run.settled.items()
            if t.startswith("L3|") and ev["status"] == "c1_candidate"
        ],
        "c1_by_label": dict(
            Counter(
                ev["payload"].get("label")
                for t, ev in run.settled.items()
                if t.startswith("L3|") and ev["status"] == "c1_candidate"
            )
        ),
        "c2": [
            {
                "task": t,
                "cumulative_spent": ev["cumulative_spent"],
                "saddle_energy_eV": ev["payload"]["saddle_energy_eV"],
                "chemical_pair": ev["payload"]["chemical_pair"],
                "event": ev["payload"]["event"],
                "class": ev["payload"]["class"],
                "audit_release": ev["payload"]["audit_release"],
            }
            for t, ev in run.settled.items()
            if t.startswith("L4|") and ev["status"] == "c2"
        ],
        "attempts": {t: n for t, n in run.attempts.items() if n > 1},
        "hosts": sorted({e["host"] for e in events if "host" in e}),
        "skips": dict(Counter(e["reason"] for e in events if e["event"] == "skip")),
        "accounts": {
            n: {"budget": a.budget, "spent": a.spent, "held": a.held}
            for n, a in run.accounts.items()
        },
        "ledger_total": sum(r["forwards"] for r in tally(run.ledger_path).values()),
        "unreconciled": [
            e for e in events if e["event"] == "reconcile" and not e["reconciled"]
        ],
    }
    if run.hf is not None:
        report.update(_hf_report(run))
    return report


def _hf_report(run):
    """Search / connection / property dimensions kept apart (HF record §6)."""
    stages = Counter()
    for (task, _attempt), row in tally(run.ledger_path).items():
        if task.startswith(("HS|", "HC|")):
            stages[task[:2]] += row["forwards"]
    connected = [
        (t, ev)
        for t, ev in run.settled.items()
        if t.startswith("HC|") and ev["status"] == "connected"
    ]
    out = {}
    if run.same_graph_ridge_last:
        tagged = [
            ev
            for t, ev in run.settled.items()
            if t.startswith("HS|") and ev["payload"].get("same_graph_ridge")
        ]
        out["same_graph_ridge"] = {
            "settled": len(tagged),
            "deferred": sum(r == "HS_deferred:same_graph_ridge" for r in run.skipped.values()),
        }
        if run.same_graph_sample:
            out["same_graph_ridge"] |= {
                "skipped": sum(
                    r == "HS_skipped:same_graph_ridge" for r in run.skipped.values()
                ),
                "sample_candidates": sorted(
                    t
                    for t, ev in run.settled.items()
                    if t.startswith("HS|")
                    and ev["payload"].get("same_graph_ridge")
                    and ev["status"] == "converged_candidate"
                ),
            }
    clusters = run.start_clusters
    if clusters is not None:
        first = run.start_cluster.first_round
        out["start_clusters"] = {
            "clusters": len({c["leader"] for c in clusters.values()}),
            "first_round": sum(c["rank"] < first for c in clusters.values()),
            "second_round": sum(c["rank"] >= first for c in clusters.values()),
            "radius_A": run.start_cluster.noise_multiple * run.hf.geometry_noise_A,
        }
    return out | {
        "hs_status": _counts(run, "HS|"),
        "hc_status": _counts(run, "HC|"),
        "property_status": "not_run",
        "hf_forwards": dict(stages),
        "hf_events": [
            {
                "task": t,
                "cumulative_spent": ev["cumulative_spent"],
                "saddle_energy_eV": ev["payload"]["saddle_energy_eV"],
                "chemical_pair": ev["payload"]["chemical_pair"],
                "label": ev["payload"]["label"],
                "broken": (ev["payload"]["fingerprint"] or {}).get("broken"),
                "formed": (ev["payload"]["fingerprint"] or {}).get("formed"),
                "class": ev["payload"]["class"],
            }
            for t, ev in connected
        ],
    }


def run(run_dir, backends):
    """Run (or resume) a frozen D2 run to completion; write and return report.json.
    Phases are strict (design §4.2): every sampled intent's L1 scan is settled or
    skipped before any L2/L3/L4 work starts."""
    state = open_run(run_dir)
    _l1(state, backends)
    _l2(state, backends)
    _promote(state, backends)
    report = _report(state)
    (state.dir / "report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report
