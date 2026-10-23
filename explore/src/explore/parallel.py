"""Independent controllers sharing one resident, batched energy/force evaluator.

Each task owns its calculator cache, guard, optimizer and ledger. Only the caller
thread touches the model. At each round it batches one pending structure from each
active controller; completed controllers are replaced from the input queue. No
timed batching window, shared ASE calculator or shared ledger segment is used.
Analytic-Hessian requests, when a backend is supplied, run one at a time on the
same coordinator thread.
"""

from __future__ import annotations

import threading
import time
from collections import deque
from concurrent.futures import Future, ThreadPoolExecutor
from dataclasses import dataclass
from typing import Callable

import numpy as np
from ase.calculators.calculator import Calculator, all_changes

from .ledger import BudgetExhausted
from .probe import constrained_scan

_STOPPED = "batch coordinator stopped"


@dataclass
class ScanJob:
    atoms: object
    intent: object
    settings: object
    config: object
    ledger: object
    task: str
    attempt: int = 0
    reserve: int | None = None
    artifact: str | None = None


@dataclass
class ScanOutcome:
    result: object = None
    error: str | None = None


@dataclass
class TaskJob:
    """``run(factory) -> result``; each ``factory()`` call returns a new client."""

    run: Callable
    ledger: object


@dataclass
class TaskOutcome:
    result: object = None
    error: str | None = None


class _ClientCalculator(Calculator):
    implemented_properties = ["energy", "forces"]

    def __init__(self, request, ledger, has_hessian=False, units=(1.0, 1.0)):
        super().__init__()
        self.request = request
        self.ledger = ledger
        self.evaluations = 0
        # The resident model's factors, so prrs._analytic_hessian scales a raw
        # Hessian from this client exactly as it would from the model itself.
        self.energy_units_to_eV, self.length_units_to_A = (float(u) for u in units)
        if has_hessian:
            # Attached only when a backend exists: prrs probes hasattr(get_hessian).
            self.get_hessian = self._get_hessian

    def get_property(self, name, atoms=None, allow_calculation=True):
        self.ledger.note_request()
        before = self.evaluations
        value = super().get_property(name, atoms, allow_calculation)
        if self.evaluations == before:
            self.ledger.note_cache_hit()
        return value

    def calculate(
        self, atoms=None, properties=("energy", "forces"), system_changes=all_changes
    ):
        super().calculate(atoms, properties, system_changes)
        self.evaluations += 1
        frame = atoms.copy()
        frame.set_constraint()
        energy, forces = self.request("forces", frame)
        self.results = {"energy": float(energy), "forces": np.asarray(forces).copy()}

    def _get_hessian(self, atoms=None):
        if atoms is None:
            atoms = self.atoms
        frame = atoms.copy()
        frame.set_constraint()
        return np.asarray(self.request("hessian", frame))


def run_many(jobs, evaluate_batch, max_active, hessian=None, units=(1.0, 1.0)):
    """Run independent tasks, returning input-ordered outcomes and round timings.

    ``job.run(factory)`` may call ``factory()`` any number of times; each call
    returns a new client calculator whose forwards are charged to ``job.ledger``.
    ``evaluate_batch(frames)`` returns one (energy_eV, forces_eV_A) per frame;
    it must preserve each frame's electronic state. Supply a resident evaluator
    object to reuse its model across successive calls. Budgets are per job; a
    batch of m structures still pays m structure forwards, including failed calls.
    ``hessian(atoms) -> (3N, 3N)`` (raw model units, see ``units``), when given,
    gives every client ``get_hessian``; each request runs alone on the coordinator
    thread and is charged one forward just before the call. General
    backend/task exceptions are retained in TaskOutcome.error.
    """
    jobs = list(jobs)
    if max_active < 1:
        raise ValueError("max_active must be positive")
    if len({id(j.ledger) for j in jobs}) != len(jobs):
        raise ValueError("concurrent tasks need independent ledgers")
    has_hessian = hessian is not None
    condition = threading.Condition()
    pending, active = {}, set()
    remaining = deque(range(len(jobs)))
    outcomes = [None] * len(jobs)
    rounds = []
    hessian_calls = []
    wait_s = [0.0] * len(jobs)
    controller_cpu_s = [0.0] * len(jobs)
    stopping = False

    def request(index, kind, atoms):
        future = Future()
        start = time.perf_counter()
        with condition:
            if stopping:
                raise RuntimeError(_STOPPED)
            pending[index] = (kind, atoms, future)
            condition.notify_all()
        try:
            return future.result()
        finally:
            wait_s[index] += time.perf_counter() - start

    def run(index):
        job = jobs[index]
        cpu_start = time.thread_time()
        try:

            def factory():
                return _ClientCalculator(
                    lambda kind, a: request(index, kind, a),
                    job.ledger,
                    has_hessian,
                    units,
                )

            outcomes[index] = TaskOutcome(result=job.run(factory))
        except Exception as exc:
            outcomes[index] = TaskOutcome(error=f"{type(exc).__name__}: {exc}")
        finally:
            controller_cpu_s[index] = time.thread_time() - cpu_start

    def controller_done(index, future):
        # Completion callback covers failures before run() reaches its try, and
        # failures in cleanup as well. Never leave a finished worker in active.
        with condition:
            if future.cancelled():
                outcomes[index] = TaskOutcome(error="controller cancelled")
            elif future.exception() is not None:
                exc = future.exception()
                outcomes[index] = TaskOutcome(error=f"{type(exc).__name__}: {exc}")
            active.discard(index)
            condition.notify_all()

    start = time.perf_counter()
    with ThreadPoolExecutor(max_workers=max_active) as executor:
        # Every request taken out of ``pending`` and not yet answered; released by
        # the finally clause if the coordinator is interrupted mid-round.
        taken = []
        try:
            while True:
                with condition:
                    while remaining and len(active) < max_active:
                        index = remaining.popleft()
                        active.add(index)
                        future = executor.submit(run, index)
                        future.add_done_callback(lambda f, i=index: controller_done(i, f))
                    if not active:
                        break
                    if not active.issubset(pending):
                        condition.wait()
                        continue
                    taken = [(i, *pending.pop(i)) for i in sorted(pending)]
                batch, indices, serial = [], [], []
                for index, kind, atoms, future in taken:
                    if kind == "hessian":
                        serial.append((index, atoms, future))
                        continue
                    try:
                        # Charge at physical dispatch, not when a controller
                        # queues work. Cancelled queued requests have no forward.
                        jobs[index].ledger.charge()
                    except BudgetExhausted as exc:
                        future.set_exception(exc)
                    else:
                        batch.append((atoms, future))
                        indices.append(index)
                if batch:
                    before = time.perf_counter()
                    outcome = "completed"
                    try:
                        values = list(evaluate_batch([atoms for atoms, _ in batch]))
                        if len(values) != len(batch):
                            raise ValueError(
                                "batch evaluator returned the wrong number of structures"
                            )
                    except Exception as exc:
                        finished = time.perf_counter()
                        outcome = f"error:{type(exc).__name__}"
                        for _, future in batch:
                            future.set_exception(exc)
                    else:
                        finished = time.perf_counter()
                        for (_, future), value in zip(batch, values):
                            future.set_result(value)
                    rounds.append(
                        {
                            "tasks": indices,
                            "size": len(batch),
                            "wall_s": finished - before,
                            "outcome": outcome,
                        }
                    )
                for index, atoms, future in serial:
                    try:
                        jobs[index].ledger.charge()
                    except BudgetExhausted as exc:
                        future.set_exception(exc)
                        continue
                    before = time.perf_counter()
                    outcome = "completed"
                    try:
                        value = np.asarray(hessian(atoms))
                    except Exception as exc:
                        outcome = f"error:{type(exc).__name__}"
                        future.set_exception(exc)
                    else:
                        future.set_result(value)
                    hessian_calls.append(
                        {
                            "task": index,
                            "wall_s": time.perf_counter() - before,
                            "outcome": outcome,
                        }
                    )
                taken = []
        finally:
            # Also release waiters when the coordinator is interrupted between
            # batches, so executor shutdown cannot wait forever on their futures.
            with condition:
                stopping = True
                for *_, future in [*taken, *pending.values()]:
                    if not future.done():
                        future.set_exception(RuntimeError(_STOPPED))
                condition.notify_all()
    wall = time.perf_counter() - start
    evaluator_s = sum(r["wall_s"] for r in rounds) + sum(h["wall_s"] for h in hessian_calls)
    return outcomes, {
        "wall_s": wall,
        "outside_evaluator_wall_s": wall - evaluator_s,
        "batch_rounds": rounds,
        "hessian_calls": hessian_calls,
        "client_wait_s": wait_s,
        "controller_cpu_s": controller_cpu_s,
    }


def scan_many(jobs, evaluate_batch, max_active, hessian=None, units=(1.0, 1.0)):
    """Run independent ``constrained_scan``s through ``run_many``.

    ``constrained_scan`` calls its factory once, so each scan has one client, as
    before. Timing keys are those of ``run_many``.
    """
    jobs = list(jobs)

    def task(job):
        def run(factory):
            return constrained_scan(
                job.atoms,
                job.intent,
                job.settings,
                factory,
                job.config,
                job.ledger,
                job.task,
                attempt=job.attempt,
                reserve=job.reserve,
                artifact=job.artifact,
            )

        return TaskJob(run, job.ledger)

    outcomes, timing = run_many(
        [task(job) for job in jobs], evaluate_batch, max_active, hessian, units
    )
    return [ScanOutcome(result=o.result, error=o.error) for o in outcomes], timing


class MACEBatchEvaluator:
    """Resident single-model, uncompiled MACE energy/force backend for scans.

    Builds each graph with the calculator's KeySpecification (including charge,
    spin and head), then combines disjoint graphs for ONE model call. Stress is
    kept enabled as in MACECalculator. This adapter supplies only energy/forces;
    certification and analytic Hessians continue through the original calculator.
    """

    def __init__(self, calculator):
        if len(calculator.models) != 1 or calculator.use_compile:
            raise ValueError("scan batching supports one uncompiled model")
        if calculator.model_type not in ("MACE", "PolarMACE", "EnergyDipoleMACE"):
            raise ValueError("scan batching requires an energy/force MACE model")
        self.calculator = calculator
        self.calls = []

    def __call__(self, atoms_list):
        import torch
        from mace import data
        from mace.tools import torch_geometric, torch_tools

        calc = self.calculator

        def sync():
            if str(calc.device).startswith("cuda"):
                torch.cuda.synchronize(calc.device)

        sync()
        start = time.perf_counter()
        keys = data.KeySpecification(
            info_keys=calc.info_keys,
            arrays_keys={**calc.arrays_keys, calc.charges_key: "charges"},
        )
        with torch_tools.default_dtype(calc.default_dtype):
            graphs = [
                data.AtomicData.from_config(
                    data.config_from_atoms(a, key_specification=keys, head_name=calc.head),
                    z_table=calc.z_table,
                    cutoff=calc.r_max,
                    heads=calc.available_heads,
                )
                for a in atoms_list
            ]
        batch = torch_geometric.Batch.from_data_list(graphs).to(calc.device)
        batch = calc._clone_batch(batch)
        dtype = next(calc.models[0].parameters()).dtype
        for key in batch.keys:
            if torch.is_tensor(batch[key]) and torch.is_floating_point(batch[key]):
                batch[key] = batch[key].to(dtype=dtype)
        sync()
        prepared = time.perf_counter()
        out = calc.models[0](
            batch.to_dict(),
            compute_stress=True,
            training=False,
            compute_edge_forces=calc.compute_atomic_stresses,
            compute_atomic_stresses=calc.compute_atomic_stresses,
        )
        sync()
        computed = time.perf_counter()
        energies = out["energy"].detach().cpu().numpy() * calc.energy_units_to_eV
        forces = out["forces"].detach().cpu().numpy() * (
            calc.energy_units_to_eV / calc.length_units_to_A
        )
        offsets = np.cumsum([0] + [len(a) for a in atoms_list])
        results = [
            (float(energies[i]), forces[offsets[i] : offsets[i + 1]].copy())
            for i in range(len(atoms_list))
        ]
        self.calls.append(
            {
                "size": len(atoms_list),
                "prepare_s": prepared - start,
                "model_s": computed - prepared,
                "return_s": time.perf_counter() - computed,
            }
        )
        return results
