"""Independent scan controllers sharing one resident, batched energy/force evaluator.

Each scan owns its calculator cache, guard, optimizer and ledger. Only the caller
thread touches the model. At each round it batches one pending structure from each
active controller; completed controllers are replaced from the input queue. No
timed batching window, shared ASE calculator or shared ledger segment is used.
"""

from __future__ import annotations

import threading
import time
from collections import deque
from concurrent.futures import Future, ThreadPoolExecutor
from dataclasses import dataclass

import numpy as np
from ase.calculators.calculator import Calculator, all_changes

from .ledger import BudgetExhausted
from .probe import constrained_scan


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


class _ClientCalculator(Calculator):
    implemented_properties = ["energy", "forces"]

    def __init__(self, request, ledger):
        super().__init__()
        self.request = request
        self.ledger = ledger
        self.evaluations = 0

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
        energy, forces = self.request(frame)
        self.results = {"energy": float(energy), "forces": np.asarray(forces).copy()}


def scan_many(jobs, evaluate_batch, max_active):
    """Run independent scans, returning input-ordered outcomes and round timings.

    ``evaluate_batch(frames)`` returns one (energy_eV, forces_eV_A) per frame;
    it must preserve each frame's electronic state. Supply a resident evaluator
    object to reuse its model across successive calls. Budgets are per job; a
    batch of m structures still pays m structure forwards, including failed calls.
    General backend/controller exceptions are retained in ScanOutcome.error.
    """
    jobs = list(jobs)
    if max_active < 1:
        raise ValueError("max_active must be positive")
    if len({id(j.ledger) for j in jobs}) != len(jobs):
        raise ValueError("concurrent scans need independent ledgers")
    condition = threading.Condition()
    pending, active = {}, set()
    remaining = deque(range(len(jobs)))
    outcomes = [None] * len(jobs)
    rounds = []
    wait_s = [0.0] * len(jobs)
    controller_cpu_s = [0.0] * len(jobs)
    stopping = False

    def request(index, atoms):
        future = Future()
        start = time.perf_counter()
        with condition:
            if stopping:
                raise RuntimeError("scan batch coordinator stopped")
            pending[index] = (atoms, future)
            condition.notify_all()
        try:
            return future.result()
        finally:
            wait_s[index] += time.perf_counter() - start

    def run(index):
        job = jobs[index]
        cpu_start = time.thread_time()
        try:
            calc = _ClientCalculator(lambda a: request(index, a), job.ledger)
            result = constrained_scan(
                job.atoms,
                job.intent,
                job.settings,
                lambda: calc,
                job.config,
                job.ledger,
                job.task,
                attempt=job.attempt,
                reserve=job.reserve,
                artifact=job.artifact,
            )
            outcomes[index] = ScanOutcome(result=result)
        except Exception as exc:
            outcomes[index] = ScanOutcome(error=f"{type(exc).__name__}: {exc}")
        finally:
            controller_cpu_s[index] = time.thread_time() - cpu_start

    def controller_done(index, future):
        # Completion callback covers failures before run() reaches its try, and
        # failures in cleanup as well. Never leave a finished worker in active.
        with condition:
            if future.cancelled():
                outcomes[index] = ScanOutcome(error="controller cancelled")
            elif future.exception() is not None:
                exc = future.exception()
                outcomes[index] = ScanOutcome(error=f"{type(exc).__name__}: {exc}")
            active.discard(index)
            condition.notify_all()

    start = time.perf_counter()
    with ThreadPoolExecutor(max_workers=max_active) as executor:
        batch = []
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
                    indices = sorted(pending)
                    batch = [pending.pop(i) for i in indices]
                selected, eligible = [], []
                for index, (atoms, future) in zip(indices, batch):
                    try:
                        # Charge at physical dispatch, not when a controller
                        # queues work. Cancelled queued requests have no forward.
                        jobs[index].ledger.charge()
                    except BudgetExhausted as exc:
                        future.set_exception(exc)
                    else:
                        selected.append((atoms, future))
                        eligible.append(index)
                batch, indices = selected, eligible
                if not batch:
                    continue
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
                batch = []
        finally:
            # Also release waiters when the coordinator is interrupted between
            # batches, so executor shutdown cannot wait forever on their futures.
            with condition:
                stopping = True
                for _, future in [*batch, *pending.values()]:
                    if not future.done():
                        future.set_exception(RuntimeError("scan batch coordinator stopped"))
                condition.notify_all()
    wall = time.perf_counter() - start
    return outcomes, {
        "wall_s": wall,
        "outside_evaluator_wall_s": wall - sum(r["wall_s"] for r in rounds),
        "batch_rounds": rounds,
        "client_wait_s": wait_s,
        "controller_cpu_s": controller_cpu_s,
    }


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
