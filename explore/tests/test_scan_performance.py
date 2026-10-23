import threading
import cProfile

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes

from explore.ledger import Ledger, counted_factory
from explore.parallel import ScanJob, scan_many
import explore.parallel as parallel
from explore.probe import Intent, ScanSettings, constrained_scan
from prrs.config import SearchConfig


class CoupledBonds(Calculator):
    """Invariant triangle with a known constrained minimum at every r01."""

    implemented_properties = ["energy", "forces"]

    def calculate(
        self, atoms=None, properties=("energy", "forces"), system_changes=all_changes
    ):
        super().calculate(atoms, properties, system_changes)
        p = atoms.positions
        q, r, s = (np.linalg.norm(p[j] - p[i]) for i, j in ((0, 1), (0, 2), (1, 2)))
        a = r - (1.7 + 0.2 * (q - 1.5))
        b = s - 1.6
        energy = a * a + 0.1 * b * b + 0.4 * (q - 1.5) ** 2
        forces = np.zeros_like(p)
        for (i, j), d in zip(
            ((0, 1), (0, 2), (1, 2)), (0.8 * (q - 1.5) - 0.4 * a, 2 * a, 0.2 * b)
        ):
            f = d * (p[j] - p[i]) / np.linalg.norm(p[j] - p[i])
            forces[i] += f
            forces[j] -= f
        self.results = {"energy": energy, "forces": forces}


def triangle():
    q, r, s = 1.5, 1.7, 1.6
    x = (r * r + q * q - s * s) / (2 * q)
    return Atoms("Ar3", positions=((0, 0, 0), (q, 0, 0), (x, np.sqrt(r * r - x * x), 0)))


CONFIG = SearchConfig(closed_shell_only=False)
INTENT = Intent("q", ((0, 1),), (1.0,))


def settings(reuse=True, n=7):
    return ScanSettings(
        tuple(0.05 * i for i in range(n)), 1e-6, 150, 3.0, reuse_curvature=reuse
    )


def check_oracle(result):
    assert result.status == "completed"
    for p in result.points:
        assert p.converged
        assert abs(p.q_actual - p.q_target) < 1e-10
        assert abs(p.energy_eV - 0.4 * p.delta_A**2) < 1e-9
        assert abs(p.reaction_force - 0.8 * p.delta_A) < 3e-6
        assert p.fmax_projected < 1e-6


def test_curvature_transfer_reduces_real_calls_and_preserves_constrained_solution(tmp_path):
    ledgers, results = [], []
    for reuse in (False, True):
        ledger = Ledger(tmp_path / f"{reuse}.jsonl", "test", budget=2000)
        result = constrained_scan(
            triangle(),
            INTENT,
            settings(reuse),
            counted_factory(CoupledBonds, ledger),
            CONFIG,
            ledger,
            str(reuse),
        )
        check_oracle(result)
        ledgers.append(ledger)
        results.append(result)
    assert ledgers[1].spent < ledgers[0].spent
    # Initial minimum is a zero-step point: no learned Hessian exists yet.
    assert results[1].points[0].steps == 0


def test_batched_controllers_keep_solutions_budgets_and_caches_independent(tmp_path):
    jobs = []
    for i in range(5):
        jobs.append(
            ScanJob(
                triangle(),
                INTENT,
                settings(n=3 + i),
                CONFIG,
                Ledger(tmp_path / f"job{i}.jsonl", "test", budget=2000),
                str(i),
            )
        )
    jobs[1].ledger.budget = 1
    batches = []
    owner = threading.get_ident()

    def evaluate(frames):
        assert threading.get_ident() == owner  # the model has one owner thread
        batches.append(len(frames))
        values = []
        for atoms in frames:
            atoms.calc = CoupledBonds()
            values.append((atoms.get_potential_energy(), atoms.get_forces()))
        return values

    outcomes, metrics = scan_many(jobs, evaluate, max_active=3)
    assert outcomes[1].result.status == "budget_exhausted"
    assert jobs[1].ledger.spent == 1
    for i, outcome in enumerate(outcomes):
        assert outcome.error is None
        if i != 1:
            check_oracle(outcome.result)
    assert sum(j.ledger.spent for j in jobs) == sum(batches)
    assert max(batches) == 3
    assert len(metrics["batch_rounds"]) == len(batches)


def test_backend_failure_releases_all_waiters_and_keeps_paid_cost(tmp_path):
    jobs = [
        ScanJob(
            triangle(),
            INTENT,
            settings(),
            CONFIG,
            Ledger(tmp_path / f"{i}.jsonl", "test", budget=100),
            str(i),
        )
        for i in range(3)
    ]

    def fail(frames):
        raise RuntimeError("failed batch")

    outcomes, _ = scan_many(jobs, fail, max_active=2)
    assert all(o.result is None and o.error == "RuntimeError: failed batch" for o in outcomes)
    assert all(j.ledger.spent == 1 for j in jobs)
    assert all(j.ledger._segment is None for j in jobs)


def test_worker_startup_exception_never_leaves_coordinator_waiting(tmp_path, monkeypatch):
    owner = threading.get_ident()
    original = parallel.time.thread_time

    def broken_clock():
        if threading.get_ident() != owner:
            raise RuntimeError("controller startup failed")
        return original()

    monkeypatch.setattr(parallel.time, "thread_time", broken_clock)
    jobs = [
        ScanJob(
            triangle(),
            INTENT,
            settings(),
            CONFIG,
            Ledger(tmp_path / f"{i}.jsonl", "test", budget=100),
            str(i),
        )
        for i in range(3)
    ]
    outcomes, _ = scan_many(jobs, lambda frames: pytest.fail("no model call expected"), 2)
    assert all(o.error == "RuntimeError: controller startup failed" for o in outcomes)
    assert all(j.ledger.spent == 0 for j in jobs)


def test_coordinator_interruption_releases_waiters_with_parent_profiler_enabled(tmp_path):
    class StopProbe(BaseException):
        pass

    jobs = [
        ScanJob(
            triangle(),
            INTENT,
            settings(),
            CONFIG,
            Ledger(tmp_path / f"{i}.jsonl", "test", budget=100),
            str(i),
        )
        for i in range(2)
    ]

    def stop(frames):
        raise StopProbe()

    profiler = cProfile.Profile()
    profiler.enable()
    try:
        with pytest.raises(StopProbe):
            scan_many(jobs, stop, 2)
    finally:
        profiler.disable()
    assert all(j.ledger._segment is None and j.ledger.spent == 1 for j in jobs)
