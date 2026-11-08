import threading

import numpy as np

from explore.ledger import Ledger
from explore.parallel import TaskJob, run_many
from test_scan_performance import CoupledBonds, triangle


def evaluator(calls):
    def batch(frames):
        calls.append(len(frames))
        out = []
        for f in frames:
            a = f.copy()
            a.calc = CoupledBonds()
            out.append((a.get_potential_energy(), a.get_forces()))
        return out

    return batch


def energy_job(ledger, positions):
    def run(factory):
        a = triangle()
        a.positions = positions
        a.calc = factory()
        return a.get_potential_energy()

    return TaskJob(run, ledger)


def test_tasks_are_batched_and_each_ledger_pays_its_own(tmp_path):
    calls = []
    jobs = [
        energy_job(
            Ledger(tmp_path / "l.jsonl", "t", replay=False), triangle().positions + 0.01 * k
        )
        for k in range(4)
    ]
    outcomes, timing = run_many(jobs, evaluator(calls), max_active=4)
    assert all(o.error is None for o in outcomes)
    assert [j.ledger.spent for j in jobs] == [1, 1, 1, 1]
    assert calls == [4]


def test_hessian_requests_run_alone_on_the_coordinator_and_cost_one(tmp_path):
    seen = []

    def hessian(atoms):
        seen.append(threading.current_thread() is threading.main_thread())
        return np.eye(3 * len(atoms))

    def run(factory):
        a = triangle()
        a.calc = factory()
        a.get_potential_energy()
        return factory().get_hessian(a)

    ledger = Ledger(tmp_path / "l.jsonl", "t", replay=False)
    (out,), _ = run_many([TaskJob(run, ledger)], evaluator([]), 1, hessian=hessian)
    assert out.error is None and out.result.shape == (9, 9)
    assert seen == [True] and ledger.spent == 2


def test_no_hessian_backend_means_no_get_hessian_attribute(tmp_path):
    def run(factory):
        return hasattr(factory(), "get_hessian")

    ledger = Ledger(tmp_path / "l.jsonl", "t", replay=False)
    (out,), _ = run_many([TaskJob(run, ledger)], evaluator([]), 1)
    assert out.result is False


def test_one_failing_task_does_not_stop_the_others(tmp_path):
    def bad(factory):
        raise ValueError("boom")

    good = energy_job(Ledger(tmp_path / "g.jsonl", "t", replay=False), triangle().positions)
    outcomes, _ = run_many(
        [TaskJob(bad, Ledger(tmp_path / "b.jsonl", "t", replay=False)), good],
        evaluator([]),
        2,
    )
    assert outcomes[0].error.startswith("ValueError") and outcomes[1].error is None
