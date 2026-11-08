import json
import threading

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.calculators.lj import LennardJones

from explore.ledger import (
    BudgetExhausted,
    Ledger,
    counted_factory,
    instrument,
    tally,
    view_factory,
)
from prrs.calculators import double_well_factory


def lines(path):
    return [json.loads(x) for x in path.read_text(encoding="utf-8").splitlines() if x.strip()]


def test_segment_writes_open_and_close_with_counts(tmp_path):
    ledger = Ledger(tmp_path / "ledger.jsonl", arm="D2", budget=10)
    with ledger.segment(task="t1", stage="probe") as seg:
        ledger.note_request()
        ledger.charge()
        ledger.note_request()
        seg.outcome = "completed"
    events = lines(tmp_path / "ledger.jsonl")
    assert [e["event"] for e in events] == ["open", "close"]
    close = events[1]
    assert close["forwards"] == 1 and close["requests"] == 2
    assert (
        close["outcome"] == "completed" and close["arm"] == "D2" and close["stage"] == "probe"
    )
    assert ledger.spent == 1 and ledger.remaining == 9


def test_charge_refused_before_forward_when_budget_is_spent(tmp_path):
    ledger = Ledger(tmp_path / "l.jsonl", arm="D2", budget=1)
    with ledger.segment(task="t", stage="probe"):
        ledger.charge()
        with pytest.raises(BudgetExhausted):
            ledger.charge()
    assert ledger.spent == 1


def test_reservation_is_held_back_and_unused_part_returned(tmp_path):
    ledger = Ledger(tmp_path / "l.jsonl", arm="D2", budget=10)
    with ledger.segment(task="a", stage="probe", reserve=6):
        assert ledger.remaining == 4
        ledger.charge()
        ledger.charge()
        assert ledger.remaining == 4  # drawn from the reservation, not the pool
    assert ledger.spent == 2 and ledger.remaining == 8


def test_reservation_larger_than_remaining_is_refused_with_no_spend(tmp_path):
    ledger = Ledger(tmp_path / "l.jsonl", arm="D2", budget=5)
    with pytest.raises(BudgetExhausted):
        with ledger.segment(task="a", stage="probe", reserve=6):
            pass
    assert ledger.spent == 0
    assert not (tmp_path / "l.jsonl").exists()


def test_exception_inside_segment_is_recorded_and_reraised(tmp_path):
    ledger = Ledger(tmp_path / "l.jsonl", arm="D2", budget=5)
    with pytest.raises(ZeroDivisionError):
        with ledger.segment(task="a", stage="probe"):
            ledger.charge()
            1 / 0
    close = lines(tmp_path / "l.jsonl")[-1]
    assert close["outcome"] == "error:ZeroDivisionError" and close["forwards"] == 1


def test_reopen_after_crash_charges_open_segment_reservation(tmp_path):
    path = tmp_path / "l.jsonl"
    ledger = Ledger(path, arm="D2", budget=100)
    with ledger.segment(task="done", stage="probe"):
        ledger.charge()
    # simulate a crash: an open line with no close
    with path.open("a", encoding="utf-8") as f:
        f.write(
            json.dumps(
                {
                    "event": "open",
                    "arm": "D2",
                    "task": "crashed",
                    "attempt": 0,
                    "stage": "probe",
                    "reserve": 7,
                }
            )
            + "\n"
        )
    again = Ledger(path, arm="D2", budget=100)
    assert again.spent == 1 + 7
    assert again.unreconciled == [
        {"task": "crashed", "attempt": 0, "stage": "probe", "reserve": 7}
    ]


def test_forwards_outside_segments_are_kept_as_unattributed(tmp_path):
    path = tmp_path / "l.jsonl"
    ledger = Ledger(path, arm="D0", budget=None)
    ledger.charge()
    ledger.close()
    close = lines(path)[-1]
    assert close["stage"] == "unattributed" and close["forwards"] == 1


def test_segments_do_not_nest(tmp_path):
    ledger = Ledger(tmp_path / "l.jsonl", arm="D2")
    with ledger.segment(task="a", stage="probe"):
        with pytest.raises(RuntimeError):
            with ledger.segment(task="b", stage="quench"):
                pass


def test_retry_with_same_key_does_not_refund_crashed_reservation(tmp_path):
    path = tmp_path / "l.jsonl"
    ledger = Ledger(path, arm="D2", budget=100)
    with ledger.segment(task="t", stage="probe"):
        ledger.charge()
    # simulate a crash: an open line for the same key with no close
    with path.open("a", encoding="utf-8") as f:
        f.write(
            json.dumps(
                {
                    "event": "open",
                    "arm": "D2",
                    "task": "t",
                    "attempt": 0,
                    "stage": "probe",
                    "reserve": 7,
                }
            )
            + "\n"
        )
    resumed = Ledger(path, arm="D2", budget=100)
    with resumed.segment(task="t", stage="probe"):
        resumed.charge()
    again = Ledger(path, arm="D2", budget=100)
    assert again.spent == 1 + 7 + 1
    assert again.unreconciled == [{"task": "t", "attempt": 0, "stage": "probe", "reserve": 7}]


def dimer(r):
    return Atoms("C2", positions=[[0, 0, 0], [r, 0, 0]])


def test_cache_hit_is_a_request_not_a_forward(tmp_path):
    ledger = Ledger(tmp_path / "l.jsonl", arm="D2", budget=10)
    atoms = dimer(1.3)
    atoms.calc = instrument(double_well_factory(), ledger)
    with ledger.segment(task="t", stage="probe"):
        atoms.get_potential_energy()
        atoms.get_forces()  # same geometry: ASE cache
        atoms.positions[1, 0] = 1.4
        atoms.get_forces()  # new geometry: a real forward
    close = lines(tmp_path / "l.jsonl")[-1]
    assert close["forwards"] == 2 and close["requests"] == 3


class Boom(Calculator):
    implemented_properties = ["energy", "forces"]

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        raise RuntimeError("model blew up")


def test_failed_forward_is_still_charged(tmp_path):
    ledger = Ledger(tmp_path / "l.jsonl", arm="D2", budget=10)
    atoms = dimer(1.3)
    atoms.calc = instrument(Boom(), ledger)
    with pytest.raises(RuntimeError):
        with ledger.segment(task="t", stage="probe"):
            atoms.get_potential_energy()
    assert ledger.spent == 1


def test_budget_refusal_happens_before_the_model_is_called(tmp_path):
    calls = []

    class Spy(Calculator):
        implemented_properties = ["energy", "forces"]

        def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
            calls.append(1)
            self.results = {"energy": 0.0, "forces": np.zeros((len(atoms), 3))}

    ledger = Ledger(tmp_path / "l.jsonl", arm="D2", budget=0)
    atoms = dimer(1.3)
    atoms.calc = instrument(Spy(), ledger)
    with ledger.segment(task="t", stage="probe"):
        with pytest.raises(BudgetExhausted):
            atoms.get_potential_energy()
    assert calls == []


def test_double_instrumentation_is_an_error(tmp_path):
    ledger = Ledger(tmp_path / "l.jsonl", arm="D2")
    calc = instrument(double_well_factory(), ledger)
    with pytest.raises(ValueError):
        instrument(calc, ledger)


def test_counted_factory_counts_through_prrs_guard(tmp_path):
    from prrs.config import SearchConfig
    from prrs.reliability import GuardedCalculator

    ledger = Ledger(tmp_path / "l.jsonl", arm="D2", budget=10)
    factory = counted_factory(double_well_factory, ledger)
    atoms = dimer(1.3)
    atoms.calc = GuardedCalculator(factory(), SearchConfig(closed_shell_only=False))
    with ledger.segment(task="t", stage="probe"):
        atoms.get_potential_energy()
        atoms.get_forces()
    assert ledger.spent == 1


def test_analytic_hessian_is_charged(tmp_path):
    class WithHessian(Calculator):
        implemented_properties = ["energy", "forces"]

        def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
            self.results = {"energy": 0.0, "forces": np.zeros((len(atoms), 3))}

        def get_hessian(self, atoms=None):
            return np.zeros((3 * len(atoms), 3 * len(atoms)))

    ledger = Ledger(tmp_path / "l.jsonl", arm="D2", budget=10)
    calc = instrument(WithHessian(), ledger)
    with ledger.segment(task="t", stage="hessian"):
        calc.get_hessian(dimer(1.3))
    assert ledger.spent == 1


def test_cache_hits_are_counted_per_segment(tmp_path):
    ledger = Ledger(tmp_path / "l.jsonl", arm="t", budget=100)
    atoms = dimer(1.0)
    atoms.calc = instrument(LennardJones(), ledger)
    with ledger.segment(task="a", stage="probe"):
        atoms.get_potential_energy()  # forward
        atoms.get_potential_energy()  # cache hit
        atoms.get_forces()  # cache hit (LJ computes both)
    ledger.close()
    close = [e for e in lines(tmp_path / "l.jsonl") if e["event"] == "close"][0]
    assert close["forwards"] == 1 and close["cache_hits"] == 2


def test_reused_factory_builds_once_and_still_counts_every_forward(tmp_path):
    from explore.ledger import reused_counted_factory

    builds = {"n": 0}

    def factory():
        builds["n"] += 1
        return LennardJones()

    ledger = Ledger(tmp_path / "l.jsonl", arm="t", budget=100)
    make = reused_counted_factory(factory, ledger)
    first, second = make(), make()
    assert first is second and builds["n"] == 1
    a, b = dimer(1.0), dimer(1.2)
    a.calc, b.calc = first, second
    a.get_potential_energy()
    b.get_potential_energy()
    a.get_potential_energy()  # positions differ from b: a real forward again
    assert ledger.spent == 3


def test_reused_factory_drops_cached_results_on_each_call(tmp_path):
    # Two callers with the same positions but different info (e.g. charge) must not
    # share a cached result: ASE's state check ignores atoms.info.
    from explore.ledger import reused_counted_factory

    ledger = Ledger(tmp_path / "l.jsonl", arm="t", budget=100)
    make = reused_counted_factory(LennardJones, ledger)
    a = dimer(1.0)
    a.calc = make()
    a.get_potential_energy()
    b = dimer(1.0)
    b.info["charge"] = -1
    b.calc = make()
    b.get_potential_energy()
    assert ledger.spent == 2


class ChargeAware(Calculator):
    implemented_properties = ["energy", "forces"]

    def __init__(self):
        super().__init__()
        self.calls = 0

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        self.calls += 1
        charge = float(self.atoms.info.get("charge", 0))
        self.results = {
            "energy": charge + float(self.atoms.positions.sum()),
            "forces": np.zeros((len(self.atoms), 3)),
        }


def test_replay_false_starts_from_zero_on_a_shared_file(tmp_path):
    path = tmp_path / "l.jsonl"
    first = Ledger(path, arm="D2", budget=5, replay=False)
    with first.segment(task="a", stage="probe"):
        first.charge()
        first.charge()
    second = Ledger(path, arm="D2", budget=5, replay=False)
    assert second.spent == 0 and second.remaining == 5
    assert Ledger(path, arm="D2").spent == 2


def test_concurrent_task_ledgers_share_one_file(tmp_path):
    path = tmp_path / "l.jsonl"

    def work(k):
        ledger = Ledger(path, arm="D2", budget=1000, replay=False)
        for _ in range(50):
            with ledger.segment(task=f"t{k}", stage="probe"):
                ledger.charge()
        ledger.close()

    threads = [threading.Thread(target=work, args=(k,)) for k in range(8)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    rows = tally(path)
    assert {key: row["forwards"] for key, row in rows.items()} == {
        (f"t{k}", 0): 50 for k in range(8)
    }
    assert all(row["open"] == 0 for row in rows.values())


def test_tally_reports_open_segments_and_ignores_a_torn_last_line(tmp_path):
    path = tmp_path / "l.jsonl"
    ledger = Ledger(path, arm="D2", replay=False)
    with ledger.segment(task="done", stage="probe"):
        ledger.charge()
    ledger._write(
        {"event": "open", "task": "cut", "attempt": 1, "stage": "probe", "reserve": None}
    )
    with path.open("a", encoding="utf-8") as fh:
        fh.write('{"event": "close", "task": "cut"')
    rows = tally(path)
    assert rows[("done", 0)] == {"forwards": 1, "open": 0}
    assert rows[("cut", 1)] == {"forwards": 0, "open": 1}


def test_counted_view_charges_real_forwards_only(tmp_path):
    inner = ChargeAware()
    ledger = Ledger(tmp_path / "l.jsonl", arm="D2", budget=10, replay=False)
    atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.74]])
    atoms.calc = view_factory(inner, ledger)()
    with ledger.segment(task="t", stage="probe") as seg:
        atoms.get_potential_energy()
        atoms.get_forces()
        atoms.positions[1, 2] += 0.1
        atoms.get_potential_energy()
    assert inner.calls == 2 and ledger.spent == 2 and seg.forwards == 2
    assert seg.cache_hits >= 1


def test_views_never_share_a_result_across_electronic_states(tmp_path):
    inner = ChargeAware()
    ledger = Ledger(tmp_path / "l.jsonl", arm="D2", replay=False)
    make = view_factory(inner, ledger)
    a = Atoms("H", positions=[[0, 0, 0]])
    a.info["charge"] = -1
    a.calc = make()
    b = a.copy()
    b.info["charge"] = 0
    b.calc = make()
    assert a.get_potential_energy() == -1.0
    assert b.get_potential_energy() == 0.0
    assert ledger.spent == 2


def test_counted_view_refuses_before_the_model_runs(tmp_path):
    inner = ChargeAware()
    ledger = Ledger(tmp_path / "l.jsonl", arm="D2", budget=1, replay=False)
    atoms = Atoms("H", positions=[[0, 0, 0]])
    atoms.calc = view_factory(inner, ledger)()
    atoms.get_potential_energy()
    atoms.positions += 0.1
    with pytest.raises(BudgetExhausted):
        atoms.get_potential_energy()
    assert inner.calls == 1


def test_view_factory_refuses_an_already_counted_calculator(tmp_path):
    ledger = Ledger(tmp_path / "l.jsonl", arm="D2", replay=False)
    inner = instrument(ChargeAware(), ledger)
    with pytest.raises(ValueError):
        view_factory(inner, ledger)
