import json

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes

from explore.ledger import BudgetExhausted, Ledger, counted_factory, instrument
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
