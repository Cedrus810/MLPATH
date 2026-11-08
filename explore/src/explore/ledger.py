"""Forward ledger for the evaluation protocol (EVALUATION_PROTOCOL_DESIGN_2026-10-13 §2).

One unit is one real model forward on one structure. ASE's property cache sits in
front of ``calculate``: a cache hit is a request, not a forward. A forward that was
started and then failed is still paid for. Reopening a ledger never refunds anything:
a segment that was opened and never closed is charged its reservation as an upper
bound and listed in ``unreconciled``.
"""

from __future__ import annotations

import contextlib
import json
import os
import socket
import threading
import time
from pathlib import Path

from ase.calculators.calculator import Calculator, all_changes

_WRITE_LOCK = threading.Lock()


class BudgetExhausted(RuntimeError):
    """Raised before a forward (or a reservation) the remaining budget cannot cover."""


class Segment:
    """Per-segment tallies, written in the close event.

    ``cache_hits`` = physical-calculator property re-reads that triggered no new forward
    (a counted ``get_property`` during which the ledger's spend did not change). Under
    ``prrs.reliability.GuardedCalculator`` this is about 1 per forward (the guard reads
    energy then forces from the physical calculator, so the second read is a hit when the
    physical calculator produced both properties in one calculate);
    it is NOT a measure of structure-level cache savings.
    """

    def __init__(self, task, attempt, stage, reserve):
        self.task, self.attempt, self.stage = task, attempt, stage
        self.reserve = reserve
        self.reserve_left = reserve
        self.forwards = 0
        self.requests = 0
        self.cache_hits = 0
        self.outcome = None
        self.artifact = None


class Ledger:
    def __init__(self, path, arm, budget=None, replay=True):
        self.path = Path(path)
        self.arm = arm
        self.budget = budget
        self.spent = 0
        self.reserved = 0
        self.unreconciled = []
        self.host = socket.gethostname()
        self._segment = None
        self._unattributed = 0
        if replay and self.path.exists():
            self._replay()

    @property
    def remaining(self):
        if self.budget is None:
            return None
        return self.budget - self.spent - self.reserved

    def _replay(self):
        pending = {}
        for line in self.path.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            event = json.loads(line)
            key = (event["task"], event["attempt"], event["stage"])
            if event["event"] == "open":
                pending.setdefault(key, []).append(event)
            elif event["event"] == "close":
                stack = pending.get(key)
                if stack:
                    stack.pop()
                self.spent += event["forwards"]
        for (task, attempt, stage), stack in pending.items():
            for event in stack:
                self.spent += event["reserve"] or 0
                self.unreconciled.append(
                    {
                        "task": task,
                        "attempt": attempt,
                        "stage": stage,
                        "reserve": event["reserve"],
                    }
                )

    def _write(self, event):
        record = {"arm": self.arm, "host": self.host, "t": time.time(), **event}
        line = json.dumps(record) + "\n"
        with _WRITE_LOCK:
            self.path.parent.mkdir(parents=True, exist_ok=True)
            with self.path.open("a", encoding="utf-8") as handle:
                handle.write(line)
                handle.flush()
                os.fsync(handle.fileno())

    @contextlib.contextmanager
    def segment(self, task, stage, attempt=0, reserve=None):
        if self._segment is not None:
            raise RuntimeError("ledger segments do not nest")
        if reserve is not None and self.budget is not None and reserve > self.remaining:
            raise BudgetExhausted(f"reserve {reserve} exceeds remaining {self.remaining}")
        seg = Segment(task, attempt, stage, reserve)
        ident = {"task": task, "attempt": attempt, "stage": stage}
        self._write({"event": "open", **ident, "reserve": reserve})
        if reserve is not None:
            self.reserved += reserve
        self._segment = seg
        started = time.monotonic()
        try:
            yield seg
        except BaseException as exc:
            if seg.outcome is None:
                seg.outcome = f"error:{type(exc).__name__}"
            raise
        finally:
            self._segment = None
            if seg.reserve_left is not None:
                self.reserved -= seg.reserve_left
            self._write(
                {
                    "event": "close",
                    **ident,
                    "forwards": seg.forwards,
                    "requests": seg.requests,
                    "cache_hits": seg.cache_hits,
                    "wall_s": time.monotonic() - started,
                    "outcome": seg.outcome,
                    "artifact": seg.artifact,
                }
            )

    def charge(self):
        seg = self._segment
        if seg is not None and seg.reserve_left:
            seg.reserve_left -= 1
            self.reserved -= 1
        elif self.budget is not None and self.remaining < 1:
            raise BudgetExhausted(f"no budget left for a forward (spent {self.spent})")
        self.spent += 1
        if seg is None:
            self._unattributed += 1
        else:
            seg.forwards += 1

    def note_request(self):
        if self._segment is not None:
            self._segment.requests += 1

    def note_cache_hit(self):
        if self._segment is not None:
            self._segment.cache_hits += 1

    def close(self):
        if self._unattributed:
            ident = {"task": "<unattributed>", "attempt": 0, "stage": "unattributed"}
            self._write(
                {
                    "event": "close",
                    **ident,
                    "forwards": self._unattributed,
                    "requests": 0,
                    "cache_hits": 0,
                    "wall_s": None,
                    "outcome": None,
                    "artifact": None,
                }
            )
            self._unattributed = 0


def instrument(calc, ledger):
    """Route a physical calculator's real evaluations through the ledger.

    ``calculate`` is where ASE does the work and is charged *before* the call, so a
    refused forward never reaches the model and a failed one is still paid for.
    ``get_property`` is counted as a request; if it triggers no ``calculate`` it is also
    a cache hit, i.e. a physical-calculator property re-read with no new forward. Under
    ``GuardedCalculator`` that is ~1 per forward and does not measure structure-level
    ASE cache reuse at the guard (that never reaches this counted calculator). A backend's own
    ``get_hessian`` bypasses ``calculate`` and is charged as one forward.
    """
    if getattr(calc, "_explore_ledger", None) is not None:
        raise ValueError("calculator is already counted by a ledger")
    inner_calculate = calc.calculate
    inner_get_property = calc.get_property

    def calculate(*args, **kwargs):
        ledger.charge()
        return inner_calculate(*args, **kwargs)

    def get_property(name, atoms=None, allow_calculation=True):
        ledger.note_request()
        before = ledger.spent
        value = inner_get_property(name, atoms, allow_calculation)
        if ledger.spent == before:
            ledger.note_cache_hit()
        return value

    calc.calculate = calculate
    calc.get_property = get_property
    if hasattr(calc, "get_hessian"):
        inner_hessian = calc.get_hessian

        def get_hessian(*args, **kwargs):
            ledger.charge()
            return inner_hessian(*args, **kwargs)

        calc.get_hessian = get_hessian
    calc._explore_ledger = ledger
    return calc


def counted_factory(factory, ledger):
    """A factory whose every calculator is counted by ``ledger``."""

    def make():
        return instrument(factory(), ledger)

    return make


def reused_counted_factory(factory, ledger):
    """Like ``counted_factory`` but the physical calculator is built once and reused.

    Measured on b04 (DONE_LOG 2026-10-22): one C2 certification calls the factory 13
    times; rebuilding POLAR-1 each time cost 8.2 s of 39.5 s. Every forward is still
    charged. Each call drops the calculator's cached results, because ASE's state check
    ignores ``atoms.info`` (charge, spin): two callers with equal positions but different
    electronic state must not share a result. Within one caller the cache works as usual.
    """
    calc = None

    def make():
        nonlocal calc
        if calc is None:
            calc = instrument(factory(), ledger)
        calc.reset()
        return calc

    return make


def tally(path):
    """Closed forwards and still-open segments per (task, attempt) in a ledger file.

    For a file shared by many per-task ``Ledger(replay=False)`` instances (discovery D2).
    ``open > 0`` means the attempt was interrupted inside a segment, so ``forwards`` is
    only a lower bound. A last line without its newline was torn by a kill and is
    ignored: it was never acknowledged."""
    path = Path(path)
    rows = {}
    if not path.exists():
        return rows
    text = path.read_text(encoding="utf-8")
    lines = text.split("\n")
    if not text.endswith("\n"):
        lines = lines[:-1]
    for line in lines:
        if not line.strip():
            continue
        event = json.loads(line)
        row = rows.setdefault((event["task"], event["attempt"]), {"forwards": 0, "open": 0})
        if event["event"] == "open":
            row["open"] += 1
        elif event["event"] == "close":
            row["forwards"] += event["forwards"]
            if event["stage"] != "unattributed":
                row["open"] -= 1
    return rows


class CountedView(Calculator):
    """A per-caller counted window onto one resident physical calculator.

    ``reused_counted_factory`` binds the resident model to one ledger; D2 needs one
    ledger per task attempt (``scan_many`` and per-task caps), so the model is shared
    and only the counting is per caller. Each view has its own ASE cache, and a fresh
    view is made per factory call, so two callers that differ only in ``atoms.info``
    (charge, spin) never share a result. The forward is charged before the model runs."""

    def __init__(self, inner, ledger):
        super().__init__()
        self.inner = inner
        self.ledger = ledger
        self.implemented_properties = list(inner.implemented_properties)
        if hasattr(inner, "get_hessian"):
            self.get_hessian = self._get_hessian

    def get_property(self, name, atoms=None, allow_calculation=True):
        self.ledger.note_request()
        before = self.ledger.spent
        value = super().get_property(name, atoms, allow_calculation)
        if self.ledger.spent == before:
            self.ledger.note_cache_hit()
        return value

    def calculate(
        self, atoms=None, properties=("energy", "forces"), system_changes=all_changes
    ):
        self.ledger.charge()
        super().calculate(atoms, properties, system_changes)
        self.inner.reset()
        self.inner.calculate(self.atoms, properties, system_changes)
        self.results = dict(self.inner.results)

    def _get_hessian(self, *args, **kwargs):
        self.ledger.charge()
        return self.inner.get_hessian(*args, **kwargs)


def view_factory(inner, ledger):
    """Factory of ``CountedView``s on ``inner`` charged to ``ledger``."""
    if getattr(inner, "_explore_ledger", None) is not None:
        raise ValueError("calculator is already counted by a ledger")

    def make():
        return CountedView(inner, ledger)

    return make
