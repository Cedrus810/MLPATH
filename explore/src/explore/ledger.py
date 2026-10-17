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
import time
from pathlib import Path


class BudgetExhausted(RuntimeError):
    """Raised before a forward (or a reservation) the remaining budget cannot cover."""


class Segment:
    def __init__(self, task, attempt, stage, reserve):
        self.task, self.attempt, self.stage = task, attempt, stage
        self.reserve = reserve
        self.reserve_left = reserve
        self.forwards = 0
        self.requests = 0
        self.outcome = None
        self.artifact = None


class Ledger:
    def __init__(self, path, arm, budget=None):
        self.path = Path(path)
        self.arm = arm
        self.budget = budget
        self.spent = 0
        self.reserved = 0
        self.unreconciled = []
        self.host = socket.gethostname()
        self._segment = None
        self._unattributed = 0
        if self.path.exists():
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
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with self.path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(record) + "\n")
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

    def close(self):
        if self._unattributed:
            ident = {"task": "<unattributed>", "attempt": 0, "stage": "unattributed"}
            self._write(
                {
                    "event": "close",
                    **ident,
                    "forwards": self._unattributed,
                    "requests": 0,
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
    ``get_property`` is counted as a request (cache hits included). A backend's own
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
        return inner_get_property(name, atoms, allow_calculation)

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
