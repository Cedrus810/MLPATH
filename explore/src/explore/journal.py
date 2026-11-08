"""D2 run journal and the two v1 budget accounts (DISCOVERY_LAYER_DESIGN §6.2 v1 row;
EVALUATION_PROTOCOL §2.2-2.3).

The journal is the coordinator's commit order: dispatch (a cap is held), settle (the
real spend replaces the hold), reconcile (an interrupted attempt is charged), skip (a
task the account could not hold). Accounts are rebuilt from it on every start."""

from __future__ import annotations

import json
import os
import socket
import time
from pathlib import Path


def repair_tail(path):
    """Drop a last line torn by a kill (no trailing newline); keep everything before it."""
    path = Path(path)
    if not path.exists():
        return False
    data = path.read_bytes()
    if not data or data.endswith(b"\n"):
        return False
    path.write_bytes(data[: data.rfind(b"\n") + 1])
    return True


class Account:
    def __init__(self, name, budget):
        self.name = name
        self.budget = int(budget)
        self.spent = 0
        self.held = 0

    @property
    def available(self):
        return self.budget - self.spent - self.held

    def hold(self, cap):
        if cap > self.available:
            return False
        self.held += cap
        return True

    def settle(self, cap, spent):
        self.held -= cap
        self.spent += spent

    def charge(self, amount):
        self.spent += amount


class Journal:
    def __init__(self, path):
        self.path = Path(path)
        repair_tail(self.path)
        self.events = []
        if self.path.exists():
            for line in self.path.read_text(encoding="utf-8").splitlines():
                if line.strip():
                    self.events.append(json.loads(line))

    def append(self, event):
        event = {**event, "host": socket.gethostname(), "t": time.time()}
        line = json.dumps(event) + "\n"
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with self.path.open("a", encoding="utf-8") as handle:
            handle.write(line)
            handle.flush()
            os.fsync(handle.fileno())
        self.events.append(event)
