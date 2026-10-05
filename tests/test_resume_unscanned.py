"""Resume keeps the reasons no scan record carries (search._rebuild_for_resume [3]).

A direction selected and then never offered -- the loop stopped at a budget before its
batch -- has no scan record, and its reason lived only in memory: every resumed run
reported it as "unknown". The double-well resume tests never reach such a budget; the
relay testbed's reactant does (conformer_trials_per_microstate stops its last
direction), which is where this was found, in discover mode and in connect mode alike.
"""

import pytest
from ase.optimize import LBFGS

from prrs.calculators import ProtonRelay, proton_relay_atoms, proton_relay_factory
from prrs.config import SearchConfig
from prrs.search import ReactionSearch


def _product():
    atoms = proton_relay_atoms((1, 1, 1))
    atoms.calc = ProtonRelay()
    LBFGS(atoms, logfile=None).run(fmax=1e-4, steps=2000)
    atoms.calc = None
    return atoms


CHILD = """
import os
import sys

sys.path.insert(0, {root!r} + "/src")

import prrs.search as search
from ase.io import read
from prrs.runner import run_trial as real_run_trial

limit = int(sys.argv[3])
seen = {{"n": 0}}


def counted(*args, **kwargs):
    seen["n"] += 1
    if seen["n"] > limit:
        os._exit(137)
    return real_run_trial(*args, **kwargs)


search.run_trial = counted

from prrs.calculators import proton_relay_atoms, proton_relay_factory
from prrs.config import SearchConfig

config = SearchConfig(**{config!r})
search.ReactionSearch(proton_relay_factory, config).run(
    proton_relay_atoms((0, 0, 0)), sys.argv[1], product=read(sys.argv[2])
)
"""

BASE = dict(
    families=("stretch", "compress"),
    closed_shell_only=False,
    max_depth=3,
    max_trials=50,
    seed=7,
)


@pytest.mark.parametrize("connect", [False, True], ids=["discover", "connect"])
def test_killed_and_resumed_on_the_relay_is_byte_identical(tmp_path, connect):
    """Killed at trial 36. In connect mode the reactant took 32 trials and the second
    I-a choice -- the d = 2 microstate out of [4, 2, 4, 4], not the first-queued -- was
    made at 32, so the kill lands after a real choice with the queue still waiting; the
    queue persists in insertion order and the rule is deterministic on it."""
    import json
    import subprocess
    import sys
    from pathlib import Path

    from ase.io import read, write

    root = Path(__file__).resolve().parents[1]
    write(str(tmp_path / "product.extxyz"), _product())
    product = read(str(tmp_path / "product.extxyz"))
    settings = dict(BASE, connect_queue_priority=connect, connect_direction_priority=connect)
    config = SearchConfig(**settings)
    full = ReactionSearch(proton_relay_factory, config).run(
        proton_relay_atoms((0, 0, 0)), tmp_path / "full", product=product
    )
    assert full["status"] == "completed"
    if connect:
        assert len(full["connect_queue_order"]) >= 2

    script = tmp_path / "child.py"
    script.write_text(CHILD.format(root=str(root), config=settings), encoding="utf-8")
    killed = tmp_path / "killed"
    subprocess.run(
        [sys.executable, str(script), str(killed), str(tmp_path / "product.extxyz"), "36"],
        check=False,
        capture_output=True,
    )
    mid = json.loads((killed / "network.json").read_text(encoding="utf-8-sig"))
    assert mid["status"] == "running"
    assert mid.get("queue"), "the kill must happen with the queue still populated"
    if connect:
        assert len(mid.get("connect_queue_order", [])) >= 2, "and after I-a has chosen"

    resumed = ReactionSearch(proton_relay_factory, config).run(
        proton_relay_atoms((0, 0, 0)), killed, resume=True, product=product
    )
    assert resumed["status"] == "completed"
    for name in ("network.json", "attempts.jsonl"):
        assert (killed / name).read_bytes() == (tmp_path / "full" / name).read_bytes(), name
    # the reason that used to come back as "unknown"
    assert "unknown" not in resumed["termination_detail"]["directions_dropped_by"]
