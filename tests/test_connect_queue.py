"""Connect mode, G2 / intervention I-a (PLAN_2026-09-28.md): closest to the product first.

Off (the default) the queue is first-in first-out exactly as before. On, every choice is
recorded with the distances it chose among, so "it always took the closest" is read off
the record rather than trusted, and nothing is ever dropped from the queue.
"""

from dataclasses import replace

import pytest
from ase.optimize import LBFGS

from prrs.calculators import ProtonRelay, proton_relay_atoms, proton_relay_factory
from prrs.config import SearchConfig
from prrs.search import ReactionSearch, _closest_first

CONFIG = SearchConfig(
    families=("stretch", "compress"),
    closed_shell_only=False,
    max_depth=3,
    max_trials=45,
    seed=7,
)


def _product():
    atoms = proton_relay_atoms((1, 1, 1))
    atoms.calc = ProtonRelay()
    LBFGS(atoms, logfile=None).run(fmax=1e-4, steps=2000)
    atoms.calc = None
    return atoms


@pytest.fixture(scope="module")
def connect_run(tmp_path_factory):
    config = replace(CONFIG, connect_queue_priority=True)
    return ReactionSearch(proton_relay_factory, config).run(
        proton_relay_atoms((0, 0, 0)), tmp_path_factory.mktemp("ia") / "run", product=_product()
    )


def test_every_choice_is_the_closest_queued_microstate_first_queued_among_equals(connect_run):
    order = connect_run["connect_queue_order"]
    assert len(order) >= 2, "the run must have chosen between queued microstates"
    for choice in order:
        closest = min(choice["queued_d"])
        assert choice["d_to_product"] == closest
        assert choice["position"] == choice["queued_d"].index(closest)


def test_the_recorded_distance_is_the_microstates_own(connect_run):
    d_of = {
        state["id"]: state["d_to_product"]
        for node in connect_run["chemical_nodes"]
        for state in node["microstates"]
    }
    for choice in connect_run["connect_queue_order"]:
        assert choice["d_to_product"] == d_of[choice["microstate"]]


def test_a_farther_microstate_is_kept_queued_not_dropped(connect_run):
    """The upper layer's guard (plan section 1 item 4): a choice made with something
    farther in the queue leaves it in the queue."""
    passed_over = [
        c for c in connect_run["connect_queue_order"] if max(c["queued_d"]) > c["d_to_product"]
    ]
    assert passed_over, "the run must have passed over a farther microstate at least once"
    for choice in passed_over:
        assert len(choice["queued_d"]) >= 2


def test_off_nothing_is_recorded_and_the_queue_is_fifo(tmp_path):
    config = replace(CONFIG, max_trials=12)
    network = ReactionSearch(proton_relay_factory, config).run(
        proton_relay_atoms((0, 0, 0)), tmp_path / "run", product=_product()
    )
    assert "connect_queue_order" not in network


def test_the_switch_without_a_product_refuses_to_start(tmp_path):
    config = replace(CONFIG, connect_queue_priority=True)
    with pytest.raises(ValueError, match="product"):
        ReactionSearch(proton_relay_factory, config).run(
            proton_relay_atoms((0, 0, 0)), tmp_path / "run"
        )
    assert not (tmp_path / "run").exists()


def test_equal_distances_are_served_first_queued_first():
    assert _closest_first([4, 2, 2, 6]) == 1
    assert _closest_first([2, 2, 2]) == 0
    assert _closest_first([6]) == 0
