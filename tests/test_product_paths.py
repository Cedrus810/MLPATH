"""Connect mode, G4 (PLAN_2026-09-28.md): the R -> P path, read off the network, recorded.

Two views, never merged: `observed` walks only directed edges a trial itself travelled
(the strongest evidence the network has), `any_evidence` also walks continuations and
the undirected edges a saddle's two-sided descent supports. Every edge of a path carries
its own evidence lists, so a path through a TS-only edge says so on that edge. A product
not reached is an answer too: the closest node and its distance.
"""

import pytest

from prrs.product import product_paths


def _node(node_id, depth, compatible, d):
    return {
        "id": node_id,
        "depth": depth,
        "product_compatible": compatible,
        "microstates": [{"id": f"{node_id}/m0000", "d_to_product": d}],
    }


def _reaction(source, target, attempts=(), continuations=(), ts=None):
    edge = {
        "source": source,
        "target": target,
        "attempts": list(attempts),
        "continuations": [
            {"seed_attempt": c, "mechanism": "boundary_continuation"} for c in continuations
        ],
    }
    if ts:
        edge["ts_support"] = list(ts)
    return edge


def _network(nodes, reactions=(), channels=()):
    return {
        "chemical_nodes": nodes,
        "reactions": list(reactions),
        "ts_connected_channels": [{"ts": ts, "ends": list(ends)} for ts, ends in channels],
    }


NODES = [
    _node("c0000", 0, False, 6),
    _node("c0001", 1, False, 4),
    _node("c0002", 1, False, 2),
    _node("c0003", 2, True, 0),
]


def test_an_observed_path_lists_every_edge_with_its_own_evidence():
    network = _network(
        NODES,
        [
            _reaction("c0000", "c0002", attempts=["t000011"], ts=["ts0001"]),
            _reaction("c0002", "c0003", attempts=["t000033"]),
            _reaction("c0000", "c0001", attempts=["t000004"]),
        ],
        [("ts0001", ("c0000", "c0002"))],
    )
    paths = product_paths(network)
    observed = paths["observed"]
    assert observed["connected"] is True
    assert [(e["from"], e["to"]) for e in observed["edges"]] == [
        ("c0000", "c0002"),
        ("c0002", "c0003"),
    ]
    first = observed["edges"][0]
    assert first["attempts"] == ["t000011"] and first["ts"] == ["ts0001"]
    assert first["continuations"] == []
    assert observed["edges"][1]["ts"] == []


def test_a_ts_only_link_connects_any_evidence_but_not_observed():
    network = _network(
        NODES,
        [_reaction("c0000", "c0002", attempts=["t000011"])],
        [("ts0005", ("c0003", "c0002"))],  # undirected, and listed the other way round
    )
    paths = product_paths(network)
    assert paths["observed"]["connected"] is False
    anyway = paths["any_evidence"]
    assert anyway["connected"] is True
    last = anyway["edges"][-1]
    assert (last["from"], last["to"]) == ("c0002", "c0003")
    assert last["attempts"] == [] and last["ts"] == ["ts0005"]


def test_a_continuation_counts_for_any_evidence_only():
    network = _network(
        NODES,
        [_reaction("c0000", "c0003", continuations=["t000020"])],
    )
    paths = product_paths(network)
    assert paths["observed"]["connected"] is False
    edge = paths["any_evidence"]["edges"][0]
    assert edge["continuations"] == [
        {"seed_attempt": "t000020", "mechanism": "boundary_continuation"}
    ]


def test_a_reverse_observation_is_not_a_forward_one():
    network = _network(NODES, [_reaction("c0003", "c0000", attempts=["t000009"])])
    assert product_paths(network)["observed"]["connected"] is False


def test_not_reached_names_the_closest_node():
    network = _network(NODES[:3], [_reaction("c0000", "c0002", attempts=["t000011"])])
    paths = product_paths(network)
    assert paths["observed"]["connected"] is False
    assert paths["closest"] == {"nodes": ["c0002"], "d_to_product": 2}


def test_the_shortest_path_is_chosen_and_ties_break_by_node_id():
    network = _network(
        NODES,
        [
            _reaction("c0000", "c0002", attempts=["t1"]),
            _reaction("c0000", "c0001", attempts=["t2"]),
            _reaction("c0001", "c0003", attempts=["t3"]),
            _reaction("c0002", "c0003", attempts=["t4"]),
        ],
    )
    edges = product_paths(network)["observed"]["edges"]
    assert [e["to"] for e in edges] == ["c0001", "c0003"]  # c0001 < c0002, both length 2


def test_a_network_without_a_root_is_refused():
    network = _network(NODES[1:])
    with pytest.raises(ValueError, match="root"):
        product_paths(network)


def test_a_product_run_publishes_its_path(tmp_path):
    """The double well: the dissociated well is one trial away from the bonded one."""
    from prrs.calculators import demo_atoms, double_well_factory
    from prrs.config import SearchConfig
    from prrs.schema import validate
    from prrs.search import ReactionSearch

    config = SearchConfig(
        families=("stretch", "compress"),
        geometry_amplitudes_A=(0.1, 0.3, 0.7, 1.0),
        timestep_fs=0.1,
        response_steps=20,
        quench_fmax_eV_A=0.001,
        closed_shell_only=False,
        max_depth=2,
        max_trials=8,
        refinement_steps=2,
    )
    product = demo_atoms()
    product.positions[:, 0] *= 2.0
    network = ReactionSearch(double_well_factory, config).run(
        demo_atoms(), tmp_path / "run", product=product
    )
    observed = network["product_paths"]["observed"]
    assert observed["connected"] is True
    assert observed["edges"][0]["from"] == "c0000"
    assert observed["edges"][-1]["to"] == observed["target"]
    assert all(edge["attempts"] for edge in observed["edges"])
    assert network["product_paths"]["closest"]["d_to_product"] == 0
    assert validate(network) == []
    del network["product_paths"]
    assert "top: product run missing product_paths" in validate(network)
