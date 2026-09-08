from dataclasses import replace
import json
import numpy as np
import pytest
from ase.calculators.calculator import all_changes
from ase.io import read
from prrs.calculators import DoubleWell, demo_atoms, double_well_factory
from prrs.config import SearchConfig
from prrs.perturbations import Probe
from prrs.runner import relax_source, run_trial
from prrs.search import ReactionSearch
from prrs.cli import export_path


@pytest.fixture
def protocol():
    return SearchConfig(
        families=("stretch", "compress"),
        geometry_amplitudes_A=(0.1, 0.3, 0.7, 1.0),
        timestep_fs=0.1,
        response_steps=20,
        quench_fmax_eV_A=0.001,
        # The synthetic double well is a software benchmark, not carbon
        # chemistry -- its own docstring says so -- and C2 has no neutral
        # closed-shell valence assignment, so the domain gate has nothing
        # meaningful to say about it.
        closed_shell_only=False,
        max_depth=2,
        max_trials=20,
        refinement_steps=2,
    )


def test_hidden_product_discovery_dedup_directionality_and_export(tmp_path, protocol):
    output = tmp_path / "search"
    network = ReactionSearch(double_well_factory, protocol).run(demo_atoms(), output)
    # The two minima are different chemistry: bonded C2 versus two separated atoms.
    assert len(network["chemical_nodes"]) == 2
    assert sorted(node["fragments"] for node in network["chemical_nodes"]) == [
        ["C", "C"],
        ["C2"],
    ]
    microstates = [state for node in network["chemical_nodes"] for state in node["microstates"]]
    distances = sorted(
        read(str(output / state["structure"])).get_distance(0, 1) for state in microstates
    )
    assert distances == pytest.approx([1.2, 2.4], abs=0.001)
    # The barrier top is a first-order saddle and belongs in the candidate pool.
    for candidate in network["ts_candidates"]:
        assert candidate["saddle_order"] == 1
        assert candidate["imaginary_wavenumbers_icm"][0] < 0
        assert read(str(output / candidate["structure"])).get_distance(0, 1) == pytest.approx(
            1.8, abs=0.01
        )
    # Reverse edge must have independently executed evidence from the other node.
    for edge in network["reactions"]:
        assert edge["source"] != edge["target"]
        assert edge["observed_only"]
        for trial_id in edge["attempts"]:
            record = json.loads((output / "trials" / trial_id / "result.json").read_text())
            assert record["source"] == edge["source"]
            assert record["target"] == edge["target"]
            assert record["quench_converged"]
            assert record["tail_persistence_steps"] >= 3
    assert network["attempts"] <= protocol.max_trials
    attempts = [
        json.loads(line) for line in (output / "attempts.jsonl").read_text().splitlines()
    ]
    assert len(attempts) == network["attempts"]
    assert any(record["stage"] == "refine" for record in attempts)
    path = tmp_path / "path.traj"
    export_path(output, network["reactions"][0]["attempts"][0], path, 7)
    frames = read(str(path), index=":")
    assert len(frames) == 7
    assert frames[0].get_distance(0, 1) == pytest.approx(1.2, abs=0.001)
    assert frames[-1].get_distance(0, 1) == pytest.approx(2.4, abs=0.001)
    assert all(frame.info["purpose"] == "unrefined_path_seed" for frame in frames)


def test_seeded_search_repeats_and_never_overwrites(tmp_path, protocol):
    config = replace(protocol, max_depth=1, max_trials=8, refinement_steps=0)
    a = ReactionSearch(double_well_factory, config).run(demo_atoms(), tmp_path / "a")
    b = ReactionSearch(double_well_factory, config).run(demo_atoms(), tmp_path / "b")
    assert a == b
    with pytest.raises(FileExistsError):
        ReactionSearch(double_well_factory, config).run(demo_atoms(), tmp_path / "a")


def test_low_amplitude_is_elastic(tmp_path, protocol):
    config = replace(protocol, geometry_amplitudes_A=(0.05, 0.1), max_depth=1)
    network = ReactionSearch(double_well_factory, config).run(
        demo_atoms(), tmp_path / "elastic"
    )
    assert len(network["chemical_nodes"]) == 1
    assert network["chemical_nodes"][0]["reservoir"]["occupied"] == 1
    assert network["reactions"] == []
    assert network["classification_counts"] == {"elastic": 2}
    assert network["admission_counts"] == {"existing_microstate": 2}


class BrokenCalculator(DoubleWell):
    def calculate(
        self, atoms=None, properties=("energy", "forces"), system_changes=all_changes
    ):
        super().calculate(atoms, properties, system_changes)
        if atoms.get_distance(0, 1) > 1.4:
            self.results["forces"][:] = np.nan


def test_unphysical_trial_is_saved_and_never_enters_network(tmp_path, protocol):
    network = ReactionSearch(BrokenCalculator, replace(protocol, max_depth=1)).run(
        demo_atoms(), tmp_path / "nan"
    )
    assert len(network["chemical_nodes"]) == 1
    assert not network["reactions"]
    attempts = [
        json.loads(line)
        for line in (tmp_path / "nan" / "attempts.jsonl").read_text().splitlines()
    ]
    assert any(r["failure"] and r["failure"]["code"] == "nonfinite" for r in attempts)
    for r in attempts:
        if r["status"] == "rejected":
            assert r["target"] is None
            assert (tmp_path / "nan" / r["rejected_structure"]).exists()


def test_failed_quench_is_not_a_basin(tmp_path, protocol):
    source = relax_source(demo_atoms(), double_well_factory, protocol)
    outcome = run_trial(
        source,
        Probe("stretch", (0, 1), 1, 0.7, 12),
        double_well_factory,
        replace(protocol, quench_steps=1, quench_fmax_eV_A=1e-12),
        tmp_path / "trial",
        "t000000",
    )
    assert outcome.record["status"] == "quench_failed"
    assert outcome.endpoint is None


def test_pulse_is_off_in_free_phase(tmp_path, protocol):
    source = relax_source(demo_atoms(), double_well_factory, protocol)
    result = run_trial(
        source,
        Probe("pulse", (0, 1), 1, 0.6, 12),
        double_well_factory,
        replace(protocol, pulse_steps=10),
        tmp_path / "pulse",
        "t000000",
    )
    assert result.record["status"] == "completed"
    observations = [
        json.loads(line)
        for line in (tmp_path / "pulse" / "observations.jsonl").read_text().splitlines()
    ]
    assert any(row["phase"] == "forced" and row["bias_energy_eV"] != 0 for row in observations)
    assert all(
        row["bias_energy_eV"] == 0 for row in observations if row["phase"].startswith("free")
    )
    assert abs(result.record["pulse_work_residual_eV"]) < 1e-5


def test_saddle_following_names_the_endpoints_without_inventing_a_direction(tmp_path, protocol):
    """A shared saddle is undirected evidence and must not create a directed edge.

    Descending both sides of the double well's barrier top has to land on the bonded pair
    and on the separated pair. That identifies the channel's ends; it says nothing about
    which direction was observed, so the directed edges still need their own trials.
    """
    output = tmp_path / "irc"
    network = ReactionSearch(double_well_factory, protocol).run(demo_atoms(), output)
    assert network["ts_candidates"]
    assert network["ts_connected_channels"]

    node_ids = {node["id"] for node in network["chemical_nodes"]}
    for channel in network["ts_connected_channels"]:
        assert channel["undirected"] is True
        assert set(channel["ends"]) <= node_ids
        assert len(set(channel["ends"])) == 2  # the two minima really differ
    for candidate in network["ts_candidates"]:
        ends = candidate["connects"]
        assert len(ends) == 2 and {e["sign"] for e in ends} == {1, -1}
        assert all("barrier" not in candidate["connects_meaning"].split("no ")[0] for _ in ends)

    # Saddle support annotates the directed edges but never substitutes for them.
    for edge in network["reactions"]:
        # Two kinds of directed evidence, and each must agree with its own trial record.
        # A trial in `attempts` observed the transition itself, so its record names the same
        # target. A boundary continuation's seed trial did not -- it stopped at the boundary,
        # which is why the continuation ran -- so it is listed apart and must NOT claim the
        # target. Collapsing the two would let a seed trial be read as an observation.
        assert edge["attempts"] or edge["continuations"], (
            "a directed edge still needs directed evidence of some kind"
        )
        assert edge["observed_only"]
        assert edge["mechanisms"]
        for trial_id in edge["attempts"]:
            record = json.loads((output / "trials" / trial_id / "result.json").read_text())
            assert record["source"] == edge["source"] and record["target"] == edge["target"]
        for continuation in edge["continuations"]:
            record = json.loads(
                (output / "trials" / continuation["seed_attempt"] / "result.json").read_text()
            )
            assert record["source"] == edge["source"]
            assert record["target"] != edge["target"], (
                "a continuation seed that had already reached the target would not "
                "need continuing, and would belong in attempts"
            )
            assert continuation["mechanism"] == "boundary_continuation"
    # Support is symmetric: both directions of one channel cite the same saddles.
    supports = {
        frozenset((e["source"], e["target"])): tuple(e.get("ts_support", ()))
        for e in network["reactions"]
    }
    assert len(set(supports.values())) == 1


def test_saddle_following_can_be_switched_off(tmp_path, protocol):
    config = replace(protocol, irc_enabled=False)
    network = ReactionSearch(double_well_factory, config).run(demo_atoms(), tmp_path / "noirc")
    assert network["ts_candidates"]
    assert network["ts_connected_channels"] == []
    assert all("connects" not in candidate for candidate in network["ts_candidates"])


class SharpWell(DoubleWell):
    """The same two minima with a narrow bump on the ridge, so the saddle is sharp.

    The double well's barrier top has lambda_1 = -0.74 and a quench started near it stalls
    there, which is how every transition state in these tests and in the ethanol benchmark
    was found. That is a property of soft reaction coordinates, not an algorithm: on
    malonaldehyde the proton-transfer saddle sits at -37.3, sixty probes quenched into
    basins every time, and the saddle was never seen. This potential reproduces that
    regime -- lambda_1 = -14.5, twenty times sharper -- with everything else unchanged.
    """

    HEIGHT, WIDTH = 0.6, 0.12

    def calculate(
        self, atoms=None, properties=("energy", "forces"), system_changes=all_changes
    ):
        super().calculate(atoms, properties, system_changes)
        vector = atoms.positions[1] - atoms.positions[0]
        r = float(np.linalg.norm(vector))
        u = (r - 1.8) / self.WIDTH
        bump = self.HEIGHT * np.exp(-u * u)
        derivative = -2 * u * bump / self.WIDTH
        self.results["energy"] += float(bump)
        force = derivative * vector / r
        self.results["forces"] = self.results["forces"] + np.array([force, -force])


@pytest.fixture
def sharp_protocol():
    return SearchConfig(
        families=("stretch", "compress"),
        geometry_amplitudes_A=(0.1, 0.3, 0.7, 1.0),
        timestep_fs=0.1,
        response_steps=20,
        quench_fmax_eV_A=0.001,
        closed_shell_only=False,  # see the fixture above
        max_depth=1,
        max_trials=20,
        refinement_steps=2,
    )


def test_a_sharper_ridge_is_still_reached_by_climbing(tmp_path, sharp_protocol):
    """With the climb switched off this potential still yields a saddle -- one internal
    coordinate makes stalling easy in a way a nine-atom molecule does not -- so the contrast
    here is not "no saddle at all" but where the saddle came from, which is recorded."""
    config = replace(sharp_protocol, saddle_search_enabled=False)
    network = ReactionSearch(SharpWell, config).run(demo_atoms(), tmp_path / "sharp_off")
    assert network["reaction_channels"]
    assert all(
        candidate["found_by"] == "quench_stalled" for candidate in network["ts_candidates"]
    )
    assert all(
        scan["saddle_search"]["attempted"] is False for scan in network["amplitude_scans"]
    )


def test_a_reacting_scan_climbs_to_its_saddle_and_resolves_both_ends(tmp_path, sharp_protocol):
    """The closed loop: observed reaction -> probe frame -> climb -> saddle -> endpoints."""
    output = tmp_path / "sharp_on"
    network = ReactionSearch(SharpWell, sharp_protocol).run(demo_atoms(), output)

    assert network["reaction_channels"]
    attempted = [
        scan["saddle_search"]
        for scan in network["amplitude_scans"]
        if scan["saddle_search"]["attempted"]
    ]
    assert attempted, "a scan that reacted must attempt the climb"

    climbed = [c for c in network["ts_candidates"] if c["found_by"] == "min_mode_following"]
    assert climbed, [s for s in attempted]
    for candidate in climbed:
        assert candidate["saddle_order"] == 1
        assert candidate["imaginary_wavenumbers_icm"][0] < 0
        # It is the ridge of this potential, not some other stationary point.
        assert read(str(output / candidate["structure"])).get_distance(0, 1) == pytest.approx(
            1.8, abs=0.02
        )
        # Both sides were descended and they land on the two different minima.
        ends = candidate["connects"]
        assert len(ends) == 2 and {e["sign"] for e in ends} == {1, -1}
        assert len({e.get("chemical_node") for e in ends}) == 2

    # The channel cites the saddle that was climbed for it.
    supported = [c for c in network["reaction_channels"] if c["ts_support"]]
    assert supported
    assert any(
        candidate["id"] in channel["ts_support"]
        for channel in supported
        for candidate in climbed
    )


def test_a_transient_product_is_checked_but_never_admitted(tmp_path, sharp_protocol):
    """A crossing that did not hold hands over its endpoint and keeps its claim out.

    The endpoint still goes through the minimum check, because a saddle is a saddle whether
    or not the product was stable, and the amplitude still brackets the dividing surface.
    What it must not do is enter the network as a product.
    """
    output = tmp_path / "transient"
    ReactionSearch(SharpWell, sharp_protocol).run(demo_atoms(), output)
    records = [
        json.loads(line) for line in (output / "attempts.jsonl").read_text().splitlines()
    ]
    transient = [r for r in records if (r.get("crossing") or {}).get("observed")]
    for record in transient:
        assert record["status"] == "unconfirmed"
        assert record["admission"] in ("withheld_untrusted_topology", "ts_candidate")
        assert record["target_microstate"] is None
        assert record["minimum_confirmation"] is not None, "the endpoint was still checked"


def test_a_run_that_dies_inside_publish_still_records_why(tmp_path, protocol):
    """The failure record must not depend on the thing that failed.

    The old handler called `publish()` and then wrote the manifest. When the exception
    CAME FROM publish() -- a non-finite value somewhere in the network, which
    `io.dumps_strict` refuses on purpose -- the handler raised the same error again and the
    manifest never got written. `status` stayed "running", which reads like the process was
    killed by a signal. Three seeds of the T=300 scan died exactly that way and it cost a
    wrong diagnosis; `io.py` records the incident, and the error message was improved, but
    the handler kept the same order.

    Smallest record first, largest last, each in its own try. `failure.json` holds only
    strings, so it survives whatever the network could not serialise.
    """
    import json
    from pathlib import Path

    import prrs.search as search_module

    output = tmp_path / "dies_in_publish"
    real = search_module.atomic_json
    state = {"armed": False}

    def explode_on_network(path, data):
        if state["armed"] and Path(path).name == "network.json":
            raise ValueError("synthetic: non-finite value in the network")
        return real(path, data)

    search_module.atomic_json = explode_on_network
    try:
        state["armed"] = True
        with pytest.raises(ValueError, match="synthetic"):
            ReactionSearch(double_well_factory, protocol).run(demo_atoms(), output)
    finally:
        search_module.atomic_json = real

    report = json.loads((output / "failure.json").read_text())
    assert report["type"] == "ValueError"
    assert "synthetic" in report["message"]
    assert report["traceback"], "a failure record without a traceback is not diagnosable"
    # The snapshot failed too, and the record has to say so rather than leaving a reader
    # to trust a stale network.json.
    assert "network_snapshot_failed" in report, report
    # The manifest is small and string-only, so it must have survived.
    manifest = json.loads((output / "config.json").read_text())
    assert manifest["status"] == "failed", manifest["status"]
