import json

import numpy as np
import pytest
from ase import Atoms
from molecules import chloride_methyl_chloride, ethanol, relabel, water_dimer
from prrs.calculators import double_well_factory
from prrs.chemistry import canonical_labels
from prrs.config import SearchConfig
from prrs.perturbations import rotatable_torsions
from prrs.state import encode

import explore.poses as poses
from explore.ledger import Ledger, counted_factory
from explore.poses import (
    Conformer,
    _clashes,
    conformer_grid,
    conformers_from_qualification,
    contact_ratio,
    fragments,
    place,
    qualify_conformers,
    sample_poses,
    single_fragment_poses,
    write_l0,
)
from explore.streams import stream

CONFIG = SearchConfig(closed_shell_only=False)


def test_stream_is_keyed_not_ordered():
    a = stream(0, "poses", "x", 3).normal(size=4)
    b = stream(0, "poses", "x", 3).normal(size=4)
    c = stream(0, "poses", "x", 4).normal(size=4)
    assert np.array_equal(a, b) and not np.array_equal(a, c)


def test_fragments_are_ordered_by_identity_not_index():
    atoms = chloride_methyl_chloride()
    frags = fragments(atoms, CONFIG)
    moved, perm = relabel(atoms, seed=3)
    frags_moved = fragments(moved, CONFIG)
    assert [f.key for f in frags] == [f.key for f in frags_moved]
    for f, g in zip(frags, frags_moved):
        assert sorted(int(perm[i]) for i in g.atoms) == list(f.atoms)


def test_identical_fragments_are_ordered_by_geometry():
    atoms = water_dimer()
    frags = fragments(atoms, CONFIG)
    moved, perm = relabel(atoms, seed=5)
    frags_moved = fragments(moved, CONFIG)
    for f, g in zip(frags, frags_moved):
        assert sorted(int(perm[i]) for i in g.atoms) == list(f.atoms)


def test_ethanol_grid_has_one_genuine_torsion_and_six_points():
    atoms = ethanol()
    (frag,) = fragments(atoms, CONFIG)
    confs, report = conformer_grid(atoms, frag, CONFIG, steps_per_turn=6, cap=1000, seed=0)
    assert report["genuine_torsions"] == 1 and report["rotors_not_enumerated"] == 1
    assert len(confs) == 6
    zero = next(c for c in confs if c.grid == (0,))
    assert np.allclose(zero.positions, atoms.positions[list(frag.atoms)], atol=1e-10)
    graph = encode(atoms, CONFIG.bond_scale)
    assert all(c.status == "generated" for c in confs)
    for c in confs:
        moved = atoms.copy()
        moved.positions[list(frag.atoms)] = c.positions
        assert encode(moved, CONFIG.bond_scale).edges == graph.edges


def test_grid_dihedrals_step_by_sixty_degrees():
    atoms = ethanol()
    (frag,) = fragments(atoms, CONFIG)
    confs, _ = conformer_grid(atoms, frag, CONFIG, steps_per_turn=6, cap=1000, seed=0)
    graph = encode(atoms, CONFIG.bond_scale)
    labels = canonical_labels(atoms.numbers, graph.edges, graph.index)
    genuine, _ = rotatable_torsions(graph, labels, {int(i): int(i) for i in graph.index})
    idx = list(genuine[0].indices)
    base = atoms.get_dihedral(*idx)
    seen = set()
    for c in confs:
        moved = atoms.copy()
        moved.positions[list(frag.atoms)] = c.positions
        seen.add(round(((moved.get_dihedral(*idx) - base) % 360.0) / 60.0) % 6)
    assert seen == set(range(6))


def test_monatomic_and_diatomic_fragments_have_only_the_input():
    atoms = chloride_methyl_chloride()
    frags = fragments(atoms, CONFIG)
    cl = next(f for f in frags if len(f.atoms) == 1)
    confs, report = conformer_grid(atoms, cl, CONFIG, steps_per_turn=6, cap=10, seed=0)
    assert len(confs) == 1 and confs[0].grid == () and report["genuine_torsions"] == 0
    assert np.isfinite(confs[0].positions).all()


def test_cap_keeps_input_and_reports_censoring():
    atoms = ethanol()
    (frag,) = fragments(atoms, CONFIG)
    confs, report = conformer_grid(atoms, frag, CONFIG, steps_per_turn=6, cap=3, seed=0)
    assert len(confs) == 3 and any(c.grid == (0,) for c in confs)
    assert report["grid_size"] == 6 and report["censored"] == 3


def test_clashing_grid_point_is_kept_with_status(monkeypatch):
    atoms = ethanol()
    (frag,) = fragments(atoms, CONFIG)
    original = atoms.positions[list(frag.atoms)].copy()

    def fake(positions, numbers, config):
        return not np.allclose(positions, original, atol=1e-8)

    monkeypatch.setattr(poses, "_clashes", fake)
    confs, report = conformer_grid(atoms, frag, CONFIG, steps_per_turn=6, cap=1000, seed=0)
    assert len(confs) == 6
    assert next(c for c in confs if c.grid == (0,)).status == "generated"
    assert all(c.status == "overlap" for c in confs if c.grid != (0,))
    assert report["status_counts"]["overlap"] == 5


def test_clashes_helper_uses_min_distance_ratio():
    numbers = np.array([6, 6])
    close = np.array([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]])
    far = np.array([[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]])
    assert _clashes(close, numbers, CONFIG) is True
    assert _clashes(far, numbers, CONFIG) is False


def test_changed_bond_graph_is_kept_with_status(monkeypatch):
    import types

    import prrs.state

    atoms = ethanol()
    (frag,) = fragments(atoms, CONFIG)
    original = atoms.positions[list(frag.atoms)].copy()
    real = prrs.state.encode

    def fake(sub, bond_scale):
        graph = real(sub, bond_scale)
        if np.allclose(sub.positions, original, atol=1e-8):
            return graph
        return types.SimpleNamespace(edges=("changed",), index=graph.index)

    monkeypatch.setattr(prrs.state, "encode", fake)
    confs, report = conformer_grid(atoms, frag, CONFIG, steps_per_turn=6, cap=1000, seed=0)
    assert len(confs) == 6
    assert next(c for c in confs if c.grid == (0,)).status == "generated"
    assert all(c.status == "graph_changed" for c in confs if c.grid != (0,))
    assert report["status_counts"] == {"generated": 1, "graph_changed": 5}


def test_diatomic_fragment_has_only_the_input():
    from ase import Atoms

    atoms = Atoms("OHCl", positions=[[0, 0, 0], [0.0, 0.0, 0.97], [0.0, 0.0, 6.0]])
    frags = fragments(atoms, CONFIG)
    oh = next(f for f in frags if len(f.atoms) == 2)
    confs, report = conformer_grid(atoms, oh, CONFIG, steps_per_turn=6, cap=10, seed=0)
    assert len(confs) == 1 and confs[0].grid == () and report["genuine_torsions"] == 0
    assert np.isfinite(confs[0].positions).all()
    assert np.allclose(confs[0].positions, atoms.positions[list(oh.atoms)], atol=1e-10)


def test_ethanol_grid_geometries_survive_relabelling():
    atoms = ethanol()
    moved, perm = relabel(atoms, seed=7)

    def geometries(a, perm_map):
        (frag,) = fragments(a, CONFIG)
        confs, _ = conformer_grid(a, frag, CONFIG, steps_per_turn=6, cap=1000, seed=0)
        out = []
        for c in confs:
            full = np.zeros_like(a.positions)
            full[list(frag.atoms)] = c.positions
            # map back to original labelling: original atom perm[k] is new atom k
            back = np.zeros_like(full)
            back[perm_map] = full
            out.append(back)
        return out

    ident = np.arange(len(atoms))
    g0 = geometries(atoms, ident)
    g1 = geometries(moved, perm)
    assert len(g0) == len(g1) == 6
    for x in g1:
        assert any(np.allclose(x, y, atol=1e-8) for y in g0)


DW = SearchConfig(closed_shell_only=False, quench_fmax_eV_A=0.001)
FS = {"C2": {"charge": 0, "multiplicity": 1}}


def _c2(r, cid):
    return Conformer(cid, "C2", (0, 1), np.array([[0.0, 0, 0], [r, 0, 0]]), (), "generated")


def test_qualification_dedupes_same_basin_and_keeps_distinct(tmp_path):
    ledger = Ledger(tmp_path / "l.jsonl", arm="P", budget=5000)
    template = Atoms("C2", positions=[[0, 0, 0], [1.2, 0, 0]])
    rows = qualify_conformers(
        template,
        [_c2(1.25, "a"), _c2(1.30, "b"), _c2(2.35, "c"), _c2(1.2, "d")],
        counted_factory(double_well_factory, ledger),
        DW,
        ledger,
        seed=0,
        fragment_states=FS,
    )
    ledger.close()
    status = {r["conformer_id"]: r["status"] for r in rows}
    assert status["a"] == "qualified"
    assert status["b"] == "duplicate" and rows[1]["duplicate_of"] == "a"
    assert status["d"] == "duplicate"
    assert ledger.spent > 0
    stages = {
        e.get("stage") for e in map(json.loads, (tmp_path / "l.jsonl").read_text().splitlines())
    }
    assert stages <= {"prepare", None}
    # c is a different basin (and bond graph), reported as qualified or at least not a duplicate of a
    assert status["c"] == "qualified"


def test_equal_energy_but_different_bond_graph_is_not_a_duplicate(tmp_path):
    ledger = Ledger(tmp_path / "l.jsonl", arm="P", budget=5000)
    rows = qualify_conformers(
        Atoms("C2", positions=[[0, 0, 0], [1.2, 0, 0]]),
        [_c2(1.2, "short"), _c2(2.4, "long")],
        counted_factory(double_well_factory, ledger),
        DW,
        ledger,
        seed=0,
        fragment_states=FS,
    )
    ledger.close()
    e = {r["conformer_id"]: r["energy_eV"] for r in rows}
    assert abs(e["short"] - e["long"]) <= DW.basin_energy_eV  # energy alone would merge them
    ends = [Atoms("C2", positions=r["positions"]) for r in rows]
    assert encode(ends[0], DW.bond_scale).edges != encode(ends[1], DW.bond_scale).edges
    assert [r["status"] for r in rows] == ["qualified", "qualified"]


def test_grid_points_with_bad_status_are_skipped_without_forwards(tmp_path):
    ledger = Ledger(tmp_path / "l.jsonl", arm="P", budget=5000)
    bad = Conformer("x", "C2", (0, 1), np.array([[0.0, 0, 0], [0.3, 0, 0]]), (), "overlap")
    rows = qualify_conformers(
        Atoms("C2"),
        [bad],
        counted_factory(double_well_factory, ledger),
        DW,
        ledger,
        seed=0,
        fragment_states=FS,
    )
    assert rows[0]["status"] == "skipped:overlap" and ledger.spent == 0


def test_budget_exhaustion_is_a_status_not_an_exception(tmp_path):
    ledger = Ledger(tmp_path / "l.jsonl", arm="P", budget=3)
    rows = qualify_conformers(
        Atoms("C2"),
        [_c2(1.3, "a")],
        counted_factory(double_well_factory, ledger),
        DW,
        ledger,
        seed=0,
        fragment_states=FS,
    )
    assert rows[0]["status"] == "budget_exhausted"


def test_fragment_state_is_stamped_and_reaches_the_calculator(tmp_path):
    seen = []

    def recording_factory():
        calc = double_well_factory()
        original = calc.calculate

        def calculate(atoms=None, properties=("energy",), system_changes=()):
            seen.append((atoms.info.get("charge"), atoms.info.get("spin")))
            return original(atoms, properties, system_changes)

        calc.calculate = calculate
        return calc

    ledger = Ledger(tmp_path / "l.jsonl", arm="P", budget=5000)
    system = Atoms("C2", positions=[[0, 0, 0], [1.2, 0, 0]])
    system.info["charge"], system.info["spin"] = -1, 2
    rows = qualify_conformers(
        system,
        [_c2(1.25, "a")],
        counted_factory(recording_factory, ledger),
        DW,
        ledger,
        seed=0,
        fragment_states=FS,
    )
    ledger.close()
    assert rows[0]["status"] == "qualified"
    assert seen and set(seen) == {(0, 1)}


def test_missing_fragment_state_is_a_status_without_forwards(tmp_path):
    ledger = Ledger(tmp_path / "l.jsonl", arm="P", budget=5000)
    rows = qualify_conformers(
        Atoms("C2"),
        [_c2(1.3, "a")],
        counted_factory(double_well_factory, ledger),
        DW,
        ledger,
        seed=0,
        fragment_states={},
    )
    assert rows[0]["status"] == "electronic_state_missing" and ledger.spent == 0


def test_hessian_gate_reason_maps_to_gate_rejected(tmp_path, monkeypatch):
    import prrs.runner as runner

    monkeypatch.setattr(
        runner, "confirm_minimum", lambda *a, **k: (False, {"reason": "gate_overlap"})
    )
    ledger = Ledger(tmp_path / "l.jsonl", arm="P", budget=5000)
    rows = qualify_conformers(
        Atoms("C2", positions=[[0, 0, 0], [1.2, 0, 0]]),
        [_c2(1.25, "a")],
        counted_factory(double_well_factory, ledger),
        DW,
        ledger,
        seed=0,
        fragment_states=FS,
    )
    assert rows[0]["status"] == "gate_rejected:overlap"


def _pair():
    atoms = chloride_methyl_chloride()
    frags = fragments(atoms, CONFIG)
    confs = [conformer_grid(atoms, f, CONFIG, 6, 100, 0)[0][0] for f in frags]
    return atoms, frags, confs


def test_placement_hits_the_shell_ratio():
    atoms, frags, (ca, cb) = _pair()
    poses = sample_poses(atoms, frags, ca, cb, shells=(1.8, 2.4), n_samples=4, seed=0)
    assert len(poses) == 8
    for p in poses:
        a = p.positions[list(frags[0].atoms)]
        b = p.positions[list(frags[1].atoms)]
        r = contact_ratio(
            a, atoms.numbers[list(frags[0].atoms)], b, atoms.numbers[list(frags[1].atoms)]
        )
        assert abs(r - p.shell) < 1e-6 and p.status == "placed"


def test_fragments_move_rigidly():
    atoms, frags, (ca, cb) = _pair()
    for p in sample_poses(atoms, frags, ca, cb, shells=(2.0,), n_samples=3, seed=1):
        for f, c in zip(frags, (ca, cb)):
            got = p.positions[list(f.atoms)]
            d0 = np.linalg.norm(c.positions[:, None] - c.positions[None], axis=-1)
            d1 = np.linalg.norm(got[:, None] - got[None], axis=-1)
            assert np.allclose(d0, d1, atol=1e-10)


def test_appending_samples_keeps_existing_poses():
    atoms, frags, (ca, cb) = _pair()
    few = sample_poses(atoms, frags, ca, cb, shells=(2.0,), n_samples=3, seed=2)
    more = sample_poses(atoms, frags, ca, cb, shells=(2.0,), n_samples=6, seed=2)
    for p, q in zip(few, more):
        assert p.pose_id == q.pose_id and np.array_equal(p.positions, q.positions)


def test_relabelled_input_gives_the_same_poses():
    atoms = water_dimer()
    frags = fragments(atoms, CONFIG)
    confs = [conformer_grid(atoms, f, CONFIG, 6, 100, 0)[0][0] for f in frags]
    poses = sample_poses(atoms, frags, *confs, shells=(2.0,), n_samples=3, seed=4)
    moved, perm = relabel(atoms, seed=7)
    mfr = fragments(moved, CONFIG)
    mconfs = [conformer_grid(moved, f, CONFIG, 6, 100, 0)[0][0] for f in mfr]
    mposes = sample_poses(moved, mfr, *mconfs, shells=(2.0,), n_samples=3, seed=4)
    for p, q in zip(poses, mposes):
        assert np.allclose(q.positions, p.positions[perm], atol=1e-8)


def test_identical_fragments_get_distinct_conformer_and_pose_ids():
    atoms = water_dimer()
    frags = fragments(atoms, CONFIG)
    assert [f.ordinal for f in frags] == [0, 1]
    ca, cb = [conformer_grid(atoms, f, CONFIG, 6, 100, 0)[0][0] for f in frags]
    assert ca.conformer_id != cb.conformer_id
    poses = sample_poses(atoms, frags, ca, cb, shells=(1.8, 2.4), n_samples=3, seed=0)
    ids = [p.pose_id for p in poses]
    assert len(ids) == 6 and len(set(ids)) == 6


def test_three_fragments_are_refused():
    atoms = water_dimer() + water_dimer()[:3]
    atoms.positions[6:] += [0, 0, 8.0]
    frags = fragments(atoms, CONFIG)
    assert len(frags) == 3
    with pytest.raises(ValueError, match="exactly 2 fragments"):
        sample_poses(atoms, frags, None, None, shells=(2.0,), n_samples=1, seed=0)


def test_single_fragment_poses_are_the_conformers():
    atoms = ethanol()
    (frag,) = fragments(atoms, CONFIG)
    confs, _ = conformer_grid(atoms, frag, CONFIG, 6, 100, 0)
    poses = single_fragment_poses(atoms, confs)
    assert [p.conformers for p in poses] == [(c.conformer_id,) for c in confs]
    assert all(p.shell is None for p in poses)


def test_l0_record_has_the_contract_fields(tmp_path):
    atoms, frags, (ca, cb) = _pair()
    poses = sample_poses(atoms, frags, ca, cb, shells=(2.0,), n_samples=1, seed=0)
    write_l0(
        tmp_path / "l0.jsonl",
        poses,
        {"charge": -1, "multiplicity": 1},
        {"rule": "L0-v1", "seed": 0},
        [f.atoms for f in frags],
    )
    rec = json.loads((tmp_path / "l0.jsonl").read_text().splitlines()[0])
    for k in (
        "source_microstate",
        "conformer",
        "pose",
        "atoms_map",
        "electronic_state",
        "generation_rule",
        "seed",
        "preparation",
        "qualification",
        "cost",
        "positions",
    ):
        assert k in rec
    assert rec["preparation"] == "raw_placement"
    assert rec["atoms_map"] == {"fragments": [list(f.atoms) for f in frags]}
    assert "t_A" in rec


def test_place_hits_shell_when_max_radius_pair_is_far():
    pos = np.array([[0.0, 0.0, 0.0], [1.61, 0.0, 0.0]])
    num = np.array([1, 53])
    _, ratio, status, _ = place(
        pos, num, pos, num, np.array([1.0, 0, 0]), np.array([0, 0, 0, 1.0]), 2.4
    )
    assert status == "placed" and abs(ratio - 2.4) < 1e-6


def test_place_cap_returns_failed_without_hanging(monkeypatch):
    monkeypatch.setattr(poses, "contact_ratio", lambda *a: 5.0)
    pos = np.array([[0.0, 0.0, 0.0]])
    num = np.array([1])
    _, _, status, _ = place(
        pos, num, pos, num, np.array([1.0, 0, 0]), np.array([0, 0, 0, 1.0]), 2.0
    )
    assert status == "placement_failed"


def test_mismatched_conformer_atoms_are_refused():
    atoms = water_dimer()
    frags = fragments(atoms, CONFIG)
    ca, cb = [conformer_grid(atoms, f, CONFIG, 6, 100, 0)[0][0] for f in frags]
    with pytest.raises(ValueError, match="do not match"):
        sample_poses(atoms, frags, cb, ca, shells=(2.0,), n_samples=1, seed=0)


def test_non_generated_conformer_marks_the_pose():
    from dataclasses import replace

    atoms, frags, (ca, cb) = _pair()
    bad = replace(ca, status="overlap")
    (p,) = sample_poses(atoms, frags, bad, cb, shells=(2.0,), n_samples=1, seed=0)
    assert p.status == "conformer_overlap"
    (q,) = sample_poses(atoms, frags, ca, cb, shells=(2.0,), n_samples=1, seed=0)
    assert q.status == "placed" and q.t_A is not None


def _qualify(tmp_path, convs, budget=5000, **kw):
    ledger = Ledger(tmp_path / "l.jsonl", arm="P", budget=budget)
    rows = qualify_conformers(
        Atoms("C2", positions=[[0, 0, 0], [1.2, 0, 0]]),
        convs,
        counted_factory(double_well_factory, ledger),
        DW,
        ledger,
        seed=0,
        fragment_states=FS,
        **kw,
    )
    return rows, ledger


def test_quench_failed_row_keeps_geometry_and_diagnostics(tmp_path, monkeypatch):
    import prrs.runner as runner

    monkeypatch.setattr(
        runner, "quench", lambda atoms, guard, config: (False, {"rounds": np.array([1.0, 2.0])})
    )
    rows, _ = _qualify(tmp_path, [_c2(1.3, "a")])
    assert rows[0]["status"] == "quench_failed"
    assert np.allclose(rows[0]["positions"], [[0, 0, 0], [1.3, 0, 0]])
    json.dumps(rows[0]["diagnostics"])
    assert rows[0]["diagnostics"]["quench"]["rounds"] == [1.0, 2.0]


def test_not_minimum_row_keeps_geometry_and_diagnostics(tmp_path, monkeypatch):
    import prrs.runner as runner

    monkeypatch.setattr(
        runner,
        "confirm_minimum",
        lambda *a, **k: (False, {"reason": "negative_mode", "lam": np.float64(-0.5)}),
    )
    rows, _ = _qualify(tmp_path, [_c2(1.25, "a")])
    assert rows[0]["status"] == "not_minimum"
    assert rows[0]["positions"] is not None
    json.dumps(rows[0]["diagnostics"])
    assert rows[0]["diagnostics"]["confirm"]["lam"] == -0.5
    assert rows[0]["diagnostics"]["quench"] is not None


def test_unexpected_exception_is_a_status_and_the_loop_continues(tmp_path, monkeypatch):
    import prrs.runner as runner

    real = runner.quench
    calls = []

    def flaky(atoms, guard, config):
        calls.append(1)
        if len(calls) == 1:
            raise RuntimeError("boom")
        return real(atoms, guard, config)

    monkeypatch.setattr(runner, "quench", flaky)
    rows, _ = _qualify(tmp_path, [_c2(1.25, "a"), _c2(2.35, "c")])
    assert rows[0]["status"] == "error:RuntimeError"
    assert rows[1]["status"] == "qualified"


def test_attempt_reaches_the_ledger_segments(tmp_path):
    _, ledger = _qualify(tmp_path, [_c2(1.25, "a")], attempt=3)
    ledger.close()
    events = [json.loads(x) for x in (tmp_path / "l.jsonl").read_text().splitlines()]
    opens = [e for e in events if e.get("event") == "open"]
    assert opens and all(e["attempt"] == 3 for e in opens)


def test_conformers_from_qualification_carry_relaxed_positions(tmp_path):
    convs = [_c2(1.25, "a"), _c2(1.30, "b"), _c2(2.35, "c")]
    rows, _ = _qualify(tmp_path, convs)
    out = conformers_from_qualification(convs, rows)
    assert [c.conformer_id for c in out] == ["a", "c"]  # b is a duplicate
    for c, src in zip(out, (convs[0], convs[2])):
        assert c.status == "qualified"
        assert (c.fragment_key, c.atoms, c.grid) == (src.fragment_key, src.atoms, src.grid)
    by_id = {r["conformer_id"]: r for r in rows}
    assert np.allclose(out[0].positions, by_id["a"]["positions"])
    assert not np.allclose(out[0].positions, convs[0].positions)


def test_qualified_conformer_is_usable_by_sample_poses():
    from dataclasses import replace

    atoms, frags, (ca, cb) = _pair()
    (p,) = sample_poses(
        atoms, frags, replace(ca, status="qualified"), cb, shells=(2.0,), n_samples=1, seed=0
    )
    assert p.status == "placed"


def test_place_normalises_the_direction():
    atoms, frags, (ca, cb) = _pair()
    num_a, num_b = atoms.numbers[list(frags[0].atoms)], atoms.numbers[list(frags[1].atoms)]
    d = np.array([1.0, 2.0, -0.5])
    q = np.array([0.0, 0.0, 0.0, 1.0])
    unit = poses.place(ca.positions, num_a, cb.positions, num_b, d / np.linalg.norm(d), q, 2.0)
    big = poses.place(ca.positions, num_a, cb.positions, num_b, 7.0 * d, q, 2.0)
    assert unit[2] == big[2] == "placed"
    assert np.allclose(unit[0], big[0]) and np.isclose(unit[3], big[3])
