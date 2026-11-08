import json

from explore.certify import certify
from explore.ledger import Ledger, counted_factory
from prrs.calculators import proton_relay_atoms, proton_relay_factory
from prrs.config import SearchConfig
from twosaddle import SADDLE_Y, local_coords, two_saddle_atoms, two_saddle_factory

CONFIG = SearchConfig(closed_shell_only=False)
FRAC = 0.9


def _run(tmp_path, guess, factory, budget=200000, known=(), name="c", config=CONFIG):
    ledger = Ledger(tmp_path / f"{name}.jsonl", arm="L4", budget=budget)
    guess = guess.copy()
    guess.calc = factory()
    e0 = guess.get_potential_energy()
    cert = certify(
        guess,
        counted_factory(factory, ledger),
        config,
        ledger,
        task=name,
        seed=0,
        ceiling_eV=e0 + 5.0,
        reference_eV=e0,
        known=known,
    )
    ledger.close()
    return cert, ledger


def test_both_saddles_of_one_endpoint_pair_are_separate_c2(tmp_path):
    certs = []
    for k, y in enumerate(SADDLE_Y):
        guess = two_saddle_atoms(0.05, FRAC * y)
        cert, ledger = _run(
            tmp_path,
            guess,
            two_saddle_factory,
            known=[(f"s{j}", c) for j, c in enumerate(certs)],
            name=f"s{k}",
        )
        assert cert.status == "c2", cert.diagnostics
        certs.append(cert)
    a, b = certs
    assert {s["chemical_key"] for s in a.sides} == {s["chemical_key"] for s in b.sides}
    ya = local_coords(a.saddle_positions)[1]
    yb = local_coords(b.saddle_positions)[1]
    assert abs(ya - yb) > 0.5  # different saddles


def test_geometry_not_energy_separates_the_two_saddles(tmp_path):
    wide = SearchConfig(closed_shell_only=False, basin_energy_eV=0.2)
    certs = []
    for k, y in enumerate(SADDLE_Y):
        cert, _ = _run(
            tmp_path,
            two_saddle_atoms(0.05, FRAC * y),
            two_saddle_factory,
            known=[(f"s{j}", c) for j, c in enumerate(certs)],
            name=f"w{k}",
            config=wide,
        )
        assert cert.status == "c2", cert.diagnostics
        certs.append(cert)
    gap = abs(certs[0].saddle_energy_eV - certs[1].saddle_energy_eV)
    assert gap < wide.basin_energy_eV  # energy cannot decide; the RMSD check did


def test_known_entry_reloaded_from_json_still_deduplicates(tmp_path):
    guess = two_saddle_atoms(0.05, FRAC * SADDLE_Y[1])
    first, _ = _run(tmp_path, guess, two_saddle_factory, name="a")
    first.chemical_pair = json.loads(json.dumps(first.chemical_pair))
    assert isinstance(first.chemical_pair, list)
    second, _ = _run(tmp_path, guess, two_saddle_factory, known=[("a", first)], name="b")
    assert second.status == "duplicate_c2"


def test_saddle_refine_records_guess_to_saddle_displacement(tmp_path):
    cert, _ = _run(tmp_path, _guess(), two_saddle_factory, name="disp")
    d = cert.diagnostics["saddle_refine"]["guess_to_saddle"]
    assert d["max_atom_displacement_A"] > 0 and d["rmsd_kabsch_A"] > 0


def test_guess_constraints_are_ignored(tmp_path):
    from ase.constraints import FixAtoms

    guess = _guess()
    guess.set_constraint(FixAtoms(indices=[3]))
    cert, _ = _run(tmp_path, guess, two_saddle_factory, name="fix")
    assert cert.status == "c2", cert.diagnostics
    assert len(guess.constraints) == 1  # caller's atoms untouched


def test_a_minimum_geometry_is_not_certified(tmp_path):
    cert, ledger = _run(tmp_path, two_saddle_atoms(1.0, 0.0), two_saddle_factory, name="min")
    assert cert.status.startswith("unresolved:saddle_"), cert.status
    assert cert.spent > 0 and cert.spent == ledger.spent


def test_the_same_saddle_twice_is_a_duplicate(tmp_path):
    guess = two_saddle_atoms(0.05, FRAC * SADDLE_Y[1])
    first, _ = _run(tmp_path, guess, two_saddle_factory, name="a")
    second, _ = _run(tmp_path, guess, two_saddle_factory, known=[("a", first)], name="b")
    assert second.status == "duplicate_c2" and second.duplicate_of == "a"


def test_proton_hop_is_c2_with_a_bond_change(tmp_path):
    left, right = proton_relay_atoms((0, 0, 0)), proton_relay_atoms((1, 0, 0))
    guess = left.copy()
    guess.positions = 0.5 * (left.positions + right.positions)
    cert, ledger = _run(tmp_path, guess, proton_relay_factory, name="hop")
    assert cert.status == "c2", cert.diagnostics
    assert cert.event["fingerprint"]["broken"] and cert.event["fingerprint"]["formed"]
    assert cert.spent == ledger.spent > 0


def test_relayed_descent_is_not_c2(tmp_path, monkeypatch):
    import prrs.runner as runner

    real = runner.descend_saddle

    def relayed(*args, **kwargs):
        endpoint, checked, relays = real(*args, **kwargs)
        return endpoint, checked, relays + [{"hop": 0, "saddle_order": 1}]

    monkeypatch.setattr(runner, "descend_saddle", relayed)
    cert, _ = _run(tmp_path, two_saddle_atoms(0.05, FRAC * SADDLE_Y[1]), two_saddle_factory)
    assert cert.status == "unresolved:relayed"
    assert any(side["relays"] for side in cert.sides)


def test_both_sides_in_one_basin_is_reported_not_counted(tmp_path, monkeypatch):
    import prrs.state as state

    monkeypatch.setattr(state, "same_basin", lambda a, b, c: True)
    cert, _ = _run(tmp_path, two_saddle_atoms(0.05, FRAC * SADDLE_Y[1]), two_saddle_factory)
    assert cert.status == "unresolved:both_sides_same_basin"


def test_identity_failure_is_unresolved(tmp_path, monkeypatch):
    import prrs.chemistry as chemistry

    def boom(*a, **k):
        raise chemistry.IdentityUnavailable("automorphism_budget")

    monkeypatch.setattr(chemistry, "chemical_key", boom)
    cert, _ = _run(tmp_path, two_saddle_atoms(0.05, FRAC * SADDLE_Y[1]), two_saddle_factory)
    assert cert.status == "unresolved:identity"


def test_budget_exhaustion_keeps_what_was_paid(tmp_path):
    cert, ledger = _run(
        tmp_path, two_saddle_atoms(0.05, FRAC * SADDLE_Y[1]), two_saddle_factory, budget=30
    )
    assert cert.status == "budget_exhausted"
    assert cert.spent == ledger.spent and ledger.spent <= 30


def test_ledger_stages_follow_the_protocol(tmp_path):
    _run(tmp_path, two_saddle_atoms(0.05, FRAC * SADDLE_Y[1]), two_saddle_factory, name="st")
    stages = {json.loads(x)["stage"] for x in (tmp_path / "st.jsonl").read_text().splitlines()}
    assert {"saddle_refine", "hessian", "descent"} <= stages


def _guess():
    return two_saddle_atoms(0.05, FRAC * SADDLE_Y[1])


def test_incomplete_automorphism_group_is_identity_unresolved(tmp_path, monkeypatch):
    import prrs.chemistry as chemistry

    first, _ = _run(tmp_path, _guess(), two_saddle_factory, name="a")
    real = chemistry.automorphisms

    def capped(atoms, *a, **k):
        _, info = real(atoms, *a, **k)
        return [list(range(len(atoms)))], {**info, "reason": "candidate_limit"}

    monkeypatch.setattr(chemistry, "automorphisms", capped)
    cert, _ = _run(tmp_path, _guess(), two_saddle_factory, known=[("a", first)], name="b")
    assert cert.status == "unresolved:identity"
    assert cert.diagnostics["automorphisms"]["reason"] == "candidate_limit"


def test_no_hessian_configuration_spends_nothing(tmp_path):
    cfg = SearchConfig(closed_shell_only=False, minimum_check="probe")
    cert, ledger = _run(tmp_path, _guess(), two_saddle_factory, config=cfg)
    assert cert.status == "unresolved:no_hessian" and cert.spent == 0 == ledger.spent


def test_hessian_without_order_names_its_reason(tmp_path, monkeypatch):
    import prrs.runner as runner

    bad = {"saddle_order": None, "reason": "trivial_mode_floor_above_tolerance"}
    monkeypatch.setattr(runner, "confirm_minimum", lambda *a, **k: (False, dict(bad)))
    cert, _ = _run(tmp_path, _guess(), two_saddle_factory)
    assert cert.status == "unresolved:hessian_trivial_mode_floor_above_tolerance"


def test_order_two_is_rejected(tmp_path, monkeypatch):
    import prrs.runner as runner

    monkeypatch.setattr(runner, "confirm_minimum", lambda *a, **k: (False, {"saddle_order": 2}))
    cert, _ = _run(tmp_path, _guess(), two_saddle_factory)
    assert cert.status == "rejected:order_2"


def test_descent_failure_is_unresolved(tmp_path, monkeypatch):
    import prrs.runner as runner

    monkeypatch.setattr(runner, "descend_saddle", lambda *a, **k: (None, {"reason": "x"}, []))
    cert, _ = _run(tmp_path, _guess(), two_saddle_factory)
    assert cert.status == "unresolved:descent_failed"


def test_budget_exhaustion_inside_descent_keeps_saddle_and_hessian(tmp_path):
    cert, ledger = _run(tmp_path, _guess(), two_saddle_factory, budget=125)
    assert cert.status == "budget_exhausted", (cert.status, ledger.spent)
    assert cert.saddle_positions is not None and cert.hessian is not None
    assert cert.spent == ledger.spent <= 125


def test_full_run_has_no_unattributed_forwards(tmp_path):
    cert, _ = _run(tmp_path, _guess(), two_saddle_factory, name="full")
    assert cert.status == "c2"
    rows = [json.loads(x) for x in (tmp_path / "full.jsonl").read_text().splitlines()]
    assert not [r for r in rows if r.get("task") == "<unattributed>"]


def test_budget_exhaustion_inside_descent_is_in_descent_stage(tmp_path):
    _run(tmp_path, _guess(), two_saddle_factory, budget=125, name="d")
    rows = [json.loads(x) for x in (tmp_path / "d.jsonl").read_text().splitlines()]
    closes = [r for r in rows if r["event"] == "close"]
    assert [r["stage"] for r in closes][-1] == "descent", closes
    assert {r["stage"] for r in closes} >= {"saddle_refine", "hessian", "descent"}


def _bridged_methyl():
    """CH3 ... O: the minimum has three equivalent H; at the "saddle" H1 bonds to both C
    and O, so the saddle graph only swaps H2 and H3. Returns (minimum, saddle A, saddle
    B), B being A with H1 and H2 relabelled -- the same saddle, transferring H2."""
    import numpy as np
    from ase import Atoms

    hs = [(1.03, 0.0, -0.36), (-0.51, 0.89, -0.36), (-0.51, -0.89, -0.36)]
    minimum = Atoms("CHHHO", positions=[(0, 0, 0), *hs, (0, 0, 3.5)])
    a = minimum.copy()
    a.positions[4] = (2.2, 0.0, -0.8)
    a.positions[1] = (1.25, 0.0, -0.45)  # 1.33 A from C, ~1.0 A from O
    b = a.copy()
    b.positions[[1, 2]] = a.positions[[2, 1]]
    rot = np.array([[np.cos(2.1), -np.sin(2.1), 0], [np.sin(2.1), np.cos(2.1), 0], [0, 0, 1]])
    b.positions = b.positions @ rot.T  # and rigidly turned, as a second encounter would be
    return minimum, a, b


def test_pair_symmetries_keep_the_minimum_symmetry_the_saddle_graph_lost():
    from explore.certify import pair_symmetries
    from prrs.chemistry import automorphisms, symmetric_rmsd

    cfg = SearchConfig(closed_shell_only=False)
    minimum, a, b = _bridged_methyl()
    saddle_only, _ = automorphisms(a, cfg.bond_scale, cfg.active_atoms, cfg.automorphism_limit)
    perms, info = pair_symmetries(a, (minimum,), cfg)
    assert info == {}
    assert len(saddle_only) < len(perms)
    assert symmetric_rmsd(a, b, saddle_only, cfg.active_atoms) > cfg.basin_rmsd_A  # old rule
    assert symmetric_rmsd(a, b, perms, cfg.active_atoms) < 1e-6
    # the union keeps every saddle-graph permutation (nothing the old rule merged is lost)
    assert {tuple(map(int, p)) for p in saddle_only} <= {tuple(map(int, p)) for p in perms}
