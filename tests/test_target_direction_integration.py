"""The whole target-direction chain, at search level, on a P2-shaped system.

The unit tests pin mode selection, overlap, refusal and crossing continuity. What they
cannot pin is the chain: topology-switch frame -> tangent -> follow_min_mode -> saddle
registration -> two-sided relayed descent -> channel. Every link of that was added or
changed for P2, and the only search-level test before this one ran on a double well with a
single internal mode, where "followed the target" and "followed the only mode there is"
are indistinguishable.

So the system here is deliberately shaped like 3-oxobutanal was: a stiff graph-changing
coordinate surrounded by much softer decoys. Measured on it, the transfer saddle sits at
lambda = -51.6 while the softest internal mode is 0.42 -- a factor of 120, the same kind of
separation that made round 2 climb eight rotor barriers at -65 to -80 cm^-1 while the
proton transfer sat at -3186.
"""

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from prrs.config import SearchConfig

SPAN, WIDTH = 2.6, 0.4


class ProtonPlusDecoy(Calculator):
    """A shared proton between two heavy atoms, plus a pendant on a very soft spring.

    The transfer changes the bond set -- (0,1) becomes (1,2) -- so it is a reaction and not
    a conformational change, and the two wells have different chemical keys (CH + CH versus
    C + CH2). The pendant contributes internal modes three to four orders of magnitude
    softer than the transfer (0.42 against -51.6 at the ridge), which is what makes a
    target direction necessary rather than decorative.
    """

    implemented_properties = ["energy", "forces"]

    def calculate(
        self, atoms=None, properties=("energy", "forces"), system_changes=all_changes
    ):
        super().calculate(atoms, properties, system_changes)
        positions = atoms.positions
        energy = 0.0
        gradient = np.zeros_like(positions)

        def spring(i, j, k, r0):
            nonlocal energy
            vector = positions[j] - positions[i]
            distance = np.linalg.norm(vector)
            energy += k * (distance - r0) ** 2
            pull = 2 * k * (distance - r0) * vector / distance
            gradient[i] -= pull
            gradient[j] += pull

        spring(0, 2, 30.0, SPAN)  # holds the heavy pair
        r01 = np.linalg.norm(positions[1] - positions[0])
        r12 = np.linalg.norm(positions[2] - positions[1])
        total, difference = r01 + r12, r01 - r12
        shape = (difference / WIDTH) ** 2 - 1
        energy += 20.0 * (total - SPAN) ** 2 + 0.5 * shape**2
        d_total = 40.0 * (total - SPAN)
        d_difference = 0.5 * 2 * shape * 2 * difference / WIDTH**2
        for (i, j), derivative in (
            ((0, 1), d_total + d_difference),
            ((1, 2), d_total - d_difference),
        ):
            vector = positions[j] - positions[i]
            pull = derivative * vector / np.linalg.norm(vector)
            gradient[i] -= pull
            gradient[j] += pull
        spring(2, 3, 0.35, 1.1)  # the soft decoy

        # Distance springs alone leave every transverse direction free on a collinear
        # arrangement -- not soft, FREE, which is a different thing: measured, those modes
        # read 0.0052 against an eigenvalue tolerance of 0.001, so finite-difference
        # truncation flipped one negative and no descent endpoint was ever confirmed a
        # minimum. That is geometry, not a badly chosen spring constant, so it needs a term
        # that actually resists bending. Written in cos(theta) rather than theta because
        # theta = pi is exactly where the angle coordinate's gradient is singular, and this
        # arrangement sits there.
        for centre, left, right, stiffness in ((1, 0, 2, 0.6), (2, 1, 3, 0.6)):

            def bend_energy(pos, centre=centre, left=left, right=right):
                u = pos[left] - pos[centre]
                v = pos[right] - pos[centre]
                cosine = float(u @ v / (np.linalg.norm(u) * np.linalg.norm(v)))
                return 1.0 + cosine

            energy += stiffness * bend_energy(positions)
            step = 1e-6
            for atom in (centre, left, right):
                for axis in range(3):
                    shifted = positions.copy()
                    shifted[atom, axis] += step
                    plus = bend_energy(shifted)
                    shifted[atom, axis] -= 2 * step
                    gradient[atom, axis] += (
                        stiffness * (plus - bend_energy(shifted)) / (2 * step)
                    )
        self.results = {"energy": float(energy), "forces": -gradient}


def factory():
    return ProtonPlusDecoy()


def source_atoms():
    r01 = (SPAN - WIDTH) / 2
    return Atoms("CHCH", positions=[[0.0, 0, 0], [r01, 0, 0], [SPAN, 0, 0], [SPAN + 1.1, 0, 0]])


@pytest.fixture
def protocol():
    return SearchConfig(
        seed=11,
        max_trials=24,
        max_depth=1,
        families=("stretch", "compress"),
        geometry_amplitudes_A=(0.1, 0.3, 0.5, 0.7),
        # Four atoms, so every pair is within a few bonds; the default four-hop floor
        # would leave no direction to probe at all.
        min_pair_hops=1,
        timestep_fs=0.25,
        response_steps=60,
        quench_fmax_eV_A=0.002,
        min_mode_fmax_eV_A=0.002,
        min_mode_max_steps=80,
        quench_steps=400,
        # A synthetic radical pair is not carbon chemistry and the valence gate has
        # nothing meaningful to say about it.
        closed_shell_only=False,
        refinement_steps=2,
    )


def _write_trial(directory, frames, switch_at):
    """A trial directory shaped like the real one, with a topology switch at `switch_at`.

    Built rather than produced by a search run on purpose. Driving this four-atom system
    through the proposal machinery does not work and should not: pulling the proton off C0
    passes through a three-centre geometry where the two fragments are momentarily one, and
    the collateral gate refuses it as `collateral_fragmentation`. That is the gate being
    right -- 3-oxobutanal is a single molecule, so its fragment count never changes, while
    this test system's does. Fixing that would mean designing a fifth atom into the
    potential to bridge the pair, which tests the geometry rather than the chain.

    So the chain is exercised on real code with a real trajectory: `_dividing_frame` reads
    these files exactly as it reads a trial's own.
    """
    from ase.io import Trajectory
    from prrs.io import append_json
    from prrs.state import encode

    directory.mkdir(parents=True, exist_ok=True)
    previous = None
    with Trajectory(str(directory / "response.traj"), "w") as trajectory:
        for index, atoms in enumerate(frames):
            graph = encode(atoms, 1.2)
            events = {"formed": [], "broken": []}
            if previous is not None:
                events = {
                    "formed": [list(e) for e in sorted(graph.edges - previous)],
                    "broken": [list(e) for e in sorted(previous - graph.edges)],
                }
            previous = graph.edges
            phase = "free" if index else "source"
            frame = atoms.copy()
            frame.calc = None
            frame.info.update({"phase": phase, "time_fs": 0.25 * index})
            trajectory.write(frame)
            append_json(
                directory / "observations.jsonl",
                {
                    "frame": index,
                    "phase": phase,
                    "step": index,
                    "events": events,
                    "potential_eV": 0.0,
                    "edges": [list(e) for e in sorted(graph.edges)],
                },
            )
    assert switch_at < len(frames)


def test_the_whole_target_direction_chain_on_a_stiff_reaction(tmp_path, protocol):
    """Topology-switch frame -> tangent -> climb -> two-sided relayed descent -> two states.

    Every link was added or changed for P2. The decoys here are four orders of magnitude
    softer than the reaction, so a climb that ignored the tangent would land on one of them
    -- which is exactly what round 2 did on the real system.
    """
    from prrs.runner import descend_saddle, follow_min_mode
    from prrs.search import _dividing_frame
    from prrs.chemistry import chemical_key

    # A trajectory that walks the proton across, so frame k is where the graph switches.
    span = np.linspace(-WIDTH, WIDTH * 0.6, 9)
    frames = []
    for q in span:
        atoms = source_atoms()
        positions = atoms.positions.copy()
        positions[1, 0] = (SPAN + q) / 2
        atoms.set_positions(positions)
        frames.append(atoms)
    _write_trial(tmp_path / "trial", frames, switch_at=4)

    seed, origin, tangent, source_name = _dividing_frame(tmp_path / "trial")
    assert seed is not None, origin
    assert source_name == "topology_switch_difference", (origin, source_name)
    assert tangent is not None

    # The tangent is dominated by the transferring proton, which is the whole point.
    magnitudes = np.linalg.norm(tangent, axis=1)
    assert magnitudes.argmax() == 1, magnitudes

    saddle, diagnostics = follow_min_mode(
        seed, factory, protocol, direction=tangent, reference=0.0
    )
    assert saddle is not None, diagnostics
    assert diagnostics["saddle_order"] == 1, diagnostics
    assert diagnostics["target_overlap"] > protocol.min_mode_overlap
    # The softest decoy sits at 0.42, the transfer at -51.6. Two orders of magnitude
    # separate "followed the tangent" from "followed whatever was softest".
    assert min(diagnostics["negative_eigenvalues"]) < -10.0, diagnostics

    # Both sides, through the shared relayed descent, must name two different states.
    mode = np.asarray(diagnostics["unstable_mode"])
    landed = {}
    for sign in (1, -1):
        endpoint, checked, relays = descend_saddle(saddle, mode, factory, protocol, sign)
        assert endpoint is not None, (sign, checked)
        assert checked.get("saddle_order") == 0, (sign, checked)
        landed[sign] = chemical_key(endpoint, protocol.bond_scale)["key"]
        for hop in relays:
            # Orientation is claimed or explicitly not claimed; never silently assumed.
            assert "travel_overlap" in hop
    assert landed[1] != landed[-1], landed
    assert set(landed.values()) == {
        chemical_key(frames[0], protocol.bond_scale)["key"],
        chemical_key(frames[-1], protocol.bond_scale)["key"],
    }, landed


def test_without_a_target_the_same_search_prefers_the_soft_decoy(tmp_path, protocol):
    """The control. `follow_min_mode` keeps its old meaning when no direction is supplied,
    and on this system that means the softest 0.46 mode rather than the -51.6 transfer.
    Run directly rather than through the search, because the search always supplies one.
    """
    from prrs.runner import curvature_spectrum, follow_min_mode

    atoms = source_atoms()
    atoms.calc = factory()  # curvature_spectrum reads the gradient off the atoms
    eigenvalues, vectors, _, _ = curvature_spectrum(atoms, factory, protocol, source="fd")
    tolerance = protocol.minimum_check_eigenvalue_tol
    internal = [i for i, v in enumerate(eigenvalues) if abs(v) > tolerance]
    assert abs(eigenvalues[internal[-1]]) > 100 * abs(eigenvalues[internal[0]]), [
        float(eigenvalues[i]) for i in internal
    ]

    _, unaimed = follow_min_mode(atoms, factory, protocol)
    assert unaimed["target_supplied"] is False
    first = (unaimed.get("history") or [{}])[0]
    assert first.get("followed_order") == internal[0], (first, internal)

    masses = atoms.get_masses()
    stiff = internal[-1]
    tangent = vectors[:, stiff].reshape(-1, 3) / np.sqrt(masses[:, None])
    _, aimed = follow_min_mode(atoms, factory, protocol, direction=tangent)
    assert aimed["target_supplied"] is True
    assert (aimed.get("history") or [{}])[0].get("followed_order") == stiff


def test_this_potentials_forces_match_its_energy():
    """The bend terms are differentiated numerically inside the calculator, so the one thing
    that could silently invalidate every assertion above is a force that is not minus this
    potential's own energy gradient. Two other test potentials in this suite had exactly
    that wrong, and it was invisible because their springs started at equilibrium."""
    atoms = source_atoms()
    atoms.calc = factory()
    analytic = atoms.get_forces()
    numeric = np.zeros_like(analytic)
    step = 1e-6
    for atom in range(len(atoms)):
        for axis in range(3):
            for sign in (1, -1):
                shifted = atoms.copy()
                shifted.calc = factory()
                positions = shifted.positions.copy()
                positions[atom, axis] += sign * step
                shifted.set_positions(positions)
                numeric[atom, axis] += -sign * shifted.get_potential_energy() / (2 * step)
    assert np.abs(analytic - numeric).max() < 1e-6, np.abs(analytic - numeric).max()
