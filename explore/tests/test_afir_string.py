"""Simplified string + climbing image (explore.afir_string) on the TwoSaddle bench: nodes
stay evenly spread, the climbing image reaches a bench saddle with E/F only, counts close.
CPU only."""

import numpy as np
from twosaddle import SADDLE_Y, TwoSaddleInvariant, local_coords, two_saddle_atoms

from explore.afir_string import reparametrize, resample, string_climb
from explore.ledger import Ledger


class Batch:
    def __init__(self):
        self.structures = 0

    def __call__(self, frames):
        out = []
        for a in frames:
            self.structures += 1
            b = a.copy()
            b.calc = TwoSaddleInvariant()
            out.append((b.get_potential_energy(), b.get_forces()))
        return out


def _template():
    a = two_saddle_atoms(-1.0, 0.0)
    a.info.update(charge=0, spin=1)
    return a


def test_reparametrize_spreads_nodes_evenly():
    xs = [two_saddle_atoms(x, 0.0).positions for x in (-1.0, -0.95, -0.9, -0.85, 1.0)]
    y = reparametrize(np.stack(xs))
    d = np.linalg.norm((y[1:] - y[:-1]).reshape(len(y) - 1, -1), axis=1)
    assert np.allclose(d, d.mean(), rtol=1e-6)
    assert len(resample(xs, 9)) == 9


def test_string_and_climbing_image_find_the_bench_saddle(tmp_path):
    """Start: all interior nodes crowded next to the reactant (the CI-NEB failure mode);
    reparametrization spreads them, the climbing image converges to a bench saddle."""
    t = _template()
    crowd = [two_saddle_atoms(x, 0.05).positions for x in np.linspace(-1.0, -0.8, 10)]
    nodes = (
        [two_saddle_atoms(-1.0, 0.0).positions] + crowd + [two_saddle_atoms(1.0, 0.0).positions]
    )
    ledger = Ledger(tmp_path / "s.jsonl", arm="T", replay=False)
    ev = Batch()
    out = string_climb(t, nodes, ev, ledger, "s", string_steps=2000, ci_steps=4000)
    assert ledger.spent == ev.structures == out["forwards"]
    assert out["string_status"] == "converged"
    assert out["ci_status"] == "converged", out["ci_status"]
    x, y, _ = local_coords(out["ci_positions"])
    assert abs(x) < 1e-2 and min(abs(y - s) for s in SADDLE_Y) < 1e-2
    xs = sorted(local_coords(n)[0] for n in out["string_nodes"])
    assert xs[0] < -0.9 and xs[-1] > 0.9 and np.diff(xs).max() < 0.5  # spread, not crowded
