"""b05: descend both sides of the three continuation saddles, and decompose their mode.

Task 2 of the 2026-09-09 ruling. `runs/b05_quench_continuation/` left three order-1
stationary points at the elimination product composition (isobutylene + chloride +
water), 1.277 eV below the reactant complex, each with one imaginary mode at -75 to
-77 cm^-1 whose displacement is 99.4 % on the water.

Two things are asked here, and the second exists because the first is not enough:

  * both sides of each saddle, through the existing descend_saddle + confirm_minimum
    protocol, six trajectories and six endpoints kept separately
  * a decomposition of the negative mode into the water's centre-of-mass TRANSLATION,
    its ROTATION about that centre, and its INTERNAL deformation. "99.4 % on the water"
    says the mode is localised there; it does NOT say the water moves rigidly, and
    H15 alone carrying 0.967 of it is a reason to look rather than to conclude

Mass-weighting convention, stated because the decomposition is meaningless without it:
`hessian_spectrum` returns eigenvectors of the MASS-WEIGHTED Hessian, so a component v
corresponds to the Cartesian displacement u_i = v_i / sqrt(m_i). The three subspaces are
built in mass-weighted coordinates, where they are mutually orthogonal, and the reported
shares are of the squared mass-weighted norm -- the same quantity the eigenvector is
normalised in.

Nothing here can make the elimination channel `observed_direct`: a saddle whose mode is
a product-side libration connects two arrangements of the products, not the reactants to
the products.
"""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

from ase.io import read, write  # noqa: E402

from prrs.chemistry import chemical_key  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402
from prrs.runner import (  # noqa: E402
    PathRecorder,
    confirm_minimum,
    descend_saddle,
    hessian_spectrum,
)

CONTINUATION = ROOT / "runs" / "b05_quench_continuation"
# the transferred hydrogen differs per sample; the water is O14 + H15 + that hydrogen
SAMPLES = {"t000002": 4, "t000058": 12, "t000080": 8}
OXYGEN, HYDROXIDE_H = 14, 15


def load_config():
    text = (ROOT / "benchmarks" / "b05_tertiary_halide" / "config.json").read_text()
    return SearchConfig(
        **json.loads(
            "\n".join(line for line in text.splitlines() if not line.lstrip().startswith("//"))
        )
    )


def decompose(atoms, mode, members):
    """Split a mass-weighted mode on `members` into translation, rotation, internal.

    In mass-weighted coordinates the rigid translations of a fragment are spanned by
    sqrt(m_i) e_a, and its rotations about the fragment's centre of mass by
    sqrt(m_i) (e_a x r_i). Both are built here, orthonormalised, and projected out in
    that order; what survives is internal deformation. Shares are fractions of the mode's
    squared mass-weighted norm restricted to the fragment.
    """
    masses = atoms.get_masses()
    positions = atoms.positions
    members = list(members)
    weights = np.sqrt(masses[members])
    centre = np.average(positions[members], axis=0, weights=masses[members])
    offsets = positions[members] - centre

    vector = mode.reshape(-1, 3)[members]  # mass-weighted components on the fragment
    total = float((vector**2).sum())

    basis = []
    for axis in range(3):
        translation = np.zeros((len(members), 3))
        translation[:, axis] = weights
        rotation = np.cross(np.eye(3)[axis], offsets) * weights[:, None]
        for candidate in (translation, rotation):
            flat = candidate.ravel().astype(float)
            for kept in basis:
                flat = flat - np.dot(flat, kept) * kept
            norm = np.linalg.norm(flat)
            if norm > 1e-9:
                basis.append(flat / norm)
    flat_mode = vector.ravel()
    translation_share = 0.0
    rotation_share = 0.0
    # the first three orthonormalised vectors came from translations, the rest rotations;
    # rebuild the split by projecting on each family separately for clarity
    families = {"translation": [], "rotation": []}
    for axis in range(3):
        translation = np.zeros((len(members), 3))
        translation[:, axis] = weights
        families["translation"].append(translation.ravel())
        families["rotation"].append(
            (np.cross(np.eye(3)[axis], offsets) * weights[:, None]).ravel()
        )
    shares = {}
    residual = flat_mode.copy()
    for name in ("translation", "rotation"):
        vectors = []
        for raw in families[name]:
            flat = raw.astype(float)
            for kept in vectors:
                flat = flat - np.dot(flat, kept) * kept
            norm = np.linalg.norm(flat)
            if norm > 1e-9:
                vectors.append(flat / norm)
        component = sum(np.dot(residual, v) * v for v in vectors) if vectors else 0.0
        shares[name] = float(np.dot(component, component)) / total if total else 0.0
        residual = residual - component
    shares["internal"] = float(np.dot(residual, residual)) / total if total else 0.0
    del translation_share, rotation_share, basis
    return shares, total


def main():
    from benchmarks.b04_sn2_chloride.calculator import polar1

    config = load_config()
    out = ROOT / "runs" / "b05_saddle_descent_budget"
    out.mkdir(parents=True, exist_ok=True)
    summary = {}

    source = read(ROOT / "runs" / "b05_probe3_budget" / "chemical" / "c0000_m0000.extxyz")
    source.info["charge"], source.info["spin"] = config.total_charge, config.multiplicity
    source.calc = polar1()
    reference = source.get_potential_energy()
    source_key = chemical_key(
        source,
        config.bond_scale,
        charge_sensitive=True,
        charge=config.total_charge,
        multiplicity=config.multiplicity,
    )

    for trial, transferred in SAMPLES.items():
        saddle = read(CONTINUATION / f"{trial}_endpoint.extxyz")
        saddle.info["charge"], saddle.info["spin"] = config.total_charge, config.multiplicity
        saddle.calc = polar1()
        values, vectors = hessian_spectrum(saddle, polar1, config, return_vectors=True)
        mode = vectors[:, 0]
        water = [OXYGEN, HYDROXIDE_H, transferred]
        shares, water_norm = decompose(saddle, mode, water)
        whole = float((mode**2).sum())
        print(f"\n{'=' * 74}\n{trial}: lowest eigenvalue {values[0]:.5f} eV/(A^2 amu)")
        print(f"  mode weight on the water fragment: {water_norm / whole:.3%}")
        print(
            "  within the water: translation {translation:.1%} · rotation {rotation:.1%}"
            " · internal {internal:.1%}".format(**shares)
        )
        print(
            f"  water geometry at the saddle: r(O-H15) {saddle.get_distance(OXYGEN, HYDROXIDE_H):.4f}"
            f" · r(O-H{transferred}) {saddle.get_distance(OXYGEN, transferred):.4f}"
            f" · angle {saddle.get_angle(HYDROXIDE_H, OXYGEN, transferred):.2f} deg"
        )
        entry = {
            "eigenvalue_eV_A2_amu": float(values[0]),
            "mode_weight_on_water": water_norm / whole,
            "within_water": shares,
            "sides": {},
        }

        for sign in (+1, -1):
            # The house recorder, so the six descents land on disk in the same shape
            # the search itself writes -- one directory per side.
            recorder = PathRecorder(out / f"{trial}_side{sign:+d}", config, "descent")
            # .copy() drops the calculator and descend_saddle reads the saddle's energy on
            # its first line, so the copy gets one back before it goes anywhere.
            start = saddle.copy()
            start.calc = polar1()
            endpoint, checked, relays = descend_saddle(
                start, mode, polar1, config, sign, seed=abs(sign), recorder=recorder
            )
            side = {"relays": relays, "reached": endpoint is not None}
            if endpoint is None:
                side.update(
                    stopped_because=checked.get("reason"),
                    saddle_order=checked.get("saddle_order"),
                )
                print(
                    f"  side {sign:+d}: no minimum -- {side['stopped_because']} "
                    f"(order {side['saddle_order']})"
                )
            else:
                endpoint.calc = polar1()
                ok, endpoint_check = confirm_minimum(endpoint, polar1, config, seed=abs(sign))
                key = chemical_key(
                    endpoint,
                    config.bond_scale,
                    charge_sensitive=True,
                    charge=config.total_charge,
                    multiplicity=config.multiplicity,
                )
                write(str(out / f"{trial}_side{sign:+d}.extxyz"), endpoint)
                side.update(
                    confirmed_minimum=bool(ok),
                    saddle_order=endpoint_check.get("saddle_order"),
                    imaginary=endpoint_check.get("imaginary_wavenumbers_icm"),
                    energy_vs_source_eV=float(endpoint.get_potential_energy() - reference),
                    fragments=key["fragments"],
                    key=key["key"][:16],
                    same_key_as_source=key["key"] == source_key["key"],
                    r_O_H15=float(endpoint.get_distance(OXYGEN, HYDROXIDE_H)),
                    r_O_Htransferred=float(endpoint.get_distance(OXYGEN, transferred)),
                )
                print(
                    f"  side {sign:+d}: confirm_minimum={ok} order={side['saddle_order']} "
                    f"dE={side['energy_vs_source_eV']:+.4f} eV fragments={side['fragments']} "
                    f"key={side['key']} relays={len(relays)}"
                )
            entry["sides"][f"{sign:+d}"] = side
        summary[trial] = entry

    (out / "summary.json").write_text(json.dumps(summary, indent=1, default=str))
    print(f"\nwritten {out}/summary.json")
    return 0


# 5.0 protocol (docs/experiments/b05_budget_rerun_protocol.md): the only two changes
# are these two budgets and the output directory. Nothing else is touched.
_base_load_config = load_config


def load_config():
    from dataclasses import replace

    return replace(_base_load_config(), quench_steps=1500, irc_relay_hops=8)


if __name__ == "__main__":
    sys.exit(main())
