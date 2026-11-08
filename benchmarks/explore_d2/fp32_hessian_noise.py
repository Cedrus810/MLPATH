"""fp32 Hessian-eigenvalue and force noise of MACE-POLAR-1-M (plan "D2 fp32 only", Task 4).

  python fp32_hessian_noise.py --out DIR [--only NAME]

Nine structures: the seven soft-mode structures of fp32_hessian_softmodes.py, b04 ts0000
and the first b05 sB C2 saddle labelled "reactive". Per structure, 5 random rigid
rotations (seed 20261023): analytic ``hessian_spectrum`` (lowest 12 INTERNAL eigenvalues,
rigid-body modes identified by overlap with the trivial basis, not by position), forces
rotated back (per-component std), then one ``confirm_minimum`` with default config.
Every forward goes through one Ledger (DIR/ledger.jsonl). hessian_noise.json is written
after each structure; the "derived" block applies the brief's rules (Task 4) and is only
REPORTED -- no PRRS default is touched. Measurement, not a gate.
"""

import argparse
import json
import socket
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT / "src"), str(ROOT / "explore" / "src")]

SEED = 20261023
ROTATIONS = 5
N_MODES = 12
N_HVP_DIRECTIONS = 200
FIRST_STRUCTURE_LIMIT_S = 120.0
OPT = ROOT / "benchmarks/carbene_ambimodal/results/2026-10-08/opt"
B04_TS = ROOT / "runs/b04_probe2_openceiling/ts_candidates/ts0000.extxyz"
B05_JOURNAL = ROOT / "benchmarks/explore_d2/results/2026-10-23/b05_d2_sB/journal.jsonl"
B05_POSES = ROOT / "benchmarks/explore_d2/results/2026-10-23/b05_d2_sB/poses.jsonl"
# name -> (path or None, total_charge); b05 is loaded from the journal
STRUCTURES = {
    "TS-2OMe": (OPT / "TS-2OMe/geometry.extxyz", 0),
    "TS-2OMe-cis": (OPT / "TS-2OMe-cis/geometry.extxyz", 0),
    "Dimethoxycarbene": (OPT / "Dimethoxycarbene/geometry.extxyz", 0),
    "Cyclohexene-carboxaldehyde": (
        OPT / "1-Cyclohexene-1-carboxaldehyde_6O/geometry.extxyz",
        0,
    ),
    "b09-complex": (ROOT / "benchmarks/carbene_ambimodal/probe_b09/input.extxyz", 0),
    "TS-CCl2-C4": (OPT / "TS-CCl2-C4/geometry.extxyz", 0),
    "TS-CCl2": (OPT / "TS-CCl2/geometry.extxyz", 0),
    "b04-ts0000": (B04_TS, -1),
    "b05-sB-E2-saddle": (None, -1),
}
CURRENT = {
    "minimum_check_eigenvalue_tol": 1e-3,
    "minimum_check_step_A": 0.01,
    "analytic_hessian_check": 1e-2,
}


def plain(x):
    if isinstance(x, dict):
        return {str(k): plain(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [plain(v) for v in x]
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, np.generic):
        return x.item()
    return x


def one_sig_fig(x):
    return float(f"{x:.0e}")


def round_up_1sf(x):
    exponent = int(np.floor(np.log10(x)))
    return float(np.ceil(x / 10**exponent - 1e-9) * 10**exponent)


def random_rotation(rng):
    q, r = np.linalg.qr(rng.normal(size=(3, 3)))
    q = q * np.sign(np.diag(r))
    if np.linalg.det(q) < 0:
        q[:, 0] = -q[:, 0]
    return q


def load_b05_saddle():
    """First settle event with status "c2" and payload event label "reactive"."""
    from ase import Atoms

    found = None
    with B05_JOURNAL.open(encoding="utf-8") as handle:
        for line in handle:
            if '"settle"' not in line or '"status": "c2"' not in line:
                continue
            event = json.loads(line)
            payload = event.get("payload") or {}
            label = (payload.get("event") or {}).get("label")
            if (
                event.get("event") == "settle"
                and event.get("status") == "c2"
                and (label == "reactive")
            ):
                found = event
                break
    if found is None:
        raise RuntimeError("no c2 settle event with label 'reactive' in the b05 journal")
    task = found["task"]
    pose = None
    with B05_POSES.open(encoding="utf-8") as handle:
        for line in handle:
            row = json.loads(line)
            if task.startswith("L4|" + row["pose_id"] + "|"):
                pose = row
                break
    if pose is None:
        raise RuntimeError(f"no pose in poses.jsonl matches task {task[:80]}")
    atoms = Atoms(
        numbers=pose["numbers"], positions=found["payload"]["saddle_positions"], pbc=False
    )
    atoms.info.update(pose["info"])
    return atoms, {"task": task, "saddle_energy_eV": found["payload"].get("saddle_energy_eV")}


def load_structure(name, path, charge):
    from ase.io import read

    from prrs.config import SearchConfig
    from prrs.search import _stamped_with_charge

    config = SearchConfig(charge_sensitive=True, total_charge=charge, multiplicity=1)
    source = {}
    if path is None:
        atoms, source = load_b05_saddle()
    else:
        atoms = read(path)
        source = {"path": str(path)}
    return _stamped_with_charge(atoms, config), config, source


def internal_eigenvalues(values, vectors, basis):
    """Drop the rigid-body modes by overlap with the trivial basis; ascending internal rest."""
    overlap = np.linalg.norm(basis @ vectors, axis=0)
    trivial = np.argsort(-overlap)[: len(basis)]
    keep = np.setdiff1d(np.arange(len(values)), trivial)
    return values[keep], float(np.min(overlap[trivial])), float(np.max(overlap[keep]))


def typical_hvp(values, vectors, masses, rng):
    """Median over random unit directions of max|H d| (raw eV/A^2), the reference scale
    used by prrs _analytic_hessian's relative check. Rebuilt from the projected spectrum."""
    root = np.sqrt(np.repeat(masses, 3))
    mass_weighted = (vectors * values) @ vectors.T
    raw = root[:, None] * mass_weighted * root[None, :]
    dirs = rng.normal(size=(N_HVP_DIRECTIONS, raw.shape[0]))
    dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)
    return float(np.median(np.abs(dirs @ raw).max(axis=1)))


def measure(name, atoms, config, source, factory, ledger, rng):
    from prrs.io import snapshot
    from prrs.reliability import GuardedCalculator
    from prrs.runner import _trivial_modes, confirm_minimum, hessian_spectrum

    entry = {"source": source, "n_atoms": len(atoms), "total_charge": config.total_charge}
    start = time.perf_counter()
    masses = atoms.get_masses()
    entry["min_mass_amu"] = float(masses.min())
    centre = np.average(atoms.positions, axis=0, weights=masses)
    eigs, forces_back, rotations_s, overlaps = [], [], [], []
    hvp_scale = None
    with ledger.segment(name, "measure"):
        for r in range(ROTATIONS):
            t0 = time.perf_counter()
            rot = random_rotation(rng)
            rotated = atoms.copy()
            rotated.positions = (atoms.positions - centre) @ rot.T + centre
            snap = snapshot(rotated, GuardedCalculator(factory(), config))
            forces_back.append(np.asarray(snap.get_forces()) @ rot)
            values, vectors = hessian_spectrum(
                snap, factory, config, return_vectors=True, source="analytic"
            )
            basis = _trivial_modes(snap.positions, masses)
            internal, min_trivial, max_internal = internal_eigenvalues(values, vectors, basis)
            overlaps.append([min_trivial, max_internal])
            eigs.append(internal[:N_MODES])
            if hvp_scale is None:
                hvp_scale = typical_hvp(values, vectors, masses, np.random.default_rng(0))
            rotations_s.append(time.perf_counter() - t0)
        eigs = np.array(eigs)
        forces_back = np.array(forces_back)
        entry["lowest12_per_rotation"] = eigs
        entry["sigma_lambda"] = eigs.std(axis=0, ddof=1)
        entry["ptp_lambda"] = np.ptp(eigs, axis=0)
        entry["max_sigma_lambda"] = float(entry["sigma_lambda"].max())
        entry["max_ptp_lambda"] = float(entry["ptp_lambda"].max())
        entry["trivial_overlap_min_vs_internal_max"] = overlaps
        sigma_f = forces_back.std(axis=0, ddof=1)
        entry["sigma_F_eV_A"] = {
            "max": float(sigma_f.max()),
            "median": float(np.median(sigma_f)),
            "components": sigma_f.size,
        }
        entry["typical_hvp_eV_A2"] = hvp_scale
        entry["rotation_wall_s"] = rotations_s
        confirm_atoms = snapshot(atoms.copy(), GuardedCalculator(factory(), config))
        t0 = time.perf_counter()
        try:
            ok, checked = confirm_minimum(confirm_atoms, factory, config, seed=0)
            recheck = checked.get("internal_mode_recheck") or {}
            entry["confirm"] = {
                "ok": bool(ok),
                "reason": checked.get("reason"),
                "saddle_order": checked.get("saddle_order"),
                "imaginary_wavenumbers_icm": checked.get("imaginary_wavenumbers_icm"),
                "lowest_eigenvalues": checked.get("lowest_eigenvalues"),
                "trivial_mode_floor_worst_residual": checked.get(
                    "trivial_mode_floor_worst_residual"
                ),
                "internal_mode_recheck": {
                    k: recheck.get(k)
                    for k in ("performed", "agree", "blocked", "resolution_stable")
                },
                "wall_s": time.perf_counter() - t0,
            }
        except Exception as exc:  # recorded, not hidden
            entry["confirm"] = {"error": repr(exc), "wall_s": time.perf_counter() - t0}
    entry["wall_s"] = time.perf_counter() - start
    return plain(entry)


def derive(structures):
    done = [v for v in structures.values() if "max_sigma_lambda" in v]
    if not done:
        return None
    sigma_lambda_max = max(v["max_sigma_lambda"] for v in done)
    raw3 = 3.0 * sigma_lambda_max
    tol_new = max(CURRENT["minimum_check_eigenvalue_tol"], one_sig_fig(raw3))
    sigma_f = max(v["sigma_F_eV_A"]["max"] for v in done)
    m_min = min(v["min_mass_amu"] for v in done)
    # FD column noise (raw eV/A^2) = sqrt(2) sigma_F / h; the eigenvalue tolerance is
    # mass-weighted (eV/A^2/amu), so the raw noise is divided by the lightest mass.
    h_needed = 3.0 * np.sqrt(2.0) * sigma_f / (m_min * tol_new)
    h_new = max(CURRENT["minimum_check_step_A"], float(h_needed))
    noise_raw = float(np.sqrt(2.0) * sigma_f / h_new)
    ratios = {
        n: noise_raw / v["typical_hvp_eV_A2"]
        for n, v in structures.items()
        if v.get("typical_hvp_eV_A2")
    }
    ratio_max = max(ratios.values())
    check_literal = 3.0 * ratio_max
    confirms = [v["confirm"] for v in done if "reason" in v.get("confirm", {})]
    floor_hits = sum(c["reason"] == "trivial_mode_floor_above_tolerance" for c in confirms)
    recheck_hits = sum(
        c["reason"] in ("internal_mode_recheck_blocked", "internal_mode_recheck_disagrees")
        for c in confirms
    )
    return {
        "note": "reported only; PRRS defaults are not changed",
        "structures_used": [n for n, v in structures.items() if "max_sigma_lambda" in v],
        "current": CURRENT,
        "minimum_check_eigenvalue_tol": {
            "value": tol_new,
            "inputs": {"max_sigma_lambda_over_structures": sigma_lambda_max, "3x": raw3},
            "formula": "max(1e-3, round_1_sig_fig(3 * max over structures and the lowest 12 "
            "internal modes of std_ddof1(lambda over 5 rotations)))",
        },
        "minimum_check_step_A": {
            "value": h_new,
            "value_rounded_up_1sf": round_up_1sf(h_new),
            "inputs": {
                "max_sigma_F_eV_A": sigma_f,
                "min_mass_amu": m_min,
                "tol_new": tol_new,
                "h_needed": float(h_needed),
            },
            "formula": "max(0.01, 3*sqrt(2)*sigma_F_max / (m_min * tol_new)); from "
            "sqrt(2)*sigma_F/(h*m_min) <= tol_new/3 (FD column noise sqrt(2)*sigma_F/h in "
            "eV/A^2, divided by the lightest mass to be in the mass-weighted units of tol)",
        },
        "analytic_hessian_check": {
            "value_literal": check_literal,
            "value_not_tighter_than_current": max(
                CURRENT["analytic_hessian_check"], check_literal
            ),
            "inputs": {
                "fd_column_noise_eV_A2_at_h_new": noise_raw,
                "noise_over_typical_hvp_per_structure": ratios,
                "max_ratio": ratio_max,
            },
            "formula": "3 * max over structures of [sqrt(2)*sigma_F_max/h_new / "
            "typical_hvp]; typical_hvp = median over 200 random unit directions of "
            "max|H d| (raw eV/A^2, the reference of prrs' own relative check)",
        },
        "reported_only_not_changed": {
            "trivial_floor_fraction_triggers": floor_hits,
            "internal_mode_recheck_triggers": recheck_hits,
            "of_confirm_minimum_runs": len(confirms),
            "confirm_reasons": {
                n: v["confirm"].get("reason") for n, v in structures.items() if "confirm" in v
            },
        },
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--only", choices=sorted(STRUCTURES), help="measure one structure only")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    from benchmarks.b04_sn2_chloride.calculator import polar1
    from explore.ledger import Ledger, view_factory

    t0 = time.perf_counter()
    calc = polar1("cuda", "float32")
    ledger = Ledger(out / "ledger.jsonl", "fp32-hessian-noise", replay=False)
    factory = view_factory(calc, ledger)
    rec = {
        "host": socket.gethostname(),
        "seed": SEED,
        "rotations": ROTATIONS,
        "dtype": "float32",
        "load_s": time.perf_counter() - t0,
        "structures": {},
        "derived": None,
    }
    rng = np.random.default_rng(SEED)
    path = out / "hessian_noise.json"
    names = [args.only] if args.only else list(STRUCTURES)
    stop = None
    for k, name in enumerate(names):
        file, charge = STRUCTURES[name]
        atoms, config, source = load_structure(name, file, charge)
        try:
            rec["structures"][name] = measure(name, atoms, config, source, factory, ledger, rng)
        except Exception as exc:  # recorded, not hidden
            rec["structures"][name] = {"error": repr(exc), "source": source}
        rec["forwards_so_far"] = ledger.spent
        path.write_text(json.dumps(rec, indent=2, default=str), encoding="utf-8")
        wall = rec["structures"][name].get("wall_s")
        print(f"{name}: {wall} s, forwards so far {ledger.spent}", flush=True)
        if k == 0 and len(names) > 1 and wall is not None and wall > FIRST_STRUCTURE_LIMIT_S:
            stop = (
                f"first structure {name} took {wall:.0f} s > {FIRST_STRUCTURE_LIMIT_S:.0f} s: "
                "stopping per plan rule; results so far are in hessian_noise.json"
            )
            break
    ledger.close()
    rec["stopped_early"] = stop
    rec["derived"] = derive(rec["structures"])
    path.write_text(json.dumps(plain(rec), indent=2, default=str), encoding="utf-8")
    if rec["derived"] is not None:
        print(json.dumps(plain(rec["derived"]), indent=2, default=str))
    if stop:
        print(stop, file=sys.stderr)
        sys.exit(2)


if __name__ == "__main__":
    main()
