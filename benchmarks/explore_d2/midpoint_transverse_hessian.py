"""Transverse curvature at bisection midpoints: failed versus converged
(docs/PLAN_HF_MIDPOINT_TRANSVERSE_2026-11-02.md; TODO item 9). Offline diagnostic; the
Hessian never enters a search and nothing is written back to a run.

A midpoint relaxation moves only in the transverse space: Cartesian displacements
orthogonal to the six rigid-body motions and to grad q. Per frame, the model's analytic
Hessian (units as in `prrs.runner._analytic_hessian`) is projected onto that space and its
eigenvalues (eV/A^2) are written, plus the same in mass weighting (eV/(A^2 amu), the
units of the 10-30 saddle spectra), and the rigid-only projected spectrum's lowest modes
with their |cos| to grad q. A finite-difference check (three seeded directions, central
difference at H_STEP_A) is recorded, never used to reject.

  F  every recorded failed midpoint: the frame of least transverse force in its last
     segment; the wander (last frame minus that frame) is decomposed on the transverse
     eigenvectors after Kabsch alignment onto the frame (fraction of its transverse norm^2
     in the lowest 1 / 3 / 5).
  S  every recorded converged midpoint: its last frame.

One JSON per frame; an existing file is skipped, so rerunning the same command resumes.
LIMIT > 0 keeps only the first LIMIT frames per group (smoke).

    python midpoint_transverse_hessian.py OUT_DIR LIMIT RUN_DIR [RUN_DIR ...]
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT / "src"), str(ROOT / "explore" / "src")]
sys.path.insert(0, str(Path(__file__).resolve().parent))

from ase import Atoms  # noqa: E402

from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.d2 import HFSettings, _intent  # noqa: E402
from explore.ts_ef import _align, _Constraint  # noqa: E402
from hs_transverse_diag import MIDPOINT_ORIGINS, segments  # noqa: E402
from prrs.runner import _hessian_providers, _trivial_modes  # noqa: E402

H_STEP_A = 0.04
CHECK_DIRECTIONS = 3
SEED = 20261102


def _lines(path):
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.strip():
                yield json.loads(line)


def frames(run_dir):
    m = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    tau = HFSettings(**m["hf"]).midpoint_transverse_fmax_eV_A
    poses = {r["pose_id"]: r for r in _lines(run_dir / "poses.jsonl")}
    rows = {r["intent_id"]: r for r in _lines(run_dir / "population.jsonl")}
    out = {"F": [], "S": []}
    for ev in _lines(run_dir / "journal.jsonl"):
        if ev.get("event") != "settle" or not ev["task"].startswith("HS|"):
            continue
        rsegs = segments((ev["payload"].get("refine") or {}).get("history") or [])
        for j, s in enumerate(rsegs):
            if s[0]["origin"] not in MIDPOINT_ORIGINS:
                continue
            base = {
                "run": run_dir.name,
                "task": ev["task"],
                "segment": j,
                "q": s[0]["q_target"],
            }
            row = rows[ev["task"][3:]]
            ctx = {"row": row, "pose": poses[row["pose_id"]]}
            if s[-1]["fmax_projected"] <= tau:
                out["S"].append(base | ctx | {"frame": len(s) - 1, "seg": s})
            elif s is rsegs[-1] and ev["status"] == "unresolved:transverse_not_converged":
                k = int(np.argmin([f["fmax_projected"] for f in s]))
                out["F"].append(base | ctx | {"frame": k, "seg": s})
    return out


def complement(fixed):
    """Orthonormal basis of the space orthogonal to the columns of ``fixed``."""
    q, _ = np.linalg.qr(fixed)
    projector = np.eye(len(q)) - q @ q.T
    w, v = np.linalg.eigh(projector)
    return v[:, w > 0.5]


def spectrum(hessian, fixed):
    basis = complement(fixed)
    w, v = np.linalg.eigh(basis.T @ hessian @ basis)
    return w, basis @ v


def analyse(item, calc):
    seg = item["seg"]
    frame = seg[item["frame"]]
    pose = item["pose"]
    atoms = Atoms(numbers=pose["numbers"], positions=np.asarray(frame["positions"]))
    atoms.info.update(pose["info"])
    atoms.calc = calc
    positions = atoms.positions.copy()
    n = len(atoms)
    masses = atoms.get_masses()
    grad_q = _Constraint(_intent(item["row"]), item["q"]).gradient(positions).ravel()

    t0 = time.monotonic()
    (provider,) = _hessian_providers(calc)
    scale = (
        float(getattr(provider, "energy_units_to_eV", 1.0))
        / float(getattr(provider, "length_units_to_A", 1.0)) ** 2
    )
    raw = scale * np.asarray(provider.get_hessian(atoms), dtype=float).reshape(3 * n, 3 * n)
    hessian_s = time.monotonic() - t0
    asymmetry = float(np.abs(raw - raw.T).max() / max(np.abs(raw).max(), 1e-30))
    hessian = 0.5 * (raw + raw.T)

    rng = np.random.default_rng(SEED)
    deviations = []
    for _ in range(CHECK_DIRECTIONS):
        d = rng.normal(size=3 * n)
        d /= np.linalg.norm(d)
        f = []
        for sign in (1, -1):
            atoms.positions = positions + sign * H_STEP_A * d.reshape(n, 3)
            f.append(atoms.get_forces().ravel())
        measured = -(f[0] - f[1]) / (2 * H_STEP_A)
        deviations.append(
            float(np.abs(hessian @ d - measured).max() / max(np.abs(measured).max(), 1e-12))
        )
    atoms.positions = positions

    centered = positions - positions.mean(axis=0)
    rigid = []
    for axis in range(3):
        t = np.zeros_like(positions)
        t[:, axis] = 1.0
        rigid += [t.ravel(), np.cross(np.eye(3)[axis], centered).ravel()]
    rigid = np.array(rigid).T
    w_t, v_t = spectrum(hessian, np.column_stack([rigid, grad_q]))
    w_r, v_r = spectrum(hessian, rigid)
    unit_q = grad_q / np.linalg.norm(grad_q)

    inv_root = 1 / np.sqrt(np.repeat(masses, 3))
    hessian_mw = hessian * np.outer(inv_root, inv_root)
    trivial = np.array(_trivial_modes(positions, masses)).reshape(-1, 3 * n).T
    w_mw, _ = spectrum(hessian_mw, np.column_stack([trivial, grad_q * inv_root]))

    record = {k: item[k] for k in ("run", "task", "segment", "frame", "q")} | {
        "group": item["group"],
        "segment_steps": seg[-1]["relax_step"],
        "fproj": frame["fmax_projected"],
        "fmax": frame["fmax_eV_A"],
        "energy_eV": frame["energy_eV"],
        "transverse_eigs": w_t.tolist(),
        "transverse_eigs_mass_weighted": w_mw.tolist(),
        "rigid_only_lowest": w_r[:5].tolist(),
        "rigid_only_lowest_cos_grad_q": np.abs(v_r[:, :5].T @ unit_q).tolist(),
        "fd_check_step_A": H_STEP_A,
        "fd_check_deviations": deviations,
        "hessian_asymmetry": asymmetry,
        "hessian_s": hessian_s,
    }
    if item["group"] == "F":
        # Remove rigid-body drift first: a large rotation is not linear in the basis.
        wander = _align(positions, np.asarray(seg[-1]["positions"])) - positions
        d = wander.ravel()
        coeff = v_t.T @ d
        total = float(coeff @ coeff)
        record |= {
            "wander_norm_A": float(np.linalg.norm(d)),
            "wander_transverse_fraction": total / max(float(d @ d), 1e-30),
            "wander_fraction_lowest": {
                str(k): float((coeff[:k] ** 2).sum() / max(total, 1e-30)) for k in (1, 3, 5)
            },
        }
    return record


def main():
    out_dir, limit = Path(sys.argv[1]), int(sys.argv[2])
    out_dir.mkdir(parents=True, exist_ok=True)
    groups = {"F": [], "S": []}
    for run_dir in map(Path, sys.argv[3:]):
        for g, items in frames(run_dir).items():
            groups[g] += [it | {"group": g} for it in items]
    print(json.dumps({g: len(v) for g, v in groups.items()}), flush=True)
    calc = polar1("cuda", "float32")
    for g, items in groups.items():
        for i, item in enumerate(items[:limit] if limit > 0 else items):
            out = out_dir / f"{g}_{i:03d}.json"
            if out.exists():
                continue
            t0 = time.monotonic()
            rec = analyse(item, calc)
            rec["wall_s"] = time.monotonic() - t0
            out.write_text(json.dumps(rec) + "\n", encoding="utf-8")
            print(
                json.dumps(
                    {
                        "file": out.name,
                        "lowest": rec["transverse_eigs"][:3],
                        "fd": max(rec["fd_check_deviations"]),
                        "wall_s": rec["wall_s"],
                    }
                ),
                flush=True,
            )


if __name__ == "__main__":
    main()
