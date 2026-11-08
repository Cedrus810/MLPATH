"""Read-only: compare the 10-18 constrained-scan peak region with the POLAR-1 b04 TS.

Reference: runs/b04_probe2_openceiling/ts_candidates/ts0000.extxyz (PRRS min-mode
following, analytic Hessian order 1; QM-refined and IRC-confirmed in
docs/experiments/b04_qm/README.md). No model is loaded; energies are the stored ones.
The TS is D3h-like, so the two Cl labels are tried both ways and the smaller RMSD kept.
"""

import json
import sys
from pathlib import Path

import numpy as np
from ase.io import read

ROOT = Path(__file__).resolve().parents[2]
TS = ROOT / "runs/b04_probe2_openceiling/ts_candidates/ts0000.extxyz"
C, CL_LEAVE, CL_NUC = 0, 1, 5


def kabsch_rmsd(a, b):
    a = a - a.mean(axis=0)
    b = b - b.mean(axis=0)
    u, _, vt = np.linalg.svd(a.T @ b)
    d = np.sign(np.linalg.det(u @ vt))
    rot = u @ np.diag([1.0, 1.0, d]) @ vt
    return float(np.sqrt(np.mean(np.sum((a @ rot - b) ** 2, axis=1))))


def dist(p, i, j):
    return float(np.linalg.norm(p[i] - p[j]))


def main(run_dir, out_path):
    ts = read(TS)
    p_ts = ts.positions
    swap = list(range(len(ts)))
    swap[CL_LEAVE], swap[CL_NUC] = CL_NUC, CL_LEAVE
    e_ts = float(ts.get_potential_energy())  # extxyz energy= goes to the attached calc
    rows = []
    for name in ("forward", "fine_fwd", "fine_rev"):
        for f in read(Path(run_dir) / f"{name}.extxyz", index=":"):
            q = float(f.info["q_target"])
            if not 2.25 <= q <= 2.31:
                continue
            p = f.positions
            rows.append(
                {
                    "scan": name,
                    "q_target": q,
                    "energy_minus_ts_eV": float(f.info["energy_eV"]) - e_ts,
                    "reaction_force": float(f.info["reaction_force"]),
                    "r_C_Clnuc": dist(p, C, CL_NUC),
                    "r_C_Clleave": dist(p, C, CL_LEAVE),
                    "rmsd_to_ts_A": min(kabsch_rmsd(p, p_ts), kabsch_rmsd(p, p_ts[swap])),
                }
            )
    record = {
        "reference": str(TS.relative_to(ROOT)),
        "ts_energy_eV": e_ts,
        "ts_r_C_Cl": [dist(p_ts, C, CL_LEAVE), dist(p_ts, C, CL_NUC)],
        "points": rows,
    }
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    Path(out_path).write_text(json.dumps(record, indent=2))
    print("TS r(C-Cl)", ["%.4f" % x for x in record["ts_r_C_Cl"]], "E", e_ts)
    print(" scan      q      E-E_TS(meV)  dE/dq   rC-Clnuc rC-Clleave RMSD")
    for r in rows:
        print(
            " %-8s %.3f %+9.2f %+7.3f  %.4f   %.4f   %.4f"
            % (
                r["scan"],
                r["q_target"],
                1e3 * r["energy_minus_ts_eV"],
                r["reaction_force"],
                r["r_C_Clnuc"],
                r["r_C_Clleave"],
                r["rmsd_to_ts_A"],
            )
        )


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
