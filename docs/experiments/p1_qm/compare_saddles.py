"""A_theta[E_m] vs T[E_m] vs T[E_QM] on malonaldehyde. Run from the repo root."""
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]      # repo root, independent of cwd
os.chdir(ROOT)
sys.path.insert(0, str(ROOT / "src"))
from ase.io import read
from prrs.chemistry import automorphisms, symmetric_rmsd
from prrs.config import SearchConfig
from prrs.runner import _wavenumbers, curvature_spectrum

M = "/home/ruigengji/MLP/mace/MACE-OFF24_medium.model"
from prrs import torch_guard
torch_guard.install(trusted=[M])
from mace.calculators import MACECalculator

calc = MACECalculator(model_paths=M, model_type="MACE", device="cuda", default_dtype="float64")
cfg = SearchConfig(quench_fmax_eV_A=0.01, minimum_check_step_A=0.01)

G = Path("docs/experiments/p1_qm")
orca = read(G / "ts_orca_on_maceoff24.xyz")
dft = read(G / "ts.xyz")
prrs = read("runs/p1_round8/ts_candidates/ts0000.extxyz")
perms, _ = automorphisms(orca, 1.2)

def internals(a):
    O = [i for i, z in enumerate(a.numbers) if z == 8]
    H = [i for i, z in enumerate(a.numbers) if z == 1]
    h = min(H, key=lambda i: a.get_distance(i, O[0]) + a.get_distance(i, O[1]))
    return a.get_distance(*O), a.get_distance(h, O[0]), a.get_distance(h, O[1])

print("对称感知 RMSD")
print("  T[E_m] (ORCA RS-P-RFO on MACE-OFF24)  vs  A_theta[E_m] (PRRS min-mode)  "
      f"{symmetric_rmsd(orca, prrs, perms):.4f} A")
print("  T[E_m]                                vs  T[E_QM] (DFT)                "
      f"{symmetric_rmsd(orca, dft, perms):.4f} A")
print("  A_theta[E_m]                          vs  T[E_QM]                      "
      f"{symmetric_rmsd(prrs, dft, perms):.4f} A")
print()
print(f"{'':38s} r(O..O)  r(O1-H)  r(O2-H)")
for name, a in (("T[E_m]      ORCA RS-P-RFO / MACE-OFF24", orca),
                ("A_theta[E_m]  PRRS min-mode / MACE-OFF24", prrs),
                ("T[E_QM]     DFT r2SCAN-3c", dft)):
    x = internals(a)
    print(f"{name:38s} {x[0]:7.4f}  {x[1]:7.4f}  {x[2]:7.4f}")
print()
print("PRRS 自己的解析 Hessian（同一套刚体投影），逐个几何：")
for name, a in (("T[E_m]      ORCA on MACE-OFF24", orca),
                ("A_theta[E_m]  PRRS on MACE-OFF24", prrs)):
    probe = a.copy()
    probe.info["charge"], probe.info["spin"] = 0, 1
    probe.calc = calc
    values, _vec, floor, prov = curvature_spectrum(probe, lambda: calc, cfg, source="analytic")
    neg = [float(v) for v in values if v < -cfg.minimum_check_eigenvalue_tol]
    print(f"  {name:34s} 负本征值 {len(neg)} 个  source={prov.get('source')}")
    print(f"  {'':34s} {[round(v, 6) for v in neg]} eV/(A^2 amu)")
    print(f"  {'':34s} {[round(w, 1) for w in _wavenumbers(neg)]} cm^-1   刚体地板 {floor}")
