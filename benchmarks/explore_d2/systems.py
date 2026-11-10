"""Per-system constants for run_d2.py (docs/PLAN_HF5_CROSS_SYSTEM_2026-11-08.md §3).

b04 / b05 values are imported from their own runners (run_b04_d2.py, run_b05_d2.py) so
there is one source for each measured number; those runners are not changed. A value that
has not been measured for a system is None, and run_d2.py refuses the option that needs it.
"""

from dataclasses import dataclass
from pathlib import Path

from run_b04_d2 import HF_B04
from run_b05_d2 import FRAGMENT_STATES as FRAGMENT_STATES_B05
from run_b05_d2 import HF_B05, HS_BRACKET_B05

from explore.d2 import HFSettings, HSBracket

ROOT = Path(__file__).resolve().parents[2]


@dataclass(frozen=True)
class System:
    name: str
    input_path: Path
    charge: int
    multiplicity: int
    # keyed by prrs.chemistry._formula of each fragment (ase Hill order)
    fragment_states: dict
    hf: HFSettings | None
    hs_bracket: HSBracket | None
    note: str


SYSTEMS = {
    "b04": System(
        name="b04",
        input_path=ROOT / "benchmarks/b04_sn2_chloride/input.extxyz",
        charge=-1,
        multiplicity=1,
        fragment_states={
            "CH3Cl": {"charge": 0, "multiplicity": 1},
            "Cl": {"charge": -1, "multiplicity": 1},
        },
        hf=HF_B04,
        # HF-5 a2 (DONE_LOG 2026-11-09): sigma_c p90 0.0069 -> B = 0.003 + 3 x 0.0069 =
        # 0.0237 (user 2026-11-09 「可以」); coarse threshold 2e-3 as b05 (b04 sigma_c not larger).
        hs_bracket=HSBracket(
            midpoint_coarse_transverse_fmax_eV_A=2e-3, midpoint_sign_bound_eV_A=0.0237
        ),
        note="Cl- + CH3Cl; HF values: DONE_LOG 2026-10-25 (run_b04_d2.HF_B04)",
    ),
    "b05": System(
        name="b05",
        input_path=ROOT / "benchmarks/b05_tertiary_halide/input.extxyz",
        charge=-1,
        multiplicity=1,
        fragment_states=FRAGMENT_STATES_B05,
        hf=HF_B05,
        hs_bracket=HS_BRACKET_B05,
        note="tBuCl + OH-; consistency entry only, production b05 runs use run_b05_d2.py",
    ),
    "n1": System(
        name="n1",
        input_path=ROOT / "benchmarks/n1_diels_alder/input.extxyz",
        charge=0,
        multiplicity=1,
        fragment_states={
            "C5H6": {"charge": 0, "multiplicity": 1},
            "C3H4O": {"charge": 0, "multiplicity": 1},
        },
        # HF-5 n1 (DONE_LOG 2026-11-09, rf_jitter_scan.py over 8 L0 orientations, bond-changing
        # ridges only, (point, tau) with repeats > branch_rmsd_A 0.05 A dropped, max taken):
        # sigma 0.0027, geometry noise 0.0104 A (user 2026-11-09「可以」). B not measured (n2).
        hf=HFSettings(
            reaction_force_sigma=0.0027,
            geometry_noise_A=0.0104,
            midpoint_transverse_fmax_eV_A=1e-3,
            hs_cap=3000,
            hc_cap=3000,
        ),
        hs_bracket=None,
        note="cyclopentadiene + acrolein; input built by benchmarks/n1_diels_alder/build.py",
    ),
}
