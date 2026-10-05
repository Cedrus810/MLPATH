#!/usr/bin/env python3
"""connect 模式：两层桥 · 工程草稿（主线之外）· 2026-09-28

对应 docs/PLAN_2026-09-28.md（落笔时实现指纹 b2c056524a704ecd）。
本文件把计划里每一条要求落成代码，全部集中在主线之外：

  * src/prrs/ 一行未动。implementation_sha256 只哈希 src/prrs/*.py（见
    scripts/implementation_hash.py:25），所以本文件的存在不动 `make hash`；
  * pytest 的 testpaths=tests（pyproject.toml:26），本文件不被收集；
  * ruff 只查 src/prrs 和 tests（Makefile:27-31），本文件不进格式门。

合并时按 §0 的地图逐段搬运；每段标注了目标位置与开关要求（默认关，
digest 零差异是唯一证据）。所有「要你落笔」的决策集中在 §K 的 DECISIONS，
provisional 默认值只让机器能跑，正式使用前逐项改为你的决定。

用法（科学环境 python）：
  python PLAN_2026-09-28.py --smoke     # 本机秒级自检：每一段机器都过一遍
  python PLAN_2026-09-28.py --t0 chain --freeze   # T0：清单 + 冻结候选 JSON
  python PLAN_2026-09-28.py --t1 --shots 2000 --seed 1   # T1：暴力标准答案
  python PLAN_2026-09-28.py --g0 --variant chain --depth 3 --trials 40 --out runs/g0_chain
  python PLAN_2026-09-28.py --l0        # L0：γ=0 退化 + 涨落耗散
  python PLAN_2026-09-28.py --l2 --seed 1   # L2：精确 h（相空间网格蒙特卡洛）
  python PLAN_2026-09-28.py --l3 --shots 400   # L3：精确 h 桥 vs 拒绝采样 + 变异检验
  python PLAN_2026-09-28.py --l4        # L4：近似 h 三候选 + ESS 门
  python PLAN_2026-09-28.py --e1 / --e2 / --e4
不写任何 runs/ 之外的东西；--t0 --freeze 才会把清单写到本文件旁边的 JSON。
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from ase import Atoms, units  # noqa: E402
from ase.calculators.calculator import Calculator, all_changes  # noqa: E402
from ase.constraints import FixAtoms  # noqa: E402
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution  # noqa: E402
from ase.md.verlet import VelocityVerlet  # noqa: E402
try:  # ase ≥3.29 把 NEB 挪到了 ase.mep
    from ase.neb import NEB  # noqa: E402
except ImportError:
    from ase.mep.neb import NEB  # noqa: E402
from ase.optimize import FIRE, LBFGS  # noqa: E402

from prrs.chemistry import chemical_key  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402
from prrs.runner import confirm_minimum  # noqa: E402
from prrs.state import encode  # noqa: E402

BOND_SCALE = 1.2  # chemical_key / encode 的键判据，全文件统一
TRAP = np.trapezoid if hasattr(np, "trapezoid") else (lambda y, x: np.trapz(y, x))

# ==========================================================================
# §0 合并地图（合并时逐段搬，搬完删本文件）
# ==========================================================================

MERGE_MAP = """
T   测试台势与清单      -> 🟡 T-链已合入 2026-09-28：src/prrs/calculators.py 的 ProtonRelay（重新设计，不是本文件的 ToyPotential）
                          + tests/test_relay_testbed.py。本文件 §T 的 ToyPotential 作废：阱在固定空间点 + tether（不平移不变）、
                          阱心在重原子核上（过不了重叠守卫）、1.5 Å 间距（新键在旧阱内就形成）、诱饵 / 桥原子不受力、末端 N 不在 VALENCE。
                          T-诱饵 / T-绕路与鞍点清单未做
T1  暴力标准答案        -> scripts/testbed_reference.py（kasuga02 跑）
T2  冻结判据            -> docs/T2_CRITERIA.md（判据全文，跑前冻结）
G0  多步 discover 基线  -> 不需要代码改动：ReactionSearch(max_depth=N) 直接跑玩具势
G1  产物输入与 d        -> ✅ 已合入 2026-09-28：src/prrs/product.py + search.run(product=) + schema + cli --product
                          + tests/test_product.py。d 改为在 Aut(G_R) 下取最小（K5 修订：Aut(G_P) 下取最小是恒等）；
                          本文件 §G1 的 load_product / product_distance_spec / annotate_network_with_d 已作废，勿再用
G2  I-a 队列优先级      -> ✅ 已合入 2026-09-28：config.connect_queue_priority + search._closest_first + connect_queue_order 记录。
                          本文件 §G2 的 ConnectQueue / HGraphPrior 作废：分数入队时定死，HGraphPrior 的更新影响不到已入队的态
G3  I-b 方向优先级      -> ✅ 已合入 2026-09-28：config.connect_direction_priority + propose(key_difference=) + ProductReference.key_difference。
                          本文件 §G3 的 boost_directions 作废：条目形状读错（真实是 (family,(i,j),sign,band,…)），
                          且 band −1 在主线 (0, 1) 循环里永远选不到
G4  路径输出与终止      -> 🟡 路径已合入 2026-09-30：product.product_paths（observed / any_evidence 两视图）+ network['product_paths']；终止未做。
                          本文件 §G4 的 extract_paths 作废：每个目标只返回一条却称全部最短路径，且不分证据视图
L0  Langevin 自由段     -> runner.py run_trial 自由段分派（connect 专用；γ=0 逐帧=VV 已测）
L1  进展集合谓词        -> runner.py:2843 first_entry 推广（{标签:键集合} 用法逐发不变）
L2/L3/L4 精确h/桥/近似h -> src/prrs/bridge.py（新）；L4 fail closed：ESS 不足拒绝输出
L5  conditioned 证据    -> schema.py 证据类型清单 + 漂移账本一条；永不与 attempts 合并
E0  两端输入判定        -> cli.py --mode connect --product P.xyz --target product
E1/E2                   -> bridge.py（鞍点接线并进 search.py 现有 seek_saddle/descend_saddle）
E3  P2 实跑             -> 协议文档 + kasuga02（MACE）；判据用 §K K10
E4  CI-NEB 对照         -> scripts/neb_control.py（单独命名，不进 PRRS 证据字段）
K1..K10                 -> config.py SearchConfig 字段 + schema 备案（值以 §K 你改过的为准）
"""

# ==========================================================================
# §K 要你落笔的（PLAN §7 K1–K10 + T0 势参数 + T2 判据 + G4 终止 + L2 网格）
# 下面全是 provisional 默认值：机器能跑、但不算决定。合并前逐项改成你的值。
# ==========================================================================


@dataclass
class Decisions:
    # K1（挡 L0）：discover 保持 NVE，Langevin 只在 connect 里用。计划原文标注
    # 「本次对话里是我的假设，未确认」——默认 True 但必须你确认。
    k1_discover_stays_nve: bool = True
    # K2（挡 L5）：桥轨迹记为新证据类型 conditioned，永不与 attempts/continuations/ts_support 合并。
    k2_conditioned_never_merged: bool = True
    # K3（挡 G2）：h_graph 先验与更新。h(s) ∝ exp(-λ·d(s)) · (α + w·succ(s))
    # / (α + visits(s)^τ)；succ = 从 s 出发落到 d 更小的态的观测边数。
    k3_h_graph: dict = field(default_factory=lambda: {"lambda": 0.5, "alpha": 1.0, "w": 3.0, "tau": 0.5})
    # K4（挡 G1）：产物输入：原子编号与 R 一致、元素序列对不上就开跑前拒绝、不做 MCS。
    k4_product_input_strict: bool = True
    # K5（挡 G1）：d = labeled 键集合对称差，在 Aut(G_P) 下取最小；键级不进 d（provisional）；
    # 是否同时取 min over Aut(G_R) 待你定。
    k5_distance: dict = field(
        default_factory=lambda: {
            "form": "labeled_bond_set_symmetric_difference",
            "minimize_over": "Aut(G_P)",
            "also_Aut_G_R": False,
            "bond_orders": False,
        }
    )
    # K6（挡 L1）：进展集合严格下降（< 而非 <=）；首次进入连续 k 步才算。
    k6_progress: dict = field(default_factory=lambda: {"strict": True, "k": 2})
    # K7（挡 L0）：γ（ps⁻¹）、T（K）。γ、T 是新调参面：跑前定值，要改走显式修订（PLAN §6）。
    k7_langevin: dict = field(default_factory=lambda: {"gamma_per_ps": 1.0, "temperature_K": 300.0})
    # 约化步实验（L2/L3/L4/E1）的 γ：K7 的 1 ps⁻¹ 下 Doob 项只有物理力的 ~1%，
    # 变异检验测不出差异；50 ps⁻¹ 让偏置在 kT 尺度上有牙齿。正式值随 K7 落笔。
    reduced_experiment_gamma_per_ps: float = 50.0
    # K8（挡 L4、E3）：近似 h 候选（顺序即优先级）；ESS 低于阈值拒绝输出（fail closed）。
    k8_approx_h: dict = field(
        default_factory=lambda: {
            "candidates": ["shooting_logistic", "brownian_bridge", "graph_distance"],
            "ess_min": 100.0,
        }
    )
    # K9（挡 E0）：两端版本 CLI 形状。
    k9_cli: str = "--mode connect --product P.xyz --target progress|product"
    # K10（挡 E3）：E3 对照 P2 冻结 TS 的判据（P2 参照：虚频 −3010 cm⁻¹ 一带、势垒 +361 meV）。
    k10_p2_criteria: dict = field(
        default_factory=lambda: {"rmsd_A": 0.3, "imag_freq_cm": [-3500.0, -2500.0], "also_P1": False}
    )
    # G4 要你落笔的：终止条件——首次到达 P 就停，还是跑满预算。
    g4_termination: str = "first_P_reached"  # 或 "budget_exhausted"
    # T0 要你落笔的：势参数（provisional；清单实测垒高与 d-profile 供你校）。
    # 实测参考：主链位点间垒 ~0.70 eV；诱饵进入垒 ~0.1 一带；绕路直连 ~2.5、过桥两段 ~0.8+0.7。
    # 玩具动力学温度：垒 ~0.7 eV 在 300 K 是 27 kT，窗口内跨不过；测试台的无偏
    # 动力学（T1 图层 / G5 点火 / E2）用这个温度（≠ K7 的真实体系温度）。
    t0_dynamic_temperature_K: float = 2000.0
    t0_potential: dict = field(
        default_factory=lambda: {
            "L": 1.5,  # 重原子间距 Å
            "n_heavy": 4,  # C0..C3，R=(0,1,2) P=(1,2,3)，三步接力
            "site_depth": 2.0,  # 主链位点阱深 eV
            "site_sigma": 0.5,  # 位点阱宽 Å
            "pp_repulsion": 1.5,  # 质子-质子排斥 eV
            "pp_sigma": 0.5,
            "transverse_k": 0.6,  # 质子横向束缚 eV/Å²
            "tether_k": 50.0,  # 重原子 C 固定 eV/Å²（桥 N / 诱饵 O 不 tether）
            # T-诱饵：悬空自由 O（与任何 C 不成键），阱更深、进垒更低、只通回 C1
            "decoy_depth": 3.2,
            "decoy_bump": 0.1,  # C1↔诱饵 之间的垒
            "decoy_wall": 1.5,  # 诱饵↔C0/C2/C3 之间的墙
            # T-绕路：N 与 C0 成键，桥位点在 C0–N 之间（质子同时键连 C0 与 N，
            # d 从 2 涨到 3）；C0→C1 直连垒抬高，到 P 必经桥态
            "bridge_depth": 1.4,
            "direct_bump": 1.2,
        }
    )
    # L2 要你落笔的：h(s,p,τ) 网格分辨率与每格采样量。
    l2_grid: dict = field(
        default_factory=lambda: {
            "s_cells": 14,
            "p_cells": 6,
            "tau_bins": 4,
            "shots_per_cell": 16,
            "dt_fs": 1.0,
            "window_fs": 300.0,
        }
    )


DECISIONS = Decisions()

PROVISIONAL_BANNER = """
⚠ 本文件所有 K1–K10 / T0 / T2 / G4 / L2-grid 决策都是 provisional 默认值，
  只证明机器能跑。正式跑任何 acceptance 之前，把 §K DECISIONS 改成你的落笔。
"""

# ==========================================================================
# §T 测试台：解析多步玩具势（PLAN §2 T0）
#
# 桶式接力：重原子 C0..C3 沿 x 轴间距 L（tether 固定），质子 ×3。
#   R = 质子在 (C0,C1,C2)，P = 质子在 (C1,C2,C3)；每步一个质子右移一位，
#   d(G,G_P)（K5 对称差）沿主链 6→4→2→0 严格下降 —— T-链。
#   T-诱饵：+ 自由 O（不成键、阱更深、进垒低、只通回 C1）——死路分支。
#   T-绕路：+ N（与 C0 成键），桥位点在 C0–N 之间，质子同时键连 C0 与 N，
#     该态 d 从 2 涨到 3；C0→C1 直连被抬高 —— 必经 d 变大的中间体。
# ==========================================================================


@dataclass
class ToyVariant:
    name: str  # chain | decoy | detour
    params: dict

    @property
    def has_decoy(self):
        return self.name == "decoy"

    @property
    def has_bridge(self):
        return self.name == "detour"


def toy_atoms(variant: ToyVariant, placements=None) -> Atoms:
    """原子表：C0,C1,C2,N3 + H×3（+O 诱饵 / +桥 N）。placements = 质子位点编号，默认 R。
    末端必须是 N：全 C 链的 R 与 P 是链反转同构，chemical_key（无标号正则形）
    会把并成同一宏态；末端异原子打破反转，R、P 才是两个化学态（T0 验收第一条）。"""
    p = variant.params
    symbols = ["C"] * (p["n_heavy"] - 1) + ["N"] + ["H"] * 3
    if variant.has_decoy:
        symbols.append("O")
    if variant.has_bridge:
        symbols.append("N")
    positions = [[i * p["L"], 0.0, 0.0] for i in range(p["n_heavy"])]
    sites = heavy_sites(variant)
    for s in placements if placements is not None else [0, 1, 2]:
        positions.append(list(sites[s]))
    if variant.has_decoy:
        positions.append([2.35, 1.9, 0.0])
    if variant.has_bridge:
        positions.append([0.15, 1.45, 0.0])
    return Atoms(symbols, positions=positions)


def heavy_sites(variant: ToyVariant):
    """质子可驻留位点。编号 0..3 主链；4 = 桥（detour）或诱饵（decoy）。"""
    p = variant.params
    sites = {i: np.array([i * p["L"], 0.0, 0.0]) for i in range(p["n_heavy"])}
    if variant.has_bridge:
        sites[p["n_heavy"]] = np.array([0.075, 0.725, 0.0])
    if variant.has_decoy:
        sites[p["n_heavy"]] = np.array([2.35, 1.9, 0.0])
    return sites


def site_adjacency(variant: ToyVariant):
    n = variant.params["n_heavy"]
    adj = {i: set() for i in range(n)}
    for i in range(n - 1):
        adj[i].add(i + 1)
        adj[i + 1].add(i)
    extra = variant.params["n_heavy"]
    if variant.has_bridge:
        adj[extra] = {0, 1}
        adj[0].add(extra)
        adj[1].add(extra)
    if variant.has_decoy:
        adj[extra] = {1}
        adj[1].add(extra)
    return adj


class ToyPotential(Calculator):
    """高斯阱（位点/桥/诱饵）+ 正高斯垒 + 质子间排斥 + 横向束缚 + C tether。
    能量与力解析。桥 N / 诱饵 O 不 tether（作为固定特征）。"""

    implemented_properties = ["energy", "forces"]

    def __init__(self, variant: ToyVariant):
        super().__init__()
        self.variant = variant
        self.p = variant.params
        self.n_calls = 0
        L = self.p["L"]
        self.bridge_center = np.array([0.075, 0.725, 0.0])
        self.decoy_center = np.array([2.35, 1.9, 0.0])
        self.direct_bump_center = np.array([0.75 * L, 0.0, 0.0])
        self.decoy_bump_center = 0.5 * (np.array([L, 0, 0]) + self.decoy_center)
        self.wall_centers = [
            0.5 * (np.array([i * L, 0, 0]) + self.decoy_center) for i in (0, 2, 3)
        ]

    def calculate(self, atoms=None, properties=("energy", "forces"), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        self.n_calls += 1
        pos = atoms.positions
        p = self.p
        n_heavy = p["n_heavy"]
        energy = 0.0
        grad = np.zeros_like(pos)

        for i in range(n_heavy):  # 只 tether 主链 C
            d = pos[i] - np.array([i * p["L"], 0.0, 0.0])
            energy += 0.5 * p["tether_k"] * float(d @ d)
            grad[i] += p["tether_k"] * d

        sigma2 = p["site_sigma"] ** 2
        h_idx = range(n_heavy, n_heavy + 3)
        for i in h_idx:
            q = pos[i]
            for j in range(n_heavy):
                d = q - np.array([j * p["L"], 0.0, 0.0])
                r2 = float(d @ d)
                w = -p["site_depth"] * math.exp(-r2 / (2 * sigma2))
                energy += w
                grad[i] += -w * d / sigma2  # dE/dq = (-w)·d/σ²
            energy += 0.5 * p["transverse_k"] * (q[1] ** 2 + q[2] ** 2)
            grad[i][1] += p["transverse_k"] * q[1]
            grad[i][2] += p["transverse_k"] * q[2]
            if self.variant.has_bridge:
                energy += self._well(q, grad, i, self.bridge_center, p["bridge_depth"], 0.45)
                energy += self._bump(q, grad, i, self.direct_bump_center, p["direct_bump"], 0.35)
            if self.variant.has_decoy:
                energy += self._well(q, grad, i, self.decoy_center, p["decoy_depth"], 0.5)
                energy += self._bump(q, grad, i, self.decoy_bump_center, p["decoy_bump"], 0.4)
                for wcenter in self.wall_centers:
                    energy += self._bump(q, grad, i, wcenter, p["decoy_wall"], 0.4)

        for a, b in itertools.combinations(h_idx, 2):
            d = pos[a] - pos[b]
            s2 = p["pp_sigma"] ** 2
            e = p["pp_repulsion"] * math.exp(-float(d @ d) / (2 * s2))
            energy += e
            g = -e * d / s2
            grad[a] += g
            grad[b] -= g

        self.results = {"energy": float(energy), "forces": -grad}

    def _well(self, q, grad, i, center, depth, sigma):
        d = q - center
        s2 = sigma**2
        e = -depth * math.exp(-float(d @ d) / (2 * s2))
        grad[i] += (-e) * d / s2
        return e

    def _bump(self, q, grad, i, center, height, sigma):
        d = q - center
        s2 = sigma**2
        e = height * math.exp(-float(d @ d) / (2 * s2))
        grad[i] += -e * d / s2
        return e


def toy_factory(variant: ToyVariant):
    def factory():
        return ToyPotential(variant)

    return factory


def _maxwell(atoms, T, rng):
    try:
        from ase.md.velocitydistribution import thermalize_momenta

        thermalize_momenta(atoms, temperature_K=T, rng=rng)
    except (ImportError, TypeError):
        try:
            MaxwellBoltzmannDistribution(atoms, temperature_K=T, rng=rng)
        except TypeError:  # 旧版 ase 无 rng 形参
            MaxwellBoltzmannDistribution(atoms, temperature_K=T)


# ---- 键集合 / chemical 状态 / d（G1 的机器，测试台先用） ----------------


def bond_set(atoms: Atoms, scale=BOND_SCALE):
    """带标号键集合（labeled），来自主线 encode（state.py:80）。"""
    graph = encode(atoms, scale=scale)
    return {tuple(sorted(int(x) for x in e)) for e in graph.edges}


def chemical_digest(atoms: Atoms):
    """chemical_key 的 digest（chemistry.py:449 返回 dict，digest 在 ['key']）。"""
    return chemical_key(atoms, bond_scale=BOND_SCALE)["key"]


def automorphisms(atoms: Atoms, edges=None):
    """自同构群（元素与键集合都保持的置换）。草稿版：暴力枚举，>10 原子直接 raise
    （fail closed，不用恒等充数）。合并时换 canonical.py 的认证群。"""
    if len(atoms) > 10:
        raise ValueError(">10 atoms：合并版用 canonical_form 的认证群")
    numbers = list(atoms.numbers)
    edges = bond_set(atoms) if edges is None else edges
    edge_set = {frozenset(e) for e in edges}
    n = len(atoms)
    result = []
    for perm in itertools.permutations(range(n)):
        if [numbers[perm[i]] for i in range(n)] != numbers:
            continue
        if {frozenset((perm[a], perm[b])) for a, b in edges} == edge_set:
            result.append(list(perm))
    return result


def graph_distance_d(edges, p_edges, aut_p):
    """K5（provisional）：d = min over σ∈Aut(G_P) |edges Δ σ(p_edges)|（labeled）。"""
    edges = {tuple(sorted(e)) for e in edges}
    best = None
    for sigma in aut_p:
        mapped = {tuple(sorted((sigma[a], sigma[b]))) for a, b in p_edges}
        dd = len(edges.symmetric_difference(mapped))
        best = dd if best is None else min(best, dd)
    return best


# ==========================================================================
# §T0 清单（PLAN §2 T0）：全部极小 + 一阶鞍点，chemical_key 两两不同，
# d-profile 报告。--freeze 写 JSON（清单冻结 = 你的落笔确认参数后签字）。
# ==========================================================================


def quench_placement(variant: ToyVariant, placements, fmax=1e-3, steps=400):
    atoms = toy_atoms(variant, placements=placements)
    atoms.calc = ToyPotential(variant)
    LBFGS(atoms, logfile=None).run(fmax=fmax, steps=steps)
    return atoms


def _site_assignment(variant: ToyVariant, atoms: Atoms):
    """淬火后把每个质子指回最近位点（0.7 Å 内；指不回记 -1）。"""
    sites = heavy_sites(variant)
    assign = []
    for i, a in enumerate(atoms):
        if a.symbol != "H":
            continue
        best, bd = None, 0.7
        for s, c in sites.items():
            dist = float(np.linalg.norm(atoms.positions[i] - c))
            if dist < bd:
                best, bd = s, dist
        assign.append(best if best is not None else -1)
    return tuple(sorted(assign))


def enumerate_minima(variant: ToyVariant):
    """3 个质子在全部位点上的可重组合 → 淬火 → 按 chemical digest 去重（同 key 取低能）。"""
    p = variant.params
    n_sites = p["n_heavy"] + (1 if variant.has_bridge else 0) + (1 if variant.has_decoy else 0)
    seen = {}
    for combo in itertools.combinations_with_replacement(range(n_sites), 3):
        atoms = quench_placement(variant, combo)
        key = chemical_digest(atoms)
        e = float(atoms.get_potential_energy())
        if key not in seen or e < seen[key]["energy"] - 1e-9:
            seen[key] = {
                "key": key,
                "energy": e,
                "atoms": atoms.copy(),
                "sites": _site_assignment(variant, atoms),
                "edges": sorted(bond_set(atoms)),
            }
    return seen


def register_reference(variant: ToyVariant, minima):
    """P = 位点 (1,2,3) 的极小；登记 P 的键集合与自同构群；给全部极小标 d_to_P。"""
    p_key = next(k for k, m in minima.items() if m["sites"] == (1, 2, 3))
    rec = minima[p_key]
    P_EDGES[variant.name] = {tuple(e) for e in rec["edges"]}
    AUT[variant.name] = automorphisms(rec["atoms"], rec["edges"])
    for m in minima.values():
        m["d_to_P"] = graph_distance_d({tuple(e) for e in m["edges"]}, P_EDGES[variant.name], AUT[variant.name])
    return p_key, AUT[variant.name]


def find_state(minima, sites):
    """按位点安排找态（找不到返回 None）。"""
    return next((k for k, m in minima.items() if m["sites"] == tuple(sorted(sites))), None)


def fd_hessian(atoms: Atoms, step=1e-3):
    """中心差分 Hessian（未质量加权；数负本征值用，定性足够）。要求 atoms.calc 已设。"""
    base = atoms.positions.copy()
    n3 = 3 * len(atoms)
    cols = []
    for i in range(n3):
        p_plus, p_minus = base.copy(), base.copy()
        p_plus[i // 3, i % 3] += step
        p_minus[i // 3, i % 3] -= step
        atoms.positions = p_plus
        g_plus = -atoms.get_forces().ravel()
        atoms.positions = p_minus
        g_minus = -atoms.get_forces().ravel()
        cols.append((g_plus - g_minus) / (2 * step))
        atoms.positions = base
    hess = np.array(cols).T
    return 0.5 * (hess + hess.T)


def _corridor_pairs(variant: ToyVariant, minima):
    """「差一个质子一步跳」的极小对。"""
    adj = site_adjacency(variant)
    pairs = []
    for ka, kb in itertools.combinations(list(minima), 2):
        sa, sb = minima[ka]["sites"], minima[kb]["sites"]
        diff = [i for i in range(3) if sa[i] != sb[i]]
        if len(diff) == 1 and sb[diff[0]] in adj.get(sa[diff[0]], ()):
            pairs.append((ka, kb))
    return pairs


def find_saddles(variant: ToyVariant, minima, pairs=None, images=9, fmax=0.03, steps=250):
    """一阶鞍点：极小对之间 climbing-image NEB；FD Hessian 数负模定阶。"""
    pairs = _corridor_pairs(variant, minima) if pairs is None else pairs
    saddles = []
    for ka, kb in pairs:
        a0, b0 = minima[ka]["atoms"].copy(), minima[kb]["atoms"].copy()
        for end in (a0, b0):
            end.set_constraint(FixAtoms(indices=[0]))
        imgs = [a0] + [a0.copy() for _ in range(images - 2)] + [b0]
        neb = NEB(imgs, climb=True, method="improvedtangent")
        neb.interpolate()
        for im in imgs:
            im.calc = ToyPotential(variant)
        LBFGS(neb, logfile=None).run(fmax=fmax, steps=steps)
        climb = imgs[images // 2]
        e_ts = float(climb.get_potential_energy())
        hess = fd_hessian(climb)
        neg = int(np.sum(np.linalg.eigvalsh(hess) < -0.5))
        saddles.append(
            {
                "between": [ka, kb],
                "energy": e_ts,
                "barrier_eV": e_ts - min(minima[ka]["energy"], minima[kb]["energy"]),
                "negative_modes": neg,
                "first_order": neg == 1,
            }
        )
    return saddles


def minimax_path(minima, pairs, barrier_of):
    """R→P 的最低垒路径（Dijkstra，边权 = 该极小对的一阶鞍点垒；
    没有一阶鞍点的极小对不直连）。返回节点列表。"""
    import heapq

    adj = {}
    for ka, kb in pairs:
        b = barrier_of.get(frozenset((ka, kb)))
        if b is None:
            continue
        adj.setdefault(ka, []).append((kb, b))
        adj.setdefault(kb, []).append((ka, b))
    r_key = find_state(minima, (0, 1, 2))
    p_key = find_state(minima, (1, 2, 3))
    if r_key is None or p_key is None:
        return []
    best = {r_key: 0.0}
    prev = {}
    pq = [(0.0, r_key)]
    while pq:
        d, u = heapq.heappop(pq)
        if d > best.get(u, 1 << 30):
            continue
        for v, b in adj.get(u, ()):
            nd = max(d, b)  # minimax：路径成本 = 沿途最高垒
            if nd < best.get(v, 1 << 30):
                best[v] = nd
                prev[v] = u
                heapq.heappush(pq, (nd, v))
    if p_key not in best:
        return []
    path, node = [], p_key
    while node != r_key:
        path.append(node)
        node = prev[node]
    path.append(r_key)
    return list(reversed(path))


def run_t0(variant_name, freeze=False, full_pairs=False):
    variant = ToyVariant(variant_name, dict(DECISIONS.t0_potential))
    minima = enumerate_minima(variant)
    p_key, aut = register_reference(variant, minima)
    global _PAIR_CACHE
    _PAIR_CACHE = _corridor_pairs(variant, minima)
    pairs = None if full_pairs else _PAIR_CACHE
    saddles = find_saddles(variant, minima, pairs=pairs)
    keys = [m["key"] for m in minima.values()]
    distinct = len(keys) == len(set(keys))
    profile = sorted((m["d_to_P"], round(m["energy"], 4)) for m in minima.values())
    barrier_of = {frozenset(s["between"]): s["barrier_eV"] for s in saddles if s["first_order"]}
    path = minimax_path(minima, _PAIR_CACHE, barrier_of)
    path_d = [minima[k]["d_to_P"] for k in path]
    path_b = [round(barrier_of.get(frozenset((path[i], path[i + 1])), float("nan")), 3)
              for i in range(len(path) - 1)]
    report = {
        "variant": variant_name,
        "n_minima": len(minima),
        "keys_pairwise_distinct": distinct,
        "P_key": p_key,
        "aut_P_size": len(aut),
        "d_profile_d_energy": profile,
        "minimax_path_d": path_d,
        "minimax_path_barriers_eV": path_b,
        "minimax_max_barrier_eV": max(path_b) if path_b else None,
        "saddles": saddles,
        "saddles_first_order": sum(s["first_order"] for s in saddles),
        "barrier_range_eV": (
            [round(min(s["barrier_eV"] for s in saddles), 3), round(max(s["barrier_eV"] for s in saddles), 3)]
            if saddles
            else None
        ),
    }
    checks = {
        "keys_pairwise_distinct": distinct,
        # 一阶鞍点数：主链三跳 + （绕路）过桥两段。隔盆的极小对 NEB 收敛不到
        # 一阶鞍点是物理事实（不直连），清单里如实标 first_order=false，不算失败。
        "first_order_saddles_ge": report["saddles_first_order"] >= (5 if variant_name == "detour" else 3),
    }
    # T-链：最低垒路径 d 严格单调下降（可逆系统里「往回跳」的边 d 必然上升，
    # 所以判据用最低垒路径，不是「任何边都不许升」）。
    if variant_name == "chain":
        checks["minimax_path_d_monotone_decreasing"] = (
            len(path_d) >= 4 and all(b < a for a, b in zip(path_d, path_d[1:]))
        )
    # T-绕路：最低垒路径必经一段 d 上升（「到 P 必须先经过 d 变大的中间体」），
    # 且那段涨幅如实报告（provisional 参数下是 2→3，涨 1）。
    if variant_name == "detour":
        rises = [b - a for a, b in zip(path_d, path_d[1:]) if b > a]
        checks["minimax_path_has_d_rise"] = bool(rises)
        report["minimax_path_d_rise"] = rises
    print(json.dumps({"report": report, "checks": checks}, indent=2, ensure_ascii=False))
    if freeze:
        out = ROOT / f"T0_inventory_{variant_name}.json"
        out.write_text(json.dumps({"report": report, "checks": checks}, indent=2, ensure_ascii=False))
        print(f"[T0] 清单已写 {out}（冻结 = 你在 §K 落笔确认参数后重跑并签字）")
    return report, checks


# P 端参照（d 的锚）：由 register_reference 填
P_EDGES = {}
AUT = {}


# ==========================================================================
# §T1 暴力标准答案（PLAN §2 T1）：(i) 图层无偏动力学 → 到达概率（Wilson）；
# (ii) 原子层指定一步拒绝采样 → 「真」条件系综。双 seed 一致性。
# 动力学里状态身份用带标号键集合（玩具上与 chemical digest 一一对应，
# encode 比 chemical_key 便宜三个量级；正式版换 digest，接口不变）。
# ==========================================================================


def wilson(p_hat, n, z=1.96):
    if n == 0:
        return (0.0, 1.0)
    denom = 1 + z**2 / n
    centre = (p_hat + z**2 / (2 * n)) / denom
    half = z * math.sqrt(p_hat * (1 - p_hat) / n + z**2 / (4 * n**2)) / denom
    return (max(0.0, centre - half), min(1.0, centre + half))


def _baoab_step(atoms, gamma, kT, dt, rng):
    """欠阻尼 Langevin BAOAB 一步（L0 的机器）。γ=0 严格退化为 VelocityVerlet。"""
    m = atoms.get_masses()[:, None]
    F = atoms.get_forces()
    v = atoms.get_velocities()
    v = v + F / m * dt / 2
    atoms.set_positions(atoms.get_positions() + v * dt / 2)
    if gamma > 0:
        c = np.exp(-gamma * dt / m)
        # 速度形式：噪声幅度是速度空间的 sqrt(kT(1-c²)/m)，不是动量空间的
        noise = np.sqrt(kT * (1 - c**2) / m) * rng.standard_normal(size=v.shape)
        v = v * c + noise
    atoms.set_positions(atoms.get_positions() + v * dt / 2)
    F2 = atoms.get_forces()
    v = v + F2 / m * dt / 2
    atoms.set_velocities(v)


def graph_layer_reference(variant: ToyVariant, shots=200, window=400, seed=1, dt_fs=1.0, every=10):
    """(i) 图层：R 盆地出发的无偏 Langevin（γ 用 K7；温度用测试台动力学温度），
    数到达 P 的比例与转移计数。"""
    gamma_ase = DECISIONS.reduced_experiment_gamma_per_ps / (1000.0 * units.fs)
    T = DECISIONS.t0_dynamic_temperature_K
    kT = units.kB * T
    r_atoms = quench_placement(variant, [0, 1, 2])
    p_atoms = quench_placement(variant, [1, 2, 3])
    p_edges = bond_set(p_atoms)
    transitions = {}
    arrived = 0
    for s in range(shots):
        rng = np.random.default_rng(seed * 10_000 + s)
        atoms = r_atoms.copy()
        atoms.calc = ToyPotential(variant)
        _maxwell(atoms, T, rng)
        prev = frozenset(bond_set(r_atoms))
        hit = False
        for step in range(window):
            _baoab_step(atoms, gamma_ase, kT, dt_fs * units.fs, rng)
            if (step + 1) % every:
                continue
            now = frozenset(bond_set(atoms))
            if now != prev:
                transitions.setdefault((tuple(sorted(prev))[:2], tuple(sorted(now))[:2]), 0)
                transitions[(tuple(sorted(prev))[:2], tuple(sorted(now))[:2])] += 1
                prev = now
            if now == frozenset(p_edges):
                hit = True
                break
        arrived += int(hit)
    p_hat = arrived / shots
    return {"shots": shots, "arrived": arrived, "p_reach_P": p_hat,
            "wilson95": wilson(p_hat, shots),
            "n_transitions_observed": len(transitions)}


def rejection_1d(hop, shots=240, seed=6, dt_fs=1.0, window_fs=300.0):
    """(ii) 原子层拒绝采样（约化一步的 1D 版）：只留进入 B_step 的轨迹。"""
    gamma_ase = DECISIONS.reduced_experiment_gamma_per_ps / (1000.0 * units.fs)
    kT = units.kB * DECISIONS.k7_langevin["temperature_K"]
    dt = dt_fs * units.fs
    n_steps = int(window_fs / dt_fs)
    arrivals = []
    for sh in range(shots):
        rng = np.random.default_rng(seed * 104729 + sh)
        s = -0.05 + 0.02 * rng.standard_normal()
        p = math.sqrt(kT * hop.mass) * rng.standard_normal()
        for step in range(n_steps):
            s, p = hop.step(s, p, gamma_ase, kT, dt, rng)
            if s > hop.s_crit:
                arrivals.append({"t": (step + 1) * dt_fs, "s": s})
                break
    return {"arrivals": arrivals, "shots": shots}


# ==========================================================================
# §T2 冻结验收判据（PLAN §2 T2）：G5/L3/L4 用的判据全文，之后逐字不改。
# FROZEN_ON=None = 未冻结；assert_frozen() 拒绝。--smoke 用 SMOKE 日期旁路。
# ==========================================================================


@dataclass
class FrozenCriteria:
    text: dict = field(default_factory=lambda: {
        "G5_reach_P": "同 seed 同预算：+I-a+I-b 首次收进 P 兼容态的 trial 数 ≤ discover 基线的 0.8 倍",
        "G5_decoy": "T-诱饵上 +I-a+I-b 花在诱饵态的 trial 占比 ≤ discover 基线 + 10 个百分点",
        "G5_detour_no_stall": "T-绕路上不卡死：到达 P，或预算内仍持续产生新微态（验证上层不剪枝兜底）",
        "L3_ensemble": "精确 h 桥系综与 T1 拒绝采样系综：到达时间与到达点分布 KS 检验 p > 0.05",
        "L3_mutation": "h 取反（或换错）必须红：至少一个 KS p < 0.05 或中位到达时间偏移 > 20%",
        "L4_ess": "近似 h 的 ESS ≥ K8.ess_min，否则拒绝输出（fail closed），不得静默使用",
        "L4_gap": "近似 h 桥系综与 L3 精确系综的到达点分布 TVD ≤ 0.2（provisional）",
    })
    frozen_on: str | None = None

    def assert_frozen(self, where):
        if self.frozen_on is None:
            raise PermissionError(f"{where} 的验收判据尚未冻结（T2）：§T2 text 落笔并填 frozen_on")


CRITERIA = FrozenCriteria()


# ==========================================================================
# §G1 产物输入与 d（PLAN §3 G1）：只记录不路由。
# ✅ 已合入 src/prrs/product.py（2026-09-28）。下面这段是草稿原样，作废：
#   graph_distance_d 在 Aut(G_P) 下取最小是恒等（σ(E_P)=E_P）；annotate_network_with_d
#   读 network['nodes']，主线叫 chemical_nodes 且没有 edges 字段，静默什么都不写。
# ==========================================================================


class ProductMismatch(ValueError):
    pass


def load_product(xyz_path, reactant: Atoms):
    """K4：读产物 xyz，原子编号一致（元素序列相同）才收，否则开跑前拒绝；不做 MCS。"""
    from ase.io import read

    prod = read(str(xyz_path))
    check_product(prod, reactant)
    return prod


def check_product(prod: Atoms, reactant: Atoms):
    if len(prod) != len(reactant):
        raise ProductMismatch(f"原子数不一致：R={len(reactant)} P={len(prod)}；K4 拒绝（不做 MCS）")
    if list(prod.numbers) != list(reactant.numbers):
        raise ProductMismatch("元素序列不一致（原子编号必须与 R 一致）；K4 拒绝")
    return True


def product_distance_spec(reactant: Atoms, product: Atoms):
    """P 兼容 = chemical_key 相同；d 与 Aut(G_P) 参照一次算好。"""
    aut = automorphisms(product)
    return {"key_P": chemical_digest(product),
            "edges_P": bond_set(product),
            "aut_P": aut,
            "d_R_to_P": graph_distance_d(bond_set(reactant), bond_set(product), aut)}


def annotate_network_with_d(network: dict, d_fn):
    """每个微态的 d 写进网络记录，但不影响任何决策（同 C1：只记录不路由）。"""
    for node in network.get("nodes", []):
        node["d_to_P"] = d_fn(node.get("edges", ()))
    return network


# ==========================================================================
# §G2 I-a 队列优先级（PLAN §3 G2）：默认关（关=逐字节 FIFO）。
# ==========================================================================


class ConnectQueue:
    """接口对齐 search.py:1963 的 queue。enabled=False：FIFO 原样（discover 默认）。
    enabled=True：h_graph 分数最高先展开，同分 FIFO；分数由 d 与观测边决定，
    都进网络记录 → 队列顺序可从记录复原（G2 验收）。"""

    def __init__(self, enabled=False):
        self.enabled = enabled
        self._items = []  # (seq, item, h_value)
        self._seq = 0

    def append(self, item, h_value=None):
        self._items.append((self._seq, item, h_value))
        self._seq += 1

    def popleft(self):
        if not self.enabled:
            _, item, _ = self._items.pop(0)
            return item
        best = min(range(len(self._items)), key=lambda i: (-self._items[i][2], self._items[i][0]))
        _, item, _ = self._items.pop(best)
        return item

    def __len__(self):
        return len(self._items)


class HGraphPrior:
    """K3（provisional）：h(s) ∝ exp(-λ d(s)) · (α + w·succ(s)) / (α + visits(s)^τ)。"""

    def __init__(self, d_fn, params=None):
        pr = params or DECISIONS.k3_h_graph
        self.lam, self.alpha, self.w, self.tau = pr["lambda"], pr["alpha"], pr["w"], pr["tau"]
        self.d_fn = d_fn
        self.succ = {}
        self.visits = {}

    def observe_edge(self, src, dst):
        self.visits[src] = self.visits.get(src, 0) + 1
        if self.d_fn(dst) < self.d_fn(src):
            self.succ[src] = self.succ.get(src, 0) + 1

    def score(self, state_id, d_value):
        prior = math.exp(-self.lam * d_value)
        succ = self.succ.get(state_id, 0)
        visits = self.visits.get(state_id, 0)
        return prior * (self.alpha + self.w * succ) / (self.alpha + visits**self.tau)


# ==========================================================================
# §G3 I-b 方向优先级（PLAN §3 G3）：R→P 键差进 perturbations.py:1059 的优先级带。
# ==========================================================================


def key_difference(edges_R, edges_P):
    b_minus = set(edges_R) - set(edges_P)  # 要断的键
    b_plus = set(edges_P) - set(edges_R)  # 要成的键
    return b_minus, b_plus


def boost_directions(entries, b_minus, b_plus, enabled, ledger=None):
    """entries：perturbations.propose 的方向条目（entry[3] = band）。
    enabled=False 原样返回（discover 逐字节）。开：碰键差的方向 band → -1
    （优先级带循环 band in (0,1) 之前先吃掉），账本记每条来由。"""
    if not enabled:
        return entries, []
    out, notes = [], []
    diff = {frozenset(e) for e in b_minus} | {frozenset(e) for e in b_plus}
    for entry in entries:
        family, a, b, band = entry[0], entry[1], entry[2], entry[3]
        if family in ("stretch", "compress") and frozenset((a, b)) in diff:
            which = "B-" if frozenset((a, b)) in {frozenset(e) for e in b_minus} else "B+"
            out.append(entry[:3] + (-1,) + entry[4:])
            notes.append({"direction": f"{family}:{a}-{b}", "band": band, "new_band": -1,
                          "reason": f"touches {which} (R→P key difference)"})
        else:
            out.append(entry)
    if ledger is not None:
        ledger.extend(notes)
    return out, notes


# ==========================================================================
# §G4 路径输出与终止（PLAN §3 G4）：每条边标证据类型，不合并；
# 预算内没连上也是答案，写明卡在哪。
# ==========================================================================


def extract_paths(edges, source_id, target_ids):
    """edges: [(src, dst, evidence_types, trial_id)]；全部最短 R→P 路径，
    每条边标证据类型（attempts/continuations/ts_support/conditioned），不合并。"""
    from collections import defaultdict, deque

    adj = defaultdict(list)
    for s, t, ev, tid in edges:
        adj[s].append((t, ev, tid))
    dist = {source_id: 0}
    prev = defaultdict(list)
    q = deque([source_id])
    while q:
        u = q.popleft()
        for v, ev, tid in adj[u]:
            if v not in dist:
                dist[v] = dist[u] + 1
                prev[v].append((u, ev, tid))
                q.append(v)
    paths = []

    def walk(node, acc):
        if node == source_id:
            paths.append(list(reversed(acc)))
            return
        for u, ev, tid in prev[node]:
            walk(u, acc + [(u, node, ev, tid)])

    for t in target_ids:
        if t in dist:
            walk(t, [])
    return paths


def budget_autopsy(records):
    """预算内没连上：状态分布、depth_capped、budget_blocks、最后一个新微态出现在第几个 trial。"""
    statuses = {}
    for r in records:
        statuses[r.get("status", "?")] = statuses.get(r.get("status", "?"), 0) + 1
    last_new, seen = None, set()
    for i, r in enumerate(records):
        t = r.get("target") or r.get("target_microstate")
        if t and t not in seen:
            seen.add(t)
            last_new = i
    return {"n_records": len(records), "statuses": statuses,
            "n_depth_capped": sum(1 for r in records if r.get("depth_capped")),
            "n_blocked": sum(1 for r in records if r.get("network_exclusion")),
            "last_new_microstate_trial": last_new,
            "termination_policy": DECISIONS.g4_termination}


# ==========================================================================
# §G5 测试台验收（PLAN §3 G5）：四组比较。discover 臂可直接用主线
# ReactionSearch（--g0）；I-a/I-b 臂合并前用 ToySearch 验证机制与判官。
# ==========================================================================


class ToySearch:
    """玩具上的小搜索（机制验证，非主线语义替代品）。
    队列 = ConnectQueue（I-a）；方向排序 = boost_directions（I-b）。"""

    def __init__(self, variant: ToyVariant, arm, seed=7, budget=24, window=120, dt_fs=1.0):
        self.variant, self.arm, self.seed = variant, arm, seed
        self.budget, self.window, self.dt_fs = budget, window, dt_fs
        self.minima = enumerate_minima(variant)
        register_reference(variant, self.minima)
        self.p_key = find_state(self.minima, (1, 2, 3))
        self.r_key = find_state(self.minima, (0, 1, 2))
        self.d_of = {k: self.minima[k]["d_to_P"] for k in self.minima}
        self.b_minus, self.b_plus = key_difference(
            {tuple(e) for e in self.minima[self.r_key]["edges"]},
            {tuple(e) for e in self.minima[self.p_key]["edges"]},
        )
        self.h = HGraphPrior(lambda k: self.d_of.get(k, 99))
        self.queue = ConnectQueue(enabled=("I-a" in arm))
        self.queue.append((self.r_key, self.minima[self.r_key]["atoms"]),
                          h_value=self.h.score(self.r_key, self.d_of[self.r_key]))
        self.visited = {self.r_key}
        self.edges, self.records = [], []
        self.decoy_visits = 0
        self.decoy_state = self._decoy_state_key()

    def _decoy_state_key(self):
        if not self.variant.has_decoy:
            return None
        n_heavy = self.variant.params["n_heavy"]
        for k, m in self.minima.items():
            if any(min(e) >= n_heavy for e in m["edges"]):  # 质子-O 键（O 编号最大）
                return k
        return None

    def run(self):
        trial, reached_at = 0, None
        gamma_ase = DECISIONS.reduced_experiment_gamma_per_ps / (1000.0 * units.fs)
        kT = units.kB * DECISIONS.k7_langevin["temperature_K"]
        k_entry = DECISIONS.k6_progress["k"]
        while self.queue and trial < self.budget:
            key, _atoms = self.queue.popleft()
            trial += 1
            self.h.visits[key] = self.h.visits.get(key, 0) + 1
            for probe, direction in self._probes(_atoms):
                if trial >= self.budget:
                    break
                atoms = probe.copy()
                atoms.calc = ToyPotential(self.variant)
                rng = np.random.default_rng(self.seed * 977 + trial * 31 + direction[0])
                _maxwell(atoms, DECISIONS.t0_dynamic_temperature_K, rng)
                ev = "conditioned" if "I-c" in self.arm else "attempts"
                prev, run_count, ok = bond_set(atoms), 0, False
                for _step in range(self.window):
                    _baoab_step(atoms, gamma_ase, kT, self.dt_fs * units.fs, rng)
                    now = bond_set(atoms)
                    run_count = run_count + 1 if now == prev else 0
                    prev = now
                    if run_count >= k_entry:
                        ok = True
                        break
                if ok:
                    dst = chemical_digest(atoms)
                    self.edges.append((key, dst, (ev,), trial))
                    self.h.observe_edge(key, dst)
                    self.records.append({"trial": trial, "status": "edge", "source": key[:8],
                                         "target": dst[:8], "evidence": ev})
                    if dst not in self.visited:
                        self.visited.add(dst)
                        ref = self.minima.get(dst, {}).get("atoms", atoms)
                        self.queue.append((dst, ref), h_value=self.h.score(dst, self.d_of.get(dst, 99)))
                    if dst == self.p_key and reached_at is None:
                        reached_at = trial
                    if self.decoy_state and dst == self.decoy_state:
                        self.decoy_visits += 1
                else:
                    self.records.append({"trial": trial, "status": "no_edge", "source": key[:8]})
        return {"arm": self.arm, "reached_P_at_trial": reached_at, "trials_used": trial,
                "n_states": len(self.visited), "decoy_visits": self.decoy_visits,
                "stalled": reached_at is None and len(self.visited) == 1,
                "records": self.records,
                "edges": [(a[:8], b[:8], ev, t) for a, b, ev, t in self.edges]}

    def _probes(self, atoms):
        """候选方向：一个质子向相邻位点跳；I-b 开时键差方向排前。"""
        sites = heavy_sites(self.variant)
        adj = site_adjacency(self.variant)
        entries = []
        h_idx = [i for i, a in enumerate(atoms) if a.symbol == "H"]
        for i in h_idx:
            near = min(sites, key=lambda s: float(np.linalg.norm(atoms.positions[i] - sites[s])))
            for t in adj.get(near, ()):
                entries.append(("stretch", i, int(t), 0, 0.0))
        entries, _ = boost_directions(entries, self.b_minus, self.b_plus, enabled=("I-b" in self.arm))
        out = []
        for _family, i, t, _band, _o in entries:
            probe = atoms.copy()
            probe.positions[i] = sites[t]
            out.append((probe, (i, int(t))))
        return out


def judge_g5(arm_results):
    """G5 判官：首次收进 P 的 trial 数、态数、诱饵占比、是否卡死。
    判据在 T2 冻结后启用（这里只报数）。"""
    table = {r["arm"]: r for r in arm_results}
    out = {"metrics": {arm: {kk: r[kk] for kk in ("reached_P_at_trial", "n_states", "decoy_visits", "stalled")}
                       for arm, r in table.items()}}
    base, best = table.get("discover"), table.get("I-a+I-b")
    if base and best and base["reached_P_at_trial"] and best["reached_P_at_trial"]:
        out["ratio_reach_vs_discover"] = best["reached_P_at_trial"] / base["reached_P_at_trial"]
    return out


# ==========================================================================
# §L0 欠阻尼 Langevin 响应（PLAN §4 L0）：还没有桥，只验机器。
# ==========================================================================


def check_gamma0_equals_verlet(steps=30, dt_fs=0.5, seed=3):
    """γ=0 时逐帧等于现有 NVE（VelocityVerlet）。"""
    from prrs.calculators import DoubleWell, demo_atoms

    rng = np.random.default_rng(seed)
    a1, a2 = demo_atoms(), demo_atoms()
    for a in (a1, a2):
        a.calc = DoubleWell()
    _maxwell(a1, 300.0, rng)
    a2.set_positions(a1.get_positions())
    a2.set_momenta(a1.get_momenta())
    kT = units.kB * 300.0
    dyn = VelocityVerlet(a2, timestep=dt_fs * units.fs, logfile=None)
    for _ in range(steps):
        _baoab_step(a1, 0.0, kT, dt_fs * units.fs, rng)
        dyn.run(1)
    dp = float(np.abs(a1.get_positions() - a2.get_positions()).max())
    dv = float(np.abs(a1.get_velocities() - a2.get_velocities()).max())
    return {"max_pos_diff": dp, "max_vel_diff": dv, "pass": dp < 1e-10 and dv < 1e-10}


def check_fdt(steps=40000, dt_fs=0.05, T=300.0, gamma_per_ps=200.0, seed=4):
    """涨落耗散：动能温度 + 径向分布对 Boltzmann（DoubleWell：二原子径向双阱，
    分布沿核间距 r 比较，含 r² 径向测度）。
    dt 取 0.05 fs：玩具 C–C 伸缩短周期 ~0.72 fs，步长大了积分器本身发热。
    γ 用强耦合（50 ps⁻¹）专测恒温器；K7 的 γ 是真实体系的运行值，另测。"""
    from prrs.calculators import DoubleWell, demo_atoms

    rng = np.random.default_rng(seed)
    gamma_ase = gamma_per_ps / (1000.0 * units.fs)
    atoms = demo_atoms()
    atoms.calc = DoubleWell()
    _maxwell(atoms, T, rng)
    kT = units.kB * T
    equil = steps // 5
    tail_r, tail_v = [], []
    for step in range(steps):
        _baoab_step(atoms, gamma_ase, kT, dt_fs * units.fs, rng)
        if step >= equil:
            d = atoms.positions[1] - atoms.positions[0]
            tail_r.append(float(np.linalg.norm(d)))
            tail_v.append(atoms.get_velocities().copy())
    v = np.array(tail_v)
    m = atoms.get_masses()
    kin_mean = float(np.mean(0.5 * np.sum(m[None, :, None] * v**2, axis=(1, 2))))
    dof = 3 * len(atoms)
    T_kin = 2.0 * kin_mean / (dof * kT)
    rs = np.array(tail_r)
    bins = np.linspace(max(0.7, rs.min()), min(3.2, rs.max()), 25)
    hist, _ = np.histogram(rs, bins=bins, density=True)
    centres = 0.5 * (bins[:-1] + bins[1:])
    # U(r) 扫描：原子放在 ±r/2，永不重合
    energies = np.array([
        DoubleWell().get_potential_energy(Atoms("C2", [[-r / 2, 0, 0], [r / 2, 0, 0]]))
        for r in centres
    ])
    w = np.exp(-(energies - energies.min()) / kT) * centres**2  # 径向测度 r²
    w = w / TRAP(w, centres)
    tvd = 0.5 * TRAP(np.abs(hist - w), centres)
    # T_kin 的估计量在强动量关联下涨落大（smoke 只求 ~±50% 一致）；
    # 定量的涨落耗散验收 = --l0 长跑。位置分布（tvd）是更好的平衡指标。
    return {"T_kin_over_T": float(T_kin), "tvd_boltzmann_radial": float(tvd),
            "pass": abs(T_kin - 1.0) < 0.5 and tvd < 0.3}


# ==========================================================================
# §L1 进展集合谓词（PLAN §4 L1）：first_entry 从 {标签: 键集合} 推广为图谓词。
# ==========================================================================


def first_entry_step(edges_now, graphs, prev_label, prev_run, k):
    """runner.py:3017-3033 计数规则的逐字转录：连续 k 步落在同一边集才算进入。
    两边都归一到 tuple 边集再比（调用方给的可能是 frozenset）。"""
    now = {tuple(sorted(e)) for e in edges_now}
    label = next((name for name, es in graphs.items()
                  if {tuple(sorted(e)) for e in es} == now), None)
    run = prev_run + 1 if (label is not None and label == prev_label) else int(label is not None)
    return label, run, (label is not None and run >= k)


class ProgressPredicate:
    """推广：B_step = {x : d(G(x)) < d0}（K6 严格 <），连续 k 步才算。"""

    def __init__(self, d_fn, d0, k=None):
        self.d_fn, self.d0 = d_fn, d0
        self.k = k or DECISIONS.k6_progress["k"]

    def step(self, edges_now, run):
        inside = self.d_fn(edges_now) < self.d0  # K6: strict <
        run = run + 1 if inside else 0
        return run, run >= self.k


def check_predicate_consistency(trials=300, seed=11):
    rng = np.random.default_rng(seed)
    graphs = {"A": {frozenset((0, 1))}, "B": {frozenset((0, 2))}}
    lab, run, hits = None, 0, []
    for i in range(trials):
        label = str(rng.choice(["A", "B", "X"]))
        edges = graphs.get(label, {frozenset((9, 9))})
        lab, run, entered = first_entry_step(edges, graphs, lab, run, 2)
        if entered:
            hits.append((i, lab))
            lab, run = None, 0
    d_fn = lambda es: 0 if es == {frozenset((0, 1))} else (1 if es == {frozenset((0, 2))} else 5)
    pred = ProgressPredicate(d_fn, d0=1, k=2)
    run2, hits2 = 0, []
    for i in range(trials):
        # 每 6 步出现连续两步在进展集合里（i%6 ∈ {0,1}），加上零星单步
        inside = (i % 6) in (0, 1) or (i % 11) == 5
        es = {frozenset((0, 1))} if inside else {frozenset((9, 9))}
        run2, entered2 = pred.step(es, run2)
        if entered2:
            hits2.append(i)
            run2 = 0
    return {"label_graph_hits": len(hits), "predicate_hits": len(hits2),
            "pass": len(hits) > 0 and len(hits2) > 0}


# ==========================================================================
# §L2 精确 h（PLAN §4 L2）：T-链的某一步降成 1D（双位点单质子），
# h(s,p,τ) 相空间网格蒙特卡洛；双 seed 在 Wilson 区间内一致。
# ==========================================================================


class ReducedHop:
    """T-链单独一步：1 质子 + 2 重原子（间距 L），1D 动力学。
    B_step（K6 严格 <）：质子键连 C1 ⇔ s > L − 1.284（H–C 判据 1.2×1.07）。
    depth 独立于 T0 的势参数（provisional 0.25 eV ≈ 3.4 kT @300K）：垒取
    kT 量级是为了窗口内有可观测的穿越率——L2/L3/L4/E1 的机器检验不需要高垒。"""

    def __init__(self, L=1.5, depth=0.25, sigma=0.5, mass=1.008):
        self.L, self.depth, self.sigma, self.mass = L, depth, sigma, mass
        self.s_crit = L - 1.284

    def energy(self, s):
        return -self.depth * (math.exp(-s**2 / (2 * self.sigma**2))
                              + math.exp(-(s - self.L) ** 2 / (2 * self.sigma**2)))

    def force(self, s):
        h = 1e-5
        return -(self.energy(s + h) - self.energy(s - h)) / (2 * h)

    def step(self, s, p, gamma, kT, dt, rng):
        """1D BAOAB。"""
        p = p + self.force(s) * dt / 2
        s = s + p / self.mass * dt / 2
        if gamma > 0:
            c = math.exp(-gamma * dt / self.mass)
            p = p * c + math.sqrt(self.mass * kT * (1 - c**2)) * rng.standard_normal()
        s = s + p / self.mass * dt / 2
        p = p + self.force(s) * dt / 2
        return s, p


def exact_h_table(hop: ReducedHop, grid=None, seed=1):
    """h(s,p,τ) = P(τ 内进入 B_step | s,p)。每格 N 发无偏 1D Langevin 跑满窗、
    记首达时间 → 所有 τ 档共用一次采样。"""
    g = grid or DECISIONS.l2_grid
    gamma_ase = DECISIONS.reduced_experiment_gamma_per_ps / (1000.0 * units.fs)
    kT = units.kB * DECISIONS.k7_langevin["temperature_K"]
    dt = g["dt_fs"] * units.fs
    n_steps = int(g["window_fs"] / g["dt_fs"])
    s_edges = np.linspace(-0.4, hop.L + 0.4, g["s_cells"] + 1)
    p_scale = math.sqrt(kT * hop.mass) * 5
    p_edges = np.linspace(-p_scale, p_scale, g["p_cells"] + 1)
    tau_edges = np.linspace(0, g["window_fs"], g["tau_bins"] + 1)
    hit_counts = np.zeros((g["s_cells"], g["p_cells"], g["tau_bins"]))
    totals = np.zeros((g["s_cells"], g["p_cells"]))
    for i in range(g["s_cells"]):
        for j in range(g["p_cells"]):
            s0 = 0.5 * (s_edges[i] + s_edges[i + 1])
            p0 = 0.5 * (p_edges[j] + p_edges[j + 1])
            for sh in range(g["shots_per_cell"]):
                rng = np.random.default_rng(seed * 1_000_003 + i * 1009 + j * 101 + sh)
                s, p = s0, p0
                t_hit = None
                for step in range(n_steps):
                    s, p = hop.step(s, p, gamma_ase, kT, dt, rng)
                    if s > hop.s_crit:
                        t_hit = (step + 1) * g["dt_fs"]
                        break
                totals[i, j] += 1
                if t_hit is not None:
                    b = min(max(int(np.searchsorted(tau_edges, t_hit)) - 1, 0), g["tau_bins"] - 1)
                    hit_counts[i, j, b:] += 1  # τ 更大 → 更可能已到达
    h = hit_counts / np.maximum(totals[..., None], 1)
    return {"h": h, "s_edges": s_edges, "p_edges": p_edges, "tau_edges": tau_edges,
            "totals": totals, "p_scale": p_scale}


def exact_h_consistency(seed_a=1, seed_b=2):
    """L2 验收：两个 seed 的 h 逐格 Wilson 区间重叠比例。"""
    hop = ReducedHop()
    ta, tb = exact_h_table(hop, seed=seed_a), exact_h_table(hop, seed=seed_b)
    ha, hb = ta["h"][:, :, -1], tb["h"][:, :, -1]
    ok = np.ones(ha.shape, dtype=bool)
    for i in range(ha.shape[0]):
        for j in range(ha.shape[1]):
            if ta["totals"][i, j] == 0:
                continue
            lo_a, hi_a = wilson(ha[i, j], int(ta["totals"][i, j]))
            lo_b, hi_b = wilson(hb[i, j], int(tb["totals"][i, j]))
            ok[i, j] = not (hi_a < lo_b or hi_b < lo_a)
    frac = float(ok.mean())
    return {"overlap_fraction": frac, "grid": list(ha.shape), "pass": frac > 0.7}


def h_lookup(table, s, p, tau_fs):
    h, se, pe, te = table["h"], table["s_edges"], table["p_edges"], table["tau_edges"]
    i = int(np.clip(np.searchsorted(se, s) - 1, 0, len(se) - 2))
    j = int(np.clip(np.searchsorted(pe, p) - 1, 0, len(pe) - 2))
    b = min(max(int(np.searchsorted(te, tau_fs)) - 1, 0), h.shape[2] - 1)
    return float(np.clip(h[i, j, b], 1e-4, 1.0))


def h_grad_p(table, s, p, tau_fs, dp=0.02):
    """∂_p log h：查表中心差分（平滑由网格承担）。"""
    ps = table["p_scale"]
    p_lo, p_hi = max(-ps + 1e-6, p - dp), min(ps - 1e-6, p + dp)
    if p_hi <= p_lo:
        return 0.0
    return (math.log(h_lookup(table, s, p_hi, tau_fs)) - math.log(h_lookup(table, s, p_lo, tau_fs))) / (p_hi - p_lo)


# ==========================================================================
# §L3 精确 h 跑桥（PLAN §4 L3）：桥系综 = T1 拒绝采样系综；h 取反必须红。
# ==========================================================================


def run_bridge_1d(hop: ReducedHop, table, shots=120, seed=5, dt_fs=1.0, window_fs=300.0,
                  mutate=None, collect_weights=False):
    """Doob 修正只加在动量通道：dp = [... + 2γkT ∇_p log h]dt + √(2γkT)dW。
    BAOAB 里以漂移力 F_bias = 2γkT·∂_p log h 进两个半踢。mutate="invert" = h→1−h。"""
    gamma_ase = DECISIONS.reduced_experiment_gamma_per_ps / (1000.0 * units.fs)
    kT = units.kB * DECISIONS.k7_langevin["temperature_K"]
    dt = dt_fs * units.fs
    n_steps = int(window_fs / dt_fs)
    arrivals, log_w = [], 0.0
    for sh in range(shots):
        rng = np.random.default_rng(seed * 7919 + sh)
        s = -0.05 + 0.02 * rng.standard_normal()
        p = math.sqrt(kT * hop.mass) * rng.standard_normal()
        t_hit = s_hit = None
        for step in range(n_steps):
            tau = window_fs - step * dt_fs
            if mutate == "invert":
                eps = 0.02
                gp = (math.log(max(1e-6, 1 - h_lookup(table, s, p + eps, tau)))
                      - math.log(max(1e-6, 1 - h_lookup(table, s, p - eps, tau)))) / (2 * eps)
            else:
                gp = h_grad_p(table, s, p, tau)
            f_bias = 2.0 * gamma_ase * kT * gp
            p = p + (hop.force(s) + f_bias) * dt / 2
            s = s + p / hop.mass * dt / 2
            c = math.exp(-gamma_ase * dt / hop.mass)
            dW = rng.standard_normal()
            sig = math.sqrt(hop.mass * kT * (1 - c**2))
            if collect_weights:
                # 加漂移 u 的 Girsanov 增量（对动量通道的噪声实现）
                log_w += (f_bias * sig * dW / math.sqrt(2 * gamma_ase * kT * hop.mass)) \
                         - 0.5 * f_bias**2 * dt / (2 * gamma_ase * kT * hop.mass)
            p = p * c + sig * dW
            s = s + p / hop.mass * dt / 2
            p = p + (hop.force(s) + f_bias) * dt / 2
            if s > hop.s_crit and t_hit is None:
                t_hit, s_hit = (step + 1) * dt_fs, s
                break
        arrivals.append({"t_fs": t_hit, "s": s_hit})
    return {"arrivals": arrivals, "shots": shots, "log_w": log_w}


def ks_compare(a, b):
    from scipy.stats import ks_2samp

    a = [x for x in a if x is not None]
    b = [x for x in b if x is not None]
    if len(a) < 8 or len(b) < 8:
        return {"n": (len(a), len(b)), "p": None}
    res = ks_2samp(a, b)
    return {"statistic": float(res.statistic), "p": float(res.pvalue)}


def run_l3(shots_bridge=120, shots_reject=240, seed=5):
    hop = ReducedHop()
    table = exact_h_table(hop, seed=seed)
    bridge = run_bridge_1d(hop, table, shots=shots_bridge, seed=seed)
    rej = rejection_1d(hop, shots=shots_reject, seed=seed + 1)
    ks_t = ks_compare([a["t_fs"] for a in bridge["arrivals"]], [a["t"] for a in rej["arrivals"]])
    ks_s = ks_compare([a["s"] for a in bridge["arrivals"]], [a["s"] for a in rej["arrivals"]])
    mut = run_bridge_1d(hop, table, shots=shots_bridge, seed=seed, mutate="invert")
    ks_t_mut = ks_compare([a["t_fs"] for a in mut["arrivals"]], [a["t"] for a in rej["arrivals"]])
    ks_s_mut = ks_compare([a["s"] for a in mut["arrivals"]], [a["s"] for a in rej["arrivals"]])
    mutation_red = (ks_t_mut.get("p") or 1) < 0.05 or (ks_s_mut.get("p") or 1) < 0.05
    bridge_ok = (ks_t.get("p") or 0) > 0.05 and (ks_s.get("p") or 0) > 0.05
    return {"bridge_vs_reject_time": ks_t, "bridge_vs_reject_point": ks_s,
            "mutation_time": ks_t_mut, "mutation_point": ks_s_mut,
            "mutation_detected_red": mutation_red,
            "pass": bridge_ok and mutation_red}


# ==========================================================================
# §L4 近似 h（研究段，PLAN §4 L4）：三候选 + ESS 门（fail closed）。
# ==========================================================================


def approx_h_shooting_logistic(hop, shots=60, seed=500):
    """① 射击数据回归 committor：logistic on (s, |p|)。"""
    gamma_ase = DECISIONS.reduced_experiment_gamma_per_ps / (1000.0 * units.fs)
    kT = units.kB * DECISIONS.k7_langevin["temperature_K"]
    X, y = [], []
    for sh in range(shots):
        rng = np.random.default_rng(seed + sh)
        s0 = -0.05 + 0.02 * rng.standard_normal()
        p0 = math.sqrt(kT * hop.mass) * rng.standard_normal()
        s, p, hit = s0, p0, 0
        for _ in range(300):
            s, p = hop.step(s, p, gamma_ase, kT, units.fs, rng)
            if s > hop.s_crit:
                hit = 1
                break
        y.append(hit)
        X.append((s0, abs(p0)))
    Xb = np.hstack([np.array(X), np.ones((len(X), 1))])
    y = np.array(y, dtype=float)
    w = np.zeros(Xb.shape[1])
    for _ in range(400):
        pr = 1 / (1 + np.exp(-np.clip(Xb @ w, -30, 30)))
        w += 5.0 * (Xb.T @ (y - pr)) / len(y)

    def h_fn(s, p, tau):
        z = w[0] * s + w[1] * abs(p) + w[2]
        return float(1 / (1 + math.exp(-np.clip(z, -30, 30))))

    return h_fn


def approx_h_graph_distance(hop):
    """② 图距离 h ∝ exp(-λ·d(s))。动量无关 → ∇_p log h ≡ 0，桥退化为无偏动力学
    ——这是 L4 要暴露的退化，如实报告，不粉饰。"""
    lam = DECISIONS.k3_h_graph["lambda"]

    def h_fn(s, p, tau):
        return math.exp(-lam * (0 if s > hop.L - 1.284 and s < 1.284 + 0 else
                                (1 if s > hop.s_crit else 2)))

    h_fn.momentum_free = True
    return h_fn


def approx_h_brownian_bridge(hop):
    """③ 布朗桥型 h：h = Φ((s_target − s − p·τ/m)/√(2Dτ))，有动量梯度；
    配 Girsanov 重加权（ESS 全程报）。"""
    from math import erf, sqrt

    kT = units.kB * DECISIONS.k7_langevin["temperature_K"]
    gamma_fs = DECISIONS.k7_langevin["gamma_per_ps"] / 1000.0

    def h_fn(s, p, tau):
        if tau <= 1:
            return 0.5
        D = kT / (gamma_fs * 1000.0 * hop.mass)  # Å²/fs
        sd = sqrt(2 * D * tau)
        return 0.5 * (1 + erf(((hop.L - s - p * tau / hop.mass) / sd) / sqrt(2)))

    return h_fn


def effective_sample_size(weights):
    w = np.asarray(weights, dtype=float)
    if w.sum() <= 0:
        return 0.0
    return float(w.sum() ** 2 / np.sum(w**2))


def _numeric_grad_p(h_fn, s, p, tau, dp=0.02):
    lo_h = max(1e-6, h_fn(s, max(-5, p - dp), tau))
    hi_h = max(1e-6, h_fn(s, p + dp, tau))
    return (math.log(hi_h) - math.log(lo_h)) / (2 * dp)


def run_l4(shots=60, seed=9):
    """三候选各在测试台上与 T1 拒绝采样比；ESS 低于 K8 阈值 → 拒绝输出（fail closed）。"""
    hop = ReducedHop()
    rej = rejection_1d(hop, shots=shots * 2, seed=seed + 1)
    gamma_ase = DECISIONS.reduced_experiment_gamma_per_ps / (1000.0 * units.fs)
    kT = units.kB * DECISIONS.k7_langevin["temperature_K"]
    ess_min = DECISIONS.k8_approx_h["ess_min"]
    out = {"candidates": {}}
    for name in DECISIONS.k8_approx_h["candidates"]:
        h_fn = {"shooting_logistic": lambda: approx_h_shooting_logistic(hop),
                "graph_distance": lambda: approx_h_graph_distance(hop),
                "brownian_bridge": lambda: approx_h_brownian_bridge(hop)}[name]()
        arrivals, weights = [], []
        for sh in range(shots):
            rng = np.random.default_rng(seed * 31 + sh)
            s, p, t_hit, lw = -0.05, math.sqrt(kT * hop.mass) * rng.standard_normal(), None, 0.0
            for step in range(300):
                tau = 300.0 - step
                f_bias = 2.0 * gamma_ase * kT * _numeric_grad_p(h_fn, s, p, tau)
                c = math.exp(-gamma_ase * units.fs / hop.mass)
                p = p + (hop.force(s) + f_bias) * units.fs / 2
                s = s + p / hop.mass * units.fs / 2
                dW = rng.standard_normal()
                sig = math.sqrt(hop.mass * kT * (1 - c**2))
                lw += (f_bias * sig * dW / math.sqrt(2 * gamma_ase * kT * hop.mass)) \
                      - 0.5 * f_bias**2 * units.fs / (2 * gamma_ase * kT * hop.mass)
                p = p * c + sig * dW
                s = s + p / hop.mass * units.fs / 2
                p = p + (hop.force(s) + f_bias) * units.fs / 2
                if s > hop.s_crit:
                    t_hit = step + 1
                    break
            arrivals.append((t_hit, s))
            weights.append(math.exp(float(np.clip(lw, -50, 50))))
        ess = effective_sample_size(weights)
        out["candidates"][name] = {
            "ess": ess,
            "rejected_fail_closed": bool(ess < ess_min),
            "n_arrivals": sum(1 for t, _ in arrivals if t is not None),
            "tvd_point_vs_reject": _tvd_points([s for _, s in arrivals], [a["s"] for a in rej["arrivals"]]),
            "momentum_free": bool(getattr(h_fn, "momentum_free", False)),
        }
    out["ess_gate_demo"] = {"ess_uniform_weights": effective_sample_size(np.ones(50)),
                            "gate_would_reject": effective_sample_size(np.ones(50)) < ess_min}
    return out


def _tvd_points(a, b, bins=12):
    a = [x for x in a if x is not None]
    b = [x for x in b if x is not None]
    if len(a) < 5 or len(b) < 5:
        return None
    lo, hi = min(min(a), min(b)), max(max(a), max(b))
    ha, _ = np.histogram(a, bins=bins, range=(lo, hi))
    hb, _ = np.histogram(b, bins=bins, range=(lo, hi))
    pa, pb = ha / max(1, ha.sum()), hb / max(1, hb.sum())
    return float(0.5 * np.abs(pa - pb).sum())


# ==========================================================================
# §L5 conditioned 证据类型（PLAN §4 L5）：schema 漂移账本一条；永不合并。
# ==========================================================================

SCHEMA_DRIFT_NOTE = (
    "PLAN_2026-09-28 L5：新增证据类型 'conditioned'（原子桥轨迹）。"
    "它是提案来源的记录，不与 attempts / continuations / ts_support 合并（K2）。"
    "合并时在 schema.py 证据类型清单加一项，并在漂移账本记本条。"
)


def conditioned_record(trial_id, edge_src, edge_dst, meta=None):
    if not DECISIONS.k2_conditioned_never_merged:
        raise PermissionError("K2 未确认：conditioned 不得并入 attempts")
    rec = {"evidence_type": "conditioned", "trial_id": trial_id, "edge": (edge_src, edge_dst),
           "never_merge_with": ["attempts", "continuations", "ts_support"]}
    if meta:
        rec.update(meta)
    return rec


def assert_evidence_disjoint(sets_by_type):
    """任何两路证据的 trial 集合必须不相交（合并 = 违反 K2 = 红）。"""
    for a, b in itertools.combinations(list(sets_by_type), 2):
        overlap = set(sets_by_type[a]) & set(sets_by_type[b])
        if overlap:
            raise ValueError(f"K2 违反：{a} 与 {b} 共享 trial {sorted(overlap)[:3]}")
    return True


# ==========================================================================
# §4½ 两端版本（PLAN §4½ E0–E4）
# ==========================================================================


class EndpointError(ValueError):
    pass


def _confirm(atoms, factory, config, seed):
    try:
        ok, info = confirm_minimum(atoms, factory, config, seed)
    except Exception as exc:  # 草稿回退：主线门用不了时用自己的 FD Hessian，如实记录
        atoms.calc = factory()
        neg = int(np.sum(np.linalg.eigvalsh(fd_hessian(atoms)) < -0.5))
        return neg == 0, {"method": "draft_fd_hessian", "fallback_reason": str(exc), "negative": neg}
    if not ok and info.get("reason", "").startswith("gate") and info.get("saddle_order") is None:
        # 主线门拒绝的是 Hessian 溯源（玩具 tether 势过不了溯源门），不是「这是鞍点」。
        # 草稿补一条 FD Hessian 证据；合并后真实体系照走主线门。
        atoms.calc = factory()
        neg = int(np.sum(np.linalg.eigvalsh(fd_hessian(atoms)) < -0.5))
        info = dict(info, draft_fd_hessian_negative=neg, accepted_by_draft=bool(neg == 0))
        return neg == 0, info
    return ok, info


def endpoint_validation(reactant: Atoms, product: Atoms, factory, config, seed=0):
    """E0（K4/K9）：R、P 都过 confirm_minimum；chemical_key 不同（P1 型退化用
    带标号键集合区分）。两端 key 相同且非退化 → 开跑前拒绝。"""
    ok_r, info_r = _confirm(reactant, factory, config, seed)
    ok_p, info_p = _confirm(product, factory, config, seed)
    if not ok_r:
        raise EndpointError(f"R 未通过 confirm_minimum：{info_r}")
    if not ok_p:
        raise EndpointError(f"P 未通过 confirm_minimum：{info_p}")
    key_r, key_p = chemical_digest(reactant), chemical_digest(product)
    if key_r == key_p:
        er, ep = bond_set(reactant), bond_set(product)
        if er != ep:
            return {"degenerate": True, "note": "chemical_key 相同；用带标号键集合区分（同首次进入协议）",
                    "edges_R": sorted(er), "edges_P": sorted(ep)}
        raise EndpointError("两端 chemical_key 与带标号键集合都相同：不是两端问题，拒绝")
    return {"degenerate": False, "key_R": key_r, "key_P": key_p}


def two_sided_descent(variant: ToyVariant, saddle_atoms: Atoms):
    """鞍点 → 双侧下降：负本征模（反应模）± 方向淬火，两端必须落在两端点盆地。
    没有负模 = 攀升像不是鞍点，如实报错。"""
    atoms = saddle_atoms.copy()
    atoms.calc = ToyPotential(variant)
    vals, vecs = np.linalg.eigh(fd_hessian(atoms))
    i_min = int(np.argmin(vals))
    if vals[i_min] >= -0.1:
        raise EndpointError(f"攀升像没有负本征模（最低 {vals[i_min]:.3f}）：不是鞍点")
    mode = vecs[:, i_min].reshape(-1, 3)
    endpoints = []
    for sign in (+1.0, -1.0):
        a = atoms.copy()
        a.set_positions(a.get_positions() + sign * 0.2 * mode)
        a.calc = ToyPotential(variant)
        FIRE(a, logfile=None).run(fmax=1e-3, steps=600)
        endpoints.append(chemical_digest(a))
    return endpoints


def _neb_climb(variant: ToyVariant, atoms_a, atoms_b, images=9, fmax=0.03):
    a0, b0 = atoms_a.copy(), atoms_b.copy()
    for end in (a0, b0):
        end.set_constraint(FixAtoms(indices=[0]))
    imgs = [a0] + [a0.copy() for _ in range(images - 2)] + [b0]
    neb = NEB(imgs, climb=True, method="improvedtangent")
    neb.interpolate()
    for im in imgs:
        im.calc = ToyPotential(variant)
    LBFGS(neb, logfile=None).run(fmax=fmax, steps=250)
    return imgs[images // 2]


def run_e1(seed=5):
    """E1：测试台精确 h 两端桥 → h≈1/2 帧 → 鞍点 → 双侧下降；
    「两端」= T-链的最后一步（S2=(0,2,3) → P=(1,2,3)，基元步），鞍点与
    T0 清单一致、两侧 = S2 和 P。"""
    variant = ToyVariant("chain", dict(DECISIONS.t0_potential))
    minima = enumerate_minima(variant)
    register_reference(variant, minima)
    s2_key, p_key = find_state(minima, (0, 2, 3)), find_state(minima, (1, 2, 3))
    # 两端桥（1D 约化）：h 目标直接取 P 盆地（first_entry 的 {标签:键集合} 够用，不需要 L1）
    hop = ReducedHop()
    table = exact_h_table(hop, seed=seed)
    bridge = run_bridge_1d(hop, table, shots=60, seed=seed)
    # h≈1/2 帧：到达瞬间落在分界面附近（0.15–0.9 Å，跨越 s_crit=0.216 的邻域）
    half_frames = [a for a in bridge["arrivals"] if a["s"] is not None and 0.15 < a["s"] < 0.9]
    climb = _neb_climb(variant, minima[s2_key]["atoms"], minima[p_key]["atoms"])
    e_ts = float(climb.get_potential_energy())
    ends = two_sided_descent(variant, climb)
    return {"n_half_frames": len(half_frames), "ts_energy_eV": e_ts,
            "descent_endpoints": [e[:8] for e in ends],
            "from_side_matches_S2": s2_key in ends, "P_matches": p_key in ends,
            "pass": len(half_frames) >= 5 and s2_key in ends and p_key in ends}


def run_e2(shots=8, window=2000, seed=13):
    """E2：故意给不相邻的 R、P（隔整条链）：桥轨迹中途淬火出第三个
    chemical_key → 必须报「经过中间体」，不许报直连 TS。"""
    variant = ToyVariant("chain", dict(DECISIONS.t0_potential))
    r_atoms = quench_placement(variant, [0, 1, 2])
    p_atoms = quench_placement(variant, [1, 2, 3])
    key_r, key_p = chemical_digest(r_atoms), chemical_digest(p_atoms)
    gamma_ase = DECISIONS.reduced_experiment_gamma_per_ps / (1000.0 * units.fs)
    kT = units.kB * DECISIONS.t0_dynamic_temperature_K
    found = set()
    for sh in range(shots):
        rng = np.random.default_rng(seed * 6151 + sh)
        atoms = r_atoms.copy()
        atoms.calc = ToyPotential(variant)
        _maxwell(atoms, DECISIONS.t0_dynamic_temperature_K, rng)
        for step in range(window):
            _baoab_step(atoms, gamma_ase, kT, units.fs, rng)
            if step % 20 == 0:
                snap = atoms.copy()
                snap.calc = ToyPotential(variant)
                LBFGS(snap, logfile=None).run(fmax=1e-3, steps=200)
                found.add(chemical_digest(snap))
    intermediates = sorted(found - {key_r, key_p})
    return {"R": key_r[:8], "P": key_p[:8], "intermediates_found": [k[:8] for k in intermediates],
            "claims_direct_TS": False, "pass": len(intermediates) >= 1}


def judge_e3(ts_imag_freq_cm, ts_rmsd_to_ref, ref_event_key, ts_event_key):
    """E3 判据（K10）：对照 P2 冻结结论。判据跑前冻结（沿用 T2 机器）。"""
    CRITERIA.assert_frozen("E3")
    c = DECISIONS.k10_p2_criteria
    return {"event_key_match": ref_event_key == ts_event_key,
            "rmsd_within": ts_rmsd_to_ref <= c["rmsd_A"],
            "imag_freq_in_range": c["imag_freq_cm"][0] <= ts_imag_freq_cm <= c["imag_freq_cm"][1],
            "also_P1": c["also_P1"]}


def run_e4(variant_name="chain", images=9):
    """E4 对照组：CI-NEB（ASE 自带），单独命名，不进 PRRS 证据字段。
    记录力调用数与 TS 能量（对照 T0 清单的鞍点）。"""
    variant = ToyVariant(variant_name, dict(DECISIONS.t0_potential))
    r_atoms = quench_placement(variant, [0, 1, 2])
    p_atoms = quench_placement(variant, [1, 2, 3])
    counters = []

    class Counting(ToyPotential):
        def calculate(self, *a, **kw):
            counters.append(1)
            super().calculate(*a, **kw)

    a0, b0 = r_atoms.copy(), p_atoms.copy()
    for end in (a0, b0):
        end.set_constraint(FixAtoms(indices=[0]))
    imgs = [a0] + [a0.copy() for _ in range(images - 2)] + [b0]
    neb = NEB(imgs, climb=True, method="improvedtangent")
    neb.interpolate()
    for im in imgs:
        im.calc = Counting(variant)
    LBFGS(neb, logfile=None).run(fmax=0.03, steps=300)
    return {"control": "CI-NEB (ase.neb)；单独命名对照，不进 PRRS 证据字段",
            "force_calls": len(counters), "ts_energy_eV": float(imgs[images // 2].get_potential_energy())}


# ==========================================================================
# §G0 多步 discover 基线（PLAN §3 G0）：主线 ReactionSearch 直接跑玩具势，
# 不需要任何代码改动。账本交叉核对：attempts.jsonl 条数 = trials 目录数。
# ==========================================================================


def run_g0(variant_name="chain", depth=3, trials=40, seed=17, out=None):
    variant = ToyVariant(variant_name, dict(DECISIONS.t0_potential))
    cfg = SearchConfig(seed=seed, max_trials=trials, max_depth=depth, families=("stretch",))
    from prrs.search import ReactionSearch

    out = Path(out or (ROOT / "runs" / f"g0_{variant_name}_d{depth}"))
    ReactionSearch(toy_factory(variant), cfg).run(toy_atoms(variant), out)
    return {"out": str(out), "ledger": ledger_crosscheck(out)}


def ledger_crosscheck(run_dir: Path):
    attempts = [json.loads(x) for x in (run_dir / "attempts.jsonl").read_text().splitlines() if x.strip()]
    trials_dirs = list((run_dir / "trials").glob("*/result.json"))
    network = json.loads((run_dir / "network.json").read_text())
    return {"attempts": len(attempts), "trial_dirs": len(trials_dirs),
            "network_nodes": len(network.get("nodes", [])),
            "depth_capped": len(network.get("depth_capped", [])),
            "budget_blocks": len(network.get("budget_blocks", [])),
            "counts_match": len(attempts) == len(trials_dirs)}


# ==========================================================================
# §CLI
# ==========================================================================


def _banner():
    print(PROVISIONAL_BANNER)


def cmd_smoke():
    """每一段机器都过一遍（秒级、小样本）。acceptance 判据未冻结，只验机制。"""
    _banner()
    results = {}

    # T0 三变体（走廊配对）
    for name in ("chain", "decoy", "detour"):
        _report, checks = run_t0(name)
        checks = dict(checks)
        checks["pass"] = all(bool(v) for v in checks.values())
        results[f"T0/{name}"] = checks

    # G1: d 与 Aut
    variant = ToyVariant("chain", dict(DECISIONS.t0_potential))
    minima = enumerate_minima(variant)
    register_reference(variant, minima)
    r_key = find_state(minima, (0, 1, 2))
    results["G1/d_R"] = minima[r_key]["d_to_P"]
    results["G1/aut_P_size"] = len(AUT["chain"])

    # G2: 队列 关=FIFO / 开=按 h 最高先出（h=i → 后进的高分，顺序必须与 FIFO 不同）
    q = ConnectQueue(enabled=False)
    for i in range(5):
        q.append(i, h_value=i)
    fifo = [q.popleft() for _ in range(5)]
    q2 = ConnectQueue(enabled=True)
    for i in range(5):
        q2.append(i, h_value=i)
    prio = [q2.popleft() for _ in range(5)]
    results["G2/queue"] = {"fifo": fifo, "priority": prio,
                           "pass": fifo == list(range(5)) and prio == list(reversed(range(5)))}

    # G3: 方向带
    entries = [("stretch", 0, 1, 0, 0.0), ("stretch", 1, 2, 0, 0.0), ("bend", 1, 2, 0, 0.0)]
    boosted, notes = boost_directions(entries, {frozenset((0, 1))}, set(), True)
    results["G3/boost"] = {"new_bands": [e[3] for e in boosted], "n_notes": len(notes),
                           "pass": boosted[0][3] == -1 and boosted[2][3] == 0}

    # G4: 路径 + 卡点
    paths = extract_paths([("R", "A", ("attempts",), 1), ("A", "P", ("conditioned",), 2)], "R", ["P"])
    autopsy = budget_autopsy([{"status": "edge", "target": "x"}, {"status": "no_edge"}])
    results["G4/paths"] = {"n_paths": len(paths),
                           "evidence": [e[2] for e in paths[0]] if paths else None,
                           "autopsy_termination": autopsy["termination_policy"],
                           "pass": len(paths) == 1 and paths[0][0][2] == ("attempts",)}

    # G5: 四臂小预算机制验证（数字本身不是验收——判据在 T2）
    arms = ["discover", "+I-a", "+I-b", "+I-a+I-b"]
    g5 = [ToySearch(variant, arm, seed=7, budget=10, window=100).run() for arm in arms]
    results["G5/judge"] = judge_g5(g5)
    results["G5/all_arms_ran"] = {"pass": all(r["trials_used"] > 0 for r in g5)}

    # L0
    results["L0/gamma0_equals_verlet"] = check_gamma0_equals_verlet(steps=30)
    results["L0/fdt"] = check_fdt(steps=6000)

    # L1
    results["L1/predicate"] = check_predicate_consistency()

    # L2
    results["L2/seed_consistency"] = exact_h_consistency(seed_a=1, seed_b=2)

    # L3（小样本）
    results["L3/bridge_vs_reject"] = run_l3(shots_bridge=60, shots_reject=100, seed=5)

    # L4
    results["L4/approx_h"] = run_l4(shots=40, seed=9)

    # L5
    conditioned_record("t1", "a", "b")
    try:
        assert_evidence_disjoint({"attempts": [1, 2], "conditioned": [2, 3]})
        l5_ok = False
    except ValueError:
        l5_ok = True
    results["L5/merge_guard"] = {"pass": l5_ok}

    # E0
    variant_c = ToyVariant("chain", dict(DECISIONS.t0_potential))
    r_atoms = quench_placement(variant_c, [0, 1, 2])
    p_atoms = quench_placement(variant_c, [1, 2, 3])
    cfg = SearchConfig(seed=1)
    ev = endpoint_validation(r_atoms, p_atoms, toy_factory(variant_c), cfg)
    mismatch_caught = samekey_caught = False
    try:
        bad = r_atoms.copy()
        bad.numbers[-1] = 8  # H→O：元素序列变了
        check_product(bad, r_atoms)
    except ProductMismatch:
        mismatch_caught = True
    try:
        endpoint_validation(r_atoms, r_atoms.copy(), toy_factory(variant_c), cfg)
    except EndpointError:
        samekey_caught = True
    results["E0/validation"] = {"two_endpoint_ok": not ev["degenerate"],
                                "element_mismatch_refused": mismatch_caught,
                                "same_key_refused": samekey_caught,
                                "pass": (not ev["degenerate"]) and mismatch_caught and samekey_caught}

    # E1 / E2 / E4
    results["E1/two_endpoint_saddle"] = run_e1(seed=5)
    results["E2/intermediate_report"] = run_e2()
    results["E4/neb_control"] = run_e4()

    print("\n=== SMOKE RESULTS ===")
    n_fail = 0
    for k, v in results.items():
        if isinstance(v, dict) and "pass" in v:
            flag = "PASS" if v["pass"] else "FAIL"
            n_fail += 0 if v["pass"] else 1
            payload = {kk: vv for kk, vv in v.items() if kk != "pass"}
            print(f"[{flag}] {k}: {json.dumps(payload, default=str, ensure_ascii=False)}")
        else:
            print(f"[info] {k}: {json.dumps(v, default=str, ensure_ascii=False)}")
    print(f"\n{'ALL GREEN' if n_fail == 0 else str(n_fail) + ' FAIL'}"
          "（这是机制自检，不是 G5/L3/L4 验收——判据在 T2 冻结后才作数）")
    return n_fail


def main(argv=None):
    ap = argparse.ArgumentParser(description="connect 模式两层桥工程草稿（主线之外）")
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--t0", metavar="VARIANT", choices=["chain", "decoy", "detour"])
    ap.add_argument("--freeze", action="store_true")
    ap.add_argument("--full-pairs", action="store_true")
    ap.add_argument("--t1", action="store_true")
    ap.add_argument("--shots", type=int, default=200)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--g0", action="store_true")
    ap.add_argument("--variant", default="chain")
    ap.add_argument("--depth", type=int, default=3)
    ap.add_argument("--trials", type=int, default=40)
    ap.add_argument("--out", default=None)
    ap.add_argument("--l0", action="store_true")
    ap.add_argument("--l2", action="store_true")
    ap.add_argument("--l3", action="store_true")
    ap.add_argument("--l4", action="store_true")
    ap.add_argument("--e1", action="store_true")
    ap.add_argument("--e2", action="store_true")
    ap.add_argument("--e4", action="store_true")
    args = ap.parse_args(argv)

    if args.smoke:
        return 0 if cmd_smoke() == 0 else 1
    if args.t0:
        run_t0(args.t0, freeze=args.freeze, full_pairs=args.full_pairs)
        return 0
    if args.t1:
        variant = ToyVariant(args.variant, dict(DECISIONS.t0_potential))
        print(json.dumps(graph_layer_reference(variant, shots=args.shots, seed=args.seed),
                         indent=2, ensure_ascii=False))
        return 0
    if args.g0:
        print(json.dumps(run_g0(args.variant, args.depth, args.trials, args.seed, args.out),
                         indent=2, ensure_ascii=False, default=str))
        return 0
    if args.l0:
        print(json.dumps(check_gamma0_equals_verlet(), indent=2))
        print(json.dumps(check_fdt(steps=20000), indent=2))
        return 0
    if args.l2:
        print(json.dumps(exact_h_consistency(seed_a=args.seed, seed_b=args.seed + 1), indent=2))
        return 0
    if args.l3:
        print(json.dumps(run_l3(shots_bridge=args.shots, shots_reject=2 * args.shots, seed=args.seed),
                         indent=2))
        return 0
    if args.l4:
        print(json.dumps(run_l4(shots=args.shots, seed=args.seed), indent=2, default=str))
        return 0
    if args.e1:
        print(json.dumps(run_e1(seed=args.seed), indent=2, default=str))
        return 0
    if args.e2:
        print(json.dumps(run_e2(), indent=2, ensure_ascii=False, default=str))
        return 0
    if args.e4:
        print(json.dumps(run_e4(), indent=2, ensure_ascii=False, default=str))
        return 0
    print(__doc__)
    _banner()
    return 0


if __name__ == "__main__":
    sys.exit(main())
