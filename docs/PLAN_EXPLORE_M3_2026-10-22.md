# 发现层 v1 · M3 姿态与意图 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 建 L0（来源构象 + 双分子姿态）与 v1 意图生成（距离原语 + 二元组合），并用纯几何数出各开发体系的“意图数 × 姿态数”，交用户据此定 L1 规模。

**Architecture:** 在独立包 `explore/` 内新增 `poses.py`、`intents.py`，单向 import `prrs`（片段、颜色、扭转、自同构、淬火与极小确认原语全部复用）。`probe.py` 只改一处：多对意图的可行域（I4）。计数脚本放 `benchmarks/explore_counts/`，只做几何，不加载模型。

**Tech Stack:** Python 3.12（`openmm_dev`，生产环境，不装不升）、ASE 3.29.0、numpy、scipy 1.18.1（`Rotation`）、RDKit 2025.09.6（仅计数脚本里给 N1 生成初始几何）、pytest；ruff 在 `.venv-tools`。

**Spec:** [发现层设计](DISCOVERY_LAYER_DESIGN_2026-10-13.md) §1.1、§2、§3、§9（按 v1 范围）；[评估协议](EVALUATION_PROTOCOL_DESIGN_2026-10-13.md) §2.1、§2.3、§5.3；上位计划 [PLAN_EXPLORE_V1](PLAN_EXPLORE_V1_2026-10-13.md) 里程碑表 M3 行。

状态：**2026-10-22 已施工完成**（用户未逐条审计划，指示“先推主线”后按 Subagent-driven 执行；施工中的裁决 R1–R11 见 DONE_LOG 末节）。U1、U3、U4、U5 仍待用户依据计数表决定。

## 必读（开工前全部读完）

| 文件 | 为什么要读 |
| --- | --- |
| `CLAUDE.md`（项目根目录） | 设计稿表；`openmm_dev` 是生产环境 |
| `docs/DISCOVERY_LAYER_DESIGN_2026-10-13.md` | §1.1 v1 范围；§2 L0 的硬约束（原始姿态与极小分开存、不从已知 TS 建起点、ID 与随机流稳定）；§3 意图规则 |
| `docs/EVALUATION_PROTOCOL_DESIGN_2026-10-13.md` §2、§5.3 | 共享准备成本 P 怎么记；搜索端不可见参考答案 |
| `docs/PLAN_EXPLORE_V1_2026-10-13.md` Global Constraints、Task 3 | 本计划沿用的约束；`Intent`/`constrained_scan` 的现状 |
| `explore/src/explore/probe.py`、`ledger.py` | `Intent(intent_id, pairs, coefs)`、`_feasible`、`Ledger.segment` |
| `src/prrs/perturbations.py`：`rotatable_torsions`、`step_torsion`、`_hop_counts`、`_bond_directions` 注 [1] | 扭转枚举与“二三跳内的非键对是伪装的弯曲”规则 |
| `src/prrs/chemistry.py`：`canonical_labels`、`_fragments`、`_formula`、`automorphisms`、`symmetric_rmsd` | 片段与对称性原语 |
| `src/prrs/runner.py`：`quench`、`confirm_minimum` docstring | 构象资格的失败语义 |

## 本计划需要用户先定的事

| # | 事项 | 计划里的建议 | 依据 / 为什么要你定 |
| --- | --- | --- | --- |
| U1 | 构象网格 | 每个 genuine 扭转在其基本域 `[0, 2π/order)` 上取 `steps_per_turn // order` 个点，`steps_per_turn = 6`（60°）；转子（甲基等对称 rotor）不枚举；环构象 v1 不覆盖，记为排除项 | 60° 是能同时含 0° 与 180°（N1 丙烯醛 s-cis/s-trans）的最粗均匀格；这是新规则，按全局规则 §8 由你定 |
| U2 | （已删）严格重复容差 | 代码不做网格点去重；资格后的去重沿用 PRRS `basin_rmsd_A` / `basin_energy_eV` | 不新设判据 |
| U3 | 姿态壳层 | 壳层用“片段间最近原子对距离 /（两原子共价半径之和）”的比值 f 定义；候选值不预设，由 Task 5 先实测 b04/b05/b09 现有反应物复合物里的 f，再由你挑 | 不同元素对的接触距离不同，比值比 Å 可迁移；具体数值没有依据前不写 |
| U4 | 二元组合的局部规则 | Task 4 实现三种（`shared_atom` / `adjacent` / `all`），Task 5 把三种的数量都报给你，你选一种冻结 | 设计 §3.1 要求局部限定并保存排除数；哪种合理要看数量 |
| U5 | 远层（> `pair_cutoff_A`）非键对 | 按颜色类对整类抽样，每姿态 `n_far` 类；`n_far` 由你看完计数再定 | 设计 §3.1“远壳层保留有界样本” |

数值参数 `pair_cutoff_A = 4.0`、`min_pair_hops = 4`、`min_distance_ratio = 0.45`、`basin_rmsd_A`、`basin_energy_eV` 都是 PRRS `SearchConfig` 现有参数，直接沿用，不新设。

## Global Constraints

- 科学环境只用 `/home/ruigengji/miniforge3/envs/openmm_dev/bin/python`；不 pip/conda 安装、卸载、升级。
- `prrs`、`pars` 不得 import `explore`；`explore` 单向 import `prrs`。
- 不从已知 TS、参考产物或手工对准的目标位点构造盲搜起点（设计 §2.1）；计数脚本可读 `probe_b09/input.extxyz` 只用于计数，并在表里注明“该几何对 TS-2OMe 有偏，只计数不作盲搜起点”。
- 原始姿态与松弛后极小分开保存；本里程碑不对姿态做模型松弛（见“明确不做”）。
- 稳定 ID：由来源/片段/网格索引/样本序号/族/坐标/符号/规则版本导出；追加样本不改已有样本的几何和 ID；随机流按 (seed, 阶段, 键) 划分（设计 §3.2）。
- 不用 hash 判断缓存或续跑能否复用；`zlib.crc32` 只用于把字符串键变成随机流种子。
- 不新增会拦住运行或结果的门；I4 可行域是设计 §3.2 明文要求的“约束不可行在求值前拒绝”。
- 输出不写 `/tmp`；结果放 `benchmarks/explore_counts/results/<YYYY-MM-DD>/`。
- NFS 上不跑 git；每个 Task 结束写 `.commit_drafts/<日期_时间>_<主题>.md`。
- 测试：只跑 `make -C explore check`（本包），不跑根目录全量。

## Review Focus

1. **构象网格点与自身碰撞**（扭转转到原子重叠或键图改变）：记 `status="overlap"` / `"graph_changed"` 并保留，不静默丢弃。Task 1 测试固定。
2. **两个相同片段**（同分子二聚）：片段先后由几何决定，原子重标号后仍得到同一组姿态。Task 3 测试固定。
3. **单原子 / 双原子片段**（Cl⁻、OH⁻）：没有扭转，只有输入构象；Kabsch 对齐在线性片段上不出 NaN。Task 1 测试固定。
4. **三个及以上片段**：v1 明确拒绝（`ValueError("v1 poses support 1 or 2 fragments")`），不猜怎么摆。Task 3 测试固定。
5. **多对意图目标落在重叠下限以下**：在任何 forward 之前判 `infeasible`。Task 4 测试固定。

## 明确不做（本里程碑）

- 姿态的模型松弛（无约束极小 / 保持间距的约束准备）：要用模型，放到 M5 编排或 M10 零号批里随成本账本一起做。
- 角/扭转/片段运动意图族、高阶组合、近似去重（v2）。
- 环构象枚举。
- 任何模型求值的实跑；Task 2 的资格函数只在解析势上测。

---

## 文件结构

```
explore/src/explore/
  streams.py      随机流：stream(seed, *parts) -> np.random.Generator
  poses.py        片段、构象网格、构象资格、双分子摆放、L0 记录读写
  intents.py      距离原语、二元组合、远层抽样、计数
  probe.py        （改）_feasible 换成按重叠下限的可行域
explore/tests/
  molecules.py    测试用小分子（乙醇、Cl⁻+CH3Cl、水二聚体）与重标号工具
  test_poses.py
  test_intents.py
  test_probe.py   （加 I4 测试）
benchmarks/explore_counts/
  count_l0.py     纯几何计数（b04、b05、N1、b09）
  results/<date>/ counts.json、contact_ratios.json、COUNTS.md
```

---

### Task 1: 随机流、片段与构象网格

**Files:**
- Create: `explore/src/explore/streams.py`
- Create: `explore/src/explore/poses.py`
- Create: `explore/tests/molecules.py`
- Test: `explore/tests/test_poses.py`

**Interfaces:**
- Consumes: `prrs.state.encode`、`prrs.chemistry.canonical_labels/_fragments/_formula/automorphisms/symmetric_rmsd`、`prrs.perturbations.rotatable_torsions/step_torsion`、`prrs.config.SearchConfig`
- Produces:
  - `stream(seed: int, *parts) -> np.random.Generator`
  - `Fragment(key: str, atoms: tuple[int, ...])`
  - `fragments(atoms, config) -> list[Fragment]`
  - `Conformer(conformer_id: str, fragment_key: str, atoms: tuple[int, ...], positions: np.ndarray, grid: tuple[int, ...], status: str)`
  - `conformer_grid(atoms, fragment, config, steps_per_turn: int, cap: int, seed: int) -> tuple[list[Conformer], dict]`
  - `RULE_VERSION = "L0-v1"`

- [ ] **Step 1: 写测试分子工具**

```python
# explore/tests/molecules.py
import numpy as np
from ase import Atoms
from ase.build import molecule


def ethanol():
    return molecule("CH3CH2OH")


def chloride_methyl_chloride():
    # b04-like: CH3Cl plus a free Cl- on the back side, charge handled by the config
    return Atoms(
        "CClHHHCl",
        positions=[
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.785],
            [1.028, 0.0, -0.373],
            [-0.514, 0.89, -0.373],
            [-0.514, -0.89, -0.373],
            [0.0, 0.0, -3.6],
        ],
    )


def water_dimer():
    a = molecule("H2O")
    b = molecule("H2O")
    b.positions += [0.3, 0.2, 3.0]
    return a + b


def relabel(atoms, seed=0):
    """A random atom permutation; returns (permuted atoms, perm) with new[k] = old[perm[k]]."""
    perm = np.random.default_rng(seed).permutation(len(atoms))
    return atoms[perm], perm
```

- [ ] **Step 2: 写失败测试**

```python
# explore/tests/test_poses.py
import numpy as np
import pytest

from explore.poses import conformer_grid, fragments
from explore.streams import stream
from molecules import chloride_methyl_chloride, ethanol, relabel, water_dimer
from prrs.config import SearchConfig
from prrs.perturbations import rotatable_torsions
from prrs.state import encode
from prrs.chemistry import canonical_labels

CONFIG = SearchConfig(closed_shell_only=False)


def test_stream_is_keyed_not_ordered():
    a = stream(0, "poses", "x", 3).normal(size=4)
    b = stream(0, "poses", "x", 3).normal(size=4)
    c = stream(0, "poses", "x", 4).normal(size=4)
    assert np.array_equal(a, b) and not np.array_equal(a, c)


def test_fragments_are_ordered_by_identity_not_index():
    atoms = chloride_methyl_chloride()
    frags = fragments(atoms, CONFIG)
    moved, perm = relabel(atoms, seed=3)
    frags_moved = fragments(moved, CONFIG)
    assert [f.key for f in frags] == [f.key for f in frags_moved]
    for f, g in zip(frags, frags_moved):
        assert sorted(int(perm[i]) for i in g.atoms) == list(f.atoms)


def test_identical_fragments_are_ordered_by_geometry():
    atoms = water_dimer()
    frags = fragments(atoms, CONFIG)
    moved, perm = relabel(atoms, seed=5)
    frags_moved = fragments(moved, CONFIG)
    for f, g in zip(frags, frags_moved):
        assert sorted(int(perm[i]) for i in g.atoms) == list(f.atoms)


def test_ethanol_grid_has_one_genuine_torsion_and_six_points():
    atoms = ethanol()
    (frag,) = fragments(atoms, CONFIG)
    confs, report = conformer_grid(atoms, frag, CONFIG, steps_per_turn=6, cap=1000, seed=0)
    assert report["genuine_torsions"] == 1 and report["rotors_not_enumerated"] == 1
    assert len(confs) == 6
    zero = next(c for c in confs if c.grid == (0,))
    assert np.allclose(zero.positions, atoms.positions[list(frag.atoms)], atol=1e-10)
    graph = encode(atoms, CONFIG.bond_scale)
    assert all(c.status == "generated" for c in confs)
    for c in confs:
        moved = atoms.copy()
        moved.positions[list(frag.atoms)] = c.positions
        assert encode(moved, CONFIG.bond_scale).edges == graph.edges


def test_grid_dihedrals_step_by_sixty_degrees():
    atoms = ethanol()
    (frag,) = fragments(atoms, CONFIG)
    confs, _ = conformer_grid(atoms, frag, CONFIG, steps_per_turn=6, cap=1000, seed=0)
    graph = encode(atoms, CONFIG.bond_scale)
    labels = canonical_labels(atoms.numbers, graph.edges, graph.index)
    genuine, _ = rotatable_torsions(graph, labels, {int(i): int(i) for i in graph.index})
    idx = list(genuine[0].indices)
    base = atoms.get_dihedral(*idx)
    seen = set()
    for c in confs:
        moved = atoms.copy()
        moved.positions[list(frag.atoms)] = c.positions
        seen.add(round(((moved.get_dihedral(*idx) - base) % 360.0) / 60.0) % 6)
    assert seen == set(range(6))


def test_monatomic_and_diatomic_fragments_have_only_the_input():
    atoms = chloride_methyl_chloride()
    frags = fragments(atoms, CONFIG)
    cl = next(f for f in frags if len(f.atoms) == 1)
    confs, report = conformer_grid(atoms, cl, CONFIG, steps_per_turn=6, cap=10, seed=0)
    assert len(confs) == 1 and confs[0].grid == () and report["genuine_torsions"] == 0
    assert np.isfinite(confs[0].positions).all()


def test_cap_keeps_input_and_reports_censoring():
    atoms = ethanol()
    (frag,) = fragments(atoms, CONFIG)
    confs, report = conformer_grid(atoms, frag, CONFIG, steps_per_turn=6, cap=3, seed=0)
    assert len(confs) == 3 and any(c.grid == (0,) for c in confs)
    assert report["grid_size"] == 6 and report["censored"] == 3


def test_clashing_grid_point_is_kept_with_status():
    atoms = ethanol()
    (frag,) = fragments(atoms, CONFIG)
    tight = SearchConfig(closed_shell_only=False, min_distance_ratio=0.99)
    confs, report = conformer_grid(atoms, frag, tight, steps_per_turn=6, cap=1000, seed=0)
    assert len(confs) == 6
    assert report["status_counts"].get("overlap", 0) >= 1
```

- [ ] **Step 3: 跑测试确认失败**

Run: `cd /home/ruigengji/MLPATH/explore && /home/ruigengji/miniforge3/envs/openmm_dev/bin/python -m pytest -q tests/test_poses.py`
Expected: FAIL，`ModuleNotFoundError: explore.poses`（`tests/` 需在 `pythonpath`：若 `pyproject.toml` 的 `pythonpath` 只有 `["src", "../src"]`，加上 `"tests"`）。

- [ ] **Step 4: 实现 `streams.py` 与 `poses.py` 的片段、构象部分**

```python
# explore/src/explore/streams.py
"""Keyed random streams (DISCOVERY_LAYER_DESIGN §3.2): a draw is addressed by its key,
so appending candidates never moves the streams of existing ones. crc32 only turns a
string key into seed words; nothing is cached or reused by it."""

import zlib

import numpy as np


def stream(seed, *parts):
    words = [int(seed)] + [zlib.crc32(str(p).encode()) for p in parts]
    return np.random.default_rng(np.random.SeedSequence(words))
```

```python
# explore/src/explore/poses.py
"""L0 for discovery v1 (DISCOVERY_LAYER_DESIGN_2026-10-13 §2): source conformers on a
torsion grid, and rigid placements of two fragments. Pure geometry except
`qualify_conformers`, which spends model forwards through the ledger (shared cost P)."""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass

import numpy as np

from .streams import stream

RULE_VERSION = "L0-v1"


@dataclass(frozen=True)
class Fragment:
    key: str
    atoms: tuple[int, ...]


@dataclass
class Conformer:
    conformer_id: str
    fragment_key: str
    atoms: tuple[int, ...]
    positions: np.ndarray
    grid: tuple[int, ...]
    status: str


def _kabsch_onto(moving, target):
    """Rigidly superpose `moving` onto `target` (same atom order); exact for n <= 2."""
    a = moving - moving.mean(axis=0)
    b = target - target.mean(axis=0)
    if len(a) < 3:
        return a + target.mean(axis=0) if len(a) == 1 else _align_line(a, b) + target.mean(axis=0)
    u, _, vt = np.linalg.svd(a.T @ b)
    d = np.sign(np.linalg.det(u @ vt)) or 1.0
    rot = u @ np.diag([1.0, 1.0, d]) @ vt
    return a @ rot + target.mean(axis=0)


def _align_line(a, b):
    # a diatomic: rotate its axis onto the target axis; the axial angle is undefined
    ua, ub = a[1] - a[0], b[1] - b[0]
    ua, ub = ua / np.linalg.norm(ua), ub / np.linalg.norm(ub)
    from scipy.spatial.transform import Rotation

    rot, _ = Rotation.align_vectors([ub], [ua])
    return rot.apply(a)


def fragments(atoms, config):
    from prrs.chemistry import _formula, _fragments, canonical_labels
    from prrs.state import encode

    graph = encode(atoms, config.bond_scale)
    labels = canonical_labels(atoms.numbers, graph.edges, graph.index)
    out = []
    for group in _fragments(graph.edges, graph.index):
        colours = ",".join(sorted(str(labels[i]) for i in group))
        out.append(Fragment(f"{_formula(atoms.numbers, group)}|{colours}", tuple(sorted(group))))
    # identical fragments tie on key: break by centroid, which relabelling cannot move
    return sorted(
        out,
        key=lambda f: (f.key, tuple(np.round(atoms.positions[list(f.atoms)].mean(axis=0), 8))),
    )


def _fragment_atoms(atoms, fragment):
    return atoms[list(fragment.atoms)]


def _ordered_torsions(sub, config):
    """Every genuine torsion of the fragment (identity orbits: two equivalent methoxy
    groups are two coordinates), in an order relabelling cannot change."""
    from prrs.chemistry import canonical_labels
    from prrs.perturbations import rotatable_torsions
    from prrs.state import encode

    graph = encode(sub, config.bond_scale)
    labels = canonical_labels(sub.numbers, graph.edges, graph.index)
    genuine, rotors = rotatable_torsions(graph, labels, {int(i): int(i) for i in graph.index})

    def key(t):
        j, k = t.indices[1], t.indices[2]
        mid = (sub.positions[j] + sub.positions[k]) / 2
        return (tuple(sorted((str(labels[j]), str(labels[k])))), tuple(np.round(mid, 8)))

    return sorted(genuine, key=key), rotors, graph


def _clashes(positions, numbers, config):
    from ase.data import covalent_radii

    r = covalent_radii[numbers]
    d = np.linalg.norm(positions[:, None] - positions[None], axis=-1)
    np.fill_diagonal(d, np.inf)
    return bool(np.any(d < config.min_distance_ratio * (r[:, None] + r[None])))


def conformer_grid(atoms, fragment, config, steps_per_turn, cap, seed):
    from prrs.perturbations import step_torsion
    from prrs.state import encode

    sub = _fragment_atoms(atoms, fragment)
    torsions, rotors, graph = _ordered_torsions(sub, config)
    axes = [range(max(1, steps_per_turn // t.symmetry_order)) for t in torsions]
    grid = list(itertools.product(*axes))
    origin = tuple(0 for _ in torsions)
    if len(grid) > cap:
        rest = [g for g in grid if g != origin]
        pick = stream(seed, "conformer_cap", fragment.key).permutation(len(rest))[: cap - 1]
        chosen = [origin] + [rest[i] for i in sorted(pick)]
    else:
        chosen = grid
    confs, counts = [], {}
    for g in chosen:
        work = sub.copy()
        for t, idx in zip(torsions, g):
            n = max(1, steps_per_turn // t.symmetry_order)
            step_torsion(work, t, idx * (2 * math.pi / t.symmetry_order) / n)
        status = "generated"
        if _clashes(work.positions, work.numbers, config):
            status = "overlap"
        elif encode(work, config.bond_scale).edges != graph.edges:
            status = "graph_changed"
        counts[status] = counts.get(status, 0) + 1
        confs.append(
            Conformer(
                conformer_id=f"{fragment.key}:t{'.'.join(map(str, g))}:{RULE_VERSION}",
                fragment_key=fragment.key,
                atoms=fragment.atoms,
                positions=_kabsch_onto(work.positions, sub.positions),
                grid=g,
                status=status,
            )
        )
    report = {
        "genuine_torsions": len(torsions),
        "rotors_not_enumerated": len(rotors),
        "rings_not_enumerated": True,
        "grid_size": len(grid),
        "censored": len(grid) - len(chosen),
        "status_counts": counts,
    }
    return confs, report
```

- [ ] **Step 5: 跑测试确认通过**

Run: `cd /home/ruigengji/MLPATH/explore && /home/ruigengji/miniforge3/envs/openmm_dev/bin/python -m pytest -q tests/test_poses.py`
Expected: PASS（8 passed）

- [ ] **Step 6: 写 commit 草稿**

新建 `.commit_drafts/<日期_时间>_explore-m3-conformers.md`：message `explore: keyed streams, fragments and torsion-grid conformers (M3)`；文件清单 `explore/src/explore/streams.py`、`poses.py`、`explore/tests/molecules.py`、`test_poses.py`（及 `pyproject.toml` 若改了 pythonpath）；附 `git add … && git commit -m "…"`。

---

### Task 2: 构象资格（淬火 + 极小确认 + 去重，记共享成本 P）

**Files:**
- Modify: `explore/src/explore/poses.py`（追加）
- Test: `explore/tests/test_poses.py`（追加）

**Interfaces:**
- Consumes: `Conformer`（Task 1）、`explore.ledger.Ledger.segment`、`prrs.runner.quench/confirm_minimum`、`prrs.reliability.GuardedCalculator/GateRejected`、`prrs.io.snapshot`、`prrs.chemistry.automorphisms/symmetric_rmsd`
- Produces: `qualify_conformers(atoms, conformers, factory, config, ledger, seed) -> list[dict]`，每项 `{"conformer_id", "status", "positions", "energy_eV", "duplicate_of"}`，`status ∈ {"qualified", "duplicate", "quench_failed", "not_minimum", "gate_rejected:<code>", "budget_exhausted", "skipped:<grid status>"}`

- [ ] **Step 1: 写失败测试**（双势阱 C2：极小在 1.2 Å 与 2.4 Å）

```python
from ase import Atoms

from explore.ledger import Ledger, counted_factory
from explore.poses import Conformer, qualify_conformers
from prrs.calculators import double_well_factory

DW = SearchConfig(closed_shell_only=False, quench_fmax_eV_A=0.001)


def _c2(r, cid):
    return Conformer(cid, "C2", (0, 1), np.array([[0.0, 0, 0], [r, 0, 0]]), (), "generated")


def test_qualification_dedupes_same_basin_and_keeps_distinct(tmp_path):
    ledger = Ledger(tmp_path / "l.jsonl", arm="P", budget=5000)
    template = Atoms("C2", positions=[[0, 0, 0], [1.2, 0, 0]])
    rows = qualify_conformers(
        template,
        [_c2(1.25, "a"), _c2(1.30, "b"), _c2(2.35, "c"), _c2(1.2, "d")],
        counted_factory(double_well_factory, ledger),
        DW,
        ledger,
        seed=0,
    )
    ledger.close()
    status = {r["conformer_id"]: r["status"] for r in rows}
    assert status["a"] == "qualified" and status["c"] == "qualified"
    assert status["b"] == "duplicate" and rows[1]["duplicate_of"] == "a"
    assert status["d"] == "duplicate"
    assert ledger.spent > 0
    stages = {e.get("stage") for e in map(json.loads, (tmp_path / "l.jsonl").read_text().splitlines())}
    assert stages <= {"prepare", None}


def test_grid_points_with_bad_status_are_skipped_without_forwards(tmp_path):
    ledger = Ledger(tmp_path / "l.jsonl", arm="P", budget=5000)
    bad = Conformer("x", "C2", (0, 1), np.array([[0.0, 0, 0], [0.3, 0, 0]]), (), "overlap")
    rows = qualify_conformers(
        Atoms("C2"), [bad], counted_factory(double_well_factory, ledger), DW, ledger, seed=0
    )
    assert rows[0]["status"] == "skipped:overlap" and ledger.spent == 0


def test_budget_exhaustion_is_a_status_not_an_exception(tmp_path):
    ledger = Ledger(tmp_path / "l.jsonl", arm="P", budget=3)
    rows = qualify_conformers(
        Atoms("C2"), [_c2(1.3, "a")], counted_factory(double_well_factory, ledger), DW, ledger, seed=0
    )
    assert rows[0]["status"] == "budget_exhausted"
```

（文件顶部补 `import json`。）

- [ ] **Step 2: 跑测试确认失败**

Run: `cd /home/ruigengji/MLPATH/explore && /home/ruigengji/miniforge3/envs/openmm_dev/bin/python -m pytest -q tests/test_poses.py -k qualif or skipped or budget`
Expected: FAIL，`ImportError: qualify_conformers`

- [ ] **Step 3: 实现**

```python
def qualify_conformers(atoms, conformers, factory, config, ledger, seed):
    """Quench each generated conformer, confirm the minimum, and merge strict duplicates
    with the existing PRRS basin rule (`basin_rmsd_A`, `basin_energy_eV`) minimised over
    the graph automorphisms. Every forward is shared preparation cost P (stage
    "prepare"); failures keep their status and are never dropped."""
    from prrs.chemistry import automorphisms, symmetric_rmsd
    from prrs.io import snapshot
    from prrs.reliability import GateRejected, GuardedCalculator
    from prrs.runner import confirm_minimum, quench

    from .ledger import BudgetExhausted

    kept, rows = [], []
    for conf in conformers:
        row = {"conformer_id": conf.conformer_id, "status": None, "positions": None,
               "energy_eV": None, "duplicate_of": None}
        rows.append(row)
        if conf.status != "generated":
            row["status"] = f"skipped:{conf.status}"
            continue
        sub = atoms[list(conf.atoms)].copy()
        sub.set_constraint()
        sub.positions = conf.positions.copy()
        try:
            guard = GuardedCalculator(factory(), config)
            sub.calc = guard
            with ledger.segment(task=conf.conformer_id, stage="prepare") as seg:
                converged, _ = quench(sub, guard, config)
                seg.outcome = "converged" if converged else "quench_failed"
                endpoint = snapshot(sub, guard) if converged else None
            if endpoint is None:
                row["status"] = "quench_failed"
                continue
            with ledger.segment(task=conf.conformer_id, stage="prepare") as seg:
                confirmed, _ = confirm_minimum(endpoint, factory, config, seed)
                seg.outcome = "minimum" if confirmed else "not_minimum"
            if not confirmed:
                row["status"] = "not_minimum"
                continue
        except BudgetExhausted:
            row["status"] = "budget_exhausted"
            continue
        except GateRejected as exc:
            row["status"] = f"gate_rejected:{exc.code}"
            continue
        energy = float(endpoint.get_potential_energy())
        row.update(positions=endpoint.positions.copy(), energy_eV=energy)
        group, _ = automorphisms(endpoint, config.bond_scale, None, config.automorphism_limit)
        for other_id, other in kept:
            if (
                abs(energy - other.get_potential_energy()) <= config.basin_energy_eV
                and symmetric_rmsd(endpoint, other, group) <= config.basin_rmsd_A
            ):
                row.update(status="duplicate", duplicate_of=other_id)
                break
        else:
            row["status"] = "qualified"
            kept.append((conf.conformer_id, endpoint))
    return rows
```

- [ ] **Step 4: 跑测试确认通过**

Run: `cd /home/ruigengji/MLPATH/explore && /home/ruigengji/miniforge3/envs/openmm_dev/bin/python -m pytest -q tests/test_poses.py`
Expected: PASS（11 passed）

- [ ] **Step 5: 写 commit 草稿** `.commit_drafts/<日期_时间>_explore-m3-qualify.md`（message `explore: conformer qualification charged as shared preparation (M3)`）。

---

### Task 3: 双分子姿态（摆放、采样、L0 记录）

**Files:**
- Modify: `explore/src/explore/poses.py`（追加）
- Test: `explore/tests/test_poses.py`（追加）

**Interfaces:**
- Consumes: `Fragment`、`Conformer`、`stream`、`fragments`
- Produces:
  - `Pose(pose_id: str, conformers: tuple[str, ...], sample: int | None, shell: float | None, direction: np.ndarray | None, quaternion: np.ndarray | None, positions: np.ndarray, min_contact_ratio: float | None, status: str)`
  - `contact_ratio(pos_a, num_a, pos_b, num_b) -> float`
  - `place(pos_a, num_a, pos_b, num_b, direction, quaternion, shell) -> tuple[np.ndarray, float, str]`（返回 B 的新坐标、实际比值、状态 `"placed"`）
  - `sample_poses(atoms, frags, conf_a: Conformer, conf_b: Conformer, shells: tuple[float, ...], n_samples: int, seed: int) -> list[Pose]`
  - `single_fragment_poses(atoms, conformers) -> list[Pose]`
  - `write_l0(path, poses, electronic_state: dict, generation: dict) -> None`（JSONL，字段见设计 §2.3）

- [ ] **Step 1: 写失败测试**

```python
from explore.poses import contact_ratio, place, sample_poses, single_fragment_poses, write_l0


def _pair():
    atoms = chloride_methyl_chloride()
    frags = fragments(atoms, CONFIG)
    confs = [conformer_grid(atoms, f, CONFIG, 6, 100, 0)[0][0] for f in frags]
    return atoms, frags, confs


def test_placement_hits_the_shell_ratio():
    atoms, frags, (ca, cb) = _pair()
    poses = sample_poses(atoms, frags, ca, cb, shells=(1.8, 2.4), n_samples=4, seed=0)
    assert len(poses) == 8
    for p in poses:
        a = p.positions[list(frags[0].atoms)]
        b = p.positions[list(frags[1].atoms)]
        r = contact_ratio(a, atoms.numbers[list(frags[0].atoms)], b, atoms.numbers[list(frags[1].atoms)])
        assert abs(r - p.shell) < 1e-6 and p.status == "placed"


def test_fragments_move_rigidly():
    atoms, frags, (ca, cb) = _pair()
    for p in sample_poses(atoms, frags, ca, cb, shells=(2.0,), n_samples=3, seed=1):
        for f, c in zip(frags, (ca, cb)):
            got = p.positions[list(f.atoms)]
            d0 = np.linalg.norm(c.positions[:, None] - c.positions[None], axis=-1)
            d1 = np.linalg.norm(got[:, None] - got[None], axis=-1)
            assert np.allclose(d0, d1, atol=1e-10)


def test_appending_samples_keeps_existing_poses():
    atoms, frags, (ca, cb) = _pair()
    few = sample_poses(atoms, frags, ca, cb, shells=(2.0,), n_samples=3, seed=2)
    more = sample_poses(atoms, frags, ca, cb, shells=(2.0,), n_samples=6, seed=2)
    for p, q in zip(few, more):
        assert p.pose_id == q.pose_id and np.array_equal(p.positions, q.positions)


def test_relabelled_input_gives_the_same_poses():
    atoms = water_dimer()
    frags = fragments(atoms, CONFIG)
    confs = [conformer_grid(atoms, f, CONFIG, 6, 100, 0)[0][0] for f in frags]
    poses = sample_poses(atoms, frags, *confs, shells=(2.0,), n_samples=3, seed=4)
    moved, perm = relabel(atoms, seed=7)
    mfr = fragments(moved, CONFIG)
    mconfs = [conformer_grid(moved, f, CONFIG, 6, 100, 0)[0][0] for f in mfr]
    mposes = sample_poses(moved, mfr, *mconfs, shells=(2.0,), n_samples=3, seed=4)
    for p, q in zip(poses, mposes):
        assert np.allclose(q.positions, p.positions[perm], atol=1e-8)


def test_three_fragments_are_refused():
    atoms = water_dimer() + water_dimer()[:3]
    atoms.positions[6:] += [0, 0, 8.0]
    frags = fragments(atoms, CONFIG)
    assert len(frags) == 3
    with pytest.raises(ValueError, match="1 or 2 fragments"):
        sample_poses(atoms, frags, None, None, shells=(2.0,), n_samples=1, seed=0)


def test_single_fragment_poses_are_the_conformers():
    atoms = ethanol()
    (frag,) = fragments(atoms, CONFIG)
    confs, _ = conformer_grid(atoms, frag, CONFIG, 6, 100, 0)
    poses = single_fragment_poses(atoms, confs)
    assert [p.conformers for p in poses] == [(c.conformer_id,) for c in confs]
    assert all(p.shell is None for p in poses)


def test_l0_record_has_the_contract_fields(tmp_path):
    atoms, frags, (ca, cb) = _pair()
    poses = sample_poses(atoms, frags, ca, cb, shells=(2.0,), n_samples=1, seed=0)
    write_l0(tmp_path / "l0.jsonl", poses, {"charge": -1, "multiplicity": 1},
             {"rule": "L0-v1", "seed": 0})
    rec = json.loads((tmp_path / "l0.jsonl").read_text().splitlines()[0])
    for k in ("source_microstate", "conformer", "pose", "atoms_map", "electronic_state",
              "generation_rule", "seed", "preparation", "qualification", "cost", "positions"):
        assert k in rec
    assert rec["preparation"] == "raw_placement"
```

- [ ] **Step 2: 跑测试确认失败**

Run: `cd /home/ruigengji/MLPATH/explore && /home/ruigengji/miniforge3/envs/openmm_dev/bin/python -m pytest -q tests/test_poses.py`
Expected: FAIL，`ImportError: contact_ratio`

- [ ] **Step 3: 实现**

```python
@dataclass
class Pose:
    pose_id: str
    conformers: tuple[str, ...]
    sample: int | None
    shell: float | None
    direction: np.ndarray | None
    quaternion: np.ndarray | None
    positions: np.ndarray
    min_contact_ratio: float | None
    status: str


def _radii(numbers):
    from ase.data import covalent_radii

    return covalent_radii[np.asarray(numbers)]


def contact_ratio(pos_a, num_a, pos_b, num_b):
    d = np.linalg.norm(pos_a[:, None] - pos_b[None], axis=-1)
    return float(np.min(d / (_radii(num_a)[:, None] + _radii(num_b)[None])))


def place(pos_a, num_a, pos_b, num_b, direction, quaternion, shell):
    """Put B's centroid on the ray from A's centroid along `direction`, B rotated about
    its centroid by `quaternion`, at the OUTERMOST distance where the closest-contact
    ratio equals `shell`. March inward with the Lipschitz-safe step
    (g - shell) * min(r_i + r_j), which can never jump over a crossing, then bisect."""
    from scipy.spatial.transform import Rotation

    centre_a = pos_a.mean(axis=0)
    body = Rotation.from_quat(quaternion).apply(pos_b - pos_b.mean(axis=0))
    rmin = float(np.min(_radii(num_a)[:, None] + _radii(num_b)[None]))
    reach = np.max(np.linalg.norm(pos_a - centre_a, axis=1)) + np.max(np.linalg.norm(body, axis=1))

    def g(t):
        return contact_ratio(pos_a, num_a, body + centre_a + t * direction, num_b)

    outer = reach + 2 * shell * rmin + 1.0
    t = outer
    while g(t) > shell:
        step = (g(t) - shell) * rmin
        if step < 1e-9:
            break
        t -= step
    lo, hi = t, outer
    while g(lo) > shell:  # the march stops at or just above the crossing
        lo -= 1e-9
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if g(mid) <= shell else (lo, mid)
        if hi - lo < 1e-10:
            break
    new_b = body + centre_a + hi * direction
    return new_b, g(hi), "placed"


def _unit(v):
    return v / np.linalg.norm(v)


def sample_poses(atoms, frags, conf_a, conf_b, shells, n_samples, seed):
    """Uniform approach directions (normalised 3-D Gaussian) and uniform orientations
    (normalised 4-D Gaussian = Haar on SO(3)), one keyed stream per sample index, so
    the k-th sample is the same however many are drawn. Each sample is placed at every
    shell."""
    if len(frags) != 2:
        raise ValueError("v1 poses support 1 or 2 fragments")
    fa, fb = frags
    num_a, num_b = atoms.numbers[list(fa.atoms)], atoms.numbers[list(fb.atoms)]
    pair_key = f"{conf_a.conformer_id}+{conf_b.conformer_id}"
    poses = []
    for k in range(n_samples):
        rng = stream(seed, "pose", pair_key, k)
        direction = _unit(rng.normal(size=3))
        quaternion = _unit(rng.normal(size=4))
        for shell in shells:
            new_b, ratio, status = place(
                conf_a.positions, num_a, conf_b.positions, num_b, direction, quaternion, shell
            )
            positions = atoms.positions.copy()
            positions[list(fa.atoms)] = conf_a.positions
            positions[list(fb.atoms)] = new_b
            poses.append(
                Pose(
                    pose_id=f"{pair_key}:k{k}:f{shell:g}:{RULE_VERSION}",
                    conformers=(conf_a.conformer_id, conf_b.conformer_id),
                    sample=k,
                    shell=float(shell),
                    direction=direction,
                    quaternion=quaternion,
                    positions=positions,
                    min_contact_ratio=ratio,
                    status=status,
                )
            )
    return poses


def single_fragment_poses(atoms, conformers):
    poses = []
    for c in conformers:
        positions = atoms.positions.copy()
        positions[list(c.atoms)] = c.positions
        poses.append(
            Pose(c.conformer_id, (c.conformer_id,), None, None, None, None, positions, None, c.status)
        )
    return poses


def write_l0(path, poses, electronic_state, generation):
    import json
    from pathlib import Path

    with Path(path).open("w", encoding="utf-8") as fh:
        for p in poses:
            fh.write(json.dumps({
                "source_microstate": None,  # set when the pose is qualified as a basin
                "conformer": list(p.conformers),
                "pose": p.pose_id,
                "atoms_map": "input order",
                "electronic_state": electronic_state,
                "generation_rule": {**generation, "shell": p.shell, "sample": p.sample},
                "seed": generation.get("seed"),
                "preparation": "raw_placement",
                "qualification": p.status,
                "cost": {"forwards": 0},
                "direction": None if p.direction is None else p.direction.tolist(),
                "quaternion": None if p.quaternion is None else p.quaternion.tolist(),
                "min_contact_ratio": p.min_contact_ratio,
                "positions": p.positions.tolist(),
            }) + "\n")
```

注意 `test_three_fragments_are_refused` 传 `None` 作构象：长度检查必须在读取构象之前，上面的实现已满足。

- [ ] **Step 4: 跑测试确认通过**

Run: `cd /home/ruigengji/MLPATH/explore && /home/ruigengji/miniforge3/envs/openmm_dev/bin/python -m pytest -q tests/test_poses.py`
Expected: PASS（18 passed）

- [ ] **Step 5: 写 commit 草稿** `.commit_drafts/<日期_时间>_explore-m3-poses.md`（message `explore: two-fragment rigid poses with keyed uniform sampling and L0 records (M3)`）。

---

### Task 4: 距离意图、二元组合与 I4 可行域

**Files:**
- Create: `explore/src/explore/intents.py`
- Modify: `explore/src/explore/probe.py`（`_feasible` 与 `constrained_scan` 调用处）
- Test: `explore/tests/test_intents.py`；`explore/tests/test_probe.py`（追加）

**Interfaces:**
- Consumes: `explore.probe.Intent`、`stream`、`fragments`、`prrs.state.encode`、`prrs.perturbations._hop_counts`、`prrs.chemistry.canonical_labels`
- Produces:
  - `IntentSpec(intent: Intent, sign: int, family: str, layer: str, scope: str)`；`intent.intent_id` 即稳定 ID
  - `distance_intents(atoms, pose_id, config, n_far: int, seed: int) -> tuple[list[IntentSpec], dict]`
  - `pair_intents(singles: list[IntentSpec], atoms, config, locality: str) -> tuple[list[IntentSpec], dict]`，`locality ∈ {"shared_atom", "adjacent", "all"}`
  - `probe._feasible(intent, q0, deltas, numbers, config) -> bool`

- [ ] **Step 1: 写失败测试**

```python
# explore/tests/test_intents.py
import numpy as np

from explore.intents import distance_intents, pair_intents
from explore.probe import Intent, _feasible
from molecules import chloride_methyl_chloride, ethanol, relabel
from prrs.config import SearchConfig

CONFIG = SearchConfig(closed_shell_only=False)


def _mapped(specs, perm):
    """Intent identity in ORIGINAL atom indices: perm[k] is the original index of new atom k."""
    out = set()
    for s in specs:
        pairs = tuple(tuple(sorted((int(perm[i]), int(perm[j])))) for i, j in s.intent.pairs)
        out.add((tuple(sorted(zip(pairs, s.intent.coefs))), s.sign, s.layer))
    return out


def _plain(specs):
    return _mapped(specs, np.arange(100))


def test_bonds_get_both_signs_and_close_intramolecular_pairs_are_excluded():
    atoms = ethanol()
    singles, report = distance_intents(atoms, "p", CONFIG, n_far=0, seed=0)
    bonded = [s for s in singles if s.layer == "bonded"]
    assert len(bonded) == 2 * 8  # ethanol has 8 bonds
    assert report["excluded"]["bend_in_disguise"] > 0
    assert all(s.family == "distance" and len(s.intent.pairs) == 1 for s in singles)


def test_cross_fragment_pairs_are_scoped_inter():
    atoms = chloride_methyl_chloride()
    singles, _ = distance_intents(atoms, "p", CONFIG, n_far=0, seed=0)
    inter = [s for s in singles if s.scope == "inter"]
    assert inter and all(5 in s.intent.pairs[0] for s in inter)


def test_far_layer_is_sampled_by_colour_class_and_counted():
    atoms = chloride_methyl_chloride()
    atoms.positions[5] = [0, 0, -9.0]  # push Cl- beyond pair_cutoff_A
    none, rep0 = distance_intents(atoms, "p", CONFIG, n_far=0, seed=0)
    some, rep1 = distance_intents(atoms, "p", CONFIG, n_far=1, seed=0)
    assert not [s for s in none if s.layer == "far"]
    far = [s for s in some if s.layer == "far"]
    assert far and rep1["far_classes_total"] >= rep1["far_classes_sampled"] == 1
    # a sampled colour class comes in whole: the three H...Cl- pairs together
    hcl = [s for s in far if {int(atoms.numbers[i]) for i in s.intent.pairs[0]} == {1, 17}]
    assert len(hcl) in (0, 6)


def test_relabelling_maps_the_intent_set():
    atoms = ethanol()
    singles, _ = distance_intents(atoms, "p", CONFIG, n_far=2, seed=0)
    moved, perm = relabel(atoms, seed=11)
    msingles, _ = distance_intents(moved, "p", CONFIG, n_far=2, seed=0)
    assert _mapped(msingles, perm) == _plain(singles)
    combos, _ = pair_intents(singles, atoms, CONFIG, "adjacent")
    mcombos, _ = pair_intents(msingles, moved, CONFIG, "adjacent")
    assert _mapped(mcombos, perm) == _plain(combos)


def test_localities_are_nested_and_ids_are_unique():
    atoms = ethanol()
    singles, _ = distance_intents(atoms, "p", CONFIG, n_far=0, seed=0)
    counts = {}
    for loc in ("shared_atom", "adjacent", "all"):
        combos, rep = pair_intents(singles, atoms, CONFIG, loc)
        ids = [c.intent.intent_id for c in combos]
        assert len(ids) == len(set(ids))
        assert all(len(c.intent.pairs) == 2 and c.family == "distance_pair" for c in combos)
        counts[loc] = len(combos)
    assert counts["shared_atom"] <= counts["adjacent"] <= counts["all"]


def test_same_sign_pair_below_overlap_floor_is_infeasible():
    numbers = np.array([6, 6, 6])
    intent = Intent("i", ((0, 1), (1, 2)), (1.0, 1.0))
    floor = 2 * CONFIG.min_distance_ratio * (2 * 0.76)
    assert not _feasible(intent, floor + 0.5, (0.0, -1.0), numbers, CONFIG)
    assert _feasible(intent, floor + 0.5, (0.0, -0.4), numbers, CONFIG)


def test_mixed_sign_pair_is_not_bounded_by_the_floor():
    numbers = np.array([6, 6, 6, 6])
    intent = Intent("i", ((0, 1), (2, 3)), (1.0, -1.0))
    assert _feasible(intent, 0.0, (0.0, -5.0, 5.0), numbers, CONFIG)
```

追加到 `explore/tests/test_probe.py`：

```python
def test_multi_pair_infeasible_target_spends_no_forward(tmp_path):
    ledger = Ledger(tmp_path / "l.jsonl", arm="D2", budget=1000)
    atoms = Atoms("C3", positions=[[0, 0, 0], [1.5, 0, 0], [3.0, 0, 0]])
    intent = Intent("i", ((0, 1), (1, 2)), (1.0, 1.0))
    settings = ScanSettings((0.0, -2.5), 1e-4, 50, 10.0)
    result = constrained_scan(atoms, intent, settings, double_well_factory, CONFIG, ledger, "t")
    assert result.status == "infeasible" and ledger.spent == 0
```

- [ ] **Step 2: 跑测试确认失败**

Run: `cd /home/ruigengji/MLPATH/explore && /home/ruigengji/miniforge3/envs/openmm_dev/bin/python -m pytest -q tests/test_intents.py tests/test_probe.py`
Expected: FAIL，`ModuleNotFoundError: explore.intents`；`_feasible` 签名不符

- [ ] **Step 3: 改 `probe._feasible`**

```python
def _feasible(intent, q0, deltas, numbers, config):
    """Design §3.2: reject an unreachable constraint before any evaluation. Every distance
    is bounded below by the existing overlap guard, r_k >= min_distance_ratio *
    (R_i + R_j), so with all coefficients positive q cannot go below sum c_k floor_k
    (all negative: above it). Mixed signs leave q unbounded both ways."""
    from ase.data import covalent_radii

    floors = [
        config.min_distance_ratio * (covalent_radii[numbers[i]] + covalent_radii[numbers[j]])
        for i, j in intent.pairs
    ]
    bound = sum(c * f for c, f in zip(intent.coefs, floors))
    if all(c > 0 for c in intent.coefs):
        return all(q0 + d > bound for d in deltas)
    if all(c < 0 for c in intent.coefs):
        return all(q0 + d < bound for d in deltas)
    return True
```

调用处改为 `if not _feasible(intent, q0, settings.deltas_A, work.numbers, config):`。原 `test_negative_target_distance_is_infeasible_before_any_forward` 仍应通过（负距离必然低于下限）。

- [ ] **Step 4: 实现 `intents.py`**

```python
"""v1 intents (DISCOVERY_LAYER_DESIGN_2026-10-13 §3.1): distance primitives and their
binary combinations. An intent says which distances move and in which sign; it never
names a product, a reference bond change or a saddle direction. Units are Å throughout
(v1 has only distances, so §3.2's per-family scale is 1 Å and nothing is mixed)."""

from __future__ import annotations

import itertools
from dataclasses import dataclass

from .poses import RULE_VERSION, fragments
from .probe import Intent
from .streams import stream


@dataclass(frozen=True)
class IntentSpec:
    intent: Intent
    sign: int
    family: str
    layer: str
    scope: str


def _ident(pose_id, family, pairs, coefs, sign):
    body = ";".join(f"{i}-{j}" for i, j in pairs)
    cs = ",".join(f"{c:+g}" for c in coefs)
    return f"{pose_id}|{family}:{body}:{cs}:{'+' if sign > 0 else '-'}:{RULE_VERSION}"


def distance_intents(atoms, pose_id, config, n_far, seed):
    from prrs.chemistry import canonical_labels
    from prrs.perturbations import _hop_counts
    from prrs.state import encode

    graph = encode(atoms, config.bond_scale)
    labels = canonical_labels(atoms.numbers, graph.edges, graph.index)
    hops = _hop_counts(graph)
    owner = {i: k for k, f in enumerate(fragments(atoms, config)) for i in f.atoms}
    near, far_classes = [], {}
    excluded = {"bend_in_disguise": 0}
    n = len(atoms)
    for i, j in itertools.combinations(range(n), 2):
        bonded = (i, j) in graph.edges
        scope = "inter" if owner[i] != owner[j] else f"intra:f{owner[i]}"
        if bonded:
            near.append(((i, j), "bonded", scope))
            continue
        sep = hops.get((i, j))
        if sep is not None and sep < config.min_pair_hops:
            excluded["bend_in_disguise"] += 1
            continue
        if atoms.get_distance(i, j) <= config.pair_cutoff_A:
            near.append(((i, j), "near", scope))
        else:
            key = tuple(sorted((str(labels[i]), str(labels[j])))) + (scope.split(":")[0],)
            far_classes.setdefault(key, []).append(((i, j), "far", scope))
    classes = sorted(far_classes)
    order = stream(seed, "far", pose_id).permutation(len(classes)) if classes else []
    chosen = [classes[k] for k in list(order)[:n_far]]
    rows = near + [m for key in chosen for m in far_classes[key]]
    specs = []
    for pair, layer, scope in rows:
        for sign in (-1, 1):
            specs.append(IntentSpec(
                Intent(_ident(pose_id, "distance", (pair,), (1.0,), sign), (pair,), (1.0,)),
                sign, "distance", layer, scope,
            ))
    report = {
        "excluded": excluded,
        "bonded_pairs": sum(1 for r in near if r[1] == "bonded"),
        "near_pairs": sum(1 for r in near if r[1] == "near"),
        "far_pairs_total": sum(len(v) for v in far_classes.values()),
        "far_classes_total": len(classes),
        "far_classes_sampled": len(chosen),
    }
    return specs, report


def _adjacent(p, q, edges):
    if set(p) & set(q):
        return True
    return any((min(a, b), max(a, b)) in edges for a in p for b in q)


def pair_intents(singles, atoms, config, locality):
    """Unordered pairs of distinct single-pair intents (bonded and near layers only),
    coefficients (1, 1) and (1, -1), both signs: 4 intents per combination. Pairs are
    put in index order so (1, -1) and (-1, 1) are not both generated."""
    from prrs.state import encode

    edges = encode(atoms, config.bond_scale).edges
    base = sorted({s.intent.pairs[0] for s in singles if s.layer in ("bonded", "near")})
    scope = {s.intent.pairs[0]: s.scope for s in singles}
    pose_id = singles[0].intent.intent_id.split("|")[0] if singles else ""
    specs, counts = [], {"shared_atom": 0, "adjacent": 0, "all": 0}
    for p, q in itertools.combinations(base, 2):
        shared = bool(set(p) & set(q))
        adj = _adjacent(p, q, edges)
        counts["all"] += 1
        counts["adjacent"] += adj
        counts["shared_atom"] += shared
        keep = {"shared_atom": shared, "adjacent": adj, "all": True}[locality]
        if not keep:
            continue
        cat = "inter" if "inter" in (scope[p], scope[q]) else (
            "intra" if scope[p] == scope[q] else "cross_intra")
        for coefs in ((1.0, 1.0), (1.0, -1.0)):
            for sign in (-1, 1):
                specs.append(IntentSpec(
                    Intent(_ident(pose_id, "distance_pair", (p, q), coefs, sign), (p, q), coefs),
                    sign, "distance_pair", locality, cat,
                ))
    return specs, {"combinations_by_locality": counts, "intents": len(specs)}
```

注意：`_hop_counts` 的键是 `(min, max)`，与 `itertools.combinations` 的 `(i, j)`（i<j）一致；`graph.edges` 也是有序对，用前先在测试里核对一次（若 `edges` 存的不是 `(min, max)`，改成 `frozenset` 比较）。

- [ ] **Step 5: 跑测试确认通过**

Run: `cd /home/ruigengji/MLPATH/explore && /home/ruigengji/miniforge3/envs/openmm_dev/bin/python -m pytest -q tests/test_intents.py tests/test_probe.py`
Expected: PASS（test_intents 7 passed；test_probe 原有全部 + 1 新增 passed）

- [ ] **Step 6: 跑本包检查**

Run: `make -C /home/ruigengji/MLPATH/explore check`
Expected: fmt-check、lint 通过；pytest 全部通过（原 31 + 本计划新增 26 左右）

- [ ] **Step 7: 写 commit 草稿** `.commit_drafts/<日期_时间>_explore-m3-intents.md`（message `explore: distance intents, binary combinations and multi-pair feasibility (M3, I4)`）。

---

### Task 5: 纯几何计数表（交用户定 L1 规模）

**Files:**
- Create: `benchmarks/explore_counts/count_l0.py`
- Output: `benchmarks/explore_counts/results/<YYYY-MM-DD>/{contact_ratios.json, counts.json, COUNTS.md}`

**Interfaces:**
- Consumes: Task 1–4 的全部公开函数；`benchmarks/b04_sn2_chloride/input.extxyz`、`benchmarks/b05_tertiary_halide/input.extxyz`、`benchmarks/carbene_ambimodal/probe_b09/input.extxyz`；RDKit（N1 几何）
- Produces: 交用户的计数表；不被代码依赖

- [ ] **Step 1: 写脚本**

脚本做四件事，全部不加载模型：

1. **现有复合物的接触比值**（给 U3 定壳层的依据）：对 b04、b05、b09 输入各算 `contact_ratio(片段 A, 片段 B)`，写 `contact_ratios.json`。
2. **N1 几何**：RDKit `Chem.AddHs(Chem.MolFromSmiles("C1=CC=CC1"))`、`"C=CC=O"` 各自 `AllChem.EmbedMolecule(mol, randomSeed=0)` + `MMFFOptimizeMolecule`，拼成一个 ASE `Atoms`（第二片段沿 x 平移 8 Å），作为计数用输入；表中注明“MMFF 几何，仅计数，不是 POLAR-1 极小”。
3. **构象计数**：每个片段 `conformer_grid(..., steps_per_turn=6, cap=10**6, seed=0)`，报 genuine 扭转数、转子数、网格大小、各状态数。
4. **姿态与意图计数**：每体系取片段 A、B 的输入构象（grid 全 0）配对，壳层用第 1 步实测的各体系比值和其 ±0.4（仅为看敏感度，不是建议值），`n_samples ∈ {8, 32}`；对前 16 个姿态逐个算 `distance_intents(n_far=0)` 与三种 `pair_intents`，报每姿态单原语数（按 bonded/near 与 intra/inter 分）、二元组合数（三种局部规则）、远层颜色类总数。总量估计 = 每姿态意图数均值 × 姿态数 × 构象配对数。
5. **L1 成本外推**：`forwards ≈ 意图数 × 每扫描点数 × 8.7`，每扫描点数取 b04 冒烟的 18、单点 8.7 forward 也来自 b04（DONE_LOG 10-13）；另给墙钟 `≈ 意图数 × 18 × 1.07 s`，换算成单张 2080 Ti 的 GPU 小时（1.07 s/点 同样是 b04 实测）。表中写明“b04（6 原子）实测外推，体系越大每点 BFGS 步数和单次 forward 都可能更慢，只作量级参考”。（用户 2026-10-22：担心主线计算太慢，定 L1 规模时要能直接看到时间。）

`COUNTS.md` 用一张表按体系列出以上数字，表头写清所用规则版本与参数。

- [ ] **Step 2: 本机跑**（纯几何；先用 `--systems b04` 跑一次看耗时，再跑全部）

Run: `cd /home/ruigengji/MLPATH && OMP_NUM_THREADS=4 /home/ruigengji/miniforge3/envs/openmm_dev/bin/python benchmarks/explore_counts/count_l0.py --out benchmarks/explore_counts/results/<日期>/ --systems b04 b05 n1 b09`
Expected: 三个输出文件；若全部体系超过 5 分钟，先停，报实测耗时再说。

- [ ] **Step 3: 更新文档**

`docs/DONE_LOG.md` 追加一节（交付、测试数、计数表路径、未做的事）；`docs/PLAN_EXPLORE_V1_2026-10-13.md` 状态行补“M3 已完成，计数表待用户定 L1 规模”。

- [ ] **Step 4: 写 commit 草稿** `.commit_drafts/<日期_时间>_explore-m3-counts.md`。

- [ ] **Step 5: 把计数表交给用户**，连同 U1–U5 中仍未定的项，请用户冻结：壳层、`n_samples`、构象上限、局部规则、`n_far`。M5（D2 编排）计划要等这一步。

---

## 自查记录

- 设计覆盖：§2.1 来源构象（Task 1 网格、Task 2 资格与去重）；§2.2 姿态（Task 3，均匀方向与 Haar 朝向、原始姿态单存、重标号不变、单分子不造姿态）；§2.3 L0 字段（Task 3 `write_l0`）；§3.1 距离与二元组合、排除规则与远层抽样（Task 4）；§3.2 稳定 ID、按键随机流、求值前拒绝不可行约束（Task 1/3/4）；§1.1 “pilot 前纯几何计数”（Task 5）。姿态的模型松弛与“完全无约束极小”工件按“明确不做”推迟到 M5/M10，理由是它要花模型成本，应随编排一起入账。
- 无占位：各步给了代码与命令；Task 5 的脚本按步骤 1 的五项写，输出字段已列。
- 类型一致：`Conformer`、`Pose`、`IntentSpec`、`_feasible(intent, q0, deltas, numbers, config)` 在各 Task 用法一致。
- Review Focus 五条都落到了具体测试：碰撞点（Task 1）、相同片段（Task 3）、单/双原子片段（Task 1）、三片段（Task 3）、多对意图不可行（Task 4）。
