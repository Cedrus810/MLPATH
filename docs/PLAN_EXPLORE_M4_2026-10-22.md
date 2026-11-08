# 发现层 v1 · M4 C2 认证 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 建共同认证后端 L4：从一个鞍点初猜出发，得到 C2 / duplicate_C2 / unresolved / rejected 之一，全程计入账本；先补完 M2 遗留的账本项 I1–I3。

**Architecture:** 新增 `explore/src/explore/certify.py`，只组合 PRRS 现有原语：`follow_min_mode`（精修到一阶驻点）→ `confirm_minimum`（独立完整 Hessian，要求恰一个负模，并给出不稳定模）→ 两侧 `descend_saddle`（各自以 `confirm_minimum` 确认极小）→ `classify` / `chemical_key`（身份）→ 与已认证项比较去重。PRRS 源码不动。测试用 PRRS 现有的 `ProtonRelay`（不同端点）和一个新写的旋转不变双鞍点势（同端点两鞍点）。

**Tech Stack:** Python 3.12（`openmm_dev`，生产环境，不装不升）、ASE 3.29.0、numpy、pytest；ruff 在 `.venv-tools`。

**Spec:** [评估协议](EVALUATION_PROTOCOL_DESIGN_2026-10-13.md) §1.1 C2 定义、§2.3 阶段名、§3.2 第 2–5 步；[发现层设计](DISCOVERY_LAYER_DESIGN_2026-10-13.md) §1（L4 是共同后端，不能把 DB 的搜索搬进 D2）、§5 L4 行、§9 失败码；上位计划 [PLAN_EXPLORE_V1](PLAN_EXPLORE_V1_2026-10-13.md) M4 行。

状态：**2026-10-22 已施工完成**（用户指示“继续往下”，Subagent-driven；裁决 M4-R1…R6 与 b04 冒烟结果见 DONE_LOG 末节）。

## 必读（开工前全部读完）

| 文件 | 为什么要读 |
| --- | --- |
| `CLAUDE.md`（项目根目录） | 设计稿表；`openmm_dev` 是生产环境 |
| `docs/EVALUATION_PROTOCOL_DESIGN_2026-10-13.md` §1.1、§2.3、§3.2 | C2 的定义：无偏 PES 驻点收敛、恰一显著负模、两侧下降到经确认的端点、事件与身份闭合；逐步过程不虚构一个直接鞍点；同端点多鞍点不能靠端点对去重 |
| `docs/DISCOVERY_LAYER_DESIGN_2026-10-13.md` §1、§5、§9 | L4 是各臂共用的后端；失败码必须分开 |
| `src/prrs/runner.py`：`follow_min_mode`、`confirm_minimum`、`descend_saddle`、`hessian_spectrum` 的 docstring | 复用原语及其返回约定 |
| `src/prrs/search.py` 约 1400–1560 行（`register_saddle`、`connect_saddle`） | PRRS 自己怎样把这几个原语串起来，照着做 |
| `explore/src/explore/ledger.py`、`probe.py` | 账本段、计数包装；I1–I3 要改的地方 |

## 与上位计划不同的地方

**与你说的不同：** PLAN_EXPLORE_V1 的 M4 出口写的是“解析双鞍点台两个鞍点都认证为 C2”。**原因：** 冻结的 `benchmarks/pars_two_saddles` 是外场势（它的 README 说明“刚体运动移除已关闭，因为势是外场”），而 PRRS 的 Hessian 原语总是投影掉六个刚体模，PRRS 里也没有关掉这一步的开关；直接拿来用会把真实模式当成刚体模删掉。所以 M4 在 `explore/tests/` 新写一个**旋转平移不变**的同构双鞍点势（同一套 `E=(x²+y²−1)²+a·y²+b·y³+z_k·z²/2`，但 x、y、z 是 H 在三个锚原子构成的分子坐标系里的坐标），出口改为在它上面两个鞍点都认证为 C2 且互不判重。冻结台原件不动。

## Global Constraints

- 科学环境只用 `/home/ruigengji/miniforge3/envs/openmm_dev/bin/python`；不 pip/conda 安装、卸载、升级。
- PRRS、PARS 源码不动；`explore` 单向 import `prrs`。
- 不新设阈值：负模判据用 `SearchConfig.minimum_check_eigenvalue_tol`，去重用 `basin_rmsd_A` / `basin_energy_eV`，接力次数用 `irc_relay_hops`，上升上限由调用方显式传入（不在 certify 里写默认值）。
- 阶段名用评估协议 §2.3 的：`saddle_refine`、`hessian`、`descent`；每个 forward 都在账本段里。
- 失败不压缩：每种失败有自己的状态，保留已得几何、诊断与已付成本。
- 证据等级不升级：约束扫描峰、负曲率、NEB 状态都不能直接当 C2；有接力（relay）的下降不是直接连接。
- 输出不写 `/tmp`；真实模型冒烟结果放 `benchmarks/explore_smoke/results/<日期>/`。
- NFS 上不跑 git；每个 Task 结束写 `.commit_drafts/<日期_时间>_<主题>.md`，**整个里程碑合成一份草稿**（M3 的教训：逐 Task 草稿会互相重叠）。
- 测试：只跑 `make -C explore check`。

## Review Focus

1. **下降经接力才到极小**（`descend_saddle` 返回的 `relays` 非空）：不能记 C2，记 `unresolved:relayed`，并保留接力记录。Task 3 测试固定。
2. **两侧落到同一个盆**：记 `unresolved:both_sides_same_basin`，单报，不当连接。Task 3 测试固定。
3. **同一对端点的第二个鞍点**：几何不同就不是重复，必须各自 C2。Task 3 在双鞍点势上固定。
4. **预算在认证中途耗尽**：状态 `budget_exhausted`，已付成本留在账本，已得几何保留。Task 3 测试固定。
5. **身份算不出**（`IdentityUnavailable`，例如自同构预算溢出）：记 `unresolved:identity`，不猜。Task 3 测试固定。

---

## 文件结构

```
explore/src/explore/
  ledger.py      （改）I2：段内记缓存命中数
  probe.py       （改）I1：attempt / reserve / artifact 透传；I3：门拒与预算耗尽的失败码写进账本，被拒几何随结果返回
  certify.py     （新）C2 认证
explore/tests/
  twosaddle.py   （新）旋转不变的同端点双鞍点势
  test_ledger.py （加 I2 测试）
  test_probe.py  （加 I1、I3 测试）
  test_twosaddle.py（新）
  test_certify.py（新）
benchmarks/explore_smoke/
  b04_certify.py （新）b04 真实模型冒烟：从 10-18 扫描峰点认证
```

---

### Task 1: 账本与探针遗留项（I1–I3）

**Files:**
- Modify: `explore/src/explore/ledger.py`（`Segment`、`segment()` 的 close 事件、`instrument.get_property`）
- Modify: `explore/src/explore/probe.py`（`constrained_scan`、`reverse_check`、`release`）
- Test: `explore/tests/test_ledger.py`、`explore/tests/test_probe.py`

**Interfaces:**
- Produces:
  - 账本 close 事件新增字段 `cache_hits`（段内：请求未触发 `calculate` 的次数）。
  - `constrained_scan(atoms, intent, settings, factory, config, ledger, task, attempt=0, reserve=None, artifact=None)`；`release(..., task, seed, attempt=0, reserve=None, artifact=None)`；`reverse_check(..., task, attempt=0, reserve=None, artifact=None)`。`artifact` 是调用方给的工件路径字符串，原样写进段的 `artifact` 字段。
  - `ScanResult` 新增 `rejected_positions: np.ndarray | None`（门拒时被拒几何，取自 `GateRejected.atoms`）；门拒时账本段 `outcome = f"gate_rejected:{code}"`，预算耗尽时 `outcome = "budget_exhausted"`（目前这两种情况段的 outcome 是 `error:<异常名>`）。

- [ ] **Step 1: 写失败测试**

```python
# test_ledger.py 追加
def test_cache_hits_are_counted_per_segment(tmp_path):
    ledger = Ledger(tmp_path / "l.jsonl", arm="t", budget=100)
    atoms = dimer(1.0)
    atoms.calc = instrument(LennardJones(), ledger)
    with ledger.segment(task="a", stage="probe"):
        atoms.get_potential_energy()   # forward
        atoms.get_potential_energy()   # cache hit
        atoms.get_forces()             # cache hit (LJ computes both)
    ledger.close()
    close = [e for e in lines(tmp_path / "l.jsonl") if e["event"] == "close"][0]
    assert close["forwards"] == 1 and close["cache_hits"] == 2
```

（`dimer`、`lines` 已在 test_ledger.py 里；如 LJ 一次 calculate 不同时给出 forces，按实测改断言并在报告里写明。）

```python
# test_probe.py 追加
def test_attempt_reserve_and_artifact_reach_the_ledger(tmp_path):
    ledger = Ledger(tmp_path / "l.jsonl", arm="D2", budget=1000)
    atoms = Atoms("C2", positions=[[0, 0, 0], [1.2, 0, 0]])
    intent = Intent("i", ((0, 1),), (1.0,))
    settings = ScanSettings((0.0, 0.1), 1e-4, 50, 10.0)
    constrained_scan(atoms, intent, settings, counted_factory(double_well_factory, ledger),
                     CONFIG, ledger, "t", attempt=3, reserve=40, artifact="scan/t.json")
    ledger.close()
    ev = [json.loads(x) for x in (tmp_path / "l.jsonl").read_text().splitlines()]
    opened = next(e for e in ev if e["event"] == "open" and e["task"] == "t")
    closed = next(e for e in ev if e["event"] == "close" and e["task"] == "t")
    assert opened["attempt"] == 3 and opened["reserve"] == 40
    assert closed["artifact"] == "scan/t.json"


def test_gate_rejection_is_named_in_the_ledger_and_geometry_kept(tmp_path, monkeypatch):
    import explore.probe as probe
    monkeypatch.setattr(probe, "_feasible", lambda *a, **k: True)
    ledger = Ledger(tmp_path / "l.jsonl", arm="D2", budget=1000)
    atoms = Atoms("C2", positions=[[0, 0, 0], [1.2, 0, 0]])
    intent = Intent("i", ((0, 1),), (1.0,))
    settings = ScanSettings(tuple(round(-0.1 * k, 10) for k in range(12)), 1e-4, 50, 10.0)
    result = constrained_scan(atoms, intent, settings,
                              counted_factory(double_well_factory, ledger), CONFIG, ledger, "t")
    ledger.close()
    assert result.status == "gate_rejected:overlap"
    assert result.rejected_positions is not None and result.rejected_positions.shape == (2, 3)
    closed = [json.loads(x) for x in (tmp_path / "l.jsonl").read_text().splitlines()
              if json.loads(x)["event"] == "close" and json.loads(x)["task"] == "t"]
    assert closed[0]["outcome"] == "gate_rejected:overlap"
```

（再加一条预算耗尽的同类测试：`budget=5`，断言段 outcome 为 `"budget_exhausted"`、`result.status == "budget_exhausted"`、已得点保留。）

- [ ] **Step 2: 跑测试确认失败**

Run: `cd /home/ruigengji/MLPATH/explore && /home/ruigengji/miniforge3/envs/openmm_dev/bin/python -m pytest -q tests/test_ledger.py tests/test_probe.py`
Expected: FAIL（`KeyError: 'cache_hits'`；`constrained_scan() got an unexpected keyword argument 'attempt'`；outcome 是 `error:GateRejected`）

- [ ] **Step 3: 实现**

`ledger.py`：

```python
class Segment:
    def __init__(self, task, attempt, stage, reserve):
        ...
        self.cache_hits = 0
```

close 事件加 `"cache_hits": seg.cache_hits`（`close()` 里 unattributed 的 close 也加 `"cache_hits": 0`）。`instrument.get_property` 改为：

```python
    def get_property(name, atoms=None, allow_calculation=True):
        ledger.note_request()
        before = ledger.spent
        value = inner_get_property(name, atoms, allow_calculation)
        if ledger.spent == before:
            ledger.note_cache_hit()
        return value
```

并加 `Ledger.note_cache_hit()`（段内 `cache_hits += 1`）。

`probe.py`：三个函数加 `attempt=0, reserve=None, artifact=None`，传给 `ledger.segment(task=..., stage=..., attempt=attempt, reserve=reserve)`，进段后 `seg.artifact = artifact`。`constrained_scan` 把 `except BudgetExhausted` / `except GateRejected` 移到 `with ledger.segment(...)` **内部**，在里面先设 `seg.outcome`（`"budget_exhausted"` / `f"gate_rejected:{exc.code}"`），门拒时 `rejected_positions = exc.atoms.positions.copy()`（`GateRejected` 的第三个参数就是 atoms，见 `src/prrs/reliability.py`）；`ScanResult` 加字段 `rejected_positions: np.ndarray | None = None`。`release` 的 quench 与 minimum_check 两段同样在段内设 outcome。

- [ ] **Step 4: 跑测试确认通过**，再跑 `make -C /home/ruigengji/MLPATH/explore check`。

- [ ] **Step 5:** 把本 Task 的文件记进本里程碑的合并草稿（见 Global Constraints）。

---

### Task 2: 旋转不变的同端点双鞍点势

**Files:**
- Create: `explore/tests/twosaddle.py`
- Test: `explore/tests/test_twosaddle.py`

**Interfaces:**
- Produces:
  - `class TwoSaddleInvariant(ase.calculators.calculator.Calculator)`，参数 `a=0.8, b=0.08, z_k=2.0, frame_k=10.0`（与冻结台 `protocol.json` 的 a、b、z_k 相同）
  - `two_saddle_factory()`、`two_saddle_atoms(x, y, z=0.0) -> Atoms`
  - 常量 `SADDLE_Y = (-0.8051774, 0.7451774)`（冻结台 README 给的两个鞍点 y 值，x=0）、极小在 `x=±1, y=0`

设计：4 个原子，锚 O0、O1、O2 组成一个刚性三角（三条边各一个谐振项，力常数 `frame_k`，平衡边长 2.5 Å，元素用 O 让键图判据不把锚连成键），H3 在锚坐标系里的坐标为
`e1 = (O1−O0)/|O1−O0|`，`e2 = (O2−O0) 去掉 e1 分量后归一`，`e3 = e1×e2`，原点 = 三锚质心，`(x, y, z) = ((H−c)·e1, (H−c)·e2, (H−c)·e3)`；能量 `E = (x²+y²−1)² + a·y² + b·y³ + z_k·z²/2 + 锚三角谐振项`。力用能量的中心差分（步长 1e-6 Å），实现简单、精度 ~1e-10 eV/Å，足够做测试（在 docstring 写明这是测试势，不是化学）。

几何要求：极小与鞍点处 H 到每个锚的距离都要在 `encode` 的成键阈值之外或之内**一致**，使两个极小的键图相同、两个鞍点的端点是同一对（同映射端点）。实现时先在测试里打印 H–O 距离，若某点落在阈值附近，调大锚三角边长（例如 3.0 Å），并在报告里写明最终参数。

- [ ] **Step 1: 写失败测试**

```python
# explore/tests/test_twosaddle.py
import numpy as np

from prrs.config import SearchConfig
from prrs.runner import confirm_minimum, hessian_spectrum
from prrs.io import snapshot
from prrs.reliability import GuardedCalculator
from twosaddle import SADDLE_Y, two_saddle_atoms, two_saddle_factory

CONFIG = SearchConfig(closed_shell_only=False)


def _snap(x, y):
    atoms = two_saddle_atoms(x, y)
    return snapshot(atoms, GuardedCalculator(two_saddle_factory(), CONFIG))


def test_energy_is_invariant_under_rigid_motion():
    from scipy.spatial.transform import Rotation

    atoms = two_saddle_atoms(0.3, 0.4, 0.1)
    atoms.calc = two_saddle_factory()
    e0 = atoms.get_potential_energy()
    moved = atoms.copy()
    moved.positions = Rotation.from_rotvec([0.3, -1.1, 0.7]).apply(moved.positions) + [1, 2, 3]
    moved.calc = two_saddle_factory()
    assert abs(moved.get_potential_energy() - e0) < 1e-10


def test_both_saddles_have_exactly_one_negative_mode():
    for y in SADDLE_Y:
        eig = hessian_spectrum(_snap(0.0, y), two_saddle_factory, CONFIG)
        assert np.sum(eig < -CONFIG.minimum_check_eigenvalue_tol) == 1


def test_minima_are_confirmed_and_share_one_bond_graph():
    from prrs.state import encode

    graphs = []
    for x in (-1.0, 1.0):
        ok, _ = confirm_minimum(_snap(x, 0.0), two_saddle_factory, CONFIG, seed=0)
        assert ok
        graphs.append(encode(two_saddle_atoms(x, 0.0), CONFIG.bond_scale).edges)
    assert graphs[0] == graphs[1]
```

- [ ] **Step 2: 跑测试确认失败**（`ModuleNotFoundError: twosaddle`）

- [ ] **Step 3: 实现 `twosaddle.py`**（按上面的设计；`two_saddle_atoms(x, y, z)` 先放好锚三角，再把 H 放到局部坐标 (x, y, z)）。注意：在中心差分里给 `atoms.positions` 加扰动时要用副本，不要改调用方的 atoms。

- [ ] **Step 4: 跑测试确认通过**；若 `confirm_minimum` 因 `trivial_mode_floor` 失败，说明力的数值精度不够，把差分步长调到 1e-5 再试，并在报告写明。

- [ ] **Step 5:** 记入合并草稿。

---

### Task 3: `certify.py`

**Files:**
- Create: `explore/src/explore/certify.py`
- Test: `explore/tests/test_certify.py`

**Interfaces:**
- Consumes: `prrs.runner.follow_min_mode/confirm_minimum/descend_saddle`、`prrs.state.classify/same_basin`、`prrs.chemistry.chemical_key/automorphisms/symmetric_rmsd/IdentityUnavailable`、`prrs.io.snapshot`、`prrs.reliability.GuardedCalculator/GateRejected`、`explore.ledger.Ledger/BudgetExhausted`
- Produces:
  - `@dataclass Certification`：`status: str`、`saddle_positions: np.ndarray | None`、`saddle_energy_eV: float | None`、`unstable_mode: list | None`、`hessian: dict | None`（`confirm_minimum` 的诊断，JSON 化）、`sides: list[dict]`（每侧 `sign, status, positions, energy_eV, relays, chemical_key`）、`event: dict | None`（`classify(minus, plus)` 的 label 与 fingerprint）、`duplicate_of: str | None`、`diagnostics: dict`、`spent: int`（本次认证花掉的 forward）
  - `certify(guess, factory, config, ledger, task, seed, ceiling_eV, reference_eV, direction=None, known=(), attempt=0) -> Certification`；`known` 是已认证的 `(id, Certification)` 序列，用于判重
  - 状态集合：`"c2"`、`"duplicate_c2"`、`"unresolved:saddle_not_found"`、`"rejected:order_<n>"`、`"unresolved:descent_failed"`、`"unresolved:relayed"`、`"unresolved:both_sides_same_basin"`、`"unresolved:identity"`、`"budget_exhausted"`、`"gate_rejected:<code>"`

流程（每步一个账本段，阶段名照协议）：

1. `saddle_refine`：`follow_min_mode(guess, factory, config, direction=direction, seed=seed, ceiling=ceiling_eV, reference=reference_eV)`。返回 `None` → `unresolved:saddle_not_found`，诊断原样保留。
2. `hessian`：对返回结构 `snapshot` 后调 `confirm_minimum(saddle, factory, config, seed)`。它在鞍点上返回 `(False, checked)`，`checked["saddle_order"]` 是负模数、`checked["unstable_mode"]` 是模向量。`saddle_order != 1` → `rejected:order_<n>`；`checked["reason"]` 以 `gate_` 开头 → `gate_rejected:<code>`；`saddle_order is None` → `unresolved:saddle_not_found`。**不要**用 `follow_min_mode` 自己的 `unstable_mode` 代替这一步：协议要求独立的完整 Hessian。
3. `descent`：`sign` 取 `+1`、`-1` 各一次 `descend_saddle(saddle, mode, factory, config, sign, seed=seed + (1 if sign > 0 else 2))`。任一侧 endpoint 为 `None` → `unresolved:descent_failed`；任一侧 `relays` 非空 → `unresolved:relayed`（接力是更弱的证据，协议 §3.2 第 4 步：逐步过程认证每条局部边，不虚构直接鞍点）。
4. 身份：两侧 `same_basin(minus, plus, config)` 为真 → `unresolved:both_sides_same_basin`；否则 `classify(minus, plus, config)` 得事件，两侧各算 `chemical_key(endpoint, config.bond_scale, config.active_atoms, charge_sensitive=config.charge_sensitive, charge=config.total_charge, multiplicity=config.multiplicity)`（实现前核对 `chemical_key` 的返回结构，取它的身份字段）；`IdentityUnavailable` → `unresolved:identity`。退化反应（两侧同 chemical_key 但键图不同）允许，记 C2。
5. 判重：对 `known` 里每个 C2，若两侧 chemical_key 无序对相同、鞍点能量差 ≤ `basin_energy_eV`、并且鞍点几何在自同构下的 `symmetric_rmsd ≤ basin_rmsd_A` → `duplicate_c2`，`duplicate_of` 填其 id。只有端点对相同而几何不同 → 不是重复。
6. `BudgetExhausted` 在任一步 → `budget_exhausted`，已得字段保留；`GateRejected` → `gate_rejected:<code>`。`spent` = 结束时 `ledger.spent` − 开始时。

- [ ] **Step 1: 写失败测试**

```python
# explore/tests/test_certify.py
import numpy as np
import pytest

from explore.certify import certify
from explore.ledger import Ledger, counted_factory
from prrs.calculators import proton_relay_atoms, proton_relay_factory
from prrs.config import SearchConfig
from twosaddle import SADDLE_Y, two_saddle_atoms, two_saddle_factory

CONFIG = SearchConfig(closed_shell_only=False)


def _run(tmp_path, guess, factory, budget=200000, known=(), name="c"):
    ledger = Ledger(tmp_path / f"{name}.jsonl", arm="L4", budget=budget)
    guess = guess.copy()
    guess.calc = factory()
    e0 = guess.get_potential_energy()
    cert = certify(guess, counted_factory(factory, ledger), CONFIG, ledger, task=name,
                   seed=0, ceiling_eV=e0 + 5.0, reference_eV=e0, known=known)
    ledger.close()
    return cert, ledger


def test_both_saddles_of_one_endpoint_pair_are_separate_c2(tmp_path):
    certs = []
    for k, y in enumerate(SADDLE_Y):
        guess = two_saddle_atoms(0.05, 0.9 * y)
        cert, _ = _run(tmp_path, guess, two_saddle_factory,
                       known=[(f"s{j}", c) for j, c in enumerate(certs)], name=f"s{k}")
        assert cert.status == "c2", cert.diagnostics
        certs.append(cert)
    a, b = certs
    assert {s["chemical_key"] for s in a.sides} == {s["chemical_key"] for s in b.sides}
    assert abs(a.saddle_positions[3, 1] - b.saddle_positions[3, 1]) > 0.5  # different saddles


def test_the_same_saddle_twice_is_a_duplicate(tmp_path):
    guess = two_saddle_atoms(0.05, 0.9 * SADDLE_Y[1])
    first, _ = _run(tmp_path, guess, two_saddle_factory, name="a")
    second, _ = _run(tmp_path, guess, two_saddle_factory, known=[("a", first)], name="b")
    assert second.status == "duplicate_c2" and second.duplicate_of == "a"


def test_proton_hop_is_c2_with_a_bond_change(tmp_path):
    left, right = proton_relay_atoms((0, 0, 0)), proton_relay_atoms((1, 0, 0))
    guess = left.copy()
    guess.positions = 0.5 * (left.positions + right.positions)
    cert, ledger = _run(tmp_path, guess, proton_relay_factory, name="hop")
    assert cert.status == "c2", cert.diagnostics
    assert cert.event["fingerprint"]["broken"] and cert.event["fingerprint"]["formed"]
    assert cert.spent == ledger.spent > 0


def test_relayed_descent_is_not_c2(tmp_path, monkeypatch):
    import prrs.runner as runner
    real = runner.descend_saddle

    def relayed(*args, **kwargs):
        endpoint, checked, relays = real(*args, **kwargs)
        return endpoint, checked, relays + [{"hop": 0, "saddle_order": 1}]

    monkeypatch.setattr(runner, "descend_saddle", relayed)
    cert, _ = _run(tmp_path, two_saddle_atoms(0.05, 0.9 * SADDLE_Y[1]), two_saddle_factory)
    assert cert.status == "unresolved:relayed"
    assert any(side["relays"] for side in cert.sides)


def test_both_sides_in_one_basin_is_reported_not_counted(tmp_path, monkeypatch):
    import prrs.state as state
    monkeypatch.setattr(state, "same_basin", lambda a, b, c: True)
    cert, _ = _run(tmp_path, two_saddle_atoms(0.05, 0.9 * SADDLE_Y[1]), two_saddle_factory)
    assert cert.status == "unresolved:both_sides_same_basin"


def test_identity_failure_is_unresolved(tmp_path, monkeypatch):
    import prrs.chemistry as chemistry

    def boom(*a, **k):
        raise chemistry.IdentityUnavailable("automorphism_budget")

    monkeypatch.setattr(chemistry, "chemical_key", boom)
    cert, _ = _run(tmp_path, two_saddle_atoms(0.05, 0.9 * SADDLE_Y[1]), two_saddle_factory)
    assert cert.status == "unresolved:identity"


def test_budget_exhaustion_keeps_what_was_paid(tmp_path):
    cert, ledger = _run(tmp_path, two_saddle_atoms(0.05, 0.9 * SADDLE_Y[1]),
                        two_saddle_factory, budget=30)
    assert cert.status == "budget_exhausted"
    assert cert.spent == ledger.spent and ledger.spent <= 30


def test_ledger_stages_follow_the_protocol(tmp_path):
    import json

    _run(tmp_path, two_saddle_atoms(0.05, 0.9 * SADDLE_Y[1]), two_saddle_factory, name="st")
    stages = {json.loads(x)["stage"] for x in (tmp_path / "st.jsonl").read_text().splitlines()}
    assert {"saddle_refine", "hessian", "descent"} <= stages
```

（`IdentityUnavailable` 的构造参数以 `src/prrs/chemistry.py` 实际签名为准；monkeypatch 的目标要是 `certify.py` 调用时真正查找的那个名字——`certify.py` 必须在函数体内 `from prrs.state import ...` / `from prrs.chemistry import ...`，与 `poses.py` 的做法一致，否则补丁不生效。初猜离鞍点太远导致 `follow_min_mode` 走到别处时，把 0.9 改成 0.97 并在报告写明。）

- [ ] **Step 2: 跑测试确认失败**（`ModuleNotFoundError: explore.certify`）

- [ ] **Step 3: 实现 `certify.py`**（按上面的 6 步流程；诊断里的 numpy 数组转成 list，复用 `poses._jsonable`）。

- [ ] **Step 4: 跑测试确认通过**，再跑 `make -C /home/ruigengji/MLPATH/explore check`。

- [ ] **Step 5:** 记入合并草稿。

---

### Task 4: b04 真实模型冒烟（认证 10-18 扫描峰点）

**Files:**
- Create: `benchmarks/explore_smoke/b04_certify.py`
- Output: `benchmarks/explore_smoke/results/<日期>/certify_b04/`

脚本：读 `benchmarks/explore_smoke/results/2026-10-18/run/fine_fwd.extxyz` 中 `q_target == 2.28` 的那一帧作初猜（10-22 已确认它与已认证 TS 的 RMSD 为 0.004 Å），按 b04 的电荷 −1 / 多重度 1 盖章，`polar1` 计数工厂，`ceiling_eV = 初猜能量 + 1.0`、`reference_eV = 初猜能量`（这两个是本次冒烟的显式参数，不是默认值），跑一次 `certify`，写 `certify.json`（状态、鞍点能量、虚频、两侧端点能量与 chemical_key、事件、`spent`、各段 forward 与墙钟）和 `saddle.extxyz`、`minus.extxyz`、`plus.extxyz`，再与 `runs/b04_probe2_openceiling/ts_candidates/ts0000.extxyz` 做只读 RMSD 对照（复用 `compare_ts.py` 的 kabsch 函数）。

- [ ] **Step 1:** 写脚本，ruff format + lint。
- [ ] **Step 2:** 冒烟：先在 183 上查负载（`nvidia-smi`、`ps`），空闲就在 183 的 GPU 上按全局规则第 3 章的标准写法启动（`setsid nohup … & echo $!` 记 PID 到结果目录的 `pids.txt`），`OMP_NUM_THREADS=MKL_NUM_THREADS=4`；预计几分钟（6 原子；依据：10-18 四段扫描 698 forward 共约 2.5 分钟）。超过 10 分钟没结束就按 PID 核对后停掉并报告。结果在 183 上读。
- [ ] **Step 3:** 期望：`c2`，两侧端点为 SN2 的两个离子–偶极复合物，事件为断 C0–Cl1、成 C0–Cl5（或反向），鞍点与 ts0000 的 RMSD 在 0.01 Å 量级。结果与期望不符就如实报告，不调参重跑。
- [ ] **Step 4:** 记入合并草稿；`docs/DONE_LOG.md` 加 M4 一节，`docs/PLAN_EXPLORE_V1_2026-10-13.md` 状态行更新。

---

## 自查记录

- 协议覆盖：C2 四要素（驻点收敛 = `follow_min_mode` 收敛到一阶驻点；恰一显著负模 = 独立 `confirm_minimum`；两侧下降到确认端点 = `descend_saddle` 内的 `confirm_minimum`；事件与身份闭合 = `classify` + `chemical_key`）各有一步；§3.2 第 4 步的“逐步过程不虚构直接鞍点”= `unresolved:relayed`；第 5 步“只看端点对不能去掉第二鞍点”= 判重要求几何也相同；§2.3 阶段名 = 三个段名。
- 与上位计划的偏离已在前面单列（解析台换成旋转不变版本）。
- 类型一致：`Certification` 字段在 Task 3 与 Task 4 用法一致；`certify(..., ceiling_eV, reference_eV, ...)` 由调用方显式给。
- Review Focus 五条各对应 Task 3 的一个测试。
