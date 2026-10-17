# 发现层 v1 与评估量尺施工计划 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 按两份已审定设计稿的 v1 范围，建出“计数账本 + 约束松弛探针 + 姿态/意图 + C2 认证 + 四臂 + 参照 + 零号批”的最小闭环，回答“同预算下能否多找到并保住经认证的鞍点”。

**Architecture:** 新建独立包 `explore/`（与 `mature/` 同形态：自带 pyproject / Makefile / tests，单向 `import prrs`，`prrs` 与 `pars` 都不 import 它）。PRRS 源码不动，因此 PRRS 指纹与 digest 不变；动力学、淬火、Hessian、爬鞍点、两侧下降全部复用 `prrs.runner` 现有原语。PARS 保真修订在 `pars/` 内单独做。

**Tech Stack:** Python 3.12（`openmm_dev`，生产环境，不装不升）、ASE 3.29.0（`FixInternals` bondcombos、`ase.mep.dimer`）、numpy、pytest；ruff 在 `.venv-tools`。

**Spec:** [评估协议设计](EVALUATION_PROTOCOL_DESIGN_2026-10-13.md)（按其“v1 范围”）、[发现层设计](DISCOVERY_LAYER_DESIGN_2026-10-13.md)（按其 §1.1）。两份均为用户 2026-10-13 审定。

状态：**用户 2026-10-13 审定并选择 Subagent-driven 执行；M1–M2（Task 1–6）已完成，b04 冒烟跑通（未做验收），见 DONE_LOG 末节。** M3 起开工前另写逐步计划给用户审。执行中的偏离：约束改为自写单标量投影 `FixBondCombo`（ASE `FixInternals` 在线性几何上报错）；冒烟改在 181 上跑（183 的 GPU 被他人占用）。

## 必读（开工前全部读完）

| 文件 | 为什么要读 |
| --- | --- |
| `CLAUDE.md`（项目根目录） | 设计稿表、`openmm_dev` 是生产环境、183 的结果要在 183 上读 |
| `docs/EVALUATION_PROTOCOL_DESIGN_2026-10-13.md` | 量尺：C1/C2 定义、成本口径（§2）、v1 范围 |
| `docs/DISCOVERY_LAYER_DESIGN_2026-10-13.md` | 被量对象：§1.1 v1 范围、§4 约束松弛、§4.3 越脊与跳支 |
| `docs/PROJECT_DIRECTION_REVIEW_2026-10-13.md` §8 | 今天每条用户决定的原话记录 |
| `mature/Makefile`、`mature/pyproject.toml` | 新包照抄的形态 |
| `src/prrs/runner.py`：`quench`、`confirm_minimum`、`follow_min_mode`、`descend_saddle`、`hessian_spectrum` 的 docstring | 复用的原语及其失败语义 |

## Global Constraints

- 科学环境只用 `/home/ruigengji/miniforge3/envs/openmm_dev/bin/python`；不 pip/conda 安装、卸载、升级。
- `prrs` 与 `pars` 不得 import `explore`；`explore` 可以单向 import `prrs`；`pars` 仍禁止 `import prrs`。
- 主预算单位：**一个结构在一个模型上真正发起一次 forward 记 1**；ASE 缓存命中不计；发起后失败或被拒仍计 1（协议 §2.1）。
- 恢复不回滚已付成本；崩溃的段按其预留上界计入并标记 unreconciled（协议 §2.3）。
- 不用 hash 判断缓存或续跑能否复用（全局规则第 8 章第 4 条）。
- 不新增会拦住运行或结果的门；协议明文规定的预算许可除外。数值阈值（能量噪声、ΔE 上限等）从实测噪声推出，测试里显式传入，不写默认魔数。
- 输出不写 `/tmp`；实验脚本放 `benchmarks/explore_<主题>/`，结果放 `benchmarks/explore_<主题>/results/<YYYY-MM-DD>/`。
- 在 NFS 工作区不跑 git；每个里程碑结束写 `.commit_drafts/<日期_时间>_<主题>.md`，由用户提交。
- 超过 5 分钟、远程、GPU 的作业先报计划表等用户说“跑”；冒烟（几分钟就停）不需批准。

## Review Focus

1. **扫描中途触发可靠性门（重叠、力上限）**：应当以 `gate_rejected:<code>` 结束本次扫描并保留已得点，而不是抛出到上层。由 Task 3 的测试固定。
2. **预算在扫描中途耗尽**：已得点保留，状态 `budget_exhausted`，账本里该段 close 行如实记录 forwards。由 Task 3 固定。
3. **约束目标不可达（原子对目标距离为负或坐标奇异）**：应在任何求值前拒绝，状态 `infeasible`，0 forward。由 Task 3 固定。
4. **同一 calculator 被计数两次**（工厂产物再被包一层）：必须报错，不能静默双计。由 Task 2 固定。
5. **崩溃后重开账本**：已关闭段照常计入，未关闭段按预留上界计入并列入 `unreconciled`。由 Task 1 固定。

---

## 文件结构

```
explore/
  pyproject.toml          包元数据；pytest pythonpath = ["src", "../src"]
  Makefile                fmt / fmt-check / lint / test / check（照抄 mature）
  README.md               一段话说明 + 指向两份设计稿
  LICENSE                 复制 mature/LICENSE
  src/explore/
    __init__.py
    ledger.py             M1：账本、预算许可、预留、重放；calculator 计数包装
    probe.py              M2：Intent、约束松弛扫描、扫描判读、释放
    poses.py              M3：来源构象 + 双分子姿态库
    intents.py            M3：通用意图生成（距离 + 二元组合）
    certify.py            M4：L4 C2 认证（组合 prrs 原语）
    d2.py                 M5：方案 A 的 v1 编排（固定分级、两账户、均匀抽查、续跑）
    arms.py               M6：D0 / D1 / DB 适配到同一 L0 与账本
    score.py              M10：R(B)、首次命中成本、错误接受与未决
  tests/
    conftest.py
    test_ledger.py        M1
    test_probe.py         M2
benchmarks/explore_smoke/   M2 真实模型冒烟脚本与结果
benchmarks/explore_refs/    M7 参照认证脚本与参照注册表
```

## 里程碑总表

本计划**只把 M1–M2 写到逐步粒度**。M3 以后要依赖 M2 冒烟的实测（单点求值数、约束扫描是否真出脊、能量噪声），所以每个里程碑开工前，再按同样格式补一份逐步计划给用户审。不提前写依赖未知数的代码。

| 里程碑 | 交付 | 依赖 | 出口（怎样算完成） | 资源 |
| --- | --- | --- | --- | --- |
| M0 环境核对 | 已完成（2026-10-13）：ASE 3.29.0，`FixInternals` 支持 bondcombos，`ase.mep.dimer` 可导入，geomeTRIC 1.1.1 | — | 结论记入本计划 | 本机只读 |
| **M1 计数账本** | `explore/ledger.py`；新包骨架 | — | Task 1–2 的测试全过；DoubleWell 上计数与外部计数器一致 | 本机 CPU |
| **M2 探针原语** | `explore/probe.py`；b04 真实模型冒烟 | M1 | Task 3–5 测试全过；冒烟给出单点 forward 数、墙钟、能量噪声实测；b04 上 Cl⁻→C 压缩扫描是否出脊如实报告 | 本机 CPU + 183 一张 GPU 冒烟 |
| M3 姿态与意图 | `poses.py`、`intents.py`；各体系意图数 × 姿态数的纯几何计数表 | M2 | 原子重标号不改变意图集合；计数表交用户，据此定 L1 规模 | 本机 CPU |
| M4 C2 认证 | `certify.py`（`follow_min_mode` → 内部 Hessian 一负模 → `descend_saddle` ± → `confirm_minimum` → 身份） | M2 | 解析双鞍点台两个鞍点都认证为 C2；非鞍点初猜得 unresolved/rejected | 本机 CPU |
| M5 D2 编排 | `d2.py`：固定分级、搜索额/认证保留额两账户、均匀抽查、kill→resume | M3、M4 | 解析台上 kill→resume 不重复计费；抽查项带 `audit_release` 标记 | 本机 CPU |
| M6 对照臂 | `arms.py`：D0（`prrs` 搜索 + 计数工厂）、D1（姿态 + pulse）、DB（ASE Dimer） | M1、M3、M4 | 四臂在解析台上同 L0、同账本跑通 | 本机 CPU |
| M7 参照认证 | b05、N1、N3 的目标 PES 参照注册表 + QM 核验 | M4 | 每项有驻点、一负模、两侧连接、QM 结构/连接判定；不合格体系标诊断 | 181/183 ORCA + 183 GPU，**先报计划** |
| M8 PARS 保真 | ALG-11/15/16（+02）修订，作为选择比较的独立臂 | — | 深挖报告的三个合成反例在修订版上不再成立；`make -C pars check` | 本机 CPU，可与 M2–M6 并行 |
| M9 选择比较 | 四臂：随机 K / 最低能量 K / 现行 PARS / 修订 PARS；冻结池后全池 oracle | M6、M8 | 解析台上跑通，口径同双鞍点台 | 本机 CPU |
| M10 零号批 | 每开发体系 6 配对块，其中 2 块同 seed 重跑 | M5–M7、M9 | 成本分布、噪声、可行性报告交用户；用户据此冻结 B*、δ、正式块数 | 183 GPU，**先报计划** |

v1 不做（v2 再议）：四账户、π_i 外推、认证反馈自适应、角/扭转/片段运动意图族、ART、committee/batching、mature A1 臂与 mature 预算预留修正（mature 进入比较时再做）。

D0 的阶段归属说明：D0 走 `prrs` 内部流程，`explore` 只能在 calculator 层计数，所以 D0 的账本阶段统一记 `d0_search`，不细分 quench/Hessian。这与协议 v1 的最小字段一致，在报告里注明。

---

### Task 1: 新包骨架与账本核心

**Files:**
- Create: `explore/pyproject.toml`、`explore/Makefile`、`explore/README.md`、`explore/LICENSE`（复制 `mature/LICENSE`）
- Create: `explore/src/explore/__init__.py`、`explore/src/explore/ledger.py`
- Test: `explore/tests/test_ledger.py`

**Interfaces:**
- Produces: `Ledger(path, arm, budget=None)`；`Ledger.segment(task, stage, attempt=0, reserve=None)` 上下文管理器，yield `Segment`（可写 `.outcome`、`.artifact`）；`Ledger.charge()`；`Ledger.note_request()`；`Ledger.remaining`；`Ledger.spent`；`Ledger.unreconciled`；`Ledger.close()`；异常 `BudgetExhausted`。

- [ ] **Step 1: 建包骨架**

`explore/pyproject.toml`：

```toml
[build-system]
requires = ["setuptools>=77"]
build-backend = "setuptools.build_meta"

[project]
name = "mlpath-explore"
version = "0.1.0"
description = "Discovery layer v1 (constrained-relaxation probe) and the evaluation ledger, on top of PRRS"
readme = "README.md"
requires-python = ">=3.10"
license = "AGPL-3.0-only"
license-files = ["LICENSE"]
authors = [{ name = "MLPATH / PRRS contributors" }]
dependencies = ["ase>=3.23,<4", "numpy>=1.24"]

[project.optional-dependencies]
test = ["pytest>=8"]

[tool.setuptools.packages.find]
where = ["src"]

[tool.pytest.ini_options]
testpaths = ["tests"]
pythonpath = ["src", "../src"]

[tool.ruff]
line-length = 96

[tool.ruff.lint]
select = ["E", "F", "W"]
ignore = ["E501"]
```

`explore/Makefile`：复制 `mature/Makefile`，把其中的 `src/mature` 全部替换为 `src/explore`。

`explore/README.md`：

```markdown
# explore — 发现层 v1 与评估账本

按 `docs/DISCOVERY_LAYER_DESIGN_2026-10-13.md` §1.1 与 `docs/EVALUATION_PROTOCOL_DESIGN_2026-10-13.md` v1 范围施工。
单向依赖 `prrs`（`../src`）；`prrs`、`pars` 不得 import 本包。测试：`make test`（需 `openmm_dev`）。
```

`explore/src/explore/__init__.py`：

```python
"""Discovery layer v1 and the evaluation ledger (see docs/*_2026-10-13.md)."""
```

- [ ] **Step 2: 写失败测试**

`explore/tests/test_ledger.py`：

```python
import json

import pytest

from explore.ledger import BudgetExhausted, Ledger


def lines(path):
    return [json.loads(x) for x in path.read_text(encoding="utf-8").splitlines() if x.strip()]


def test_segment_writes_open_and_close_with_counts(tmp_path):
    ledger = Ledger(tmp_path / "ledger.jsonl", arm="D2", budget=10)
    with ledger.segment(task="t1", stage="probe") as seg:
        ledger.note_request()
        ledger.charge()
        ledger.note_request()
        seg.outcome = "completed"
    events = lines(tmp_path / "ledger.jsonl")
    assert [e["event"] for e in events] == ["open", "close"]
    close = events[1]
    assert close["forwards"] == 1 and close["requests"] == 2
    assert close["outcome"] == "completed" and close["arm"] == "D2" and close["stage"] == "probe"
    assert ledger.spent == 1 and ledger.remaining == 9


def test_charge_refused_before_forward_when_budget_is_spent(tmp_path):
    ledger = Ledger(tmp_path / "l.jsonl", arm="D2", budget=1)
    with ledger.segment(task="t", stage="probe"):
        ledger.charge()
        with pytest.raises(BudgetExhausted):
            ledger.charge()
    assert ledger.spent == 1


def test_reservation_is_held_back_and_unused_part_returned(tmp_path):
    ledger = Ledger(tmp_path / "l.jsonl", arm="D2", budget=10)
    with ledger.segment(task="a", stage="probe", reserve=6):
        assert ledger.remaining == 4
        ledger.charge()
        ledger.charge()
        assert ledger.remaining == 4  # drawn from the reservation, not the pool
    assert ledger.spent == 2 and ledger.remaining == 8


def test_reservation_larger_than_remaining_is_refused_with_no_spend(tmp_path):
    ledger = Ledger(tmp_path / "l.jsonl", arm="D2", budget=5)
    with pytest.raises(BudgetExhausted):
        with ledger.segment(task="a", stage="probe", reserve=6):
            pass
    assert ledger.spent == 0
    assert lines(tmp_path / "l.jsonl") == []


def test_exception_inside_segment_is_recorded_and_reraised(tmp_path):
    ledger = Ledger(tmp_path / "l.jsonl", arm="D2", budget=5)
    with pytest.raises(ZeroDivisionError):
        with ledger.segment(task="a", stage="probe"):
            ledger.charge()
            1 / 0
    close = lines(tmp_path / "l.jsonl")[-1]
    assert close["outcome"] == "error:ZeroDivisionError" and close["forwards"] == 1


def test_reopen_after_crash_charges_open_segment_reservation(tmp_path):
    path = tmp_path / "l.jsonl"
    ledger = Ledger(path, arm="D2", budget=100)
    with ledger.segment(task="done", stage="probe"):
        ledger.charge()
    # simulate a crash: an open line with no close
    with path.open("a", encoding="utf-8") as f:
        f.write(json.dumps({"event": "open", "arm": "D2", "task": "crashed", "attempt": 0,
                            "stage": "probe", "reserve": 7}) + "\n")
    again = Ledger(path, arm="D2", budget=100)
    assert again.spent == 1 + 7
    assert again.unreconciled == [{"task": "crashed", "attempt": 0, "stage": "probe", "reserve": 7}]


def test_forwards_outside_segments_are_kept_as_unattributed(tmp_path):
    path = tmp_path / "l.jsonl"
    ledger = Ledger(path, arm="D0", budget=None)
    ledger.charge()
    ledger.close()
    close = lines(path)[-1]
    assert close["stage"] == "unattributed" and close["forwards"] == 1


def test_segments_do_not_nest(tmp_path):
    ledger = Ledger(tmp_path / "l.jsonl", arm="D2")
    with ledger.segment(task="a", stage="probe"):
        with pytest.raises(RuntimeError):
            with ledger.segment(task="b", stage="quench"):
                pass
```

- [ ] **Step 3: 跑测试确认失败**

Run: `cd /home/ruigengji/MLPATH/explore && /home/ruigengji/miniforge3/envs/openmm_dev/bin/python -m pytest tests/test_ledger.py -q`
Expected: FAIL，`ModuleNotFoundError: No module named 'explore.ledger'`

- [ ] **Step 4: 实现 `ledger.py` 核心**

```python
"""Forward ledger for the evaluation protocol (EVALUATION_PROTOCOL_DESIGN_2026-10-13 §2).

One unit is one real model forward on one structure. ASE's property cache sits in
front of ``calculate``: a cache hit is a request, not a forward. A forward that was
started and then failed is still paid for. Reopening a ledger never refunds anything:
a segment that was opened and never closed is charged its reservation as an upper
bound and listed in ``unreconciled``.
"""

from __future__ import annotations

import contextlib
import json
import os
import socket
import time
from pathlib import Path


class BudgetExhausted(RuntimeError):
    """Raised before a forward (or a reservation) the remaining budget cannot cover."""


class Segment:
    def __init__(self, task, attempt, stage, reserve):
        self.task, self.attempt, self.stage = task, attempt, stage
        self.reserve = reserve
        self.reserve_left = reserve
        self.forwards = 0
        self.requests = 0
        self.outcome = None
        self.artifact = None


class Ledger:
    def __init__(self, path, arm, budget=None):
        self.path = Path(path)
        self.arm = arm
        self.budget = budget
        self.spent = 0
        self.reserved = 0
        self.unreconciled = []
        self.host = socket.gethostname()
        self._segment = None
        self._unattributed = 0
        if self.path.exists():
            self._replay()

    @property
    def remaining(self):
        if self.budget is None:
            return None
        return self.budget - self.spent - self.reserved

    def _replay(self):
        pending = {}
        for line in self.path.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            event = json.loads(line)
            key = (event["task"], event["attempt"], event["stage"])
            if event["event"] == "open":
                pending[key] = event
            elif event["event"] == "close":
                pending.pop(key, None)
                self.spent += event["forwards"]
        for (task, attempt, stage), event in pending.items():
            self.spent += event["reserve"] or 0
            self.unreconciled.append(
                {"task": task, "attempt": attempt, "stage": stage, "reserve": event["reserve"]}
            )

    def _write(self, event):
        record = {"arm": self.arm, "host": self.host, "t": time.time(), **event}
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with self.path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(record) + "\n")
            handle.flush()
            os.fsync(handle.fileno())

    @contextlib.contextmanager
    def segment(self, task, stage, attempt=0, reserve=None):
        if self._segment is not None:
            raise RuntimeError("ledger segments do not nest")
        if reserve is not None and self.budget is not None and reserve > self.remaining:
            raise BudgetExhausted(f"reserve {reserve} exceeds remaining {self.remaining}")
        seg = Segment(task, attempt, stage, reserve)
        ident = {"task": task, "attempt": attempt, "stage": stage}
        self._write({"event": "open", **ident, "reserve": reserve})
        if reserve is not None:
            self.reserved += reserve
        self._segment = seg
        started = time.monotonic()
        try:
            yield seg
        except BaseException as exc:
            if seg.outcome is None:
                seg.outcome = f"error:{type(exc).__name__}"
            raise
        finally:
            self._segment = None
            if seg.reserve_left is not None:
                self.reserved -= seg.reserve_left
            self._write(
                {
                    "event": "close",
                    **ident,
                    "forwards": seg.forwards,
                    "requests": seg.requests,
                    "wall_s": time.monotonic() - started,
                    "outcome": seg.outcome,
                    "artifact": seg.artifact,
                }
            )

    def charge(self):
        seg = self._segment
        if seg is not None and seg.reserve_left:
            seg.reserve_left -= 1
            self.reserved -= 1
        elif self.budget is not None and self.remaining < 1:
            raise BudgetExhausted(f"no budget left for a forward (spent {self.spent})")
        self.spent += 1
        if seg is None:
            self._unattributed += 1
        else:
            seg.forwards += 1

    def note_request(self):
        if self._segment is not None:
            self._segment.requests += 1

    def close(self):
        if self._unattributed:
            ident = {"task": "<unattributed>", "attempt": 0, "stage": "unattributed"}
            self._write(
                {"event": "close", **ident, "forwards": self._unattributed, "requests": 0,
                 "wall_s": None, "outcome": None, "artifact": None}
            )
            self._unattributed = 0
```

- [ ] **Step 5: 跑测试确认通过**

Run: `cd /home/ruigengji/MLPATH/explore && /home/ruigengji/miniforge3/envs/openmm_dev/bin/python -m pytest tests/test_ledger.py -q`
Expected: 8 passed

- [ ] **Step 6: 格式与 lint**

Run: `cd /home/ruigengji/MLPATH/explore && make fmt && make lint`
Expected: lint 无报错

### Task 2: calculator 计数包装

**Files:**
- Modify: `explore/src/explore/ledger.py`（文件末尾追加）
- Test: `explore/tests/test_ledger.py`（追加）

**Interfaces:**
- Consumes: Task 1 的 `Ledger.charge()`、`Ledger.note_request()`
- Produces: `instrument(calc, ledger) -> calc`；`counted_factory(factory, ledger) -> callable`。被包装的 calculator 的 `calculate`、`get_property`、`get_hessian`（若有）都经过账本；`get_hessian` 计 1 forward。

- [ ] **Step 1: 写失败测试**

追加到 `explore/tests/test_ledger.py`（新增的 import 行并入文件顶部的 import 区：ruff 规则 E402 不允许文件中部 import）：

```python
import numpy as np
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes

from explore.ledger import counted_factory, instrument
from prrs.calculators import double_well_factory


def dimer(r):
    return Atoms("C2", positions=[[0, 0, 0], [r, 0, 0]])


def test_cache_hit_is_a_request_not_a_forward(tmp_path):
    ledger = Ledger(tmp_path / "l.jsonl", arm="D2", budget=10)
    atoms = dimer(1.3)
    atoms.calc = instrument(double_well_factory(), ledger)
    with ledger.segment(task="t", stage="probe"):
        atoms.get_potential_energy()
        atoms.get_forces()  # same geometry: ASE cache
        atoms.positions[1, 0] = 1.4
        atoms.get_forces()  # new geometry: a real forward
    close = lines(tmp_path / "l.jsonl")[-1]
    assert close["forwards"] == 2 and close["requests"] == 3


class Boom(Calculator):
    implemented_properties = ["energy", "forces"]

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        raise RuntimeError("model blew up")


def test_failed_forward_is_still_charged(tmp_path):
    ledger = Ledger(tmp_path / "l.jsonl", arm="D2", budget=10)
    atoms = dimer(1.3)
    atoms.calc = instrument(Boom(), ledger)
    with pytest.raises(RuntimeError):
        with ledger.segment(task="t", stage="probe"):
            atoms.get_potential_energy()
    assert ledger.spent == 1


def test_budget_refusal_happens_before_the_model_is_called(tmp_path):
    calls = []

    class Spy(Calculator):
        implemented_properties = ["energy", "forces"]

        def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
            calls.append(1)
            self.results = {"energy": 0.0, "forces": np.zeros((len(atoms), 3))}

    ledger = Ledger(tmp_path / "l.jsonl", arm="D2", budget=0)
    atoms = dimer(1.3)
    atoms.calc = instrument(Spy(), ledger)
    with ledger.segment(task="t", stage="probe"):
        with pytest.raises(BudgetExhausted):
            atoms.get_potential_energy()
    assert calls == []


def test_double_instrumentation_is_an_error(tmp_path):
    ledger = Ledger(tmp_path / "l.jsonl", arm="D2")
    calc = instrument(double_well_factory(), ledger)
    with pytest.raises(ValueError):
        instrument(calc, ledger)


def test_counted_factory_counts_through_prrs_guard(tmp_path):
    from prrs.config import SearchConfig
    from prrs.reliability import GuardedCalculator

    ledger = Ledger(tmp_path / "l.jsonl", arm="D2", budget=10)
    factory = counted_factory(double_well_factory, ledger)
    atoms = dimer(1.3)
    atoms.calc = GuardedCalculator(factory(), SearchConfig(closed_shell_only=False))
    with ledger.segment(task="t", stage="probe"):
        atoms.get_potential_energy()
        atoms.get_forces()
    assert ledger.spent == 1


def test_analytic_hessian_is_charged(tmp_path):
    class WithHessian(Calculator):
        implemented_properties = ["energy", "forces"]

        def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
            self.results = {"energy": 0.0, "forces": np.zeros((len(atoms), 3))}

        def get_hessian(self, atoms=None):
            return np.zeros((3 * len(atoms), 3 * len(atoms)))

    ledger = Ledger(tmp_path / "l.jsonl", arm="D2", budget=10)
    calc = instrument(WithHessian(), ledger)
    with ledger.segment(task="t", stage="hessian"):
        calc.get_hessian(dimer(1.3))
    assert ledger.spent == 1
```

- [ ] **Step 2: 跑测试确认失败**

Run: `cd /home/ruigengji/MLPATH/explore && /home/ruigengji/miniforge3/envs/openmm_dev/bin/python -m pytest tests/test_ledger.py -q`
Expected: FAIL，`ImportError: cannot import name 'counted_factory'`

- [ ] **Step 3: 实现包装**

追加到 `explore/src/explore/ledger.py`：

```python
def instrument(calc, ledger):
    """Route a physical calculator's real evaluations through the ledger.

    ``calculate`` is where ASE does the work and is charged *before* the call, so a
    refused forward never reaches the model and a failed one is still paid for.
    ``get_property`` is counted as a request (cache hits included). A backend's own
    ``get_hessian`` bypasses ``calculate`` and is charged as one forward.
    """
    if getattr(calc, "_explore_ledger", None) is not None:
        raise ValueError("calculator is already counted by a ledger")
    inner_calculate = calc.calculate
    inner_get_property = calc.get_property

    def calculate(*args, **kwargs):
        ledger.charge()
        return inner_calculate(*args, **kwargs)

    def get_property(name, atoms=None, allow_calculation=True):
        ledger.note_request()
        return inner_get_property(name, atoms, allow_calculation)

    calc.calculate = calculate
    calc.get_property = get_property
    if hasattr(calc, "get_hessian"):
        inner_hessian = calc.get_hessian

        def get_hessian(*args, **kwargs):
            ledger.charge()
            return inner_hessian(*args, **kwargs)

        calc.get_hessian = get_hessian
    calc._explore_ledger = ledger
    return calc


def counted_factory(factory, ledger):
    """A factory whose every calculator is counted by ``ledger``."""

    def make():
        return instrument(factory(), ledger)

    return make
```

- [ ] **Step 4: 跑测试确认通过**

Run: `cd /home/ruigengji/MLPATH/explore && /home/ruigengji/miniforge3/envs/openmm_dev/bin/python -m pytest tests/test_ledger.py -q`
Expected: 14 passed

- [ ] **Step 5: 格式与 lint**

Run: `cd /home/ruigengji/MLPATH/explore && make fmt && make lint`
Expected: lint 无报错

- [ ] **Step 6: M1 commit 草稿**

写 `.commit_drafts/<YYYY-MM-DD_HHMM>_explore-ledger.md`，内容三项：commit message `explore: forward ledger and counted calculator factory (protocol §2 v1)`；文件清单（`explore/` 下新建的 8 个文件，每个一句话）；可直接复制的 `git add explore/ && git commit -m "…"`。

### Task 3: Intent 与约束松弛扫描

**Files:**
- Create: `explore/src/explore/probe.py`
- Test: `explore/tests/test_probe.py`

**Interfaces:**
- Consumes: `Ledger`、`Ledger.segment`、`BudgetExhausted`、`counted_factory`（Task 1–2）；`prrs.reliability.GuardedCalculator`、`GateRejected`；`prrs.state.encode`。
- Produces:
  - `Intent(intent_id: str, pairs: tuple[tuple[int, int], ...], coefs: tuple[float, ...])`，方法 `value(positions) -> float`、`gradient(positions) -> np.ndarray`。q(x) = Σ c_k r_k，单位 Å。
  - `ScanSettings(deltas_A: tuple[float, ...], fmax_eV_A: float, max_steps_per_point: int, max_rise_eV: float)`；`deltas_A[0]` 必须是 0.0。
  - `ScanPoint`（字段：`delta_A, q_target, q_actual, energy_eV, reaction_force, fmax_projected, converged, steps, edges, positions`）。
  - `ScanResult(intent_id: str, q0: float, points: list[ScanPoint], status: str)`；status ∈ `completed / censored_energy / budget_exhausted / gate_rejected:<code> / infeasible`。
  - `constrained_scan(atoms, intent, settings, factory, config, ledger, task) -> ScanResult`；`factory` 应是 `counted_factory` 的产物。

- [ ] **Step 1: 写失败测试**

`explore/tests/test_probe.py`：

```python
import numpy as np
import pytest
from ase import Atoms
from ase.calculators.lj import LennardJones

from explore.ledger import Ledger, counted_factory
from explore.probe import Intent, ScanSettings, constrained_scan
from prrs.calculators import double_well_factory
from prrs.config import SearchConfig

CONFIG = SearchConfig(closed_shell_only=False, quench_fmax_eV_A=0.001)


def dw_energy(r):
    x = (r - 1.8) / 0.6
    return 0.4 * (x * x - 1) ** 2


def dw_slope(r):
    x = (r - 1.8) / 0.6
    return 1.6 * x * (x * x - 1) / 0.6


def test_intent_gradient_matches_finite_difference():
    rng = np.random.default_rng(0)
    pos = rng.normal(size=(4, 3)) * 1.5
    intent = Intent("i", ((0, 1), (2, 3)), (1.0, -0.5))
    grad = intent.gradient(pos)
    h = 1e-6
    for a in range(4):
        for k in range(3):
            plus, minus = pos.copy(), pos.copy()
            plus[a, k] += h
            minus[a, k] -= h
            fd = (intent.value(plus) - intent.value(minus)) / (2 * h)
            assert abs(fd - grad[a, k]) < 1e-6


def scan_double_well(tmp_path, budget=1000, deltas=None):
    ledger = Ledger(tmp_path / "l.jsonl", arm="D2", budget=budget)
    atoms = Atoms("C2", positions=[[0, 0, 0], [1.2, 0, 0]])
    settings = ScanSettings(
        deltas_A=deltas or tuple(round(0.1 * k, 10) for k in range(11)),
        fmax_eV_A=1e-4,
        max_steps_per_point=50,
        max_rise_eV=10.0,
    )
    result = constrained_scan(
        atoms, Intent("cc", ((0, 1),), (1.0,)), settings,
        counted_factory(double_well_factory, ledger), CONFIG, ledger, task="dw",
    )
    return result, ledger


def test_double_well_scan_reproduces_the_analytic_profile(tmp_path):
    result, ledger = scan_double_well(tmp_path)
    assert result.status == "completed" and len(result.points) == 11
    for p in result.points:
        r = 1.2 + p.delta_A
        assert abs(p.q_actual - r) < 1e-6
        assert abs((p.energy_eV - result.points[0].energy_eV) - (dw_energy(r) - dw_energy(1.2))) < 1e-8
        assert abs(p.reaction_force - dw_slope(r)) < 1e-6
    assert ledger.spent > 0


def test_unconstrained_degrees_of_freedom_relax(tmp_path):
    ledger = Ledger(tmp_path / "l.jsonl", arm="D2", budget=5000)
    atoms = Atoms("Ar3", positions=[[0, 0, 0], [1.12, 0, 0], [0.3, 1.6, 0]])
    start = atoms.positions[2].copy()
    settings = ScanSettings(deltas_A=(0.0, 0.1), fmax_eV_A=1e-3, max_steps_per_point=500, max_rise_eV=10.0)
    result = constrained_scan(
        atoms, Intent("ar", ((0, 1),), (1.0,)), settings,
        counted_factory(lambda: LennardJones(sigma=1.0, epsilon=0.01, rc=10.0), ledger),
        CONFIG, ledger, task="lj",
    )
    last = result.points[-1]
    assert last.converged and last.fmax_projected < 1e-3
    assert np.linalg.norm(last.positions[2] - start) > 0.01


def test_energy_cap_censors_and_keeps_points(tmp_path):
    ledger = Ledger(tmp_path / "l.jsonl", arm="D2", budget=1000)
    atoms = Atoms("C2", positions=[[0, 0, 0], [1.2, 0, 0]])
    settings = ScanSettings(tuple(round(0.1 * k, 10) for k in range(11)), 1e-4, 50, max_rise_eV=0.2)
    result = constrained_scan(
        atoms, Intent("cc", ((0, 1),), (1.0,)), settings,
        counted_factory(double_well_factory, ledger), CONFIG, ledger, task="cap",
    )
    assert result.status == "censored_energy"
    assert result.points[-1].energy_eV - result.points[0].energy_eV > 0.2
    assert all(p.energy_eV - result.points[0].energy_eV <= 0.2 for p in result.points[:-1])


def test_budget_exhaustion_mid_scan_keeps_points(tmp_path):
    result, ledger = scan_double_well(tmp_path, budget=3)
    assert result.status == "budget_exhausted"
    assert ledger.spent == 3


def test_overlap_gate_ends_scan_without_raising(tmp_path):
    deltas = tuple(round(-0.1 * k, 10) for k in range(12))  # compress 1.2 A towards 0
    result, _ = scan_double_well(tmp_path, deltas=deltas)
    assert result.status == "gate_rejected:overlap"
    assert len(result.points) >= 1


def test_negative_target_distance_is_infeasible_before_any_forward(tmp_path):
    result, ledger = scan_double_well(tmp_path, deltas=(0.0, -1.5))
    assert result.status == "infeasible" and result.points == [] and ledger.spent == 0


def test_first_delta_must_be_zero():
    with pytest.raises(ValueError):
        ScanSettings(deltas_A=(0.1, 0.2), fmax_eV_A=1e-3, max_steps_per_point=10, max_rise_eV=1.0)
```

- [ ] **Step 2: 跑测试确认失败**

Run: `cd /home/ruigengji/MLPATH/explore && /home/ruigengji/miniforge3/envs/openmm_dev/bin/python -m pytest tests/test_probe.py -q`
Expected: FAIL，`ModuleNotFoundError: No module named 'explore.probe'`

- [ ] **Step 3: 实现**

`explore/src/explore/probe.py`：

```python
"""Bounded constrained-relaxation probe (DISCOVERY_LAYER_DESIGN_2026-10-13 §4, v1).

At each target value of the intent coordinate q(x) = sum_k c_k r_k the structure is
relaxed on the unbiased PES with q held exactly (ASE FixInternals bondcombos); every
other internal degree of freedom is free. The energy recorded is the PES energy --
there is no penalty term -- and ``reaction_force`` is dE*/dq read from the
unconstrained forces, which at a constrained minimum is the Lagrange multiplier.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from ase.constraints import FixInternals
from ase.optimize import BFGS

from .ledger import BudgetExhausted


@dataclass(frozen=True)
class Intent:
    intent_id: str
    pairs: tuple[tuple[int, int], ...]
    coefs: tuple[float, ...]

    def value(self, positions) -> float:
        return float(
            sum(c * np.linalg.norm(positions[j] - positions[i]) for (i, j), c in zip(self.pairs, self.coefs))
        )

    def gradient(self, positions) -> np.ndarray:
        grad = np.zeros_like(positions, dtype=float)
        for (i, j), c in zip(self.pairs, self.coefs):
            v = positions[j] - positions[i]
            u = v / np.linalg.norm(v)
            grad[j] += c * u
            grad[i] -= c * u
        return grad

    def combo(self):
        return [[i, j, c] for (i, j), c in zip(self.pairs, self.coefs)]


@dataclass(frozen=True)
class ScanSettings:
    deltas_A: tuple[float, ...]
    fmax_eV_A: float
    max_steps_per_point: int
    max_rise_eV: float

    def __post_init__(self):
        if not self.deltas_A or self.deltas_A[0] != 0.0:
            raise ValueError("deltas_A must start at 0.0 (the intent's own start)")


@dataclass
class ScanPoint:
    delta_A: float
    q_target: float
    q_actual: float
    energy_eV: float
    reaction_force: float
    fmax_projected: float
    converged: bool
    steps: int
    edges: list
    positions: np.ndarray


@dataclass
class ScanResult:
    intent_id: str
    q0: float
    points: list
    status: str


def _feasible(intent, q0, deltas):
    # single-pair intents: a target distance must stay positive
    if len(intent.pairs) == 1:
        c = intent.coefs[0]
        return all((q0 + d) / c > 0 for d in deltas)
    return True


def constrained_scan(atoms, intent, settings, factory, config, ledger, task):
    from prrs.reliability import GateRejected, GuardedCalculator
    from prrs.state import encode

    work = atoms.copy()
    work.set_constraint()
    q0 = intent.value(work.positions)
    points = []
    if not _feasible(intent, q0, settings.deltas_A):
        return ScanResult(intent.intent_id, q0, points, "infeasible")
    work.calc = GuardedCalculator(factory(), config)
    status = "completed"
    e_start = None
    try:
        with ledger.segment(task=task, stage="probe") as seg:
            for delta in settings.deltas_A:
                target = q0 + delta
                work.set_constraint(FixInternals(bondcombos=[[target, intent.combo()]]))
                opt = BFGS(work, logfile=None)
                converged = bool(opt.run(fmax=settings.fmax_eV_A, steps=settings.max_steps_per_point))
                energy = float(work.get_potential_energy())
                f_pes = work.get_forces(apply_constraint=False)
                g = intent.gradient(work.positions)
                point = ScanPoint(
                    delta_A=float(delta),
                    q_target=float(target),
                    q_actual=intent.value(work.positions),
                    energy_eV=energy,
                    reaction_force=-float(np.vdot(f_pes, g) / np.vdot(g, g)),
                    fmax_projected=float(np.max(np.linalg.norm(work.get_forces(), axis=1))),
                    converged=converged,
                    steps=int(opt.nsteps),
                    edges=sorted(map(list, encode(work, config.bond_scale).edges)),
                    positions=work.positions.copy(),
                )
                points.append(point)
                e_start = energy if e_start is None else e_start
                if energy - e_start > settings.max_rise_eV:
                    status = "censored_energy"
                    break
            seg.outcome = status
    except BudgetExhausted:
        status = "budget_exhausted"
    except GateRejected as exc:
        status = f"gate_rejected:{exc.code}"
    finally:
        work.set_constraint()
    return ScanResult(intent.intent_id, q0, points, status)
```

注：门代码取 `GateRejected.code`（`src/prrs/reliability.py:8-12`；`args[0]` 是 message，不是代码）。

- [ ] **Step 4: 跑测试确认通过**

Run: `cd /home/ruigengji/MLPATH/explore && /home/ruigengji/miniforge3/envs/openmm_dev/bin/python -m pytest tests/test_probe.py -q`
Expected: 8 passed

- [ ] **Step 5: 格式与 lint**

Run: `cd /home/ruigengji/MLPATH/explore && make fmt && make lint`

### Task 4: 扫描判读——`ridge_suspected` 与 `branch_bracketed`

**Files:**
- Modify: `explore/src/explore/probe.py`（追加）
- Test: `explore/tests/test_probe.py`（追加）

**Interfaces:**
- Consumes: `ScanPoint`、`ScanResult`、`constrained_scan`（Task 3）
- Produces:
  - `read_scan(points, noise_eV) -> list[dict]`：纯函数。返回触发列表，每项 `{"kind": "ridge_suspected", "peak": k, "release": k+1}` 或 `{"kind": "jump", "before": k-1, "after": k}`。
  - `bracket_jump(points, jump, reverse_point) -> dict | None`：纯函数。反向扫描在 `before` 处仍停留在跳后的键图上，就返回 `{"kind": "branch_bracketed", "reactant_side": before, "product_side_positions": reverse_point.positions}`，否则返回 None。
  - `reverse_check(points, jump, intent, settings, factory, config, ledger, task) -> dict | None`：从 `points[after]` 的几何出发，用 `constrained_scan` 只扫一点（`deltas_A=(0.0, points[before].delta_A - points[after].delta_A)`），再交给 `bracket_jump`。

判读规则（发现层 §4.3，v1）：
- `ridge_suspected`：存在内部点 k，前 k+2 个点都收敛；E[k] − E[0] > noise；E[k] − min(E[k+1:]) > noise；`reaction_force[k-1]` 与 `reaction_force[k+1]` 异号。取第一个满足的 k。
- `jump`：相邻收敛点 k−1、k 的键图不同，且 E[k] < E[k−1] − noise（跨过键图变化时能量下降）。
- `noise_eV` 由调用方传入；它的实测来源见 Task 6，测试里显式给值。

- [ ] **Step 1: 写失败测试**

追加到 `explore/tests/test_probe.py`（新增的 import 行并入文件顶部的 import 区：ruff 规则 E402 不允许文件中部 import）：

```python
from explore.probe import ScanPoint, bracket_jump, read_scan


def fake_point(delta, energy, rf, edges, converged=True):
    return ScanPoint(delta, 0.0, 0.0, energy, rf, 0.0, converged, 0, edges, np.zeros((2, 3)))


def test_double_well_scan_reads_one_ridge_at_the_barrier(tmp_path):
    result, _ = scan_double_well(tmp_path)
    triggers = read_scan(result.points, noise_eV=1e-3)
    ridges = [t for t in triggers if t["kind"] == "ridge_suspected"]
    assert len(ridges) == 1
    assert abs(result.points[ridges[0]["peak"]].q_actual - 1.8) < 1e-6
    assert ridges[0]["release"] == ridges[0]["peak"] + 1


def test_monotone_rise_reads_nothing():
    pts = [fake_point(0.1 * k, 0.1 * k, 1.0, [[0, 1]]) for k in range(6)]
    assert read_scan(pts, noise_eV=1e-3) == []


def test_peak_below_noise_is_not_a_ridge():
    es = [0.0, 0.0005, 0.0008, 0.0005, 0.0]
    rfs = [1.0, 1.0, 0.0, -1.0, -1.0]
    pts = [fake_point(0.1 * k, e, f, [[0, 1]]) for k, (e, f) in enumerate(zip(es, rfs))]
    assert read_scan(pts, noise_eV=1e-3) == []


def test_unconverged_point_before_the_peak_blocks_the_ridge():
    es, rfs = [0.0, 0.2, 0.4, 0.2, 0.0], [1.0, 1.0, 0.0, -1.0, -1.0]
    pts = [fake_point(0.1 * k, e, f, [[0, 1]], converged=(k != 1)) for k, (e, f) in enumerate(zip(es, rfs))]
    assert [t for t in read_scan(pts, 1e-3) if t["kind"] == "ridge_suspected"] == []


def test_graph_change_with_energy_drop_is_a_jump():
    pts = [
        fake_point(0.0, 0.0, 1.0, [[0, 1]]),
        fake_point(0.1, 0.5, 1.0, [[0, 1]]),
        fake_point(0.2, -0.3, -1.0, [[0, 1], [1, 2]]),
    ]
    assert {"kind": "jump", "before": 1, "after": 2} in read_scan(pts, noise_eV=1e-3)


def test_reverse_point_staying_on_the_new_graph_brackets_the_jump():
    pts = [fake_point(0.0, 0.0, 1.0, [[0, 1]]), fake_point(0.1, -0.3, -1.0, [[0, 1], [1, 2]])]
    jump = {"kind": "jump", "before": 0, "after": 1}
    stayed = fake_point(0.0, -0.1, -1.0, [[0, 1], [1, 2]])
    went_back = fake_point(0.0, 0.0, 1.0, [[0, 1]])
    assert bracket_jump(pts, jump, stayed)["kind"] == "branch_bracketed"
    assert bracket_jump(pts, jump, went_back) is None
```

- [ ] **Step 2: 跑测试确认失败**

Run: `cd /home/ruigengji/MLPATH/explore && /home/ruigengji/miniforge3/envs/openmm_dev/bin/python -m pytest tests/test_probe.py -q`
Expected: FAIL，`ImportError: cannot import name 'read_scan'`

- [ ] **Step 3: 实现**

追加到 `explore/src/explore/probe.py`：

```python
def read_scan(points, noise_eV):
    """Release triggers from a finished scan (§4.3). Pure: no evaluations."""
    triggers = []
    energies = [p.energy_eV for p in points]
    for k in range(1, len(points) - 1):
        if not all(p.converged for p in points[: k + 2]):
            break
        rises = energies[k] - energies[0] > noise_eV
        falls = energies[k] - min(energies[k + 1 :]) > noise_eV
        flips = points[k - 1].reaction_force * points[k + 1].reaction_force < 0
        if rises and falls and flips and energies[k] >= energies[k - 1] and energies[k] >= energies[k + 1]:
            triggers.append({"kind": "ridge_suspected", "peak": k, "release": k + 1})
            break
    for k in range(1, len(points)):
        a, b = points[k - 1], points[k]
        if a.converged and b.converged and a.edges != b.edges and b.energy_eV < a.energy_eV - noise_eV:
            triggers.append({"kind": "jump", "before": k - 1, "after": k})
    return triggers


def bracket_jump(points, jump, reverse_point):
    """A jump is bracketed when, scanned back to ``before``, the structure stays on the
    graph it jumped to: two branches coexist at that target (hysteresis)."""
    before, after = points[jump["before"]], points[jump["after"]]
    if reverse_point.converged and reverse_point.edges == after.edges != before.edges:
        return {
            "kind": "branch_bracketed",
            "reactant_side": jump["before"],
            "product_side_positions": reverse_point.positions,
        }
    return None


def reverse_check(atoms_template, points, jump, intent, settings, factory, config, ledger, task):
    before, after = points[jump["before"]], points[jump["after"]]
    start = atoms_template.copy()
    start.positions = after.positions.copy()
    back = ScanSettings(
        deltas_A=(0.0, before.delta_A - after.delta_A),
        fmax_eV_A=settings.fmax_eV_A,
        max_steps_per_point=settings.max_steps_per_point,
        max_rise_eV=settings.max_rise_eV,
    )
    result = constrained_scan(start, intent, back, factory, config, ledger, task)
    if result.status != "completed" or len(result.points) < 2:
        return None
    return bracket_jump(points, jump, result.points[-1])
```

`reverse_check` 的第一个参数是任一带正确元素与 info（电荷等）的 `Atoms` 模板，只取其元素与 info，坐标来自扫描点。上面 Interfaces 段里 `reverse_check` 的签名以这里为准（多一个 `atoms_template`）。

- [ ] **Step 4: 跑测试确认通过**

Run: `cd /home/ruigengji/MLPATH/explore && /home/ruigengji/miniforge3/envs/openmm_dev/bin/python -m pytest tests/test_probe.py -q`
Expected: 14 passed

- [ ] **Step 5: 格式与 lint**

Run: `cd /home/ruigengji/MLPATH/explore && make fmt && make lint`

### Task 5: 释放——解约束、淬火、极小确认、分类

**Files:**
- Modify: `explore/src/explore/probe.py`（追加）
- Test: `explore/tests/test_probe.py`（追加）

**Interfaces:**
- Consumes: `prrs.runner.quench(atoms, guard, config)`（返回 `(converged, diagnostics)`）、`prrs.runner.confirm_minimum(candidate, factory, config, seed)`（返回 `(confirmed, diagnostics)`）、`prrs.io.snapshot(atoms, physical)`、`prrs.state.classify(source, endpoint, config)`；`Ledger.segment`。
- Produces: `release(positions, template, source, factory, config, ledger, task, seed) -> dict`，键：`status`（`c1_candidate / quench_failed / not_minimum / budget_exhausted / gate_rejected:<code>`）、`label`（classify 的 label 或 None）、`fingerprint`（classify 的 fingerprint 或 None）、`endpoint`（`Atoms` 或 None）。`source` 是带能量的来源快照（`snapshot` 产物）。

- [ ] **Step 1: 写失败测试**

追加到 `explore/tests/test_probe.py`（新增的 import 行并入文件顶部的 import 区：ruff 规则 E402 不允许文件中部 import）：

```python
from explore.probe import release
from prrs.io import snapshot
from prrs.reliability import GuardedCalculator


def dw_source():
    atoms = Atoms("C2", positions=[[0, 0, 0], [1.2, 0, 0]])
    return snapshot(atoms, GuardedCalculator(double_well_factory(), CONFIG))


def test_release_past_the_ridge_lands_in_the_far_basin(tmp_path):
    ledger = Ledger(tmp_path / "l.jsonl", arm="D2", budget=5000)
    template = Atoms("C2", positions=np.zeros((2, 3)))
    out = release(np.array([[0, 0, 0], [2.0, 0, 0]], float), template, dw_source(),
                  counted_factory(double_well_factory, ledger), CONFIG, ledger, task="rel", seed=0)
    assert out["status"] == "c1_candidate"
    assert abs(out["endpoint"].get_distance(0, 1) - 2.4) < 1e-2
    assert out["label"] == "dissociation"
    stages = [json.loads(x)["stage"] for x in (tmp_path / "l.jsonl").read_text().splitlines()]
    assert "quench" in stages and "minimum_check" in stages


def test_release_before_the_ridge_returns_to_source(tmp_path):
    ledger = Ledger(tmp_path / "l.jsonl", arm="D2", budget=5000)
    template = Atoms("C2", positions=np.zeros((2, 3)))
    out = release(np.array([[0, 0, 0], [1.6, 0, 0]], float), template, dw_source(),
                  counted_factory(double_well_factory, ledger), CONFIG, ledger, task="rel", seed=0)
    assert out["label"] == "elastic"
```

并在 `test_probe.py` 顶部 import 区加 `import json`。

- [ ] **Step 2: 跑测试确认失败**

Run: `cd /home/ruigengji/MLPATH/explore && /home/ruigengji/miniforge3/envs/openmm_dev/bin/python -m pytest tests/test_probe.py -q`
Expected: FAIL，`ImportError: cannot import name 'release'`

- [ ] **Step 3: 实现**

追加到 `explore/src/explore/probe.py`：

```python
def release(positions, template, source, factory, config, ledger, task, seed):
    """Drop every external constraint, quench on the unbiased PES, confirm the minimum,
    classify against the source. Static: no momenta, no claim of a thermal crossing."""
    from prrs.io import snapshot
    from prrs.reliability import GateRejected, GuardedCalculator
    from prrs.runner import confirm_minimum, quench
    from prrs.state import classify

    atoms = template.copy()
    atoms.set_constraint()
    atoms.positions = np.asarray(positions, float).copy()
    atoms.set_momenta(np.zeros((len(atoms), 3)))
    out = {"status": None, "label": None, "fingerprint": None, "endpoint": None}
    try:
        guard = GuardedCalculator(factory(), config)
        atoms.calc = guard
        with ledger.segment(task=task, stage="quench") as seg:
            converged, _ = quench(atoms, guard, config)
            seg.outcome = "converged" if converged else "quench_failed"
        if not converged:
            out["status"] = "quench_failed"
            return out
        endpoint = snapshot(atoms, guard)
        with ledger.segment(task=task, stage="minimum_check") as seg:
            confirmed, _ = confirm_minimum(endpoint, factory, config, seed)
            seg.outcome = "minimum" if confirmed else "not_minimum"
        if not confirmed:
            out["status"] = "not_minimum"
            return out
        label, fingerprint = classify(source, endpoint, config)
        out.update(status="c1_candidate", label=label, fingerprint=fingerprint, endpoint=endpoint)
    except BudgetExhausted:
        out["status"] = "budget_exhausted"
    except GateRejected as exc:
        out["status"] = f"gate_rejected:{exc.code}"
    return out
```

注：`prrs.state.classify` 返回 `(label, fingerprint)`（`src/prrs/state.py:142-175`，已核实）。`c1_candidate` 只表示“释放到了确认极小”；能否标 C1 还要看来源盆地是否确认（协议 §1.1），由 M5 的编排判断。

- [ ] **Step 4: 跑测试确认通过**

Run: `cd /home/ruigengji/MLPATH/explore && /home/ruigengji/miniforge3/envs/openmm_dev/bin/python -m pytest tests/test_probe.py -q`
Expected: 16 passed

- [ ] **Step 5: 全包检查**

Run: `cd /home/ruigengji/MLPATH/explore && make check`
Expected: fmt-check、lint、test 全过

### Task 6: b04 真实模型冒烟（测单点成本与能量噪声）

**Files:**
- Create: `benchmarks/explore_smoke/b04_scan.py`
- Create（运行产物）: `benchmarks/explore_smoke/results/<YYYY-MM-DD>/`（`ledger.jsonl`、`scan.json`、`noise.json`、`run.log`、`pids.txt`）

**Interfaces:**
- Consumes: `Intent`、`ScanSettings`、`constrained_scan`、`read_scan`、`release`、`Ledger`、`counted_factory`；`benchmarks.b04_sn2_chloride.calculator:polar1`；`prrs.search._stamped_with_charge`。

这是冒烟：一个意图、十几个点、几分钟就停，不需要批准。目的只有三个：单点 forward 数与墙钟、能量噪声、这个原语在真实 PES 上能不能出脊。b04 的已知 TS 只用于事后对照，不进入扫描。

- [ ] **Step 1: 写脚本**

`benchmarks/explore_smoke/b04_scan.py`：

```python
"""Smoke: one constrained-relaxation scan on b04 (Cl- + CH3Cl, POLAR-1, charge -1).

Measures per-point forwards and wall time, the energy noise of a relaxed point, and
whether the probe reads a ridge. Not an acceptance run; the known TS is only compared
after the fact.
"""

import json
import sys
from dataclasses import asdict
from pathlib import Path

import numpy as np
from ase.io import read

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "explore" / "src"), str(ROOT)]

from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.ledger import Ledger, counted_factory  # noqa: E402
from explore.probe import Intent, ScanSettings, constrained_scan, read_scan, release  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402
from prrs.io import snapshot  # noqa: E402
from prrs.reliability import GuardedCalculator  # noqa: E402
from prrs.search import _stamped_with_charge  # noqa: E402
from prrs.state import encode  # noqa: E402


def main(out_dir):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=False)
    config = SearchConfig(charge_sensitive=True, total_charge=-1, multiplicity=1)
    atoms = _stamped_with_charge(read(ROOT / "benchmarks/b04_sn2_chloride/input.extxyz"), config)
    graph = encode(atoms, config.bond_scale)
    bonded = {i for e in graph.edges for i in e}
    nucleophile = next(i for i, s in enumerate(atoms.get_chemical_symbols()) if s == "Cl" and i not in bonded)
    carbon = next(i for i, s in enumerate(atoms.get_chemical_symbols()) if s == "C")
    ledger = Ledger(out / "ledger.jsonl", arm="smoke", budget=20000)
    factory = counted_factory(polar1, ledger)
    q0 = float(np.linalg.norm(atoms.positions[nucleophile] - atoms.positions[carbon]))
    steps = int(round((q0 - 1.9) / 0.1))
    settings = ScanSettings(
        deltas_A=tuple(round(-0.1 * k, 10) for k in range(steps + 1)),
        fmax_eV_A=0.01,
        max_steps_per_point=300,
        max_rise_eV=3.0,
    )
    intent = Intent(f"cmp_{nucleophile}_{carbon}", ((nucleophile, carbon),), (1.0,))
    result = constrained_scan(atoms, intent, settings, factory, config, ledger, task="scan")

    # noise: re-relax one mid-scan point from a 0.02 A random perturbation, 3 times
    mid = result.points[len(result.points) // 2]
    rng = np.random.default_rng(0)
    repeats = []
    for r in range(3):
        start = atoms.copy()
        start.positions = mid.positions + rng.normal(scale=0.02, size=mid.positions.shape)
        rep = constrained_scan(start, intent, ScanSettings((0.0,), settings.fmax_eV_A,
                               settings.max_steps_per_point, settings.max_rise_eV),
                               factory, config, ledger, task=f"noise{r}")
        repeats.append(rep.points[0].energy_eV if rep.points else None)
    noise = {"mid_delta_A": mid.delta_A, "energies_eV": repeats, "reference_eV": mid.energy_eV}
    (out / "noise.json").write_text(json.dumps(noise, indent=2))

    triggers = read_scan(result.points, noise_eV=max(1e-4, float(np.ptp([e for e in repeats if e is not None]))))
    released = None
    ridge = next((t for t in triggers if t["kind"] == "ridge_suspected"), None)
    if ridge:
        source = snapshot(atoms, GuardedCalculator(polar1(), config))
        rel = release(result.points[ridge["release"]].positions, atoms, source, factory, config,
                      ledger, task="release", seed=0)
        released = {k: rel[k] for k in ("status", "label", "fingerprint")}
    ledger.close()
    points = [{k: v for k, v in asdict(p).items() if k != "positions"} for p in result.points]
    (out / "scan.json").write_text(json.dumps(
        {"status": result.status, "q0": result.q0, "points": points, "triggers": triggers,
         "release": released, "spent": ledger.spent}, indent=2, default=str))


if __name__ == "__main__":
    main(sys.argv[1])
```

注：噪声重复里随机扰动会改变 q 的起点；`constrained_scan` 以扰动后的 q 为 0 点，所以三次重复测的是“同一 PES 区域、不同起点的约束极小能量离散”，它是能量噪声的上界估计。若实测离散被 q 偏移主导，改为在扰动后先把 q 拉回 mid.q_target（`deltas_A=(0.0, mid.q_target - q_perturbed)`），取第二点能量。冒烟报告里写明用的是哪种。

- [ ] **Step 2: 本地先跑语法与导入检查（不加载模型）**

Run: `cd /home/ruigengji/MLPATH && /home/ruigengji/miniforge3/envs/openmm_dev/bin/python -c "import ast,sys; ast.parse(open('benchmarks/explore_smoke/b04_scan.py').read())"`
Expected: 无输出

- [ ] **Step 3: 查 183 负载**

Run: `ssh korakuen@192.168.0.183 'bash -c "nvidia-smi --query-gpu=index,utilization.gpu,memory.used --format=csv && uptime"'`
Expected: 有一张 GPU 空闲；若都忙，停下报告用户，不抢卡。

- [ ] **Step 4: 在 183 上启动冒烟并记 PID**

```bash
D=benchmarks/explore_smoke/results/$(date +%F)
ssh korakuen@192.168.0.183 "bash -c 'cd /home/ruigengji/MLPATH && mkdir -p $D && CUDA_VISIBLE_DEVICES=<空闲卡号> OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 setsid nohup /home/ruigengji/miniforge3/envs/openmm_dev/bin/python benchmarks/explore_smoke/b04_scan.py $D/run > $D/run.log 2>&1 < /dev/null & echo \$!'" >> $D/pids.txt
```

- [ ] **Step 5: 在 183 上读结果**

Run: `ssh korakuen@192.168.0.183 'bash -c "cd /home/ruigengji/MLPATH/benchmarks/explore_smoke/results/<日期> && tail -5 run.log && python3 -c \"import json;d=json.load(open(\\\"run/scan.json\\\"));print(d[\\\"status\\\"],d[\\\"spent\\\"],d[\\\"triggers\\\"],d[\\\"release\\\"])\""'`
跑超过 15 分钟就按 PID 核对后停下（`ps -o pid,etime,cmd -p <PID> && kill <PID>`），把“冒烟本身不止几分钟”作为结果报告用户。

- [ ] **Step 6: 写冒烟记录与 M2 commit 草稿**

在 `docs/DONE_LOG.md` 末尾追加一节：单点 forward 数（ledger 的 probe 段 forwards ÷ 点数）、单点墙钟、能量噪声实测与测法、是否出脊、释放结果、与 b04 已知 TS 距离的事后对照；写明“冒烟，未做验收”。写 `.commit_drafts/<YYYY-MM-DD_HHMM>_explore-probe.md`（`explore/src/explore/probe.py`、`explore/tests/test_probe.py`、`benchmarks/explore_smoke/b04_scan.py`、`docs/DONE_LOG.md`）。

**M2 之后的决策点（交用户）：** 根据冒烟实测的单点成本和是否出脊，决定 M3 的意图规模与扫描网格，再写 M3 的逐步计划。若 b04 上不出脊，先回到发现层设计 §4 查原因（网格、fmax、求解器换支），不往下叠功能。

---

## 自查记录

- **设计覆盖：** 协议 §2（账本）→ Task 1–2；发现层 §4.1 约束松弛 → Task 3；§4.2 能量截断 → Task 3；§4.3 越脊、跳支反向扫描 → Task 4；§4.3 释放 → Task 5；§4.2 “阈值来自开发噪声” → Task 6。L0、意图、C2、编排、对照臂、参照、PARS 保真、选择比较、零号批 → M3–M10，各自开工前另写逐步计划。
- **类型一致：** `Intent / ScanSettings / ScanPoint / ScanResult / constrained_scan / read_scan / bracket_jump / reverse_check / release`、`Ledger / segment / charge / note_request / instrument / counted_factory` 在各 Task 中一致；`reverse_check` 以 Task 4 Step 3 的签名为准（含 `atoms_template`）。
- **已按源码核实：** `GateRejected.code` 是门代码；`prrs.state.classify` 返回 `(label, fingerprint)`。
