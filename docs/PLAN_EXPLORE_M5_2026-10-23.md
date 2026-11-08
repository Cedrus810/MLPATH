# 发现层 v1 · M5 D2 编排 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 把 M1–M4 的原语和 10-23 的性能原语接成方案 A 的 v1 编排：固定分级（先完成 L1 基础覆盖，再升级）、搜索额与 C2 认证保留额两个账户、暂缓候选的均匀抽查（标 `audit_release`）、kill→resume 不回滚已付成本、不重复计费。

**Architecture:** 新增 `explore/src/explore/journal.py`（运行日志 + 两账户）和 `explore/src/explore/d2.py`（冻结人口表、按波次调度、升级、报告）。一次 run 是一个目录：`manifest.json`、`poses.jsonl`、`population.jsonl`、`journal.jsonl`、`ledger.jsonl`、`report.json`。`run()` 只由这些文件和两个后端决定：同一条命令重跑，就跳过已结算任务，对被打断的任务按账本对账（无法核对就按上限计费），再用新的 attempt 号重做。L1/L2 扫描走 `scan_many` + 常驻 float32 批量后端；释放与 L4 认证走常驻 float64 calculator，经新加的按任务计数视图 `view_factory` 计费。PRRS 源码不动。

**Tech Stack:** Python 3.12（`openmm_dev`，生产环境，不装不升）、ASE 3.29.0、numpy、pytest；ruff 在 `.venv-tools`。

**Spec:** [发现层设计](DISCOVERY_LAYER_DESIGN_2026-10-13.md) §1.1（v1 范围）、§4.2（先粗扫后细化）、§4.3（`ridge_suspected` / `branch_bracketed` / `audit_release`）、§5（L0–L4）、§6.1–6.3 的 v1 行、§9（预算、恢复、失败码）；[评估协议](EVALUATION_PROTOCOL_DESIGN_2026-10-13.md) §2.1–2.3（每结构计费、预留与结算、崩溃按预留上界计入）、§4.1（首次命中成本 = C2 闭合时的累计臂内支出）；上位计划 [PLAN_EXPLORE_V1](PLAN_EXPLORE_V1_2026-10-13.md) M5 行与“M5 起的执行约束”。

状态：**用户 2026-10-23 批准（默认 M5-D1…D6），已施工完成。** 解析台 kill→resume 通过；b04 真实模型冒烟得 c2（与 ts0000 RMSD 3.2e-4 Å）；见 DONE_LOG 末节。noise_eV(b04, fp32) = 1.953125e-3 eV（M5-D5 默认规则）。

## 必读（开工前全部读完）

| 文件 | 为什么要读 |
| --- | --- |
| `CLAUDE.md`（项目根目录） | 设计稿表；`openmm_dev` 是生产环境；183 的结果在 183 上读 |
| `docs/DISCOVERY_LAYER_DESIGN_2026-10-13.md` 全文 | 被实现的对象；v1 做什么、不做什么（§1.1 表） |
| `docs/EVALUATION_PROTOCOL_DESIGN_2026-10-13.md` §2、§4.1 | 计费口径、预留/结算、崩溃与 unreconciled、首次命中成本的定义 |
| `docs/PLAN_EXPLORE_V1_2026-10-13.md` 里程碑表、“M5 起的执行约束”、“性能执行层例外” | 一进程守一张 GPU、calculator 复用、搜索 fp32 / 认证 fp64；`scan_many` 只是执行原语 |
| `docs/DONE_LOG.md` 2026-10-13 至 10-23 各节 | b04 实测：单点 forward、噪声 4.7e-5 eV（fp64）、认证 200 forward、fp32 误差、软模 Hessian 在 fp32 下的未决 |
| `docs/SCAN_THROUGHPUT_2026-10-23.md` | `scan_many` / `MACEBatchEvaluator` 的契约：按结构计费、每个扫描独立账本 |
| `explore/src/explore/{ledger,probe,parallel,certify,intents,poses,streams}.py` | 要组合的全部原语 |

## 与上位计划 / 设计稿不同的地方

1. **出口在旋转不变双鞍点势上验，不在冻结外场台上验。** 原因与 M4 相同：冻结台是外场势，PRRS 的 Hessian 总是投影掉刚体模。复用 M4 的 `explore/tests/twosaddle.py`。
2. **M5 不做 L3 的 MD。** 设计 §1.1 写着“L3 MD 只给抽查项与静态信息不足项”，但 MD 的时窗 t_free 属于 R8（[DECISIONS_PENDING](DECISIONS_PENDING.md) R8）未裁决，评估协议 §1.1 要求跑前冻结。M5 的 L3 只做静态释放（解约束 → 淬火 → 极小确认），抽查项也只走静态释放 + L4。MD 等 R8 定了再另立任务。
3. **§6.1 第 3 条（反复回同盆 / 同鞍点者减少重复额度）M5 不做。** 它要用认证结果反过来调额度，和 v1 的“固定分级、不做认证反馈”边界很近；先做无反馈的固定分级，作为后续消融的基线。重复认证照常付费并单报 `duplicate_c2`。
4. **L2 只含“跳支的反向扫描”，不含“峰附近幅度细化”。** 10-22 b04 冒烟中，认证从 0.01 Å 细扫的峰点出发；M5 直接从 0.1 Å 粗扫的峰点出发认证，细化不做（`follow_min_mode` 本来就是从初猜精修）。若 b04 冒烟（Task 8）出现粗扫峰点精修失败，再报用户加细化。

## 需要你定（每条附默认建议；不定就按默认施工）

| # | 问题 | 默认建议 | 理由 |
| --- | --- | --- | --- |
| M5-D1 | 升级项从哪个点释放、从哪个点认证 | `ridge_suspected`：释放点 = `read_scan` 给的 `release`（峰后一点），认证初猜 = 峰点。`branch_bracketed`：释放点 = 反扫末点（产物支），认证初猜 = 跳前点。抽查项：释放点 = 最后一个收敛点，认证初猜 = 除起点外能量最高的收敛点 | 前两条沿用设计 §4.3 和现有 `read_scan` / `bracket_jump` 的输出；抽查规则是新的冻结规则，只选点、不设门 |
| M5-D2 | 幅度格有一部分越过重叠下限时怎么办 | **截断定义域**：保留从 0 开始、仍在可行域内的最长前缀；剩不到 2 个点才记 `infeasible`；截断如实写进人口表 `domain_truncated` | 10-22 终审提醒：现行 `_feasible` 只要有一个 delta 越界就整条拒，大压缩幅度会让全部短键压缩意图都判 infeasible。截断不新设阈值（下限还是现有 `min_distance_ratio·ΣR`） |
| M5-D3 | `read_scan` 要不要加连续性检查 | **不加**。疑似脊只是触发条件，C2 由 L4 判 | 10-18 复查表明“能量–反力积分差”在陡而连续的脊上会误报（0.1 Å 步长差 0.045 eV，细扫降到 1e-3）；按全局规则 §8 不新增判据 |
| M5-D4 | 释放用哪个精度 | ~~**float64**（与认证共用常驻 calculator）~~ **作废（2026-10-23）**：POLAR-1 权重是 float32，float64 只是升格，见 PLAN_EXPLORE_V1 精度一条 | 释放含极小确认（Hessian）；10-22 实测 fp32 下 2/7 软模结构过不了 PRRS 数值自检，会变成未决 |
| M5-D5 | 正式运行的 `noise_eV` | Task 7 测完后由你定；建议沿用 10-13 的定义：同一点 5 次扰动重复的能量极差，4 个点取最大 | fp32 能量的最小可分辨量约 2 meV（−26132 eV 处 float32 的 ulp ≈ 1.95e-3 eV），fp64 下测得的 4.7e-5 eV 不能沿用 |
| M5-D6 | 先完成全部 L1 再升级（严格分相） | **是** | 设计 §4.2 原文：“在其已选集合的 L1 基础试探覆盖完成后才开放 L2”。代价是首次命中成本偏后，但可归因 |

不阻塞 M5、留到计分（M10）前再定：U1/U3/U4/U5 与 b09 抽样（M5 把它们当输入参数，不写死）；M4 留下的 `same_basin` 按固定标号判还是按对称商判。

## Global Constraints

- 科学环境只用 `/home/ruigengji/miniforge3/envs/openmm_dev/bin/python`；不 pip/conda 安装、卸载、升级。
- PRRS、PARS 源码不动；`explore` 单向 import `prrs`。
- 主预算单位：一个结构在一个模型上真正发起一次 forward 记 1；批量求值一批 m 个结构记 m（协议 §2.1）。
- 恢复不回滚已付成本：被打断的任务若账本里还有没关闭的段，按该任务的上限（cap）计入并标 `reconciled: false`；没有未关闭段则按账本实数计入（协议 §2.3）。
- 不用 hash 判断缓存或续跑能否复用；续跑以任务 ID + journal 为准。
- 不新增会拦住运行或结果的门。所有数值（预算、cap、幅度、fmax、noise_eV、抽查比例）由调用方显式传入，`D2Settings` 没有默认值。
- 随机都走 `explore.streams.stream` 的按键随机流：追加意图不改变已有意图的顺序键和抽查签。
- 输出不写 `/tmp`；真实模型脚本放 `benchmarks/explore_d2/`，结果放 `benchmarks/explore_d2/results/<YYYY-MM-DD>/`。
- 在 NFS 工作区不跑 git；整个里程碑合成一份 `.commit_drafts/<日期_时间>_explore-m5.md`。
- 测试：只跑 `make -C explore check`（本机 CPU）。GPU 命令里不夹带测试。
- Task 7、Task 8 是 GPU 作业：开跑前给计划表（命令、节点、当前 `nvidia-smi`/`ps`、预计耗时及依据、能否续跑、输出目录），等用户说“跑”。

（测试文件约定：各 Task 里“追加”的测试代码所带的 import 一律并到文件顶部，避免 ruff E402 和重复 import。）

## Review Focus

1. **kill 时有多个扫描同时在跑**：每个没结算的任务各自对账；有未关闭段的按 cap 计入，一个都不退款。Task 2、Task 6 测试固定。
2. **续跑重做已结算的任务，或复用旧 attempt 号**：不允许。每个任务最多结算一次，重做的 attempt 号严格递增。Task 4、Task 6 测试固定。
3. **搜索额被扫描吃光后，认证保留额照样能用**：释放被跳过（记 `L3_skipped:search_budget`），认证照常完成。Task 5 测试固定。
4. **kill 留下半行 JSON**（journal 或 ledger 的最后一行没写完）：丢掉这半行后继续；如果丢的是某段的 close，该任务就当未关闭处理，按 cap 计入。Task 1、Task 2 测试固定。
5. **姿态没带电子态**（`atoms.info` 缺 `charge` 或 `spin`）：在 `init_run` 时拒绝，0 forward。Task 3 测试固定。

---

## 文件结构

```
explore/src/explore/
  ledger.py     （改）Ledger(replay=...)、写文件加锁、tally()、CountedView / view_factory
  journal.py    （新）repair_tail、Account、Journal
  probe.py      （改）feasible_prefix()、reverse_plan()（reverse_check 改为调用它，行为不变）
  d2.py         （新）D2Settings、Backends、init_run、open_run、run、报告
explore/tests/
  test_ledger.py    （改）Task 1 的新测试
  test_journal.py   （新）
  test_d2.py        （新）Task 3–5
  test_d2_resume.py （新）Task 6：子进程 kill→resume
  d2_bench.py       （新）双鞍点势上的 D2 测试台（外部计数器、两后端）
  d2_driver.py      （新）kill→resume 子进程入口
benchmarks/explore_d2/
  fp32_noise.py     （新）Task 7
  run_b04_d2.py     （新）Task 8
```

---

### Task 1: 账本——共享文件、按任务对账、常驻模型计数视图

**Files:**
- Modify: `explore/src/explore/ledger.py`
- Test: `explore/tests/test_ledger.py`

**Interfaces:**
- Produces:
  - `Ledger(path, arm, budget=None, replay=True)`：`replay=False` 时不重放已有文件（多个按任务的 Ledger 共用一个文件时用）。所有 Ledger 写文件经同一把模块锁。
  - `tally(path) -> dict[(task, attempt), {"forwards": int, "open": int}]`：已关闭段的 forward 之和，以及还没关闭的段数；文件最后半行（没有换行结尾）忽略。
  - `CountedView(inner, ledger)`：ASE Calculator；每次真实 `calculate` 先 `ledger.charge()`，再调用 `inner.reset()` + `inner.calculate(...)`，复制其 results；有自己的缓存。`inner` 有 `get_hessian` 时也提供 `get_hessian`，每次计 1。
  - `view_factory(inner, ledger) -> Callable[[], CountedView]`：`inner` 已被 `instrument` 过就报 `ValueError`。

- [ ] **Step 1: 写失败的测试**（追加到 `explore/tests/test_ledger.py`）

```python
import threading

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes

from explore.ledger import BudgetExhausted, Ledger, instrument, tally, view_factory


class ChargeAware(Calculator):
    implemented_properties = ["energy", "forces"]

    def __init__(self):
        super().__init__()
        self.calls = 0

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        self.calls += 1
        charge = float(self.atoms.info.get("charge", 0))
        self.results = {
            "energy": charge + float(self.atoms.positions.sum()),
            "forces": np.zeros((len(self.atoms), 3)),
        }


def test_replay_false_starts_from_zero_on_a_shared_file(tmp_path):
    path = tmp_path / "l.jsonl"
    first = Ledger(path, arm="D2", budget=5, replay=False)
    with first.segment(task="a", stage="probe"):
        first.charge()
        first.charge()
    second = Ledger(path, arm="D2", budget=5, replay=False)
    assert second.spent == 0 and second.remaining == 5
    assert Ledger(path, arm="D2").spent == 2


def test_concurrent_task_ledgers_share_one_file(tmp_path):
    path = tmp_path / "l.jsonl"

    def work(k):
        ledger = Ledger(path, arm="D2", budget=1000, replay=False)
        for _ in range(50):
            with ledger.segment(task=f"t{k}", stage="probe"):
                ledger.charge()
        ledger.close()

    threads = [threading.Thread(target=work, args=(k,)) for k in range(8)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    rows = tally(path)
    assert {key: row["forwards"] for key, row in rows.items()} == {
        (f"t{k}", 0): 50 for k in range(8)
    }
    assert all(row["open"] == 0 for row in rows.values())


def test_tally_reports_open_segments_and_ignores_a_torn_last_line(tmp_path):
    path = tmp_path / "l.jsonl"
    ledger = Ledger(path, arm="D2", replay=False)
    with ledger.segment(task="done", stage="probe"):
        ledger.charge()
    ledger._write({"event": "open", "task": "cut", "attempt": 1, "stage": "probe", "reserve": None})
    with path.open("a", encoding="utf-8") as fh:
        fh.write('{"event": "close", "task": "cut"')
    rows = tally(path)
    assert rows[("done", 0)] == {"forwards": 1, "open": 0}
    assert rows[("cut", 1)] == {"forwards": 0, "open": 1}


def test_counted_view_charges_real_forwards_only(tmp_path):
    inner = ChargeAware()
    ledger = Ledger(tmp_path / "l.jsonl", arm="D2", budget=10, replay=False)
    atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.74]])
    atoms.calc = view_factory(inner, ledger)()
    with ledger.segment(task="t", stage="probe") as seg:
        atoms.get_potential_energy()
        atoms.get_forces()
        atoms.positions[1, 2] += 0.1
        atoms.get_potential_energy()
    assert inner.calls == 2 and ledger.spent == 2 and seg.forwards == 2
    assert seg.cache_hits >= 1


def test_views_never_share_a_result_across_electronic_states(tmp_path):
    inner = ChargeAware()
    ledger = Ledger(tmp_path / "l.jsonl", arm="D2", replay=False)
    make = view_factory(inner, ledger)
    a = Atoms("H", positions=[[0, 0, 0]])
    a.info["charge"] = -1
    a.calc = make()
    b = a.copy()
    b.info["charge"] = 0
    b.calc = make()
    assert a.get_potential_energy() == -1.0
    assert b.get_potential_energy() == 0.0
    assert ledger.spent == 2


def test_counted_view_refuses_before_the_model_runs(tmp_path):
    inner = ChargeAware()
    ledger = Ledger(tmp_path / "l.jsonl", arm="D2", budget=1, replay=False)
    atoms = Atoms("H", positions=[[0, 0, 0]])
    atoms.calc = view_factory(inner, ledger)()
    atoms.get_potential_energy()
    atoms.positions += 0.1
    with pytest.raises(BudgetExhausted):
        atoms.get_potential_energy()
    assert inner.calls == 1


def test_view_factory_refuses_an_already_counted_calculator(tmp_path):
    ledger = Ledger(tmp_path / "l.jsonl", arm="D2", replay=False)
    inner = instrument(ChargeAware(), ledger)
    with pytest.raises(ValueError):
        view_factory(inner, ledger)
```

- [ ] **Step 2: 跑测试确认失败**

Run: `make -C explore test`
Expected: FAIL，`ImportError: cannot import name 'tally'`

- [ ] **Step 3: 实现**

`ledger.py` 顶部 import 改为：

```python
import contextlib
import json
import os
import socket
import threading
import time
from pathlib import Path

from ase.calculators.calculator import Calculator, all_changes

_WRITE_LOCK = threading.Lock()
```

`Ledger.__init__` 改签名并按 `replay` 决定是否重放：

```python
    def __init__(self, path, arm, budget=None, replay=True):
        self.path = Path(path)
        self.arm = arm
        self.budget = budget
        self.spent = 0
        self.reserved = 0
        self.unreconciled = []
        self.host = socket.gethostname()
        self._segment = None
        self._unattributed = 0
        if replay and self.path.exists():
            self._replay()
```

`Ledger._write` 加锁：

```python
    def _write(self, event):
        record = {"arm": self.arm, "host": self.host, "t": time.time(), **event}
        line = json.dumps(record) + "\n"
        with _WRITE_LOCK:
            self.path.parent.mkdir(parents=True, exist_ok=True)
            with self.path.open("a", encoding="utf-8") as handle:
                handle.write(line)
                handle.flush()
                os.fsync(handle.fileno())
```

文件末尾追加：

```python
def tally(path):
    """Closed forwards and still-open segments per (task, attempt) in a ledger file.

    For a file shared by many per-task ``Ledger(replay=False)`` instances (discovery D2).
    ``open > 0`` means the attempt was interrupted inside a segment, so ``forwards`` is
    only a lower bound. A last line without its newline was torn by a kill and is
    ignored: it was never acknowledged."""
    path = Path(path)
    rows = {}
    if not path.exists():
        return rows
    text = path.read_text(encoding="utf-8")
    lines = text.split("\n")
    if not text.endswith("\n"):
        lines = lines[:-1]
    for line in lines:
        if not line.strip():
            continue
        event = json.loads(line)
        row = rows.setdefault((event["task"], event["attempt"]), {"forwards": 0, "open": 0})
        if event["event"] == "open":
            row["open"] += 1
        elif event["event"] == "close":
            row["forwards"] += event["forwards"]
            if event["stage"] != "unattributed":
                row["open"] -= 1
    return rows


class CountedView(Calculator):
    """A per-caller counted window onto one resident physical calculator.

    ``reused_counted_factory`` binds the resident model to one ledger; D2 needs one
    ledger per task attempt (``scan_many`` and per-task caps), so the model is shared
    and only the counting is per caller. Each view has its own ASE cache, and a fresh
    view is made per factory call, so two callers that differ only in ``atoms.info``
    (charge, spin) never share a result. The forward is charged before the model runs."""

    def __init__(self, inner, ledger):
        super().__init__()
        self.inner = inner
        self.ledger = ledger
        self.implemented_properties = list(inner.implemented_properties)
        if hasattr(inner, "get_hessian"):
            self.get_hessian = self._get_hessian

    def get_property(self, name, atoms=None, allow_calculation=True):
        self.ledger.note_request()
        before = self.ledger.spent
        value = super().get_property(name, atoms, allow_calculation)
        if self.ledger.spent == before:
            self.ledger.note_cache_hit()
        return value

    def calculate(self, atoms=None, properties=("energy", "forces"), system_changes=all_changes):
        self.ledger.charge()
        super().calculate(atoms, properties, system_changes)
        self.inner.reset()
        self.inner.calculate(self.atoms, properties, system_changes)
        self.results = dict(self.inner.results)

    def _get_hessian(self, *args, **kwargs):
        self.ledger.charge()
        return self.inner.get_hessian(*args, **kwargs)


def view_factory(inner, ledger):
    """Factory of ``CountedView``s on ``inner`` charged to ``ledger``."""
    if getattr(inner, "_explore_ledger", None) is not None:
        raise ValueError("calculator is already counted by a ledger")

    def make():
        return CountedView(inner, ledger)

    return make
```

- [ ] **Step 4: 跑测试确认通过**

Run: `make -C explore check`
Expected: fmt-check、lint 通过，全部测试 PASS（原 118 个 + 本 Task 7 个）

---

### Task 2: 运行日志与两账户

**Files:**
- Create: `explore/src/explore/journal.py`
- Test: `explore/tests/test_journal.py`

**Interfaces:**
- Produces:
  - `repair_tail(path) -> bool`：文件最后一行没有换行结尾（kill 写了一半）就截掉这半行，返回是否截过。
  - `Account(name, budget)`：`spent`（已结算实数 + 对账计入）、`held`（未结算的预留）、`available = budget − spent − held`、`hold(cap) -> bool`、`settle(cap, spent)`、`charge(amount)`。
  - `Journal(path)`：打开时先 `repair_tail`，再把全部事件读进 `events`；`append(event)` 写一行、fsync、追加到 `events`。

- [ ] **Step 1: 写失败的测试**（`explore/tests/test_journal.py`）

```python
from explore.journal import Account, Journal, repair_tail


def test_account_holds_caps_and_settles_actual():
    acc = Account("search", 100)
    assert acc.hold(60)
    assert not acc.hold(50)
    acc.settle(60, 7)
    assert (acc.spent, acc.held, acc.available) == (7, 0, 93)
    acc.charge(10)
    assert acc.available == 83


def test_journal_drops_a_torn_last_line_and_appends_cleanly(tmp_path):
    path = tmp_path / "journal.jsonl"
    Journal(path).append({"event": "dispatch", "task": "a"})
    with path.open("a", encoding="utf-8") as fh:
        fh.write('{"event": "sett')
    again = Journal(path)
    assert again.events == [{"event": "dispatch", "task": "a"}]
    again.append({"event": "settle", "task": "a"})
    assert [e["event"] for e in Journal(path).events] == ["dispatch", "settle"]


def test_repair_tail_leaves_a_complete_file_alone(tmp_path):
    path = tmp_path / "x.jsonl"
    path.write_text('{"a": 1}\n', encoding="utf-8")
    assert repair_tail(path) is False
    assert path.read_text(encoding="utf-8") == '{"a": 1}\n'
    assert repair_tail(tmp_path / "missing.jsonl") is False
```

- [ ] **Step 2: 跑测试确认失败**

Run: `make -C explore test`
Expected: FAIL，`ModuleNotFoundError: No module named 'explore.journal'`

- [ ] **Step 3: 实现** `explore/src/explore/journal.py`

```python
"""D2 run journal and the two v1 budget accounts (DISCOVERY_LAYER_DESIGN §6.2 v1 row;
EVALUATION_PROTOCOL §2.2-2.3).

The journal is the coordinator's commit order: dispatch (a cap is held), settle (the
real spend replaces the hold), reconcile (an interrupted attempt is charged), skip (a
task the account could not hold). Accounts are rebuilt from it on every start."""

from __future__ import annotations

import json
import os
from pathlib import Path


def repair_tail(path):
    """Drop a last line torn by a kill (no trailing newline); keep everything before it."""
    path = Path(path)
    if not path.exists():
        return False
    data = path.read_bytes()
    if not data or data.endswith(b"\n"):
        return False
    path.write_bytes(data[: data.rfind(b"\n") + 1])
    return True


class Account:
    def __init__(self, name, budget):
        self.name = name
        self.budget = int(budget)
        self.spent = 0
        self.held = 0

    @property
    def available(self):
        return self.budget - self.spent - self.held

    def hold(self, cap):
        if cap > self.available:
            return False
        self.held += cap
        return True

    def settle(self, cap, spent):
        self.held -= cap
        self.spent += spent

    def charge(self, amount):
        self.spent += amount


class Journal:
    def __init__(self, path):
        self.path = Path(path)
        repair_tail(self.path)
        self.events = []
        if self.path.exists():
            for line in self.path.read_text(encoding="utf-8").splitlines():
                if line.strip():
                    self.events.append(json.loads(line))

    def append(self, event):
        line = json.dumps(event) + "\n"
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with self.path.open("a", encoding="utf-8") as handle:
            handle.write(line)
            handle.flush()
            os.fsync(handle.fileno())
        self.events.append(event)
```

- [ ] **Step 4: 跑测试确认通过**

Run: `make -C explore check`
Expected: 全部 PASS

---

### Task 3: D2 冻结——人口表与定义域截断

**Files:**
- Modify: `explore/src/explore/probe.py`（加 `feasible_prefix`、`reverse_plan`；`reverse_check` 改用 `reverse_plan`）
- Create: `explore/src/explore/d2.py`（本 Task 只写 `D2Settings`、`Backends`、`init_run`）
- Create: `explore/tests/d2_bench.py`
- Test: `explore/tests/test_d2.py`

**Interfaces:**
- Consumes: `IntentSpec`（`intents.py`）、`Intent`、`_feasible`（`probe.py`）、`stream`（`streams.py`）、`SearchConfig.to_dict()`。
- Produces:
  - `feasible_prefix(intent, q0, deltas, numbers, config) -> tuple[float, ...]`
  - `reverse_plan(points, jump, settings) -> (start_positions: np.ndarray, back: ScanSettings)`
  - `D2Settings`（全部字段必填，见下）
  - `Backends(search_batch, certify_calculator)`
  - `init_run(run_dir, poses, specs_by_pose, settings, config, identity) -> None`；`poses` 是 `[(pose_id, Atoms)]`，`atoms.info` 必须有 `charge`、`spin`；`specs_by_pose[pose_id]` 是 `IntentSpec` 列表（由调用方按 U4/U5 规则生成，规则参数写进 `identity`）。
  - `population.jsonl` 每行：`intent_id, pose_id, pairs, coefs, sign, family, layer, scope, deltas, domain_truncated, sampled, status ∈ {ready, infeasible, not_sampled}, order=[position, pose_key], audit`；行按 `order` 排序（轮转：先每个姿态的第 0 个，再第 1 个…）。

- [ ] **Step 1: 测试台** `explore/tests/d2_bench.py`

```python
"""TwoSaddle bench for the D2 tests and the kill->resume driver.

Pose: H at the left minimum (-1, 0). H-O1 compression and H-O0 stretch both drive H
over to the right minimum; O0-O1 stretch only loads the anchor frame (no trigger)."""

import os

from ase.calculators.calculator import all_changes

from explore.d2 import Backends, D2Settings, init_run
from explore.intents import IntentSpec
from explore.probe import Intent
from prrs.config import SearchConfig
from twosaddle import TwoSaddleInvariant, two_saddle_atoms

CONFIG = SearchConfig(closed_shell_only=False)
POSE = "p0"


def spec(pair, sign, pose=POSE, layer="near"):
    iid = f"{pose}|distance:{pair[0]}-{pair[1]}:{'+' if sign > 0 else '-'}"
    return IntentSpec(Intent(iid, (pair,), (1.0,)), sign, "distance", layer, "inter")


SPECS = (spec((1, 3), -1), spec((0, 3), +1), spec((0, 1), +1))


def settings(**over):
    base = dict(
        seed=7,
        search_budget=20000,
        certify_budget=20000,
        scan_cap=600,
        reverse_cap=200,
        release_cap=600,
        certify_cap=3000,
        per_pose_quota=10,
        audit_fraction=0.0,
        noise_eV=1e-6,
        certify_max_rise_eV=5.0,
        magnitudes_A=tuple(round(0.1 * k, 10) for k in range(1, 18)),
        fmax_eV_A=1e-3,
        max_steps_per_point=300,
        max_rise_eV=3.0,
        max_active=2,
        wave_size=4,
    )
    base.update(over)
    return D2Settings(**base)


def pose():
    atoms = two_saddle_atoms(-1.0, 0.0)
    atoms.info.update(charge=0, spin=1)
    return atoms


def init(run_dir, specs=SPECS, **over):
    init_run(
        run_dir,
        [(POSE, pose())],
        {POSE: list(specs)},
        settings(**over),
        CONFIG,
        identity={"bench": "twosaddle"},
    )


class Tick:
    """External evaluation counter: one byte per model evaluation, written before it
    runs; evaluation number ``crash_after`` (in this process) kills it like SIGKILL."""

    def __init__(self, path=None, crash_after=-1):
        self.path, self.crash_after, self.count = path, crash_after, 0

    def __call__(self):
        self.count += 1
        if self.path is not None:
            with open(self.path, "ab") as fh:
                fh.write(b".")
        if self.count == self.crash_after:
            os._exit(137)


class CountedTwoSaddle(TwoSaddleInvariant):
    def __init__(self, tick):
        super().__init__()
        self.tick = tick

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        self.tick()
        super().calculate(atoms, properties, system_changes)


def backends(tick):
    def batch(frames):
        out = []
        for frame in frames:
            tick()
            atoms = frame.copy()
            atoms.calc = TwoSaddleInvariant()
            out.append((atoms.get_potential_energy(), atoms.get_forces()))
        return out

    return Backends(search_batch=batch, certify_calculator=CountedTwoSaddle(tick))
```

- [ ] **Step 2: 写失败的测试**（`explore/tests/test_d2.py`）

```python
import json

import pytest
from ase.data import covalent_radii

from d2_bench import CONFIG, POSE, SPECS, init, pose, settings, spec
from explore.d2 import init_run


def population(run_dir):
    lines = (run_dir / "population.jsonl").read_text(encoding="utf-8").splitlines()
    return [json.loads(line) for line in lines if line.strip()]


def two_pose_specs():
    pairs = (((1, 3), -1), ((0, 3), +1), ((0, 1), +1))
    return {p: [spec(pair, sign, pose=p) for pair, sign in pairs] for p in ("p0", "p1")}


def test_population_is_round_robin_keyed_and_quota_limited(tmp_path):
    poses = [("p0", pose()), ("p1", pose())]
    specs = two_pose_specs()
    cfg = settings(per_pose_quota=2, audit_fraction=0.5)
    init_run(tmp_path / "a", poses, specs, cfg, CONFIG, identity={})
    rows = population(tmp_path / "a")
    sampled = [r for r in rows if r["sampled"]]
    assert [r["order"][0] for r in sampled] == [0, 0, 1, 1]
    assert {r["pose_id"] for r in sampled[:2]} == {"p0", "p1"}
    assert sum(r["status"] == "not_sampled" for r in rows) == 2
    # appending an intent moves no existing order key and no audit draw
    more = {p: v + [spec((0, 2), +1, pose=p)] for p, v in specs.items()}
    init_run(tmp_path / "b", poses, more, cfg, CONFIG, identity={})
    old = {r["intent_id"]: r for r in rows}
    for r in population(tmp_path / "b"):
        if r["intent_id"] in old:
            assert r["audit"] == old[r["intent_id"]]["audit"]
            assert r["order"][1] == old[r["intent_id"]]["order"][1]
    old_ids = [r["intent_id"] for r in rows if r["pose_id"] == "p0"]
    new_ids = [r["intent_id"] for r in population(tmp_path / "b") if r["pose_id"] == "p0"]
    assert [i for i in new_ids if i in old_ids] == old_ids


def test_domain_is_truncated_at_the_overlap_floor(tmp_path):
    # H-O1 starts at 3.2146 A; the floor is min_distance_ratio * (R_O + R_H)
    floor = CONFIG.min_distance_ratio * (covalent_radii[8] + covalent_radii[1])
    init(tmp_path / "r", specs=(spec((1, 3), -1),), magnitudes_A=(0.5, 1.0, 2.0, 3.0))
    (row,) = population(tmp_path / "r")
    q0 = 3.2146
    assert row["status"] == "ready" and row["domain_truncated"]
    assert row["deltas"] == [0.0, -0.5, -1.0, -2.0]
    assert q0 - 3.0 <= floor < q0 - 2.0


def test_a_domain_with_no_feasible_step_is_infeasible(tmp_path):
    init(tmp_path / "r", specs=(spec((1, 3), -1),), magnitudes_A=(3.0, 3.1))
    (row,) = population(tmp_path / "r")
    assert row["status"] == "infeasible" and row["deltas"] == [0.0]


def test_init_refuses_a_pose_without_electronic_state_and_an_existing_run(tmp_path):
    bare = pose()
    del bare.info["spin"]
    with pytest.raises(ValueError, match="spin"):
        init_run(tmp_path / "x", [(POSE, bare)], {POSE: list(SPECS)}, settings(), CONFIG, {})
    assert not (tmp_path / "x" / "manifest.json").exists()
    init(tmp_path / "r")
    with pytest.raises(FileExistsError):
        init(tmp_path / "r")


def test_settings_have_no_silent_values():
    with pytest.raises(ValueError):
        settings(audit_fraction=1.5)
    with pytest.raises(ValueError):
        settings(magnitudes_A=(0.2, 0.1))
    with pytest.raises(ValueError):
        settings(scan_cap=0)
```

- [ ] **Step 3: 跑测试确认失败**

Run: `make -C explore test`
Expected: FAIL，`ModuleNotFoundError: No module named 'explore.d2'`

- [ ] **Step 4: `probe.py` 加两个纯函数**（放在 `_feasible` 之后；`reverse_check` 改为调用 `reverse_plan`）

```python
def feasible_prefix(intent, q0, deltas, numbers, config):
    """The longest leading part of ``deltas`` inside the overlap-floor domain (design
    §3.2). ``deltas`` start at 0 and move away monotonically, so once a target is
    infeasible every later one is too. Fewer than two points left means no scan."""
    k = len(deltas)
    while k > 0 and not _feasible(intent, q0, deltas[:k], numbers, config):
        k -= 1
    return tuple(deltas[:k])


def reverse_plan(points, jump, settings):
    """Start positions and settings of the reverse scan that checks a jump for
    hysteresis: from the post-jump point back to the pre-jump target."""
    before, after = points[jump["before"]], points[jump["after"]]
    back = ScanSettings(
        deltas_A=(0.0, before.delta_A - after.delta_A),
        fmax_eV_A=settings.fmax_eV_A,
        max_steps_per_point=settings.max_steps_per_point,
        max_rise_eV=settings.max_rise_eV,
        reuse_curvature=settings.reuse_curvature,
    )
    return after.positions.copy(), back
```

`reverse_check` 函数体前半改为：

```python
    positions, back = reverse_plan(points, jump, settings)
    start = atoms_template.copy()
    start.positions = positions
    result = constrained_scan(
        start,
        intent,
        back,
        factory,
        config,
        ledger,
        task,
        attempt=attempt,
        reserve=reserve,
        artifact=artifact,
    )
```

（其余不变；`before`、`after` 两个局部变量随之删掉。）

- [ ] **Step 5: 写 `explore/src/explore/d2.py`（本 Task 部分）**

```python
"""D2 orchestration, discovery v1 (DISCOVERY_LAYER_DESIGN_2026-10-13 §1.1, §4.2-4.3,
§5, §6.1-6.3 v1 rows, §9; EVALUATION_PROTOCOL_DESIGN_2026-10-13 §2.2-2.3, §4.1).

A run directory holds manifest.json (frozen settings, config, identity), poses.jsonl,
population.jsonl (the intent table: order keys, audit draws, scan domains),
journal.jsonl (dispatch / settle / reconcile / skip in coordinator order) and
ledger.jsonl (every forward, one Ledger per task attempt). ``run`` is a function of
these files and the two backends: rerunning it after a crash skips settled tasks,
charges interrupted attempts (their cap when a segment was left open) and redoes them
under a new attempt number. Nothing is reused by hash."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from .probe import feasible_prefix
from .streams import stream

ARM = "D2"
RULE_VERSION = "D2-v1"


@dataclass(frozen=True)
class D2Settings:
    seed: int
    search_budget: int
    certify_budget: int
    scan_cap: int
    reverse_cap: int
    release_cap: int
    certify_cap: int
    per_pose_quota: int
    audit_fraction: float
    noise_eV: float
    certify_max_rise_eV: float
    magnitudes_A: tuple[float, ...]
    fmax_eV_A: float
    max_steps_per_point: int
    max_rise_eV: float
    max_active: int
    wave_size: int

    def __post_init__(self):
        for name in (
            "search_budget",
            "certify_budget",
            "scan_cap",
            "reverse_cap",
            "release_cap",
            "certify_cap",
            "per_pose_quota",
            "max_active",
            "wave_size",
        ):
            if getattr(self, name) < 1:
                raise ValueError(f"{name} must be positive")
        if not 0.0 <= self.audit_fraction <= 1.0:
            raise ValueError("audit_fraction must be in [0, 1]")
        mags = list(self.magnitudes_A)
        if not mags or any(m <= 0 for m in mags) or mags != sorted(mags):
            raise ValueError("magnitudes_A must be positive and ascending")


@dataclass
class Backends:
    """``search_batch``: ``evaluate_batch`` for ``scan_many`` (resident float32 model in
    production). ``certify_calculator``: one resident physical calculator for release
    and L4 (float64 in production)."""

    search_batch: object
    certify_calculator: object


def _key(seed, *parts):
    return float(stream(seed, *parts).random())


def init_run(run_dir, poses, specs_by_pose, settings, config, identity):
    """Freeze a D2 run: write poses, the intent table and the manifest (last, so its
    presence marks a complete freeze). Refuses before any evaluation when a pose lacks
    its electronic state."""
    run_dir = Path(run_dir)
    if (run_dir / "manifest.json").exists():
        raise FileExistsError(f"{run_dir} already holds a D2 run; resume it with run()")
    for pose_id, atoms in poses:
        missing = [k for k in ("charge", "spin") if k not in atoms.info]
        if missing:
            raise ValueError(f"pose {pose_id} lacks atoms.info {missing}")
    rows = []
    for pose_id, atoms in poses:
        pose_key = _key(settings.seed, "pose_order", pose_id)
        specs = sorted(
            specs_by_pose.get(pose_id, ()),
            key=lambda sp: _key(settings.seed, "population", sp.intent.intent_id),
        )
        for position, spec in enumerate(specs):
            intent = spec.intent
            wanted = spec.deltas(settings.magnitudes_A)
            q0 = intent.value(atoms.positions)
            deltas = feasible_prefix(intent, q0, wanted, atoms.numbers, config)
            sampled = position < settings.per_pose_quota
            if not sampled:
                status = "not_sampled"
            else:
                status = "ready" if len(deltas) >= 2 else "infeasible"
            rows.append(
                {
                    "intent_id": intent.intent_id,
                    "pose_id": pose_id,
                    "pairs": [list(p) for p in intent.pairs],
                    "coefs": list(intent.coefs),
                    "sign": spec.sign,
                    "family": spec.family,
                    "layer": spec.layer,
                    "scope": spec.scope,
                    "deltas": [float(d) for d in deltas],
                    "domain_truncated": len(deltas) < len(wanted),
                    "sampled": sampled,
                    "status": status,
                    "order": [position, pose_key],
                    "audit": _key(settings.seed, "audit", intent.intent_id)
                    < settings.audit_fraction,
                }
            )
    rows.sort(key=lambda r: tuple(r["order"]))
    run_dir.mkdir(parents=True, exist_ok=True)
    with (run_dir / "poses.jsonl").open("w", encoding="utf-8") as fh:
        for pose_id, atoms in poses:
            record = {
                "pose_id": pose_id,
                "numbers": atoms.numbers.tolist(),
                "positions": atoms.positions.tolist(),
                "info": {"charge": int(atoms.info["charge"]), "spin": int(atoms.info["spin"])},
            }
            fh.write(json.dumps(record) + "\n")
    with (run_dir / "population.jsonl").open("w", encoding="utf-8") as fh:
        for row in rows:
            fh.write(json.dumps(row) + "\n")
    manifest = {
        "rule_version": RULE_VERSION,
        "arm": ARM,
        "settings": asdict(settings),
        "config": config.to_dict(),
        "identity": identity,
    }
    (run_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
```

- [ ] **Step 6: 跑测试确认通过**

Run: `make -C explore check`
Expected: 全部 PASS（`test_probe.py` 里 `reverse_check` 的原有测试也必须照旧通过）

---

### Task 4: D2 执行——L1/L2 扫描波次与对账恢复

**Files:**
- Modify: `explore/src/explore/d2.py`
- Test: `explore/tests/test_d2.py`

**Interfaces:**
- Consumes: Task 1 `Ledger(replay=False)`、`tally`；Task 2 `Account`、`Journal`、`repair_tail`；Task 3 `reverse_plan`；`scan_many`、`ScanJob`、`read_scan`、`ScanPoint`。
- Produces:
  - `open_run(run_dir) -> _Run`：读全部文件，重建账户，**对每个“已 dispatch 未 settle”的 attempt 写一条 reconcile 事件**（有副作用）。
  - `run(run_dir, backends) -> dict`：本 Task 中依次做 L1、L2，写 `report.json` 并返回报告（Task 5 再加升级）。
  - journal 事件字段：
    - `dispatch {task, attempt, account, cap, kind}`
    - `settle {task, attempt, account, cap, spent, status, payload, cumulative_spent}`
    - `reconcile {task, attempt, account, cap, spent, reconciled, known_forwards}`
    - `skip {task, reason, cumulative_spent}`
  - 任务 ID：`L1|<intent_id>`、`L2|<intent_id>`、`L3|<intent_id>`、`L4|<intent_id>`。

- [ ] **Step 1: 写失败的测试**（追加到 `test_d2.py`）

```python
from d2_bench import Tick, backends
from explore.d2 import open_run, run
from explore.journal import Journal
from explore.ledger import Ledger


def journal(run_dir):
    return Journal(run_dir / "journal.jsonl").events


def test_l1_scans_every_ready_intent_once_and_pays_exactly(tmp_path):
    init(tmp_path / "r")
    tick = Tick()
    rep = run(tmp_path / "r", backends(tick))
    settles = [e for e in journal(tmp_path / "r") if e["event"] == "settle"]
    l1 = [e for e in settles if e["task"].startswith("L1|")]
    assert len(l1) == 3 and all(e["attempt"] == 0 for e in l1)
    assert all(e["spent"] <= 600 for e in l1)
    assert rep["accounts"]["search"]["held"] == 0
    # D2 adds nothing to the count: model calls == ledger == accounts
    paid = rep["accounts"]["search"]["spent"] + rep["accounts"]["certify"]["spent"]
    assert tick.count == rep["ledger_total"] == paid


def test_coverage_truncation_is_recorded_not_hidden(tmp_path):
    init(tmp_path / "r", search_budget=600)
    rep = run(tmp_path / "r", backends(Tick()))
    events = journal(tmp_path / "r")
    assert sum(e["event"] == "settle" and e["task"].startswith("L1|") for e in events) == 1
    assert rep["skips"].get("L1_skipped:search_budget") == 2
    assert rep["accounts"]["search"]["spent"] <= 600


def test_interrupted_attempts_are_charged_and_retried(tmp_path):
    run_dir = tmp_path / "r"
    init(run_dir, specs=SPECS[:2])
    rows = population(run_dir)
    a, b = (f"L1|{r['intent_id']}" for r in rows[:2])
    j = Journal(run_dir / "journal.jsonl")
    for task in (a, b):
        j.append(
            {"event": "dispatch", "task": task, "attempt": 0, "account": "search",
             "cap": 600, "kind": "L1"}
        )
    led = Ledger(run_dir / "ledger.jsonl", arm="D2", replay=False)
    with led.segment(task=a, stage="probe"):
        for _ in range(37):
            led.charge()
    led._write({"event": "open", "task": b, "attempt": 0, "stage": "probe", "reserve": None})
    opened = open_run(run_dir)
    rec = {e["task"]: e for e in opened.journal.events if e["event"] == "reconcile"}
    assert (rec[a]["spent"], rec[a]["reconciled"]) == (37, True)
    assert (rec[b]["spent"], rec[b]["reconciled"]) == (600, False)
    assert opened.accounts["search"].spent == 637
    rep = run(run_dir, backends(Tick()))
    settled = {(e["task"], e["attempt"]) for e in journal(run_dir) if e["event"] == "settle"}
    assert (a, 1) in settled and (b, 1) in settled
    assert len([e for e in journal(run_dir) if e["event"] == "reconcile"]) == 2
    assert rep["unreconciled"][0]["task"] == b
```

- [ ] **Step 2: 跑测试确认失败**

Run: `make -C explore test`
Expected: FAIL，`ImportError: cannot import name 'open_run'`

- [ ] **Step 3: 实现**

`d2.py` 的 import 块替换为：

```python
from __future__ import annotations

import json
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from ase import Atoms

from .journal import Account, Journal, repair_tail
from .ledger import Ledger, tally
from .parallel import ScanJob, scan_many
from .poses import _jsonable
from .probe import (
    Intent,
    ScanPoint,
    ScanSettings,
    feasible_prefix,
    read_scan,
    reverse_plan,
)
from .streams import stream
```

在 `init_run` 之后追加：

```python
def _task(kind, row):
    return f"{kind}|{row['intent_id']}"


def _seed(settings, task):
    return int(stream(settings.seed, "task", task).integers(2**31 - 1))


def _intent(row):
    return Intent(row["intent_id"], tuple(tuple(p) for p in row["pairs"]), tuple(row["coefs"]))


def _scan_settings(settings, deltas):
    return ScanSettings(
        tuple(float(d) for d in deltas),
        settings.fmax_eV_A,
        settings.max_steps_per_point,
        settings.max_rise_eV,
        reuse_curvature=True,
    )


def _scan_payload(outcome):
    if outcome.error is not None:
        return f"error:{outcome.error.split(':')[0]}", {"error": outcome.error}
    result = outcome.result
    rejected = result.rejected_positions
    return result.status, {
        "q0": result.q0,
        "points": [_jsonable(vars(p)) for p in result.points],
        "rejected_positions": None if rejected is None else rejected.tolist(),
    }


def _points(payload):
    out = []
    for d in payload.get("points", []):
        forces = d["forces_pes"]
        out.append(
            ScanPoint(
                **{
                    **d,
                    "positions": np.asarray(d["positions"], float),
                    "forces_pes": None if forces is None else np.asarray(forces, float),
                }
            )
        )
    return out


class _Run:
    """Run state rebuilt from the run directory. Construction reconciles every
    dispatched-but-unsettled attempt and journals it (see ``open_run``)."""

    def __init__(self, run_dir):
        from prrs.config import SearchConfig

        self.dir = Path(run_dir)
        manifest = json.loads((self.dir / "manifest.json").read_text(encoding="utf-8"))
        raw = dict(manifest["settings"])
        raw["magnitudes_A"] = tuple(raw["magnitudes_A"])
        self.settings = D2Settings(**raw)
        self.config = SearchConfig.from_dict(manifest["config"])
        self.poses = {}
        for line in (self.dir / "poses.jsonl").read_text(encoding="utf-8").splitlines():
            if line.strip():
                r = json.loads(line)
                atoms = Atoms(numbers=r["numbers"], positions=r["positions"])
                atoms.info.update(r["info"])
                self.poses[r["pose_id"]] = atoms
        self.population = [
            json.loads(line)
            for line in (self.dir / "population.jsonl").read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
        self.ledger_path = self.dir / "ledger.jsonl"
        repair_tail(self.ledger_path)
        self.journal = Journal(self.dir / "journal.jsonl")
        self.accounts = {
            "search": Account("search", self.settings.search_budget),
            "certify": Account("certify", self.settings.certify_budget),
        }
        self.settled, self.skipped, self.attempts = {}, {}, {}
        unsettled = {}
        for ev in self.journal.events:
            kind = ev["event"]
            if kind == "dispatch":
                self.attempts[ev["task"]] = max(self.attempts.get(ev["task"], 0), ev["attempt"] + 1)
                unsettled[(ev["task"], ev["attempt"])] = ev
            elif kind in ("settle", "reconcile"):
                unsettled.pop((ev["task"], ev["attempt"]), None)
                self.accounts[ev["account"]].charge(ev["spent"])
                if kind == "settle":
                    self.settled[ev["task"]] = ev
            elif kind == "skip":
                self.skipped[ev["task"]] = ev["reason"]
        self._reconcile(unsettled)

    def _reconcile(self, unsettled):
        counts = tally(self.ledger_path)
        for (task, attempt), ev in unsettled.items():
            row = counts.get((task, attempt), {"forwards": 0, "open": 0})
            exact = row["open"] == 0
            spent = row["forwards"] if exact else ev["cap"]
            self.journal.append(
                {
                    "event": "reconcile",
                    "task": task,
                    "attempt": attempt,
                    "account": ev["account"],
                    "cap": ev["cap"],
                    "spent": spent,
                    "reconciled": exact,
                    "known_forwards": row["forwards"],
                }
            )
            self.accounts[ev["account"]].charge(spent)

    def total_spent(self):
        return sum(a.spent for a in self.accounts.values())

    def task_ledger(self, cap):
        return Ledger(self.ledger_path, arm=ARM, budget=cap, replay=False)

    def dispatch(self, task, account, cap, kind):
        if not self.accounts[account].hold(cap):
            return None
        attempt = self.attempts.get(task, 0)
        self.attempts[task] = attempt + 1
        self.journal.append(
            {"event": "dispatch", "task": task, "attempt": attempt, "account": account,
             "cap": cap, "kind": kind}
        )
        return attempt

    def settle(self, task, attempt, account, cap, ledger, status, payload):
        ledger.close()
        self.accounts[account].settle(cap, ledger.spent)
        event = {
            "event": "settle",
            "task": task,
            "attempt": attempt,
            "account": account,
            "cap": cap,
            "spent": ledger.spent,
            "status": status,
            "payload": payload,
            "cumulative_spent": self.total_spent(),
        }
        self.journal.append(event)
        self.settled[task] = event

    def skip(self, task, reason):
        if task in self.skipped:
            return
        self.journal.append(
            {"event": "skip", "task": task, "reason": reason,
             "cumulative_spent": self.total_spent()}
        )
        self.skipped[task] = reason

    def pending(self, task):
        return task not in self.settled and task not in self.skipped


def open_run(run_dir):
    """Load a frozen run and reconcile interrupted attempts (journals reconcile events)."""
    return _Run(run_dir)


def _run_scans(run, backends, items, cap, kind):
    """Batched scans in waves of ``wave_size``; each job holds its cap from the search
    account before dispatch and settles its real spend after the wave."""
    s = run.settings
    todo = list(items)
    while todo:
        wave = []
        for task, atoms, intent, scan_settings in todo[: s.wave_size]:
            attempt = run.dispatch(task, "search", cap, kind)
            if attempt is None:
                break
            job = ScanJob(
                atoms, intent, scan_settings, run.config, run.task_ledger(cap), task,
                attempt=attempt,
            )
            wave.append((task, attempt, job))
        if not wave:
            for task, *_ in todo:
                run.skip(task, f"{kind}_skipped:search_budget")
            return
        outcomes, _ = scan_many([job for _, _, job in wave], backends.search_batch, s.max_active)
        for (task, attempt, job), outcome in zip(wave, outcomes):
            status, payload = _scan_payload(outcome)
            run.settle(task, attempt, "search", cap, job.ledger, status, payload)
        todo = todo[len(wave) :]


def _l1(run, backends):
    s = run.settings
    items = []
    for row in run.population:
        task = _task("L1", row)
        if row["status"] == "ready" and run.pending(task):
            atoms = run.poses[row["pose_id"]].copy()
            items.append((task, atoms, _intent(row), _scan_settings(s, row["deltas"])))
    _run_scans(run, backends, items, s.scan_cap, "L1")


def _first(triggers, kind):
    return next((t for t in triggers if t["kind"] == kind), None)


def _l2(run, backends):
    s = run.settings
    items = []
    for row in run.population:
        ev = run.settled.get(_task("L1", row))
        task = _task("L2", row)
        if ev is None or not run.pending(task):
            continue
        points = _points(ev["payload"])
        triggers = read_scan(points, s.noise_eV)
        jump = _first(triggers, "jump")
        if _first(triggers, "ridge_suspected") is not None or jump is None:
            continue
        start, back = reverse_plan(points, jump, _scan_settings(s, row["deltas"]))
        atoms = run.poses[row["pose_id"]].copy()
        atoms.positions = start
        items.append((task, atoms, _intent(row), back))
    _run_scans(run, backends, items, s.reverse_cap, "L2")


def _counts(run, prefix):
    return dict(Counter(ev["status"] for t, ev in run.settled.items() if t.startswith(prefix)))


def _report(run):
    rows = run.population
    population = dict(Counter(r["status"] for r in rows))
    population["domain_truncated"] = sum(1 for r in rows if r["sampled"] and r["domain_truncated"])
    events = run.journal.events
    return {
        "rule_version": RULE_VERSION,
        "population": population,
        "l1_status": _counts(run, "L1|"),
        "l2_status": _counts(run, "L2|"),
        "skips": dict(Counter(e["reason"] for e in events if e["event"] == "skip")),
        "accounts": {
            n: {"budget": a.budget, "spent": a.spent, "held": a.held}
            for n, a in run.accounts.items()
        },
        "ledger_total": sum(r["forwards"] for r in tally(run.ledger_path).values()),
        "unreconciled": [
            e for e in events if e["event"] == "reconcile" and not e["reconciled"]
        ],
    }


def run(run_dir, backends):
    """Run (or resume) a frozen D2 run to completion; write and return report.json."""
    state = open_run(run_dir)
    _l1(state, backends)
    _l2(state, backends)
    report = _report(state)
    (state.dir / "report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report
```

（`feasible_prefix` 仍供 `init_run` 使用。本 Task 的 `run()` 还没有 Task 5 的升级相。）

- [ ] **Step 4: 跑测试确认通过**

Run: `make -C explore check`
Expected: 全部 PASS

---

### Task 5: D2 升级——释放、认证、抽查与报告

**Files:**
- Modify: `explore/src/explore/d2.py`
- Test: `explore/tests/test_d2.py`

**Interfaces:**
- Consumes: Task 1 `view_factory`；`release`、`bracket_jump`（`probe.py`）；`certify`、`Certification`（`certify.py`）。
- Produces:
  - `_classify(run, row) -> {"class": "ridge_suspected" | "branch_bracketed" | "audit" | None, "info": dict | None, "reason": str | None}`
  - 升级顺序：类别（脊 0、夹住的跳支 1、抽查 2），同类内按 `stream(seed, "promote", intent_id)` 排序。
  - L3（搜索额，`release_cap`）与 L4（认证额，`certify_cap`）的 settle payload 都带 `class`、`audit_release`、`deferred_reason`。
  - 报告新增 `promotion`、`deferred`、`release_status`、`certify_status`、`c2`（每项：`task, cumulative_spent, saddle_energy_eV, chemical_pair, event, class, audit_release`）、`c1`。

- [ ] **Step 1: 写失败的测试**（追加到 `test_d2.py`）

```python
from twosaddle import SADDLE_Y, local_coords


def l4(run_dir):
    return [e for e in journal(run_dir) if e["event"] == "settle" and e["task"].startswith("L4|")]


def test_ridges_are_released_and_certified_on_two_saddles(tmp_path):
    init(tmp_path / "r")
    rep = run(tmp_path / "r", backends(Tick()))
    assert rep["promotion"].get("ridge_suspected", 0) >= 1
    assert rep["deferred"].get("no_trigger", 0) >= 1  # O0-O1 stretch
    assert len(rep["c2"]) >= 1
    for e in l4(tmp_path / "r"):
        assert e["payload"]["class"] == "ridge_suspected"
        assert e["status"] in ("c2", "duplicate_c2")
    seen = []
    for e in l4(tmp_path / "r"):
        if e["status"] == "c2":
            _, y, _ = local_coords(e["payload"]["saddle_positions"])
            k = min(range(2), key=lambda i: abs(SADDLE_Y[i] - y))
            assert abs(SADDLE_Y[k] - y) < 1e-3
            assert k not in seen
            seen.append(k)
    c1 = rep["release_status"].get("c1_candidate", 0)
    assert c1 == rep["promotion"]["ridge_suspected"]


def test_audit_items_are_flagged_and_paid_by_d2(tmp_path):
    init(tmp_path / "r", specs=(SPECS[2],), audit_fraction=1.0)
    tick = Tick()
    rep = run(tmp_path / "r", backends(tick))
    assert rep["promotion"] == {"audit": 1}
    done = [e for e in journal(tmp_path / "r") if e["event"] == "settle"]
    tags = [e["payload"]["audit_release"] for e in done if e["task"][:2] in ("L3", "L4")]
    assert tags == [True, True]
    paid = rep["accounts"]["search"]["spent"] + rep["accounts"]["certify"]["spent"]
    assert tick.count == paid


def test_search_exhaustion_leaves_the_certification_reserve(tmp_path):
    init(tmp_path / "r", search_budget=1800, release_cap=1800)
    rep = run(tmp_path / "r", backends(Tick()))
    assert rep["skips"].get("L3_skipped:search_budget", 0) == rep["promotion"]["ridge_suspected"]
    assert len(l4(tmp_path / "r")) == rep["promotion"]["ridge_suspected"]


def test_certification_budget_exhaustion_is_a_skip_not_a_loss(tmp_path):
    init(tmp_path / "r", certify_budget=100)
    rep = run(tmp_path / "r", backends(Tick()))
    assert l4(tmp_path / "r") == []
    assert rep["skips"]["L4_skipped:certify_budget"] == rep["promotion"]["ridge_suspected"]
    assert rep["accounts"]["certify"]["spent"] == 0
    assert rep["release_status"].get("c1_candidate", 0) >= 1


def test_resume_after_completion_dispatches_nothing(tmp_path):
    init(tmp_path / "r")
    first = run(tmp_path / "r", backends(Tick()))
    n = len(journal(tmp_path / "r"))
    tick = Tick()
    second = run(tmp_path / "r", backends(tick))
    assert tick.count == 0 and len(journal(tmp_path / "r")) == n
    assert second["c2"] == first["c2"]
```

**注意（写给实现者）：** 第一个测试假设两条 H 意图的约束扫描都在 0.1 Å 格上出现 `ridge_suspected`，且从粗扫峰点出发 `follow_min_mode` 能精修到双鞍点之一。这是按势的几何推出来的（两条路径都从左极小 (−1, 0) 越到右极小 (1, 0)），没有运行过。实测不符合就停下来，把扫描能量曲线和认证状态报回来，**不要改断言去凑结果**。

- [ ] **Step 2: 跑测试确认失败**

Run: `make -C explore test`
Expected: FAIL，`KeyError: 'promotion'`

- [ ] **Step 3: 实现**

import 块加上：

```python
from ase.calculators.singlepoint import SinglePointCalculator

from .certify import Certification, certify
from .ledger import Ledger, tally, view_factory
from .probe import (
    Intent,
    ScanPoint,
    ScanSettings,
    bracket_jump,
    feasible_prefix,
    read_scan,
    release,
    reverse_plan,
)
```

（替换原来的 `from .ledger import …` 和 `from .probe import …` 两行。）

`_Run` 加一个方法：

```python
    def known_c2(self):
        """Certified saddles so far, in journal order, for certify's duplicate rule."""
        out = []
        for task, ev in self.settled.items():
            if task.startswith("L4|") and ev["status"] == "c2":
                p = ev["payload"]
                cert = Certification(
                    status="c2",
                    saddle_positions=np.asarray(p["saddle_positions"], float),
                    saddle_energy_eV=p["saddle_energy_eV"],
                    chemical_pair=tuple(p["chemical_pair"]),
                )
                out.append((task, cert))
        return out
```

在 `_l2` 之后追加：

```python
_CLASS_RANK = {"ridge_suspected": 0, "branch_bracketed": 1, "audit": 2}


def _classify(run, row):
    """Promotion class or deferral reason of one scanned intent (pure, from the journal).
    Audit eligibility = deferred with an L1 response (at least two points), drawn
    once at freeze time (``row["audit"]``)."""
    ev = run.settled[_task("L1", row)]
    points = _points(ev["payload"])
    if len(points) < 2:
        return {"class": None, "info": None, "reason": f"no_response:{ev['status']}"}
    triggers = read_scan(points, run.settings.noise_eV)
    ridge = _first(triggers, "ridge_suspected")
    if ridge is not None:
        return {"class": "ridge_suspected", "info": ridge, "reason": None}
    jump = _first(triggers, "jump")
    reason = "no_trigger"
    if jump is not None:
        rev = run.settled.get(_task("L2", row))
        if rev is None:
            reason = "jump_unchecked"
        else:
            back = _points(rev["payload"]) if rev["status"] == "completed" else []
            bracket = bracket_jump(points, jump, back[-1]) if len(back) >= 2 else None
            if bracket is not None:
                return {"class": "branch_bracketed", "info": bracket, "reason": None}
            reason = "jump_not_bracketed"
    if row["audit"]:
        return {"class": "audit", "info": None, "reason": reason}
    return {"class": None, "info": None, "reason": reason}


def _promotion_points(classified, points):
    """(release positions, certification guess) by the frozen rule M5-D1."""
    info = classified["info"]
    if classified["class"] == "ridge_suspected":
        return points[info["release"]].positions, points[info["peak"]].positions
    if classified["class"] == "branch_bracketed":
        return (
            np.asarray(info["product_side_positions"], float),
            points[info["reactant_side"]].positions,
        )
    good = [k for k, p in enumerate(points) if p.converged] or list(range(len(points)))
    interior = [k for k in good if k >= 1] or [len(points) - 1]
    peak = max(interior, key=lambda k: points[k].energy_eV)
    return points[good[-1]].positions, points[peak].positions


def _source(template, point):
    atoms = template.copy()
    atoms.positions = point.positions
    atoms.calc = SinglePointCalculator(atoms, energy=point.energy_eV, forces=point.forces_pes)
    return atoms


def _release_json(out):
    end = out.get("endpoint")
    data = {k: v for k, v in out.items() if k != "endpoint"}
    data["endpoint_positions"] = None if end is None else end.positions.tolist()
    data["endpoint_energy_eV"] = None if end is None else float(end.get_potential_energy())
    return _jsonable(data)


def _release_task(run, backends, row, start, positions, tags):
    s = run.settings
    task = _task("L3", row)
    if not run.pending(task):
        return
    attempt = run.dispatch(task, "search", s.release_cap, "L3")
    if attempt is None:
        run.skip(task, "L3_skipped:search_budget")
        return
    ledger = run.task_ledger(s.release_cap)
    template = run.poses[row["pose_id"]]
    try:
        out = release(
            positions,
            template,
            _source(template, start),
            view_factory(backends.certify_calculator, ledger),
            run.config,
            ledger,
            task,
            seed=_seed(s, task),
            attempt=attempt,
        )
        status, payload = out["status"], _release_json(out)
    except Exception as exc:
        status, payload = f"error:{type(exc).__name__}", {"error": str(exc)}
    run.settle(task, attempt, "search", s.release_cap, ledger, status, {**payload, **tags})


def _certify_task(run, backends, row, start, positions, tags):
    s = run.settings
    task = _task("L4", row)
    if not run.pending(task):
        return
    attempt = run.dispatch(task, "certify", s.certify_cap, "L4")
    if attempt is None:
        run.skip(task, "L4_skipped:certify_budget")
        return
    ledger = run.task_ledger(s.certify_cap)
    guess = run.poses[row["pose_id"]].copy()
    guess.positions = np.asarray(positions, float)
    try:
        cert = certify(
            guess,
            view_factory(backends.certify_calculator, ledger),
            run.config,
            ledger,
            task,
            seed=_seed(s, task),
            ceiling_eV=start.energy_eV + s.certify_max_rise_eV,
            reference_eV=start.energy_eV,
            known=run.known_c2(),
            attempt=attempt,
        )
        status, payload = cert.status, _jsonable(vars(cert))
    except Exception as exc:
        status, payload = f"error:{type(exc).__name__}", {"error": str(exc)}
    run.settle(task, attempt, "certify", s.certify_cap, ledger, status, {**payload, **tags})


def _promote(run, backends):
    s = run.settings
    chosen = []
    for row in run.population:
        if _task("L1", row) in run.settled:
            classified = _classify(run, row)
            if classified["class"] is not None:
                chosen.append((row, classified))
    chosen.sort(
        key=lambda rc: (
            _CLASS_RANK[rc[1]["class"]],
            _key(s.seed, "promote", rc[0]["intent_id"]),
        )
    )
    for row, classified in chosen:
        points = _points(run.settled[_task("L1", row)]["payload"])
        release_at, guess_at = _promotion_points(classified, points)
        tags = {
            "class": classified["class"],
            "audit_release": classified["class"] == "audit",
            "deferred_reason": classified["reason"],
        }
        _release_task(run, backends, row, points[0], release_at, tags)
        _certify_task(run, backends, row, points[0], guess_at, tags)
```

`_report` 末尾的 return 字典里加入下列键（在 `"l2_status"` 之后）：

```python
        "promotion": dict(promotion),
        "deferred": dict(deferred),
        "release_status": _counts(run, "L3|"),
        "certify_status": _counts(run, "L4|"),
        "c1": [
            {"task": t, "label": ev["payload"].get("label"), "class": ev["payload"]["class"],
             "audit_release": ev["payload"]["audit_release"]}
            for t, ev in run.settled.items()
            if t.startswith("L3|") and ev["status"] == "c1_candidate"
        ],
        "c2": [
            {
                "task": t,
                "cumulative_spent": ev["cumulative_spent"],
                "saddle_energy_eV": ev["payload"]["saddle_energy_eV"],
                "chemical_pair": ev["payload"]["chemical_pair"],
                "event": ev["payload"]["event"],
                "class": ev["payload"]["class"],
                "audit_release": ev["payload"]["audit_release"],
            }
            for t, ev in run.settled.items()
            if t.startswith("L4|") and ev["status"] == "c2"
        ],
```

并在 `_report` 开头（`rows = run.population` 之后）算出：

```python
    promotion, deferred = Counter(), Counter()
    for row in rows:
        if _task("L1", row) in run.settled:
            classified = _classify(run, row)
            if classified["class"] is not None:
                promotion[classified["class"]] += 1
            else:
                deferred[classified["reason"]] += 1
```

`run()` 改为：

```python
def run(run_dir, backends):
    """Run (or resume) a frozen D2 run to completion; write and return report.json.
    Phases are strict (design §4.2): every sampled intent's L1 scan is settled or
    skipped before any L2/L3/L4 work starts."""
    state = open_run(run_dir)
    _l1(state, backends)
    _l2(state, backends)
    _promote(state, backends)
    report = _report(state)
    (state.dir / "report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report
```

- [ ] **Step 4: 跑测试确认通过**

Run: `make -C explore check`
Expected: 全部 PASS

---

### Task 6: kill→resume 出口验收（解析台，子进程真 kill）

**Files:**
- Create: `explore/tests/d2_driver.py`
- Test: `explore/tests/test_d2_resume.py`

**Interfaces:**
- Consumes: `d2_bench.init`、`d2_bench.backends`、`d2_bench.Tick`；`explore.d2.run`。
- Produces: M5 出口（上位计划 M5 行）：解析台上 kill→resume 不重复计费；抽查项带 `audit_release`（Task 5 已固定）；外部计数器与分阶段总账一致。

- [ ] **Step 1: 子进程入口** `explore/tests/d2_driver.py`

```python
"""python d2_driver.py RUN_DIR COUNTER CRASH_AFTER: run or resume the TwoSaddle D2
bench. Each model evaluation appends one byte to COUNTER first; evaluation number
CRASH_AFTER of this process calls os._exit(137) (no finally blocks, like SIGKILL).
CRASH_AFTER < 0 never crashes."""

import sys
from pathlib import Path

from d2_bench import Tick, backends, init
from explore.d2 import run


def main():
    run_dir, counter, crash_after = Path(sys.argv[1]), Path(sys.argv[2]), int(sys.argv[3])
    if not (run_dir / "manifest.json").exists():
        init(run_dir)
    run(run_dir, backends(Tick(counter, crash_after)))


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: 写测试** `explore/tests/test_d2_resume.py`

```python
import json
import os
import subprocess
import sys
from collections import Counter
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def drive(run_dir, counter, crash_after):
    path = os.pathsep.join(str(ROOT / p) for p in ("src", "../src", "tests"))
    return subprocess.run(
        [sys.executable, str(ROOT / "tests" / "d2_driver.py"), str(run_dir), str(counter),
         str(crash_after)],
        env={**os.environ, "PYTHONPATH": path},
        capture_output=True,
        text=True,
        timeout=900,
    )


def events(run_dir):
    text = (run_dir / "journal.jsonl").read_text(encoding="utf-8")
    return [json.loads(line) for line in text.splitlines() if line.strip()]


def report(run_dir):
    return json.loads((run_dir / "report.json").read_text(encoding="utf-8"))


def paid(rep):
    return rep["accounts"]["search"]["spent"] + rep["accounts"]["certify"]["spent"]


@pytest.fixture(scope="module")
def reference(tmp_path_factory):
    base = tmp_path_factory.mktemp("ref")
    done = drive(base / "run", base / "count", -1)
    assert done.returncode == 0, done.stderr
    return base / "run", (base / "count").stat().st_size


def test_uninterrupted_run_external_count_matches_every_ledger(reference):
    run_dir, external = reference
    rep = report(run_dir)
    assert external == rep["ledger_total"] == paid(rep)
    assert rep["unreconciled"] == []


def crash_point(ref_dir, phase):
    if phase == "L1":
        return 5  # scans of the first wave are in flight
    first = next(e for e in events(ref_dir) if e["event"] == "settle" and e["task"].startswith("L4|"))
    return first["cumulative_spent"] - first["spent"] // 2  # inside the first certification


@pytest.mark.parametrize("phase", ["L1", "L4"])
def test_kill_and_resume_pays_once_and_finds_the_same_saddles(reference, tmp_path, phase):
    ref_dir, ref_external = reference
    run_dir, counter = tmp_path / "run", tmp_path / "count"
    first = drive(run_dir, counter, crash_point(ref_dir, phase))
    assert first.returncode == 137, first.stderr
    cut = events(run_dir)
    settled = {(e["task"], e["attempt"]) for e in cut if e["event"] == "settle"}
    open_kinds = {
        e["kind"] for e in cut
        if e["event"] == "dispatch" and (e["task"], e["attempt"]) not in settled
    }
    assert phase in open_kinds

    second = drive(run_dir, counter, -1)
    assert second.returncode == 0, second.stderr
    log, rep, ref = events(run_dir), report(run_dir), report(ref_dir)

    # same science as the uninterrupted run
    def energies(r):
        return sorted(round(c["saddle_energy_eV"], 8) for c in r["c2"])

    assert energies(rep) == energies(ref)
    assert rep["release_status"] == ref["release_status"]
    assert rep["certify_status"] == ref["certify_status"]

    # every task settles once; an interrupted attempt is retried under a larger number
    assert max(Counter(e["task"] for e in log if e["event"] == "settle").values()) == 1
    reconciles = [e for e in log if e["event"] == "reconcile"]
    assert reconciles
    for r in reconciles:
        assert any(
            e["event"] == "settle" and e["task"] == r["task"] and e["attempt"] > r["attempt"]
            for e in log
        )

    # cost: nothing refunded, nothing paid twice
    wasted = sum(r["spent"] for r in reconciles)
    assert paid(rep) == paid(ref) + wasted
    external = counter.stat().st_size
    assert external <= paid(rep)
    assert external - ref_external <= wasted
```

- [ ] **Step 3: 跑测试**

Run: `make -C explore check`
Expected: 全部 PASS。若 `paid(rep) == paid(ref) + wasted` 不成立，说明重做的任务花费和参考不同（不确定性来源），停下来报告，不放宽断言。

- [ ] **Step 4: 写里程碑提交草稿** `.commit_drafts/<日期_时间>_explore-m5.md`（commit message、文件清单逐个一句话、`git add … && git commit -m "…"`；Task 7–9 完成后补上 benchmarks 与 docs 的文件）。

---

### Task 7: b04 float32 扫描噪声与 fp32/fp64 同进程共存（183 GPU；先报计划，等“跑”）

**Files:**
- Create: `benchmarks/explore_d2/fp32_noise.py`
- Output: `benchmarks/explore_d2/results/<日期>/fp32_noise/`（在 183 上读）

**为什么做：** 上位计划要求 M5 先在 float32 下重测扫描噪声，再定 `read_scan` 的 `noise_eV`（M5-D5）。同时要确认一个前提：10-22 的脚本按“MACE 会改 torch 全局默认 dtype”一进程只放一种精度；而 M5 要在一个进程里同时常驻 fp32（搜索）和 fp64（释放、认证）。我只读核对过 `openmm_dev` 里的 `mace/calculators/mace.py`：`MACECalculator` 构造时把模型 `.float()` / `.double()`，构图时用 `torch_tools.default_dtype(...)` 上下文管理器，没有调用全局 `set_default_dtype`。这只是读代码，**【未实测】**，所以要在 GPU 上验一次。

**做什么：**
1. (a) 正向扫描：b04 输入、意图 C0–Cl5 压缩 3.6→1.9 Å、0.1 Å 步长、fmax 0.01、每点 300 步、max_rise 3.0（10-13 协议），用 fp32 `MACEBatchEvaluator` + `scan_many`，打开曲率复用。
2. (b) 噪声：取 q = 3.2、2.7、2.4、2.3 Å 四个点，每点 5 次重复。每次给每个笛卡尔分量加 N(0, 0.02 Å) 扰动，再拉回该点的目标 q 并弛豫（同上设置）。20 个重复任务作为一批 `scan_many` 跑。每点报告收敛能量的极差、标准差，以及与基点能量的差。
3. (c) 共存：同一进程里在 fp32 之后再建 fp64 POLAR-1，对 3 个几何（b04 输入、`ts0000.extxyz`、扫描峰点）分别用 fp32（建 fp64 之前和之后各一次）和 fp64 求值。另起一个只有 fp64 的进程（`--mode fp64-only`）对同样 3 个几何求值，作为对照。报告：同进程 fp64 与单精度进程 fp64 的能量、力之差；fp32 在建 fp64 前后之差。不设通过门，数值交用户判断。

- [ ] **Step 1: 写脚本**

```python
"""b04 float32 scan noise and float32/float64 coexistence (PLAN_EXPLORE_M5 Task 7).

  python fp32_noise.py --out DIR                 # (a) scan, (b) noise, (c) coexistence
  python fp32_noise.py --out DIR --mode fp64-only

Smoke on one GPU. Scan and noise forwards go through ledger.jsonl; the coexistence
evaluations (3 geometries x 4 readings) are plain calls, listed in coexist.npz. Not a
gate: numbers for the user to freeze noise_eV (M5-D5) and to judge one-process mixed
precision."""

import argparse
import json
import socket
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT / "src"), str(ROOT / "explore" / "src")]

from ase.io import read  # noqa: E402

from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.ledger import Ledger  # noqa: E402
from explore.parallel import MACEBatchEvaluator, ScanJob, scan_many  # noqa: E402
from explore.probe import Intent, ScanSettings  # noqa: E402
from prrs.config import SearchConfig  # noqa: E402
from prrs.search import _stamped_with_charge  # noqa: E402

INPUT = ROOT / "benchmarks/b04_sn2_chloride/input.extxyz"
TS = ROOT / "runs/b04_probe2_openceiling/ts_candidates/ts0000.extxyz"
INTENT = Intent("C0-Cl5", ((0, 5),), (1.0,))
NOISE_Q = (3.2, 2.7, 2.4, 2.3)
REPEATS = 5
SIGMA_A = 0.02
SEED = 20261023
FMAX, STEPS, RISE = 0.01, 300, 3.0


def config():
    return SearchConfig(charge_sensitive=True, total_charge=-1, multiplicity=1)


def source():
    atoms = _stamped_with_charge(read(INPUT), config())
    atoms.info["spin"] = 1
    return atoms


def settings(deltas):
    return ScanSettings(tuple(deltas), FMAX, STEPS, RISE, reuse_curvature=True)


def evaluate(calc, geometries):
    out = []
    for atoms in geometries:
        a = atoms.copy()
        a.calc = calc
        out.append((float(a.get_potential_energy()), a.get_forces().copy()))
    return out


def geometries(peak_positions=None):
    src = source()
    ts = read(TS)
    ts.info.update(src.info)
    geos = [src, ts]
    if peak_positions is not None:
        peak = src.copy()
        peak.positions = peak_positions
        geos.append(peak)
    return geos


def run_main(out):
    rec = {"host": socket.gethostname(), "seed": SEED, "sigma_A": SIGMA_A}
    t0 = time.perf_counter()
    calc32 = polar1("cuda", "float32")
    rec["load32_s"] = time.perf_counter() - t0
    batch = MACEBatchEvaluator(calc32)
    cfg = config()
    src = source()
    q0 = INTENT.value(src.positions)
    deltas = [0.0] + [round(-0.1 * k, 10) for k in range(1, int(round((q0 - 1.9) / 0.1)) + 1)]
    ledger_path = out / "ledger.jsonl"
    scan_job = ScanJob(src, INTENT, settings(deltas), cfg,
                       Ledger(ledger_path, "fp32-scan", replay=False), "scan")
    t0 = time.perf_counter()
    (scan,), timing = scan_many([scan_job], batch, 1)
    rec["scan_wall_s"] = time.perf_counter() - t0
    if scan.error is not None:
        rec["scan_error"] = scan.error
        (out / "noise.json").write_text(json.dumps(rec, indent=2), encoding="utf-8")
        return
    points = scan.result.points
    rec["scan"] = {
        "status": scan.result.status,
        "forwards": scan_job.ledger.spent,
        "q": [p.q_actual for p in points],
        "energy_eV": [p.energy_eV for p in points],
        "converged": [p.converged for p in points],
    }
    rng = np.random.default_rng(SEED)
    jobs, keys = [], []
    for q in NOISE_Q:
        base = min(points, key=lambda p: abs(p.q_target - q))
        for r in range(REPEATS):
            a = src.copy()
            a.positions = base.positions + rng.normal(0.0, SIGMA_A, base.positions.shape)
            back = settings((0.0, base.q_target - INTENT.value(a.positions)))
            jobs.append(ScanJob(a, INTENT, back, cfg,
                                Ledger(ledger_path, "fp32-noise", replay=False), f"noise|{q}|{r}"))
            keys.append((q, base))
    t0 = time.perf_counter()
    outcomes, _ = scan_many(jobs, batch, 4)
    rec["noise_wall_s"] = time.perf_counter() - t0
    noise = {}
    for (q, base), job, o in zip(keys, jobs, outcomes):
        row = noise.setdefault(str(q), {"q_target": base.q_target, "base_eV": base.energy_eV,
                                        "energies_eV": [], "failed": [], "forwards": 0})
        row["forwards"] += job.ledger.spent
        last = None if o.result is None or not o.result.points else o.result.points[-1]
        if last is not None and last.converged:
            row["energies_eV"].append(last.energy_eV)
        else:
            row["failed"].append(o.error or (o.result.status if o.result else None))
    for row in noise.values():
        e = np.asarray(row["energies_eV"])
        row["ptp_eV"] = float(np.ptp(e)) if len(e) >= 2 else None
        row["std_eV"] = float(e.std(ddof=1)) if len(e) >= 2 else None
        row["max_dev_from_base_eV"] = float(np.max(np.abs(e - row["base_eV"]))) if len(e) else None
    rec["noise"] = noise
    peak = max(points, key=lambda p: p.energy_eV)
    geos = geometries(peak.positions)
    before = evaluate(calc32, geos)
    t0 = time.perf_counter()
    calc64 = polar1("cuda", "float64")
    rec["load64_s"] = time.perf_counter() - t0
    after32 = evaluate(calc32, geos)
    same64 = evaluate(calc64, geos)
    np.savez(out / "coexist.npz", peak=peak.positions,
             e32_before=[v[0] for v in before], e32_after=[v[0] for v in after32],
             e64=[v[0] for v in same64],
             f32_before=np.stack([v[1] for v in before]),
             f32_after=np.stack([v[1] for v in after32]),
             f64=np.stack([v[1] for v in same64]))
    (out / "noise.json").write_text(json.dumps(rec, indent=2), encoding="utf-8")


def run_fp64_only(out):
    peak = np.load(out / "coexist.npz")["peak"]
    calc64 = polar1("cuda", "float64")
    values = evaluate(calc64, geometries(peak))
    np.savez(out / "fp64_only.npz", e64=[v[0] for v in values],
             f64=np.stack([v[1] for v in values]))
    a, b = np.load(out / "coexist.npz"), np.load(out / "fp64_only.npz")
    compare = {
        "fp64_same_process_minus_alone_eV": (a["e64"] - b["e64"]).tolist(),
        "fp64_force_maxdiff_eV_A": float(np.abs(a["f64"] - b["f64"]).max()),
        "fp32_after_minus_before_eV": (a["e32_after"] - a["e32_before"]).tolist(),
        "fp32_force_maxdiff_after_before_eV_A": float(np.abs(a["f32_after"] - a["f32_before"]).max()),
        "fp32_minus_fp64_eV": (a["e32_after"] - a["e64"]).tolist(),
    }
    (out / "coexist.json").write_text(json.dumps(compare, indent=2), encoding="utf-8")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--mode", choices=("main", "fp64-only"), default="main")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (run_main if args.mode == "main" else run_fp64_only)(out)


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: 静态检查**（本机，不跑模型）

Run: `/home/ruigengji/MLPATH/.venv-tools/bin/ruff check benchmarks/explore_d2/fp32_noise.py && /home/ruigengji/MLPATH/.venv-tools/bin/ruff format --check benchmarks/explore_d2/fp32_noise.py`
Expected: 通过

- [ ] **Step 3: 报计划表，等用户说“跑”**（模板见下方“GPU 作业计划”），然后按全局规则的标准写法在 183 上 `setsid nohup` 启动 `main`，跑完再启动 `fp64-only`。PID 记进结果目录的 `pids.txt`；进度和结果在 183 上读。

- [ ] **Step 4: 结果写进 DONE_LOG，`noise_eV` 交用户定（M5-D5）。**

---

### Task 8: b04 D2 端到端冒烟（183 GPU；先报计划，等“跑”）

**Files:**
- Create: `benchmarks/explore_d2/run_b04_d2.py`
- Output: `benchmarks/explore_d2/results/<日期>/b04_d2/`（run 目录本身即结果；在 183 上读）

**做什么：** 一个进程、一张 GPU：建 b04 L0（CH₃Cl 构象经资格确认，其费用记入单独的共享准备账本 `prepare.jsonl`；Cl⁻ 单原子不做资格确认）→ 2 个姿态 → 每个姿态的距离单原语意图 → `init_run` → `run`。L1/L2 用常驻 fp32 批量后端，释放和认证用常驻 fp64。中途按 PID 杀一次再重跑同一条命令，在真实模型上看续跑。冒烟取值都写在脚本常量里，**不是 U1–U5 的冻结值**。

- [ ] **Step 1: 写脚本**

```python
"""b04 end-to-end D2 smoke (PLAN_EXPLORE_M5 Task 8).

  python run_b04_d2.py --run DIR --noise-eV X

One process, one GPU: a resident float32 batch model for L1/L2 and a resident float64
POLAR-1 for release and certification. Rerunning the same command resumes. Smoke
values below are NOT the frozen U1-U5 / budget choices."""

import argparse
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT / "src"), str(ROOT / "explore" / "src")]

from ase.io import read  # noqa: E402

from benchmarks.b04_sn2_chloride.calculator import polar1  # noqa: E402
from explore.d2 import Backends, D2Settings, init_run, run  # noqa: E402
from explore.intents import distance_intents  # noqa: E402
from explore.ledger import Ledger, view_factory  # noqa: E402
from explore.parallel import MACEBatchEvaluator  # noqa: E402
from explore.poses import (  # noqa: E402
    conformer_grid,
    conformers_from_qualification,
    fragments,
    qualify_conformers,
    sample_poses,
)
from prrs.config import SearchConfig  # noqa: E402
from prrs.search import _stamped_with_charge  # noqa: E402

INPUT = ROOT / "benchmarks/b04_sn2_chloride/input.extxyz"
SEED = 20261023
SMOKE = {"shells": [2.0], "n_samples": 2, "steps_per_turn": 6, "conformer_cap": 1000, "n_far": 2}


def settings(noise_eV):
    return D2Settings(
        seed=SEED, search_budget=8000, certify_budget=3000, scan_cap=400, reverse_cap=100,
        release_cap=300, certify_cap=800, per_pose_quota=100, audit_fraction=0.1,
        noise_eV=noise_eV, certify_max_rise_eV=1.0,
        magnitudes_A=tuple(round(0.1 * k, 10) for k in range(1, 18)),
        fmax_eV_A=0.01, max_steps_per_point=300, max_rise_eV=3.0, max_active=4, wave_size=16,
    )


def build_l0(atoms, config, calc64, run_dir):
    ledger = Ledger(run_dir.parent / f"{run_dir.name}_prepare.jsonl", arm="P", replay=False)
    frags = fragments(atoms, config)
    states = {f.key: {"charge": -1 if len(f.atoms) == 1 else 0, "multiplicity": 1} for f in frags}
    chosen = []
    for f in frags:
        confs, _ = conformer_grid(atoms, f, config, SMOKE["steps_per_turn"],
                                  SMOKE["conformer_cap"], SEED)
        if len(f.atoms) > 1:
            rows = qualify_conformers(atoms, confs, view_factory(calc64, ledger), config,
                                      ledger, SEED, states)
            confs = conformers_from_qualification(confs, rows)
        chosen.append(confs[0])
    ledger.close()
    poses = sample_poses(atoms, frags, chosen[0], chosen[1], SMOKE["shells"],
                         SMOKE["n_samples"], SEED)
    out = []
    for p in poses:
        if p.status != "placed":
            continue
        a = atoms.copy()
        a.positions = p.positions
        out.append((p.pose_id, a))
    return out, ledger.spent


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--noise-eV", type=float, required=True)
    args = ap.parse_args()
    run_dir = Path(args.run)
    config = SearchConfig(charge_sensitive=True, total_charge=-1, multiplicity=1)
    t0 = time.perf_counter()
    calc32 = polar1("cuda", "float32")
    calc64 = polar1("cuda", "float64")
    load_s = time.perf_counter() - t0
    if not (run_dir / "manifest.json").exists():
        atoms = _stamped_with_charge(read(INPUT), config)
        atoms.info["spin"] = 1
        poses, prepare = build_l0(atoms, config, calc64, run_dir)
        specs = {pid: distance_intents(a, pid, config, SMOKE["n_far"], SEED)[0] for pid, a in poses}
        identity = {"system": "b04", "model": "MACE-POLAR-1-M", "search_dtype": "float32",
                    "certify_dtype": "float64", "l0": SMOKE, "prepare_forwards": prepare}
        init_run(run_dir, poses, specs, settings(args.noise_eV), config, identity)
    report = run(run_dir, Backends(MACEBatchEvaluator(calc32), calc64))
    report["load_s_this_process"] = load_s
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: 静态检查**（同 Task 7 Step 2）

- [ ] **Step 3: 报计划表，等用户说“跑”。** 用户同意后：启动 → 在第一个 L4 dispatch 出现后（在 183 上 `tail` journal 看到）按 PID 杀掉 → 用同一条命令重跑到完成。

- [ ] **Step 4: 结果写进 DONE_LOG**：forward 分账（P / 搜索 / 认证）、各阶段墙钟、升级与暂缓分布、C1/C2、重复、unreconciled、续跑前后对账。C2 鞍点与 `ts0000` 对齐 RMSD 做事后诊断（只读）。

---

### Task 9: 文档

**Files:**
- Modify: `docs/PLAN_EXPLORE_V1_2026-10-13.md`（状态行）、`docs/DONE_LOG.md`、`docs/README.md`（10-13 那一行）、`explore/README.md`（模块列表）
- Modify: `.commit_drafts/<日期_时间>_explore-m5.md`（补全文件清单）

- [ ] **Step 1:** 上位计划状态行写明 M5 已施工，以及哪些已验证（解析台 kill→resume）、哪些没有（多通道收益、MD、正式预算）。
- [ ] **Step 2:** DONE_LOG 新增一节：本计划的偏离、M5-D1…D6 的实际裁决、测试数（亲自跑 `make -C explore check` 的输出）、Task 7/8 的数字及读取位置。
- [ ] **Step 3:** 回复里写“文档已更新：…”。

---

## GPU 作业计划（Task 7、8；报给用户时把 `nvidia-smi` / `ps` 的现场输出填进去）

| 项 | Task 7 噪声 + 共存 | Task 8 b04 D2 冒烟 |
| --- | --- | --- |
| 命令 | `python benchmarks/explore_d2/fp32_noise.py --out benchmarks/explore_d2/results/<日期>/fp32_noise`，跑完再加 `--mode fp64-only` | `python benchmarks/explore_d2/run_b04_d2.py --run benchmarks/explore_d2/results/<日期>/b04_d2 --noise-eV <M5-D5 的值>` |
| 节点 | 183 一张 2080 Ti（以 `nvidia-smi` 为准，卡被占就不抢） | 同左 |
| CPU 线程 | `OMP_NUM_THREADS=MKL_NUM_THREADS=4` | 同左 |
| 预计耗时及依据 | 加载 30–105 s/次（10-22、10-23 实测），共 3 次；扫描约 84 forward × 0.1 s；噪声 20 个短任务，4 路合批（10-23 实测 4 路 323 forward / 16.6 s）。合计约 5–8 分钟 | 每姿态约 11 条意图（10-22 计数表），2 个姿态约 22 条扫描 × 约 84–158 forward；4 路合批约 0.05 s/forward → 约 2–3 分钟；释放约 70 forward、认证约 200 forward/个，fp64 约 0.16 s/forward（10-22 实测）。升级个数事先不知道：按 2–6 个估，约 1–4 分钟。加 2 次加载，再加一次 kill→resume。合计约 8–15 分钟 |
| 能否续跑 | 不可续跑（短测量，失败就整次重跑） | 可续跑（重跑同一条命令即可） |
| 输出目录 | `benchmarks/explore_d2/results/<日期>/fp32_noise/` | `benchmarks/explore_d2/results/<日期>/b04_d2/` |

## 自查记录

- 覆盖设计稿：两账户（Task 2、5）；固定分级与严格分相（Task 4、5，M5-D6）；均匀抽查 + `audit_release`（Task 3 抽签、Task 5 执行）；kill→resume 不回滚、不重复计费（Task 2、4、6）；失败码不压缩（scan 状态、`error:*`、各 skip 原因、L3/L4 原样保留状态）；首次命中成本 = C2 结算时的 `cumulative_spent`（Task 5 报告）。设计 §1.1 中 L3 MD 与 §6.1 第 3 条明确不做，原因见“不同的地方”。
- 占位扫描：没有 TBD。Task 5、Task 6 中两处依赖未运行过的结果（双鞍点势上出脊、重做花费确定），都写了“不符就停下来报告，不改断言”。
- 类型一致：任务 ID 前缀 `L1|…L4|`、journal 字段、`D2Settings` 字段在 Task 3–8 中一致；`view_factory(inner, ledger)` 在 Task 1、5、8 中一致；`reverse_plan(points, jump, settings)` 在 Task 3、4 中一致。
- Review Focus 五条：1 → Task 6；2 → Task 4（`test_interrupted_attempts_are_charged_and_retried`）+ Task 6；3 → Task 5（`test_search_exhaustion_leaves_the_certification_reserve`）；4 → Task 1、Task 2；5 → Task 3。
