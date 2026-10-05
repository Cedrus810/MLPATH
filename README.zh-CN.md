# PRRS — 扰动–响应反应搜索（Perturbation–Response Reaction Search）

[![License: AGPL v3](https://img.shields.io/badge/License-AGPL_v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)
![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue)
![Lint: ruff](https://img.shields.io/badge/lint-ruff-261230)

[English](README.md) | 简体中文

**MLPATH 项目入口**：[整体逻辑与当前进度](docs/PROJECT_OVERVIEW.md) · [文档导航](docs/README.md) · [PARS 路径模块](pars/README.md)。
MLPATH 包含 PRRS 反应发现、PARS 路径选择与精修，以及规划中的射击成熟层。本文介绍 PRRS；整体关系与 connect 进度见项目总览。

**PRRS 是一个面向机器学习势（MLP）的反应发现方法。**

> 给分子一个受控扰动 → 撤去扰动，让它在 MLP 势能面上自由响应 → 淬火 → 看它落到哪个化学态。

不预设反应坐标，不预设产物，不从过渡态猜起。先发现连通性（可能通道与产物），
之后才做过渡态的定位与验证。

---

## 机理

```
源状态 → 受控扰动（内坐标） → 短 MLP 轨迹 → 淬火 → 分类
                                                │
    源构象 ◄── chemical_key 相同 ──► 新化学态
                                                │
                   一阶鞍点（n₋ = 1）──► TS pool
```

每次试验是一个两段式的小数值实验：先是把扰动**增量交付**到内坐标上的**交付段**，
然后是撤去扰动后无偏松弛的**自由响应段**。三种可能结局：

1. **回到源极小** —— 弹性证据：在测试幅度下，这个方向不易通向别的盆地。
2. **松弛到另一个极小** —— 盆地发现。构象与新化学态被显式区分：裁决依据是化学
   连通性，不是能量或距离。
3. **恰好一个负 Hessian 本征值的驻点** —— TS 候选，进入 TS pool，随后沿两侧下降
   方向验证它连接哪两个态。

### 撑起这个方法的几个设计决策

- **扰动作用在 primitive internal coordinates 上**（键长、键角、二面角、集体
  fragment 运动、分子间接近坐标），而不是笛卡尔噪声。这是实测结论不是偏好：乙醇
  早期只用径向扰动的试验全部弹回，而三个真实构象盆地就在 4 meV 处 —— 成对径向
  位移在轨迹长度内耦合不进扭转。**扰动集合必须张开化学上相关的自由度。**
- **交付是增量的**（每步 ≤ 0.05 Å/rad，每步重算梯度），因为大扰动 = 一串局部
  有效的小扰动。请求的目标位移与附带损伤分开报告；共价键损伤与原子重叠是硬门。
- **对称性显式处理。** 化学等价原子绝不产生人为的"新反应"：身份用自同构不变的
  图不变量哈希，构象去重用对称感知 RMSD，立体化学签名只取构型手性。对称转子
  **降级但不删除** —— 其中一个贡献过真实的甲基旋转鞍点。
- **三层身份，永不互相折叠**：*化学态*（连通性 + 构型立体化学）、*构象水库*
  （一个态多个几何 —— 反应可及性依赖起始构象）、*反应事件*（键差 key 在两端
  自同构群下正则化）。退化反应（如等价 H 的质子转移）是**自环**，不是新节点。
- **收敛判据按模分级。** 固定力容差只控制"力有多小"，不控制"离驻点多远"：
  乙醇扭转比键伸缩软约 150 倍，软模有独立的 |g/k| 判据，真正自由的模只记录、
  不强行判"收敛"。
- **失败必须关闭（fail closed）。** 判不出来就拒绝并说明理由，绝不用默认值替代；
  存在性必须由执行证实，不能由记录证实。

---

## PRRS 负责——和不负责——什么

| 问题 | 比较对象 | 检验的是 |
|---|---|---|
| **搜索**好不好？ | PRRS vs **同一张 MLP 面**上的参照搜索 | 方法本身 |
| **势**准不准？ | MLP 过渡态 vs QM | 能量模型 |

PRRS 负责发现通道、识别产物、定位并验证 TS 候选、构建网络；**不**负责给出最终
的定量势垒：势垒高度与虚频是所选势能面的属性，交下游 QM 检验。把两种比较混在
一起，误差就无法归因 —— 苯甲酰丙酮三模型对照中，各模型的化学身份字段逐字段一致，
而势垒差 3.7 倍，正是这条边界要切分的东西。

---

## 当前状态

基准套件（`benchmarks/`）的规矩：一格 = 一个反应体系 + **跑前冻结**的判据 +
一个机械按判据评定的脚本。

| 阶段 | 体系 | 新打的轴 | 状态 |
|---|---|---|---|
| P0 | 乙醇 | 构象发现；框架 | 完成 —— 13/13 |
| P1 | 丙二醛 | 反应事件 —— **退化**自环 | 完成 —— 13/13，QM 复核通过 |
| P2 | 3-氧代丁醛烯醇 | 网络扩张 A→B | 完成 —— 16/16 seed |
| P3 | 苯甲酰丙酮 | 芳香性、螯合、竞争异构体 | 完成 —— 判据冻结并判定通过 |
| b04 | Cl⁻ + CH₃Cl | **总电荷**、2 碎片、重原子坐标 | 完成 —— 判据 2026-09-13 冻结，TS 已 QM 确认 |
| b05 | 叔卤代烷 + 亲核试剂 | 取代/消除竞争 | 进行中 |
| b06 | 亲电加成底物 | 极性加成 | 计划中 |
| b07 | 环己烯（逆 DA）· 丁二烯 + 乙烯（DA） | **碎片数 1↔2** | 计划中（探针已绿） |

P0–P3 已转为冻结的回归夹具：结果仍有效，正式验收职责移交 P4+ 的独立 QM 协议。
这条序列刻意从反复的质子转移走向带电体系、分子间反应、竞争通道、重原子
断键/成键、碎片数变化。

### 独立 QM 复核

协议：**r2SCAN-3c** 出几何与频率 → **ωB97M-V/def2-TZVPD** 在同一几何上补单点。
（r2SCAN-3c 单独用会把 SN2 中心势垒错约 11 kcal —— 几何好，能量必须补。）

- **丙二醛（P1）**：从 PRRS 鞍点出发的双向 IRC 收敛到预期极小，两端均 0 虚频。
  MLP 与 QM 鞍点几何差 **0.0151 Å**（对称感知 RMSD）。MLP 势垒 419 meV 对 QM
  134 meV 差得很大 —— 这是能量模型的误差，不是没找到反应几何。
- **SN2：Cl⁻ + CH₃Cl（b04）**：PRRS 找到的过渡态精化到 2.3645/2.3645 Å，与独立
  求得的 QM 鞍点小数点后四位相同；双向 IRC 到达预期反应物与产物；反应事件 key
  与 PRRS 逐字一致。

> **即使 MLP 不能准确复现 QM 势垒能量，PRRS 也能给出正确的反应事件与过渡态几何。**

### 机器学习势的选择

P0–P3 冻结在 **MACE-OFF24** 上，不重跑、不改标签。从 b04 起默认模型是
**MACE-POLAR-1-M**，依据是 QM 对照而不是偏好：它是唯一在"中性质子转移化学"和
"离子 SN2 基准"两边都成立的候选。MACE-OFF24 对总分子电荷完全不响应（电荷
0/−1/+1 能量逐位相同），离子化学因此不能用。冻结结果不会自动成为后续阶段的
参照。

---

## 安装

要求：**Python ≥ 3.10**；核心依赖只有 ASE 与 NumPy。

```bash
pip install -e .            # 核心：prrs 命令行；可在任何 ASE calculator 上运行
pip install -e ".[mace]"    # MACE 支持（mace-torch）
pip install -e ".[test]"    # pytest
```

冻结 P0–P3 夹具用的 MACE-OFF24 medium checkpoint 随仓库发布在 `models/`
（sha256 前缀 `e5ccf5837f685899`）。每个基准在 `benchmarks/b0*/model.json`
里声明自己的模型 —— **绝不隐式下载任何模型**。

### 两个环境，有意为之

**科学环境**（conda 规格在 `env/openmm_dev.yml`；基线与沿革在 `env/BASELINE.md`）
产出全部对外数字，其包版本是结果的一部分。开发工具（ruff）放在独立的
`.venv-tools`（`make tools`），让格式化永远不会推动科学依赖。backend parity 曾因
openmm 的一个 patch 版本变红 —— 绝不往科学环境里装开发工具。

## 快速上手

可选的局部极小／TS 候选精修入口为 `prrs optimize`，通过 geomeTRIC 支持 ASE calculator
与 OpenMM/MACE Context。命令、输出与验证边界见 [精修说明](docs/GEOMETRIC_REFINEMENT.md)。

```bash
# 最小自带演示
prrs demo --output runs/demo --seed 42

# 一次搜索运行
prrs search input.extxyz --output runs/my_run \
    --config config.json --calculator mace --calculator-kwargs calc.json

# 导出未精化的路径种子，供下游 NEB / dimer / string 精化；
# 淬火迭代被排除 —— 它们不是动力学
prrs path --run runs/my_run --attempt <id> --output seed.extxyz --images 9
```

每次运行写一个 manifest（含 `implementation_sha256`，见[开发](#开发)）和一个
`network.json` 反应网络。`runs/README.md` 记录磁盘上哪些运行计入证据 ——
`network.json` **存在**不等于那次运行跑完了。

### 实时观察一次运行

```bash
python scripts/prrs_watch.py runs/my_run
```

只读终端 UI：按字节偏移增量 tail 每条试验的 `observations.jsonl`、轮询
`network.json`，把演化的势能景观直接画在终端里（真彩色）。它**绝不写 `runs/`**；
全分辨率 PNG 写到 `.scratch/prrs_watch/<run>/`。

---

## 文件地图

```
src/prrs/                 搜索核心 —— 19 个模块，约 12k 行
  config.py               全部协议参数 + 校验；不合法配置在启动前拒绝
  internal.py             primitive internals、精确刚体 fragment 旋转、continuation 位移
  chemistry.py            三层身份、自同构群、对称感知 RMSD
  canonical.py            individualization–refinement 的可认证典范标记
  perturbations.py        六个扰动家族、对称商、分层方向预算
  network.py              Registry：接纳、构象水库、域门
  state.py                活性子集图编码、Kabsch RMSD、响应分类
  reliability.py          每次力评估的失败关闭门控
  runner.py               单次试验：交付、mode-resolved 淬火、Hessian/极小检验、鞍点下降
  search.py               预算调度、幅度细化、TS pool、run manifest
  completeness.py         盆地归属的条件完备性证书
  product.py              connect 模式的产物目标：只记录距离，绝不路由
  schema.py               显式的 network.json schema（schema 2）
  calculators.py          calculator 工厂
  openmm_backend.py       经实测验证的 openmm-ml/MACE 原语
  torch_guard.py          唯一一个经审计的 torch.load 包装
  validation.py           独立交叉校验路径（ase.vibrations）—— 搜索侧从不 import
tests/                    31 个文件、349 个测试函数
benchmarks/               基准套件 —— 每个反应体系一个目录（判据 + 评定脚本）
docs/                     工作记录：交接、验收文档、计划、实验日志
scripts/                  溯源、实现哈希、实时观察器
runs/                     实验输出 —— runs/README.md 说了哪些计入证据
models/                   MACE-OFF24_medium.model（冻结 P0–P3 夹具的模型）
env/                      科学环境基线：BASELINE.md、conda 规格、pip freeze
```

---

## 开发

```bash
make help      # 各 target 的作用
make check     # fmt-check + lint + test + digest，按此顺序 —— 完整验证门
make digest    # 黄金样本回归门
make hash      # run manifest 会记的 implementation_sha256
```

- **pytest 问的是"代码对不对"，digest 问的是另一个问题** —— 代码是否仍然产出与
  磁盘上已有运行相同的数字（332 文件 / 2136 帧，每个差异都被分类为
  `graph_invariant` / `parity_only` / `key_only`）。**digest 零差异是"零行为
  改动"的唯一被承认的证据。** 有意的改动用 `--allow parity_only` 之类放行 ——
  但放行只确认差异的*形状*，不确认新值是对的；那需要化学判断。
- **`implementation_sha256`**（`make hash`）是每个 run manifest 的溯源锚：对
  `src/prrs/*.py` 的内容哈希，不依赖 VCS 状态。绝不手抄；现行指纹的记录在
  `docs/HANDOFF.md` 与 `env/BASELINE.md` 维护。
- 最近一次记录的完整验证门：**476 passed, 1 skipped，digest 零差异**（2026-09-26）。

### 项目的规矩

这些是工作规矩，不是建议；违反它们产出的结果以前被作废过。

- **判据跑前冻结**，之后绝不原地修改；要改必须作为显式修订，记录旧文本原文、
  改的理由、以及修订发生在看到哪些数据之后。
- **干预分别命名**，改进才能归因到具体某一项。
- **未完成 ≠ 失败。** 分母是**发起数**；未完成与崩溃按原样记录并附原因
  （`scripts/scan_status.sh` 机械区分三态）。
- **失败必须关闭（fail closed）。** 判不出来就拒绝并说明理由，绝不用默认值替代。
- **门不会因为挡住了想要的结果而移动。** 先量清楚它挡掉的是什么。
- **测量优先于推断。** 被记录的证伪是资产，不是耻辱（`docs/HANDOFF.md` §6）。

---

## 文档

| | |
|---|---|
| [docs/README.md](docs/README.md) | 文档索引 —— 哪些是当前的、哪些只是记录 |
| [docs/HANDOFF.md](docs/HANDOFF.md) | 交接总览；冻结结论以此为准 |
| [docs/DECISIONS_PENDING.md](docs/DECISIONS_PENDING.md) | 只管没做完的 |
| [docs/PRRS_STATUS.md](docs/PRRS_STATUS.md) | 实现状态，逐条实测 |
| [docs/onboarding/README.md](docs/onboarding/README.md) | 新人导读 |
| [benchmarks/README.md](benchmarks/README.md) | 基准套件的规矩与队列 |
| [env/BASELINE.md](env/BASELINE.md) | 环境基线与指纹沿革 |
| [docs/PLAN_2026-09-28.md](docs/PLAN_2026-09-28.md) | 现行计划 —— connect 模式两层桥 |

截至 2026-09-28 的活跃工作：**connect 模式两层桥计划**（已知反应物与产物，
求其间通道）—— 代码草稿在仓库根 `PLAN_2026-09-28.py`，待分段并入 `src/prrs`；
**首次进入 committor**（主线 A）2026-09-26 完成并确认；静态性能审查
（`docs/REVIEW_2026-09-28.md`，PERF-01…28）待按"先量后建"分诊。

方法原稿：
[`mlp_native_perturbation_response_reaction_search.md`](mlp_native_perturbation_response_reaction_search.md)。

---

## 许可协议

本项目以 **GNU Affero General Public License v3.0 only**（`SPDX: AGPL-3.0-only`）
发布，完整协议文本见 [LICENSE](LICENSE)。

简言之：可以自由使用、研究、分享和修改本代码，但任何衍生作品 —— 无论是分发还是
仅通过网络对外提供服务 —— 都必须以同一协议公开完整对应源码。`-only` 表示仅按
AGPL-3.0 这一版本授权，不自动延伸到未来的 GPL/AGPL 版本。
