# PRRS — 扰动–响应反应搜索

[![License: AGPL v3](https://img.shields.io/badge/License-AGPL_v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)
![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue)
![Lint: ruff](https://img.shields.io/badge/lint-ruff-261230)

[English](README.md) | **简体中文**

**MLPATH 项目入口:**[架构与现状(中文)](docs/PROJECT_OVERVIEW.md) · [文档索引](docs/README.md) · [PARS 路径模块](pars/README.md) · [路径成熟层](mature/README.md)。
MLPATH 包含 PRRS 发现、PARS 路径选择/精修(MVP-1..5 已落地)和基于 shooting 的成熟层(`mature/`,M0+M1 已落地)。本 README 讲 PRRS;项目总览跟踪三者的边界与 connect 模式进展。

**PRRS 是为机器学习原子间势(MLP)设计的反应发现方法。**

> 以受控方式扰动一个分子,撤掉扰动,让体系在 MLP 势能面上自由响应,淬灭,
> 观察它落到哪个化学态。

PRRS **不需要**预定义反应坐标、预定义产物、或过渡态初猜。它先发现可能的反应
通道与产物;过渡态随后再定位与验证。

---

## 工作原理

```
种子态 → 受控扰动 → 短 MLP 轨迹 → 淬灭 → 分类
                                                  │
     源构象异构体 ◄── 同一化学身份 ──► 新化学态
                                                  │
                       一阶鞍点(n₋ = 1)──► TS 池
```

每次试验是一个小型数值实验,分两相:**施加相**,扰动以小步长加到内坐标上;
**自由响应相**,撤掉扰动,体系无偏弛豫。三种结局:

1. **回到源极小**——弹性证据。该方向在测试幅度下不易通向另一个盆地。
2. **弛豫到另一个极小**——发现一个盆地。构象异构体与新化学态显式区分:
   裁决的是化学连通性,不是能量或距离。
3. **恰好一个负 Hessian 本征值的驻点**——过渡态候选,进入 TS 池,随后沿两侧
   下坡方向验证它连接哪两个态。

### 支撑方法的设计选择

- **扰动作用在原始内坐标上**——键长、键角、二面角、片段集体运动、分子间接近
  ——而不是笛卡尔噪声。理由是实测的:早期乙醇试验里纯径向扰动全部回到起点,
  而三个真实构象盆地就在 4 meV 之外,因为成对径向移动在轨迹长度内根本耦合不进
  二面角。**扰动集必须张满化学上相关的自由度。**
- **施加是增量的**(每步 ≤ 0.05 Å/rad,每步重算梯度),因为大扰动是一串局部
  有效的小扰动。请求位移与附带损伤分开报告;共价键损伤与原子重叠是硬门。
- **对称性显式处理。**化学等价原子不会制造人工"新反应":身份按自同构不变量
  哈希,构象去重用对称感知 RMSD,立体签名只用构型手性。对称转子*降权而非删除*
  ——其中一个贡献了真实的甲基转动鞍点。
- **三层身份永不合并**:*化学态*(连通性 + 构型立体化学)、*构象库*(每态多个
  几何,因为反应可达性依赖起始构象)、*反应事件*(键差键名在两端点自同构群下
  规范化)。简并反应——如等价氢质子转移——是自环,不是新节点。
- **收敛判据按模式分辨。**固定力容差控制的是力有多小,不是离驻点有多近;乙醇
  二面角比键伸缩软约 150 倍,所以软模有自己的梯度/刚度判据,真自由模被如实记录,
  不被逼着"收敛"。
- **Fail closed。**任何无法判定的东西都带理由拒绝,绝不用默认值顶替;存在性
  主张必须由执行证明,不能由记录证明。

---

## PRRS 主张什么、不主张什么

| 问题 | 比较 | 检验的是什么 |
|---|---|---|
| **搜索**好不好? | PRRS 对参照搜索,在*同一* MLP 面上 | 方法 |
| **势**准不准? | MLP 过渡态对 QM | 能量模型 |

PRRS 负责发现通道、识别产物、定位并验证过渡态候选、构建网络。它**不**负责
最终定量势垒:势垒高度与虚频是所选势面的性质,由下游 QM 复核。把两种比较混在
一起会让误差无法归因——苯甲酰丙酮三模型对比里每个化学身份字段跨模型一致而
势垒高度差 3.7 倍,这正是本设计强制分开的分离。

---

## 状态

基准套件(`benchmarks/`)遵循一条规则:一格 = 一个反应体系 + 跑**之前**冻结的
验收判据 + 一个按判据机械判分的脚本。

| 阶段 | 体系 | 新检验轴 | 状态 |
|---|---|---|---|
| P0 | 乙醇 | 构象发现;框架 | 完成——13/13 |
| P1 | 丙二醛 | 反应事件——**简并**自环 | 完成——13/13,QM 复核 |
| P2 | 3-氧代丁醛烯醇 | 反应网络扩张 A→B | 完成——16/16 种子 |
| P3 | 苯甲酰丙酮 | 芳香性、螯合、竞争异构体 | 完成——判据冻结并通过 |
| b04 | Cl⁻ + CH₃Cl | **总电荷**、2 片段、重原子坐标 | 完成——判据 2026-09-13 冻结,TS 经 QM 证实 |
| b05 | 叔卤代烃 + 亲核试剂 | 取代对消除竞争 | 进行中——D2-HF v1.4 双种子全程已跑(2026-11-08);HS 括号规则冻结于 v7 |
| b06 | 亲电加成底物 | 极性加成 | 计划中 |
| b07 | 环己烯(逆 DA)· 丁二烯 + 乙烯(DA) | **片段数 1↔2** | 计划中——n1 Diels–Alder 探针脚本已落地 |

P0–P3 现已冻结为回归基准:其结果持续有效,正式验收转入 P4+ 独立 QM 协议。系列
刻意从反复的质子转移例子走向带电体系、分子间反应、竞争路径、重原子成断键和
片段数变化。

### 独立 QM 交叉核对

协议:几何与频率用 **r2SCAN-3c** → 同几何上 **ωB97M-V/def2-TZVPD** 单点。
(r2SCAN-3c 单独用会把 SN2 中心势垒低估约 11 kcal——几何是好的,能量必须复核。)

- **丙二醛(P1)。** 从 PRRS 鞍点出发的双向 IRC 收敛到预期极小,两端虚频均为零。
  MLP 与 QM 鞍点几何差 **0.0151 Å**(对称感知 RMSD)。MLP 势垒(419 meV)与
  QM(134 meV)差异显著——这是能量模型的误差,不是没找到反应几何。
- **SN2:Cl⁻ + CH₃Cl(b04)。** PRRS 找到的过渡态精修到 2.3645/2.3645 Å,与独立
  定位的 QM 鞍点四位小数一致。双向 IRC 到达预期反应物与产物,反应事件键名与
  PRRS 完全一致。

> **即使 MLP 不能准确复现 QM 势垒能量,PRRS 也能识别反应事件与过渡态几何。**

### 机器学习势的选择

P0–P3 冻结在 **MACE-OFF24** 上,不回溯重跑或重标注。b04 起默认 **MACE-POLAR-1-M**,
依据是 QM 交叉对比而非偏好:它是唯一在测试中同时扛住中性质子转移化学与离子
SN2 基准的候选。MACE-OFF24 对总电荷完全不响应(电荷 0/−1/+1 能量逐位相同),
这在离子化学里不合格。冻结结果从不自动充当后续阶段的参照。

---

## 安装

要求:**Python ≥ 3.10**;核心依赖仅 ASE 与 NumPy。

```bash
pip install -e .            # 核心:prrs 命令行;可在任何 ASE calculator 上跑
pip install -e ".[mace]"    # MACE 支持(mace-torch)
pip install -e ".[test]"    # pytest
```

冻结 P0–P3 基准所用的 MACE-OFF24 medium 检查点随仓库分发(`models/`,sha256
前缀 `e5ccf5837f685899`)。每个基准在 `benchmarks/b0*/model.json` 里声明自己的
模型——**从不隐式下载任何模型**。

### 两个环境,有意为之

**科学环境**(conda 规格在 `env/openmm_dev.yml`;基线与历史在 `env/BASELINE.md`)
产出全部发表数字,其包版本是每个结果的一部分。开发工具(ruff)放在单独的
`.venv-tools`(`make tools`),格式化永远动不了科学依赖。后端一致性曾在一次
openmm *patch* 发布上翻红——永远不要把开发工具装进科学环境。

## 快速开始

可选的局部极小/TS 候选精修:`prrs optimize`。支持 ASE calculator 与既有
OpenMM/MACE context,经 geomeTRIC;见[命令、输出与验证边界](docs/GEOMETRIC_REFINEMENT.md)。

```bash
# 最小自包含演示
prrs demo --output runs/demo --seed 42

# 一次搜索
prrs search input.extxyz --output runs/my_run \
    --config config.json --calculator mace --calculator-kwargs calc.json

# connect 模式:给定产物,每个微态记录到它的距离(只记录,不路由)
prrs search input.extxyz --product product.extxyz --output runs/my_connect \
    --config config.json --calculator mace

# 导出未精修的路径种子供下游 NEB / dimer / string 精修;淬灭迭代被排除——
# 它们不是动力学
prrs path --run runs/my_run --attempt <id> --output seed.extxyz --images 9
```

每次运行写 manifest(含 `implementation_sha256`,见[开发](#开发))和
`network.json` 反应网络。`runs/README.md` 记录哪些在盘运行算证据——
`network.json` *存在*不等于运行完成。

### 观看实时运行

```bash
python scripts/prrs_watch.py runs/my_run
```

只读终端 UI:增量 tail 每个试验的 `observations.jsonl`,轮询 `network.json`,
直接在终端里绘制演化中的势能面(真彩色)。它从不写 `runs/`;全分辨率 PNG 到
`.scratch/prrs_watch/<run>/`。

---

## 仓库布局

```
src/prrs/                 搜索核心——25 个模块,约 1.3 万行
  config.py               全部协议参数 + 校验;非法配置在启动前拒绝
  internal.py             原始内坐标、精确刚体片段转动、延拓位移
  chemistry.py            三层身份、自同构群、对称感知 RMSD
  canonical.py            经认证的 canonical 标号(individualization–refinement)
  perturbations.py        六个扰动族、对称商、分层方向预算
  network.py              注册表:准入、构象库、域门
  state.py                活跃子集图编码、Kabsch RMSD、响应分类
  reliability.py          每次力调用的 fail-closed 门控
  runner.py               一次试验:施加、按模式分辨的淬灭、Hessian/极小检查、鞍点下坡
  search.py               预算调度、幅度精修、TS 池、运行 manifest
  completeness.py         盆地成员资格的条件完备性证书
  product.py              connect 模式产物目标:只记录距离,从不路由
  window_audit.py         匹配响应视界的重放审计(从不为静态边界背书)
  geometric_backend.py    可选的 geomeTRIC 精修(ASE 或既有 OpenMM Context)
  geometric_cli.py        精修的 CLI 接线(惰性导入)
  provenance.py           calculator 身份,供续跑核验(不执行)
  schema.py               显式 network.json schema(schema 2)
  calculators.py          calculator 工厂
  openmm_backend.py       已测试的 openmm-ml/MACE 原语
  torch_guard.py          torch.load 的单一受审计包装
  validation.py           独立交叉核对路径(ase.vibrations)——搜索永不导入它
tests/                    36 个文件,388 个测试函数
benchmarks/               基准套件——每个反应体系一个目录(判据 + 判分器)
docs/                     工作记录:交接、验收文档、计划、实验日志
scripts/                  provenance、实现哈希、实时观察器
runs/                     实验输出——runs/README.md 说明哪些算证据
models/                   MACE-OFF24_medium.model(冻结的 P0–P3 基准模型)
env/                      科学环境基线:BASELINE.md、conda 规格、pip freeze
```

---

## 开发

```bash
make help      # 每个目标做什么
make check     # fmt-check + lint + test + digest,按此顺序——完整门
make digest    # 金样回归门
make hash      # 运行 manifest 会记录的 implementation_sha256
```

- **pytest 问代码对不对;digest 问另一个问题**——代码是否仍产生与盘上运行
  相同的数字(332 文件 / 2136 帧,每个差异归类为 `graph_invariant` /
  `parity_only` / `key_only`)。**零差异是零行为变化的唯一可接受证据。**有意的
  变更用 `--allow parity_only` 等放行——但放行确认的是差异的*形状*,不是新值
  正确;那需要化学判断。
- **`implementation_sha256`**(`make hash`)是每个运行 manifest 的溯源锚:
  `src/prrs/*.py` 的内容哈希,与 VCS 状态无关。永不手抄;记录在案的指纹跟踪在
  `docs/HANDOFF.md` 与 `env/BASELINE.md`。
- 最近一次完整门:**537 passed,7 skipped,digest 零差异**(2026-10-01;跳过是
  环境性的——门机器上无 MACE 检查点 / CUDA)。

### 项目规矩

这些是工作规则,不是建议;违反它们产出的结果曾被作废。

- **验收判据在数据存在之前冻结**,此后绝不就地修改;修改是显式修订,记录旧文、
  理由、以及修订时已看过哪些数据。
- **干预各自命名**,一次改进只归因于一个变化。
- **未完成 ≠ 失败。**分母是*发起的*试验数;未完成与崩溃如实记录并注明原因
  (`scripts/scan_status.sh` 机械区分三态)。
- **Fail closed。**判不了就带理由拒绝——绝不用默认值顶替。
- **门不因为挡住了想要的结果而搬动。**先量化门挡住了什么。
- **先测量再推断。**被记录的否证是资产,不是难堪(`docs/HANDOFF.md` §6)。

---

## 文档

| | |
|---|---|
| [docs/README.md](docs/README.md) | 文档索引——哪些现行,哪些只是记录 |
| [docs/HANDOFF.md](docs/HANDOFF.md) | 交接总览;冻结结论在这里 |
| [docs/DECISIONS_PENDING.md](docs/DECISIONS_PENDING.md) | 只有未决项 |
| [docs/PRRS_STATUS.md](docs/PRRS_STATUS.md) | 实现状态,逐条主张对应实测 |
| [docs/onboarding/README.md](docs/onboarding/README.md) | 新人指南 |
| [benchmarks/README.md](benchmarks/README.md) | 基准套件规则与队列 |
| [env/BASELINE.md](env/BASELINE.md) | 环境基线与指纹历史 |
| [docs/HF5_SUMMARY_2026-11-10.md](docs/HF5_SUMMARY_2026-11-10.md) | 当前线——HF5/HF6 小结与计划 |
| [docs/PLAN_2026-09-28.md](docs/PLAN_2026-09-28.md) | connect 模式双层桥——G1–G4 已落地,到产物即停的终止待做 |

截至 2026-11-11 的活跃工作:connect 模式、PARS(MVP-1..5)与成熟层 M0+M1 维持
原状。**D2-HF 线**已跑 v1.1 → v1.4:同图山脊被采样、HS / HC 并行波、HS 括号
规则冻结于 v7、v1.4 以双种子全程跑完并带丢失模式诊断;fp32-Hessian 诊断与
F3 假山脊审计已关闭,救回的中点候选是新的 E2 TS。下一批构建块:HF5 跨体系
(`docs/PLAN_HF5_CROSS_SYSTEM_2026-11-08.md`,小结在
`docs/HF5_SUMMARY_2026-11-10.md`)与 HF6 dimer(`docs/PLAN_HF6_DIMER_2026-11-10.md`);
n1 Diels–Alder 探针脚本在 `benchmarks/n1_diels_alder/`。

方法原始笔记:
[`mlp_native_perturbation_response_reaction_search.md`](mlp_native_perturbation_response_reaction_search.md)。

---

## 许可证

本项目仅按 **GNU Affero General Public License v3.0**(`SPDX: AGPL-3.0-only`)
授权。全文见 [LICENSE](LICENSE)。

简言之:可以使用、学习、分享和修改本代码,但任何衍生作品——无论分发还是
通过网络提供服务——必须在同一许可下发布,并附完整对应源码。"only"表示
代码恰按 AGPL-3.0 授权,不自动适用未来的 GPL/AGPL 版本。
