# b05 — 叔丁基氯 + OH⁻

**这一格问什么**（按 `../README.md` 规矩 5，标签在输出侧）：

> **取代通道**上，反应物与产物之间是否存在一个 `confirm_minimum` 通过的极小？
> **消除通道**上，同一个问题，独立回答。
> 两条通道独立报告正证据：`observed_stepwise` 或 `observed_direct`；其余为 `unanswered` 并注明原因。
> 有限搜索未找到极小不能推断为协同反应。正式评定以本格 `CRITERIA.md` 的冻结/修订状态为准。

教科书会说叔丁基氯走离子化（SN1/E1，分步，经碳正离子），强碱下走 E2。
**这些话都不进判据。** 它们在 `reference.json` 里，作为可被本基准证伪的预期。

## 为什么它不能被 b04 代替

| | b04 | b05 |
| --- | --- | --- |
| 通道 | 一条（取代） | **两条，竞争同一个底物** |
| 反应中心 | 甲基碳（一级） | **叔碳** |
| 亲核试剂 | Cl⁻（只有亲核性） | **OH⁻（同时是碱）** |
| 产物碎片数 | 2 | 取代 **2**、消除 **3** |
| 断的键 | 一根 C–Cl | C–Cl **与** C–H（哪一根是输出） |

**取代与消除必须在同一格里同时报**（`../TAXONOMY.md` B 组）：它们争夺同一个底物，
分开做会把"这个体系走哪条路"变成两次独立的单选题，而那不是化学里的问题形状。

**为什么碱选 OH⁻ 而不是 Cl⁻。** Cl⁻ 碱性弱，只给取代通道 ——
那样 b05 就只是"b04 换个更胖的底物"，**违反规矩 3（一格一条轴）**。
OH⁻ 两条通道都开，且元素 H C O Cl 全在 `chemistry.VALENCE` 的七个之内，
总电荷 −1 与 b04 相同，**b04 那套启用工作原样继承**（`_stamped_with_charge`、
`--calculator` 钩子、`closed_shell_only=False` 的模型域声明）。

## 输入几何：放置规则是机械的，而且被强制执行

规则写在 `build.py` 的 docstring 里，三条，任何一条不满足 `build.py` 直接退出：

1. 背面 C3 轴上（三个甲基保持等价，不偏向任一面）
2. 氧在所有底物原子的成键距离之外，且键图报 **2 个碎片**
3. 在 `pair_cutoff_A` 内既能看到中心碳、也能看到至少一个 β-H

落点是 **b04 自己的 3.6 Å 和 b04 自己的 4.0 cutoff，没有任何旋钮被调**。
实测：O···C 3.600 Å、三个反式共平面 H 各 2.463 Å、2 个碎片、最近接触 2.463 Å
对成键阈值 1.16 Å。

**这一格的草稿违反过这条规则两次，两次都是规则自己抓住的**，写在这里因为它比结论有用：

| 草稿 | 错在哪 | 谁抓住的 |
| --- | --- | --- |
| 甲基极角写成 `π − 109.47°` | 三个甲基与 Cl 同侧，Cl···C = 1.93 Å，键图多长三根键、Cl 度数 4 | **身份层**：gate 7 报 `degree_exceeds_valence` |
| 为"避开氢键"移到 4.30 Å 并把 cutoff 提到 4.6 | 用错了阈值（`hbond_acceptor_A` 只管 N/O/F–H 供体，甲基 C–H 不是供体） | 读代码 + 读 b04 落盘 |
| 把上一条与"提高 cutoff"合成一次改动，并据此断言 cutoff 是 gate 6 失败的原因 | **归因错误。** 同一几何下 4.0 与 4.6 选中数都是 24、都 cap 绑定，旧版 gate 6 两者都 FAIL；那次 PASS 属于坏几何 | 用户 2026-09-09 只读核查 |

基覆盖的实测（同一修正后几何，唯一变量是 cutoff）：

```
                        候选   选中   覆盖   cap_bound   截断者
b04 probe3  ( 6 原子)     28     14    50%    false      direction_quota_fraction = 0.5
b05        (16 原子, 4.0) 152     24    16%    TRUE       max_directions
b05        (16 原子, 4.6) 191     24    13%    TRUE       max_directions
```

**两个 cutoff 下选中数都是 24、都被 cap 绑定，旧版 gate 6 在两者下都 FAIL。**
所以提高 cutoff **不是**那次 FAIL 的原因 —— 原因是方向预算，而 cutoff 只把覆盖从
16% 稀释到 13%。（本文件早先一版写成"cutoff 25 → 37 使门 PASS → FAIL"，
那是把**几何修复**与**cutoff 变更**混成了一次干预：25 是甲基极角写错时的候选数，
那次 PASS 来自坏几何。已更正。）

**b04 的 50% 更要紧**：6 个原子、cap 从未绑定，被截断的是配额分数本身。
覆盖截断不是大分子才有的问题，它在每个尺寸上都在。

`perturbations.propose` 的 docstring 早写了这件事（"A fixed total budget cannot hold a
fixed hit rate while candidates grow with the molecule"）—— b05 是第一次把它量出来，
而量出来的东西触发了**方向合同裁决**（用户 2026-09-09，见 `docs/DECISIONS_PENDING.md`）：
**这一格因此暂停。**

## 状态

> **2026-09-22 更新：gate 9 已 PASS（1149 ms），必需门 7/7。** 下面这一段是
> 2026-09-09 的状态，原样保留作为记录 —— 它描述的枚举器已被同年 09-09 的逐点回溯
> 剪枝换掉，本节此后一直是过期文本。

~~`blocked`（2026-09-09）。**必需门 5/6，gate 9 FAIL**~~ —— 自同构群在
`automorphism_limit` 内枚举不完：`|Aut| = 1296`（三个等价甲基），需要 >2×10⁶ 个候选与
**19.8 秒**，默认是 20 000。第一次探针因此在第 1 个 trial 就以 `IdentityUnavailable`
停下（fail-closed，按设计）。

**这是枚举器的 scaling 问题，不是旋钮设小了**：代价按最大轨道的阶乘走
（轨道 3/6/9 → 0.8 ms / 30 ms / 19.8 s），而 b05 的九氢轨道来自叔碳本身。
两条出路（提限并每结构付 20 秒 / 改进剪枝）的代价不同，是设计决定，见 `preflight.md`。

**判据未冻结，冻结是用户的行为。**

前置继承 b04 的三处缺口 —— **都已解决**，所以这一格没有自己的阻塞项：
域门那条按层归位到 `model.json` + preflight，`model_type` 走 `--calculator` 钩子
（**共用 b04 的工厂，不复制一份**），电荷由 `_stamped_with_charge` 盖章、冲突拒开跑。

## 已知的、尚未测到的

- **`r_C_H_beta` 是 `candidate_not_selected`**：基里有 **5 个**候选方向能移动那根 C–H，
  seed 42 下一个都没被选中（gate 6 逐坐标打印基与选中两列）。
  **这不是基的缺口**（那会是 `basis_uncovered`，要改输入或家族），
  是方向预算的产物 —— exhaustive 或更大预算就能覆盖。
  两个词分开是因为修法不同，而这一格是最可能撞上它的地方：
  按裁决，负结论在覆盖被截断时一律 `unanswered`。
  （本文件前一版写成"没有任何方向直接驱动它"，已更正。）
- **三碎片产物态从未被任何模型验过**（`model.json` 的 `unverified_for_this_grid`）。
  消除产物是烯 + 水 + 卤离子，而 −1 必须落在氯上而不是抹开。gate 8 量的是这件事。
- **POLAR-1 的离子验证只有一个数、一个阴离子**（b04 的 Cl⁻ SN2，中心势垒 +2.09 kcal）。
  羟离子、C–H 抽提、烯产物**都没有**。任何来自这一格的能量陈述都继承这一句。
- **`collateral_gate` 的"碎片数减少"分支**（两碎片并一碎片）到 b04 probe3 为止仍未触发。
  这一格的取代通道正是 2 → 1，所以它可能第一次被走到。
- **成本未测。** b04 是 6 原子 39.5 s/trial；这里 16 原子，MD 与解析 Hessian 都更贵，
  探针预算与墙钟要先用一次小规模运行量出来，不从 b04 外推。

## 怎么跑

计算在 183（kasuga02）或本机（kasuga01），见 `env/BASELINE.md`。

```bash
python benchmarks/b05_tertiary_halide/build.py benchmarks/b05_tertiary_halide/input.extxyz
python benchmarks/b05_tertiary_halide/preflight.py            # 1/2/3/6/7，CPU
python benchmarks/b05_tertiary_halide/preflight.py --models   # 4/5/8，一块 GPU

# config.json 带 // 注释，load_json 不吃，先剥掉（同 b04 的 RUN.md）
scripts/prrs_remote.sh search benchmarks/b05_tertiary_halide/input.extxyz \
    --output runs/b05_probe1 --config .scratch/b05_config.json \
    --calculator benchmarks.b04_sn2_chloride.calculator:polar1
```
