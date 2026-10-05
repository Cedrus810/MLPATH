# b04 判据

**判据文本依据用户 2026-09-09 裁决逐条写定。** 冻结＝逐字不改；改动走文末修订记录，
并注明**修订发生在看到哪些数据之后**。

> **已冻结 · 2026-09-13**（用户落笔）。此刻起本文件逐字不改，任何改动走文末修订记录。
> 冻结时 `check.py` 逐条对得上本文文本，自检 5 例通过；`§5` 第三条的实现口径见该节末的实现说明。
> 冻结**不追认**任何既有 run：`probe1–4` 都在冻结之前跑的，它们是探针不是验收。

配置：`config.json`（同日对齐历史事实：`direction_mode: "sampled"`、
`saddle_search_max_rise_eV: 10.0` 为非 binding 探针包络）· 模型：`model.json`
（MACE-POLAR-1-M）· 实现树：由每次运行的 `implementation_sha256` 记录，不在此手抄。

---

## 0. 这一格问什么

> **取代通道**上，反应物与产物之间的连通性，有没有**正证据**？

`Cl⁻ + CH₃Cl` 是同一性取代，两端在自同构下等价，所以这条通道是**自环**。
问法与 b01 的退化质子转移同形（`docs/HANDOFF.md` §3）。

教科书说 SN2 是协同的。**那句话不在这里**，它在 `reference.json` 里，
作为一条可被本基准证伪的预期。

## 1. 只接受正证据（裁决的核心）

三个判定值，**只有前两个是正证据**：

| 判定 | 需要什么 |
| --- | --- |
| `observed_stepwise` | 该通道上确认了一个极小，**并且**有证据把它接入这条通道 |
| `observed_direct` | 一个 order-1 鞍点，其**两侧下降 / IRC 直接到达该通道的两端** |
| `unanswered` | **其余一律**，并附细分理由（§4） |

**"没找到中间体"不是"协同"。** 有限抽样只给"未观测到"，不给"不存在" ——
这是本轮停线的理由，写在 `docs/DECISIONS_PENDING.md` 首节。
所以 `concerted` / `stepwise` 作为**结论标签从判据里退休**：它们只能是对
`observed_*` 的解读，不能是 pass/fail 的条目。

**两条都成立时两条都报。** 一条通道上同时存在直接路径与分步路径不是矛盾，
是这条轴上最有信息量的结果，不得压成单一机理标签。

### 1.1 `exhaustive` 也推不出协同

即便 `direction_mode: "exhaustive"`，`searched_no_hit` 能证明的只有
**"在被枚举的提案基内没有命中"**。它**不自动**推出协同 ——
协同仍然必须由 `observed_direct` 给出。

（b04 的三次探针都是 `sampled`，见 `config.json`。）

## 2. b04 的 `observed_direct` 逐条（退化自环的正确写法）

**不得要求两个不同的 chemical node。** 同一性取代两端是同一化学态，
要求两个节点等于要求这个反应不是它自己。四条一起成立才算：

| # | 条目 | 读哪个字段 |
| --- | --- | --- |
| D1 | 存在 `saddle_order == 1` 的鞍点候选 | `ts_candidates[*].saddle_order` |
| D2 | 两侧下降到达**两个不同的 microstate** | `ts_candidates[*].connects[*].microstate`，两条记录的值不同 |
| D3 | 两侧端点都是**确认过的极小** | `connects[*].status` ∈ {`existing_microstate`, `new_microstate`, `new_chemical_node`}。这三者只在 `descend_saddle` 的 `confirm_minimum` 通过后才可能出现（`runner.py`；未通过则记 `not_a_minimum` 或 `quench_failed`） |
| D4 | 该通道有**明确的键变事件**，且被这个鞍点支持 | `reaction_channels[*].broken` / `.formed` 非空，`event_key` 存在，且 `.ts_support` 含该 TS 的 id |

**D2 用 microstate 而不是 chemical node，是这一格的关键。** `chemical_key` 在自环两端相同，
所以节点数扛不动这个证据；扛得动的是"两个不同微态 + 一个具名的键变事件"。

**不进判据的**（规矩 2）：势垒、虚频数值、与 QM 的偏差。
两个候选模型在这个反应上过渡区能量差 2.9 倍、最低本征值差 12.8 倍
（`reference.json`；波数写法下 √12.83 = 3.58 倍，同一个事实）。
写进去就是把模型的分歧记到方法头上。那两个数**本身**也不能按字面读：
D3h 点的残余梯度从未测量（规矩 6）。

## 3. `observed_stepwise` 逐条

| # | 条目 | 读哪个字段 |
| --- | --- | --- |
| S1 | 通道上存在一个 `confirm_minimum` 通过的极小，且它**不是**该通道的任一端点 | 微态记录 + `connects[*].microstate` |
| S2 | 有证据把它接入这条通道：从它出发的鞍点或下降记录连到该通道的端点 | `ts_candidates[*].connects`、`reaction_channels[*].observed_by` |

**S1/S2 不成立不等于协同**（§1）。b04 若两条都不成立，该问法的分步侧记 `unanswered`。

## 4. `unanswered` 的细分理由（必须给出，不得只写 unanswered）

| 理由 | 什么时候用 |
| --- | --- |
| `basis_uncovered` | 候选基里就没有能移动相关坐标的方向 |
| `candidate_not_selected` | 是候选，但没被选中（`direction_ledger.uncovered_ids`） |
| `proposed_not_executed` | 被选中但预算没轮到执行（`direction_ledger.not_executed_ids`） |
| `gate_refused` | 被门拒绝，**附门名与理由**。包括"极小已确认但被化学标签挡在网络外"那一类 —— 结构在 `network["withheld_candidates"]` 里可检索，QM 复核照样能取用 |
| `searched_no_hit` | 真的搜过且没命中。**只有此项才是"搜过"** |

## 5. 负记录何时可采纳（按微态裁决）

**逐微态判定，不按整次运行。** 某个微态的负记录可采纳，当且仅当该微态同时满足：

```
direction_ledger_per_microstate[<id>].uncovered_ids      == []
direction_ledger_per_microstate[<id>].not_executed_ids   == []
且该微态上没有与该通道相关的 gate refusal
```

三条齐了，该微态上的 `searched_no_hit` 才是**一条关于该微态、在被枚举的提案基内的、
可采纳的负记录** —— 与 §1.1 同一口径。**它仍然不是关于势能面的陈述**：
提案基是有限的、`sampled` 下还是被截断的，所以它不支持任何全局不存在断言。

**一个深层微态被截断，只让该微态 `unanswered`，不抹掉别处的正面发现。**
正证据是局部的、可累加的；负结论不是。

**b04 probe3 在这一条下的实际状态**：基覆盖 14/28，`cap_bound false`，
截断者是 `direction_quota_fraction` —— 所以 **b04 的负侧一律 `unanswered`**，
而 `observed_direct`（§2 四条）不受影响，它是正证据。

**实现说明（2026-09-13 冻结时记，属修订记录范围）**：第三条"该微态上没有与该通道相关的
gate refusal"要求**逐微态**归属，而 `network.json` 只写 run 级的 `outcome_counts`，
不把 refusal 归到微态。所以 `check.py` 保守实现它：**run 内有任何 refusal 就 censor 全部微态**；
**`outcome_counts` 字段缺失同样 censor**（"没有记录"不是"没有发生"，与 no-ledger 分支同一口径）。
这是实现比文本**更严**，不是更松 —— 它只会让可采纳的负记录更少。
要拿到文本原本的逐微态口径，需要运行侧把 refusal 归属写进记录，那是一次代码改动，不是判据改动。

## 6. 分母与记账（规矩 4）

一份报告里，分母是**发起数**；被门拒绝的记被拒并附门名与理由，**不记成失败**。
`direction_ledger` 的两个单位不得相加：`directions_*` 回答"这根坐标试没试过"，
`amplitudes_*` 回答"沿它走到了哪一步"。

## 修订记录

| 日期 | 改了什么 | 看到哪些数据之后 |
| --- | --- | --- |
| 2026-09-09 | 判据文本首次写定（本文件此前只有流程） | `runs/b04_probe3_rigidbond`（0 拒绝 · 1 个 bracket · ts0000 复现）· `docs/experiments/b04_qm/`（双向 IRC、两端 0 虚频、`event_key` 逐字一致）· 基覆盖 14/28 的实测 |
| 2026-09-13 | **冻结生效。** 文本本身只加了两处：文首的冻结戳、§5 末的实现说明（记 `check.py` 把第三条实现得比文本更严，及其原因）。判据条目本身一字未改 | 三路审计核实 `check.py` 与本文的两处不一致并修掉：① clean-ledger 分支曾打印 "this IS a statement about the surface"，与 §5 "它仍然不是关于势能面的陈述"直接冲突；② §5 的三条只机械实现了两条，第三条（gate refusal）缺失。两处都已修，自检 3 例扩到 5 例。另：`benchmarks/README.md` 规矩 5 与 `TAXONOMY.md` 同日改写（"不存在 ⇒ 协同"是无效推理形式），本格判据的引用源头因此才一致 |
