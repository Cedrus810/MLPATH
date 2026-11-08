# 施工计划：延续试步松弛失败时的有界折半重试（TODO 第 10 条）

> **10-29 更新：** b05 yayoi08 整跑中折半重试未观察到候选收益，用户指示改为“未收敛试步首次失败，就结束该侧”（HF-EF-branches-v4，DONE_LOG 10-29）。下文 §3 第 2、5 条的折半重试已不是现行策略；同机 v3 / v4 对照用 `benchmarks/explore_d2/hs_policy_compare.py`。

日期：2026-10-28。状态：用户「是的」，已施工为 v3；10-29 改为 v4（首次失败即终止该侧）；同机对照已做。结果见 [HF_WALK_POLICY_RESULTS_2026-10-29.md](HF_WALK_POLICY_RESULTS_2026-10-29.md)。
依据：[HF 设计记录](TS_HESSIAN_FREE_DESIGN_RECORD_2026-10-24.md) §5.2（先缩小步长分别延续两侧；不能形成连续区间就保存双侧几何、返回明确的分支未决）；[跳支定位计划](PLAN_HF_JUMP_2026-10-25.md) J1（`continue_branches_ef` / `_walk` 的折半与终止）；DONE_LOG 2026-10-28「HS `transverse_not_converged` 离线诊断」。

## 必读（开工前全部读完）

- `docs/TS_HESSIAN_FREE_DESIGN_RECORD_2026-10-24.md` §5.2、§5.3、§7：延续与终止的原始要求、数值约束、验收方式。
- `docs/PLAN_HF_JUMP_2026-10-25.md` §3 J1：`_walk` 的“疑似跳变 → 折半 → 低于 `min_step_q` 才终止”是本计划要扩展的机制。
- `docs/DONE_LOG.md` 2026-10-28 诊断一节（含更正）：178 项的证据，以及哪些说法未确认。
- `docs/TODO.md` 第 9、10 条：本计划只做第 10 条；第 9 条（中点阈值）在本计划之后。
- `explore/src/explore/ts_branch.py`（`_walk`、`continue_branches_ef`）、`explore/src/explore/ts_ef.py`（`_relax`、`_RelaxState`、`_Stop`）、`explore/tests/test_ts_branch.py`。

## 1. 用户原话（2026-10-28）

- “我建议先处理第 10 条的退出语义，再处理第 9 条的数值判据：延续试步松弛失败后，允许有界折半重试，保留最后已接受点并回滚试探曲率；未收敛帧只作诊断，不能用于反力变号括区。端点失败则需单独处理，因为还没有可回退的已接受延续点。”
- “178 项支持‘疑似失支’，尚不能全部确认换盆。大位移、能量下降也可能来自同一分支上的缓慢松弛。应通过缩小延续步长后的响应判断。”

## 2. 现状

`_walk` 每一步从当前已接受点出发，在 `q + sign·step` 做约束松弛（`_relax`，最多 `max_relax_steps` 200 步，横向阈值 `transverse_fmax_eV_A` 0.003）。返回后按位移判断是否疑似跳变：疑似就把 step 折半重试，step 低于 `min_step_q` 仍跳才终止该侧（`terminated`）。

`_relax` 在 200 步内没收敛时抛 `_Stop("transverse_not_converged")`。`_walk` 不接这个异常，它一路冒到 `continue_branches_ef`，整个 HS 以 `unresolved:transverse_not_converged` 结束，折半机制没机会起作用。b05 HF 四块中 166 项属于这种情况（178 项减去端点失败的 12 项）。

已有的、本计划直接依赖的性质：试步用 `state.copy()`，只有被接受才提交；`_relax` 只在收敛时写 state。所以失败的试步不会污染已接受点的几何和曲率。

## 3. 改动

只改 `explore/src/explore/ts_branch.py`：

1. `_walk` 中把试步的 `relax(...)` 包在 `try` 里，只接 `_Stop("transverse_not_converged")`；其他 `_Stop`（`above_ceiling`、`constraint_not_converged`、`degenerate_coordinate`）、`BudgetExhausted`、`GateRejected` 照旧向上抛。
2. 接到时：当前点、`state`、`rates` 都不变（丢弃 `trial_state`）；若 `step > min_step_q`，`dq = step / 2` 继续；否则该侧终止：`status = "terminated"`，`termination = "relax_unconverged"`，并记录诊断（试步 q、step、最后一帧的横向力 / 完整力 / 能量变化 / 相对当前点的对齐 RMSD、位置）。疑似跳变终止的一侧记 `termination = "jump"`，原有 `jump` 字段不变。
3. 失败试步的帧照旧进入 `result.history`（`_relax` 逐帧追加，计费不变），但**不进** `side["points"]`，不参与 `rates`，也不参与反力变号判断；所以不可能用于括区或交给 `refine_ts_ef`。
4. 总状态：任一侧找到最大值仍是 `continuous_maximum`（之后照旧交给 `refine_ts_ef`）。都没找到时：两侧都是 `reached_end` 或 `termination = "jump"`，记 `unresolved:hysteresis`（原语义）；有一侧因 `relax_unconverged` 终止，记新的 `unresolved:walk_relax_unconverged`，不冒充 hysteresis。
5. 有界：每次失败最多折半到 `min_step_q`。默认 `step_q` 0.01、`min_step_q` 1e-3，即同一位置最多再试 4 次，每次最多 200 步松弛，最坏多花约 800 forward；`max_points` 400 照旧兜底。不新增设置项。
6. 端点（`continue_branches_ef` 开头两端的松弛）不改，仍是 `unresolved:transverse_not_converged`。原因是此时还没有已接受的延续点可以回退（用户原话）。处理方式另议。
7. 算法版本 `HF-EF-branches-v2` → `v3` 写进 provenance。D2 规则版本不变（与 10-28 括区修复同一做法），在 DONE_LOG 注明从哪天起的 HS 用 v3。

## 4. 测试（`explore/tests/test_ts_branch.py`）

`_walk` 接受 `relax` 可调用对象，所以前四项可以用假的 `relax` 直接测，不依赖解析势的巧合：

1. 大步失败、小步收敛：该侧折半后继续，走到端点；失败帧不在 `points` 里；失败前后的 `state`（几何与曲率）完全不变。
2. 某个 q 附近无论步长多小都失败：该侧 `terminated`，`termination = "relax_unconverged"`，诊断字段齐全；试探次数符合 `step_q → min_step_q` 的折半次数。
3. 失败试步的最后一帧反力与当前点异号，但已接受点不变号：不得返回 `maximum_bracketed`。
4. 其他 `_Stop`（如 `above_ceiling`）与 `BudgetExhausted` 仍向上抛，行为不变。
5. 解析台（现有折叠势）把 `max_relax_steps` 设得很小，迫使部分试步失败：`continue_branches_ef` 端到端得到 `walk_relax_unconverged` 或折半后继续；外部计数 == 账本（失败试步的 forward 照常计费）。
6. 原有 `test_ts_branch.py`、`test_ts_ef.py`、`test_d2_hf.py` 全过（无失败试步时行为逐字节不变）。

变异检验：把第 3 条改成“失败帧也计入反力判断”，以及把第 2 条改成“失败即终止不折半”，对应测试必须失败。

## 5. 验证（真实模型）

| 步骤 | 内容 | 资源 | 预计耗时 | 依据 |
| --- | --- | --- | --- | --- |
| 冒烟 | 从 v1.1 sA 的 166 项中按 journal 顺序取前 3 项，用新代码重跑 HS（`_hs_task` 同设置，只读 journal，每项上限 `hs_cap` 3000） | 本机 181，2080 Ti | 约 5 分钟 | 10-28 括区冒烟实测 0.09–0.12 s/forward；每项最坏约 3000 forward |
| 有界对照（另报计划表、等“跑”） | 166 项中按固定种子抽 30 项，新旧代码同输入**在同一台机器上**各跑一次 HS（冒烟显示跨硬件不可复现）；看新状态分布（继续走到端点 / 找到最大值 / `walk_relax_unconverged`）、折半后在多小步长上收敛、每项 forward | 180/181/183 之一，1 GPU | 约 30 × 2 × 1500 × 0.1 s ≈ 2.5 h（粗估，冒烟后按实测单项耗时重估） | 同上 |

判读只按“缩小步长后的响应”下结论：折半后以小位移收敛、继续前进，说明原来是步长或慢松弛；降到 `min_step_q` 仍失败，才记为该侧在此终止。两者都不等于“已确认换盆”，也都不等于找到了 TS。

## 6. 不做

- 第 9 条（中点阈值、σ 重测）。
- 端点失败的处理。
- 不改 `_relax`、`refine_ts_ef`、阈值、`max_relax_steps`，不新增门或设置项。
- 不改 D2 调度。

## 未确认（Claude 推断）

- 166 项中有多少在折半后能继续前进，未知；也可能大多数降到 `min_step_q` 仍失败，那样只是把“整个 HS 结束”变成“该侧终止、另一侧继续”，额外花约 800 forward。
- refine 额度用满时，每项多花的 forward 会挤掉后面的项；净效果要看有界对照。
- `walk_relax_unconverged` 这个状态名是我拟的。
