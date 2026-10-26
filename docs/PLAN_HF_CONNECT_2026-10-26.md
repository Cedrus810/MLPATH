# 施工计划：只用 E/F 的两侧连接追踪（HF 记录 F7；`trace_sides_ef`）

日期：2026-10-26。用户 10-26：要我按 HF 记录 F7 写只用 E/F 的连接追踪施工计划，回答「可以」。
依据：[HF 设计记录](TS_HESSIAN_FREE_DESIGN_RECORD_2026-10-24.md) §2（两侧下降是连接追踪，保存启动方向、位移尺度、下降记录与端点；不等同完整 IRC）、§4 F7、§6（`trace_sides_ef`、结果维度分开）；[HF-2/3 原语](HF_EF_LOCATOR_2026-10-25.md)；[跳支定位计划](PLAN_HF_JUMP_2026-10-25.md)；DONE_LOG 2026-10-25→26 各节。
本文只是计划。施工进展（2026-10-26）：K1–K3 完成，5 个热启动候选 5/5 E/F 连通、键变化与旧 C2 一致；第一轮后改了起点条件与端点资格判据（见 DONE_LOG 同日）。K4 未做。

## 必读（开工前全部读完）

- `docs/TS_HESSIAN_FREE_DESIGN_RECORD_2026-10-24.md` §2、§4 F7、§6：本计划的原文要求与结果维度。
- `docs/DONE_LOG.md` 2026-10-25→26 的 τ 1e-3 与曲率复用两节：现有 5 个候选、旧核验的成本与失败原因。
- `src/prrs/runner.py` 的 `descend_saddle`、`follow_unstable_mode`、`quench`、`confirm_minimum`：被替代的旧链。
- `explore/src/explore/ts_ef.py`、`ts_branch.py`：搜索结果里有什么（`bracket`、`history`、`direction_cartesian`、分支点）。

## 1. 现状（实测）

旧连接核验 = `confirm_minimum`（完整 Hessian）→ 每侧 `descend_saddle`：

| 环节 | 依赖 | 本批实测 |
| --- | --- | --- |
| 起步方向 | Hessian 不稳定模 × `irc_step_A` 0.15 Å | — |
| 下降 | FIRE + 扭转 Newton 打磨（`quench`，E/F） | b05 三个候选中两个因打磨不收敛失败（`soft_mode_not_converged`、`soft_mode_moved_on_the_last_round`，后者未满足放行条件） |
| 端点判定 | `confirm_minimum` 完整 Hessian；判为鞍点就沿其不稳定模接力 | 一个候选接力 4 次后失败 |
| 成本 | — | b04 173、b05 `6-9` 715–775、其余 1144–2120 forward；搜索本身 215–284 |

`same_basin`、`classify`、`chemical_key` 只用键图、能量与对齐几何，可直接复用。

## 2. 方法

搜索结果本身已经给出两侧：候选 x* 是 q 约束极小支 E*(q) 上的极大点，分支在 x* 两侧的点都是这条支上更低的结构。

1. **起点取自搜索，不用 Hessian 模（F7 原文：「TS 两侧的启动方向优先来自搜索本身」）。** 每侧取同一支上离 x* 最近、且反力已分辨的点：|λ| ≥ 3σ（σ 是该体系实测的同支反力噪声）。来源是 `EFResult.history` 与 `BranchResult` 的分支点，都是已付费的约束极小结构。
   - 推断（未确认）：E*(q) 在 x* 处是极大，所以分支切向在不稳定方向上有分量，两侧点局部位于分界面的两侧。这一条要靠 C3 的实测核对，不当作前提。
   - 若某侧没有满足条件的点，则沿 `direction_cartesian` 以 ±a 位移起步。幅度 a 从小到大取，直到能量下降超过 3 个能量量子，取值与历史都记录。
2. **自由下降只用 E/F**：撤掉约束，用 BFGS 下降到完整力 ≤ `quench_fmax_eV_A`（0.03，与旧链同一判据），并启用 `ts_ef` 已有的曲率复用。不做扭转打磨：扭转打磨只影响软转子是否完全静止，不改变键图，而连接判定只看键图与盆。逐帧记录能量，下降中能量单调非增。
3. **端点资格用 E/F 有限扰动**：对端点做 k 次 σ_x 扰动后重新下降，全部回到同一盆（`same_basin`，同一键图）才记为 `ef_stable`。这只是 E/F 等级的稳定性证据（F7 原文：「E/F 松弛和有限扰动检查仅提供相应等级的稳定性证据」），不写成“极小已认证”。k 与 σ_x 先在解析台上定，写明依据。
4. **连接判定**：两侧都得到 `ef_stable` 端点，且不同盆，记 `connected`，再用 `classify` 给出事件标签与指纹；同盆记 `same_basin`；任一侧失败记 `failed:<原因>`；结果说不清的记未决。
5. **性质检验（Hessian）作为独立可选维度**，不作门：`validate_ts` 可单独对候选和端点求完整 Hessian，结果另存，不改写连接结果（HF 记录 §6 表「MLP Hessian 与独立 QM 的结论也不互相覆盖」）。

不做：不调用 `confirm_minimum`、`curvature_spectrum`、`descend_saddle`；不加接力；不改 PRRS 共享原语；不改 C1/C2；不接 D2（属 HF-4）。

## 3. 施工步骤

| 步 | 内容 | 完成条件 |
| --- | --- | --- |
| K1 | `explore/src/explore/ts_connect.py`：`trace_sides_ef(search_result, atoms, factory, config, ledger, task, *, settings)` → `ConnectionResult`（`to_dict()`）。逐侧保存起点来源（分支点 / 切向位移 a）、下降帧、端点 E/F、扰动复核、键图；分段计费 `connect_ef`。 | — |
| K2 | 解析台 `explore/tests/test_ts_connect.py`：二维双盆势（含鞍点），覆盖：两侧到达两个极小；起点在两侧；模型 Hessian 禁用；端点为鞍点时扰动复核能识别；预算 / 能量上限截断；计费与账本一致。 | 全部通过；`test_ts_ef` / `test_ts_branch` 不变 |
| K3 | 真实模型有界对照 `benchmarks/explore_d2/hf_connect_bounded.py`：对 `results/2026-10-26/hf_warm_start/` 的 5 个候选（只读）跑 `trace_sides_ef`。对照量：端点键图与旧 C2 两侧 `chemical_key` / 指纹是否一致，以及连接成本与旧核验 forward 的比较。 | b04 连通 SN2；b05 `6-9` 连通 E2（与旧指纹相同）；另三个 b05 如实报告；成本逐项列出 |
| K4 | 可选：对 K3 端点各求一次完整 Hessian，只作独立性质对照，看 `ef_stable` 端点是否都是极小。 | 一致率如实报告；不一致的逐个列出，不回改 K3 结论 |

K3 的成本依据：b05 自由下降每侧预计 100–300 forward，这一点没有实测，旧 quench 的 forward 数可从旧核验账本里拆出来核对；扰动复核 k × 约 30。冒烟级，在本机或 180/181 上跑。

## 4. 与旧链的关系

旧 `certify` / `descend_saddle` 保留作基线与历史复现。新连接结果是 HF 结果的独立维度（HF 记录 §6），不重标旧 C2。

## 未确认（Claude 推断）

- 两侧分支点局部位于分界面两侧（§2 第 1 点）：依据只是 E*(q) 在候选处为极大；K3 用 b04 / b05 `6-9`（旧链已连通的两例）核对。
- 不做扭转打磨不影响键图判定：依据是打磨只调整软转子；K3 若出现同键图不同盆，再回来看。
- `quench_fmax_eV_A` 0.03 沿用旧链；它在 fp32 下是否合适，未单独测过。
