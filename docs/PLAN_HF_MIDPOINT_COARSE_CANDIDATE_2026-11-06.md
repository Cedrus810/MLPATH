# 施工计划：粗阶段出候选（P1，HF-EF-bracket-v6，TODO 第 9 条 9c）

日期：2026-11-06。依据：[改法选项](HF_MIDPOINT_FINE_STAGE_OPTIONS_2026-11-06.md) §3 P1 与下面的用户决定。

## 必读（开工前全部读完）

- `docs/HF_MIDPOINT_FINE_STAGE_OPTIONS_2026-11-06.md`：P1 的依据、数据与误判代价。
- `docs/PLAN_HF_MIDPOINT_COARSE_2026-11-02.md` §1 与"未确认"第 1 条：被撤回的那条 v4 附加规则。
- `explore/src/explore/ts_ef.py` 的 `coarse_sign`、`stationary`、`candidate` 与二分循环。

## 用户决定（2026-11-06，原文）

问"P1 可以定下来吗？"，回答「可以」。Q2（精阶段噪声受限出口）未定。

## 1. 规则

在 v5（两阶段中点 + 精阶段 FIRE）上：粗阶段两侧都收敛后，`coarse_sign` 依次检查括区宽度、几何（`branch_rmsd_A`、`branch_step_A`）、双侧反力差门。三项都通过、且任一侧满足 `stationary`（完整力 ≤ `fmax_eV_A`、横向力 ≤ `transverse_fmax_eV_A`、约束残差）时，取满足的一侧中完整力最小者作为候选，HS 结束。其余情况与 v5 完全相同（包括三项检查任一不过时转精阶段）。

不变：候选判据数值、几何门、反力差门、B、粗阶段 BFGS、精阶段 FIRE 与 1e-3、200 步上限、括区端点与延续。

## 2. 开关与版本

- `EFSettings.midpoint_coarse_candidate: bool = False`。只允许在 v5 上开（要求两阶段规则且 `midpoint_fine_optimizer="fire"`），开启时 provenance 为 `HF-EF-bracket-v6`；不开时 v3 / v4 / v5 行为不变。
- 粗阶段出候选的那两帧记 `coarse_candidate = True`，`coarse_sign` 仍为 0（汇总脚本按 `coarse_sign` 统计提前收括区，不受影响）。
- D2 / `HFSettings` / manifest 不接，与 v4 / v5 相同，验收通过后再定。

## 3. 测试（CPU，解析台）

- 设置校验：不是 bool、没有两阶段规则、精阶段不是 FIRE 时报错。
- 粗阶段一侧已满足候选判据、门都通过：v6 在粗阶段出候选（origin 为 `*_coarse`，之后没有精阶段帧），候选满足候选判据；同一输入下 v5 仍转精阶段（现有测试）。
- 粗阶段一侧满足候选判据但反力差超门：v6 不在粗阶段出候选，转精阶段。
- 不开 P1 时现有测试全部不变。
- 变异检验（临时目录 + `-o pythonpath=`）：把候选检查挪到反力差门之前；忽略开关（总在粗阶段出候选）。各应使至少 1 个测试失败。

## 4. 验收（GPU，正式运行前报计划表）

照 11-04 v5 验收（`midpoint_fire_validate.py`）：22 个失败意图 + 15 个对照，v5 / v6 同机配对各跑完整 HS 一次；结局或候选不同（有无，或超过运行判重标准）的项两版各复核 2 次。统计：救回数、对照候选是否保留、候选几何身份（对称 RMSD，`candidate_identity.py` 同一判据）、后续卡点、forward 与墙钟。新增候选做 HC 另报。

## 5. 进展

- 2026-11-07：代码、测试（113 passed，两个变异均被测出）与冒烟完成（失败、对照各 1 项，v6 都在粗阶段出候选，与 v5 候选 ΔE ≤ 0.002 eV、对称 RMSD ≤ 0.0018 Å）。数字见 DONE_LOG 同日。
- 2026-11-07 验收（PBS yayoi40，作业 24581）：对照 15 项候选全部保留（ΔE ≤ 0.0039 eV、对称 RMSD ≤ 0.0049 Å）；失败组首轮出候选 15 → 16 项，复核中 `k0·2-4` 1/3 → 3/3、`k0·8-14` 0/3 → 2/3、`k2·6-8` 0/3 → 3/3；v6 候选 38 / 39 出自粗阶段；全部候选属已知 3 类，失败组都是新 E2 鞍点，新增事件 0；精阶段失败 8 → 1；forward 失败组 −18%、对照组 −8%。剩余失败：第一个中点之前 5、粗阶段 BFGS 2、精阶段 1。

## 未确认（Claude 推断）

- P1 只在 v5 上开放，是为了与验收对象一致、版本号不交叉；P1 本身不依赖精阶段用哪种优化器。
