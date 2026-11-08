# 施工计划：refine 端点出候选（R1，HF-EF-bracket-v7，TODO 第 9 条 9c''）

日期：2026-11-07。依据：DONE_LOG 2026-11-07"第一个中点之前的 HS 停止普查"与下面的用户决定。

## 必读（开工前全部读完）

- `docs/DONE_LOG.md` 2026-11-07 普查一节：`no_resolved_maximum_bracket` 3 项的数字（端点帧完整力 0.0028–0.0030、横向力 ≤ 0.0029，能量在新 E2 鞍点那一段）。
- `docs/PLAN_HF_MIDPOINT_COARSE_CANDIDATE_2026-11-06.md`：v6（P1），R1 建在它上面。
- `explore/src/explore/ts_ef.py` 的 `refine_ts_ef`：两端重松弛后的 `no_resolved_maximum_bracket` 检查、`stationary`、`candidate`。

## 用户决定（2026-11-07，原文）

问"R1 可以定吗？"（R1：refine 端点已满足候选判据就直接出候选），回答「可以」。

## 1. 规则

`refine_ts_ef` 两端重松弛后，现行：左端反力 > tol 且右端 < −tol 才进入二分，否则 `no_resolved_maximum_bracket`。R1 只改"否则"这一支：

- 右端满足 `stationary` 且左端反力 > tol → 右端出候选；
- 左端满足 `stationary` 且右端反力 < −tol → 左端出候选；
- 两者都成立时取完整力小的一端；
- 其余照旧 `no_resolved_maximum_bracket`。

已分辨的括区（两端都过 tol）照旧二分，不在端点出候选。候选判据数值、二分、v6 其余规则都不变。

## 2. 开关与版本

`EFSettings.endpoint_candidate: bool = False`，要求 v6（`midpoint_coarse_candidate=True`），开启时 provenance `HF-EF-bracket-v7`。端点候选帧记 `endpoint_candidate = True`。D2 不接。

## 3. 测试（CPU，解析台）

- 设置校验：不是 bool、没有 v6 时报错。
- 极大剖面，右端恰在 TS：v6 为 `no_resolved_maximum_bracket`，v7 由右端出候选（origin `right`），满足候选判据，没有中点帧。
- 极小剖面，右端恰在极小（左端反力 < −tol）：v7 仍为 `no_resolved_maximum_bracket`。
- 不开时现有测试不变。
- 变异检验：去掉"另一端须指向最大值"的条件；忽略开关。各应使至少 1 个测试失败。

## 4. 冒烟与验收

- 冒烟（定向）：从 10-29 v4 sB journal 取 3 个 `no_resolved_maximum_bracket` 项的延续括区（`maximum_bracketed` 一侧的最后两点），v6 / v7 各跑一次 `refine_ts_ef`（冷启动），看 v7 是否由端点出候选、候选能量与几何身份。
- 回归验收（正式运行前报计划表）：22 个失败意图 + 15 个对照，v6 / v7 同机配对，照 11-07 v6 验收做法；看对照候选是否保留、端点候选出现几次、几何身份。`no_resolved_maximum_bracket` 本身少见（v6 验收 32 次中 3 次），配对验收里可能碰不到，所以定向冒烟是主要证据。

## 5. 进展

- 2026-11-07：代码、测试（118 passed，两个变异均被测出）、定向冒烟完成：3 项冷启动重跑后 1 项仍为未分辨括区，v7 由右端出候选（新 E2 鞍点，对称 RMSD 0.024 Å）；另 2 项括区已分辨，R1 不起作用。配对回归验收未跑。数字见 DONE_LOG 同日。

## 未确认（Claude 推断）

- "另一端的反力须指向最大值"是我加的限制，用户原文没有。理由：没有这一条，在 q 方向上是极小的完全松弛结构（整体力也为零）也会被当成候选。普查中的 3 项都满足这一条。
- R1 只在 v6 上开放，是为了版本号不交叉；R1 本身不依赖 P1。
- 端点候选只有一次松弛，没有中点那种两侧延续的分支一致性检查；端点与另一端的连续性由延续的跳变检查保证（两者是同一侧相邻的两个延续点）。
