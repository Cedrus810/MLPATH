# 施工计划：同键图的脊最后派发（R3a-2，D2-HF-v1.2，TODO 第 9 条 R3 / F4）

日期：2026-11-07。依据：[R3 改法选项](HF_HS_NO_MAXIMUM_OPTIONS_2026-11-07.md) §3 R3a 与下面的用户决定。

## 必读（开工前全部读完）

- `docs/HF_HS_NO_MAXIMUM_OPTIONS_2026-11-07.md`：测量、R3a 的依据与误判代价。
- `docs/PLAN_START_CLUSTER_2026-10-28.md`：D2-HF-v1.1 的两轮派发；本规则在它之上再加一组。
- `explore/src/explore/d2.py` 的 `_promote`、`_start_clusters`、`_hf_tasks`、`_hs_task`、`_hf_report`。

## 用户决定（2026-11-07，原文）

问"R3a 用 R3a-2（降优先级、不拒跑），可以吗？"，回答「可以」。

## 1. 规则

在 D2-HF-v1.1 上：

- 被升级的项中，`_classify` 为 `ridge_suspected`、且 L1 扫描峰点 k 与 k+1 的 `edges`（扫描自带，与 `jump` 触发同一定义）相同的，记为同键图的脊（`_same_graph_ridge`）。
- 派发顺序：第一轮 → 第二轮（`HS_deferred:same_start`）→ 同键图的脊。最后这组只用前两轮剩下的额度，额度不够的记 `HS_deferred:same_graph_ridge`，不拒跑、不记成收敛或无反应。
- 起点分簇时，同键图的脊排在最后参与，所以它们不会成为簇首、把成键改变的脊挤到第二轮。
- release（L3）照旧对所有升级项执行，顺序不变。
- HS 结算的 payload 加 `same_graph_ridge`（bool）标签；报告加 `same_graph_ridge: {settled, deferred}`。

不变：`read_scan` 的触发条件本身、`_classify` 的类别、升级排序、HS / HC 的算法与上限。

## 2. 开关与版本

`init_run(..., same_graph_ridge_last=True)`，需要 `hf` 与 `start_cluster`，冻结为 `D2-HF-v1.2`；manifest 写 `"same_graph_ridge_last": true`（v1.1 及以前的 manifest 不写这个键，保持原样）。旧运行不追认。HS 的 E/F 设置（v6 / v7）没有接进 D2（9d，另议）。

## 3. 测试（已施工）

- `_same_graph_ridge`：峰前后成键相同 → True；改变 → False；非脊类别 → False。
- 冻结：缺 `start_cluster` 或非 bool 报错；v1.1 manifest 无此键；v1.2 规则版本与键。
- 解析台（TwoSaddle，把第一个升级项人为标成同键图的脊）：额度充足时结局与参照相同、该项最后派发、不再当簇首、外部计数 = 账本；额度只够一个 HS 时该项记 `HS_deferred:same_graph_ridge`。
- 变异检验：起点分簇不重排、同键图的脊不放到最后、暂缓原因写错，分别使 2 / 2 / 1 个新测试失败。
- `test_ts_ef test_ts_branch test_ts_connect test_ts_warm_start test_d2_hf test_d2 test_d2_resume test_probe` 147 passed；全包 `make check` 未跑。

## 4. 离线重放（已做）与整跑验收（待定）

离线重放（`benchmarks/explore_d2/same_graph_ridge_replay.py`，八块 b05 HF，用新代码的 `_classify` 与 `_same_graph_ridge`）：

- 2 634 个升级项中，1 663 个被排到最后（ridge 类共 1 932 个）。它们原来花了 HS 总费用的 87.3%（1.64M / 1.88M forward）。
- 这些项里候选 0 个（全部 52 个候选都不受影响）。
- 22 个新 E2 意图与 15 个对照意图都不受影响。对照中 `k5·3-14` 的标签在 v4 sA 命中一次，但那是同一缩写下的另一个起始姿态（以 `hysteresis` 结束），出候选的那个不受影响。
- HC 状态没有读入（重放用的 `load` 只保留 L1 / L2 / HS），以候选代替。
- v4 两块中，各有 32–34 个非同键图的项原来因额度用完从没跑 HS（`branch_bracketed` 18–21、`ridge_suspected` 3–4、`audit` 8–12）。v1.2 下它们排在同键图的脊之前。

整跑验收：b05 两个种子（PBS），与 10-29 yayoi32 v4 整跑对比事件、候选、forward、各组跑到的项数。建议与 9d（v6 / v7 接入 D2）合并做一次，正式运行前报计划表。

## 未确认（Claude 推断）

- 在 b05 的额度（refine 240k）下，非同键图的项跑完后额度还剩很多，同键图的脊照样会把它花掉。所以 v1.2 本身在 b05 上省不了多少 forward，主要作用是让原来跑不到的 32–34 项先跑。它们能出什么，没有数据。
- 只有 b05 一个体系的证据（见选项文档"未确认"）。
