# explore — 发现层 v1 与评估账本

按 `docs/DISCOVERY_LAYER_DESIGN_2026-10-13.md` §1.1 与 `docs/EVALUATION_PROTOCOL_DESIGN_2026-10-13.md` v1 范围施工。
单向依赖 `prrs`（`../src`）；`prrs`、`pars` 不得 import 本包。测试：`make test`（需 `openmm_dev`）。

约束扫描默认复用相邻已收敛点的曲率作为下一点的初值，并重置优化历史；
`ScanSettings(reuse_curvature=False, ...)` 保留旧基线。

`explore.parallel.scan_many(jobs, backend, max_active=4)` 可让多个独立扫描共享
一个常驻的 `MACEBatchEvaluator(calc)`，模型只建立一次。各 `ScanJob` 必须有
独立的 `Ledger`；一批 m 个结构计 m 次求值。返回顺序与输入一致，失败单独保留。
该接口只执行扫描，不包含发现调度、恢复或 C2 认证；认证继续使用原 calculator。
执行合同和验证范围见 [扫描吞吐记录](../docs/SCAN_THROUGHPUT_2026-10-23.md)。

`explore.d2`：方案 A 的 v1 编排（[M5 计划](../docs/PLAN_EXPLORE_M5_2026-10-23.md)）。`init_run(run_dir, poses, specs_by_pose,
settings, config, identity)` 冻结人口表；`run(run_dir, Backends(search_batch, certify_calculator))` 按 L1 → L2 → 释放/认证
执行并写 `report.json`。重跑同一条命令即续跑（已结算任务跳过，被打断的按账本或上限计费后用新 attempt 重做）。
**只在写出该 run 的主机上续跑**：NFS 属性缓存下别的主机可能看到旧 journal。`D2Settings` 没有默认值。

`explore.ts_ef.refine_ts_ef`：独立的连续约束分支 E/F 反力零点定位入口，
返回可序列化的 `EFResult`（驻点候选、完整力与残差、括区、方向、历史、成本）。
它与旧 C2/D2 链分开，尚未接入 D2 或恢复；候选性质与连接仍待独立核验。
接口、解析验证与限制见 [HF-2/3 实现记录](../docs/HF_EF_LOCATOR_2026-10-25.md)。
