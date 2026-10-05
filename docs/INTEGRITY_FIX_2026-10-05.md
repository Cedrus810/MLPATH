# 端点身份、结果完整性与导入对齐修复 · 2026-10-05

修复 [审查记录](INTEGRITY_REVIEW_2026-10-05.md) 中的两项 P1 和一项 P2。
基线证据保留；修复后记录为
`benchmarks/pars_integrity/results/2026-10-05/after-fix.json`。
记录中的 inspected_commit 是施工前 HEAD，实际执行源码由
`pars_implementation_sha256=57fcad4d8785761bcdf862337e972fc9a2688cb87b735b8527adf1e30dc6b3d0` 固定。

## 端点身份

下降落点及独立源/产物盆地见证，都要求整体 RMSD **和每原子最大位移**
不超过原有身份阈值 0.05 Å；使用同一映射与可选 Kabsch 对齐。
共享距离计算，不在 PARS 中重建 PRRS 的化学节点身份。
端点塌缩也同时检查局部最大位移，避免把稀释后的不同端点判成同一极小值。

500 原子完整 CI-NEB 反例现在返回 `failed`，`barrier_eV=null`。
落点仍为中间极小值和外侧极小值；即便 RMSD 仅 0.0447 Å，约 1 Å 的局部偏差
明确拒绝错误端点匹配。正确的相邻盆地连接仍能通过。
回归同时覆盖单原子对照、大体系正确连接及独立盆地见证错配。
这是数值几何验收，不能替代化学身份、连续 IRC 或独立 QM 验证。

## 结果工件

extxyz 在目标目录写临时文件，完整写入、fsync 后再原子替换；失败清理临时文件。
随后提交的候选 JSON 保存 `geometry_artifact`：文件名、SHA-256 和帧数。
续跑及聚合在进入计算/写汇总之前核验既有工件；还检查候选身份、帧编号及状态。
聚合后的每条记录保留同一绑定，消费汇总时也应核对对应几何工件。

- 删除/截断工件、错误帧数/候选身份、没有绑定的旧记录：拒绝静默跳过。
- 写轨迹中途失败：旧 JSON 与旧轨迹均完整保留，可正常恢复。
- 替换轨迹后、提交 JSON 前中断：旧哈希不匹配，恢复明确拒绝。
- `--only` 不能绕过另一个损坏结果，也不能把它混进汇总。
- `--redo` 按被指定的记录文件名替换，损坏的内部 candidate_id 不能越过作用域。
- calculator 失败且没有 band 时，保存 `failed`、空工件和零帧绑定，仍可续跑。

旧记录不会自动升级成新验收结果。需要用原 selection、显式 calculator factory
执行 `pars refine ... --redo` 重做；源 run 不改写。
JSON 与 extxyz 仍是两个文件，保证是保留可核验旧组或检出不完整新组，
并不声称两文件具有单次联合原子提交。

## PRRS 导入

逐行要求 frame 为整数，且精确等于从零开始的轨迹索引；拒绝倒序、重复、
缺失索引和布尔索引。轨迹存在 phase 时必须与观察记录一致。
所有帧及端点必须保持原子数和元素序列。
result 与 attempts 的 attempt_id/status/classification 冲突也拒绝。
同元素原子的置换不能仅靠元素序列检出；生产者仍须保留映射顺序。

## 验收

- `OPENBLAS_NUM_THREADS=1 make -C pars check`：格式、lint、模块边界通过，190 项测试通过。
- 根独立鞍点参考/硬预算回归：6 项通过。
- 独立复制 PARS 到临时目录运行：189 项通过，1 项历史 run smoke 因缺少 run 跳过。
- 本轮三项实际复现脚本重跑：错误连接与倒序输入被拒绝，断写保护旧结果。
- 同一冻结双鞍点池、选择计划和每臂 6000 次力预算重跑：PARS 与 channel-pruning
  均覆盖 2/2，各使用 5464 次力计算；新增几何检查没有额外 force evaluation。

双鞍点记录另存于 `benchmarks/pars_two_saddles/results/2026-10-05-integrity/`。
原基准结果与本轮所有历史源 run 均保留。
该解析基准不证明一般优越性或真实分子/QM 有效性。
此前 factory/model 行为参数两项小缺口仍按既定顺序保留于
`docs/DEFERRED_FIXES_2026-10-04.md`。
