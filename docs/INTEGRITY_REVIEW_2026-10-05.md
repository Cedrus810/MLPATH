# 新确认的完整性缺口 · 2026-10-05

此页保留 `cf9c90d` 基线的复现证据，原报告不改写。
三项缺口已修复，见 [修复与验收记录](INTEGRITY_FIX_2026-10-05.md)。
脚本：`benchmarks/pars_integrity/audit.py`；冻结证据：
`benchmarks/pars_integrity/results/2026-10-05/report.json`。
本次只验证解析外场和合成 PRRS 文件，未证明任何历史分子结果已受影响。

## P1：全体系 RMSD 稀释局部差异，第三盆地被误认成端点

`pars/refinement/ase_neb.py` 的下降落点匹配只限制平均每原子 RMSD ≤ 0.05 Å。
499 个谐振固定旁观原子加 1 个活性原子；活性坐标势能为
`E=x²(x²−1)²`，极小值为 −1、0、+1，两个一级鞍点为 ±1/√3。
端点设 −1/+1；它们的整体 RMSD 为 `2/√500≈0.08944 Å`，不会被端点塌缩检查拒绝。
中间极小值到任一端点的整体 RMSD 仅 `1/√500≈0.04472 Å`，尽管关键原子相差 1 Å。

真实 `AseNebBackend.refine` 返回 `ci_neb_converged` 和约 0.148148 eV 能垒；
实际下降落点为中间极小值与一个外侧极小值。
粗、细下降均稳定到同一对真实落点，并把中间极小值误认成另一个外侧端点。
因此上一轮的下降分辨率复核不能解决身份判据错误。
解析外场使用 `refinement_remove_rigid_motion=False`，并非分子验证。

修复验收：身份判据不能因加入不参与反应的原子而变宽；
下降落点匹配与源/产物见证匹配应共享明确的映射和局部偏差约束。
保持 PRRS 身份权威，PARS 不自行重建化学节点身份。
本反例应拒绝报告端点间一级鞍点连接/能垒，正常同盆地精修仍应通过。

## P1：结果 JSON 与轨迹没有完整性绑定，断写后续跑跳过坏结果

`pars/cli.py` 的续跑只比较 JSON `conditions`，并不检查配套 extxyz 存在或完整。
轨迹直接覆盖原文件；`--redo` 写轨迹失败时，旧 JSON 提交标记仍在。
相同条件下再次调用，会输出 `0 refined, 1 already done` 并返回 0。

两条独立复现：删除已经记录的 extxyz；注入写轨迹中途异常并留下截断文本。
两次正常续跑均跳过。此处初次精修记录为有意的 index-three `failed`，
没有将 CLI 返回 0 等同于科学成功；缺陷是把缺失/损坏的交付工件当成已完成记录。

修复验收：轨迹写入临时文件并原子替换；JSON 提交标记绑定内容哈希、
候选身份和帧数量，恢复与聚合都核验。
两文件切换中途崩溃必须保持一组可核验旧结果，或拒绝/重做不完整新结果；
同条件正常恢复仍不得重复计算。

## P2：PRRS 观察记录与轨迹仅检查行数，重排后静默错配

`pars/adapters/prrs_run.py` 文档明确 `row.frame == trajectory index`，
实际读取只检查长度，然后按文件顺序把 phase/energy/edges/uncertainty 配到构型。
保持 8 行但倒序 `[7,...,0]`，读取成功：源构型配上 `quenched` 阶段和 −0.8 eV，
正确值应是 `source` 与 −1.0 eV。这会污染反应段裁切和后续目标计算。

修复验收：在生成 DiscoveryTrace 前验证 frame 索引连续唯一且顺序一致，
核对存在的逐帧阶段信息和原子映射；错误工件应明确拒绝。
还应审查 result 与 attempts 分类一致性、端点原子映射；后两项本次未单独复现。
当前属于损坏输入的校验缺口，排在前两项之后。

## 重跑与顺序

```bash
OPENBLAS_NUM_THREADS=1 /home/ruigengji/miniforge3/envs/openmm_dev/bin/python \
  benchmarks/pars_integrity/audit.py --output /tmp/pars-integrity-report.json
```

顺序：端点身份 → 结果事务/完整性 → 导入对齐。
此前两项 factory/model 行为参数的小缺口仍见 `docs/DEFERRED_FIXES_2026-10-04.md`，
本轮未修。真实分子、独立 QM、完整网络预算与优越性证据仍未补齐。
