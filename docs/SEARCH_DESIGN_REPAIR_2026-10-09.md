# PRRS 搜索设计修复 · 2026-10-09

本轮依据 b09 probe1 及随后用户提供的排序重算，选择方案 B，并修复有可重复工程证据的一般缺陷。原设计的 [§3.1](../mlp_native_perturbation_response_reaction_search.md) 将非键接近列为可能的成键探针，§8 要求先广泛探索再使用响应自适应。软运动的方向重叠本身不足以承担“哪个原子对会反应”的先验。

## 状态与证据边界

| 问题 | 本轮处理 | 证据及剩余边界 |
| --- | --- | --- |
| 排序方向与执行方向不同 | 保留用户已完成的解析 `delivered_direction` 修复；环内新交付也同步更新切向量 | 用户报告 b04 排名保持第一、旧打分重现 b09 24/24；本轮没有重跑该 POLAR-1 谱 |
| 投影残量冒充软模式 | 重新对角化去掉刚体与实际扭转运动后的 Hessian；SVD 按真实秩去掉依赖约束 | 小矩阵近纯扭转反例、投影算子残差、正交性及秩回归；不是 b09 的已发生原因 |
| 细化挤占方向覆盖 / exhaustive 漏算细化 | 每个源微态先完成所有选中方向的粗扫，再用余量细化；预算预检保留粗扫计数 | 三层预算分别恰好等于粗扫需求、剩余预算细化、粗扫/两阶段交界/细化中真实进程退出后的逐字节续跑回归 |
| 环内交付拖动其他共价键 | 无刚体分割的位移改为质量度量下的约束延拓，每步恢复非目标共价键；门限不变 | 苯环三族目标精度与保键、切向一致性、不相容约束拒绝；b09 原五个拒绝探针的纯几何重放全部过原门 |
| 双分子体系的软模成键先验 | **尚未建立有效的替代先验**。先用全候选粗扫取得响应基线，再比较固定预算的候选策略 | 新预算协议已准备，GPU 动力学尚未执行；不能声称 b09 已抓到反应 |

此前审查的 0.368 是“某个实际位移向量与旧端点梯度的重叠”，用户重算的 0.515 是“新执行方向与四个软模式的最大重叠”；两者不是同一观测量。用户测得 b09 四模投影保留 99.4–99.9%，排除了近纯扭转残量放大作为此次失败原因。

## B 的精确范围

`search.py::finish_scan` 复用原有窗口细化、样本整理、committor、鞍点精修和边界延拓逻辑。调度分成源微态内的两遍：

1. 按已有提案顺序执行所有 coarse 批次，继续记录三层预算截断原因。
2. 在同一源微态剩余预算内细化已有返回/离开窗口，再生成原来的扫描分析记录。

“先全部粗扫”指当前源微态的选中集合，不意味着整个动态网络所有后续微态都先粗扫。若 sampled 没选中成键方向，B 本身不会增加候选覆盖；若粗扫本身已经超预算，也仍会截断并留账。

恢复事务仍以微态开始前的 `resume_point.json` 为边界。任一遍被杀都回滚该微态的在途 trial / path / TS 并重放两遍，避免持久化不完整的待细化内存队列。原运行保持原源码身份；新源码不能续接旧指纹的运行。

## 环内交付的物理含义

`internal.py::constrained_displace` 在每个小目标增量上解质量加权最小位移，约束为目标坐标变化与所有非目标原有共价键长度。有限步后的非线性约束用几何迭代恢复。无模型能量或力参与，目标步长和单次原子校正都有限幅。

它改变的是扰动初态的交付路径，不是自由响应动力学：约束只用于交付，后续 MD 和 quench 保持原物理势能面。无法兼容的几何约束显式失败；clash、fragment、能量、力等门继续判定。此交付仍不能代表最小能量路径或协同反应坐标。

[b09_constrained_delivery.json](../benchmarks/carbene_ambimodal/results/2026-10-09/b09_constrained_delivery.json) 记录原五个 `collateral_bond` 拒绝探针的重放：非目标键变化约 2.2e−16–1.1e−15 Å、目标误差约 8.3e−17–4.4e−16，五个全部通过原门限。每个包括审计的交付约 3.6–3.8 秒。**只验证几何，未评估 POLAR-1 能量或运行 MD。**

## 回到设计稿：成键探索的下一份可执行协议

先不添加卡宾、烯碳或已知产物提示。第一步用现有 exhaustive 候选基绕开软模筛选，获得覆盖明确的响应资料：

- 输入、模型、电荷/自旋、seed、幅度、响应时窗、门限均沿用 b09 probe1。
- [config_coarse_coverage.json](../benchmarks/carbene_ambimodal/probe_b09/config_coarse_coverage.json)：432 个源候选、1408 个粗扫样本；`max_trials` / `chemical_trials_per_node` / `conformer_trials_per_microstate` 同为 1408。
- `refinement_steps` 保留 3，但方案 B 在粗扫用满预算后没有 trial 可细化。这份协议验证粗扫覆盖；另开含余量的运行验证细化。
- 使用新的 `runs/b09_coarse_coverage`；执行机器记录当前源码、模型与输入身份。本轮没有开跑。

```sh
PYTHONPATH=.:src python -m prrs search benchmarks/carbene_ambimodal/probe_b09/input.extxyz \
  --output runs/b09_coarse_coverage \
  --config benchmarks/carbene_ambimodal/probe_b09/config_coarse_coverage.json \
  --calculator benchmarks.b04_sn2_chloride.calculator:polar1
```

验收先看源微态：432 个方向 / 1408 个 coarse 是否交付、是否有拒绝或预算丢失、关键接近方向是否执行，再读化学变化。不能由 `status=completed` 推断完整覆盖。深微态、有限时窗和原始起点偏向 TS-2OMe 的限制保持明确。

得到响应基线后，在相同物理 trial 预算上比较旧谱排序、逐族无放回均匀采样、按原子对/运动覆盖分层的首轮策略。比较的是各类实际交付覆盖、有效响应、化学新节点、确认连接及实际成本；已知 TS 坐标只做事后审计。对强相关的跨片段平移要同时记“原子对覆盖”与“实际运动覆盖”，不能误把多个近共线对当成独立探索，也不能仅因运动相近就合并不同成键位点。

不通过这些对照就不把新的排序启发式变成默认，不宣称代码修复建立了双分子反应发现能力。

## 本轮验证

固定源码指纹：`d6a6986e147aa853574597dd10f463c2c286162491cf62949b67f960aba70128`。

- `make check` 完成：fmt-check / lint 通过；**592 passed, 2 skipped**；历史 digest **332 files / 2136 frames / 4 event pairs 零差异**。digest 验证身份层与历史帧，不验证新搜索编排的行为等价。
- 定向测试 **64 passed**，包括粗扫预算、两遍之间及细化内真实进程退出后的逐字节恢复、受限 Hessian、保键交付与原有内坐标/解析打分。
- 变异检验使用 `/tmp` 中的独立源码副本（源码副本没有拷进仓库）。恢复边扫边细化 → 3 个预算回归失败；恢复投影残量归一化 → 3 个谱回归失败；恢复环内无约束交付 → 3 个保键回归失败。记录 [mechanism-mutations.json](../benchmarks/carbene_ambimodal/results/2026-10-09/search_repair/mechanism-mutations.json)、[mutations/](../benchmarks/carbene_ambimodal/results/2026-10-09/search_repair/mutations/)。
- 新全粗扫配置通过现有 exhaustive 预检：1408 次需求与三个上限一致。
- PARS 文件接口检查：新合成运行的 12 条账本均通过 result/ledger 校验；默认过滤不收 C2 解离，显式包含 `dissociation` 后导入 `t000002` / `t000003` 两条路径。没有改变 PARS 默认科学筛选规则。记录 [coarse-first-pars-check.json](../benchmarks/carbene_ambimodal/results/2026-10-09/search_repair/coarse-first-pars-check.json)。

中途曾在一轮续跑测试尚未结束时继续编辑源码，触发了正确的实现指纹拒绝；该轮不能作为验收。上述全量与定向结果在源码固定后取得。

全量日志：[mechanisms-full-check.log](../benchmarks/carbene_ambimodal/results/2026-10-09/search_repair/mechanisms-full-check.log)；定向日志：[mechanisms-final-focused.log](../benchmarks/carbene_ambimodal/results/2026-10-09/search_repair/mechanisms-final-focused.log)。全部验证记录（原在 `/tmp`，10-09 拷入仓库）见 [search_repair/](../benchmarks/carbene_ambimodal/results/2026-10-09/search_repair/README.md)。没有运行新的 POLAR-1 GPU 搜索，也没有把已有 b09 运行改成新协议结果。
