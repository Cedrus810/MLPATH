# HF-2/3：独立 E/F 驻点候选定位原语

日期：2026-10-25。源码：[`explore/ts_ef.py`](../explore/src/explore/ts_ef.py)。

## b05 复核结论

F6 单变量结果不支持修复完成：反应性 E2 从 1 变 0，新增命中主要为重复构象鞍点。
下一步应验证新的定位内核。完整 Hessian 依赖与旧更新策略是已确认的实现问题；
但本轮结果不能单独证明 E2 丢失的具体物理原因，也未确认 −18880.50 eV 鞍点的运动类型。

只读脚本 [`hf_f1_audit.py`](../benchmarks/explore_d2/hf_f1_audit.py) 生成
[`hf_f1_audit.json`](../benchmarks/explore_d2/results/2026-10-25/hf_f1_audit.json)：

| 诊断 | sA | sB |
| --- | --- | --- |
| lost_target_mode 总数 | 122 | 111 |
| 第 0 步退出、初始重叠记录为空 | 30 | 29 |
| 成功初选后才丢模 | 92 | 82 |

空值原因已复现：`follow_min_mode` 返回字典先写刚计算的 `target_overlap` 和
`selected_mode_index`，再展开仍是初始值的 `**bounds`，覆盖了两个字段。
本轮只调换展开顺序，搜索数值与阈值不变；历史 journal 不修改、不补造缺失值。

另需区分：第一次选择对照扫描切向，之后对照的是上一次选中的本征向量。
所以 `lost_target_mode` 不是“某个本征模与最初扫描方向不再重叠”的直接观测，
也不能凭它区分真实方向漂移、模混合、刷新间隔或谱的数值问题。
fp32 sA 原反应性 E2 对应的意图在 v1.1 第 **80** 步停止，初始重叠 **0.7030**、
末次重叠 **0.01818**；25–35 步是全部失败项的中位数。

## 本轮实现

`refine_ts_ef(atoms, intent, bracket, factory, config, ledger, task, ...)`：

1. 输入同一原子顺序的两个 `ScanPoint`，按实际坐标目标 q 排序，重新约束松弛并读取 E/F；
   不信任保存的收敛标记或反力，局部峰可以低于扫描起点。
2. 横向残余力与约束残差合格后，要求反力按 q 递增为正→负，且超过声明的反力容差。
3. 每次二分 q，分别从左右已松弛结构延续到同一中点；不插值两套几何。
   两次延续的对齐 RMSD、步间 RMSD 和反力差超限时返回分支未决。
4. 双侧分支一致性检查通过后，分别检查完整无偏力、横向力和约束残差；任一结构达标，
   即从达标者中选完整力最小的返回 `converged_candidate`。`reaction_force_tol` 只用于括区
   符号分辨，不再作为候选的独立反力绝对值门，也不要求两侧同时驻点收敛。
   只缩窄区间不算收敛；反力噪声底无法消除时返回 `force_resolution_unresolved`。
5. 分别保留候选 E/F、两侧括区、所有已评估帧、局部切向和相对初始方向漂移、
   设置、坐标、原子序号、电子态、停止原因及成本。`to_dict()` 可直接写 JSON。

只调用 E/F 和已有可靠性守卫。BFGS 用轨迹 secant 学习近似曲率，未调用模型完整 Hessian、
`curvature_spectrum` 或逐坐标差分重建；每次实际模型求值仍由 counted factory 计费。
能量上限逐个已评估结构检查，越界保留该次费用与结构并停止。

独立结果为 `EFResult`，其 `property_status` 和 `connection_status` 初始均为 `not_run`。
这不是 C2，也不自动声明一级鞍点。旧 `certify` 与 D2-v1.1 继续走原链；本轮没有默认替换。

使用约定与 `constrained_scan` 一致：

```python
result = refine_ts_ef(
    atoms, intent, (before, after), counted_factory(factory, ledger),
    config, ledger, task="ts-local-0", settings=EFSettings(), ceiling_eV=ceiling,
)
record = result.to_dict()  # 包含未决与截断记录；由调用方持久化
```

## 已验证与限制

解析台采用平移／旋转不变的耦合双键距离势，驻点可独立求得。
测试覆盖完整谱与模型 Hessian 禁用、完整力／几何误差、刚体变换、反向扫描、
负坐标系数、极小、平坦、同图双分支、力噪声底、步数／预算／能量截断和坐标退化。
外部 calculator 计数与账本一致；结果严格 JSON 序列化通过。

最终 `make -C explore check`：格式、lint、**177 passed**（新增 19 项）。
PRRS 相关 **96 passed**；修改前后历史摘要检查均为 332 文件、2136 帧、4 事件对不变。
检查日志与源码指纹见 [`hf_ef_checks.json`](../benchmarks/explore_d2/results/2026-10-25/hf_ef_checks.json)。

分支阈值是显式数值诊断，不能证明分支连续或横向稳定；两侧都跳到同一个分支仍可能漏检。
当前仅支持非周期分子。两侧反力不能组成局部最大括区时明确未决，不回退到旧精修器。
原子映射由输入顺序保证，尚无自动重排验证。

本轮尚未接入 D2 的调度、原子落盘及 kill→resume，也没有运行真实模型或 QM。
下一关：选 b04 及 b05 原 E2 的已保存扫描括区做有界局部对照，分别保存搜索与独立性质／
连接结果；确认新内核的真实 E/F 尺度和方向后再接 HF-4，不直接重复整批 b05。

## 补充（2026-10-25，Claude）：真实模型对照与分支力门

有界对照与反力抖动实测见 `DONE_LOG.md` 同日三节。`EFSettings.reaction_force_sigma`（新增，可选）：调用方按体系实测的同支反力 1σ，给定时 `branch_force_disagreement` 门为 `3·√2·σ`；不给定保持原 `2 × reaction_force_tol`。b04 fp32 实测 σ ≈ 0.0034 eV/Å，`reaction_force_tol` 0.003 约为 1σ。

## 候选判据修复（2026-10-25）

`HF-EF-bracket-v2` 将双侧分支一致性与单结构驻点收敛分开：先检查双侧几何、步间几何
与反力差，再按第 4 项选候选。完整笛卡尔力已包含反应坐标分量；另外要求
`|dE/dq| <= reaction_force_tol` 会引入依赖 q 归一化的额外限制。
实测 σ 仍只决定双侧反力差门，不放宽完整力或横向力容差。实际设置、分支力门与
`candidate_rule` 写入 provenance，历史 v1 结果不覆盖。

只读检查历史 `hf_ef_bounded_sigma/b04.json`：第 7 个中点 q = 2.2798038977 的反力差
0.00835765 < 0.01442498 eV/Å，双侧对齐 RMSD 0.00036956 Å；右延续完整力
0.00580621、横向力 0.00257815 eV/Å、约束残差 0，按新规则合格。
左延续完整力 0.01407712 eV/Å，不合格，不妨碍返回已通过分支检查的右侧结构。
这是对保存帧的离线判据复核，**不是模型重跑或独立 TS 认证**。

第 8 个中点的反力差 0.04318588 eV/Å 仍超过原分支门，新规则不会绕过它；
此前测得的 σ 来自另一 q 附近，不能直接解释此差异。当前环境 `nvidia-smi` 无法连接驱动，
未做新的 GPU 测量。后续需在同一坐标系固定这两套几何，比较 fp32／fp64 的反力，
并对照各自进一步横向松弛后的反力，区分求值精度和残余几何的影响。

新增解析回归覆盖两种延续顺序下的单侧达标、达标侧不能绕过反力分歧、σ 不放宽完整力。
`test_ts_ef.py` **24 passed**；`make -C explore check` 格式、lint、全包 **182 passed**。
离线逐中点复核、检查结果与源码／历史输入 SHA256 见
[`hf_ef_candidate_v2_checks.json`](../benchmarks/explore_d2/results/2026-10-25/hf_ef_candidate_v2_checks.json)。

## 松弛曲率复用与切向预测（2026-10-26）

`HF-EF-bracket-v3` 与 `HF-EF-branches-v2` 增加两个 `EFSettings` 开关，默认均开启：
`reuse_curvature`、`predict_tangent`。两者均设为 `False` 可恢复每个目标点冷启动 BFGS
的数值策略，用于同预算对照。完整力、横向力、反力差及几何分支判据不变。

每条分支独占运行期 `_RelaxState`，保留最近两个已接受点与 E/F 学得的 BFGS 曲率。
换 q 时仍创建新优化器，只传 `P H P + alpha n nᵀ` 到新约束切空间，不传跨约束
secant 历史。零步收敛保留已有曲率。几何坐标系不匹配时丢弃提示。
预测步取该支两个已接受点的对齐割线外推，每个原子位移限制在 `max_step_A` 内，
再投影到目标约束；不会混合两支几何。

每个试步使用状态副本，只有分支检查接受后才写回。跳变拒绝、步长折半、能量／预算／
松弛失败均不污染已接受状态。二分两侧各用独立副本，延续转交二分时保留对应端点的曲率。
运行状态不序列化为结果或恢复证书；设置与算法版本写入 provenance，逐帧记录
`curvature_reused` 和 `predictor_step_A`。新状态未接入 D2 调度或 kill→resume。

新增解析测试覆盖真实 E/F 计费下降、解析驻点一致、真实折叠滞回保留、拒绝试步回滚、
预算／能量失败、零步曲率保留与新约束法向投影；相关测试 **39 passed**。
真实模型有界 A/B 入口为
[`hf_warm_start_bounded.py`](../benchmarks/explore_d2/hf_warm_start_bounded.py)，
各次搜索单独限制为 3000 forward，独立核验另记账。

本机 2080 Ti 的同轮 A/B：b04 冷启动 413 forward / 37.9 s → 合用 215 / 18.9 s；
b05 `6-9` 859 / 76.1 s → 284 / 25.5 s。优化后 5/5 出候选；b04 与 b05 `6-9` 的新候选
独立验为一阶、两侧连通 SN2／E2，另三项本轮只复核搜索。全包 `make -C explore check`
**197 passed**，格式／lint 通过。逐例性能、核验、账本与源代码快照入口：
[`SUMMARY.json`](../benchmarks/explore_d2/results/2026-10-26/hf_warm_start/SUMMARY.json)。
