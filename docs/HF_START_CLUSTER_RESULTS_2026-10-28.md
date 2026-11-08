# D2-HF-v1.1 两轮派发：b05 正式运行结果（2026-10-28）

结论：同额度下覆盖了更多不同起点，并保住旧的两类 E/F 连通事件；未找到新事件，也未减少本次总 forward。sB 的第二类事件由第二轮补回，证明不能把“每簇先跑一次”变成永久删除其余成员。出现 11 项括区接口异常，完整工程验收仍有缺口。

## 输入与核查范围

- 旧：`benchmarks/explore_d2/results/2026-10-26/b05_hf_s{A,B}/`，`D2-HF-v1`，写入节点 yayoi28。
- 新：`benchmarks/explore_d2/results/2026-10-28/b05_hf_v11_s{A,B}/`，`D2-HF-v1.1`，写入节点 yayoi04。
- 四份最终 report 均存在。逐项核对 journal 结算费用、ledger 的 close forward、report 两账户，全部一致；held 为 0，unreconciled 为空。没有调用模型，也没有启动新运行。
- 同种子的 settings、HFSettings、SearchConfig 与 L0 元数据相同；新 manifest 增加 start_cluster。实际 L1/L2 分类已有少量差别，因此不是相同轨迹上的纯顺序重放，也没有运行时源码哈希证明只有调度代码不同。
- 本次读取的是 kasuga01 上可见的 NFS 落盘文件。尝试直接连接写入节点时，yayoi04.cluster.local 的 DNS 解析失败；没有完成写入节点实时复核，不能据此断言远程进程退出码。下面是这些文件对应的结果。

证据与分析脚本快照：[`hf_v11_comparison/`](../benchmarks/explore_d2/results/2026-10-28/hf_v11_comparison/)。入口 `summary.json`；`audit.json` 保存逐项费用/状态，`event_match.json` 保存事件几何匹配，`brackets.json` 保存异常括区。summary 记录输入文件大小/时间及 manifest/report 的 SHA256；未给大体积 journal 补造摘要。当前分析代码哈希仅标识重放代码，不代表远程运行源码证明。

## 同额度的覆盖与产出

| 指标（旧 → 新） | sA | sB |
| --- | --- | --- |
| refine forward | 237127 → 237635 | 238176 → 237193 |
| HS 已结算项数 | 236 → 290 | 240 → 275 |
| HS 有实际 forward 的项数 | 236 → 283 | 240 → 271 |
| 有实际 forward 的起点簇 | 141 → 196 | 140 → 201 |
| 新规则第一轮派发簇数 | 203 | 205 |
| 新规则第二轮实际派发项数 | 87 | 70 |
| 新规则第二轮暂缓项数 | 42 | 55 |
| 对称合并后的 E/F 连通事件 | 2 → 2 | 2 → 2 |
| journal 首末事件间隔 / h | 8.53 → 9.27 | 8.81 → 9.71 |

旧起点簇覆盖数来自已保存的 `repeat_response_census/replay_start_clusters.json`（同一 r=0.014 Å、m=1 规则）；新起点簇数直接核对 HS payload 的 leader/rank。新派发簇中 7/4 项在首次 forward 前异常，所以不能把 203/205 全部算成实际完成搜索的起点。排除这些项后，覆盖分别增加约 39%/44%。这仍是两个各自运行的描述性比较。

这次 refine 几乎用满预算，第二轮继续消耗余额；此前离线重放的“省 34–62%”不能当成这次实际总 forward 的节省。墙钟也没有变快；写入节点不同，未做硬件/负载控制，不能把耗时差归因于调度。

## 事件匹配与第二轮的必要性

按姿态图自同构群（1296 个置换）对旧 connected 几何与新 connected 几何做对称 RMSD 匹配，同时检查 chemical_pair 和能量：旧事件全部能匹配到新事件，最大 RMSD 约 0.00585 Å、最大能量差 0.001953125 eV，均在原规则的 0.15 Å / 0.05 eV 内。两类之间约 0.57 Å。保住的是两类鞍点几何，化学反应对仍是同一种 E2；没有新增反应类型。

**新 sB 的反例：** `k3 / distance:10-13` 是某簇首，其 HS 为 converged_candidate，但 HC 为 `failed:tangent_no_resolved_drop`；该簇 rank=1 的 `k5 / distance:6-9` 在第二轮 connected，才补齐第二类事件。若硬删第二轮，本次 sB 只剩一类已连通事件。此前旧批次“含事件簇首都连通”的观察不能推广到新运行。

两块的 property_status 都仍为 not_run。本记录只核对 E/F 搜索与连通、几何事件身份和费用，不构成 index-one、IRC 或 QM 验收。

## 新暴露的括区接口缺陷

HS 出现 `error:ValueError` 7/4 项，全部属于 branch_bracketed，全部 0 forward。原始异常是 `bracket must have positive width`，从 `_hs_task` 调用 `continue_branches_ef` 时触发。

对已保存的 L1/L2 重建输入括区：10 项两端 q 完全相等，1 项只差 2.22e-16 Å（重算浮点差异）。`bracket_jump` 的合同是在同一个 q 上确认两分支共存；`_hf_bracket` 却把反应物点和回扫产物点作为 q 区间交给要求正宽度的延续器。问题在分支证据与搜索括区的交接，不能靠缩放聚类半径或任意加一个 q 偏移解决。

下一步先明确并修复这道交接：保留同 q 的两支证据，区分无有效区间与合法搜索输入；添加正向/反向扫描、完全同 q 和机器精度宽度的回归，再重放现有 branch_bracketed 记录。修复前不据本次结果宣布完整验收，也不启动下一轮昂贵对比。本轮没有修改搜索/调度代码或历史运行结果。

## 补充（Claude，2026-10-28；只读 report / journal）

- 上表关键数字已用四份 report 独立复核：refine 花费、HS 结算 290/275、`HS_deferred:same_start` 42/55、`error:ValueError` 7/4、第二轮实际派发 129−42 = 87、125−55 = 70，全部一致。
- **括区缺陷自 HF-4 起就存在，只是被预算挡住了。** 10-26 两块 v1 中没有任何 `branch_bracketed` 项结算过 HS：升级排序里脊在前、跳支在后，refine 额度在轮到它们之前就用完了（`HS_skipped:refine_budget` 89/96）。v1.1 第一轮覆盖每个簇，才第一次把这类项送进 HS（sA 11 项：7 error + 4 unresolved；sB 7 项：4 + 3）。所以这不是两轮派发引入的回归。HF-4 的解析台测试也没有覆盖跳支括区（【未核实】TwoSaddle 台是否会产生 `branch_bracketed`）。
- 新两块各有 7 / 3 个 `hf_duplicate`，`connected` 各 2 个：10-28 修复的去重规则（`pair_symmetries`）第一次在真实运行里生效，对称等价的 E2 被合并，没有再报成 4–6 个“不同事件”。
