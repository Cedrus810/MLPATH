# MLPATH 文档入口

**先读 [项目总览](PROJECT_OVERVIEW.md)**：目标、数据流、模块边界、当前进度与后续工作，按逻辑顺序组织。状态核对日期：2026-10-08。

本索引负责告诉你去哪找资料。日期计划和交接记录保留原路径，避免打断脚本与记录的引用；其中“当前”“下一步”“未实施”描述的是相应时期。最新实现状态见总览及各工作线后续记录，冻结判据与用户裁决仍以原件为依据。

## 按问题进入

| 你要了解什么 | 入口 |
| --- | --- |
| 整个项目怎样连接 | [PROJECT_OVERVIEW.md](PROJECT_OVERVIEW.md) |
| 10-13 评估协议与发现层重设计 | [评估协议设计](EVALUATION_PROTOCOL_DESIGN_2026-10-13.md)、[发现层设计](DISCOVERY_LAYER_DESIGN_2026-10-13.md)；用户 10-13 审定（按各自 v1 范围）。[施工计划](PLAN_EXPLORE_V1_2026-10-13.md)：M1–M2 已完成（`explore/` 账本 + 约束松弛探针，b04 冒烟跑通、未验收），M3 起开工前另写逐步计划 |
| PRRS 数学与使用方式 | [onboarding/README.md](onboarding/README.md)；[中文 README](../README.zh-CN.md) |
| 哪些命令现在能跑 | [PRRS README](../README.md)、[PARS README](../pars/README.md)；各自 `--help` |
| 历史结论与项目规矩 | [HANDOFF.md](HANDOFF.md)；规矩在其 §2 |
| 冻结判据与基准结果 | [benchmarks/README.md](../benchmarks/README.md)、各体系判据和运行记录 |
| 新旧结果为何不同 | 对应日期报告、run manifest、[DONE_LOG.md](DONE_LOG.md) |
| 10-03 审查问题与修复 | [大问题复核](REVIEW_2026-10-03.md)、[10-04 修复与验证](FIXES_2026-10-04.md)；10-05 收口见 [完整性复审](INTEGRITY_REVIEW_2026-10-05.md)、[INTEGRITY_FIX](INTEGRITY_FIX_2026-10-05.md)、[复审的复审](REREVIEW_2026-10-05.md)；复审三项已修复，真实 b04 核验见 [DONE_LOG 末节](DONE_LOG.md) |
| 10-09 b09 探针与搜索设计修复 | [搜索设计修复](SEARCH_DESIGN_REPAIR_2026-10-09.md)（谱排序按执行运动打分、受限软模、先粗扫后细化、环内保键交付）；b09 结果见 [probe_b09/README](../benchmarks/carbene_ambimodal/probe_b09/README.md) 与 [FOLLOWUP](../benchmarks/carbene_ambimodal/probe_b09/FOLLOWUP_2026-10-09.md) |
| b09 差距测量现状（10-10 交接） | [HANDOFF_B09_GAP_2026-10-10.md](HANDOFF_B09_GAP_2026-10-10.md)：b09 定位为终极验收、当前只测差距；全部实验、已识别差距与待决定项 |
| PRRS 离 b09 验收差多远 | [B09_GAP_ANALYSIS_2026-10-10.md](B09_GAP_ANALYSIS_2026-10-10.md)：按目标逐级对照、差距分参数级 / 设计稿预见未实现 / 超出发现层三类；用户 10-10：不着急做 |
| 哪些决策尚未收口 | 项目总览 §9 → [DECISIONS_PENDING.md](DECISIONS_PENDING.md) 及对应工作线计划 |
| 如何复跑与验证改动 | 两套 Makefile；[env/BASELINE.md](../env/BASELINE.md) 是环境与历史验证记录 |

## 按工作线进入

| 工作线 | 主要资料 | 怎么读 |
| --- | --- | --- |
| PRRS 核心 | [HANDOFF](HANDOFF.md)、[09-20 计划](PLAN_2026-09-20.md)、[09-24 审查](REVIEW_2026-09-24.md) | 核心已有多轮实现；日期计划不是当前未完成清单 |
| 首次进入 committor | [首次进入协议](P1_COMMITTOR_FIRST_ENTRY_PROTOCOL.md)、[射击协议与执行记录](P1_COMMITTOR_PROTOCOL.md)、[09-25 计划](PLAN_2026-09-25.md) | R6/A5 已结清；k=160、cap=1280 等属于相应 P1 协议，不是普适默认值 |
| connect | [09-28 计划及进度](PLAN_2026-09-28.md)、`tests/test_connect_*`、`tests/test_relay_*` | G1–G3 已落地，G4 部分落地；G5 是探索，L/E 尚未实施 |
| PARS | [独立包说明](../pars/README.md)、[10-02 审查与 P1 实测](PARS_REVIEW_2026-10-02.md)、[原始设计](../Pareto-Adaptive%20Reaction%20Path%20Search%20%28PARS%29%20%E2%80%94%20MLPATH%20Integration%20Design.md) | select/refine 可用；10-05 电子态交接修复后真实 b04（POLAR-1）精修通过；P1、b04 都是单鞍点，不能证明多样性选择有效 |
| geomeTRIC 局部精修 | [GEOMETRIC_REFINEMENT.md](GEOMETRIC_REFINEMENT.md) | `prrs optimize` 支持极小／TS 候选优化，ASE 与 OpenMM；驻点与连接认证仍需另做 |
| 路径成熟层 | [独立模块 `mature/`](../mature/README.md)、[施工计划](PLAN_GEN_COMPAS_MATURATION_2026-10-03.md)、[Gen-COMPAS 可行性](GEN_COMPAS_FEASIBILITY_2026-10-03.md) | **M0+M1 已落地**（preflight 全 supported、schema/io + 53 测试，见 [preflight.json](experiments/maturation_preflight/preflight.json) 与 DONE_LOG 末节）；M2+M3 未实施，Q/B 为后续研究 |
| 性能 | [09-28 审查](REVIEW_2026-09-28.md)、[HANDOFF 实现沿革](HANDOFF.md) | PERF 候选与后续已落地项分开看；先测量，再决定改动 |
| 基准/QM | [基准队列](../benchmarks/MANIFEST.md)、[P1 QM](experiments/p1_qm/README.md)、[b04 QM](experiments/b04_qm/README.md) | 先读判据与来源；旧队列快照不代替后来协议 |

## 冻结判据与协议

以下原件保留；变更另写修订或新批次协议，不修改旧验收条件：

- [P1_MALONALDEHYDE_ACCEPTANCE.md](P1_MALONALDEHYDE_ACCEPTANCE.md)
- [P2_OXOBUTANAL_ACCEPTANCE.md](P2_OXOBUTANAL_ACCEPTANCE.md)
- [P2_MULTISEED_PROTOCOL.md](P2_MULTISEED_PROTOCOL.md)
- [P2_TEMPERATURE_SCAN_PROTOCOL.md](P2_TEMPERATURE_SCAN_PROTOCOL.md)
- [P3_BENZOYLACETONE_ACCEPTANCE.md](P3_BENZOYLACETONE_ACCEPTANCE.md)
- [P1_COMMITTOR_FIRST_ENTRY_PROTOCOL.md](P1_COMMITTOR_FIRST_ENTRY_PROTOCOL.md)
- [b04 CRITERIA.md](../benchmarks/b04_sn2_chloride/CRITERIA.md)

`P1_COMMITTOR_PROTOCOL.md` 含历史批次与后续确认结果，引用时注明批次和终止规则。09-22 二次扰动批次不可作为有效射击统计引用；后来修复/确认见后续章节。

## 历史、完成记录与理论资料

| 文件 | 用途与限制 |
| --- | --- |
| [DONE_LOG.md](DONE_LOG.md) | 完成事项、验证方式与实现沿革；按需追溯 |
| [ROUNDS.md](ROUNDS.md) | 09-04 重构与 QM 参照记录，承接 HANDOFF §10–§11 |
| [DECISIONS_PENDING.md](DECISIONS_PENDING.md) | 裁决、修订与追加的历史链；顶部旧队列不代替后来收口 |
| [PRRS_STATUS.md](PRRS_STATUS.md) | 08-31 实现快照；部分结论被后续测量推翻 |
| [PRRS_CURVATURE_ROUTES.md](PRRS_CURVATURE_ROUTES.md) | 曲率理论与部分测量；TS 判据勘误见 09-24 审查 |
| [PRRS_ENGINEERING_PLAN.md](PRRS_ENGINEERING_PLAN.md) | 初始工程设计，已有实现取代其任务状态 |
| [PRRS_CORRECTNESS_ROUND.md](PRRS_CORRECTNESS_ROUND.md) | 正确性审查与框架改动记录 |
| [P3_BENZOYLACETONE_PROBES.md](P3_BENZOYLACETONE_PROBES.md) | 探针记录；验收依据在 P3 acceptance |
| [09-05 审查](REVIEW_2026-09-05.md)、[09-04 快照](SESSION_2026-09-04.md)、[09-13 快照](SESSION_2026-09-13.md) | 对应日期工作记录 |
| [原始搜索方案](../mlp_native_perturbation_response_reaction_search.md) | 最初的方法意图 |

`design/` 是设计资料，`experiments/` 是测量/评定脚本及参照材料，`figures/` 是作图工具。通过引用它们的协议定位用途，不把每个脚本都当主流程入口。
