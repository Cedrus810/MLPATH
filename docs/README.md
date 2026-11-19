# MLPATH 文档入口

**先读 [项目总览](PROJECT_OVERVIEW.md)**：目标、数据流、模块边界与后续顺序。2026-10-24 已按用户明确的“面向 MLP 的 Hessian-free TS 搜索”整理主线；细节见 [设计记录与接手路线](TS_HESSIAN_FREE_DESIGN_RECORD_2026-10-24.md)。各工作线的科学结论仍以其协议和工件为准。

本索引负责告诉你去哪找资料。日期计划和交接记录保留原路径，避免打断脚本与记录的引用；其中“当前”“下一步”“未实施”描述的是相应时期。最新实现状态见总览及各工作线后续记录，冻结判据与用户裁决仍以原件为依据。

## 按问题进入

| 你要了解什么 | 入口 |
| --- | --- |
| 粗路径到 TS 的新颖性研究（2026-11-19） | [方向充分性与 E/F 定位](TS_DIRECTION_NOVELTY_2026-11-19.md)：String / SPM / BITSS / OCI-NEB 近邻、动作分歧指标、理论边界及最小消融；研究建议，未施工 |
| HF TS 论文初稿（2026-11-13，开发批次 v0.1） | [写作入口](../paper/README.md)、[英文正文](../paper/MLPATH_MainText.md)、[SI](../paper/MLPATH_SI.md)；含从现有汇总生成的两张图与证据表，采用一阶且匹配参照口径并记录镜像修正 |
| 首篇横向验证设计（2026-11-13） | [HORIZONTAL_VALIDATION_DESIGN](../paper/HORIZONTAL_VALIDATION_DESIGN_2026-11-13.md)：三个开发体系＋三个公开数据留出槽位，六体系局部、四体系端到端；共同 C2 与预算内完成计分，未实施 |
| 整个项目怎样连接 | [PROJECT_OVERVIEW.md](PROJECT_OVERVIEW.md) |
| 核心问题、当前缺陷、怎么修与未来方向 | [HF TS 设计记录](TS_HESSIAN_FREE_DESIGN_RECORD_2026-10-24.md)：八项缺口与 HF-1…HF-6；[HF-2/3 实现记录](HF_EF_LOCATOR_2026-10-25.md)：独立 E/F 原语、b05 丢模审计、解析测试及下一关；[跳支定位施工计划](PLAN_HF_JUMP_2026-10-25.md)：b05 E2 的 §5.2 双侧延续（J0/J1 已施工）；[E/F 连接追踪施工计划](PLAN_HF_CONNECT_2026-10-26.md)：F7 `trace_sides_ef`（K1–K3 已施工）；[HF-4 施工计划](PLAN_HF4_D2_2026-10-26.md)：E/F 链接入 D2 调度、账本与续跑（H1–H4 已施工，H5 待节点） |
| 10-13 评估协议与发现层重设计 | [评估协议设计](EVALUATION_PROTOCOL_DESIGN_2026-10-13.md)、[发现层设计](DISCOVERY_LAYER_DESIGN_2026-10-13.md)；用户 10-13 审定（按各自 v1 范围）。[施工计划](PLAN_EXPLORE_V1_2026-10-13.md)：M1–M4 已完成（账本、约束松弛探针、姿态与意图、C2 认证；b04 冒烟跑通、未验收）；[M5 D2 编排](PLAN_EXPLORE_M5_2026-10-23.md)已完成（两账户、分级、抽查、kill→resume；b04 冒烟得 c2，未验收） |
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
| 算力节约待办（重复响应、构象入口、暂缓状态） | [TODO.md](TODO.md)（2026-10-26 起）；[同一起点降优先级计划](PLAN_START_CLUSTER_2026-10-28.md)（已施工为 D2-HF-v1.1）；[延续试步失败处理计划](PLAN_HF_WALK_RETRY_2026-10-28.md)（现行为 v4：首次失败即终止该侧）；[结果记录](HF_WALK_POLICY_RESULTS_2026-10-29.md) |
| 二分中点卡住（TODO 第 9 条，11-01 至 11-02）：诊断、冷启动对照、横向曲率、阈值选项、分级终止（v4）与验收 | [冷启动对照计划](PLAN_HF_MIDPOINT_COLD_2026-11-01.md)；[横向曲率计划](PLAN_HF_MIDPOINT_TRANSVERSE_2026-11-02.md)；[改法选项](HF_MIDPOINT_THRESHOLD_OPTIONS_2026-11-02.md)；[分级终止施工计划（含用户决定与验收结果）](PLAN_HF_MIDPOINT_COARSE_2026-11-02.md)；[O5 优化器对照](PLAN_HF_MIDPOINT_O5_2026-11-02.md)；[精阶段 FIRE（v5）施工与验收](PLAN_HF_MIDPOINT_FINE_FIRE_2026-11-04.md)；[剩余精阶段卡点的改法选项（11-06）](HF_MIDPOINT_FINE_STAGE_OPTIONS_2026-11-06.md)；[粗阶段出候选（v6）施工计划](PLAN_HF_MIDPOINT_COARSE_CANDIDATE_2026-11-06.md)；[refine 端点出候选（v7）施工计划](PLAN_HF_ENDPOINT_CANDIDATE_2026-11-07.md)；[HS 两侧无最大值（R3）的改法选项（11-07）](HF_HS_NO_MAXIMUM_OPTIONS_2026-11-07.md)；[同键图的脊最后派发（D2-HF-v1.2）施工计划](PLAN_HF_SAME_GRAPH_RIDGE_2026-11-07.md)；[v7 接入 D2（D2-HF-v1.3）与 b05 整跑](PLAN_HF_D2_V13_2026-11-07.md)；[同键图的脊剩余花费的改法选项（11-07）](HF_SAME_GRAPH_RIDGE_SPEND_OPTIONS_2026-11-07.md)；[S3 + HS / HC 并行（D2-HF-v1.4）施工计划](PLAN_HF_D2_V14_2026-11-07.md) |
| 二分中点工作线交接（11-01 至 11-06：诊断、v4 分级终止、O5、v5 精阶段 FIRE、候选身份与 HC、待办） | [HANDOFF_HF_MIDPOINT_2026-11-06.md](HANDOFF_HF_MIDPOINT_2026-11-06.md)（接替 [11-04 版](HANDOFF_HF_MIDPOINT_2026-11-04.md)） |
| HF-5 / HF-6 工作线交接（11-08 至 11-12：b04 / N1 跨体系检验、HF-6 局部与端到端、检验脚本修正、待定事项；接替 11-08 版） | [HANDOFF_HF5_HF6_2026-11-12.md](HANDOFF_HF5_HF6_2026-11-12.md) |
| HF-6 小结（11-12：局部 HF 命中多、成本低约一个数量级；端到端 DB 0 命中；N1 两臂都漏协同 DA） | [HF6_SUMMARY_2026-11-12.md](HF6_SUMMARY_2026-11-12.md) |
| 三条通用规则的选项（11-12：同键图脊判定 G0–G3、σ / 几何噪声 / B 的测量程序 M1–M3、noise_eV 定义 N1–N3；待用户定） | [HF_GENERAL_RULES_OPTIONS_2026-11-12.md](HF_GENERAL_RULES_OPTIONS_2026-11-12.md) |
| D2-HF-v1.5 施工计划（11-13：G3 / M2 / N1 施工、冒烟、三体系 × 两种子整跑，参照覆盖与上一轮相同） | [PLAN_HF_D2_V15_2026-11-13.md](PLAN_HF_D2_V15_2026-11-13.md) |
| F8 探索 / 精修分波次交替施工计划（11-14：离线重排估算 W = 8 时 b04 / b05 前两类首次命中由约 0.9 B* 提前到 0.2–0.35 B*；搁置，未施工） | [PLAN_HF_D2_F8_WAVES_2026-11-14.md](PLAN_HF_D2_F8_WAVES_2026-11-14.md) |
| 算法结构缺陷清单（11-14 起，挖掘中：A 发射点按几何枚举、不看反应附近——三体系 84% / 65% / 84% 的发射无升级，far 层 380 次发射 0 事件） | [ALGORITHM_DEFECTS_2026-11-14.md](ALGORITHM_DEFECTS_2026-11-14.md) |
| Non-Hessian AFIR 设计稿（11-14 第 3 版，待用户审：反应逻辑照搬 AFIR——函数、γ、MC / SC-AFIR、导航与停止；TS 用 HF 反力零点内核在路径上重建括区，IRC 用 E/F 下降；不用 Hessian） | [NONHESSIAN_AFIR_DESIGN_2026-11-14.md](NONHESSIAN_AFIR_DESIGN_2026-11-14.md) |
| Non-Hessian AFIR 线交接（11-15：设计、代码、b04 结果、b05 / N1 PBS 运行与待办；接替 HF5 / HF6 交接） | [HANDOFF_AFIR_2026-11-15.md](HANDOFF_AFIR_2026-11-15.md) |
| HF-6 第二阶段：端到端 D2 对 DB（11-11 计划：同一 L0、同一总预算，R(b) 与首次命中成本） | [PLAN_HF6_DIMER_E2E_2026-11-11.md](PLAN_HF6_DIMER_E2E_2026-11-11.md) |
| HF-6 与 Dimer 比较（11-10 计划：局部定位先行，端到端另立） | [PLAN_HF6_DIMER_2026-11-10.md](PLAN_HF6_DIMER_2026-11-10.md) |
| HF-5 小结（11-10：b04 召回成立、S3 前提不成立；N1 单原语 D2 漏掉协同 DA；需统一重定的规则） | [HF5_SUMMARY_2026-11-10.md](HF5_SUMMARY_2026-11-10.md) |
| HF-5 跨体系检验（11-08 计划：b04 有界端到端 → N1 协同体系，D2-HF-v1.4 不改；b09 留作最终验收） | [PLAN_HF5_CROSS_SYSTEM_2026-11-08.md](PLAN_HF5_CROSS_SYSTEM_2026-11-08.md) |
| HS 精修与 D2 工作线交接（11-06 至 11-08：9c 诊断、v6 / v7、R3 普查、D2-HF-v1.2–v1.4、并行测量、v1.4 整跑与同硬件对照、数据被删、待办；接替 11-06 版） | [HANDOFF_HF_D2_2026-11-08.md](HANDOFF_HF_D2_2026-11-08.md) |
| b05 HF D2 工作线交接（10-26 至 11-01：去重、两轮派发、跳支括区、延续策略 v4、新 E2 鞍点与 QM 核验） | [HANDOFF_HF_D2_2026-11-01.md](HANDOFF_HF_D2_2026-11-01.md) |
| 如何复跑与验证改动 | 两套 Makefile；[env/BASELINE.md](../env/BASELINE.md) 是环境与历史验证记录 |

## 按工作线进入

| 工作线 | 主要资料 | 怎么读 |
| --- | --- | --- |
| 当前 HF TS 主线 | [设计记录与接手路线](TS_HESSIAN_FREE_DESIGN_RECORD_2026-10-24.md)、[HF-2/3 原语](HF_EF_LOCATOR_2026-10-25.md)、`explore/` | M1–M5 已有；D2 仍接完整 Hessian 精修。独立 E/F 定位通过解析测试，待真实模型验证与 HF-4 集成；HF 编号不与既有 M 编号混用 |
| PRRS 核心 | [HANDOFF](HANDOFF.md)、[09-20 计划](PLAN_2026-09-20.md)、[09-24 审查](REVIEW_2026-09-24.md) | 核心已有多轮实现；日期计划不是当前未完成清单 |
| 首次进入 committor | [首次进入协议](P1_COMMITTOR_FIRST_ENTRY_PROTOCOL.md)、[射击协议与执行记录](P1_COMMITTOR_PROTOCOL.md)、[09-25 计划](PLAN_2026-09-25.md) | R6/A5 已结清；k=160、cap=1280 等属于相应 P1 协议，不是普适默认值 |
| connect | [09-28 计划及进度](PLAN_2026-09-28.md)、`tests/test_connect_*`、`tests/test_relay_*` | G1–G3 已落地，G4 部分落地；G5 是探索，L/E 尚未实施 |
| PARS | [独立包说明](../pars/README.md)、[10-02 审查与 P1 实测](PARS_REVIEW_2026-10-02.md)、[原始设计](../Pareto-Adaptive%20Reaction%20Path%20Search%20%28PARS%29%20%E2%80%94%20MLPATH%20Integration%20Design.md) | select/refine 可用；10-05 电子态交接修复后真实 b04（POLAR-1）精修通过；P1、b04 都是单鞍点，不能证明多样性选择有效 |
| geomeTRIC 局部精修 | [GEOMETRIC_REFINEMENT.md](GEOMETRIC_REFINEMENT.md) | `prrs optimize` 支持极小／TS 候选优化，ASE 与 OpenMM；驻点与连接认证仍需另做 |
| 路径成熟层 | [独立模块 `mature/`](../mature/README.md)、[施工计划](PLAN_GEN_COMPAS_MATURATION_2026-10-03.md)、[Gen-COMPAS 可行性](GEN_COMPAS_FEASIBILITY_2026-10-03.md) | M0–M4 已实现，M2–M4 在解析势上验证；M5 真实 MACE 工作流实跑与 M6 比较待做，Q/B 为后续研究；不作为 HF TS 首轮定位的前置任务 |
| 性能 | [曲率复用与批量控制器（10-23）](SCAN_THROUGHPUT_2026-10-23.md)、[约束扫描性能挖掘（10-23）](SCAN_PERFORMANCE_2026-10-23.md)、[09-28 审查](REVIEW_2026-09-28.md)、[HANDOFF 实现沿革](HANDOFF.md) | 扫描的每次求值成本、每点优化步数与组合规模分开；先测量，再决定改动 |
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
