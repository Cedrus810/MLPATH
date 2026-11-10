# 施工计划：HF-5 跨体系检验（b04 有界端到端 → N1 协同体系；D2-HF-v1.4 规则不变）

日期：2026-11-08。状态：b04（a0–a5）完成；N1 n0、n2 完成，n3（v1.4）不做；小结见 [HF5_SUMMARY_2026-11-10.md](HF5_SUMMARY_2026-11-10.md)（用户 11-10：不再逐体系调参数）。

## 必读（开工前全部读完）

- `docs/TS_HESSIAN_FREE_DESIGN_RECORD_2026-10-24.md` §7 HF-5（本计划的出处与完成条件）、§5.3（阈值按实测分辨率）、§6（搜索 / 性质 / 连接 / QM 分开报告）。
- `docs/PLAN_HF_D2_V14_2026-11-07.md`：要检验的规则 D2-HF-v1.4 及 b05 结果。
- `docs/PLAN_HF_MIDPOINT_COARSE_2026-11-02.md` §2：符号界线 B 在 b05 上怎么测、公式 B = 0.003 + 3σ_c 的来历（本计划在新体系上照做）。
- `docs/EVALUATION_PROTOCOL_DESIGN_2026-10-13.md` §3.1–3.2：b04 / N1 / b09 各自的角色、参照认证流程；为什么 b09 不进开发。
- `docs/HANDOFF_B09_GAP_2026-10-10.md` §1：b09 的定位。
- `benchmarks/explore_d2/run_b05_d2.py`、`run_b04_d2.py`、`l0.py`：现有 runner；本计划的新 runner 由它们合成。
- `benchmarks/explore_d2/midpoint_sign_resolution.py`、`hf_ef_rf_jitter.py`：B 与 σ / 几何噪声的测量脚本。
- `benchmarks/explore_counts/results/2026-10-22/COUNTS.md`：各体系意图数（定规模用）。

## 用户原话与依据

- 2026-11-08：问"要我写 HF-5（b04 → b09）的施工计划吗？"，回答「可以」。
- b09 的定位记在 `docs/HANDOFF_B09_GAP_2026-10-10.md` §1（标"用户 2026-10-10 确定"）：b09 是终极验收项目，当前只测差距、不在当前方法上调通。评估协议 §3.1（用户 10-13 审定）把 b09 列为"最终压力测试"，开发用解析台、b04、b05 / N1。
- 2026-10-22（DONE_LOG 记录的原话）："的确 09 太慢了……同意吧"——b09 全覆盖不可行，L1 必须抽样；U1 / U3 / U4 / U5 取值未给。

**与你说的不同：** 上一轮我说"b04 → b09"；本计划改为 **b04 → N1（环戊二烯 + 丙烯醛 DA）**，b09 不进 HF-5。原因：① 按评估协议和 b09 交接，b09 是最终验收 / 压力测试，HF-5 要按测量结果定每个体系的参数（σ、B）并看规则在新体系上是否成立——这是开发用途，在 b09 上做会污染最终验收；② N1 正是协议里的"新双分子参照"，两根键协同生成、有 endo / exo 两个鞍点，符合设计 §7"更软 / 协同体系"；③ b09 的 L0 取值（U1 / U3 / U4 / U5 与抽样方式）至今没定。b09 留到规则冻结之后作最终检验，另立计划。

## 1. 目标与完成条件（设计 §7 HF-5）

在 b05 以外的体系上用**不改动的** D2-HF-v1.4 规则端到端跑，回答：

1. v1.4 在新体系上能不能找到参照鞍点（召回），找到的候选经独立性质检验是什么结果；
2. "同键图的脊不出候选"（S3 的前提）在新体系上是否成立——看抽查样本里有没有候选；
3. 成本：各阶段 forward、墙钟；
4. 设计 §7 要求的敏感性：旋转、步长、初始化（种子）、重复运行（b04 上做全，N1 上只做种子与重复）。

搜索、性质检验、连接、参照匹配四项分开报告，检验未决不抹掉搜索几何（设计 §6）。

## 2. 不做

- 不改 `explore/src/`（`d2.py`、`ts_ef.py` 等）、不改 v1.4 规则与任何门。新体系上若出问题，只报告，改法另议。
- 不改 `run_b04_d2.py`、`run_b05_d2.py`（保留 b04 / b05 历史运行可复现）。
- 不做 b09、Dimer 基线（HF-6）、QM 核验（N1 参照的 QM 轴另立计划）。
- 不评定 b04 的冻结 `CRITERIA.md`（这里只是 HF 搜索的有界核对）。

## 3. 代码（全部在 `benchmarks/` 下）

| 文件 | 内容 |
| --- | --- |
| `benchmarks/explore_d2/systems.py`（新） | 每个体系一条：输入（路径或构建函数）、片段电子态、`HFSettings`（σ、几何噪声、中点 1e-3、hs_cap / hc_cap 3000）、`HSBracket`（粗阈值、B）、噪声。b04 的 σ 0.0042 / 几何 6.2e-4 Å 照抄 `run_b04_d2.py`（10-25 实测）；b05 一条照抄 `run_b05_d2.py`，只作一致性对照；N1 的数值在 §5 测出后填，测出前为空，用到时报错 |
| `benchmarks/explore_d2/run_d2.py`（新） | `run_b05_d2.py` 的 CLI 加 `--system`；`--measure-noise` 与 `--run` 两种模式；v1.1–v1.4 开关同 b05；新增 `--rotate SEED`（L0 之前对输入做一次随机刚体转动，转动矩阵写进 manifest 的 identity）与 `--magnitude-step`（扫描幅度网格步长，默认 0.1 Å，即现行 `MAGNITUDES_A`；随 `D2Settings.magnitudes_A` 冻结进 manifest） |
| `benchmarks/n1_diels_alder/`（新） | `build.py`：RDKit MMFF 构建环戊二烯、丙烯醛（s-cis / s-trans），POLAR-1 优化成片段极小与反应前复合物，写 `input.extxyz`；`reference.py`：endo / exo × s-cis / s-trans 四个 TS 初猜 → POLAR-1 上 geomeTRIC 鞍点优化（照 `carbene_ambimodal/optimize_b08_b09.sh` 的做法）→ Hessian 阶数 → 两侧下降与端点身份 → `reference.json`。这是参照一侧，可以用已知结构和 Hessian；搜索一侧看不到它 |
| `benchmarks/explore_d2/hf_event_validate.py`（新） | 读一次 D2 运行，对每个 HC 连通的不重复事件：① 鞍点的 POLAR-1 Hessian 阶数（性质检验，只在搜索之后做）；② 与 `reference.json` 匹配（`candidate_identity.py` 的对称 RMSD ≤ 0.15 Å 且 ΔE ≤ 0.05 eV，与运行判重同一规则）；③ 输出召回表：每个参照鞍点被几个意图 / 事件找到、首次命中时的累计 forward |

一致性检查（冒烟）：`run_d2.py --system b04 --hf --start-cluster` 用 10-28 `b04_hf_v11_smoke` 的参数，结局类别应与那次相同（HS converged 2 / hysteresis 2 / above_ceiling 2，HC 连通 SN2 1 + 判重 1）；forward 有 fp32 小差异。

## 4. b04（HF-5a）

| 步 | 内容 | 产出 |
| --- | --- | --- |
| a0 | 上面的一致性冒烟（183） | runner 可用；测单步耗时 |
| a1 | **v1.2 基线整跑**（`--hf --start-cluster --same-graph-ridge`，不加括区规则），种子 sA | 对照臂；同时提供 a2 测 B 所需的中点 |
| a2 | `midpoint_sign_resolution.py` 读 a1 的中点，测 σ_c 与粗精偏差 → B_b04 = 0.003 + 3 × p90 σ_c | **B_b04 交你定**（b05 的 B = 0.025 由你 11-02 定，公式本身标为"待验证"）；粗阈值沿用 2e-3，除非实测 σ_c 显示 2e-3 时符号分不清 |
| a3 | **v1.4 整跑** 种子 sA、sB | 主结果 |
| a4 | v1.4 敏感性：sA 同命令重复 1 次（fp32 重复噪声）、sA 加 `--rotate`、sA 加 `--magnitude-step 0.05` | 设计 §7 的旋转 / 步长 / 重复；种子差异由 a3 的 sA / sB 给出 |
| a5 | 每个运行跑 `hf_event_validate.py`（b04 参照：SN2 鞍点 −26132.2012 eV，QM 已核验，`benchmarks/b04_sn2_chloride/reference.json`、`docs/experiments/b04_qm/`） | 召回、性质、成本表 |

L0 与预算：
- 壳层 2.0（b04 实测接触比 f = 2.02，与 b04 冒烟和 count 表同一规则）；`n_samples`、`n_far`、`per_pose_quota`、各 cap、`audit_fraction` 0.1、`max_active` 8、`wave_size` 16 沿用 b05 v1.4 整跑的值，使两体系可比。
- search 额度 = 实际意图数 × a0 实测的每意图平均 forward × 2；refine 额度 = 升级项数 × 6000（每项 HS + HC 上限），即 b04 上不让额度成为瓶颈。实际花费如实报告。

## 5. N1（HF-5b）

| 步 | 内容 | 产出 |
| --- | --- | --- |
| n0 | `build.py` + `reference.py`（183，CPU / GPU 都可） | `input.extxyz`、`reference.json`：四个 TS 初猜各自是否收敛成一阶鞍点、是否互不相同、两端是否是反应物 / 降冰片烯类产物；不合格的项照评估协议 §3.1 记录、不进分母 |
| n1 | `run_d2.py --system n1 --measure-noise`（定 `noise_eV`，定义同 b05）；`hf_ef_rf_jitter.py` 在噪声测量扫描里的近峰点上测反力 σ 与几何噪声（同 b05 J0 的定义；点取自 L1 扫描，不用参照结构） | `systems.py` 的 N1 一条：noise、σ、几何噪声 |
| n2 | v1.2 基线整跑 sA → `midpoint_sign_resolution.py` → B_N1 | B_N1 交你定 |
| n3 | v1.4 整跑 sA、sB；sA 同命令重复 1 次 | 主结果 |
| n4 | `hf_event_validate.py` 对所有运行 | 召回（endo / exo 各被找到几次）、性质、S3 样本、成本 |

L0：壳层取 n0 反应前复合物的实测接触比 f（同 count 表规则）；其余同 b04。N1 每姿态 43 个单原语 + 二元组合，6 个构象配对（count 表），规模可全覆盖。

## 6. 计划表（每次正式运行前再按当时负载报一次，等你说"跑"）

| 运行 | 机器 | 资源 | 预计耗时与依据 | 续跑 | 输出 |
| --- | --- | --- | --- | --- | --- |
| a0 冒烟 | 183 | 1 GPU，8 线程 | 几分钟；10-28 同参数冒烟 refine 740 forward，HF 段 0.089 s / forward（2080 Ti） | 能 | `benchmarks/explore_d2/results/<日期>/b04_hf5_smoke/` |
| a1、a3、a4（6 个运行） | PBS 双卡节点或 183 | 每运行 1 GPU、16 线程 | 不知道，a0 后按实测每意图 forward × 意图数 × 每 forward 墙钟外推 | 能（同主机重跑同命令） | `results/<日期>/b04_hf5_*/` |
| a2 | 183 | 1 GPU，8 线程 | b05 那次（11-02）的规模与耗时作上限参考；b04 每 forward 更快 | 能（逐文件跳过） | `results/<日期>/b04_sign_resolution/` |
| n0 | 183 | 1 GPU 或 CPU 16 核 | 不知道；b09（29 原子）6 个鞍点重优化有先例，N1 19 原子、4 个初猜 | 能（逐项跳过） | `benchmarks/n1_diels_alder/results/<日期>/` |
| n1–n3（5 个运行） | PBS 双卡节点 | 每运行 1 GPU、16 线程 | 不知道，n1 后外推；count 表的旧外推（逐点串行 fp64）为 44 GPU h（单原语），v1.4 批量 fp32 下会少得多，只作上限量级 | 能 | `results/<日期>/n1_hf5_*/` |

两块各占一张卡的布局照 `b05_hf_v14_pbs.zsh`。只交冒烟跑通过的代码（a0 覆盖 b04 路径；N1 路径由 n1 的噪声测量先跑通）。

## 7. 报告（每个体系一张表，v1.2 基线与 v1.4 并列）

- 参照召回：每个参照鞍点是否找到、被几个意图找到、首次命中的累计 forward；未找到的留在分母。
- 不在参照里的连通事件：列出，做性质检验，不当成错误、也不当成新发现，交你看。
- 性质检验：阶数通过 / 不通过 / 未决（TODO 第 12 条的交叉核验门限问题可能导致未决，照实记）。
- S3：同键图的脊抽中 / 跳过数，样本 HS 结局，**样本出候选就醒目列出**。
- 成本：search / refine / HS / HC forward，墙钟。
- 敏感性（b04）：重复、旋转、步长、种子各一行，比召回与事件身份。

## 8. 待你定的事

1. 本计划（含"b09 换成 N1"）；
2. B_b04、B_N1（a2 / n2 测出后）；
3. 各正式运行的节点与开跑时间。

## 未确认（Claude 推断）

- search 额度"× 2"的余量是我定的，没有推导；目的只是不让 L1 被额度截断，实际花费会如实报告。
- N1 的 endo / exo 鞍点在 POLAR-1 上是否都存在、是否是协同一步（而不是分步）不知道，要 n0 才知道；如果只剩一个，N1 仍可检验 S3 与成本，但召回只有一项。
- "b04 上同键图的脊很少"是从它只有 6 原子推的，没有数过；S3 在 b04 上可能几乎不被触发，所以 b04 主要检验召回、成本与敏感性，S3 主要看 N1。
- 粗阈值 2e-3 沿用到新体系，依据只是 b05；a2 / n2 的 σ_c 若显著大于 b05（中位 0.0040），要回来报告。
- 运行目录里的 `journal.jsonl` 可能很大；提交草稿只列脚本、文档和汇总文件（依据是你 11-08 说的"块太大了 我git会爆炸"，具体改法你未定）。
