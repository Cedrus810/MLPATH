# PRRS — 交接文档

**阅读入口更新 · 2026-10-03**：先读 [MLPATH 项目总览](PROJECT_OVERVIEW.md) 和 [文档索引](README.md)。本文保留项目规矩（§2）、历史结论与交接沿革；各节“现在/下一步”对应写作时期。当前跨模块进度统一看总览，实验协议与裁决仍引用原件。

**最后更新** 2026-09-30（connect 测试台 T-绕路变体落地：`ProtonRelayGated`，指纹 → **`ebff4ec0`**；同日较早：connect 模式 G4 落地：R→P 路径分 observed / any_evidence 两视图记录，指纹 → `2cbcaec8`；此前 2026-09-28：续跑丢失方向停止原因的主线 bug 修复，指纹 → `bb525621`；同日较早：connect 模式 G2 = 干预 I-a 落地：队列按到产物的 d 出队，指纹 → `5487b5f1`；同日较早：connect 模式 G3 = 干预 I-b 落地：方向按到产物的键差提前，指纹 → `47803359`；同日较早：性能审查 PERF-22 落地：`state.encode` 向量化，指纹 → `22d90bd9`；同日较早：connect 测试台 T0 的 T-链变体落地：`calculators.ProtonRelay`，指纹 → `1098bcc8`；同日较早：性能审查 PERF-03a 落地：`symmetric_rmsd` 不再逐群元素拷贝 Atoms，b05 上 159 → 79 ms/次，指纹 → `f03de1bf`；同日较早：connect 模式 G1 落地：产物只记录不路由，指纹 → `06d6b071`，计划见 `PLAN_2026-09-28.md`；同日较早：性能审查 PERF-02 落地：`audit_basin` 只建一个 calculator，指纹 → `a13da932`；同日 PERF-11 `prrs_watch` 不进指纹；清单与协议见 `REVIEW_2026-09-28.md`；此前 2026-09-27：首次进入 committor 主线 A 完成：确认批四窗口全过 + 生产实现逐发等价，指纹 → `b601b124`；随后 C1 样本分类 + D1 下降几何，指纹 → `99edddde`；2026-09-27 电荷 / 多重度整数校验，指纹 → `b2c05652`；支线 B0 / B2 / C2 / C3 已跑，新决策 R8 待定；此前 2026-09-25：外部审查核对 + 续跑回滚重做 + 四处小修，指纹 → `ae47c4d9`；**正文见 `REVIEW_2026-09-24.md`，交接见 §15**；此前 2026-09-24：committor 射击二次扰动修复 + `runs/p1_committor2` 重跑：09-22 的「ridge 代理有偏」被推翻；**正文见 `DONE_LOG.md` 末节**；此前 2026-09-22：kasuga02 实跑队列七步全执行 + committor 判归修复；接手请先读 `DONE_LOG.md` 末两节与 `DECISIONS_PENDING.md` 的 2026-09-22 收口；文档索引见 `docs/README.md`）· 环境 `openmm_dev` · 工作目录 `/home/ruigengji/MLPATH`
（~~不是 git 仓库~~ 是 git 仓库，`main`；但 2026-09-25 之前工作区长期大量未提交，所以 `implementation_sha256` —— 对 `src/prrs/*.py` 的哈希 —— 仍是**判定用**的代码指纹）

**当前实现 `ebff4ec05fbc0c3d`**（2026-09-30）。别手抄这个数 —— 用 `make hash`
（`scripts/implementation_hash.py`）现算，它和 manifest 用的是同一个算法。
这一行以前是手抄的，抄完就过期了：它写着 `9c13fdb5` 时树已经是 `dba26be3`，中间隔了两轮修复。

沿革（每一步都是一次 `src/prrs/*.py` 的改动）：

| 哈希 | 日期 | 内容 |
| --- | --- | --- |
| `0058158d` | 2026-09-02 | **P2 冻结点**，`runs/p2_round7` 与 16 个 `p2_T300_seed*` 都在它上面 |
| `158f2f3b` | 2026-09-02 | P3 全部探针记录在它下面产生 |
| `9c13fdb5` | 2026-09-03 | 四项确定性缺陷修复 |
| `2f02fef9` | 2026-09-03 | `internal_mode_recheck` 的 fail-open 修复（`runs/p3_descent_recheck`） |
| `dba26be3` | 2026-09-03 | `quench` 两处收敛缺陷 + `follow_unstable_mode` 软模接纳（`runs/p3_descent_recheck_fixed`） |
| `23179dea` | 2026-09-04 | **工程化重构，零行为改动**，见 §10 |
| `cf29f66b` | 2026-09-04 | `_stamped_with_charge`：声明的总电荷必须到达势能，冲突则拒绝开跑 |
| `686a5c4f` | 2026-09-04 | `side_seed`：修负 seed（只在 `minimum_check="probe"` 下可达，未打到任何冻结结果） |
| `5d905aaa` | 2026-09-04 | `IdentityUnavailable`：枚举预算溢出不再成为身份；`admit` 的三层分离 |
| `6b372359` | 2026-09-04 | 未知记为未知（`locked_bond_parity` 为 `None` 而非 `[]`）—— **本轮唯一的身份迁移** |
| `69d53501` | 2026-09-04 | 失败原因移出哈希；`withheld_candidates` 检索路径 |
| `5b69dfed` | 2026-09-05 | **三项正确性修复**（见 §12）：验证不完整不再当通过 · 候选出口统一 · 失败落盘独立于网络快照 |
| `d7d53219` | 2026-09-09 | **扰动交付**：bond 位移改刚体片段平移（`internal.translate_fragments`），torsion 的 `rotate_fragment` 终于有了 bond 的对应物；环上键与同碎片 1–3 对仍走 continuation |
| `8418b82b` | 2026-09-09 | **bracket 高端派生字段**（REVIEW §7.3）：细化换掉高端后 `target_microstate` 不再留着 coarse 高端的值；`_bracket_high` 一处定义，不改搜索决策 |
| `093bf311` | 2026-09-09 | **`free_bonds` 并集 fail-closed**（REVIEW §7.4）：另一侧测过同一根扭转而判的不是 free 时不再商掉它（`measured_bonds_of`）。该分支在语料里从未执行过，是预防不是修复 |
| `5dc83ed4` | 2026-09-09 | **方向合同第 2/3 条**（用户裁决）：`candidate_id` 稳定 ID、`coverage_censored` / `censored_by` / `uncovered_ids` 进提案账本。只加记录字段，行为零改动 |
| `fa15a829` | 2026-09-09 | **方向合同第 4/5 条**：`DirectionSchedule` 从 `run()` 抽出（只搬调度这一片，策略未改）· `direction_ledger` 两个单位各自自洽 · `direction_mode = sampled/exhaustive`，exhaustive 预算不足则**开跑前**拒绝 |
| `7636cd82` | 2026-09-09 | 内层预算名归位（`inner_blocked`）：批中途被 `conformer_trials_per_microstate` 截断不再错记成泛化 `node_budget` · b04 判据 §5 口径收窄为「该微态在枚举提案基内的可采纳负记录」 |
| `1b96d400` | 2026-09-09 | **bend 的刚体交付**（`rotate_about_vertex`）：两侧片段绕顶点各转 δ/2，非目标键变化 2e-16 Å。三种坐标至此都有精确刚体实现；同碎片非键对与 quench 预算两项**未做**，各自另议 |
| `6b293d16` | 2026-09-09 | exhaustive 开跑前检查覆盖全部三道预算（b05 probe2 发现 `conformer_trials_per_microstate` 会绕过它） |
| `460dd519` | 2026-09-09 | **自同构枚举改逐点回溯**（用户裁决：走剪枝不提限）：b05 的 1296 阶群从 2e6 候选 / 19.8 s 拿不到，变成 3889 nodes / 28 ms；permutation 集合与旧实现逐位相同、digest 零差异；预算单位改为 `search_nodes` |
| `89f41c80` | 2026-09-09 | 运行级账本改为按 (微态, 方向) 对**求和**而非按 ID 取并（b04 probe4 抓到 28 ≠ 17+7）· 新增 `termination_detail.directions_dropped_by`：`termination` 只回答 trial 预算，别的预算按名另记 |
| `c393da78` | 2026-09-13 | **三路审计落地**（P0 生成群闭包、P1-6 (0,π) 守卫、check.py 对齐、规矩 5 改写）· 商层置换不变性 fuzz · 谱提案 · 碎片刚体转动 · b04 冻结 + `runs/b04_frozen1` 首评 `observed_direct` |
| `b4e92642` | 2026-09-20 | **PLAN 七决策落笔**（`DECISIONS_PENDING.md` 2026-09-20 收口）· **第 4 件 schema 显式化**（`schema.py` + `publish()` 记录 `schema_warnings` + CLI `--strict-schema` + runs/ 全量锚测试与漂移清单）· **P1-4 溢出拒建**（决策 7：`admit` 拒 `automorphism_budget`，不再拿 `[identity]` 当群） |
| `1a6a2971` | 2026-09-20 | **第 6 件收尾**：resume 全件——三条硬前提、盘上重建（queue 落盘 / 计数器数 paths 目录 / schedules 逐向恢复）、scan 边界回滚记 resume.json、**npz sidecar 保位精确**；验收 = kill(os._exit)→续跑→与不中断孪生 network.json/attempts.jsonl 逐字节相同 |
| `b1c4d440` | 2026-09-20 | **第 1 件正则标号全件**（`canonical.py` IR 引擎、`structural_key`、轨道取代色类：P1-3 关、P1-7 关、xfail 同次摘；先量 2136 帧零差异 ⇒ 免费）· 第 2 件 committor 机器+接线（T=0 拒绝、Wilson、bracket 记录；实跑排队 MACE 机）· 第 3 件证书机器+解析势链条（命题待批故恒盖 DRAFT 章）· 第 5 件 5.1/5.4/5.2(惰性)+测量脚本+5.0 协议 · 第 6 件 6.1/6.2（from_network 逐字段等价）· logging/max_depth 随手件 · **digest 全程零差异** |
| `ec32bc3d` | 2026-09-21 | manifest 记 `host` · `make fmt` 全量（15 个文件自 09-13 漂出格式而每道门都绿）· `make check` 加 `fmt-check` |
| `63be7064` | 2026-09-22 | **committor 判归修复**：`Registry.identify`（只读的 `key_of()["key"] → find_node → _match_microstate`）取代读 `record["target_microstate"]` —— 那个键只有扫描路径在 `admit` 之后才写，所以射击的 A/B 两个分支**在代码里不可达**，首次实跑四个 bracket 全报 `A 0 / B 0 / other 16`。另加 `landings`（每发落在哪个微态）。digest 零差异（committor 默认关，语料无此路径）· 正文见 `DONE_LOG.md` 末节 |
| `4014c4b2` | 2026-09-23 | **committor 射击不再二次扰动**：`shoot` 把 ridge 帧交给 `run_trial`，后者又施加了一遍扫描的 `probes[0]`（0.15 Å 拉伸 +0.405/+0.491 eV 势能，或 0.3 eV kick），对 419 meV 的势垒。`run_trial(deliver=False)` 只留热动量，射击记录带 `probe_delivery`。digest 零差异 · 重跑 `runs/p1_committor2` 推翻 09-22 的「ridge 代理有偏」· 正文见 `DONE_LOG.md` 末节 |
| `ae47c4d9` | 2026-09-24 | **续跑回滚重做**：旧回滚以 network.json 为边界、保留它引用的 attempt，而 `publish()` 每个 trial 后都跑 —— scan 中途 kill 后记录有洞、trial id 撞目录（kill@4 → `FileExistsError`），且在飞微态已出队、被静默丢掉。改为微态出队前写 `resume_point.json` 快照、按计数回滚（含 climb / continuation / descent / ts / withheld）· resume 核对 `input_structure_sha256` · `completeness.cover` 要求正交归一基底 · committor 记录加 `p_B_conditioned_on` / `other_fraction` / `max_nve_drift_eV_atom`。476 passed / 1 skipped，digest 零差异，变异检验五次全红（**本机 kasuga01**，见 §15）· 正文 `REVIEW_2026-09-24.md` |
| `b601b124` | 2026-09-26 | **committor 首次进入（A5）**：`run_trial(first_entry=…)` 逐步计数、连续 k 步在 A / B 图上即停，之后原有 quench 就是入口核对；`shoot` 透传；配置 `committor_entry_steps`（k，0 = 关，默认）/ `committor_max_steps`（1280）；search 的 committor 分支按进入侧判，核对不一致只计数不改判，没有 B 或两端同图时退回固定规则并在 `termination` 写明；R1：`p_B` → `p_B_given_AB`。 kasuga02：`make check` 绿、digest 零差异 332 / 2136 / 4；变异检验两条都红（`>=`→`>` 进入步差一；删 `break` 跑满上限）；**真实数据核对** `runs/p1_first_entry_production_check`：生产路径重打 A3 的 200 发（同帧同 seed、提前停），**进入侧与进入步 200/200 相同、入口 quench 落点 200/200 相同**，约 23 min。 正文 `DONE_LOG.md` 末节 |
| `99edddde` | 2026-09-26 | D1：鞍点双侧下降改用**一个** `PathRecorder` 管两支，帧计数连续，`descent.extxyz` 不再被第二支覆盖（此前历来只剩 minus 支的几何）。 C1：scan 每个样本加 `outcome`（returned / departed / ridge / transient / unresolved）与 `outcome_reason`（未定时的失败码或 `admission:…`），只记录不路由。 测试 `tests/test_scan_records.py`（双阱真实搜索里两支帧都在且与日志一一对应；各分类单测；端到端字段一致）。变异检验两条都红（D1：记录器挪回循环 → 帧数 ≠ 行数；C1：returned / departed 对调 → 两条测试红），脚本还原并核指纹 `99edddde`；变异检验在 **yayoi**（PBS 作业）上跑。未变异全量：**用户报告通过**（未贴 digest 行）。 |
| `b2c05652` | 2026-09-27 | `SearchConfig`：`total_charge` 设了值必须是 `int`，`multiplicity` 设了值必须是 >= 1 的 `int`（`type(v) is int`，`bool` 也拒）。此前身份层（`network.Registry` → `chemical_key`）原样收下 `total_charge`，势能侧（`_stamped_with_charge`）取 `int()`：1.5 会以 1.5 记身份、以 1 算能量，`True` 静默成 +1。 测试 `test_charge_and_multiplicity_must_be_integers_when_set`，修复前先红（`DID NOT RAISE`）。本机（用户授权）：全量 484 passed / 1 skipped、fmt-check / lint 绿、digest 零差异 332 / 2136 / 4。仓内配置（b04、b05、preflight）都是整数，不受影响。 |
| `a13da932` | 2026-09-28 | **PERF-02**（`REVIEW_2026-09-28.md`）：`completeness.audit_basin` 整个审计只调一次 `factory()`，能量循环和 `membership` quench 共用同一个 calculator。此前每个覆盖点调两次，b04 审计（619 点）会重载模型约 1237 次。测试 `test_audit_builds_one_calculator`；变异检验：恢复旧实现后变红（`assert 14 == 1`）；新旧 audit 记录逐字节一致（解析三阱，d1 9 个样本、d2 37 个样本）。本机（用户授权，并已和同目录会话 mlpath-35 商量好先后）：`make check` 全绿，488 passed / 1 skipped，digest 零差异 332 / 2136 / 4。488 里另含 PERF-11 的 3 条新测试（`tests/test_prrs_watch.py`；`scripts/prrs_watch.py` 不进指纹；三个变异体全红）。 |
| `06d6b071` | 2026-09-28 | **connect 模式 G1**（`PLAN_2026-09-28.md` §3，只记录不路由）：新模块 `product.py`；`ReactionSearch.run(product=…)` / CLI `search --product`。产物原子数或元素序列与反应物不一致、或反应物自同构群溢出，**开跑前拒绝、什么都不写**（K4，不做 MCS）；给了产物才写：manifest `product_structure_sha256`、`product.extxyz`、network 顶层 `product` 块、每节点 `product_compatible`（同 `Registry.key_of`）、每微态 `d_to_product`；续跑时产物不一致即拒（回滚之前）。**d 在 Aut(G_R) 下取最小**（K5 修订：Aut(G_P) 下取最小是恒等）。schema：出现 `product` 块则所有节点 / 微态字段必须齐。不给产物时逐字节不变。测试 `tests/test_product.py` 13 条（甲醇 → 甲醛 + H₂：换一个等价甲基 H 做反应 d 仍为 0；双阱端到端；续跑三种拒绝）。变异检验五条全红（去掉群最小 / 元素不一致放行 / 续跑不查产物 / schema 不查 d / 群溢出放行）。本机（用户授权）：全量 501 passed / 1 skipped、fmt-check / lint 绿、digest 零差异 332 / 2136 / 4 |
| `f03de1bf` | 2026-09-28 | **PERF-03a**（`REVIEW_2026-09-28.md`）：`chemistry.symmetric_rmsd` 直接置换坐标数组，不再对每个群元素做 `Atoms.copy` + `set_positions`；Kabsch 从 `state.aligned_rmsd` 抽成共用的 `centred_kabsch_rmsd`，算术逐位相同。先测后改：b05（|Aut| = 1296）上 159 → 79 ms/次，`free_aligned_symmetric_rmsd` 两次调用 7.66 → 4.05 s；改前冻结的判据是「每次调用降一半且结果逐位相同」，基准的全部结果（全体、active 子集、自由转子路径）逐位一致。余下的开销是每个群元素一次 SVD：批量 SVD 不保证逐位相同，没有做。测试 `test_symmetric_rmsd_is_the_relabelled_aligned_rmsd_without_copies`；变异检验：换回旧循环后变红（多出 24 次 copy）。本机：`make check` 全绿，502 passed / 1 skipped，digest 零差异 332 / 2136 / 4。 |
| `1098bcc8` | 2026-09-28 | **connect 测试台 T0（T-链变体）**（`PLAN_2026-09-28.md` §2）：`calculators.ProtonRelay` / `proton_relay_factory` / `proton_relay_atoms(sides)`。F–C–O–Cl 平面反式锯齿、间距 3.0 Å（重原子互不成键），每个空隙一个质子，2³ 种排布即 8 个极小；R = 全左、P = 全右，沿单跳路径 d = 6→4→2→0。全部项只依赖距离 / 角 / 二面角（平移转动不变）；阱半径 = 0.8 × H–X 键截断（不重叠）。针对草稿测试台审出的四处缺陷设计：固定空间点阱 + tether、阱心在重原子核上、1.5 Å 间距使新键在旧阱内就形成、诱饵 / 桥原子不受力。原型实测：直链有 2 / 8 极小过不了 `confirm_minimum`（`trivial_mode_floor_above_tolerance`），锯齿全过。测试 `tests/test_relay_testbed.py` 9 条：力 = 解析梯度（有限差分 1e-7）、刚体不变、8 个极小全过 `confirm_minimum` 且 `chemical_key` 两两不同、零模恰为 6 个（二面角项由它守，`confirm_minimum` 容忍平扭转所以守不住）、d 数剩余跳数、每个空隙新键在势垒顶之后才形成。变异检验五条全红（去二面角 / 梯度漏顶点项 / 四原子同元素 / 间距 2.2 Å / 直骨架）。**T-诱饵、T-绕路两个变体与鞍点清单未做**。本机：全量 511 passed / 1 skipped、fmt-check / lint 绿、digest 零差异 332 / 2136 / 4 |
| `22d90bd9` | 2026-09-28 | **PERF-22**（`REVIEW_2026-09-28.md`）：`state.encode` 的逐对 Python 双循环改成向量化。hysteresis 因子仍是 1 ± h，数值和原来一样；`np.nonzero` 按行优先扫上三角，边的插入顺序与原循环相同，所以下游看到的集合迭代顺序也不变。先测后改：N=16 时每次 112 → 56 µs，N=22 时 188 → 68 µs，带 hysteresis 的 N=16 为 118 → 75 µs（只快 1.57 倍，没达到改前冻结的「快一倍」；plain 路径达到了）。在 2138 帧上（digest 语料加各 benchmark 的输入），plain、hysteresis、active 子集三种调用下的边集连同迭代顺序与旧实现完全一样。测试 `test_vectorised_encode_is_the_pairwise_rule_in_the_same_order`（以旧的逐对循环为参考实现，hysteresis 带两侧各一个原子对，并包含反向对和子集外的对）；三个变异体（去掉 1−h、1+h 写反、`triu` 取 k=0）全红。本机：`make check` 全绿，512 passed / 1 skipped，digest 零差异 332 / 2136 / 4。 |
| `47803359` | 2026-09-28 | **connect 模式 G3 = 干预 I-b**（`PLAN_2026-09-28.md` §3）：配置 `connect_direction_priority`（默认关；开而无产物则开跑前拒绝）。开时每个微态按它到产物的键差（`ProductReference.key_difference`：在所有 d 最小的产物像上取并集）给方向重排：拉伸 B⁻ 里的键、压缩 B⁺ 里的对进 band −1，`propose` 的分配循环只在有键差时才走 (−1, 0, 1)，否则仍是 (0, 1) —— 这正是草稿的 bug（放进 band −1 而循环只走 0 / 1，想提前的方向一个都选不到）。**只在各族配额内重排、不增不减**：sampled 模式下键差方向仍可能被配额截掉（测试台 R 上 3 条 B⁻ 拉伸只选 2 条），此时 `uncovered_ids` 点名。提案账本记 `connect_priority_selected` 与 `connect_key_difference`，scan 的 `proposal.priority_reason` 记 `connect: B-/B+`。测试 `tests/test_connect_directions.py` 8 条（含 exhaustive 模式全部提前、甲醇三个等价甲基 H 的平局并集）。变异检验五条全红（草稿的 band 循环 / B⁻B⁺ 对调 / 无产物放行 / search 不传键差 / 平局只取第一像）。本机：全量 520 passed / 1 skipped、fmt-check / lint 绿、digest 零差异 332 / 2136 / 4 |
| `5487b5f1` | 2026-09-28 | **connect 模式 G2 = 干预 I-a**（`PLAN_2026-09-28.md` §3）：配置 `connect_queue_priority`（默认关；开而无产物则开跑前拒绝）。出队时取队列里 `d_to_product` 最小的微态、同 d 先入先出（`search._closest_first`）；**不剪枝**，更远的留在队列里。K3 未落笔，所以**无参数**：只按 d 排，相当于 K3 先验 λ→∞、不学习；「随观测更新」没做。每次选择记进 network 的 `connect_queue_order`（选中的微态、d、位置、当时全队列的 d），「每次都取最近」可从记录核对。队列仍按入队顺序落盘，规则在其上确定，续跑复现同样的选择（续跑测试见下一行的修复）。测试 `tests/test_connect_queue.py` 6 条；中途测试台上一次选择是从 [4, 2, 4, 4] 里取了 d = 2 的双跳态。变异检验四条全红（最远先 / 同 d 后入先出 / 只记录不改序 / 无产物放行）。本机：全量 526 passed / 1 skipped、fmt-check / lint 绿、digest 零差异 332 / 2136 / 4。**同时发现主线原有续跑 bug**：已选中但从未交付的方向的停止原因只在内存里，续跑后记成 `unknown`（纯 discover 下在测试台上复现），另作一条干预修 |
| `bb525621` | 2026-09-28 | **续跑丢失方向停止原因（主线原有 bug）**：方向调度在续跑时只从 scan 记录重建，而「已选中、预算截断前从未交付」的方向没有 scan 记录，它的原因（`never_offered`）只在内存里 —— 续跑后一律记成 `unknown`，`directions_dropped_by` 与两级 `direction_ledger` 都错。双阱续跑测试从来碰不到这类预算；connect 测试台的反应物会碰到（`conformer_trials_per_microstate` 截掉最后一个方向），纯 discover 下复现。修法：`resume_point.json` 里加 `unscanned_directions`（这些原因 + 从未交付任何方向的调度的提案），`_rebuild_for_resume` 注 [3] 读回；在 network.json 之外，所以任何运行发布的记录都不变。测试 `tests/test_resume_unscanned.py`：测试台上 discover / connect 两种各杀在第 36 个 trial 再续跑，network.json 与 attempts.jsonl 逐字节等于不中断的孪生运行。变异检验两条全红（不存 / 不读）。本机：全量 528 passed / 1 skipped、fmt-check / lint 绿、digest 零差异 332 / 2136 / 4 |
| `2cbcaec8` | 2026-09-30 | **connect 模式 G4：R→P 路径，只记录**（`PLAN_2026-09-28.md` §3）：`product.product_paths(network)`，给了产物的运行每次 publish 写 network 顶层 `product_paths`。两种视图**不合并**：`observed` 只走试验自己从源走到目标的有向边（`reactions[].attempts`），`any_evidence` 另加续行（`continuations`）与鞍点双侧下降的无向边（`ts_connected_channels`）；路径上每条边各自列出 attempts / continuations / ts。从 depth-0 根节点广度优先、邻居按节点 id 排序，同一网络同一路径；根节点不唯一则拒绝。没连上也是答案：`closest` 给 d 最小的节点。**终止条件（到达 P 即停）未做**，等 G4 落笔，仍跑满预算。schema：product 运行必须有 `product_paths`。测试 `tests/test_product_paths.py` 8 条（7 条纯网络 + 双阱端到端）。变异检验五条全红（两视图合并 / TS 边单向 / 平局反序 / publish 不写 / closest 取最远）。本机：全量 536 passed / 1 skipped、fmt-check / lint 绿、digest 零差异 332 / 2136 / 4 |
| `ebff4ec0` | 2026-09-30 | **connect 测试台 T0 的 T-绕路变体**：`calculators.ProtonRelayGated` / `proton_relay_gated_factory` / `proton_relay_gated_atoms(sides, gatekeeper_home)`。导线原样（先把 `ProtonRelay` 的三个内嵌项抽成模块级 `_spring` / `_angle` / `_dihedral` / `_well`，导线的能量与力**逐位不变**，对 8 个随机几何核过），再加 Br（接在 Cl 后面续锯齿）与第四个质子「看门者」（家在 Cl、外阱在 Br）；唯一的耦合项在看门者在家时把质子 5 穿越 C–O 空隙中点的势垒抬高 A = 1.5 eV。刚性扫描：门关 2.37 eV、门开 0.90 eV、看门者出 / 回 0.60–0.70 eV；最低势垒路径出去 → 质子 5 过 → 回来，**d 沿路 6, 4, 6, 4, 2, 0**，最高垒 0.92 eV。**第一版把看门者的家放在 C，d 沿路 6, 4, 4, 4, 2, 0 —— 绕路是 d 中性的**：反应物里看门者与质子 5 是同一个 C 上的两个 H，Aut(G_R) 交换它们，G1 的对称感知 d 把「一个 H 离开 C」也读成质子 5 的动作；家改到反应物里没有质子的 Cl，群平凡（阶 1）。测试 `tests/test_relay_detour.py` 7 条：力、刚体不变、16 个极小全过 `confirm_minimum` 且 key 两两不同、零模恰 6 个、门只在看门者在家时关、最便宜路径上 d 有上升且比强行过门便宜、看门者在外的态能量都更高。变异检验四条全红（去门 / 去 Br 二面角 / 门梯度漏 Cl 项 / 看门者偏好 Br）。本机：全量 543 passed / 1 skipped、fmt-check / lint 绿、digest 零差异 332 / 2136 / 4 |

改 `src/prrs/` 前先跑 `make digest`（= `python docs/experiments/digest_regression.py`）
（332 文件 / 2136 帧的 digest 黄金样本，不匹配即非零退出）。
它把差异**分类**：`graph_invariant`（图变了，P1/P2 冻结在上面）/
`parity_only`（只有立体化学 parity 动了）/ `key_only`（组件没变而 key 变了）。
有意改变 digest 的工作用 `--allow parity_only` 之类放行 ——
**但放行只确认差异的形状，不确认新值是对的**，那需要化学判断。

这份文档给接手的人。它不重复细节，只说**现在在哪、什么已冻结、什么还开着、以及这个项目的规矩**。
细节全部在 §7 列出的文档里，每一条结论都能追到一次运行。

---

## 1. 这是什么

**PRRS = Perturbation-Response Reaction Search.** 一个 MLP-native 的反应发现方法：
给分子一个受控扰动 → 无偏 NVE 自由响应 → 淬火 → 看它落到哪个化学态。
它**不**预设反应坐标，**不**预设产物，**不**从 TS 猜起。

核心设计取向（这些是设计决定，不是实现细节，改它们要重新验证一批结论）：

1. **Hessian 是验证极限，不是主观测量。** 有限幅度响应是 PES 上的尺度滤波器，
   只有 $a\to0$ 那个极限才有谱定理。有限幅度的本征值**不得**称为"频率"、
   **不得**报 cm⁻¹。见 `PRRS_CURVATURE_ROUTES.md` §E。
2. **三层身份**，不可互相折叠：
   - **化学态身份** `chemical_key` —— 图 + 碎片 + 配置立体化学，**构象被商掉**
   - **反应事件身份** `reaction_event_key` —— 键差 $B^-=E_0\setminus E_1$、$B^+=E_1\setminus E_0$，
     在 $\mathrm{Aut}(G_0)\cup\mathrm{Aut}(G_1)$ 下正则化、方向商掉
   - **通道** `reaction_channels` —— 允许自环（退化反应就是自环，不是构象变化）
   这一层被折叠回第一层过一次（P1 的退化转移），后果见 `PRRS_STATUS.md` §8.19.5。
3. **失败必须关闭（fail closed）。** 判不出来就拒绝并说明理由，绝不用默认值替代。
4. **有向观测与结构证据分开记。** `attempts` = 这次试验自己从源走到了目标；
   `continuations` = probe 被继续交付后到达；`ts_support` = 鞍点两侧下降（**无向**）。
   三者强度不同，网络里分三个字段。

---

## 2. 项目的规矩（这些是用户定的，不是建议）

**接手之前先读这一节。**违反它们产出的结果会被作废，本项目已经作废过一轮。

- **判据跑前冻结，逐字不改。** 要改必须作为**显式修订**记录：旧文本原样保留、
  改的理由、以及**修订发生在看到哪些数据之后**。已有两次修订（B1/B3、A4），
  格式照抄 `P2_OXOBUTANAL_ACCEPTANCE.md` 和 `P2_TEMPERATURE_SCAN_PROTOCOL.md`。
- **不许"边跑边改"。** 绝不为了得到绿色结果而换 seed 或调参。
- **干预必须分别命名**，这样改进能归因到具体某一项。
- **未完成 ≠ 失败。** 分母是**发起数**；未完成记未完成（0 分因为没回答），
  崩溃记崩溃并附原因。`scripts/scan_status.sh` 机械地区分这三态。
- **测量优先于推断。** 本项目有大量"我推断 X，检验杀了它"的记录，
  那些记录是资产不是耻辱 —— 见 §6。
- **旋钮先验证再设计实验。** 在确认被改的变量能推动被测机制之前设计方差实验，
  等于用算力买一个已经注定的答案。这条是 $T=0$ 扫描 8/8 那次教训的产物。
- **门（gate）不能因为它挡住了想要的结果就调。** 要动必须先量清楚它挡掉的是什么。

---

## 3. 已冻结的结论

### P0 乙醇 — 逐字段复现
频率经独立的 `ase.vibrations` 交叉验证（`src/prrs/validation.py`，**搜索侧从不 import 它**）。

### P1 丙二醛（分子内质子转移，退化）
`runs/p1_round8`，**13/13 判据**，与 `runs/p1_round7` 逐字段相同。

| | |
| --- | --- |
| 通道 | `degenerate_reaction`，自环，`event_key = a04410b76e16b117` |
| TS | 4 个，$-3189.6 / -3188.2 / -3187.3 / -3186.7$ cm⁻¹（参考区间 $-3186$） |
| 势垒 | **419.14 meV** |
| `reactions` | **0 —— 这是对的**：退化转移两端是同一化学态，所以是自环不是有向边 |

> **2026-09-04 追加，冻结文字未动**：上面这个 419.14 meV 第一次有了 QM 参照 ——
> `T[E_QM] = 134.4 meV`（ωB97M-V/def2-TZVPD//r2SCAN-3c）。
> **13/13 没有被推翻**：两张面的鞍点位置几乎重合（对称感知 RMSD **0.0151 Å**），
> 所以那 285 meV 的差是**纯能量的、落在模型上**，不是几何、不是方法。详见 §11 与
> `docs/experiments/p1_qm/README.md`。

### P2 3-氧代丁醛烯醇（不对称，跨节点 A→B）
判据由 `docs/experiments/p2_check.py` 机械评定（16 条 + 3 条记录完整性检查）。

| | |
| --- | --- |
| 最终 $T=0$ 轮 | `runs/p2_round7`，**16/16**，实现 `0058158d` |
| 最终 $T=300$ 批 | `runs/p2_T300_seed{16 个}`，单一实现 `0058158d`，16/16 完成，**0 失败** |
| $Q_1$ 有向边落在转移通道上 | **16/16**，95% Wilson $[0.81,1.00]$ |
| $Q_2$ 转移通道存在 | **16/16**，$[0.81,1.00]$ |
| $Q_4$ 有 min-mode TS 支持它 | **15/16**，$[0.72,0.99]$ |
| 全量评分 | 12 个 seed 16/16；3 个 15/16（仅 B1b）；1 个 12/16（seed 71） |
| 虚频 | $-3010.9$ 到 $-2998.8$ cm⁻¹（n = 21） |
| $\Delta E(B-A)$ | $-11.567$ 到 $-10.758$ meV（紧参考 $-11.721$） |
| 势垒 | **+361 meV** |

### 冻结时**明确不成立**的七条（不要引用成结论）

1. $n=16$ 的 16/16 **不是"接近 1"** —— 95% 下界 0.81。区间必须与点估计同行。
2. **B1b 记"不评"**（4 个含 E/Z 通道的 seed）。标签早于分类器修复；只重跑那四个
   会重新制造混合实现批次。
3. $Q_1$ 从 12/16 升到 16/16 **归于 relay 定向是旁证，不是隔离归因** ——
   三处改动同批进入，没做关掉定向的对照运行。
4. **0.81 meV 是本协议下的端点散布，不是能量分辨率** ——
   要成为分辨率必须把所有端点用同一紧容差重新松弛。没做。
5. 直接观测与续行是**两类强度不同的证据**，不可合并。
6. **seed 71 的 C4/D1/D2 是真实失败**（它没有任何 TS），已测明与 ceiling 无关。
7. **openmm 8.6.0.dev 下 backend parity 为何变红仍未查明** ——
   现在绿是因为环境回滚到 8.5.2。

---

## 4. 关键机制，以及它们为什么长这样

每一条都是从一次失败里长出来的。改之前先读对应章节。

| 机制 | 解决什么 | 在哪 |
| --- | --- | --- |
| 解析 Hessian（`hessian_source="auto"`） | MACE 的力本身是 autograd 梯度，所以二阶导已经在图里；比 FD 快 2 倍且到机器精度 | `PRRS_STATUS.md` §2.6b |
| 刚体地板 + 梯度扣除 | $\dot x^\top H\dot x=\omega^2(g\cdot x_\perp)$，转动曲率**正比于残余梯度**，扣掉后残差 $\sim10^{-16}$ | §2.6c / §2.6g |
| 目标模跟踪 | 不给目标方向就爬**全局最软**模：3-氧代丁醛上返回 8 个 $-65\sim-80$ cm⁻¹ 的转子鞍点，而质子转移在 $-3186$ | §8.16 第 1 条 |
| 切向的质量度量 $M^{-1}g_x$ | `internal.gradient` 返回协向量；直接用会给 O 比 H 大 4 倍的权重（实测 H/O = 15.87 = 质量比） | §8.16.1 |
| `flat_biased` 档 | 曲率在地板下**而梯度不为零**的坐标不是平的。一个甲基转子因此冻住 9.97 meV | §8.18 |
| boundary continuation | 鞍点两侧下降证明"两极小相连"，**不**证明"交付的运动到得了产物" | §8.19.2 |
| `arrival_through_saddle` | 一阶鞍点在恰好两个盆地的边界上；两侧同名 ⇒ 该点邻域属于那个态 | §8.20.1 |
| relay 定向对齐 | 本征向量符号任意；不对齐会让下降在 relay 后倒回来 | `P2_TEMPERATURE_SCAN_PROTOCOL.md` |
| `descend_saddle` | 两侧下降的**唯一**实现，`connect_saddle` 与续行共用 | `runner.py` |

---

## 5. 现在的状态与下一步

### P3：苯甲酰丙酮 + 离子化学（**判据仍未冻结**）

**探针已跑过一轮，记录在 `docs/P3_BENZOYLACETONE_PROBES.md`（378 行）。**
那份文档里没有一条是验收结果 —— 它测的是"判据能被写出来之前必须先知道的东西"。
下面只列结论，数字与出处全在那里。

已测掉的：身份层能处理芳香环（Kekulé 对返回恰好 2 个解，环上键正确地不被 lock）；
22 原子成本 42.5 ms/点、解析 Hessian 4.2–4.5 s、~1350 MiB/worker；
**闭合螯合盆地只有 1 个构象**（否掉了"更多软构象"这个预期，`conformer_max_per_node = 8` 够用）；
TS 存在且 order 1，势垒 **+316 meV**、虚频 **−3086 cm⁻¹**，与 P1（419 / −3186）、
P2（361 / −3010）同一带；TS 处 r(O···O) 从 2.5289 收缩到 2.3234 Å（**0.2055 Å**），**重原子门控有实测值**。

**2D 投影已跑通**（`docs/figures/p3_projection.png`，48×48，relaxed 2304/2304 全收敛，
曲面容差与路径同为 2e-3）：双阱 + 中间鞍谷，阱在 r_OO ≈ 2.52、通道底部 2.32–2.36。
**投影间隙最大 +32.9 meV、平均 +1.4，仅 41/490 帧为负**（容差修正前是 375/491）——
按 §0 这就是 (q_PT, r_OO) 的成绩单：位置层面够，能量层面在 TS 附近缺约 33 meV。
**注意该剖面不含势垒顶点**：`descend_saddle` 先走一整步 `irc_step_A` 再淬火，
首帧已比 TS 低 184.9 meV。刚性减松弛 451–4187 meV 处处为正。

**职责边界（用户 2026-09-02 界定）：PRRS 只负责 TS 在哪里、产物是什么。**
势垒和虚频的绝对值是势能面的属性，不是搜索方法的属性，交下游 QM。
三模型对照支持这条界定：**scope 内的量三模型稳健，不稳健的恰好是声明不负责的部分** ——
`chemical_key` / `event_key` / `broken` / `formed` / saddle order **逐字段一致**，
r(O···O)@TS 差 0.036 Å；而**势垒差 3.7 倍**（+316.3 / +84.8 / +96 meV）、
**虚频差 3.3 倍**（−3086 / −926 / −1436 cm⁻¹）。与 SN2 探针同一模式：分歧在过渡态本身。

**两个必须先解决的**（2026-09-04：两条都已处理，见 `benchmarks/b03_benzoylacetone/README.md`；
下面的原文保留不动。**注意 §5 开头那句「判据仍未冻结」指的是本节写作当时** ——
判据已于 2026-09-03 冻结并判定通过，见 `P3_BENZOYLACETONE_ACCEPTANCE.md`；
现在待办的是**标准化收编**：尺度指标要写成带分母/阈值/pass-fail 的逐条判据，
以及补一份机械 `check.py`，P3 的判定至今是人读出来的）：
- **双侧连接只在 POLAR-1 上拿到。** `connects_two_nodes = True`，鞍点两侧分别落回
  A（`8710064d…`）和 B（`b1a4338a…`）—— P3 第一次拿到完整的产物身份证据。
  但在 **OFF24（P0–P2 的参照模型）上 `side +1` 仍然到不了 A**，
  即使 `quench_steps = 5000`，所以不是预算问题。
  跨模型对照指向那个近自由甲基转子：**OFF24 在 A 处的最软内坐标模 0.00262，
  POLAR-1 是 0.00470**，软将近一倍。所以这个失败是**那个势能面的属性**，不是方法的 ——
  但两个模型差的不只是这一个数，仍非隔离归因。
- **必须是 Z 异构体。** E 异构体在 300 个构象下 r(O···O) 全在 4.3–4.7 Å，**根本不能螯合**，
  而它会返回一整套自洽全绿的数字 —— 身份层看图，E 和 Z 的图一样。
  **冻结的 preflight 必须含一条几何门 `r(O···O)`。**

  > **2026-09-04 已建**：`benchmarks/b03_benzoylacetone/preflight.py`。
  > 实测 E **0/300** closed（r 4.396–4.559 Å）、Z 169/300（2.328–4.073 Å）、
  > 两个源结构 2.5193 / 2.5287 Å。数值与上面原记录的差别来自 ETKDG 种子与版本，结论一致。
  > 那条 `side +1` 的失败**记为 MACE-OFF24 这张势能面的已知限制**，
  > 保留"两个模型差的不只是这一个数，仍非隔离归因"这句限定；
  > 要重开需要在同一模型上单独改那个转子软度的隔离实验。

用户定的三类指标（跑前要冻结）：
- **正确性**：仍看节点 / 通道 / TS / 有向边
- **尺度**：单独看候选数、方向 inclusion、Hessian 成本
- **不再要求所有 seed 16/16**

**Track A（中性闭壳，现有管道，MACE-OFF24）**
1. 苯甲酰丙酮质子转移 —— 不对称 β-二酮 + 苯环 + 更多软构象 + 芳香图 + 竞争异构体
2. Diels–Alder（丁二烯 + 乙烯）—— **碎片数 2→1**，PRRS 至今只做过单分子内重排

**Track B（离子化学，需要 charge-aware 模型）**
SN2 → SN1 → E2 → EA。**MACE-OFF24 做不了**：它对总电荷完全不响应
（实测 charge 0/−1/+1 逐位相同，模型里没有任何 charge/spin 模块），
而域门会按构造拒绝每一个离子端点。

### 2D 势能面 + 路径叠加（用户要的可视化方向）

坐标对建议：苯甲酰丙酮用 $(q_{\rm PT},\ r_{\rm O\cdots O})$（显示重原子门控）、
Diels–Alder 用两条成键 C–C 长 $(r_1,r_2)$（对角线上是同步 TS）。

**成本**（实测 12 原子 52.2 ms/点，OFF24）：刚性 40×40 = 1.4 min；
松弛（~100 calls/点）40×40 ≈ 2.3 h、60×60 ≈ 5.2 h。

**必须处理的方法问题**：PRRS 路径是 $3N$ 维的。投影到 2 个坐标后，
**画位置没问题，画能量不能当成曲面的能量** —— 松弛曲面在同一点的能量是其余
$3N-2$ 个坐标都松弛过的最低值，路径一般高于它。正确画法是两联图：
等高线上画路径位置，旁边画路径自身能量 vs 弧长**并叠加曲面沿同一投影线的能量**。
**两条线之间的间隙本身是诊断量** —— 它回答"选的两个坐标有没有张开这个反应"。

### SN2 结果（`runs/sn2_probe.json`，`docs/experiments/sn2_probe.py`）

气相同一性 SN2，Cl⁻ + CH₃Cl。参考态是**一次计算、总电荷 −1、离子放 20 Å**，
所以全局电荷嵌入的常数偏移在每个数里抵消。D3h 鞍点由**在二参数对称子空间内取能量极小**
得到 —— 在 D3h 几何上所有非对称方向的梯度按对称性自动为零，所以这是完整势能面上的
真驻点，不是近似。**内部判据是 Hessian 必须给出恰好一个负本征值，两个模型都通过。**

| | OMOL-0-4M | POLAR-1-M |
| --- | --- | --- |
| CH₃Cl $r$(C–Cl) | 1.7816 Å | 1.7816 Å |
| 复合物 $r$(C···Cl⁻) | 3.1355 Å | 3.1767 Å |
| **复合物深度** | $-9.57$ kcal/mol | $-9.23$ kcal/mol |
| `confirm_minimum` | True, order 0 | True, order 0 |
| **D3h 鞍点** $r$(C–Cl)×2 | 2.4236 Å | **2.2815 Å** |
| **势垒 vs 参考** | $+13.24$ | $+4.53$ |
| **势垒 vs 复合物** | $+22.81$ | $+13.76$ |
| 虚频 | $-149.5$ cm⁻¹ | $-535.6$ cm⁻¹ |
| Hessian 判定 | first-order saddle | first-order saddle |

**不依赖任何文献值的结论**：两个模型在**复合物上只差 0.34 kcal/mol**，
在**鞍点上差 8.71 kcal/mol**，鞍点键长差 0.14 Å。所以它们的分歧**不在长程静电**
（那在复合物距离上才起作用，而那里它们几乎一致），而在**五配位碳的 TS 本身**。
虚频差 3.6 倍，是第二个独立的判别量。

**依赖文献值的结论（暂不成立）**：我记忆里的气相值是复合物 $-10.5$、
势垒相对参考 $+3.0$（相对复合物约 $13.5$）。按这组数，POLAR-1-M 的
$+13.76$ 只差 0.26 kcal/mol，OMOL-0 的 $+22.81$ 差 9.3。
**但本次会话没有核实这两个数** —— 一次网络检索没有找到 Cl⁻+CH₃Cl 的具体值
（只找到 F⁻+CH₃Cl 的 15.8 / 9.6 / 3.3±0.3 kcal/mol，那是**另一个体系**）。

所以"POLAR-1-M 复现了 SN2 势垒"目前是**条件命题**。要坐实它，读这两篇里的数字：

- Gas-Phase Identity SN2 Reactions of Halide Anions and Methyl Halides,
  JACS 1996 — <https://pubs.acs.org/doi/abs/10.1021/ja9620191>
- Model Identity SN2 Reactions CH3X + X⁻ (X = F, Cl, CN, OH, SH, NH2, PH2),
  J. Phys. Chem. A — <https://pubs.acs.org/doi/10.1021/jp054734f>

两篇都在付费墙后（抓取工具进不去），有机构访问权的人一分钟就能读到
复合物深度、势垒、D3h 键长和虚频。**在那之前不要把 POLAR-1 的势垒当作已验证。**

### 模型现状（`/home/ruigengji/MLP/mace/`，1.1 GB）

| 模型 | 电荷 | 原子电荷输出 | Cl/Br/F EA 偏差 (eV) | 尺寸一致性 | 备注 |
| --- | --- | --- | --- | --- | --- |
| MACE-OFF24_medium | **无响应** | 无 | — | — | P0–P2 全部用它，**不要换** |
| mace-omol-0-xl-4M | 有 | **无** | $-0.021$ / $+0.013$ / $\mathbf{-0.858}$ | $-0.403$ kcal/mol 逐位恒定 | 78.6 ms/点（OFF24 的 1.51×） |
| MACE-POLAR-1-M | 有 | **本轮取不到，见下** | $+0.085$ / $+0.078$ / $+0.119$ | $-0.80$ + 距离依赖项 | 需 `graph_longrange 0.4.0`（已装） |

**`Qs` 在 P3 这一轮取不到**（本地 CPU 与节点 GPU 两侧一致）：只拿到 `dipole`，
`'Qs' in atoms.arrays` 是 `False`。键名 / 触发条件 / 尺寸依赖未查明。
所以本表"有，且正确定域"那一格来自更早的一次测量，**P3 没能复现它**，
而"SN1 只有 POLAR-1 能检验电荷去了哪"这条推论目前**不可用**。
偶极不能替代原子电荷 —— 偶极不变不排除电荷在内部重排。
详见 `P3_BENZOYLACETONE_PROBES.md` §3.14 与 §5.7。

**POLAR-1 的推理**（`model_type` 必传，否则走不到极化场路径）：

```python
calc = MACECalculator(model_paths=".../MACE-POLAR-1-M.model",
                      model_type="PolarMACE", device="cuda", default_dtype="float64")
atoms.info["charge"] = 0        # 总电荷
atoms.info["spin"] = 1          # 总自旋，约定未验，闭壳一律用 1
atoms.info["external_field"] = np.zeros(3)   # V/Å，可省
# 输出：calc.results["dipole"]（Debye）、atoms.arrays["Qs"]（Fukui 均衡后的原子电荷）
```

**已验证的三件（都有数）**
- 电荷完整穿过 PRRS 管道：`snapshot` 用 `atoms.copy()` 保留 `.info`，quench 正常
- POLAR-1 电荷正确定域：20 Å、总电荷 −1 → 氯上 $-0.987$、CH₃Cl 上 $-0.013$
  （**SN1 是异裂，只有 POLAR-1 能被检验电荷去了哪；OMOL-0 不输出原子电荷**）
- 离子–偶极复合物 Cl⁻···CH₃Cl（OMOL-0，完整松弛，`confirm_minimum` 确认极小）：
  **$-9.57$ kcal/mol**，文献约 $-10.5$

**未解决的三件**
- **`info["spin"]` 的约定不明**：OMOL-0 上 H₂ 在 `spin=3` 时比 `spin=1` **还低** 0.9 meV，
  物理上该高约 10 eV。**开壳体系在搞清之前不要碰。**
- **POLAR-1 的长程项机制不明**：`×r²` 单调增长（103→320），与用模型自己偶极算的
  离子–偶极理想值之比从 3.5 涨到 10.8 —— **不是 $1/r^2$**。减掉常数后大致 $1/r^{2.6\text{–}3.3}$。
  它在复合物距离（3–8 Å）上贡献 1–3 kcal/mol，任何结合能都带这个系统量。
- **SN2 已跑通，但文献值未核实。**见下节。

---

## 6. 本轮被检验杀掉的推断（保留，因为它们是方法记录）

不要重复这些。每一条都是"我有一个自洽的解释，然后一个测量把它否了"。

| 我推断的 | 杀掉它的测量 |
| --- | --- |
| 更快的 GPU 会抬高并发天花板 | RTX 5080 单点快一倍，天花板从 3.2× **降到** 1.7× |
| 争抢的是 CPU 线程 | 每进程一线程 + 绑一核：1.70× → 1.72×，**没变** |
| CUDA MPS 能解除序列化 | 交错 A/B：**0.15×**，慢 6.8 倍 |
| 那批 seed 被主机 OOM killer 杀了 | 125 GB 空闲、内核日志无 kill；日志里一直写着 `ValueError` |
| 八个 seed 被外部中断了 | **它们还在跑**。`pgrep` 跑在了共享文件系统但不共享进程表的机器上 |
| $Q_1$ 与 $Q_4$ 会互相交换 | 边界帧 vs 直接观测 $r=-0.891\,[-0.96,-0.71]$ 真的竞争；传导到 TS 的是 $-0.062\,[-0.54,+0.45]$ |
| POLAR-1 的距离漂移是真实离子–偶极物理 | 与理想 $-q\mu/r^2$ 之比从 3.5 涨到 10.8，**不是那个律** |
| 苯甲酰丙酮"更多软构象"，reservoir 容量要加 | 闭合盆地只有 **1** 个构象，80 个确认极小同一 key、跨度 0.003 meV |
| 我那条 SMILES 是螯合烯醇 | rdkit 判 **STEREOE**，300 个构象 r(O···O) 全在 4.3–4.7 Å，**0** 个闭合 |
| A 与 conf265 是相距 1.255 meV 的两个极小 | 紧容差下**收敛到同一点**，ΔE = 0.0000 meV |
| A 的负模是苯环扭转（22% 在环上） | 读错了本征向量（`vec` 是**列**索引且质量加权）。改对后 **99.83% 在甲基** |
| `climb_seed` 是个能推动走法的旋钮 | `seed` 在 `follow_min_mode` 函数体里从未被使用（`inspect` 证实） |
| 那 8 对重复会逐字段相同 | **6/8 相同**；amp 0.50 差 1.2e-7 eV，**amp 0.10 结论完全不同**（GPU 非确定性） |
| `TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD` 是 shell 继承的，pop 掉即可 | **mace 自己设**（`mace/__init__.py:5`），`import mace` 那一刻就设回来 |

**最严重的一次不是判断偏差**：我在一台看不到那些进程的机器上读 `pgrep`，
把"本机不可见"读成"没在跑"，据此建议 `rm -rf` 七个**正在跑**的作业目录。
教训写死在 `scripts/scan_status.sh` 的注释里：**跨主机时 liveness 不可知就不要声称。**

---

## 7. 文件地图

```
src/prrs/                 17 个模块，6008 行
  runner.py               单次试验、淬火、Hessian、min-mode、续行、路径记录
  search.py               ReactionSearch 主循环、通道/边登记、鞍点搜索
  chemistry.py            三层身份、classify_transition、stereo_only_difference
  network.py              Registry：接纳、微观态水库、域门
  perturbations.py        方向提案、优先级带、家族配额
  config.py               全部参数 + 校验（跑前拒绝不合法配置）
  validation.py           独立交叉校验（ase.vibrations）—— 搜索侧从不 import
tests/                    13 个文件，172 passed / 1 skipped / 0 failed
docs/
  HANDOFF.md              本文件
  PRRS_STATUS.md          1978 行，主叙事，§8.x 是每一批落地的干预与其证据
  PRRS_CURVATURE_ROUTES.md 1270 行，路线 A/B/C 与 §E 架构结论
  P1_MALONALDEHYDE_ACCEPTANCE.md   507 行，冻结判据 + 八轮记录
  P2_OXOBUTANAL_ACCEPTANCE.md      882 行，冻结判据 + 七轮 + 两次判据修订
  P2_MULTISEED_PROTOCOL.md         241 行，T=0 扫描（**已作废**，理由：旋钮不动机制）
  P2_TEMPERATURE_SCAN_PROTOCOL.md  593 行，T=300 协议 + 结果 + ceiling 诊断 + P2 冻结
  P3_BENZOYLACETONE_PROBES.md      378 行，P3 探针记录（**判据未冻结，无一条是验收结果**）
  experiments/            40 个脚本，全部只读运行产物或做单点测量
  figures/                5 张图 × SVG/PDF/PNG + make_figures.py（读盘，不算新东西）
scripts/
  p3_node_bundle/              P3 的重活，全部分片 + PER_GPU=8 工作池
    run_all.sh                 ts / quench-ab / models / pes / landscape / merge / status
    _env_guard.py              float64 断言、TF32 关闭、模型按 hash 认、torch.load 静音、hostname
    p3_ts_barrier.py           TS / 势垒 / 虚频 / 两侧下降
    p3_three_model.py          OFF24 / omol-0 / POLAR-1 对照，只报模型内部的差
    p3_pes2d.py                (q_PT, r_OO) 刚性 + 松弛曲面，含路径投影
    p3_landscape_full.py       300 构象松弛 + Hessian 确认
    p3_merge.py                合并分片；**缺片会报并以非零退出**；拒绝混容差
  run_p2_temperature_scan.sh   check / probe / run / provenance / report
  scan_status.sh               status + 进程存活 + trial 数三项交叉判断
  _saturation_probe.py         并发天花板测量（含两块卡的实测值）
  _provenance.py               混合实现批次会 exit 1
runs/                     442 MB；哪些计入结论见 runs/README.md
models/MACE-OFF24_medium.model   随仓库走，sha e5ccf5837f685899
```

## 8. 环境与硬件

```
python 3.12.13   Linux-7.2.0-1-cachyos   /home/ruigengji/miniforge3/envs/openmm_dev
mace-torch 0.3.16   graph_longrange 0.4.0   torch 2.12.0   ase 3.29.0
numpy 2.4.3   openmm 8.5.2   e3nn 0.4.4
本机 GPU  RTX 2080 Ti 11 GB    单点 52.2 ms（12 原子，OFF24）
节点 GPU  RTX 5080   16 GB     单点 23.8 ms
```

**并发实测（`scripts/run_p2_temperature_scan.sh probe`）**

```
RTX 2080 Ti   天花板 3.2×   （7 并发 3.7×、20 并发 3.2× —— 过 7 就不涨）
RTX 5080      天花板 1.7×   （1→24 并发：单点 23.8→328 ms，聚合 1.72×）
                            32 并发在 16 GB 上 OOM
```

**更快的卡天花板更低**，这排除了 GPU 计算是争抢源；绑核 + 单线程毫无变化，
排除了 CPU 线程争抢。剩下的嫌疑是多进程共享一卡的 CUDA 上下文切换，**未测**。
**不要开 MPS**（实测 0.15×）。杠杆是**多卡**，不是单卡堆并发。

## 8.1 P2 的 relay 方向判定：19 条不可判定，350 条已证明稳健

全仓 369 条 `travel_overlap` **全部来自 P2**（P3 一条没有）。
2026-09-03 由另一个 session 提出、本 session 独立复核。

**这件事不用 GPU 就判完了**，方法是换一个提法：不问"重叠有多小"，而问
**"投影误差是否大于重叠本身"**。两种投影的差异完全由质量加权质心位移决定，
那是纯几何量，逐帧存在 `runs/p2_*/paths/ts*/descent.extxyz` 里，不需要任何势能评估。

### 数据是双峰的，而且分界不是人挑的

| | 条数 | |
| --- | --- | --- |
| 下簇 | **19（5.1%）** | max \|overlap\| = **0.027051** |
| 上簇 | **350（94.9%）** | min \|overlap\| = **0.935910** |
| 间隙 | **0 条** | [0.0271, 0.9359]，宽 **0.9089** |

本 session 用了一个不预设阈值的检验：把 369 个值排序后找最大相邻跳变。
**它恰好落在 19/350 的分界，宽 0.9089，而第二大跳变只有 0.0148（小 61 倍）。**
所以这个划分是数据自身的结构，不是选出来的。

### 投影误差完全落在间隙内部

348 条下降轨迹上独立复现新旧两种投影并比较单位方向：

| | 本 session（24676 个 travel 向量） | 对方（49982 个） |
| --- | --- | --- |
| p0 | 0.0043 | 0.0043 |
| p50 | 0.1157 | 0.1029 |
| p100 | **0.3220** | **0.4420** |

**向量计数差约 2 倍、p100 差 37%，但结论对此不敏感**（两组数字都满足下面两条）：

- **350 条：裕度 0.9359 远高于误差上界（0.32 / 0.44）** → 符号不可能被这个缺陷翻转。
  **这一条现在是被证明的，不再是被假设的。**
- **19 条：裕度 ≤ 0.0271 低于误差的 5 分位（0.029 / 0.051）** → 记录下来的符号是由一个
  小于伪影的量决定的。

### "不可判定" ≠ "符号错了"

要真正定下这 19 条的符号需要 relay 的本征向量，**那个没有存盘**，所以那一步才需要势能。
复核脚本：`docs/experiments/p2_relay_overlap_margin.py`（只报告，不 gate）。

**对复算的意义：需要复算的是 19 条而不是 369 条**，`p2_round7` 的 3 条优先
（round7 是 P2 最终 T=0 轮，16/16 判据就是在它上面评的）。P2 冻结在 `0058158d`，
要动得走显式修订流程。

与 §8.16.1（切向的质量度量要用 $M^{-1}g$ 而非 $M^{1/2}g$，实测 H/O = 15.87 恰是质量比）
并排看：**这是同一条线上的第二个质量加权缺陷。**

## 9. 技术债

- ~~**manifest 不记 hostname**，只有内核串。跨主机判断运行状态靠它碰巧能用。
  改它会改变 `implementation_sha256` —— P2 已冻结，**现在可以做了**。~~
  **已做 2026-09-21**：`search.py` 的 manifest 新增 `"host": platform.node()`。
  指纹如期移位（`1a6a297143deb645…` → `ec32bc3df5d2bf14…`），digest 零差异，
  不影响续跑（`_rebuild_for_resume` 比的是 `config` 与 `implementation_sha256`）。
  （2026-09-24 起还比 `input_structure_sha256`，并要求 `resume_point.json`，见 §15。）
  记在 `env/BASELINE.md` 文首。
- **扫描不可续跑。** `run` 拒绝非空目录，一次外部中断丢掉在飞的全部工作。
  要做需重建 registry 状态，而**接纳顺序会影响结果**，做错比不做更糟。
- ~~**这个目录不是 git 仓库。**~~（已作废，见本节下文 2026-09-13 那条）所以"两份实现只差记录层"这类判断**无法独立核验** ——
  P2 的混合实现批次就是因此只能作废重跑。**建议接手后第一件事是 `git init`。**

  2026-09-04 又收到一次证据，而且比 P2 那次更难处理：`make hash --check` 报出
  `runs/` 下 **16 个实现指纹 / 107 个 run**，分布极偏 —— `d7726173` 占 26 个、
  `0058158d` 占 19 个，**其余 14 个实现分摊 62 个 run**，也就是一大半 run 的实现
  只出现过一两次。`scripts/_provenance.py` 只能在**一批之内**发现混批，
  跨批的实现漂移没有任何东西在看。

  同一天的逆 DA 探针把这件事推到了极限：扫描**开始时**树是 `0a7aa55dd2aa556e`，
  **跑完时**是 `23179deadf657452` —— 另一个 session 的注释轮次落在了运行期间。
  这次无害（那一轮只改注释；探针只 import `prrs.calculators` 与 `prrs.state`，
  两个文件在整轮重构中 **AST 逐字段相同**，在该运行窗口内**逐字节相同**），
  但它说明**"批"这个单位本身已经不成立**：漂移可以发生在一次运行的中间。
  一个 sha256 只能说"树变了"，说不出变的是注释还是判据 —— 那正是 git 提供的东西。
- **openmm 8.6.0.dev 下 4 个 backend parity 测试为何变红**，未查明。
- **`torch.load` 的包裹已在 `src/prrs/` 内收敛成一处**（`src/prrs/torch_guard.py`，
  2026-09-03）：单一 wrapper、安装一次、**没有 uninstall 路径**，device-move 改成同一个
  wrapper 上的作用域深度计数，所以退出作用域时无物可还原 ——
  原先那个"finally 静默吃掉另一处"的失效模式**结构上不再存在**。

  **两份实现现在通过一个共享信任表协作**（`_prrs_trusted` 属性），三条规则，
  每一条都是先量到失败才成为规则：

  1. 已有 guard 且暴露共享表 ⇒ **把自己的信任路径加进去，不再重复包裹**。
     `_env_guard` 原先看到 `_prrs_wrapped` 就 early-return 而不贡献任何东西，
     所以 `torch_guard` 先装时 bundle 的信任列表**完全不生效** ——
     而三模型探针要加载 omol-0 与 POLAR-1，那两个 `torch_guard` 没有理由信任过。
     用假 torch 的反事实实测：旧逻辑下两者都拿到 `weights_only=True` 而加载失败，
     **是静默失败不是告警**。
  2. **下面还有 guard 时，绝不强加限制性默认** —— 它可能知道那个文件可信。
     反向的同一条曾是真实的 `UnpicklingError`：外层设了 `weights_only=True`，
     内层按设计看到参数已存在就不再干预，于是一个内层信任表允许的 checkpoint 加载失败。
  3. 暴露 `_prrs_trusted`，让**后装**的 guard 能对自己做规则 1。

  两种安装顺序实测（三个 checkpoint：OFF24 / omol-0-xl-4M / POLAR-1-M，
  `-W error::UserWarning` armed）：

  | 顺序 | 共享表 | 三个 checkpoint | warning |
  | --- | --- | --- | --- |
  | bundle 先装，`torch_guard` 后 | 3 条 | 全部 OK | 零 |
  | **`torch_guard` 先装（信任表空），bundle 后** | 3 条 | 全部 OK | 零 |

  第二种顺序在规则 1 之前会静默失败。

  **`intact()` 一类的自检只能测到失效的一半。** 两个方向：
  *REMOVED*（第三方 finally 还原掉 guard）与 *NESTED*（第三方包在外层、guard 仍在链上工作）。
  比较 `torch.load is <自己的 wrapper>` 对两者都返回 False，后者是假警。
  `_env_guard` 的那个函数已按此更名为 `torch_load_outermost()`，
  env 字段改为 `torch_load_wrapper_outermost`，**False 是"去看一眼"的提示，不是判决**。

- **`tests/test_backend_parity.py` 有一个"静默不运行"的失效模式，比变红危险。**
  该文件顶部有四个 `pytest.importorskip`（`openmm` / `openmmml` / `mace` / `torch`）
  加两道 gate（`models/MACE-OFF24_medium.model` 存在、CUDA 可用）。
  **在只装了 numpy/ase/pytest 的环境里整个模块收集到 0 条**，pytest 照样报"全绿"。

  实测对照（2026-09-03，两个 session 各跑一次同一份代码）：

  | 环境 | 结果 | backend parity |
  | --- | --- | --- |
  | 只有 numpy/ase/pytest | 174 passed, 1 skipped | **0 条被收集** |
  | 全栈 + CUDA + 模型文件 | **180 passed**, 1 skipped | 6 passed, 1 skipped |

  差的恰好是这 6 条。**所以"测试通过"在不同环境下含义不同**，而变红至少看得见。
  交接或 CI 时 pytest 应带 `-rs`（列出 skip 原因），或显式断言收集到的用例数。
  唯一被跳过的第 7 条是 `test_report_backend_throughput`，需 `PRRS_RUN_BENCHMARKS=1`，
  那是有意的。

---

## 10–11. 2026-09-04 两轮记录（正文见 [ROUNDS.md](ROUNDS.md)）

- **10. 2026-09-04 工程化重构（`dba26be3` → `23179dea`，零行为改动）**
- **11. 2026-09-04 QM 参照：冻结结论第一次被外部尺子量过**

两节都是**已结账的**：重构零行为改动，QM 参照没有推翻任何冻结结论
（推翻了什么、量到多少都在 `ROUNDS.md`）。当前状态看 §5，接手看 §12。

## 12. 2026-09-05 交接：三项正确性修复，与接手要先知道的事

树 `69d53501` → **`5b69dfed41c4bdb2`** · `233 passed, 1 skipped` · digest 绿（一条已审查迁移）
· 复跑基线在 `env/BASELINE.md`

### 12.1 接手先读这三份

| | |
| --- | --- |
| `docs/REVIEW_2026-09-05.md` | 用户的独立代码审查，逐项验收条件。**本节是对它前三项的执行。** |
| `docs/DECISIONS_PENDING.md` | 待决事项与迁移计划。**没做完的都在这里，别去别处找。** |
| `env/BASELINE.md` | 指纹、`versions()` 九项、复跑命令 |

### 12.2 修了什么（每项都有变异检验）

**① 验证不完整不再算通过**（`runner._probe_minimum`）

`_curvature_along` 被 gate 拒绝的方向原先被 `continue` 跳过，而
`confirmed = bool(measured) and not negative` —— **四个方向被拒、一个非负，就返回通过**。
反向也错：调用方读 `not confirmed` 后写 `not_a_minimum`，
于是"没测出来"被记成"不是极小"，几何被丢。

现在三态：`complete`（所有方向都读到）/ `confirmed`（且无负值）/
`verdict ∈ {minimum, not_a_minimum, unverified}`。
`unverified` 走 `status="unverified"`、`admission="withheld"`，几何进候选清单。

**② 候选出口统一**（`search.ReactionSearch.run`）

三处 `admit` 里，`execute` 会列候选，而**两侧下降与边界延续直接丢几何**。
现在一个 `withhold(structure, source, attempt, reason, detail, confirmed)`，四处调用，
结构写 `withheld/<序号>_<来源>.extxyz`，`source ∈
{trial_endpoint, saddle_descent, boundary_continuation}`。
**几何与出处先落盘，再决定进不进网络。**

**③ 失败落盘独立于网络快照**（异常处理）

旧顺序是先 `publish()` 再写 manifest。当异常**来自** `publish()`（网络里有非有限值，
`io.dumps_strict` 按设计拒绝）时，处理器再抛一次，manifest 根本写不下去，
`status` 停在 `running` —— 看起来像被信号杀掉。
**T=300 的三个 seed 就是这么死的**，`io.py` 记着现场，但当时只改了错误消息、没改顺序。

现在最小的先写：`failure.json`（纯字符串，含 traceback）→ manifest →
`publish()` 尽力而为，失败则在 `failure.json` 里记 `network_snapshot_failed`。

### 12.3 一条被推翻的分析，接手不要再引

**"b04 没有 bracket 是因为只有最高一级反应、没有下一级可比"——错的。**

```
0.1 → m0000 返回 · 0.3 → m0000 返回 · 0.6 rejected(collateral) · 1.0 → m0001 离开
```

细化只比相邻样本（`search.py`），(0.3,0.6) 与 (0.6,1.0) 都因 0.6 是拒绝点而不成对。
**是 collateral 拒绝隔断了唯一有效的 (0.3, 1.0) 对。**
若 0.6 有效返回，现有实现就能成 bracket —— 不用加档、不用改阶梯。

**所以 bracket 不独立于 collateral**，次序是先修扰动交付再评估扫描策略。

### 12.4 QM 复核这条线已经通了两例

| | |
| --- | --- |
| **P1**（`docs/experiments/p1_qm/`） | 双向 IRC 收敛、两端 0 虚频、正向端点对极小 RMSD 0.0004 Å；`event_key a04410b76e16b117` **复现冻结值** |
| **b04**（`docs/experiments/b04_qm/`） | 从 PRRS 的 TS 出发精化到 2.3645/2.3645，与独立求得的 QM 鞍点小数点后四位相同；双向 IRC 收敛、两端 0 虚频互为镜像；`event_key 847306d3f8c57bd4` 与 PRRS 一致；中心势垒 +13.02 kcal 两组几何相同 |

**IRC 在本项目里第一次真正执行**（§3.7 那条"由执行证实"此前只有文档）。

协议：**r2SCAN-3c 出几何与频率 → ωB97M-V/def2-TZVPD 在同一几何上补单点**。
r2SCAN-3c 单独用把 SN2 中心势垒错 11 kcal，**几何好、能量必须补**。
校准：SN2 上复现文献 0.4–0.9 kcal。

**两例同一形状：几何与事件身份属于方法且正确，能量与曲率属于模型。**

### 12.5 没做完的，按用户定的顺序

```
1 可复跑环境与版本记录、原版回归        ✅ env/BASELINE.md，本轮实跑
2 三项正确性问题                       ✅ 本节
3 单独修扰动交付，比较有效扰动比例与 bracket   ← 下一步
4 扫描策略 · b04 正式验收 · 扩展 B 组
5 效率用「每个经 QM 确认的独立通道成本」衡量
```

**第 3 步的具体建议**（用户）：对适用的桥键按分子图切分片段、整体移动，保持片段内部几何；
环内键与跨碎片接近分别处理。之后把扫描显式分为
**返回 / 离开 / 被阻断 / 未知** —— 不得把拒绝当返回，也不得无记录跳过。

诊断已完成：拉 C–Cl 0.600 Å 时只有 C（0.4482 Å）与 Cl（0.1518 Å）动，
三个 H 移动 **0.0000 Å**，每根 C–H 缩短 0.0628 Å；CH₃ 若随 C 刚体平移则变化 **0.000000 Å**。
**是位移实现，不是必要耦合。** torsion 已用 `rotate_fragment` 避开，bond 没有等价物。
契约已钉成测试：被驱动的键豁免、其余键的损伤仍拒绝。

### 12.6 接手要知道的三条规矩之外的事

- ~~**这个目录仍然不是 git 仓库。**~~ **已作废（2026-09-13 核实：它是了，`main` 上有提交）。**
  当日那条观察本身仍然成立并且仍然要紧：2026-09-04 实测一次扫描的**中途**树哈希变了
  （`0a7aa55d` → `23179dea`），连"批"这个单位都不成立；
  `runs/` 下 16 个实现 / 107 个 run，14 个实现只分摊 62 个。
  **扛住可追溯性的仍然是 `implementation_sha256`，不是 git。**
- **`benchmarks/` 是新建的反应基准套件**，六条规矩在它的 `README.md`。
  b00–b03 已转 `fixture`（判据冻结过、判定通过过，验收职责移交 P4+ 的 QM 协议）。
  ~~**b04 仍缺正式判据与机械评定**（`CRITERIA.md` / `check.py` 是空壳）。~~
  **已作废（2026-09-13）**：b04 的判据与 `config.json` 同日冻结，`check.py` 逐条对得上、
  自检 5 例，冻结后首次运行 `runs/b04_frozen1` 评定为 `observed_direct`（分步侧 `unanswered`），
  状态词已是 `done`。见 §13。
- **默认模型从 P4 起是 `MACE-POLAR-1-M`**（用户 2026-09-04）。
  P0–P3 冻结在 MACE-OFF24 上，不重跑、不改标签。
  依据是两组 QM 对照，不是偏好：POLAR-1 在中性有机质子转移上偏 +12 meV、
  在离子 SN2 中心势垒上偏 +2.1 kcal，是唯一两边都成立的候选。


## 13. 2026-09-13 交接：三路审计的执行、b04 冻结、两条脑洞落地

树 `fa15a829869cdc29` → **`c393da782a1350a1`** · `416 passed, 1 skipped, 9 xfailed`
· digest **零差异**（332 files / 2136 frames / 4 event pairs）· `make lint` 绿

### 13.1 接手先读这三份

| | |
| --- | --- |
| [`SESSION_2026-09-13.md`](SESSION_2026-09-13.md) | **当日记录。** 三路审计的逐条核实（P0/P1/P2 各自成立到什么程度、证据在哪），以及两轮脑洞里**被实测推翻的前提**。这两样别处没有 |
| [`DECISIONS_PENDING.md`](DECISIONS_PENDING.md) | 「未决清单 · 2026-09-13 收口」。**没做完的都在这里** |
| [`DONE_LOG.md`](DONE_LOG.md) | 今天四条落地的细节、证据与运行目录 |

### 13.2 一句话摘要

**修了**：`reaction_event_key` 改在生成群闭包上正则化（并集不是群）·
`rotate_about_vertex` 补 `(0, π)` 守卫（唯一会把错数字写进命运-幅度记录的路径）·
`b04/check.py` 与 `CRITERIA.md` §5 对齐两处 · 规矩 5 与 `TAXONOMY.md` 的
「不存在 ⇒ 协同」改写（走修订记录）。

**建了**：`tests/test_quotient_invariance.py`（六个商点统一的置换不变性 fuzz，**首跑即抓到一条**）·
谱提案（Hessian 已算白扔，投掉平凡模与扭转后最软；b04 上把实际反应的方向从抽签变成 1/28 第一名）·
碎片刚体转动（`internal.rotate_rigid`，第三种精确实现）。

**冻了**：b04 的 `CRITERIA.md` 与 `config.json`，同日。首次评定通过。

### 13.3 接手要先知道的三件事

- **`reaction_event_key` 只在同一套原子编号内是身份。** fuzzer 首跑抓到的，已钉成
  `xfail(strict=True)`。一次 run 内共享编号所以**已记录结果不受影响**，坏的是跨 run 比 key。
  修法会改掉历史上每一个 `event_key` 摘要 —— 这是此刻**唯一挡着队列**的决定。
- **B（碎片刚体转动）关闭了一个真实的定义域缺口，但没解决 b05。** 我预测它会，实测否掉了。
  `quench_steps` / `irc_relay_hops` 那条原诊断依然成立。
- **`frozen1` vs `probe4` 不是受控对照**（ts 1→2、reactive 2→4）。两次运行实现树不同，
  且同一轮动了七个文件。要隔离谱提案的效果需要一次只切 `spectral_rank_modes` 的 A/B，**没做**。

---

## 14. 2026-09-20 交接：七决策落笔，第 4 件与 P1-4 落地

树 `c393da782a1350a1` → **`b4e92642c66c9bab`** · `420 passed, 7 skipped, 9 xfailed`
（yayoi；skip 全是环境性，本机无 MACE checkpoint / CUDA）· digest **零差异** · lint 绿

### 14.1 接手先读这三份

| | |
| --- | --- |
| [`PLAN_2026-09-20.md`](PLAN_2026-09-20.md) | **六件大事的实施计划**（现状带行号、验收、量） |
| [`DECISIONS_PENDING.md`](DECISIONS_PENDING.md) | 「未决清单 · 2026-09-20 收口」——七决策的正文与队列状态 |
| [`DONE_LOG.md`](DONE_LOG.md) | 当日三条的细节、漂移清单、两次变异检验 |

### 14.2 一句话摘要

**落笔了**：PLAN 文末的七个决策点（用户授权按推荐项落笔），逐条见 DECISIONS_PENDING
2026-09-20 节。**从此「不落笔不开写」的门是开的**，队列顺序 4→1→2→5.0→5→6→3，
下一件是 1（正则标号，唯一挡队列的）。

**建了**：`src/prrs/schema.py`——schema 2 首次被定义（此前 `schema_version: 2` 盖在一个
未定义对象上，`runs/` 70 个文件实测两种旧版本、五档字段漂移）。锚测试把整个 `runs/`
钉在「清单内的偏差逐字复现、清单外的零问题」上；漂移清单里 65 个文件的偏差全部可用
写入年代解释。`publish()` 末尾记录 `schema_warnings`（干净时不写该键，文件逐字节不变），
CLI `--strict-schema` 才把非空当失败。

**修了**：P1-4（决策 7）——`admit` 建节点遇自同构溢出改拒 `automorphism_budget`，
不再拿 `[identity]` 当群去量 RMSD。语料里该分支从未执行，digest 零差异是证据。

### 14.3 接手要先知道的三件事

- **决策已全部落笔但决策 4（证书命题）是草案**，显式批准前不发任何证书；它只挡第 3 件
  （排最后），不挡队列。
- **B′ 挂起是有触发条件的**（决策 6）：5.0 放开 `quench_steps`/`irc_relay_hops` 重跑之后
  若负曲率仍在精确周期坐标上，才单独命名启用。先跑 5.0，别直接写 B′。
- **committor（决策 3）是独立探针**：T=300 K、N=16/窗口、只进模型报告不进判据、
  默认 `committor_shots = 0`。分辨率限定（只答"接近 0/1"）必须开跑前写进协议。

### 14.4 下午追加（指纹移到 b1c4d440，正文见 DONE_LOG）

六件里**本机可做的全部做完**：第 1 件全件（引擎由 flash 写、验收全过）、第 2 件机器+接线、
第 3 件机器+解析势链条（恒盖 DRAFT 章等决策 4 批准）、第 5 件算子/分流/惰性枚举/测量脚本、
第 6 件 6.1+6.2。**排队给 kasuga02（MACE 机）**：5.0 重跑、平动阈值测量、committor 实跑。
**另起会话做**：6.3–6.5 resume 接线（运行态重建清单在 DONE_LOG）。GLM-5.3-Flash 配额事件
与"flash 写+主会话门"分工的结论如实记录在 DONE_LOG 末节。

## 15. 2026-09-24/25 交接：外部审查核对、续跑回滚重做

**下一步看 [`PLAN_2026-09-25.md`](PLAN_2026-09-25.md)**。A0（零 MD 离线重分析）**已做**；A1 **已落笔、协议已冻结**（2026-09-25）（冻结时 sha256 `b685120fe25665a4…`）；落笔值：k = 160 步（40 fs）、L_max_run = 1280 步（320 fs）、V3–V5 阈值 10% / 5% / 20%、键图判归 + `q_PT` 交叉检查。

树 `4014c4b2`（按本文件记载，会话开始时未现算）→ **`ae47c4d9142203dc`** ·
`476 passed, 1 skipped` · digest **零差异** 332 / 2136 / 4 · fmt + lint 绿 · 变异检验五次全红。
先在 kasuga01 本机跑（当时会话环境里 ssh kasuga02 被拒，用户一次性授权），会话重开后在 **kasuga02 复跑，结果一致**。

- **正文**：[`REVIEW_2026-09-24.md`](REVIEW_2026-09-24.md) —— 外部审查逐条核对（数学两条成立，
  不波及现有 TS；committor 那条只是键名）、修了什么、committor 规划器草稿与 PRRS-2 原型对照现码、
  轨迹终止规则（固定 50 fs + quench vs 首次进入即停）、计划 F0–F5、**七条待落笔**。
- **最要紧的一条**：PLAN 第 6 件（续跑）09-21 的验收只在 kill@6 测，那恰是微态之间的缝。
  scan 中途 kill 的续跑是坏的（记录有洞撞目录 + 在飞微态被丢）。已改快照回滚，验收参数化
  kill@4 / kill@6。**代价：kill 会重做在飞的那个微态。**之前的 run 不能跨指纹续跑，本来也不能。
- **没动的**：`PRRS_CURVATURE_ROUTES.md` 只在 B.7 与 C 节加勘误引用块，原文一字未删；
  `p_B` 键名、PRRS-2 新增部分都等落笔。
- **t_free 敏感性已量**（`runs/p1_tfree_sensitivity`，R7）：P1 committor 的 p_B 随自由段长度动 —— scan 8 在 40 / 80 / 160 fs 为 0.06 / 0.57 / 0.67，轨迹在 160 fs 仍回越。**「scan 8 接近 0」不可再引**；R6（终止规则）变成必答。
- **首次进入第 0 段已做**（`runs/p1_first_entry_stage0`，kasuga02，零 MD）：k ≤ 20 步是假进入（质子在脊上来回），k = 40 仍不稳，**k ≥ 80 步是平台**（k = 80 / 120 / 160 逐发判归相同，再离开 ≤ 1/16，quench 核对 0 例不一致）。回越**不是无限的**：待满 20–40 fs 之后几乎不再离开，固定 160 fs 与平台接近 —— 此前「孤立 NVE 会一直回越」的说法已更正。只用于选参数，不可引为 p_B。
- **首次进入确认批跑完并判定（2026-09-26）**：V1 过，scan 0 / 2 / 6 / 8 四窗口 V2–V5 全过，`p_B_given_AB` = **0.76 [0.63, 0.86] / 0.12 [0.06, 0.24] / 1.00 [0.93, 1.00] / 0.72 [0.58, 0.83]**（N = 50），已替换撤回的 40 fs 值（`P1_COMMITTOR_PROTOCOL.md` 末节）。主线 A 剩 A5。原记录：**协议已冻结，确认批（A3）在跑**：**2026-09-25 19:5x 起在 kasuga02 跑**（`setsid nohup`，关窗口不影响）：`docs/experiments/p1_first_entry_confirm.py` → `runs/p1_first_entry_confirm/`，日志 `.scratch/first_entry/confirm.log`。每步存帧使每步约 0.2 s（`sample_interval = 20` 时约 0.07 s），实测外推约 **15 h**（协议里「约 5 h」是冻结前的估算，不是判据）。被杀了就原命令重跑，已完成的发跳过、只丢在飞的一发。跑完后 `python docs/experiments/p1_first_entry_confirm.py --analyse` 出 `verdict.json`（V1 整批 + V2–V5 逐窗口）。
- **A5 完成（2026-09-26）**：`run_trial(first_entry=…)` 逐步计数、连续 k 步在 A / B 图上即停，之后原有 quench 就是入口核对；`shoot` 透传；配置 `committor_entry_steps`（k，0 = 关，默认）/ `committor_max_steps`（1280）；search 的 committor 分支按进入侧判，核对不一致只计数不改判，没有 B 或两端同图时退回固定规则并在 `termination` 写明；R1：`p_B` → `p_B_given_AB`。 kasuga02：`make check` 绿、digest 零差异 332 / 2136 / 4；变异检验两条都红（`>=`→`>` 进入步差一；删 `break` 跑满上限）；**真实数据核对** `runs/p1_first_entry_production_check`：生产路径重打 A3 的 200 发（同帧同 seed、提前停），**进入侧与进入步 200/200 相同、入口 quench 落点 200/200 相同**，约 23 min。 指纹 → `b601b124`。**主线 A 全部完成**；下一步见 `PLAN_2026-09-25.md` 支线 B / C。
- **支线 B0 / C2 / C3 已跑（2026-09-26）**：搜索的 bracket 也随 t_free 移动（→ 新决策 R8）；窗口普查显示最大缺口是 unresolved（154）而非窗口泛化（+24）；沿鞍点不稳定模 p_B 的变化全在 |q_PT| < 0.1（探索批，不可引用）。另发现 D1：鞍点下降 `descent.extxyz` 只剩 minus 分支。
- **C1 + D1 已改（2026-09-26，指纹 → `99edddde`）**：D1：鞍点双侧下降改用**一个** `PathRecorder` 管两支，帧计数连续，`descent.extxyz` 不再被第二支覆盖（此前历来只剩 minus 支的几何）。 C1：scan 每个样本加 `outcome`（returned / departed / ridge / transient / unresolved）与 `outcome_reason`（未定时的失败码或 `admission:…`），只记录不路由。 变异检验两条红（yayoi）；未变异全量用户报告通过。
- **所有改动未提交**，与工作区原有的大量未提交改动混在一起；git 由用户自己管。


## 16. 2026-10-01 交接：PARS 独立模块落地（MVP-1..5）

**新模块在 `pars/`，`src/prrs` 一行未动**。依据根目录《Pareto-Adaptive Reaction Path Search (PARS) — MLPATH Integration Design.md》（§6/§25/§26/§31/§43）实现 MVP-1..5：候选表示 → 静态 frontier → margin/archive → 机制感知 cover → ASE NEB 适配器。QM 反馈实验（MVP-6）未做，选择记录里已留 audit 槽位（§20）。

- **独立性契约**：`pars/src/pars` 内禁止 `import prrs`（AST 边界测试钉死）；与 PRRS 的唯一集成面是**文件契约**——run 目录的 `network.json`（`reaction_channels` 即 channel 身份）、`config.json` manifest、`attempts.jsonl`、`trials/*/{response.traj,observations.jsonl,result.json,endpoint.extxyz}`。读取器唯一入口 `pars/src/pars/adapters/prrs_run.py`。旧格式 run（无 reaction_channels）fail-closed 拒收。整个 `pars/` 文件夹可整体摘出独立成库（自带 pyproject/Makefile/README/LICENSE，依赖仅 numpy+ase）。
- **权威留在 PRRS 侧**：channel/对称/键感知直接消费产物，不重写；不确定性只用已落盘的逐帧 `uncertainty_eV_A`；缺失目标量记 None 不造假（§7），全通道缺失时该维度整体退出归一化空间。
- **关键语义决策**（如与设计字面有出入，以此为准）：(a) CERTAINLY_DOMINATED 需要**真实存在的不确定性数据**（或有配置容差）才可判——σ 未知不当 0 用，否则 margin 机制死代码、每次硬删都"伪确定"；(b) `dominance_mode: nominal` = 严格 Pareto 剪枝（§38 消融的对照组）；(c) NEB 分两段：普通 NEB 收敛后再开 CI，CI 段不收敛/势垒反升/最高像不在内侧 → 判病理 band，回退报告普通 NEB 结果并写明（§18 监控要求）；(d) 事件窗口按产物里权威的 `phase` 标签裁剪（source/free_start/free），强制帧丢弃但保留在 provenance（含 `kept_frames` 原始帧号，窗口不连续）。
- **测试**：`pars/` 自有 78 个测试全绿（openmm_dev，`cd pars && make test`），含设计 §31 不变量逐条（排列不变、尺度不变、duplicate stability、channel isolation、γ rescue 单调、cover 确定性）+ 边界测试 + 解析势上的精修契约（含 UNRESOLVED_MLP 门与 backend 失败报告）。测试不依赖 `runs/`，合成 run 目录 + 自建解析势（triple well / per-atom double well）。
- **门禁**：父仓库 `make check` 绿：fmt + lint + **537 passed, 7 skipped** + **digest 零差异** 332 / 2136 / 4（2026-10-01 实测，PARS 不在父 lint/test 路径内，属预期）。
- **真实 run 冒烟（只读）**：`pars select --run-dir runs/connect_g5_explore_detour/seed1/pIapIb` → 120 attempts → 16 候选 → 5 通道 → **5 条代表路径**（压缩 16:1），每通道带冷存档与 audit 样本；代表种子 extxyz（9 图）可直接送 NEB。`pars refine` 已在解析势上端到端跑通；真实体系精修需 `--calculator-factory` 给一个包一层 MACE 的工厂（运行时 importlib 解析，静态依赖仍为零）。
- **下一步**（等用户决策）：① MACE 工厂包装与真实体系 refine 试跑；② MVP-6 QM 反馈实验协议（ranking inversion / margin rescue 度量）；③ 是否给 pars 单独建 git 仓库（现与父仓库同样未提交，git 由用户自管）。

## 17. 2026-10-03 交接：路径成熟层独立模块（骨架 + M0 + M1）

**下一步看 [`PLAN_GEN_COMPAS_MATURATION_2026-10-03.md`](PLAN_GEN_COMPAS_MATURATION_2026-10-03.md)（§9 修订 + §4 M0–M6）与 [`docs/README.md`](README.md) 路径成熟层行**。正文见 `DONE_LOG.md` 末节。

- **新模块 [`mature/`](../mature/README.md)**（用户决定：独立模块形态，与 `pars/` 同级；修订记录在施工计划 §9）：单向 `import prrs`（复用 `runner.shoot` 与只读 `Registry.identify`，不复制动力学）；`src/prrs` 与 `pars/src/pars` 禁止反向 import，AST 边界测试钉死并做过变异检验。**不碰 `src/prrs/`，指纹仍 `ebff4ec0`**（实测 `make hash` 与 preflight 记录一致）。
- **M0 已做（只读）**：`docs/experiments/maturation_preflight.py` → `docs/experiments/maturation_preflight/preflight.json`。四个指定 run 全 supported（探索批标 reference-only）；P1 区域 = 无序对 {m0000, m0001}，chemical_key 相同而键图不同 → 首次进入 eligible；27/27 精化记录同一鞍点（0.419128–0.419208 eV）；PARS selection 缺失，M2 前须重 select。基线：537 passed / 7 skipped + digest 零差异 332/2136/4（存 `baseline_2026-10-03.json`）。
- **M1 已做**：`mature/src/mature/schema.py`（协议/区域/提案/射击/统计/诊断全件，同图区域定义期即拒、pulse 拒、漂移容差只有零号测量批可留空、A/B=0 时 p_B 为 null、词表无"已验证 TS"）+ `io.py`（原子写、提交标记定完成、同配置任务表可重建、换协议/内容/指纹拒绝续跑、半写隔离不删除）。`mature/` 自有 53 测试全绿（含已发布 Wilson 值的两条数字锚）。
- **未做（下一施工块 M2+M3）**：候选提取（等 pars select 重跑）、射击接线（`shoot(..., first_entry=...)` + `Registry.from_network` 离线判归已勘察可行）、CLI 实装。M4 前不做 console script 的正式解析方式决定。
