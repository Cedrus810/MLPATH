# 交接：二分中点卡住（TODO 第 9 条）工作线（2026-11-01 至 11-06）

写给新开的窗口。本文件接替 [HANDOFF_HF_MIDPOINT_2026-11-04.md](HANDOFF_HF_MIDPOINT_2026-11-04.md)，补上 11-04 之后的 v5（精阶段 FIRE）、候选身份和 HC。只汇总现状和去处，数字以链接的原文为准。更早的 b05 HF D2 交接见 [HANDOFF_HF_D2_2026-11-01.md](HANDOFF_HF_D2_2026-11-01.md)。

## 必读（开工前全部读完）

- `docs/TS_HESSIAN_FREE_DESIGN_RECORD_2026-10-24.md`：设计稿；§5.1 反力零点二分、§5.3"以实测分辨率决定终止"是这条线的依据。
- `docs/HF_MIDPOINT_THRESHOLD_OPTIONS_2026-11-02.md`：问题、代码里中点阈值管什么、离线重放、选项 O1–O5。
- `docs/PLAN_HF_MIDPOINT_COARSE_2026-11-02.md`：用户 11-02 对 O3 的原文决定、界线 B 的测量与决定（B = 0.025）、v4 验收。
- `docs/PLAN_HF_MIDPOINT_O5_2026-11-02.md`：优化器对照（FIRE / BFGS）的设计与结果。
- `docs/PLAN_HF_MIDPOINT_FINE_FIRE_2026-11-04.md`：用户 11-04 的原文决定（只换精阶段、状态交接、验收方式、之后的分流），v5 的实现与验收结果。
- `docs/DONE_LOG.md` 2026-11-01 至 11-06 各节：每一步的命令、数字、更正。
- `docs/TODO.md` 第 9 条及"未确认"一节：进展与仍未验证的推断。
- 全局与项目 `CLAUDE.md`：NFS 缓存（跨主机结果在写入那台读）、PBS 规矩、只按 PID 杀进程。

## 1. 问题

二分中点的约束松弛要求横向力 ≤ 1e-3 eV/Å，200 步内达不到，HS 就以 `unresolved:transverse_not_converged` 结束。b05 八块 HF 运行里，268 次中点松弛失败 71 次，来自 22 个意图。新 E2 鞍点（QM 确认，比已知 E2 TS 低 11.75 kcal/mol）就藏在这类项里。

## 2. 用户决定（原文，按日期；全文见各计划）

- 2026-11-02（O3，`PLAN_HF_MIDPOINT_COARSE_2026-11-02.md`）："粗松弛先取 2e-3；两侧反力同号、且各自都超过实测的符号分辨界线，才允许提前收括区。其余继续到 1e-3。""保留几何、约束残差和双侧反力差检查。""候选的完整力、横向力判据保持现行值。若剩下的失败仍集中在近根松弛，再做 O5 的优化器对照。"
- 2026-11-02：问"B 取 0.025，可以吗"，回答「可以」。
- 2026-11-04（v5，`PLAN_HF_MIDPOINT_FINE_FIRE_2026-11-04.md`）："先把 FIRE 只接进中点精松弛，再验收完整 HS；粗阶段和 D2 留到后面。""这轮之后再分流：若剩余失败转到粗阶段，再研究粗、精都换 FIRE，并核查粗阶段反力分辨率；若仍卡在精阶段，就集中诊断剩余帧。已有两个救回候选也应先比较几何身份，再做 HC——多个意图可能找到同一个 TS，新增候选和新增事件要分别计数。""下一步最有价值的问题是：FIRE 能否把完整搜索稳定地送到候选。 软模与噪声的机制解释可以随后补"

## 3. 已查清的事（按时间）

| 步骤 | 结论 | 记录 |
| --- | --- | --- |
| 按阶段诊断 | 不是噪声底的平台：失败段末 50 步横向力约为阈值 10 倍，步长贴近 0.05 Å 上限 | DONE_LOG 11-01 |
| 同帧冷启动对照 | 不是带入曲率的问题 | DONE_LOG 11-01 |
| 横向曲率（PBS yayoi35） | 失败帧有一个单独的软横向模（0.013–0.024 eV/Å²），游走中位 96% 沿此模 | DONE_LOG 11-02 |
| 离线重放 | 71 个失败中点中 65 个的括区方向在 2e-3 下已定 | 选项文档 §3 |
| 界线测量 | σ_c 中位 0.0040、p90 0.0073，未见一致偏差；用户定 B = 0.025 | COARSE 计划 §2 |
| v4 验收 | 对照 15 项不变；失败 22 项救回 2；剩余卡点在近根精阶段 | COARSE 计划 §5 |
| O5 对照 | 同帧冷启动：FIRE 近根卡点 33 / 38（BFGS 3 / 38），收敛反力与 BFGS 差 ≤ 0.004 | O5 计划 §6 |
| **v5 验收**（11-04） | 失败 22 项中首轮 14 项出候选（v4 1 项），复核中 v5 12 / 13 项 3 / 3 次出候选；对照 15 项候选全部保留（ΔE ≤ 0.002 eV、RMSD ≤ 0.0046 Å）；forward −8%，墙钟 +36%（外部负载在变，不归因） | FINE_FIRE 计划 §2；DONE_LOG 11-04 |
| **候选身份**（11-06，只读） | 两次验收的 104 个候选全部落在 22 个已连通事件的 3 类里；失败组 44 个（15 意图）全部是新 E2 鞍点（对称 RMSD ≤ 0.038 Å） | DONE_LOG 11-06；`results/2026-11-06/candidate_identity/` |
| **HC**（11-06） | v5 重跑 22 个失败意图：15 出候选，13 HC 连通且化学对全部同新 E2 事件，2 个 `endpoint_unstable`（高能侧没离开鞍点）；**新增候选 15 意图、新增连通 13、新增事件 0** | DONE_LOG 11-06；`results/2026-11-06/midpoint_fire_hc/` |

## 4. 代码状态（`explore/src/explore/ts_ef.py`，均默认关闭，未接入 D2）

| 规则 | 开关 | 内容 |
| --- | --- | --- |
| v3（现行默认） | — | 单阶段中点，1e-3 |
| v4 | `midpoint_coarse_transverse_fmax_eV_A` + `midpoint_sign_bound_eV_A`（两者同给） | 两侧先松弛到粗阈值，几何、反力差门、候选让位、同号且超界线全过才按符号收括区，否则继续到 1e-3 走 v3 逻辑 |
| v5 | v4 加上 `midpoint_fine_optimizer="fire"` | 只有中点精阶段用 FIRE（从粗阶段结束处开始，速度归零，不用带入曲率；接受后保存几何、清 BFGS 曲率；失败不改 state）；粗阶段、括区端点、延续仍为 BFGS |

- D2 / `HFSettings` / manifest / `_ef_settings` 都没改。所有 v4 / v5 运行都是在 benchmark 脚本里用 `midpoint_coarse_validate.settings_for` 构造的设置。
- 测试：`test_ts_ef test_ts_branch test_ts_connect test_ts_warm_start test_d2_hf test_d2 test_d2_resume` 108 passed（11-04）；v4、v5 各有变异检验；**全包 `make check` 一直没跑**。

## 5. 待办（方向需用户定）

| 编号 | 内容 | 状态 |
| --- | --- | --- |
| 9c | v5 下仍卡在中点精阶段的项（11-04：`k0·8-14`、`k1·8-14`、`k2·6-8`、`k2·8-14`；11-06 重跑：3 项，组成有 fp32 波动），以及 O5 中 FIRE 未收敛的帧 | 诊断已完成（11-06）：精阶段 1e-3 在这些帧上低于实测噪声水平；改法选项见 `docs/HF_MIDPOINT_FINE_STAGE_OPTIONS_2026-11-06.md`；P1 用户已定，施工为 v6（`docs/PLAN_HF_MIDPOINT_COARSE_CANDIDATE_2026-11-06.md`），11-07 验收通过（精阶段失败 8 → 1，对照候选保留），未接入 D2；Q2 未定；见 DONE_LOG 11-06 末节"已完成：v5 剩余精阶段卡点的诊断"。用户 11-04："若仍卡在精阶段，就集中诊断剩余帧" |
| 9c' | HC `endpoint_unstable` 的 2 项（`k2·2-4`、`k2·6-8`：高能侧终点能量等于鞍点） | 未查 |
| 9c'' | 延续阶段 `walk_relax_unconverged`（每次 3–4 项，v4 / v5 相同代码，fp32 波动下时有时无） | 已普查（11-07，只读）：与 `hysteresis` 同属两侧都无最大值，不是数值收敛问题；另 `no_resolved_maximum_bracket` 的端点帧已满足候选判据，用户定 R1，施工为 v7（`docs/PLAN_HF_ENDPOINT_CANDIDATE_2026-11-07.md`），定向冒烟通过；见 DONE_LOG 11-07 |
| 9d | v4 / v5 接入 D2：`HFSettings` 字段、新规则版本、b05 整跑 | 11-07 已施工为 D2-HF-v1.3（v1.2 + v7，`docs/PLAN_HF_D2_V13_2026-11-07.md`），冒烟通过；11-07 b05 整跑完成（候选 21 / 18，sB 首次找到新 E2，同键图的脊约 190k forward / 块无产出），见 DONE_LOG |
| 9e | 粗阶段反力差门（冻结 σ 0.0025 → 0.0106；粗阶段实测 σ_c 更大，v4 验收中挡下 2 次） | 只报告；用户决定保留门 |
| 机制 | 软模与 fp32 噪声的解释 | 用户 11-04："可以随后补" |
| S3 + 并行（11-07） | 同键图的脊按 `audit_fraction` 抽查 + HS / HC 按波次并行（D2-HF-v1.4） | 代码与测试完成，测量冒烟已做（`docs/PLAN_HF_D2_V14_2026-11-07.md`）|
| R3（11-07） | `hysteresis` / `walk_relax_unconverged` 占 HS 约 95% forward；选项文档 `docs/HF_HS_NO_MAXIMUM_OPTIONS_2026-11-07.md` | R3a-2 已施工为 D2-HF-v1.2（`docs/PLAN_HF_SAME_GRAPH_RIDGE_2026-11-07.md`），离线重放通过、整跑未做；R3b、R3c 未定 |
| TODO 11、12、10 端点、3、4、5、7 | 见 `docs/TODO.md` | 未动；12 待用户定改法 |

## 6. 提交草稿（用户提交，按顺序；都在 `.commit_drafts/`）

上一份 b05 交接所列各草稿 → `2026-11-01_1400_b05-qm-and-saddle-checks` → `2026-11-02_2000_hf-midpoint-coarse` → `2026-11-04_1944_hf-midpoint-fine-fire` → `2026-11-06_1752_hf-midpoint-identity-hc`。

几份草稿的 `git add` 有重叠文件（`ts_ef.py`、`test_ts_ef.py`、DONE_LOG、TODO、交接文档、README 等）。`git add` 加进去的是工作区当前的内容，所以先提交的草稿会带上后面的改动；要让历史分清，就合成一次提交。各草稿的 Note 里都写了这一点。

## 7. 环境与操作要点（都已实测）

- **结果在写入那台读**：本线远程结果在 183 写入，用 `ssh korakuen@192.168.0.183` 在 183 上读；本机读要等文件写完 6 h 以后（NFS 属性缓存）。
- **183 列目录会重复**：11-02 起 183 上 `ls` / `glob` 结果目录时，每个文件名会重复出现（多为 8 次），文件只有一份；在 183 上汇总一律按文件名去重（各 `*_summary.py` 已这样做）。
- **远程后台启动**：`ssh … "bash -c 'cd X && { setsid nohup CMD > log 2>&1 < /dev/null & echo \$!; }'"`。不加花括号，ssh 不返回，`$!` 拿到的是子 shell 的 PID。
- **183**：2080 Ti，40 核。别的项目（THERMO 等）有时占 CPU 或 GPU，启动前先看 `uptime` 和 `nvidia-smi`。用户要求"慢慢跑"时只开 1 进程、8 线程。v5 验收时实测 0.12–0.19 s/forward。
- **PBS**：yayoi35 / yayoi41 是 RTX 5080，`openmm_dev` 的 torch 2.12.0 / CUDA 12.9 含 sm_120，已在 yayoi35 实跑通过。脚本照模板 `/home/kasuga/.claude/templates/pbs.zsh`，范例 `benchmarks/explore_d2/midpoint_transverse_hessian_pbs.zsh`。只能交冒烟跑通过的代码。
- **变异检验**：变异代码放临时目录，`PYTHONPATH=<临时目录>:src` 配合 `-o pythonpath=` 才真的被测。
- **脚本索引**（`benchmarks/explore_d2/`）：
  - `hs_transverse_diag.py`（按阶段诊断）
  - `midpoint_cold_replay.py`（冷启动对照；也提供 `load` / `context` / `cold` / `rerun`）
  - `midpoint_transverse_hessian.py`（横向曲率）
  - `midpoint_threshold_replay.py`（离线放宽重放）
  - `midpoint_sign_resolution.py`（界线测量）
  - `midpoint_coarse_validate.py`（v3 / v4 验收；`settings_for` 构造 v3 / v4 / v5）
  - `midpoint_optimizer_compare.py`（O5）
  - `midpoint_fire_validate.py`（v4 / v5 验收）
  - `candidate_identity.py`（候选身份）
  - `midpoint_fire_hc.py`（v5 HS + HC）
  - 各自配 `_run.sh` 或 `_summary.py`。

## 8. 当前没有运行中的作业

2026-11-06 18:05 核对：183 与 181 上没有本线启动的进程。

## 未确认（Claude 推断）

- （11-06 全量诊断后：阈值在噪声内得到支持，但卡住帧的残余力只有 4–35% 在最软模上，"同一件事"只部分成立。）软模与 fp32 力噪声可能是同一件事：λ ≈ 0.019 eV/Å² 时，要让该模的力 ≤ 1e-3，位置须落在模极小 ±约 0.05 Å 内，此处回复力与测得的力抖动同量级。逐原子分配与噪声投影没算。
- v4 / v5 之间、以及 11-04 与 11-06 两次 v5 之间的个别结局差别，我归为 fp32 波动，依据是这些差别出现在两版代码相同的阶段，或者同一版本重复运行时结局也不同。没有逐项核查。
- FIRE 只用了默认参数，b05 以外的体系没有测过。
- HC `endpoint_unstable` 两项的原因没查。可能和这个鞍点的软模有关，也可能是 `branch_starts` 选的起点离鞍点太近。
- 11-06 HC 重跑得到的候选没有再用对称 RMSD 比几何，身份是按 HC 化学对和能量判的。
