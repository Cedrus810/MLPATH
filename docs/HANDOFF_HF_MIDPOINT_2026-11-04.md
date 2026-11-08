# 交接：二分中点卡住（TODO 第 9 条）工作线（2026-11-01 至 11-04）

写给接手这条线的窗口。本文件只汇总现状和去处，数字以链接的原文为准。上一份交接（b05 HF D2，10-26 至 11-01）见 [HANDOFF_HF_D2_2026-11-01.md](HANDOFF_HF_D2_2026-11-01.md)。

## 必读（开工前全部读完）

- `docs/TS_HESSIAN_FREE_DESIGN_RECORD_2026-10-24.md`：设计稿；§5.1 反力零点二分、§5.3 "以实测分辨率决定终止"是这条线的依据。
- `docs/HF_MIDPOINT_THRESHOLD_OPTIONS_2026-11-02.md`：问题、代码里中点阈值管什么、离线重放、选项 O1–O5。
- `docs/PLAN_HF_MIDPOINT_COARSE_2026-11-02.md`：用户 11-02 对 O3 的原文决定、界线 B 的测量与决定、同机验收结果。
- `docs/PLAN_HF_MIDPOINT_O5_2026-11-02.md`：优化器对照的设计与结果。
- `docs/DONE_LOG.md` 2026-11-01、11-02 各节：每一步的命令、数字、更正。
- `docs/TODO.md` 第 9 条及"未确认"一节：进展与仍未验证的推断。
- 全局与项目 `CLAUDE.md`：NFS 缓存（跨主机结果在写入那台读）、PBS 规矩、只按 PID 杀进程。

## 1. 问题

二分中点的约束松弛要求横向力 ≤ 1e-3 eV/Å（`hf.midpoint_transverse_fmax_eV_A`），200 步内达不到，HS 以 `unresolved:transverse_not_converged` 结束。b05 八块 HF 运行里 268 次中点松弛失败 71 次，来自 22 个意图。新 E2 鞍点（QM 确认，比已知 E2 TS 低 11.75 kcal/mol）就藏在这类项里。

## 2. 已查清的事（按时间）

| 步骤 | 结论 | 记录 |
| --- | --- | --- |
| 按阶段诊断 | 不是噪声底的平台：失败段末 50 步横向力约为阈值 10 倍，步长贴近 0.05 Å 上限，同一 q 上反力极差中位 0.065 | DONE_LOG 11-01 |
| 同帧冷启动对照 | 不是带入曲率的问题：失败帧冷启动 24 / 213 收敛，正常帧 70 / 90 | DONE_LOG 11-01；`PLAN_HF_MIDPOINT_COLD_2026-11-01.md` |
| 横向曲率（PBS yayoi35） | 失败帧无横向负曲率，有一个单独的软横向模（0.013–0.024 eV/Å²，次低约 0.12），游走中位 96% 沿此模；同样软的中点也有收敛的 | DONE_LOG 11-02；`PLAN_HF_MIDPOINT_TRANSVERSE_2026-11-02.md` |
| 离线重放 | 71 个失败中点中 65 个的括区方向在 2e-3 下已定（反力同号、|RF| 约 0.1） | `HF_MIDPOINT_THRESHOLD_OPTIONS_2026-11-02.md` §3 |
| 界线测量 | σ_c 中位 0.0040、p90 0.0073；粗、精一致偏差未见；**用户定 B = 0.025** | `PLAN_HF_MIDPOINT_COARSE_2026-11-02.md` §2 |
| v4 同机验收 | 对照 15 项候选不变；失败 22 项救回 2 项；成本 +1.8%；剩余卡点在近根中点的精松弛 | 同上 §5 |
| O5 优化器对照 | 同帧冷启动：FIRE 近根卡点 33 / 38（BFGS 3 / 38），远根失败帧 54 / 71（3 / 71），对照 28 / 30（23 / 30），forward 更少；收敛反力与 BFGS 差 ≤ 0.004；缩小 BFGS 步长无效 | `PLAN_HF_MIDPOINT_O5_2026-11-02.md` §6 |

## 3. 代码状态

| 部分 | 现状 |
| --- | --- |
| `explore/src/explore/ts_ef.py` | v5（11-04）：`EFSettings.midpoint_fine_optimizer="fire"` 时精阶段用 FIRE（速度归零、接受后清 BFGS 曲率、失败不改 state），粗阶段及其余仍为 BFGS；以下为 v4。v4 分级终止：`EFSettings.midpoint_coarse_transverse_fmax_eV_A` + `midpoint_sign_bound_eV_A`，两者都给才启用，provenance `HF-EF-bracket-v4`；不给则与 v3 完全相同。**默认关闭，未接入 D2**（`HFSettings`、manifest、`_ef_settings` 都没改） |
| FIRE | 只在 `benchmarks/explore_d2/midpoint_optimizer_compare.py` 里临时替换 `ts_ef.BFGS` 测过，`explore/` 里没有 |
| 测试 | `test_ts_ef.py` 新增 9 项；`test_ts_ef test_ts_branch test_ts_connect test_ts_warm_start test_d2_hf test_d2 test_d2_resume` 103 passed；两个变异检验被测出；全包 `make check` 没跑 |
| 脚本 | `benchmarks/explore_d2/` 下：`hs_transverse_diag.py`（按阶段）、`midpoint_cold_replay.py`、`midpoint_transverse_hessian.py`、`midpoint_threshold_replay.py`、`midpoint_sign_resolution.py`、`midpoint_coarse_validate.py`、`midpoint_optimizer_compare.py`，各带 `_run.sh` 或 `_summary.py` |

## 4. 待办（按建议顺序；方向需用户定）

| 编号 | 内容 | 状态 |
| --- | --- | --- |
| 9a | FIRE 接进中点松弛后做同机验收：救回数、候选一致性、成本 | **已完成（11-04，用户定只换精阶段）**：HF-EF-bracket-v5；失败 22 项中 14 项出候选并稳定复现，对照 15 项不变；见 [PLAN_HF_MIDPOINT_FINE_FIRE_2026-11-04.md](PLAN_HF_MIDPOINT_FINE_FIRE_2026-11-04.md) |
| 9b | v4 救回的 2 个候选与 v5 的 14 个候选（能量都在新 E2 鞍点一段）：先比几何身份，再做 HC；新增候选与新增事件分别计数（用户 11-04） | **几何身份已比（11-06）**：全部是新 E2 鞍点，新增候选 15 个意图、新增事件 0（DONE_LOG 11-06）；**HC 已做（11-06）**：v5 重跑 15 个出候选，13 个 HC 连通、化学对全部与新 E2 事件相同，2 个 `endpoint_unstable`（高能侧没离开鞍点，原因未查） |
| 9c | v5 下仍卡在中点精阶段的 4 项（`k0·8-14`、`k1·8-14`、`k2·6-8`、`k2·8-14`），以及 O5 里 FIRE 未收敛的帧 | 未查；用户 11-04：仍卡在精阶段就集中诊断剩余帧 |
| 9d | v4（及以后的 FIRE 改法）接入 D2：新规则版本、`HFSettings` 字段、b05 整跑 | 未开始 |
| 9e | 粗阶段反力差门：冻结 σ 0.0025 → 0.0106，粗阶段实测 σ_c 更大，验收中挡下 2 次 | 只报告；用户决定保留门 |
| 11、12、10 端点、3、4、5、7 | 见 `docs/TODO.md` | 未动；12 待用户定改法 |

## 5. 提交草稿（用户提交，按顺序；都在 `.commit_drafts/`）

上一份交接所列各草稿 → `2026-11-01_1400_b05-qm-and-saddle-checks` → `2026-11-02_2000_hf-midpoint-coarse`（本线全部：v4 代码与测试、各测量与验收、O5 脚本与结果、文档）。

## 6. 环境与操作要点（都已实测）

- **结果在 183 上读**：这条线的远程结果都在 183 写入，必须 `ssh korakuen@192.168.0.183` 在 183 上读（NFS 属性缓存 6 h）。
- **183 列目录会重复**：11-02 起 183 上 `ls` / `glob` 这些结果目录时每个文件名重复出现（多为 8 次），文件只有一份。在 183 上汇总一律按文件名去重（各 `_summary.py` 已这样做）。
- **远程后台启动的写法**：`ssh … "bash -c 'cd X && { setsid nohup CMD > log 2>&1 < /dev/null & echo \$!; }'"`。不加花括号时 `&` 作用于整个 `cd && …`，派生的子 shell 拿着 ssh 的输出，ssh 不返回，`$!` 也是子 shell 的 PID（11-01 踩过）。
- **183 负载**：别的项目（THERMO）会占 10 核左右；GPU 常空。用户 11-02 要求"慢慢跑"时只开 1 进程 8 线程。183 在 11-03 下午重启过（11-04 uptime 19 h）。
- **PBS**：GPU 节点 yayoi35 / yayoi41 为 RTX 5080；`openmm_dev` 的 torch 2.12.0 / CUDA 12.9 含 sm_120，已在 yayoi35 实跑通过。PBS 脚本照模板 `/home/kasuga/.claude/templates/pbs.zsh`，范例 `benchmarks/explore_d2/midpoint_transverse_hessian_pbs.zsh`。
- **变异检验**：变异代码放临时目录，`PYTHONPATH=<临时目录>:src` 配合 `-o pythonpath=` 才真的被测（pyproject 的 pythonpath 会盖掉 PYTHONPATH）。

## 7. 当前没有运行中的作业

2026-11-04 13:09 核对：183 与 181 上没有本线启动的进程。

## 未确认（Claude 推断）

- 软模 + fp32 力噪声可能是同一件事：λ ≈ 0.019 eV/Å² 时，要让该模的力 ≤ 1e-3，位置须落在模极小 ±约 0.05 Å 内，回复力与测得的力抖动同量级。逐原子分配与噪声投影没算。
- FIRE 的优势来自冷启动单次松弛；实际运行带入曲率、两阶段交接、以及 FIRE 默认参数在别的体系上的表现都没测。
- v4 救回的 2 个候选能量与新 E2 鞍点同段，只是能量接近，身份没查。
- σ_c 用冷启动加 0.01 Å 扰动测得，与实际运行带入曲率时的抖动可能不同。
