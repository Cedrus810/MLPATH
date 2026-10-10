# b09 PRRS 冒烟探针（2026-10-08，不是验收）

- `input.extxyz`：由 `../build_b09_input.py` 从 POLAR-1 上的 TS-2OMe 沿虚频退到反应物一侧、再做 BFGS 极小得到的复合物
  （卡宾 C 到 Cb 3.12 Å、到 Ca 3.61 Å、到醛氧 4.12 Å；无虚频；比 TS-2OMe 低 0.477 eV）。
  **它就是通向 TS-2OMe 的那个复合物，对其他通道有偏**；正式探针要另想无偏起点。
- `config_smoke.json`：b04 冻结配置去掉注释后改 `total_charge=0`、`multiplicity=1`、`max_trials=4`，
  `tracked_coordinates` 换成卡宾 C（2）到 Cb（15）、Ca（13）、醛氧（1）。其他值沿用 b04，**没有为 b09 调过**。
- 命令（183）：`PYTHONPATH=.:src python -m prrs search benchmarks/carbene_ambimodal/probe_b09/input.extxyz
  --output runs/b09_smoke0 --config benchmarks/carbene_ambimodal/probe_b09/config_smoke.json
  --calculator benchmarks.b04_sn2_chloride.calculator:polar1`，外套 `timeout 1200`。
- 目的：能不能跑通、单个 trial 多久、碰到了哪些鞍点。

## 冒烟结果：`runs/b09_smoke0`（183，PID 949899，20:45–20:48）

- 跑通了。起点在 POLAR-1 上弛豫确认为极小（25 s），作为 c0000 收进网络（10 s），4 个 trial 共 2 分 16 秒（**约 34 s/trial**）。
- **4 个 trial 全花在同一个方向上**：环己烯环上 C14–C16（CH₂–CH₂）拉伸的幅度阶梯 0.1 / 0.3 / 0.6 / 1.0 Å。
  0.1 Å 是 elastic（回到原处）；其余 3 个被 `collateral_bond` 门拒绝：相邻的 C16–C20 变化 0.052 / 0.091 / 0.119 Å，
  超过 b04 配置里的 `collateral_bond_A = 0.05`。
- 没有任何 trial 碰到卡宾的接近方向，0 反应、0 鞍点。**这次冒烟回答不了"能不能找到鞍点"**，只说明能跑、多快。
- 两件要注意的事：(1) 一个方向的幅度阶梯就是 4 个 trial，`max_trials=4` 只够扫一个方向；
  (2) `collateral_bond_A = 0.05` 是 b04（6 个原子）上定的门，在 29 原子的柔性环上拉伸会连带邻键，可能大量误拒。
  按规矩不自己改门，要和用户商量。

## 探针 1：`runs/b09_probe1`（183，PID 950185，10-08 21:08 → 10-09 00:56；log `../results/2026-10-08/b09_probe1.log`）

`config_probe1.json` = 冒烟配置把 `max_trials` 改成 96。**我的配置错误**：`conformer_trials_per_microstate` 仍是 32，
只有一个微观态时 96 根本到不了，run 在 32 个 trial 处停（`termination_detail.directions_dropped_by.conformer_trials_per_microstate = 17`）。
实测约 7.7 min/trial（会爬鞍点、做中继下降），不是冒烟的 34 s。

- **0 个反应，0 个真鞍点**。32 个 trial：16 elastic、10 `ts_candidate`、5 被 `collateral_bond` 拒（C16–C20 三次、C14–H27 两次，门 0.05 Å）、1 淬火不收敛。
- 10 个 `ts_candidate` 全是复合物上的软模：虚频 −18…−51 cm⁻¹，比起点高 1–6 meV，到起点的重原子 RMSD 0.04–0.08 Å，
  卡宾 C 到 Ca / Cb / 醛氧 3.57 / 3.08 / 4.13 Å（真 TS-2OMe 的 C:–Cb 是 1.98 Å）。
- 9/32 个 trial 花在细化上（compress 5-16、kick 5-16、compress 4-15 各 3 次），触发细化的"离开"就是这些软模假候选。
- **反应坐标从没被选中**：432 个候选方向里选 24、执行 7。`compress:2-15` / `kick:2-15`（卡宾 C → Cb）、`2-13`（→ Ca）
  都在 `uncovered_ids`（`max_directions` 截断）。实际推的是甲基 C5 → 环 C16、甲氧基 O4 → Cb、环内 C14–C16 / C14–H17 拉伸和一个弯曲。
  账本自己写明了这不能当阴性结果读。

## 方向排序诊断（10-09，`../rank_b09.py`，log `../results/2026-10-09/rank_b09.log`、`rank_b04.log`）

在 run 自己的源微观态和配置上重算排序。旧打分（端点梯度 g/√m）重现了 run 实际选中的 24 个方向（24/24），诊断可信。

- **4 个软模的非扭转分量保留 0.994–0.999**：b09 上不存在"几乎全是扭转的模被放大占名额"的情况。
- **打分修复（`src/prrs/perturbations.py` `delivered_direction`）**：位移类（compress / stretch / bend / torsion）改为用
  `apply` 实际执行的刚体运动打分，kick / pulse 不变（它们本来就沿 g/√m）。
  - b09：`compress:2-15` 分数 0.186 → 0.515，**名次 21/84 → 36/84，仍选不上**；`2-13` 39 → 47/84。
  - b04：C–Cl 拉伸仍是 stretch 族第一（0.965 → 0.990，第二名 0.184），没有退步。
- **结论**：分数现在对了，但对 b09 没用。原因在先验本身：同一对分子之间的每个 compress，执行的都是两个整分子沿原子对轴线的刚体平移，
  分数只反映"这根轴和最软的分子间运动对得多齐"，与哪两个原子会成键无关（b04 的四个 Cl⁻ 接近方向全在 0.994–1.000）。

> 同日 Codex 的独立核查和预算对照见 [FOLLOWUP_2026-10-09.md](FOLLOWUP_2026-10-09.md)，之后的源码修复（受限软模、方案 B、环内保键交付）
> 见 [docs/SEARCH_DESIGN_REPAIR_2026-10-09.md](../../../docs/SEARCH_DESIGN_REPAIR_2026-10-09.md)。上面的排序诊断是在这些修复之前做的。

## 全候选粗扫冒烟：`runs/b09_coarse_smoke`（10-09，183，PID 953421，`timeout 900`；log `../results/2026-10-09/b09_coarse_smoke.log`）

配置就是正式的 `config_coarse_coverage.json`（exhaustive，三层上限 1408），源码为 Codex 修复后的版本。15 分钟时被 timeout 停掉。

- 启动 15:40:16 → 源微态收进网络 15:41:34（约 1.3 min）。
- 完成 5 个 trial，第 6 个跑到 1 min 42 s 时被停：

| trial | 方向 | 振幅 | 结果 | 耗时 |
| --- | --- | --- | --- | --- |
| t000000 | stretch 2-4（卡宾 C–甲氧基 O） | 0.1 Å | 返回 | 2 min 21 s |
| t000001 | 同上 | 0.3 Å | 返回 | 1 min 53 s |
| t000002 | 同上 | 0.6 Å | 返回 | 2 min 03 s |
| t000003 | 同上 | 1.0 Å | ts_candidate ts0000 | 4 min 22 s |
| t000004 | **compress 2-13（卡宾 C → Ca）** | 0.1 Å | 返回 | 1 min 36 s |

- 平均约 **2.4 min/trial**（5 个的均值 147 s）。冒烟 0 的"34 s/trial"是被门拒掉的 trial 拉低的（拒掉的不跑动力学），不能用来估算。
- 新排序下 exhaustive 的第二个方向就是卡宾 C → Ca。
- **中途杀掉等于白跑**：`resume_point.json` 停在 15:41:38、`attempts: 0`，5 个已完成 trial 不在检查点里。续跑以微态为边界，
  源微态的全部粗扫是一个事务，杀掉就整体回滚。
- 日志里仍有一次 `runner.py:2279` 的 Logging error（中继下降那行），不影响计算。

## 方案 B：只扫 compress 族的全候选粗扫（10-09 15:57 起，183，PID 953745；log `../results/2026-10-09/b09_compress_coverage.log`）

- 用户 10-09 选"先跑 B"。配置 `config_compress_coverage.json` = `config_coarse_coverage.json` 只改 `families: ["compress"]`、三层上限 336。
  exhaustive 预检通过（84 个方向 × 4 个振幅 = 336）。
- 与 Codex 协议的差别：只覆盖分子间 / 非键接近一族，不含 kick（672）、stretch（112）、bend（288）。
- 预计约 13.5 h（冒烟实测 2.4 min/trial × 336），范围 9–25 h。**不可续跑**：中途停掉，源微态整批回滚。
- 输出 `runs/b09_compress_coverage`。

### 结果（10-10 在 183 上读；`../analyze_coverage_b09.py`，log/json `../results/2026-10-10/b09_compress_coverage_analysis.*`）

- 10-10 06:50:48 正常结束，`status: completed`，336/336 trial，**实测 14 h 53 min ≈ 2.66 min/trial**。
- **覆盖完整**：84 个方向全部执行，`not_executed 0`、`uncovered 0`、`coverage_censored false`、`directions_dropped_by {}`。
  预算正好等于粗扫量，所以没有细化 trial。
- **0 个反应**。273 completed、26 淬火不收敛、37 `ts_candidate`；**0 个被门拒**（所有 trial 都跑了动力学）。
  13 个 `rearrangement` 都没有键变化（formed/broken 为空），是构象变化；新收了 2 个微态，其中 `m0002` 比源微态低 214 meV
  （源微态不是最低构象）。
- **37 个 `ts_candidate` 全是复合物上的软模**：虚频 −17…−118 cm⁻¹，比源微态高 0.4–23 meV，到源微态的重原子 RMSD 0.03–0.23 Å，
  卡宾 C 到 Ca / Cb / 醛氧 ≥ 3.43 / 3.05 / 3.99 Å。"最近的 DFT 鞍点"一栏全是 TS-2OMe（0.32–0.40 Å），这只是因为反应物复合物
  本来就和 TS-2OMe 最像，没有意义。
- **反应坐标这次被推了，但推得不够远**。整条轨迹上卡宾 C 离得最近的距离：

| 方向 | 振幅 1.0 Å 时轨迹最近距离 | POLAR-1 鞍点上的值（上文表） | 结果 |
| --- | --- | --- | --- |
| compress 2-15（→ Cb） | C:–Cb 2.12 Å | TS-2OMe 1.98 Å | 弹回（构象变化） |
| compress 2-13（→ Ca） | C:–Ca 2.60 Å | TS-2OMe-C2 的 DFT 值 1.81 Å | 淬火不收敛 |
| 全部 336 个 trial 中 C:–O_ald 最近 | 3.18 Å（compress 1-4） | TS-2OMe-2 的 DFT 值 2.84 Å | — |

  **C:–O_ald 方向本身不在候选集里**（10-10 外部复核指出，已核实）：源微态上卡宾 C 到醛氧 4.123 Å，超过 `pair_cutoff_A = 4.0`，
  所以 84 个 compress 方向里没有 2-1。单凭这一点不足以解释为什么没出现 (4+1)。

  最近点就是刚推完的那一帧（源 3.124 − 1.0 = 2.124），之后体系往回走。粗扫振幅表 `geometry_amplitudes_A` 最大是 1.0 Å，
  而源微态到 POLAR-1 TS-2OMe 的 C:–Cb 要缩短 1.15 Å。细化只在已有振幅之间二分，不会超过 1.0 Å。
  所以这次的阴性结论只适用于"振幅 ≤ 1.0 Å 的刚体压缩"，不能说明 PRRS 在 b09 上找不到第二条路径。
  但也不能反过来说"再推近一点就能过"：下面的冒烟里推 1.2 Å 已经把 C:–Cb 推到 1.924 Å（比 TS 的 1.979 Å 还短），仍然弹回；
  接近姿态和其他坐标同样起作用。
- 没跑的族：kick（672）、stretch（112）、bend（288）。kick 是带动能沿 g/√m 推，不受几何振幅上限限制。

## 加大力度冒烟：`../push_b09.py`（10-10 15:25–15:40，183，PID 956356；结果 `../results/2026-10-10/push_b09/`、log `push_b09.log`）

从覆盖扫描的源微态 `c0000_m0000` 出发，用同一配置，每个探针跑一个 `run_trial`（search 用的同一个调用），终点用 `confirm_minimum` 判断。

| 探针 | 推完那一帧比源高 | 轨迹上最近 C:–Ca / C:–Cb | 结果 |
| --- | --- | --- | --- |
| compress 2-15，1.2 Å | +2.80 eV | 2.53 / 1.92 Å | 弹回，终点是复合物上的软模鞍点（+22 meV） |
| **compress 2-15，1.5 Å** | **+6.16 eV** | 1.47 / 1.44 Å | **成环丙烷（(2+1) 产物）**：推完即成 C2–C15，自由段 step 131 成 C2–C13；末帧 C:–Ca 1.52、C:–Cb 1.48、Ca–Cb 1.52、C:–O_ald 3.11 Å，比源低 1.17 eV。**但淬火 500 步没到 fmax 0.01，记为 `quench_failed`，search 里会被丢掉，不记反应** |
| kick 2-15，0.5 / 1.0 eV | — | ≥ 2.57 / 2.36 Å | 弹回（kick 只给 C2、C15 两个原子动能，能量很快进了分子内振动） |
| compress 2-13，1.2 Å | — | 2.40 / 2.05 Å | 淬火不收敛，无键变化 |
| compress 2-13，1.5 Å | — | — | 被 `collateral_fragmentation` 门拒（推完分子数 2 → 1，即非目标原子对已成键） |
| kick 2-13，0.5 / 1.0 eV | — | ≥ 2.91 / 2.57 Å | 弹回，终点是源附近的极小 |

- 结论：**刚体压缩推到 1.5 Å 能撞出 (2+1) 环丙烷**，但这是 +6 eV 的硬撞，说明不了走的是哪个鞍点；(4+1) 产物（C–O_ald 成键）没出现。
  能量参照（10-10 在 183 上用 POLAR-1 重算，`results/2026-10-08/opt/*/geometry.extxyz` 相对本次源微态 `c0000_m0000`）：
  TS-2OMe +0.477、TS-2OMe-cis +0.642、TS-2OMe-C2 +1.171、(2+1) 产物 −1.181、(4+1) 产物 −2.012 eV。
  （此前这里写的"TS-2OMe 能垒约 +3.8 kcal/mol ≈ 0.16 eV"是错的：+3.8 是 TS-2OMe-cis 相对 TS-2OMe 的能差。）
- 淬火问题：这个 trial 的产物明明到了，却因 500 步淬火没收敛被判失败。
- 覆盖扫描里 26 个 `quench_failed`（25 个首次 FIRE 未到阈值，1 个 polish 后重新优化失败；最后一帧 fmax 0.01049–0.22014 eV/Å）：
  `raw_endpoint` 上卡宾 C 到 Ca / Cb / 醛氧都 ≥ 2.94 Å，没有卡宾–环新键。**但这不等于没有漏收化学变化**：`t000327`
  （compress 4-6，1.0 Å）推完先成 O3–H6，free 第 9 步断开，第 13 步 C5–H6 断裂，末帧 C5–H6 10.22 Å——甲氧基上一个 H 被撞飞。
  它没过终点验收，也没被收录。（此前写的"没有漏掉反应"说得太满，10-10 更正。）

### p001 为什么淬火失败（10-10，`../results/2026-10-10/quench_p001/continue_quench.py`、`.log`）

- 淬火用 FIRE，上限 `quench_steps = 500`，判据 fmax 0.01 eV/Å。p001 的淬火轨迹（每 5 步采样）：能量在采样点上持续下降，
  fmax 从 4.8 降到 0.02 附近，最低 0.0169（第 495 步），第 500 步是 0.1091；采样中每隔约 50 步有一次 fmax 尖峰，
  是不是 FIRE 内部重置，单凭这些记录确认不了。
- 从同一个淬火起点（`raw_endpoint.extxyz`，即自由段末帧，fmax 4.85，最大力在 O3、C15、C2、C13）接着优化：
  - FIRE **834 步**收敛；
  - BFGS **139 步**收敛。
  
  两者落在同一个环丙烷（C:–Ca 1.52、C:–Cb 1.48、Ca–Cb 1.52、C:–O_ald 3.12–3.13 Å），能量 −16770.27 eV，比源微态低 1.18 eV。
- 结论：普通力优化下，这次失败是 FIRE 步数上限不够。注意 834 / 139 步只是单纯的力优化，没走 runner 的软模 polish、
  `confirm_minimum` 和分类，"BFGS 快六倍"也只是步数比；完整验收见下一节重放。
- 10-10 更正：最初导出的 `fire_end.extxyz` / `bfgs_end.extxyz` 带的是起点的能量和力（共享 calculator，打印能差时又算了一次起点）。
  脚本已改为最后对终点求值再写出，重跑后文件里的能量 −16770.2700 / −16770.2737 eV，fmax 0.0100 / 0.0096，与日志一致。

### 淬火重放冒烟：p001（10-10 16:42–16:45，183，PID 957756；`../replay_quench_b09.py`，结果 `../results/2026-10-10/replay_quench_smoke/`）

从 p001 的 `raw_endpoint` 出发，走 runner 的完整链：`quench`（FIRE + 软模 polish）→ `state.classify` → `confirm_minimum`。
只把 `quench_steps` 从 500 改为 1500，其余是覆盖扫描的配置。

- FIRE 834 步收敛，polish 1 轮即通过，fmax 0.00998；**分类 `reactive`，成键 C2–C13、C2–C15，`confirm_minimum` 通过（0 阶）**。
- 比源微态低 1.177 eV；与 POLAR-1 优化的 (2+1) 产物重原子 RMSD **0.029 Å**（与 (4+1) 产物 1.08 Å）。
- 耗时 134 s（淬火 118 s，即约 0.14 s/FIRE 步）。
- 结论：p001 在完整验收链下是一个被 500 步上限漏收的 (2+1) 反应。
