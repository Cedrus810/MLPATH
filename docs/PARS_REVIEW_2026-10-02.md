# PARS 审查、修复与 ground truth 实验 · 2026-10-02

对新模块 `pars/`（Pareto-Adaptive Reaction Path Search，设计稿见根目录
`Pareto-Adaptive Reaction Path Search (PARS) — MLPATH Integration Design.md`）的一次完整过程记录：
通读审查 → 逐条修复（每条先写会失败的测试）→ 在真实 run 上用 MACE 全量精修做 ground truth。

**状态**：pars 测试 116 个全部通过，`make check`（fmt-check + lint + test）通过。**未提交**
（`pars/` 整体仍是 untracked）。只动了 `pars/` 和 `docs/experiments/` 下的两个新脚本，**没有动
`src/prrs/`**，PRRS 指纹不受影响。

**结论先说**：

1. 精修后端原来大多是假成功。在 p1 两个 run 共 27 个候选上，原代码只有 4/27 找到真实鞍点，
   22/27 报 `NEB_CONVERGED`、势垒约 0.01 eV（真实值 0.419 eV）。三处修复后为 27/27，单个候选耗时
   中位数从 2601 s 降到 311 s。
2. p1 数据检验不了 PARS 的选择：27 条路径全落在同一个鞍点上，任何选择规则（包括随机抽取）都是
   满分。要检验 PARS 的核心假设（压缩精修成本又不丢机理多样性），需要多鞍点 / 多通道的数据。
3. 以前用旧代码跑出的 PARS selection 和 refine 结果都不可信，需要重跑。

---

## 1. 审查中确认并已修复的 bug

每条都有数值复现，修复前先写测试并看着它失败。

| # | 位置 | 问题 | 怎么确认的 | 修复 |
| --- | --- | --- | --- | --- |
| B1 | `frontier/objectives.py` Kabsch | 旋转矩阵用成了转置 | 随机旋转后对齐，RMSD 仍有 1.84 Å；纯刚体转动的曲率算成 0.17。原测试每帧用同一个旋转，测不出来 | `p @ rotation`；新增"每帧旋转不同"的测试 |
| B2 | `refinement/ase_neb.py` CI 接受判据 | 要求 CI 势垒 ≤ 普通 NEB 势垒，但 CI 本来就会把势垒抬到鞍点上 | 解析势上真实势垒 2.116 eV：CI 给出 2.12，普通 NEB 给出 2.10，CI 被拒 | 改成几何判据 `_climbing_stage_accepted`：收敛、顶点在内部、CI 顶点相对普通 NEB 的位移不超过一个像点间距 |
| B3 | `tests/conftest.py` `TripleWellCalculator` | 高斯项的力符号写反 | 有限差分得到 −dE/dr = (−3.33, −2.67)，解析力给 (+2.93, +2.37)。精修测试一直在能量和力对不上的势上跑，断言写的是"两种状态都算对"，所以 B2 没暴露 | 改符号；新增力等于 −梯度 的测试 |
| B4 | `frontier/dominance.py` / `margin.py` | ε-支配不反对称：`dominance_tolerance > 0` 时候选可以互相支配，前沿变空，最后零代表；方向也和 RATIONALE 写的"放宽前沿"相反 | (0, 0.05) 与 (0.05, 0) 在 tol = 0.1 时前沿为空 | 改为"每个目标都不差、且至少一个好出 tol 以上"（严格偏序，tol 越大前沿只会越宽）；硬删除要求每个目标上都隔开 tol |
| B5 | `frontier/objectives.py` 曲率 | 注释写 ∫\|d²R/ds²\|² ds，代码除以 h 而不是 h³；施力段之后的跳变也被算进去 | 路径整体放大 2 倍，旧值翻倍（正确应减半） | 除以 h³；只在施力间隙之后的连续自由段上算 |
| B6 | 曲率（B5 修完后暴露） | 在原始 MD 帧上算，实际测的是"自由段有多长 × 振动多少"，和 recrossing 重复 | ρ(曲率, 帧数) = 0.92 / 0.86，ρ(曲率, recrossing) = 0.80 / 0.89 | 在 NEB 像点分辨率（`refinement_images` 个等弧长点）上算。修后 ρ(帧数) = 0.29 / 0.23（n = 15 时 p ≈ 0.29，不显著），ρ(recrossing) = 0.06 |

## 2. 其它修改（设计层面，不改变默认行为的标了"默认不变"）

- **不确定度接线**（默认不变）：archive 把候选上的 `objective_uncertainty`（原始单位）换算到归一化
  单位后传给分层。V1 没有任何目标带实测不确定度，所以默认配置下不会出现 `CERTAINLY_DOMINATED`，
  README 里写明了。
- **`preserve_mechanism_diversity`** 原来没有被任何代码读取；现在为 False 时对整个活跃集只做一次覆盖。
- **ingest 统计**（过滤规则不变）：PRRS 源码确认 `ts_candidate` 是 quench 停在鞍点上，`unconfirmed`
  是新拓扑没能在自由段末尾保持住，两者都不是确认的 A→B，继续跳过。selection 里新增
  `counts.skipped_attempts`（按 status/classification 计数）和 `counts.channel_members_not_ingested`；
  没有候选时，报错会带上这份统计。
- **`--channels`** 改为在构建候选之前筛选，不在所选通道里的坏 trace 不会再让整个 select 失败。
- **精修端点检查**：端点弛豫没收敛、或两个端点落进同一个极小（RMSD < `refinement_endpoint_min_rmsd`
  = 0.05 Å），都报 `FAILED` 并给出原因。
- **refine 续跑**：每个候选跑完立刻写 `refined-<id>.json`，里面记录 conditions（种子 sha256、factory
  字符串、全部 `refinement_*` 设置、PARS 版本）。重跑时，条件一致的候选跳过；任何已有结果的条件
  不一致或没有记录条件，整次运行在开算之前就拒绝，并列出是哪个候选、哪一项不同；`--redo` 才覆盖。
  已知漏洞：factory 字符串不变但内部换了模型，检测不到。
- **provenance**：精修结果保留原始 run 的 provenance（种子文件路径另记在 `extra.seed_path`）；
  一个 attempt 被多个 reaction 列出时直接报错（真实 run 里没有这种情况）。
- 删掉了 `cover.py` 里的死分支；`kabsch_align` / `align_frames` 移到 `paths/resample.py`。

## 3. Ground truth 实验：全量精修 p1 的 27 个候选

**做法**：`p1_round7`（15 个）和 `p1_committor1`（12 个）的全部 reactive 候选，用 run 自己的模型
（MACE-OFF24 medium，float64）精修，种子和设置与 `pars refine` 相同。`p1_round8` / `p1_committor2`
与这两个逐字节相同，跳过。在本机（181，kasuga01）CPU 上并行，每个进程 2 线程，同时 14 个。

脚本：`docs/experiments/pars_ground_truth.py`（全量精修，可续跑，可按候选单独跑），
`docs/experiments/pars_ground_truth_score.py`（对照打分）。

### 3.1 三轮结果

| 轮次 | 设置 | 找到 0.419 eV 鞍点 | 单个耗时中位数 / 最大 | band 弧长中位数 / 最大 | 结果目录 |
| --- | --- | --- | --- | --- | --- |
| 1 | `eb`，k = 0.1（原代码） | 4/27；22 个报 `NEB_CONVERGED`，势垒 0.001–0.018 eV | 2601 / 3150 s | 2.42 / 3.0 Å | `runs/pars_ground_truth_eb_k0.1/` |
| 2 | improved tangent，k = 5，顶点分辨率检查 | 25/27；1 个如实 `FAILED`（NEB 1000 步未收敛），1 个新的假成功 | 535 / 2962 s | 1.69 / 183 Å | `runs/pars_ground_truth_it_k5_unaligned/` |
| 3 | 再加去除刚体运动 | **27/27**，势垒 0.4191–0.4192 eV | **311 / 418 s** | **0.94 / 0.9 Å** | `runs/pars_ground_truth/` |

### 3.2 每轮查出的原因

**第 1 轮：弹簧太弱。** 27 条 band 的两端都确实是两种互变异构体（O4–H 0.99 Å → O0–H 0.99 Å，端点
能量对称），端点检查没问题。问题在 band：假成功的 band 里，q_PT 在相邻两个像点之间从 +0.73 直接
跳到 −0.71（间距约 0.7 Å），鞍点附近（q ≈ 0）一个像点都没有。种子本身没问题，每个种子都有 1–5
个像点落在过渡区（|q| < 0.4），像点是在普通 NEB 里滑进盆地的。k = 0.1 时剩余的弹簧力
k·Δgap ≈ 0.05 eV/Å 已经满足 fmax，band 不会再被重新分布。

- 修复：`refinement_spring_k` = 5 eV/Å²（VTST Tools 默认值，SPRING = −5），配投影弹簧
  （`improvedtangent`，Henkelman–Jónsson），弹簧只沿切线方向，k 只影响间距、不切角。
- 修复：退回普通 NEB（不开 CI，或 CI 被拒）时，若最高像点两侧的间距超过 band 中位间距的
  `refinement_max_top_gap_ratio`（= 2）倍，报 `FAILED`，不给势垒。判据只看顶点附近，两端间距大
  不算（端点弛豫离原始帧较远时两端本来就大）。真实 band `t000029`、`t000002` 的间距序列作为测试
  的字面值。
- 合成复现：平坦盆地 + 0.42 eV 窄势垒，种子大部分弧长花在盆地里来回振动、穿越很短。k = 0.1 时
  普通 NEB "收敛"在 0.291 eV（顶点间距比 3.0），k = 5 时为 0.418 eV（1.0）。

**第 2 轮：种子里的刚体运动。** `t000002` 报 `CI_NEB_CONVERGED`、势垒 0.0002 eV：9 个像点全在盆地里，
q 在两种互变异构体之间来回翻（A A A B B A A A B），相邻像点间距均匀 1.85 Å（总弧长约 15 Å，正常
band 约 2.5 Å），间距检查和 CI 接受判据都拦不住。逐原子看：Kabsch 对齐后相邻像点要么完全相同
（0.00），要么差一次质子跳跃（0.63）；质心沿 band 漂了约 1 Å。也就是 MD 帧在种子中间发生了平移和
转动（扰动会给分子动量），而两个弛豫后的端点之间几乎没有刚体差（对齐前后 RMSD 0.26 / 0.24 Å）。
能量对刚体运动不变，所以笛卡尔坐标下的 NEB 把弧长花在刚体运动上。`t000038` 的 183 Å 弧长是同一
问题的极端情况。

- 修复：`refinement_remove_rigid_motion`（默认 True）：重采样前把所有帧对齐到第一帧（写种子和
  后端里都做），NEB 开 `remove_rotation_and_translation`，端点"同一极小"的 RMSD 只在开启时对齐。
- **限制：只对能量在刚体运动下不变的体系成立（孤立分子；PRRS 这些 run 的 scope 是
  `finite nonperiodic`）。核心嵌在固定环境里（QM/MM）或有任何外场时，必须设为 False。**
  RATIONALE 和 README 里都写了；测试里的单原子势（外场）都显式关掉了。
- 合成复现：只依赖原子间距离的三原子 A–H–B 质子转移势，种子带"出去又回来"的转动和平移。修复前
  CI "收敛"在 0.0003 eV，修复后 0.41947 eV（解析值）。

**第 3 轮**：27/27 收敛到同一个鞍点，没有新问题。

### 3.3 选择规则打分：这份数据分不出好坏

27 条路径的势垒相差只有 0.05–0.06 meV（同一个对称质子转移鞍点），distinct saddles = 1。默认规则、
四个去掉一个目标的变体、随机抽取，regret 全是 0、鞍点覆盖全是 100%。打分脚本把差值低于
`--barrier-tolerance`（0.005 eV）的势垒视为相等，并在势垒差低于这个容差时不输出 ρ（是噪声）。

| run | 默认规则选出的代表 | 实际够用 |
| --- | --- | --- |
| `p1_round7` | 8 / 15 | 1 |
| `p1_committor1` | 9 / 12 | 1 |

注意：这些 run 在审查 PARS 时已经被看过，在它们上面比较任何规则都是事后的，不是判据。

## 4. 前沿大小：从 5 变 12 再到 10，不是 bug

- 只修 Kabsch（B1）时 `p1_round7` 前沿还是 5；修 h³ 与自由段（B5）后变成 12，原因是曲率当时和
  帧数、recrossing 重复（B6）。
- 改到像点分辨率后为 10/15（`p1_committor1` 8/12）。蒙特卡洛：4 个相互独立的目标、随机数据，
  前沿平均 9.1/15、7.8/12。真实数据正好在这个基线上，所以这是 4 个独立目标、十几个点时 Pareto
  的正常表现，压不下去不是实现错了。
- 中途有一个错误判断已更正：原型里两次重采样得到 "12 → 5"、ρ ≈ 0，实际单次重采样的实现是 10，
  没有为了数字好看去换估计方法。
- 曾经怀疑是酶环境导致的：不是。`p1_round7` 是气相丙二醛（9 个原子，`active_atoms` 为空）；QM/MM
  下 MLP 只负责核心，环境原子不在这些坐标里。

## 5. 没改、需要你决定的

1. **反应管聚类（原第 5 条）**：真实数据里每个通道都只有 1 个管。9 原子、8 根键的体系，反应物和
   产物之间的 Jaccard 距离只有 0.222，低于阈值 0.30，管根本分不开。可行的方向是只在"发生过变化的
   键"上算距离，但抖动的键会把管切得过碎。config 规定默认值要在看 benchmark 之前冻结，所以没在
   这些 run 上调。
2. **压缩靠什么**：前沿（目标独立，n 小）和管（分不开）两层现在都不压缩。可选方向：重新考虑哪些
   目标进前沿（比如 `trigger_severity` 是 PRRS 从离散网格里选的扰动幅度，不太算路径本身的性质；
   `p1_committor1` 去掉它前沿从 8 降到 3），或者让 cover 承担压缩并定它的容差。都应在设计层面定好、
   冻结后再上数据。
3. **多鞍点数据**：要检验选择，需要找有多个鞍点或多个通道的 run（b05、connect 系列等），先查清
   哪些 run 真有多个鞍点。

## 6. 文件

- 修改：`pars/` 下的 `config.py`、`schema.py`、`cli.py`、`adapters/prrs_run.py`、
  `paths/{candidate,ingest,resample}.py`、`frontier/{archive,cover,dominance,margin,objectives}.py`、
  `refinement/ase_neb.py`、`README.md`，以及 `tests/` 下的 conftest 和各测试文件。
- 新增：`docs/experiments/pars_ground_truth.py`、`docs/experiments/pars_ground_truth_score.py`、本文。
- 结果：`runs/pars_ground_truth*/`（gitignored）。
- 新的 config 键：`refinement_endpoint_min_rmsd`、`refinement_spring_k`、
  `refinement_max_top_gap_ratio`、`refinement_remove_rigid_motion`，各带 RATIONALE。
