# 施工计划：Non-Hessian AFIR 第一阶段（MC-AFIR + 本项目 MIN / TS / IRC）

日期：2026-11-14。依据：[设计稿第 3 版](NONHESSIAN_AFIR_DESIGN_2026-11-14.md)（用户 2026-11-14 审定）。状态：施工中。

## 必读（开工前全部读完）

- `docs/NONHESSIAN_AFIR_DESIGN_2026-11-14.md` 全文：分工（§1）、公式（§2）、路径（§3）、MIN / TS / IRC（§4）、MC-AFIR（§5.1）、参数（§7）。
- `explore/src/explore/ts_ef.py`（`_Constraint`、`refine_ts_ef`）、`ts_branch.py`（`continue_branches_ef`）、`ts_connect.py`（`branch_starts`、`trace_sides_ef`）：TS 与连接直接调用。
- `explore/src/explore/d2.py` `_hc_run`：HS 结果 → 连接的现成接法（照抄其逻辑）。
- `benchmarks/explore_d2/l0.py`、`explore/src/explore/poses.py`：随机起始结构的来源。

## 用户原话

- 2026-11-14：问"设计稿可以按这一版定稿、我接着写第一阶段（MC-AFIR）的施工计划吗？"，回答"可以 就这么施工 （然后 和我们原先的 路线基本都是高度耦合的是吧"。

## 1. 范围

第一阶段只做 MC-AFIR（两片段之间推近，ρ = +1）：随机起始结构 → AFIR 路径 → 路径上的 MIN → 路径最高点处的 TS（HF 内核）→ 两侧连接；"最近 n = 100 条路径无新产物"或预算用完即停。SC-AFIR、网络扩展、拉开（ρ = −1）放第二阶段。旧 D2 不动。

## 2. 代码

| 文件 | 内容 |
| --- | --- |
| `explore/src/explore/afir.py` | `WeightedDistance(A, B)`：ρ_AB 的 `value` / `gradient`（解析）/ `pairs`，接口同 `probe.Intent`；`alpha_from_gamma`、`gamma_from_T_t`；`AFIRCalculator`：包住无偏计算器，给出 E + ρα·ρ_AB 及其力，同时保存无偏 E / F |
| `explore/src/explore/afir_path.py` | `run_afir_path`：在偏置能面上做 E/F 极小化（ASE BFGS），逐步记录无偏 E、ρ_AB、键图、几何；终止：偏置能面极小、经过两个能量极大、解离；10 步早弃；定义域守卫；步数与 forward 上限 |
| `explore/src/explore/afir_local.py` | 路径 → MIN（无偏极小化 + 扰动确认）；路径 → TS（最高点前后各取若干帧，固定 ρ_AB 做无偏横向松弛得 ScanPoint，找反力 + → − 的相邻对，交 `continue_branches_ef`）；TS → 连接（`branch_starts` + `trace_sides_ef`，同 `d2._hc_run`） |
| `benchmarks/afir_mc/run_mc_afir.py` | 运行器：体系表沿用 `benchmarks/explore_d2/systems.py`；随机起始结构由 L0 姿态生成器逐条生成；每条路径一个 JSON + 一个 ledger（已存在即跳过，可续跑）；停止判据；汇总 |
| `explore/tests/test_afir.py` | 见 §4 |

ρ_AB 的梯度：w = (c/r)^6、c = R_i + R_j，∂(w r)/∂r = −5w、∂w/∂r = −6w/r，所以 ∂ρ/∂r_ij = w_ij (6ρ/r_ij − 5)/Σw，再乘键向量的单位向量。

## 3. 第一阶段我定的细节（设计稿未规定，列出来，你可以改）

| 项 | 取值 | 理由 |
| --- | --- | --- |
| 路径优化器 | ASE BFGS，maxstep 0.1 Å | 只用 E/F；GRRM 的 RFO 用 Hessian |
| 偏置能面极小判据 | 偏置力 fmax ≤ 0.03 eV/Å（`config.quench_fmax_eV_A`） | 沿用项目现有的下降判据 |
| 路径步数上限 | 1000 步 | 防止不收敛的路径吃掉预算 |
| "能量极大"的判定 | 轨迹上 E 的局部极大，两侧下降都超过 noise_eV（v1.5 规则：2 × fp32 能量间距） | 不把 fp32 抖动当成极大 |
| 解离 | 键图连通分量数比起点增加，且新分量之间最近原子对 r/(R_i+R_j) > 3.0 | 与远距原子对上限 3.0 一致 |
| 定义域守卫 | 任一原子对 r/(R_i+R_j) < 0.5，或单步无偏能量下降 > 5 eV，记 `domain_exit` 并停 | AFIR 力可能把原子压到模型训练域外 |
| TS 括区取帧 | 路径最高点前后各 4 帧（按路径顺序） | 括区要覆盖最高点两侧 |
| TS / 连接的数值设置 | v1.5 的模型级常数（σ 0.0042、粗阈值 2e-3、B 0.025、FIRE 精阶段等）、hs_cap / hc_cap 3000 | 沿用已有的 HF 内核设置 |
| MIN 确认 | 3 次 0.03 Å 扰动后重新下降，键图相同且 ΔE、对齐 RMSD 在判重容差内 | 同 HC 端点检查 |
| 判重 | 键图 + 对称 RMSD（结构自身图群 ∪ 反应物图群）+ 镜像 | 设计稿 §12 |

## 4. 测试（解析势，CPU）

1. ρ_AB 的梯度对有限差分；平移、转动不变。
2. `alpha_from_gamma(100 kJ/mol)` ≈ 118.8 kJ/mol/Å；`gamma_from_T_t(298.15 K, 10 天)` ≈ 106.9 kJ/mol（Chem. Rec. 原文数值）。
3. `AFIRCalculator` 的力对有限差分；无偏 E / F 与底层计算器一致。
4. TwoSaddle 解析台：H 与 O1 之间推近（A = {H}、B = {O1}），路径到达产物盆；最高点处的 TS 括区 + HF 内核找到解析台的鞍点（`SADDLE_Y`），完整力收敛；连接两端是两个极小。
5. 早弃、定义域守卫、步数上限、预算截断都有明确结局码；外部计数 = 账本。

## 5. 运行

1. 测试通过后，b04 冒烟一次（181，几条路径），确认 AFIR 路径、TS、连接与账目。
2. 冒烟后给计划表（三个体系，预算与 D2 v1.5 的 B* 相同以便对照），等你给节点。

## 未确认（Claude 推断）

- §3 全部是我定的；其中路径优化器与 maxstep 会影响路径形状，进而影响找到哪些 TS。
- 式 (3) 仍未对照原文 PDF；实现按设计稿 §2 的推导形式。
- L0 姿态生成器作随机起始结构，与原文的"随机朝向与接近位置"不完全相同（先做构象资格与壳层摆放）。
