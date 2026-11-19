# 交接：Non-Hessian AFIR（2026-11-14 至 11-15）

只写这条新线。旧的 D2 / HF-5 / HF-6 路线已不是主线，只作对照数字的来源，不在这里展开（需要时看 `docs/HANDOFF_HF5_HF6_2026-11-12.md`）。

## 必读（开工前全部读完）

- `docs/NONHESSIAN_AFIR_DESIGN_2026-11-14.md`：设计稿（用户 11-14 审定）。§0 为什么绕开 Hessian、§1 照搬 AFIR 与本项目自己做的分工，是这条线的全部前提。
- `docs/PLAN_AFIR_MC_2026-11-14.md`：第一阶段施工计划，§3 列了我定的细节。
- `docs/ALGORITHM_DEFECTS_2026-11-14.md`：为什么放弃 D2 的发射层（发射点按几何枚举，84% / 65% / 84% 的发射白花）。
- `docs/DONE_LOG.md` 2026-11-14 末起各节：每一步的命令、数字、修正。
- 全局与项目 `CLAUDE.md`：NFS 缓存（结果在写入机读）、PBS 规矩、只按 PID 杀进程。
- `docs/AFIR_PUSH_OPTIONS_2026-11-15.md`：N1 推不出 DA 的机制（AFIR 加权 ρ 塌到一对 H···H、把成键 C–C 拉开）与三个改法选项，**待用户定**。
- 本文 §8（11-15 下午更新）：§4–§6 的部分内容已被它取代。

## 1. 用户的方向（原话，2026-11-14 至 11-15）

- "我说了我们要通用性 算法"；"通用性三个字 你听懂 啥意思吗"——一套算法、一套参数，b04 / b05 / N1 只作检验，不拿来调参。
- "发射 的 地方 根本不在  反应附近的 没用发射"；"减少  这个花销是 可以通过算法约束的是吧  你看 典型的 就是 AFIR 那套"。
- "但是我们要做  non-hessian 的  afir"；"我们IRC 有自己一套  TS 也有自己一套。。。  但是其他 能 吸收的  对于反应 逻辑  都是 一样的"；"对面是 高贵的  hessian  贵 我mlp便宜"。
- "因为 mlp 给的F 和E 靠谱但是  hessian 不靠谱  所以只能绕开  就是mlp 特化的"；"你这叫防守 写法"——MLP 的 Hessian 不可靠是前提，不必另证。
- "直接用原文参数"（γ = 106.9 kJ/mol）；"106.9 kJ/mol 相当于 25 kcal   意思是 基本反应不了"（无催化室温上限）；远距原子对上限"往下压一压  至少压倒3倍"（→ 3.0）。
- "mlp 探索的时候  是不同于 afir  我们可以大量并发 因为 分子量太小了"。
- 不升 γ，改用 SC-AFIR 原子对分组："b   走起来"。

## 2. 方法现状（一句话版）

AFIR 的反应逻辑照搬（加权距离 ρ、α(γ)、MC / SC-AFIR 分组、路径停止判据）；路径用批量 FIRE（只用 E/F）；极小用 E/F 扰动确认；TS：在路径上取帧按固定 ρ 松弛找反力 + → − 的括区 → 批量收窄 → HF 反力零点内核；IRC：E/F 两侧下降。全程不调用 Hessian（阶数检验只在事后、单独记）。

## 3. 代码（旧 D2 未改）

| 文件 | 内容 |
| --- | --- |
| `explore/src/explore/afir.py` | ρ_AB（`WeightedDistance`，接口同 `Intent`）、α(γ)、γ(T, t)、`AFIRCalculator` |
| `explore/src/explore/afir_batch.py` | 批量路径 `run_paths`、批量（约束）极小化 `minimize`、`confirm_minima`、括区 `bracket_frames` / `scan_points` / `pick_bracket`、批量收窄 `narrow_brackets`；`AFIR_MAX_ATOMS_PER_CALL` |
| `explore/src/explore/afir_fragments.py` | SC-AFIR：扰动结构、A / B 分组、候选原子对（≤ 3.0，按对称去重） |
| `explore/src/explore/afir_path.py`、`afir_local.py` | 逐条版路径与 MIN / TS / 连接（解析台测试用；运行器已改用批量版） |
| `benchmarks/afir_mc/run_mc_afir_batch.py` | 运行器：`--mode mc|sc`、`--sc-poses`、`--early-discard-steps`（默认关）；按波次、可续跑；事件按 `reaction_event_key` 计（超限时退回摘要键） |
| `benchmarks/afir_mc/afir_validate.py` | 事后检验：阶数 + 参照匹配（含镜像）+ 不同 TS / 事件 + 与 D2 v1.5 的 R(b) |
| `benchmarks/afir_mc/throughput.py` | 吞吐随并发数 |
| 测试 | `explore/tests/test_afir.py`、`test_afir_batch.py`、`test_afir_fragments.py`：20 passed |

## 4. 结果

| 体系 | 结果 | 出处 |
| --- | --- | --- |
| b04 | SC 两个种子都找到参照 SN2 鞍点（−26132.2012 eV，一阶，连接连通）；首次命中约 9.7k / 9.8k forward，D2 v1.5 为 29.6k / 31.3k；R(½B*) = 1 | `benchmarks/afir_mc/results/2026-11-14/overnight/validate.json`（fulda） |
| b04（整片段 MC 推） | 36 条路径 0 次 SN2 → 改用 SC 分组的原因 | DONE_LOG 11-14 |
| b05 | 过夜三项因事件键 bug 崩溃（已修）；11-15 冒烟：27 条路径 6 条键图改变、1 个事件，括区阶段太贵、预算用完 | `results/2026-11-15/b05_smoke_fix/` |
| N1 | 过夜两个种子因早弃误杀（约 30/32 条每波）0 个事件（早弃现已关）；MC 64 条 0 次 | `results/2026-11-14/overnight/` |
| 吞吐 | POLAR-1 每结构 4–7 ms，GPU 在并发 8–32 时饱和；约 2 000 原子 / 调用 OOM（32 GB 卡） | `results/2026-11-14/throughput_*.json` |

## 5. 在跑 / 待提交

- `benchmarks/afir_mc/afir_sc_b05_n1_pbs.zsh`（yayoi42）：b05 / N1 SC 两个种子 + 检验，输出 `results/2026-11-15/pbs/`。脚本已交给用户 `qsub`；**提交与否、作业号本交接写作时未知**，接手先问用户或查 `afir_sc.o*`。预计 3–5 h（无实测依据，按 fulda 估计外推）。
- fulda 上没有本线进程在跑（11-15 冒烟已结束）。

## 6. 待办（按我认为的先后）

1. 读 PBS 结果（fulda 上读），看 b05 / N1 在关闭早弃后能否推出反应与 TS；N1 的协同 DA 是这条线的关键检验。
2. 括区扫描太贵（b05 冒烟约 3 万 forward 只为括区）：现在对整条路径逐帧松弛；可改为先稀疏后加密或只扫最高点附近一段。改前先量。
3. SC 只做了推近（ρ = +1）；拉开、网络扩展（从新极小继续，导航 μ_i）是设计第二阶段，未做。
4. 式 (3)（α 由 γ 换算）仍未对照原文 PDF；式 (5) 在 403.15 K 例子上与原文不符（145.6 对 153.0）。
5. TS 的成本按波末累计记（上界）；要精确的首次命中成本，需要在波内记录。

## 7. 机器

- fulda（`ssh cyanus@192.168.4.199`）：RTX 5090 32 GB，24 核；本线冒烟与过夜都在这里，结果在这里读。远程启动用脚本文件 + `setsid nohup`，命令末尾加 `< /dev/null` 并给 ssh 设 `timeout`，否则 ssh 会挂住。
- yayoi42（PBS）：RTX 4080 SUPER 16 GB → `AFIR_MAX_ATOMS_PER_CALL=512`。

## 未确认（Claude 推断）

- b04 上整片段 MC 推不出 SN2 的原因（力按权重分摊、落在 Cl–C 上不足）没有专门检验。
- 早弃误杀的原因（FIRE 从静止起步 vs GRRM 按 RFO 步长定的 10 步）是推断；关闭早弃后成本会升多少未测。
- 4080 SUPER 上 512 原子 / 调用不会 OOM 是按 5090 实测按显存比例推的。
- 退回的事件键可能把不同通道合并，b05 的事件数因此可能偏少。

## 8. 2026-11-15 下午更新（接手先读这一节；§4–§6 中与此冲突处以此为准）

### 用户原话（11-15）

- 我问"要我在 fulda 上量一下 b05 括区点的横向松弛收敛曲线吗？"，答"跑"。
- 我问"要我把上限改成 2000，然后在 fulda 上用 b05 冒烟一次吗？"，答"可以 然后我们 继续 看看 整体overview  路走通了吗"。
- "可以 总而言之 大改思路 通的剩下就是  怎么精细设计算法了"（同意做 N1 只读分析）。
- "对面跑完了   然后一起重放"；"做 然后我出去吃饭"（QM 单点）。
- 我问"要我把……几种通用改法列成选项文档给你定吗？"，答"是的  我们要改一改 afir 同意吧"。
- 选项文档写完、问是否按建议施工后："你先落盘 然后我想想"——**选项尚未选定，不要施工**。

### 已完成（数字与命令见 DONE_LOG 2026-11-15 各节）

| 项 | 结果 | 出处（读的机器） |
| --- | --- | --- |
| PBS 28756（yayoi42，09:35 → 12:05） | b05 SC sA / sB 各 2 个事件（E2、甲基 C–H 去质子，fallback 键）、0 个 TS；N1 SC sA / sB 0 个事件（1 379 条 / 种子，24.9 万 / 27.2 万 forward）；validate R 全 0 | `benchmarks/afir_mc/results/2026-11-15/pbs/`（fulda） |
| b05 0 个 TS 的原因 | 括区扫描的固定 ρ 横向松弛 200 步内只 6% 收敛，`pick_bracket` 要求两点都收敛；放宽到 2000 步 93/95 收敛（中位 539 步），第 200 步的反力符号 16/95 会翻转，不能凑合用 | `results/2026-11-15/bracket_convergence/`（fulda） |
| 修复 | `run_mc_afir_batch.py` `SCAN_MAX_STEPS = 2000`（只用于括区扫描与收窄；HF 内核仍 200） | 代码 |
| b05 冒烟（修复后） | 12 次括区尝试 5 次成括区（改前 0）；内核 2 次 `converged_candidate`，其中 E2 一条 `connected`；2 次内核自己的 `transverse_not_converged`（同类 200 步上限，未改）、1 次 `hysteresis`；括区阶段约 11k forward / 次（贵） | `results/2026-11-15/b05_smoke_scan2000/`（fulda） |
| N1 机制 | 推 DA 位置的路径把一对片段间 H···H 压到 0.93–0.99 Å（权重份额 2–4% → 42–46%），能量 +1.0–1.3 eV，成键 C–C 仍 ≥ 3.2 Å；末端帧 AFIR 力把四个 DA 成键 C–C 对拉开 | `results/2026-11-15/n1_da_pairs_sA.txt`、`n1_da_replay/`（fulda） |
| QM 核对 | ωB97M-V/def2-TZVPD 单点：上述压近的能量升高 QM 1.289 / 1.014 / 1.349 eV，MLP 1.315 / 1.017 / 1.335 eV——模型没问题，"MLP 偏软"的推断已否定 | `results/2026-11-15/n1_hh_qm/`（**183**） |

### 现在的状态

- 没有进程在跑（fulda、183、PBS 都已结束）。
- 选项文档 `docs/AFIR_PUSH_OPTIONS_2026-11-15.md` 待用户定：选项 1 SC 只推单对（括区坐标需同时改，否则回到缺陷 C）；选项 2 ρ 中不计 H（k / l 为 H 时除外）；选项 3 双分子先按原版 MC-AFIR、停止按 100 条计。我建议 2 + 3 的停止判据。

### 待办（按我认为的先后；第 1 条等用户定）

1. 按用户选定的选项施工 → 三个体系各一个种子同条件重跑（看 N1 是否在 DA 位置成键、b04 / b05 的事件与 TS 是否保持）。
2. HF 内核横向松弛的 200 步上限在 b05 冒烟中造成 2/5 括区 `transverse_not_converged`；是否同样放宽，未定、未测。
3. 括区成本：改为从最高点往回逐点松弛、找到括区即停（原 §6 第 2 条）。
4. 冒烟中的 E2 TS 未做参照匹配与阶数检验。
5. 原 §6 第 3–5 条不变（第二阶段拉开 / 网络扩展；式 (3)、(5) 原文核对；首次命中成本按波内记）。

### 本节新增的脚本

`benchmarks/afir_mc/bracket_convergence.py`、`n1_da_pairs.py`、`n1_da_replay.py`，启动脚本 `fulda_bracket_convergence.sh`、`fulda_smoke_scan2000.sh`、`fulda_n1_da_replay.sh`；QM 脚本 `results/2026-11-15/n1_hh_qm/n1_hh_qm.sh`（在 183 上写入）。

### 未确认（Claude 推断）

- 原版 AFIR（GRRM，RFO + QM）在环戊二烯 + 丙烯醛上从随机朝向能否推出 DA，项目内无证据。
- b04 整片段 MC 推不出 SN2 是否也是同一塌缩机制，未查。
- 冒烟里 E2 两条路径 ΔE 差 0.65 eV 却落在同一 fallback 键下，是否同一通道未核实。

## 9. 2026-11-16 至 11-19 更新（接手先读这一节与 DONE_LOG 11-16 起各节）

- 已完成：PBS 29729（yayoi08）三体系两种子，11-15 定稿默认（单对 + 双对层、TS 坐标解耦、不停两个极大、括区逐轮展开）。N1 第一次由搜索找到参照协同 DA 鞍点（s-trans 两个，全部来自双对层）；b05 两种子三类参照全找到；b04 两种子找到 SN2。数字见 DONE_LOG 11-18。
- N1 两个 s-cis 参照没找到：反应物里只有 s-trans 丙烯醛（只读核实）；第一阶段推 O9···C6 已得到 s-cis 极小但未进入搜索。第二阶段网络扩展计划 `docs/PLAN_AFIR_NETWORK_2026-11-18.md`，**待用户定：展开的 MIN 是否都跑双对层**；未施工。
- 进行中（用户"你直接disown明早看"）：CI-NEB 作 TS 定位的对照，180 上冒烟 1 条路径，见 DONE_LOG 2026-11-19 末节；结果在 180 上读。fulda 用户在调试模型，**不要占 fulda**。
