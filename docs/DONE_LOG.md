
## PBS 脚本：`--locator window` 回归（yayoi15，2026-11-25；用户给节点）

`benchmarks/afir_mc/window_regress_pbs.zsh`（照模板；ncpus 32 / ngpus 1 / mem 118gb）。单卡 4090，**串行**（用户"串行"）：b04 sA → b09 sB（第一阶段 A/C/D + window）→ N1 sA → b05 sA → b04 sB；16 线程、`AFIR_MAX_ATOMS_PER_CALL=768`。b04 / b05 / N1 与 PBS 29729 同种子 / 壳层 / 预算 / 32 朝向，只多开 `--locator window`。最后 `afir_validate.py`（b04 / b05 / N1）。输出 `results/2026-11-25/window_regress/`。`qsub` 由用户执行。
- 更新（同日，用户"串行"）：改为单进程串行 b04 sA → b09 sB → N1 sA → b05 sA → b04 sB，16 线程、`AFIR_MAX_ATOMS_PER_CALL=768`；上一条"两路共卡"的设置作废。

## 已完成：`--locator window` 回归（PBS，yayoi15，串行，2026-11-25 18:1x → 11-26 13:30；结果在 183 上读）

`results/2026-11-25/window_regress/`；对照 PBS 29729（旧 HF 定位器，b04 / b05 / N1）与 b09 第一阶段（A/C/D + HF）。

| 运行 | 首次命中（forward）旧 → 新 | 全部参照（新） | 总花费 旧 → 新 | 墙钟（新，4090 单进程） |
| --- | --- | --- | --- | --- |
| b04 sA | SN2 23 518 → **13 336** | 1/1 | 159 985 → 47 509（条目跑完即停） | 867 s |
| b04 sB | SN2 76 424 → **19 024** | 1/1 | 165 919 → 54 142 | 475 s |
| b05 sA | event0 73 263 / event2 163 069 / event1 429 297 → **20 298 / 34 083 / 66 520** | **3/3** | 706 090 → 208 912 | 2 589 s |
| N1 sA | endo s-trans 894 511 → **656 105**；exo s-trans（镜像）1 677 450 → 1 572 846 | 2/4（含镜像；s-cis 两个仍缺，构象原因） | 1 676 490（同预算） | 17 300 s |
| b09 sB | 无鞍点 → **"反应物 → (2+1)" 一阶鞍点**（第 7 波结束累计 101 271） | 参照 0/4（≤ 0.2 Å） | 1 699 974 → 1 699 975 | 47 506 s |

- b09 sB（`benchmarks/afir_mc/b09_ts_summary.py`，新，只算几何，183 CPU）：收敛候选聚成 2 个不同鞍点：① `wave_0007:n28|2-15`，两侧 (2+1) / 反应物、`connected`，即 11-25 离线找到的同一个点（到 TS-2OMe 对称 RMSD 0.34 Å，C:–Cb 1.98）——这次由运行器自动得到；② 卡宾碳–原子 0 成键的鞍点（7 次，`n3|0-2` 一类，非参照）。TS-2OMe 本身、-cis、-C2 未找到。事件 4 个（第一阶段 3 个）。
- 读法：新定位器在 b04 / b05 上首次命中提前 1.8–6 倍、b05 三个参照全中，总花费降到约 1/3；N1 召回不变、首次命中提前；b09 首次自动定位到 (2+1) 通道鞍点。N1 / b09 的预算大头仍在双对层（N1 1.32M、b09 1.59M）。

## 决定 + 整理：成键窗口定位器设为默认（2026-11-26；用户原话"整体整理 整理  可以 可以"）

- `run_mc_afir_batch.py`：`--locator` 默认改为 `window`；`hf`（括区 + HF 内核）为对照臂。11-26 之前的运行都用 hf，复现或续跑须加 `--locator hf`；已给三个旧正式 PBS 脚本（`afir_pair_double_pbs.zsh`、`b09_phase1_pbs.zsh`、`rules_regress_pbs.zsh`）显式补上 `--locator hf`（不改任何结果）。
- 文档：[TS 定位设计](TS_LOCATOR_DESIGN_2026-11-19.md) 状态行改为"已定"、新增 §2.6（诊断、改法、b05 / b09 证据、正式回归）、§7 第 1 条标为已定；[Non-Hessian AFIR 设计稿](NONHESSIAN_AFIR_DESIGN_2026-11-14.md) 开头加指针（§4.2 定位器被取代，原文保留）；[b09 距离](B09_DISTANCE_2026-11-20.md) §7 第 2 级改为"部分"并列出新瓶颈；[总览](PROJECT_OVERVIEW.md) 加 11-26 现状；[文档入口](README.md) 改行；项目 `CLAUDE.md` 设计稿表加一行（TS 定位器）。
- 仍待用户定：构象回流 B、朝向数 E、双对层规则（是否放开共用原子 / 是否只在片段间单对里组合 / 预算占比）、b09 验收判据草稿、TS 定位设计 §7 第 2–4 项（k 的默认已取 10，窗口规则不再用固定最高节点）。

## 施工 + 冒烟：单分子输入支持（2026-11-26；用户"可以   做起来"；180，PID 1617，结果在 180 上读）

[计划](PLAN_FREEZE_AND_HELDOUT_2026-11-26.md) §3.0。另写计划本身，与双对层普查（`benchmarks/afir_mc/double_layer_census.py`，`results/2026-11-26/double_layer_census/census.json`，只读：XX 组合六次运行 0 事件 0 参照，II 产出 N1 全部 DA 参照，IX 产出 b05 全部消除参照）。
- `run_mc_afir.prepare` 允许 1 个片段；`run_mc_afir_batch.py` 片段数为 1 时不摆朝向（`--mode mc` 报错），输入结构是唯一 SC 节点。`systems.py` 登记冒烟专用 `s01`（丙二醛烯醇，RDKit ETKDG + MMFF，`benchmarks/afir_mc/inputs/s01_malonaldehyde.extxyz`，9 原子；不是基准）。
- 冒烟 `benchmarks/afir_mc/k180_k1_smoke.sh`（预算 30 000，默认 window 定位器；输出 `results/2026-11-26/k1_smoke180/`；638 s）：准备、单对层 24 条、双对层 74 条都跑通、无报错；只有 1 条双对路径成键改变，但其终点极小化未收敛（`window_no_ends`），**定位段在单分子上没有走到**。冒烟通过（准备 / 路径 / 极小化段），定位段未覆盖。输入是 RDKit 的默认构象，未必是有分子内氢键的那种烯醇式，这只影响能否推出反应，不影响代码路径。

## 施工 + 冒烟：AFIR-W1 候选规则 D1 / D2 / 构象回流（2026-11-26；用户"做起来"；180，PID 1912，结果在 180 上读）

[计划](PLAN_FREEZE_AND_HELDOUT_2026-11-26.md) §2。`run_mc_afir_batch.py` 加三个开关（默认 = 原规则）：`--double-xx {all,first}`（D1：两对都在片段内的双对组合在复合物上只在第一个节点推；回流节点从不推）、`--double-order {shuffle,class}`（D2：双对层按 II → IX → XX 排序，类内仍按种子打乱；条目记 `pair_class`）、`--reflow N`（B：单对层后，从键图不变的终点里按片段对称 RMSD（片段自同构 + 镜像，> basin_rmsd_A）找新构象，与另一片段的输入构象摆 N 个朝向作新节点 `r*`，只推片段间单对（层 `rwave_*`）与 II / IX 双对；只一层；记录 `reflow.json`）。`sc_items` 加 `node_prefix` / `include_input` / `intra` / `ledger_name` / `k_offset` / `n_poses` 参数供回流调用。
- 冒烟 `benchmarks/afir_mc/k180_w1_smoke.sh`：N1 sA，window + A/C/D + D1 + D2 + `--reflow 2`，4 个朝向，预算 80 000；输出 `results/2026-11-26/w1_smoke180/`；2 177 s。全程无报错。单对层 3 波（23 166 forward），76 个终点键图不变，**回流没有找到新构象**（0 个新节点，`rwave` 未运行：回流的摆放 / 节点路径未覆盖）；双对层按类排序，第 1 波（II）出 DA 事件，`afir_validate.py`：**ts_exo_s-trans（非镜像）首次命中 66 218 forward**，2 个不同一阶鞍点。
- 读法：只有 5 个节点（4 朝向 + 输入），单对层条目远少于 33 节点的正式运行，首次命中数字不能直接与 PBS 29729（endo 894 511）比；排序让 II 双对先付费这一点在这里起了作用。回流的摆放路径要靠更多节点（正式运行）或单独构造的测试来覆盖。
- 回流摆放路径的单独检查（CPU，不调用模型，scratch 脚本）：把 N1 输入中丙烯醛绕 C7–C8 转 180°（键图不变）构造 s-cis 终点，`reflow_poses` 识别为 1 个新构象（片段 0，来源 `n0|9-14`），未转的终点未被误判，摆出 3 个朝向（形状正确）。
- 预备：`benchmarks/afir_mc/w1_regress_pbs.zsh`（照模板，节点待填）：开发集 b04 sA / b05 sA / N1 sA / b09 sA / b09 sB 串行，W1 候选规则全开（64 朝向、回流 32），对照 11-25 window 回归与 PBS 29729。

## 文档：TPS + MSM 设计草稿（独立工作线，暂存；2026-11-26；用户"这个草稿丢到独立工作线路 留着"）

新增 [DESIGN_DRAFT_TPS_MSM_2026-11-26.md](DESIGN_DRAFT_TPS_MSM_2026-11-26.md)，`docs/README.md` 加索引行。不施工，不进当前主线计划。
