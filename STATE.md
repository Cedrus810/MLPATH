# MLPATH 当前状态（交接用，2026-11-26）

写给压缩上下文之后的自己 / 下一个窗口。细节与数字的出处都在 `docs/DONE_LOG.md` 2026-11-20 → 11-26 各节。

## 必读（开工前全部读完）

- `CLAUDE.md`（项目）：设计稿表，新增"TS 定位器（成键窗口版）"一行；环境与跨机读结果规矩。
- `docs/PLAN_FREEZE_AND_HELDOUT_2026-11-26.md`：**当前主线计划**（冻结 AFIR-W1 → 开发集回归 → 留出验证 → b09 第 3–6 级）。
- `docs/TS_LOCATOR_DESIGN_2026-11-19.md` §2.6、§7：现行默认定位器（成键窗口）的诊断、证据、回归；剩余待定项。
- `docs/B09_DISTANCE_2026-11-20.md` §2、§4、§7：b09 终极目标的分级距离、验收判据草稿（待审）、现状与瓶颈。
- `docs/NONHESSIAN_AFIR_DESIGN_2026-11-14.md` §0、§2、§7、§11：E/F only、AFIR 推力规则、参数全体系同一套（§4.2 定位器已被取代，原文保留）。
- `docs/DONE_LOG.md` 2026-11-20 → 11-26：所有实测数字、PID、结果目录、更正。
- `docs/DESIGN_DRAFT_TPS_MSM_2026-11-26.md`：独立工作线，**暂存不施工**，只需知道它存在。

## 用户决定（逐字，均为 2026-11 本对话）

- 11-20："我们项目撞车了  需要加速了"（Popcornn，arXiv:2611.22272，全文在 `docs/refs/popcornn_2611.22272/`）；"而且我们不做 几千个反应 就把那个比较神经病的最终goal 做出来"；"主要是我们最后有一个 B09 还是11 哪个奇怪的 反应类型" → 目标是 b09。
- 11-22："2  你别  浪费时间 180 继续冒烟"（未冒烟代码先交 PBS）。
- 11-26："整体整理 整理  可以 可以" → **成键窗口定位器设为默认**，HF 内核作对照臂。
- 11-26："可以   做起来"（冻结 + 留出计划）；"做起来"（D1 / D2 / 构象回流按建议施工）。
- 11-26："这个草稿丢到独立工作线路 留着  然后 回到当前 情况 我们继续看"。

## 代码现状（`benchmarks/afir_mc/run_mc_afir_batch.py` 等）

| 开关 | 默认 | 含义 |
| --- | --- | --- |
| `--locator {window,hf}` | **window** | window = `explore/src/explore/afir_window.py`（String k=`--window-k` 10 → 键图变化处窗口 [j−1, j+2] → 16 节点 String 20 次 → 固定方向爬坡 → ±0.05 Å 两侧 E/F 下降）；hf = 旧括区 + HF 内核。11-26 前的运行都是 hf，复现 / 续跑须加 `--locator hf`（三个旧正式 PBS 脚本已补） |
| `--sc-split {0,1,2}` | 0 | 1 = A：片段内原子对只在第一个节点推；2 = 另加"伙伴分子在接触范围"的节点（无证据需要，未跑） |
| `--skip-start-min` | 0 | C：扰动起点键图未变时不极小化起点 |
| `--confirm-dedup` | 0 | D：已知事件不确认；同新事件一波只确认一个代表 |
| `--peak-dedup D` | 0 | 只对 hf 定位器；d* = 0.25 Å（用户 11-22 定）；window 模式不生效 |
| `--double-xx {all,first}` | all | D1：XX 双对（两对都在片段内）在复合物上只在第一个节点推 |
| `--double-order {shuffle,class}` | shuffle | D2：双对按 II → IX → XX 排序 |
| `--reflow N` | 0 | B：一层构象回流，每个新片段构象摆 N 个朝向，只推片段间单对 + II/IX 双对；记 `reflow.json` |

其它：NFS 列目录重复计数已修（`sorted(set(glob))`）；单分子输入（片段数 1）已支持；`systems.py` 有 `b09`（无偏输入）与冒烟专用 `s01`；`classify_b09.py` 已改对称 RMSD。新脚本：`peak_window_refine.py`（任意体系成键窗口，输出交 `afir_validate.py`）、`b09_ts_summary.py`、`b09_saddle_check.py`、`b09_chain_diag.py`、`peak_ts_distance.py`、`double_layer_census.py`。测试：`explore/tests/test_afir_window.py`（+ `test_afir_string.py`，4 passed）。

## 关键结果

- **window 回归**（PBS yayoi15，11-25→26，`benchmarks/afir_mc/results/2026-11-25/window_regress/`）对 PBS 29729：b04 首次命中 23 518→13 336（sA）、76 424→19 024（sB）；b05 三参照 7.3/16.3/42.9 万→2.0/3.4/6.7 万，全中，总花费 70.6 万→20.9 万；N1 召回不变（2/4 含镜像），endo 89.5 万→65.6 万；b09 sB 自动定位"反应物→(2+1)"一阶鞍点（事后 Hessian 1 个负模 −520 cm⁻¹；与 TS-2OMe 同类、另一面进攻；POLAR-1 上比 TS-2OMe 低 0.25 eV；对称 RMSD 0.34 Å）。TS-2OMe 本身、-cis、-C2 仍未找到。
- **双对层普查**（`results/2026-11-26/double_layer_census/`）：XX 六次运行 0 事件 0 参照；II 出 N1 全部 DA 参照；IX 出 b05 全部消除参照。N1 / b09 预算大头在双对层。
- **W1 候选冒烟**（N1，4 朝向，`results/2026-11-26/w1_smoke180/`）：无报错；II 优先排序后 ts_exo_s-trans 在 66 218 forward 命中（节点少，不可与正式运行直接比）；回流未触发，回流摆放路径用构造的 s-cis 单独检查通过。
- **单分子冒烟**（s01）：准备 / 路径 / 极小化段通过，定位段未覆盖。
- 更正记录：b09 卡宾"构象不对"是原子编号假象（已更正）；b05 event2 丢失**不能**归因规则 A（换硬件后 SC 节点集合不同）。

## 下一步（按序）

1. **W1 开发集回归**：`benchmarks/afir_mc/w1_regress_pbs.zsh` 已写好，`HOST/NCPUS/NGPUS/MEM/THREADS/W1_DATE` 占位待填；串行 b04 sA → b05 sA → N1 sA → b09 sA → b09 sB，W1 全开（64 朝向、回流 32）。耗时估 30–40 h（按 window 回归外推，无实测）。11-26 用户给 yayoi31（9950X + RTX 4080 SUPER 16 GB，32 / 64 GB）：已填 ncpus 32、ngpus 1、mem 54gb、线程 16、`AFIR_MAX_ATOMS_PER_CALL=512`，输出 `benchmarks/afir_mc/results/2026-11-26/w1_regress/`；qsub 由用户执行（作业号待记）。11-27 用户给 yayoi10（EPYC 7452 + 2×3090，32 / 128 GB）：`w1_regress_sB_pbs.zsh` 跑 sB 种子（n1 → b04 在卡 0，b05 在卡 1），输出 `results/2026-11-27/w1_regress_sB/`。
2. 回归通过 → 冻结 AFIR-W1，写入项目 `CLAUDE.md` 与计划文档。
3. 留出：RGD1（2023-04-24 更新）H2 / H3 单分子槽位，按 `paper/HORIZONTAL_VALIDATION_DESIGN_2026-11-13.md` §3 选取与隔离；下载前报数据量。
4. b09 第 3–6 级：-cis（构象）、网络扩展（(2+1)→(4+1) 重排）、射击分支比（动量系综问题待用户定）。

## 待用户定

- b09 验收判据草稿（`B09_DISTANCE` §4）。
- TS 定位设计 §7 第 2–4 项（k 默认已取 10）。
- 射击阶段能否用 Hessian（准经典采样）还是只用热动量。
- 双对层是否放开共用原子（D3，建议不做）。

## 运行与机器注意

- 现在**没有我启动的进程在跑**。
- PBS 结果在 fulda（`ssh cyanus@192.168.4.199`）或 183 上读；fulda 11-25 曾 ssh 连接被重置，改用 183 读成功。180 / 183 的结果在写出它的机器上读。
- 11 GB 卡（180、183 的 2080 Ti）必须设 `AFIR_MAX_ATOMS_PER_CALL`（512 或更小；183 曾因默认 1024 OOM）；24 GB 卡 768；两进程共卡要减半。
- 远程长时间命令一律 `setsid nohup … &` 脱离 ssh，PID 写进结果目录 `pids_*.txt`；本地 `timeout` 会断开前台 ssh（11-24 出过一次，按 PID 清理过）。
- 登录节点有人扫全盘时 NFS 极慢（11-22），远程进程会卡在 D 状态，等它恢复即可。
- MLPATH 不 commit、不建分支（memory）。

## 未确认（Claude 推断）

- W1 回归耗时 30–40 h 是外推。
- 回流能否在正式运行里推出 N1 s-cis / b09 -cis 构象并定位到对应鞍点：没有证据。
- 64 朝向能否让 b09 sA 推到卡宾碳–Cb：几何估算单朝向约 11.5%，sA 实际 32 个 0 次，原因未查。
