# 交接：b05 HF D2 工作线（2026-10-26 至 11-01）

写给接手这条线的窗口。本文件只汇总现状和去处，细节与数字以链接的原文为准。

## 必读（开工前全部读完）

- `docs/TS_HESSIAN_FREE_DESIGN_RECORD_2026-10-24.md`：这条线的设计稿（HF TS 搜索的问题、F1–F8 缺陷、§5.2 跳支与滞回、§7 施工顺序）；所有改动都以它为准。
- `docs/PLAN_HF4_D2_2026-10-26.md`：D2-HF-v1 的调度、HS / HC 任务、账户与续跑合同；本线的改动都建在它上面。
- `docs/TODO.md`：用户原话、待办编号（1–12）与每条的现状；接手后先从这里挑下一步。
- `docs/HF_START_CLUSTER_RESULTS_2026-10-28.md`：D2-HF-v1.1（同一起点两轮派发）的 b05 正式结果，含另一窗口的核查与我的补充。
- `docs/HF_WALK_POLICY_RESULTS_2026-10-29.md`：延续试步失败的处理（v3 折半重试 → v4 首次失败即终止）、同机对照、v4 整跑、新 E2 鞍点的发现与检验；现行策略在此。
- `docs/DONE_LOG.md` 2026-10-26 之后各节：每一步的命令、数字、更正，以及 QM 核验全过程。
- 全局与项目 `CLAUDE.md`：机器、NFS 缓存（跨主机结果在写入那台读）、PBS 规矩、只按 PID 杀进程。

## 1. 现行代码状态（`explore/`）

| 部分 | 现状 | 记录 |
| --- | --- | --- |
| 去重 | `certify.pair_symmetries`：鞍点 + 端点图自同构群并集；certify 第 5 步与 `d2._hf_duplicate` 都用它 | DONE_LOG 10-28 |
| 调度 | `D2-HF-v1.1`：`StartCluster(noise_multiple, first_round)`，同一起点的 HS 两轮派发，第二轮额度不够记 `HS_deferred:same_start`；`run_b0{4,5}_d2.py --start-cluster` | `PLAN_START_CLUSTER_2026-10-28.md` |
| 跳支括区 | `_hf_bracket` 用正向扫描的跳变区间（反应物侧最后一点、产物侧第一点），不再用 L2 反扫点（旧括区宽度为 0） | DONE_LOG 10-28 |
| 延续试步失败 | `HF-EF-branches-v4`：`_walk` 中试步 `transverse_not_converged` 首次失败即终止该侧（`termination = "relax_unconverged"`），另一侧继续；无最大值且有一侧如此终止时总状态 `unresolved:walk_relax_unconverged`；几何跳变的折半定位保留 | `HF_WALK_POLICY_RESULTS_2026-10-29.md` §5 |
| 规则版本 | 上面三处修复都是原地修复，没有新 rule_version；引用 10-28 及以前的 HF 运行要注明用的是旧括区 / 旧延续 | DONE_LOG |

相关测试最近一次：`test_ts_branch test_ts_ef test_ts_connect test_ts_warm_start test_d2_hf test_d2 test_d2_resume` 94 passed；全包 `make check` 一直没跑。

## 2. 关键结果

- b05 HF 整跑三次（yayoi04 v1.1 / yayoi08 v3 / yayoi32 v4，同种子，硬件各不相同，L1/L2 与分簇都有差别，不是只差代码的对照）：v4 那次 HS 结算 267 / 272、候选 8 / 6、连通事件类 3 / 2。
- 冻结输入同机对照（183）：30 个失败组项上 v3 / v4 = 52.6k / 26.8k forward，候选 0 / 0；同机噪声约 ±400 forward / 项。
- **新 E2 鞍点**（v4 sA 意图 `k1 / distance:10-12`）：POLAR-1 面上比已知两类低 0.39 eV；内部模复核判一阶；QM（r2SCAN-3c，OptTS + Freq + IRC + 端点优化 + ωB97M-V/def2-TZVPD 单点）确认为一阶 E2 TS（虚频 −127.1 cm⁻¹），IRC 两端为收敛的 E2 反应物 / 产物复合物；比已知 E2 TS 低 11.75 kcal/mol（ωB97M-V），其中约 8.2 来自它的反应物复合物更低，按各自复合物的势垒 +4.08 对 +7.59 kcal/mol（电子能，参照态不同）。来源不是 v4 策略，是二分中点在该硬件上碰巧收敛（属 TODO 第 9 条那类）。
- 一阶性检验的门（`src/prrs/runner.py` `_analytic_hessian` 的交叉核验门限 0.01）在 fp32 下贴着噪声：同一几何 20–60% 概率被判“解析不可信”并回退到更差的 fd Hessian。未改（TODO 第 12 条）。

## 3. 待办（按 `docs/TODO.md` 编号）

| 编号 | 内容 | 状态 |
| --- | --- | --- |
| 9 | 二分中点阈值（1e-3）与 fp32 噪声：43 项已夹住最大值却卡在中点；新鞍点就出自这类 | 未做；先修诊断脚本（`hs_transverse_diag.py` 统一用了 0.003，中点实为 1e-3），再定是否调阈值并在该阈值重测 σ |
| 11 | HC 失败后同类备份在第二轮被暂缓“饿死” | 未处理 |
| 12 | 解析 Hessian 交叉核验门限贴噪声（PRRS 共享代码） | 测量已做，等用户决定改法 |
| 10 端点部分 | 端点松弛失败（没有可回退的已接受点） | 未处理 |
| 3、4、5、7 | 筛选 vs 直接 TS 预算对比、构象入口复用、暂缓状态、单项 HS 早停 | 未开始 |

## 4. 提交草稿（用户提交，按此顺序；都在 `.commit_drafts/`）

`2026-10-26_1800_hf4-d2` → `2026-10-28_0245_dedup-pair-symmetries` → `2026-10-28_1200_d2-hf-v11-start-cluster` → `2026-10-28_1800_hf-jump-bracket` → `2026-10-28_2000_hs-transverse-diag` → `2026-10-28_2200_hf-walk-retry` → `2026-10-29_0230_hf-walk-first-failure` → `2026-10-29_1700_hf-walk-policy-docs` → `2026-11-01_1400_b05-qm-and-saddle-checks`。后面的草稿改的文件与前面的重叠（`d2.py`、`ts_branch.py`、`DONE_LOG.md` 等），顺序不能乱。

## 5. 环境与操作要点（都已实测）

- **ORCA**：183 本地 `/home/korakuen/orca_6_1_1_linux_x86-64_shared_openmpi418_avx2/`（PATH 与 LD_LIBRARY_PATH 加该目录，系统 Open MPI 5.0.11 可用）；PBS 节点用 `module load orca/6.1.1`（带出 `openmpi/5.0.7_gcc13.3.0`，ORCA 在 `/home/apps/orca_6_1_1_avx2/`）。181 上的 `orca` 是编辑器 orca-ide，不是量化软件。
- **PBS 上跑 MPI**：select 必须写 `mpiprocs=N`，否则只给 1 个槽位，`mpirun -np >1` 报 “not enough slots”。纯 CPU 作业放 GroupB（yayoi13 / 14，48 核 192 GB）；QM 脚本模板见 `docs/experiments/b05_qm/b05_qm_ep2_pbs.zsh`。
- **ORCA 优化“正常结束”不等于收敛**：默认 50 个循环上限，松散复合物常常走不完；链脚本要检查 “OPTIMIZATION RUN DONE”（`run_endpoint2.sh` 已这样做，`run_chain*.sh` 没有）。
- **变异检验**：`explore/pyproject.toml` 的 pytest `pythonpath = ["src", ...]` 会盖掉 `PYTHONPATH`；放在临时目录的变异代码要配合 `-o pythonpath=` 才会真的被测。进程内替换模块属性的做法不受影响。
- **NFS 缓存**：181 / 183 属性与目录缓存 6 小时；PBS 或 183 写的文件，本机可能看到旧内容或看不到新文件，读结果到写入那台读（PBS 结果可在 183 读，前提是 183 之前没缓存过该目录）。
- **zsh**：本机与远程默认 shell 是 zsh，`$VAR` 不分词；多文件参数用数组或 bash。

## 6. 当前没有运行中的作业

2026-11-01 核对：183 与 181 上没有本线启动的进程；PBS 作业 22053（yayoi13）已结束。

## 未确认（Claude 推断）

- 新鞍点 QM 虚频偏软（127i），模中混有 OH 转动与甲基运动；它对应的反应物复合物是否是该 OH⁻ 位置下的全局最低，没有查。
- “v4 让第二轮覆盖回升”的方向与冻结输入对照一致，但三次整跑硬件不同，不能把差别全部归因于策略。
- 第 12 条门限问题在其他体系 / 模型上的噪声水平没有测，改法要按体系重新测量。
