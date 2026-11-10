# 交接：HS 精修规则 v6 / v7、D2-HF-v1.2–v1.4 与并行（2026-11-06 至 11-08）

写给新开的窗口。本文件接替 [HANDOFF_HF_MIDPOINT_2026-11-06.md](HANDOFF_HF_MIDPOINT_2026-11-06.md)（那份留着没动，是否标"已作废"未问用户）。这里只汇总现状和去处，数字以链接的原文和 DONE_LOG 为准。11-08 18 时前后整理，**当前没有本线的作业在跑**。

## 必读（开工前全部读完）

- `docs/TS_HESSIAN_FREE_DESIGN_RECORD_2026-10-24.md`：设计稿。§4 F4（同键图跳支没有连续性证据）、§5.1（反力零点二分 / 受保护割线）、§5.2（跳支与滞回）、§5.3（以实测分辨率决定终止）、§7（HF-5：更软 / 协同体系）是本线的依据与去向。
- `docs/HANDOFF_HF_MIDPOINT_2026-11-06.md`：11-06 之前的背景（v4 / v5、新 E2 鞍点、候选身份与 HC）。
- `docs/HF_MIDPOINT_FINE_STAGE_OPTIONS_2026-11-06.md`：9c 剩余精阶段卡点的测量与选项。P1 已施工为 v6，Q2 未定。
- `docs/PLAN_HF_MIDPOINT_COARSE_CANDIDATE_2026-11-06.md`、`docs/PLAN_HF_ENDPOINT_CANDIDATE_2026-11-07.md`：HS 括区规则 v6（粗阶段出候选）、v7（端点出候选）的规则、用户决定、验收。
- `docs/HF_HS_NO_MAXIMUM_OPTIONS_2026-11-07.md`：R3 的测量。HS 的大部分 forward 花在同键图的脊上（F4）。R3b、R3c 未定。
- `docs/PLAN_HF_SAME_GRAPH_RIDGE_2026-11-07.md`、`docs/PLAN_HF_D2_V13_2026-11-07.md`、`docs/HF_SAME_GRAPH_RIDGE_SPEND_OPTIONS_2026-11-07.md`、`docs/PLAN_HF_D2_V14_2026-11-07.md`：D2-HF-v1.2 / v1.3 / v1.4 的规则、用户原话、测试与整跑结果（v1.4 计划 §6 是最新结果）。
- `docs/DONE_LOG.md` 2026-11-06 至 11-08 各节：每一步的命令、数字、更正；**11-08 末节有数据被删的记录**。
- `docs/TODO.md` 第 9 条及"未确认"：进展与仍未验证的推断。
- 全局与项目 `CLAUDE.md`：NFS 缓存（跨主机结果在写入那台读）、PBS 规矩、只按 PID 杀进程、zsh 不分词。

## 1. 用户原话（按日期；各计划里有上下文）

- 2026-11-06："顺着handoff往下做的 同时 审视 这几个性能优化问题"（附 6 条建议：割线、FIRE 速度投影、预条件、近根联合修正、按已访问分支复用、换坐标）。
- 2026-11-06：9c 全量诊断「跑」；选项文档「可以」；问"P1 可以定下来吗？"「可以」。
- 2026-11-07：9c'' 诊断「可以」；问"R1 可以定吗？"「可以」；问 R3 选项文档「可」；问"R3a 用 R3a-2（降优先级、不拒跑），可以吗？"「可以」；9d「准备」；v1.3 整跑「跑完了」；问同键图的脊剩余花费选项页「可」。
- 2026-11-07（并行与少跑）："你看啊 我们是这么想的 我们是做 小分子 的  反应设计是吧""不是  是 这样的 你看 我们显卡 本身就很空 你同意吧""是的 并行起来（ 并行异构都行 你说呢？""并行异构  +少跑  不就行了吗""S3 是不是 合理啊（ 为啥我看 s3 合理一些？"；问"要我按'S3 + HS / HC 并行'写施工计划吗？"「可以」；施工「可以 （你看181 空的吗」。
- 2026-11-08：整跑并发「8」、yayoi05「装得下」；"可以继续跑 几组计算 现在节点空"；"A 跑  然后 BC  我全给你3090 跑"；"a 你先落盘"；"是的准备交接文档"；"然后计算完了一个是吧"；"跑完了吧"；"处理完 整理完 写handoff"。
- **没有回答的**：问"配对回归验收（v7）可以先不跑吗？"；问"这两份草稿（v1.3、v1.4 结果）的运行目录只加小文件、不加 journal.jsonl，可以吗？"。
- 2026-11-08（结果目录被拿掉）："是这样 的影响1  和2 是因为块太大了  我git会爆炸 所以没动"；"我是因为推完了所以 才给她拆了好了"。

## 2. 代码状态

### 2.1 HS 括区规则（`explore/src/explore/ts_ef.py`，均默认关闭）

| 版本 | 开关（`EFSettings`） | 内容 |
| --- | --- | --- |
| v3 | — | 单阶段中点 1e-3（默认） |
| v4 | `midpoint_coarse_transverse_fmax_eV_A` + `midpoint_sign_bound_eV_A` | 两阶段中点（b05：2e-3 / B 0.025） |
| v5 | v4 + `midpoint_fine_optimizer="fire"` | 精阶段 FIRE |
| v6 | v5 + `midpoint_coarse_candidate` | 粗阶段几何门、反力差门都过且一侧满足候选判据 → 出候选 |
| v7 | v6 + `endpoint_candidate` | refine 两端未分辨括区时，一端满足候选判据且另一端反力指向最大值 → 出候选 |

### 2.2 D2 规则（`explore/src/explore/d2.py`，在 `init_run` 冻结；`run_b05_d2.py` 开关在括号里）

| 版本 | 内容 |
| --- | --- |
| D2-HF-v1.1（`--hf --start-cluster`） | 起点分簇、两轮派发；10-26 至 10-29 的整跑 |
| D2-HF-v1.2（`--same-graph-ridge`） | 同键图的脊（扫描峰点与下一点 `edges` 相同）在第一、二轮之后派发 |
| D2-HF-v1.3（`--hs-bracket`） | v1.2 + 每个 HS 用 v7 括区规则（`HSBracket` 冻结 b05 的 2e-3 / 0.025） |
| D2-HF-v1.4（`--sample-parallel`） | v1.3 + S3：同键图的脊按 `audit_fraction` 抽样，其余 `HS_skipped:same_graph_ridge`；HS / HC 按组分波次 `run_many` 并行（共用批量求值器，`max_active` / `wave_size` 沿用冻结值），按输入顺序结算，HC 判重在结算时按顺序做 |

- 测试：`test_ts_ef test_ts_branch test_ts_connect test_ts_warm_start test_d2_hf test_d2 test_d2_resume test_probe test_run_many test_scan_performance` 166 passed（11-07）。每版都有变异检验。**全包 `make check` 这条线一直没跑。**
- 已知测试缺口：v1.4 的"HC 判重被放进并行阶段"测试台上测不出（测试台没有重复事件）。

## 3. 已查清的事（按时间）

| 步骤 | 结论 | 记录 |
| --- | --- | --- |
| 9c 剩余精阶段卡点 | 两类：粗阶段已有候选仍被要求精松弛（v4 里 Claude 自加的"候选让位"）；精阶段 1e-3 在这些帧上低于实测单次求值噪声，FIRE 已停止移动。FIRE 速度法向分量可忽略（排除建议第 2 条） | DONE_LOG 11-06 |
| v6 验收（PBS yayoi40） | 对照 15 项候选保留；精阶段失败 8 → 1 / 32；救回的都是新 E2；forward −18% / −8% | DONE_LOG 11-07 |
| 9c'' 普查（v4 两块） | `hysteresis` + `walk_relax_unconverged` 占 HS forward 约 95%，两侧反力都不翻转；`no_resolved_maximum_bracket` 的端点帧已满足候选判据 | DONE_LOG 11-07 |
| v7 定向冒烟 | 1 项由端点出候选（新 E2）；配对回归没跑；三次整跑中端点出候选一次也没触发 | DONE_LOG 11-07 |
| R3 测量（八块） | 成键不变的"脊"占 ridge 类 HS forward 92%，0 个候选；52 个候选与 22 个新 E2 意图的脊都有成键变化 | DONE_LOG 11-07；选项文档 |
| v1.3 整跑（yayoi32，与 v4 同机） | 候选 8 / 6 → 21 / 18；事件 sA 3 → 3、sB 2 → 3（sB 首次找到新 E2）；同键图的脊 363 个 HS 0 候选、约 190k forward / 块 | DONE_LOG 11-07 |
| 并行测量（三次：4080S 一次、3090 两次） | 串行 GPU 利用率约 17–23%；p8 比串行快 4.3–4.9 倍；并行与串行的结局差异和串行自身的 fp32 波动同量级（串行对串行 4–7 项不同） | DONE_LOG 11-08；`hs_parallel_bench/repeats_summary.json` |
| 批量求值噪声（181） | 批量与单结构的力噪声、平均值看不出差别（4 个对照候选结构） | DONE_LOG 11-08 |
| **v1.4 整跑 vs v1.3 同硬件对照**（yayoi05 / yayoi08，3090，布局与参数相同） | 事件两版两块都是已知三类，无新类、无丢失；候选 21 / 17 对 21 / 19；refine forward 约少 72%（69k / 65k 对 238k / 238k）；墙钟 3.6 h 对 9.0 h；同键图的脊样本 16 / 16 个，候选 0 | DONE_LOG 11-08 末节；v1.4 计划 §6 |

## 4. 作业与机器

- 本线没有作业在跑（11-08 17:57 核对：三个 PBS 作业都有 done 行）。
- 183 的 GPU：11-08 08:15 起被别人的作业（`experiments/2026-11-08_ppi_r4/`，PID 520460）占了约 8.3 GB 显存，不是本线的，不要动。之后是否还在，没再核对。
- PBS 结果按惯例在 183 上读，作业日志 `<作业名>.o<作业号>` 在项目根目录。

## 5. 结果目录被拿掉（11-08 12:02–12:04；用户知情，与 git 体积有关）

- 用户 11-08 原话："是这样 的影响1  和2 是因为块太大了  我git会爆炸 所以没动"。下面"影响"一段里的 1、2 即指此事。又："我是因为推完了所以 才给她拆了好了"——推送之后才从工作区拿掉，不是丢失。

- 不见了：`benchmarks/explore_d2/results/2026-10-26/` 整个目录（10-26 的 b05 HF 两块整跑、冒烟、`repeat_response_census` 等），以及 `results/2026-11-08/b05_hf_v14_smoke`（目录、日志、pid 文件，只剩 `_prepare.jsonl`）。
- 同一时间段，`benchmarks/explore_d2/` 下多个脚本与结果目录的修改时间被改；核对过本线的代码、测试、计划、TODO，内容未变。
  - **更正（11-08 18:59）**：上面这句不对。`explore/src/explore/ts_branch.py`、`ts_connect.py` 和 `explore/tests/test_ts_branch.py`、`test_ts_connect.py` 也在这段时间不见了，`explore.d2` 现在导入失败。见 DONE_LOG 11-08 末节。
- 没有做任何恢复。
- 影响：
  - 参照集从八块变为六块（已知三类都还在）；
  - 文档里引用这些路径的地方（DONE_LOG 10-26 至 11-07 的若干节、各选项文档的"八块"统计）仍是删除前的结果，数字不受影响，但文件已无法复查；
  - 4 份旧提交草稿（`2026-10-26_1200`、`2026-10-26_1800`、`2026-10-28_0245`、`2026-10-28_1200`）的 `git add` 引用了已删路径，照原样执行会失败（未改这些旧草稿）；本线 `2026-11-08_0800` 草稿已去掉被删的冒烟路径。

## 6. 接下来（都需要用户定）

| 事项 | 现状 / 依据 |
| --- | --- |
| v1.4 是否作为今后 D2 的默认规则 | b05 上事件不变、refine forward 少约 72%、墙钟快约 2.5 倍；只有 b05 一个体系的证据 |
| 换体系检验 | 设计 §7 HF-5（更软 / 协同体系）；也是检验 S3（"同键图的脊不出候选"）是否普遍成立的唯一办法 |
| 提交草稿里的大文件 | `2026-11-07_1530_d2-hf-v13`、`2026-11-08_1800_v14-results-handoff` 的 `git add` 含整块运行目录（`journal.jsonl` 127–379 MB / 个，合计约 1.8 GB），会撑大 git；改法待用户定 |
| R3b 单项 HS 早停、R3c 跳支两侧局部鞍点搜索（`branch_bracketed` 跳支跑了也没产出） | 未定 |
| Q2 精阶段"噪声受限"出口 | 未定；v6 之后精阶段失败已很少 |
| 9c' HC `endpoint_unstable` / `tangent_no_resolved_drop` | 每次整跑 1–4 个，未查 |
| 并行的下一步 | 求值外约 100 s 在 p4–p16 不降，瓶颈像是线程下的 Python 侧；下一步更可能是多进程而非流水线 |
| 用户 11-06 所贴建议剩余：第 1 条割线、第 3 条预条件（延续阶段）、第 5 条按已访问分支复用 | 未做；第 2 条已排除；第 4 条暂缓；第 6 条即 R3c |
| 9e 反力差门 | 用户决定保留；仍偶尔挡下 |
| TODO 第 11 条（备份被饿死） | v1.3 / v1.4 下第二轮不再被暂缓，b05 上可能已不成问题，未专门查 |
| TODO 第 12 条（Hessian 交叉核验门限）、第 3、4、5、7 条 | 未动；12 待用户定改法 |
| 全包 `make check` | 没跑 |

## 7. 提交草稿（用户提交，按顺序；都在 `.commit_drafts/`）

上一份交接所列各草稿 → `2026-11-06_1752_hf-midpoint-identity-hc` → `2026-11-07_0010_hf-midpoint-coarse-candidate` → `2026-11-07_0930_hf-endpoint-candidate` → `2026-11-07_1200_d2-same-graph-ridge` → `2026-11-07_1530_d2-hf-v13` → `2026-11-08_0800_d2-hf-v14` → `2026-11-08_0840_handoff` → `2026-11-08_1800_v14-results-handoff`。

- 几份草稿的 `git add` 有重叠文件（`ts_ef.py`、`d2.py`、测试、DONE_LOG、TODO、交接、README 等）。`git add` 加的是工作区当前内容，所以先提交的草稿会带上后面的改动；要分清历史，就一次合并提交或按顺序连续提交、中间不改文件。
- 见 §5：4 份 10-26 / 10-28 的旧草稿引用了已删路径。

## 8. 环境与操作要点（都已实测）

- **结果在写入那台读**：183 写的在 183 上读；本机读要等写完 6 h 以后（NFS 属性缓存）。PBS 节点写的结果按惯例在 183 上读。
- **183 列目录会重复**：汇总一律按文件名去重（各 `*_summary.py` 已这样做）。
- **zsh 不对 `$var` 分词**：把多个文件名放进一个变量传给 pytest 会变成一个参数，结果是"no tests ran"（11-07 踩过）。
- **远程后台启动**：`ssh … "bash -c 'cd X && { setsid nohup CMD > log 2>&1 < /dev/null & echo \$!; }'"`，PID 记进同目录的 `pids.txt`。
- **PBS**：照模板 `/home/kasuga/.claude/templates/pbs.zsh` 只改 ← 行。双卡节点（yayoi03 / 05 / 08，2×3090，EPYC 7452 32 核）上两个进程分别设 `CUDA_VISIBLE_DEVICES=0 / 1`、各 16 线程。只交冒烟跑通过的代码。v1.4 整跑现成脚本：`benchmarks/explore_d2/b05_hf_v14_pbs.zsh`。
- **事件身份核对**：`candidate_identity.py` 的 `references` / `cluster`（对称 RMSD ≤ 0.15 Å 且 ΔE ≤ 0.05 eV），参照用剩下的六块（10-28 v1.1 / v1.1r、10-29 v4）；读 6 块 journal 在 183 上约 10 分钟。
- **变异检验**：变异代码放临时目录，`PYTHONPATH=<临时目录>:../src` 配合 `-o pythonpath= -p no:cacheprovider`。
- **脚本索引**（`benchmarks/explore_d2/`）：
  - 9c 诊断：`midpoint_fine_diag.py`（+ `_run.sh`、`_summary.py`）；
  - v6 验收：`midpoint_coarse_candidate_validate.py`（+ `_run.sh`、`_summary.py`、`_pbs.zsh`）；
  - v7 冒烟：`endpoint_candidate_smoke.py`；
  - 普查：`hs_prewalk_census.py`、`hs_walk_profile.py`、`hs_ridge_continuity.py`；
  - v1.2 重放：`same_graph_ridge_replay.py`；
  - 候选身份：`candidate_identity.py`（文件名规则 v3–v7）；
  - 并行测量：`hs_parallel_bench.py`（+ `_pbs.zsh`、`_rep_pbs.zsh`）；
  - 批量噪声：`batch_noise.py`；
  - 整跑：`run_b05_d2.py`，`b05_hf_v13_pbs.zsh`、`b05_hf_v14_pbs.zsh`、`b05_hf_v13r_pbs.zsh`。

## 未确认（Claude 推断）

- 并行的瓶颈是"线程下 Python 侧串行化"：依据只是求值外时间在 p4–p16 不降，没有做 profile。
- 批量噪声只测了 4 个对照候选结构，失败组卡住的帧没测。三次测量重复的结局差异与 fp32 波动同量级，是"并行不额外改变结局"的主要证据。
- 同键图的脊"不出候选"只有 b05 一个体系的证据（八块回放 0 / 1 663，两次整跑 0 / 363 + 0 / 32 样本）。
- v1.4 的 refine 节省与墙钟缩短中，S3（少跑）与并行各占多少分不开（v1.3r → v1.4 一次改了两样）；按测量推算，HS 并发 8 单独约快 4.5–4.9 倍。
- 端点出候选（v7）在三次整跑中一次也没触发；它的价值目前只有 1 项定向冒烟作证据。
- 11-08 v1.4 与 v1.3r 的候选数差（sB 17 对 19）归为 fp32 波动，依据是串行重复之间候选数也差 1–2 个，未逐项核查。
