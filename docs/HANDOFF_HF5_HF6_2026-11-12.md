# 交接：HF-5 跨体系检验与 HF-6 Dimer 比较（2026-11-08 至 11-12）

写给新开的窗口。本文件接替 [HANDOFF_HF_D2_2026-11-08.md](HANDOFF_HF_D2_2026-11-08.md)（那份未改、未标作废）。只汇总现状和去处，数字以链接的原文与 DONE_LOG 为准。2026-11-12 13:30 前后整理，**当前没有本线的作业或进程在跑**（181 上本线进程已全部退出，PBS 作业均有 done 行）。

## 必读（开工前全部读完）

- `docs/TS_HESSIAN_FREE_DESIGN_RECORD_2026-10-24.md`：设计稿。§7（HF-5 / HF-6 的完成条件）、§8.1（首篇两层证据）、§4 F4 / F8 是本线的依据与去向。
- `docs/DISCOVERY_LAYER_DESIGN_2026-10-13.md` §3.1：意图规则对所有体系一次冻结，不为 DA / E2 写专门组合——用户 11-10 纠正方向的依据。
- `docs/EVALUATION_PROTOCOL_DESIGN_2026-10-13.md` §3（参照认证）与 DB 臂的约定：HF-6 的公平性约定出处。
- `docs/HF5_SUMMARY_2026-11-10.md`：HF-5 结论与需要一次性通用重定的规则。
- `docs/HF6_SUMMARY_2026-11-12.md`：HF-6 结论（局部与端到端）。
- `docs/PLAN_HF5_CROSS_SYSTEM_2026-11-08.md`、`docs/PLAN_HF6_DIMER_2026-11-10.md`、`docs/PLAN_HF6_DIMER_E2E_2026-11-11.md`：三份施工计划（设定与"未确认"）。
- `docs/DONE_LOG.md` 2026-11-08 末至 11-12 各节：每一步的命令、数字、更正（含 11-12"更正与完成"）。
- `docs/TODO.md` 第 12 条：检验脚本里改了 Hessian 交叉核验门限。
- 全局与项目 `CLAUDE.md`：NFS 缓存、PBS 规矩、只按 PID 杀进程、zsh。

## 1. 用户原话（按日期，关键决定）

- 2026-11-08："我们继续搞吧。。。"；问"要我写 HF-5（b04 → b09）的施工计划吗？"「可以」；问"可以按这个计划开工吗？"「可以」；ts_branch / ts_connect 被拿掉后「好了恢复了」；"停 / 183在跑别的"；"180"。
- 2026-11-09：问"B_b04 取 0.0237，可以吗？"「可以」；"可以 180 +181"；问 N1 σ / 几何噪声「可以」；"ahh nfs 爆炸了"；"别浪费时间了求求你了"；"你浪费时间没节点了 全部 别给我一个冒烟"；"去180 冒烟"。
- 2026-11-10："1 接二元意图 先查范围"；**"但是 是 我们为啥变成了 一个反应 优化一次参数 没意义了 啊"**；**"我们在干什么 / 你是不是完全搞糊涂了 方向"**；问"要我先写 HF-5 小结吗？"「对的  写」；「可以  继续做」；「开工」；"冒烟后 给我 我节点跑"。
- 2026-11-11：第二阶段计划「可以」，开工「可以」；"要不先把第二份 DB D2 做了 ？"；"3 做了"。
- 2026-11-12：「写 然后写 handoff」。
- **方向（用户 11-10 纠正后）**：不再逐体系调参数、补规则来让某个体系"跑通"；HF-5 / HF-6 是检验，负结果照实写。

## 2. 结论（详见两份小结）

| 项 | 结论 |
| --- | --- |
| HF-5 b04 | 6 个运行都找回参照 SN2；S3 的前提（同键图的脊不出候选）不成立，且依赖扫描步长 |
| HF-5 N1 | N1 参照：4 个协同 DA 一阶鞍点（endo / exo × s-cis / s-trans）；单原语 D2 只找到单键分步鞍点，DA 0 / 4 |
| HF-5 数值参数 | σ、几何噪声、B、noise_eV 各体系取法不统一；noise_eV 现定义在 N1 上测到构象盆能量差 |
| HF-6 局部 | 同一冻结括区：HF 命中参照 b04 11 / b05 19、Dimer 3 / 0，HF forward 少约 6–10 倍；N1 上 Dimer 2 项找到 DA，HF 0 |
| HF-6 端到端 | 同一 L0、同一预算、两个种子：D2 在 b04 / b05 找到参照（sB b05 2/3 严格口径），DB 三个体系都是 0；N1 两臂都是 0；D2 命中都在预算末尾 5–15% |

## 3. 代码（全部在 `benchmarks/` 下；`explore/src/` 与 D2 规则未改）

| 文件 | 内容 |
| --- | --- |
| `benchmarks/explore_d2/systems.py`、`run_d2.py` | 体系表（b04 / b05 / n1 的 HFSettings、HSBracket）与通用 runner（`--system`、`--rotate`、`--magnitude-step`、`--certify-max-rise-eV`） |
| `benchmarks/explore_d2/hf_event_validate.py` | 事件检验：阶数 + 参照匹配；`validation_config`（交叉核验门 0.1）、`match_group`（鞍点图群 ∪ 姿态图群） |
| `benchmarks/explore_d2/rf_jitter_scan.py` | σ / 几何噪声的通用测量（只取成键改变的脊） |
| `benchmarks/explore_d2/hf6_local.py`、`hf6_local_summary.py` | HF-6 局部定位两臂与汇总 |
| `benchmarks/explore_d2/hf6_e2e_db.py`、`hf6_e2e_summary.py` | HF-6 端到端 DB 臂与汇总（R(b)、首次命中成本） |
| `benchmarks/explore_d2/hf5_b04_summary.py`、`hf5_smoke_compare.py` | HF-5 b04 汇总、runner 一致性对照 |
| `benchmarks/n1_diels_alder/build.py`、`reference.py`、`input.extxyz` | N1 输入与 POLAR-1 参照鞍点 |
| `benchmarks/explore_d2/*_pbs.zsh`（`b04_hf5_*`、`b04_sign_resolution`、`n1_*`、`hf6_*`） | 本线各次 PBS 脚本 |

- `systems.py` 中 N1 的 `hs_bracket` 仍为 None（B_N1 测了两份但没有定）；N1 的 noise_eV 未定（并跑了 0.0527 / 0.00195 两种）。
- 测试：本线没有改 `explore/src/`，没有新增单元测试；`make check` 未跑。

## 4. 结果目录

- HF-5：`benchmarks/explore_d2/results/2026-11-08/`（b04 a0–a2）、`2026-11-09/`（b04 v1.4 五个运行、a5、N1 v1.2、σ 测量）、`benchmarks/n1_diels_alder/results/2026-11-09/`（N1 参照）。
- HF-6：`results/2026-11-10/hf6_local/`（局部，汇总以 `hf6_local_summary_v3.json` 为准）、`results/2026-11-11/hf6_e2e_db/`、`hf6_e2e_db_sB/`、`n1_hf5_v12_sB_n0020/`（端到端，汇总以 `hf6_e2e_summary_sA_v3.json`、`hf6_e2e_summary_sB_v2.json` 为准）。较早版本的汇总保留，b05 匹配不可信。

## 5. 接下来（都需要用户定）

| 事项 | 现状 |
| --- | --- |
| HF-5 小结 §4 的通用规则：S3 与同键图判定、数值参数的固定测量程序、noise_eV 定义 | 11-12 写成选项 [HF_GENERAL_RULES_OPTIONS_2026-11-12.md](HF_GENERAL_RULES_OPTIONS_2026-11-12.md)，待用户定；定了以后三个体系统一重跑 D2 |
| PRRS 默认的解析 Hessian 交叉核验门限 | 只在检验脚本里改为 0.1；PRRS 默认 0.01 是否改未定 |
| "内部模复核不一致"这道门 | 未查；使一阶数在重算之间摆动 |
| Dimer 找到的参照外一阶鞍点（局部 b05 8、N1 27；端到端 b05 9 / 12、N1 4 / 7） | 11-12 已查（DONE_LOG 同日末节、HF-6 小结 §6）：b05 有一组取代通道不在参照集里；N1 一组是参照 exo s-cis 的对映体（匹配规则漏了镜像，汇总未改） |
| N1 参照的 QM 一轴 | 未做（评估协议 §3.2 第 6 步） |
| 更多种子 / 置信区间 | 端到端只有 2 个种子 |
| D2 先扫描后精修，命中都在预算末尾 | 设计 §4 F8；未改 |
| b09 | 留作最终验收，HF-5 / HF-6 都没有用 |

## 6. 机器与操作要点（本轮实测）

- **180**：11-11 16:00 前后重启后 NFS 未挂上（`cd /home/ruigengji/MLPATH: No such file or directory`）；`nvidia-smi` 报驱动 / 库版本不匹配，但 torch 能用 CUDA；11-09 一次模型加载 1 480 s、NFS 炸时进程长时间 D 状态。用前先核对。
- **181（本机）**：汇总类小计算放这里；GPU 有时有他人作业。读远程写入的结果前，核对文件数与写入机或 fulda 一致。
- **fulda（`ssh cyanus@192.168.4.199`）**：用户 11-09 让在这里读 PBS 结果；能看到项目目录与 `openmm_dev`。
- **183**：11-08 用户说"183在跑别的"；本轮之后没有用。
- **PBS**：照模板只改 ← 行；qsub 由用户执行；只交冒烟跑通过的代码。同一张 16 GB 卡放 5 个进程时出现过一次 CUDA OOM（b04 v1.4 sA 重复的 1 个 HC），之后每卡 ≤ 3 个进程。
- **后台等待**：用 `kill -0 <PID>` 循环等进程退出；本轮两次把 PID 抄错，启动后先 `ps -o pid,etime,cmd -p <PID>` 核对再挂等待。zsh 下 `echo =====` 会报错。
- **NFS 新目录**：在远程机器上新建的目录，本机可能一时看不到（pid 文件写失败过一次），pid 文件在写入机上写。

## 7. 提交

本轮没有写提交草稿。本线新增与修改的文件见 §3 与各计划 / 小结；运行目录里的 `journal.jsonl` 很大（用户 11-08："块太大了 我git会爆炸"），提交时只宜加脚本、文档与汇总 JSON。上一份交接所列草稿的状态未变。

## 未确认（Claude 推断）

- N1 漏掉协同 DA 的原因是单根距离原语推出分步路径；HF-6 中 Dimer 从扫描峰点能找到 DA 与此一致，但没有专门检验。
- DB 端到端找不到参照，可能是原样 L0 姿态离反应区太远、默认 0.1 Å 位移太小；未区分。
- "S3 不通用""noise_eV 定义失效"各只有一个反例体系。
- 180 的 NFS 何时恢复不知道。
