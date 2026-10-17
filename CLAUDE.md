# MLPATH 项目规矩（项目级，补充全局 CLAUDE.md）

入口文档：[docs/PROJECT_OVERVIEW.md](docs/PROJECT_OVERVIEW.md)（整体逻辑与进度）、[docs/README.md](docs/README.md)（找文档）。
完成记录写进 `docs/DONE_LOG.md`；历史规矩在 `docs/HANDOFF.md` §2。

## 设计稿

| 工作线 | 设计稿 | 状态 |
| --- | --- | --- |
| PARS | `Pareto-Adaptive Reaction Path Search (PARS) — MLPATH Integration Design.md`（项目根目录；实现覆盖 §43 的 MVP-1..5） | 用户确认（2026-10-08） |
| PRRS 核心 | `mlp_native_perturbation_response_reaction_search.md`（项目根目录，原始方法意图） | 待用户确认 |
| 路径成熟层 `mature/` | `docs/PLAN_GEN_COMPAS_MATURATION_2026-10-03.md`（施工计划） | 用户 2026-10-11 指示施工 M2–M4（探索批序贯停止）；全文仍待确认 |
| connect | `docs/PLAN_2026-09-28.md` | 待用户确认 |
| 评估协议（量尺，先冻结） | `docs/EVALUATION_PROTOCOL_DESIGN_2026-10-13.md`（按其 v1 范围） | 用户审定（2026-10-13）；执行参数待零号批后冻结 |
| 发现层重设计（方案 A 约束松弛） | `docs/DISCOVERY_LAYER_DESIGN_2026-10-13.md`（按其 §1.1 v1 范围）；施工计划 `docs/PLAN_EXPLORE_V1_2026-10-13.md`，代码在独立包 `explore/` | 用户审定（2026-10-13）；M1–M2 已施工，M3 起待逐步计划 |

## 环境

| 环境 | 路径 | 用途 | 生产环境？ |
| --- | --- | --- | --- |
| `openmm_dev`（mamba） | `/home/ruigengji/miniforge3/envs/openmm_dev/bin/python` | PRRS / PARS / mature 的科学计算与测试（torch、MACE、OpenMM、ASE） | **是**（用户 2026-10-08 确认）。不 pip/conda 安装、卸载、升级，不拷贝 |
| `.venv-tools` | `/home/ruigengji/MLPATH/.venv-tools` | ruff 格式化 / lint（两套 Makefile 引用） | 否（用户 2026-10-08）。仍只为 ruff 服务，不往里装科学计算的包 |
| `.venv` | 项目根目录 | Windows 遗留，Linux 上不可用 | — |

183（kasuga02）与本机共用 NFS，上述路径两边相同，不要往 183 拷贝环境。

**跨主机读结果要在写结果的那台机器上读。** 181 和 183 的 NFS 都挂成 `acregmin=acdirmin=21600`（属性缓存 6 小时）+ `nocto`：
183 写的 log 和结果，181 上可能几个小时都看到旧长度、旧目录。查远程作业进度、读远程结果，一律
`ssh korakuen@192.168.0.183 '...'` 在 183 上读。2026-10-08 我因此把一批正常跑完的作业误判成"无声死亡"。

## 结果与脚本位置

- PARS / 集成实验：脚本放 `benchmarks/<主题>/`，结果放 `benchmarks/<主题>/results/<YYYY-MM-DD>/`。
- `runs/` 只放 PRRS run。
- 需要 `git` 的脚本（如 `benchmarks/pars_integrity/rereview.py` 里的 `git rev-parse`）在 NFS 工作区上跑不起来，这是预期的情况，不要绕。
