# mature — 路径成熟层(PRRS 通道的射击诊断与精化交接)

一个**独立模块**,与 [PRRS](../README.md)、[PARS](../pars/README.md) 并排。
PRRS 发现反应通道;PARS 组织、选择并精化代表路径;mature 补上选择与精化之间的
那一层:从**已发现的路径**取候选几何 → 固定几何、独立热动量的诊断射击 →
首次进入式判归 → 过渡区域诊断 → 可追溯的精化交接清单。

```text
PRRS 发现;  PARS 选择;  mature 诊断(射击);  NEB/dimer 精化;  驻点/QM 验证。
```

设计依据:[施工计划](../docs/PLAN_GEN_COMPAS_MATURATION_2026-10-03.md)
(§1 交付链、§3 文件契约、§4 M0–M6、§7 验收门)与
[可行性评估](../docs/GEN_COMPAS_FEASIBILITY_2026-10-03.md)。

---

## 模块契约

与 `pars/` 同级同形态:自带 `pyproject.toml`、`Makefile`、`LICENSE`、`tests/`,
整个文件夹可摘出。与 `pars/` 的唯一形态差别是依赖方向:

1. **`mature` 允许且必须单向 `import prrs`。** 动力学与判归复用 PRRS 现有实现
   (`runner.shoot`、只读 `Registry.identify`),**不复制**;可行性评估明文要求
   "沿用现有 shoot,不要复制动力学实现"。解析路径沿袭父仓库现行约定(`../src`
   上路径,`PYTHONPATH=src`),不经包索引;M4 落地时再决定 console script 的正式
   解析方式(未决,见计划 §3.1)。
2. **反向永不成立:** `src/prrs/` 与 `pars/src/pars/` 都不得 `import mature`
   (静态 AST 扫描,`tests/test_boundary.py`);PARS 自己的"禁止 import prrs"
   文件契约边界不受影响,也从本侧再钉一份。
3. **新增数据只写新运行目录,** 不向源 `network.json` admit 射击端点,不改旧 run
   (计划 §1);输出目录契约见计划 §3.2。

## 状态(2026-10-11)

| 项 | 状态 |
| --- | --- |
| 模块骨架(pyproject / Makefile / README / LICENSE / 边界测试) | ✅ 2026-10-03 |
| M0 清点与基线(preflight) / M1 协议与落盘 | ✅ 见 `docs/DONE_LOG.md` |
| M2 候选提取 / M3 射击诊断 / M4 CLI 与交接 | ✅ 2026-10-11,`src/mature/maturation.py`;只在解析势测试台上验证过 |
| M5 P1 工程实跑 / M6 多通道比较 / Q 学习层 / B 生成桥 | 未开始,入口条件见计划 §5–§6 |

用法:

```bash
mature prepare --protocol protocol.json --out BATCH [--dry-run]   # M2:候选 + 任务表
mature shoot   --out BATCH [--max-new-shots N] [--dry-run]       # M3:逐发射击,可续跑
mature report  --out BATCH                                       # M4:诊断、预算、精化交接
```

- 候选几何取自源 PRRS `response.traj` 的自由段(`free_start`/`free`),按弧长位置;
  PARS selection 只决定取哪几条路径(`refinement_order`,反应管优先)。插值候选标
  `path_interpolation` 并记录来源帧与权重。
- 判归:首次进入侧为 A/B;入口与淬火落点不一致只计数,不改判;未进入时淬火落到
  第三盆记 `other`,否则 `unabsorbed`;门拒与失败各自计数,都不进 `p_B_given_AB` 的分母。
- 探索批可在协议里设 `shooting.interim_look`:前 n 发的 Wilson 区间整段在诊断带外
  就停该候选;排除项已经多到全样本必然超 `other_fraction_max` 时也停。确认批、抽查批固定发数。
- 预算在每发之前检查;力评估与墙钟按已见单发最大值预留。
- 显式候选清单(`candidates.proposal_ids`):只准备列出的候选,列了却产不出的 id
  在写盘前拒绝(零号批固定 1 个候选;id 用 `prepare --dry-run` 看)。再加
  `candidates.from_batch`:候选从那一批原样拷几何(确认批、抽查批打探索批诊断过的同一几何),
  候选记录 `origin_batch`;区域键图不同、几何文件被改过都拒绝。怎么选由
  `shooting.selection_rule` 写明,选了哪些由清单写明。

## 开发

```bash
cd mature && make check   # fmt-check + lint + test(openmm_dev 科学环境)
```

提取清单:拷走文件夹、改 `Makefile` 里的 `RUFF ?=`(默认指向父仓库的
`.venv-tools`)、保证 `prrs` 可导入(见模块契约第 1 条),即可独立运行。
