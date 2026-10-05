# 轮次记录 · 2026-09-04

从 `HANDOFF.md` §10–§11 移出的正文，一字未改。交接主文档只留两行摘要。

## 10. 2026-09-04 工程化重构（`dba26be3` → `23179dea`，零行为改动）

改的是代码长什么样，不是代码做什么。**验收标准是两条同时绿**：
`pytest 224 passed, 1 skipped`（重构前 223，多的一条是新加的 `RATIONALE` 校验）
与 `digest_regression.py` 的 `332 files, 2136 frames, 4 event pairs` **零差异**。
每一步落盘前后都跑过这两条。

改动的机械性由 AST 层面证明，不是靠人读 diff：

| 改动 | 怎么证明它是零行为的 |
| --- | --- |
| `config.py` 的 30 段字段说明搬到模块底部的 `RATIONALE` 表 | 去掉注释后 token 流逐个相同（1883 个）；注释里的每个词都还在文件里 |
| 26 段 `"meaning"` 散文提成模块级常量 | 把常量代回去之后 **AST 完全相同** |
| 39 段 ≥5 行的函数体内注释搬进 docstring 的 `Notes` 脚注 | 抹掉 docstring 之后 AST 完全相同；注释词零丢失 |
| `ruff format`（line-length 96）+ `ruff check` 清零 | pytest + digest |

**注意 digest 绿证明不了什么。** 它比对的是 `runs/` 里已有帧重算出来的 digest，
所以它覆盖身份层（`chemical_key` / `event_key`）和几何，**不覆盖搜索编排**。
上面四项改动全部在 AST 层面被证明等价，这才是它们的依据；digest 只是第二道。

另有一条**只对 digest 沉默的分支**：`follow_unstable_mode` 在
「quench 未收敛 + 力容差满足 + polish 收缩」三条同时成立时才接纳端点，
而 digest 语料里没有这种几何（该现象只在苯甲酰丙酮的近自由甲基转子上出现过）。
digest 对它是**沉默**，不是**通过**；覆盖它的是 `tests/test_convergence.py` 的两条 monkeypatch 测试。

### 全树 reformat 之后，「测试通过」不足以说明测试还有牙

一个被排版打断的 assert 也会安静地通过。所以这一轮的验收还包含一次**变异检验**
（由 `mlpath-e5` 独立执行）：把 `internal_mode_recheck` 的 `"blocked": True` 翻回 `False`，
把 `follow_unstable_mode` 的三条件接纳判据换成 `if True:` —— 两条安全性测试**如期失败**，
恢复后 65 passed。**以后任何一次全量格式化都应照做**：绿本身不是证据，
能被杀死的绿才是。判据选那些 fail-closed 的分支，因为它们正是排版最容易变成空壳的地方。

### 具体做了什么

1. **`config.py`**：dataclass 字段体 215 行 → 101 行，字段列表恢复成可以一眼扫完的表。
   30 段"这个默认值为什么是这个值"的说明**一字未删**，移到文件末尾的 `RATIONALE` 字典。
   `tests/test_validation.py::test_every_rationale_entry_names_a_live_field` 钉住
   每个 key 都是真字段 —— 搬出函数体的代价是它会悄悄过期，这条测试是补偿。
   顺带 `total_charge: int = None` / `multiplicity: int = None` → `int | None`（原来是笔误）。

2. **散文常量化**：`runner.py` 18 处、`search.py` 6 处、`network.py` / `chemistry.py` 各 1 处。
   这些字符串**是输出契约的一部分**（`network.json` 里的 `meaning` 字段），
   提成常量之后改它们是一次 schema 变更，不再是一次注释编辑。

3. **注释搬进 docstring 脚注**：`runner.py` 17 段、`search.py` 12 段、`perturbations.py` 5 段、
   `chemistry.py` 3 段、`network.py` / `openmm_backend.py` 各 1 段。
   原地留 `# [n] <短标签>`，正文在函数 docstring 的 `Notes` 段按序号排列。
   `run_trial` / `Search.run` / `departed` / `_bond_directions` 原本没有 docstring，现在有了。

4. **`engine.py` 移出包**（→ `docs/design/engine_contract.py`）。
   它 215 行，全仓**零 import**，却又定义了第二个 `GateRejected`
   （和 `reliability.py` 里真正在用的那个结构不同）。
   包里放一个没人用的抽象，读起来像是"runner 通过契约和后端对话"，
   而事实是 `runner.py` 直接 import `ase.md.verlet` / `ase.optimize.FIRE`。
   文件原样保留在 `docs/design/`，加了一段说明它从未接线以及为什么。
   **`reliability.GateRejected` 一行未动** —— fail-closed 语义全靠它。

5. **`ruff` 落地**，装在 `.venv-tools/`，**刻意不装进 `openmm_dev`**。
   那个环境的版本号是每一个结果的一部分（backend parity 在一次 openmm 小版本上翻过红），
   所以任何纯开发便利都不许进它。配置在 `pyproject.toml` 的 `[tool.ruff]`。
   `ruff check` 报的 20 条已清零：13 个未使用 import、2 个死变量、
   3 处风格（`i,j,k,l` 的 `l` 是二面角自己的记法，用 `# noqa` 保留并写明理由）。
   `src/` + `tests/` 全量格式化，**代价是 `src/` + `tests/` 一共长了 1606 行（+13%）** ——
   black 系格式化器不保留紧凑多行集合，所以每个 `INTEGER_POSITIVE` 这样的元组
   都从打包排列变成一行一个元素。这是引入 formatter 换来一致性所付的价。

6. **`Makefile`**：`make check` = lint + test + digest。
   `make hash` 现算 `implementation_sha256`，算法与 `search.py` 一致但**不 import 它**
   （import 会让这个脚本和一个坏掉的实现达成一致）。两者若不一致，manifest 是对的。

### 没做，但应该做

- **`Search.run` 仍然是 837 行 / 14 个嵌套闭包 / 最深 8 层嵌套**，占 `search.py` 的 76%。
  `record_channel` / `record_reaction` / `register_saddle` / `publish` 都带着 10–25 行
  解释化学语义的 docstring —— 能写出那种 docstring 的东西应该是一等公民，
  现在它们捕获六个可变变量，无法单独测试（tests 里 `from prrs.search import` 只有 6 次，
  对比 `from prrs.runner import` 54 次）。
  **没做的原因是没有安全网**：digest 不覆盖编排层，pytest 在这一层也薄。
  正确顺序是先补 search 层的测试（或先进 git），再拆。
- **`network.json` 的 schema 仍然是隐式的**。它写着 `schema_version: 2`，
  但没有任何地方定义 schema 2 是什么 —— 它是 `search.py` 里 16 处 `network[...]`
  赋值涌现的结果，而 8 个分析脚本各自手工解析它。
- **全包零 `logging`**，一个五小时的 run 中途完全静默。
- **225 个函数里只有 3 个有参数类型注解。**

---

## 11. 2026-09-04 QM 参照：冻结结论第一次被外部尺子量过

**没有任何冻结判据被改动。** 这一节加的是一把尺子，以及用它做的一次归因。

### 11.1 尺子

ORCA 6.1.1 on `kasuga02`（`scripts/orca_remote.sh`，40 物理核 + NUMA 绑定）。两段协议：

| | 方法 |
| --- | --- |
| 几何 + 频率 | r2SCAN-3c / TightSCF |
| 能量 | ωB97M-V / def2-TZVPD / DEFGRID2，**在上一段的几何上** |

**校准**：气相 SN2 上复现文献到 0.4–0.9 kcal/mol。
**r2SCAN-3c 单独用不行** —— SN2 势垒偏低 11 kcal、丙二醛偏低 72 meV，
meta-GGA 自相互作用误差。几何好，能量必须补单点。

两个鞍点的驻点性都是**测量**的（收敛 + 恰好一个虚模），
不是从对称性论证的 —— 这正是 `runs/sn2_probe.json` 缺的那个前提。

### 11.2 归因：285 meV 是模型的，不是方法的

| | 值 |
| --- | --- |
| `A_θ[E_OFF24]` —— PRRS 在自己找到的几何上（P1 冻结，13/13） | **419.1 meV** |
| `E_OFF24` 在 DFT 几何上 | **410.4 meV** |
| `T[E_QM]` | **134.4 meV** |

**两张势能面的鞍点位置几乎重合**：PRRS 在 MACE-OFF24 上找到的四个 TS 对 DFT TS 的
对称感知 RMSD 全部是 **0.0151 Å**（自同构群 2 元，即两个氧的对换）；
内坐标上 r(O···O) 差 0.042 Å、转移中的 r(O–H) 差 0.020 Å。

所以那 285 meV 的分歧是**纯能量的，不是几何的**。

**注意 419.1 − 410.4 = 8.7 meV 不是"搜索误差"。** 鞍点在正交的 3N−1 维里是极小、
只沿反应坐标是极大，所以在一个偏离几何上取值可高可低；410.4 < 419.1 只说明
DFT 几何相对 OFF24 的鞍点主要偏在反应坐标方向上。
**"PRRS 找到的确实是 OFF24 自己的鞍点"这一条，证据是解析 Hessian**
（saddle order 1 + 力容差收敛），不是这个能量差。

（这个更窄的措辞由 `mlpath-e5` 2026-09-04 提出并纠正；原稿把 8.7 meV 读成了搜索误差界。）
虚频同理：−3186 cm⁻¹ 对 DFT 的 −909，是势垒顶太陡的另一面。

**这是"方法 / 模型"这条分界线第一次真正产出归因，而方向有利于方法。**

### 11.3 三个模型，两个化学，同一组几何

| | 丙二醛质子转移 | 离子 SN2 中心势垒 |
| --- | --- | --- |
| MACE-OFF24_medium | **+276 meV** | 做不了（对电荷逐位无响应） |
| MACE-omol-0-4M | −1.8 meV | **+10.46 kcal/mol** |
| **MACE-POLAR-1-M** | **+12.1 meV** | **+2.09 kcal/mol** |

> **两个数不要混。** 上表的 `+276 meV` 是**固定几何**下的模型误差（同一模型在 DFT 几何上的单点 410.4 对 QM 的 134.4）。
> **P1 的势垒是 419 meV，相对 QM 的偏差约 285 meV** —— 那是 PRRS 在**自己找到的几何**上的值。
> 涉及 P1 势垒时用后者。


**用户 2026-09-04 定：从 P4 起统一用 `MACE-POLAR-1-M`。**
P0–P3 冻结在 MACE-OFF24 上，**不重跑、不改标签**。

**一个被推翻的先验**：MACE-OFF 专门训练在中性有机分子上，
却在中性有机分子的质子转移上比两个 82/83 元素的通用模型差一个数量级。
"专用模型在自己域内更好"在这里不成立。

### 11.4 不能说的

1. **n = 2 个反应。** 两次测量，不是基准。
2. **这不是 `T[E_m]`。** 三个模型只在 DFT 几何上做单点，
   没有一个被允许优化出自己的鞍点。要 ORCA `ExtOpt` 驱动 MLP，~~**尚未做**~~
   —— **2026-09-04 当天稍后就做了**（`docs/experiments/p1_qm/orca/p1_T_Em/`），
   b04 在 POLAR-1 上的那份补于 2026-09-22（`b04_qm/T_Em/`）。本节这一句写就时为真，
   后来被同日的工作超过；引用时以那两个目录为准。
3. **P0–P3 的方法验收在另一个模型上**，所以它们不自动成为 P4+ 的参照。
   要继承需跑 conformance 第 2 阶（冻结配置重跑 b00–b02，只比身份字段）。
4. 文献值（SN2 −10.5/+3.0 kcal；丙二醛 3–4 kcal 电子势垒）**本项目未核验**，
   只用来让比较可证伪。

### 11.5 落盘位置

- `docs/experiments/sn2_qm/`、`docs/experiments/p1_qm/` —— 两份完整记录 + 几何
- `benchmarks/REFERENCES.md` 末节 —— 套件级的协议描述
- `benchmarks/MANIFEST.md` —— 默认模型决定与它的代价
- `scripts/orca_remote.sh` —— 可复现的启动器

---
