# 已完成条目 · 明细

从 `DECISIONS_PENDING.md` 移出的正文，一字未改。那份文档只留标题清单。
接手不必读这里 —— 只在要查某条**当时是怎么验的**时候进来。

## 已完成：第 2 项（`docs/REVIEW_2026-09-05.md` 的三项正确性问题）

`69d53501` → **`5b69dfed41c4bdb2`** · `233 passed, 1 skipped` · digest 绿。三项都有变异检验。

| | 修了什么 | 检验 |
| --- | --- | --- |
| ① | `_probe_minimum` 三态 `complete`/`confirmed`/`verdict`：被 gate 拒绝的方向不再算通过，"没测出来"不再记成"不是极小" | 改回两态 → 失败 ✅ |
| ② | 一个 `withhold()`，四处调用；两侧下降与边界延续不再直接丢几何 | 契约测试钉住「只有一处 append」「三种 source 齐全」 |
| ③ | `failure.json` 最先写，`publish()` 尽力而为；异常来自 `publish()` 时 manifest 不再写不下去（T=300 三个 seed 就是这么死的） | 改回「先 publish」→ 失败 ✅ |

## 已完成：第 3 项（扰动交付）· 2026-09-09

`5b69dfed41c4bdb2` → **`d7d53219baef17e2`** · `236 passed, 1 skipped` · digest 零差异。

4c 的三选一**选了第一条：改位移实现**。`collateral_bond_A` 一字未动 —— 那个阈值不是问题所在。

| | |
| --- | --- |
| 改了什么 | `internal.translate_fragments` + `perturbations._bond_sides`：bond 位移在能刚体切分时按图切成两侧，沿键轴按质量反比平移。**torsion 的 `rotate_fragment` 终于有了 bond 的对应物** |
| 适用范围 | 桥键、跨碎片接近（各端自成一侧）。环上的键两端仍连通 → 无切分 → 落回 continuation；旁观碎片不动 |
| 实测 | b04 那 9 次 collateral 拒绝里 **8 次消失**（`max_bond_change_A` 精确 0），目标误差 ~1e-16 |
| 仍拒绝的 4 次 | 同碎片 1–3 对（`compress [1,3]` / `[1,4]`，Cl–C–H 两跳）。**这不是同一个缺陷** —— 压缩 1–3 距离在刚体约束下无解，能做这件事的坐标是弯曲，属 `bend` 家族。按 README 规矩 4 记为被门拒绝并附理由 |
| 变异检验 | 一条反事实断言钉住旧实现的损伤（0.054 Å on CH₃Cl / 0.0628 Å on b04 输入，都过 0.05 Å 的门）；把交付改回 `displace` 它立刻红 |
| 未做 | 通用复合驱动（第二条路）—— 仍未被需要 |

细节与逐条数字：`benchmarks/b04_sn2_chloride/PROBE_2026-09-04.md` 末两节。

## 已完成：REVIEW §7.3（bracket 高端派生字段不刷新）· 2026-09-09

`d7d53219baef17e2` → **`8418b82b1f6197d6`** · `237 passed, 1 skipped`（在 183 上跑，
`PRRS_MACE_MODEL` 指向 NFS 上的 checkpoint，backend parity 那 6 条不再 skip）· digest 零差异。

细化把高端换掉之后，只有 `high_kind` 被刷新，`target_microstate` 留着 coarse 高端的值。
合成体系上直接复现出来（`SharpWell`，coarse 对是 0.3 返回 / 0.7 到 `c0001/m0000`，
细化到 0.6 停在脊上）：

```
旧：bracket high=0.6  high_kind=boundary_unattributed  target_microstate=c0001/m0000
新：bracket high=0.6  high_kind=boundary_unattributed  target_microstate=None
```

**旧记录说的是"在 0.6 这一档到达了 c0001/m0000"，而那是 0.7 的测量结果。**

修法是一个 `_bracket_high(record)`，初始化与刷新共用它，循环里那条只改 `high_kind` 的
分支删掉 —— 派生字段只有一处定义，不可能再半刷新。
**这一项不改搜索决策**：`bracket` 只进 `scan["local_brackets"]`，`src/prrs/` 里没有任何
地方把它读回去（已 grep 确认），控制流用的是循环里的 `low/high/high_record`。

变异检验：把刷新改回只刷 `high_kind` → `assert 'c0001/m0000' == None` 报出。
断言写成不变式（高端字段必须与 `samples` 里同幅度那条一致）而不是抄这一次的数，
另加一条"必须真有 bracket 移动过高端，否则这条测试什么也没测"。

---

## 已完成：第 4 项（bracket）· 2026-09-09 · `runs/b04_probe3_rigidbond`

在 **183（kasuga02）** 跑，config 与 probe2 **逐字段相同**（脚本核对，差异集合为空），
唯一变量是实现树 `8418b82b1f6197d6`。

```
拒绝 9 → 0      brackets 0 → 1      reactive 1 → 2      成本 36.1 → 39.5 s/trial
bracket  stretch [0,1]  [0.75, 0.80]  high_kind=confirmed  reached=c0000/m0001
```

**结论：现有扫描策略够用，一个字都不用改。** 0.6 一旦有效返回，(0.6, 1.0) 自己成对、
自己细化到 0.05 Å 窗口 —— 与 4b 那条更正的预测一致。**幅度阶梯不加档、上限不扩展、
`refinement_steps` 不动。** §3.1 的重测要求也过了：1.0 Å 仍然反应，0.8 Å 是新增的第二个反应样本。

**所以"缺测区间"这一项作废了，不是被绕过的**：probe3 的 20 个 `stretch [0,1]` 样本里
没有 rejected、没有 unknown，四种状态里只出现了返回与离开两种。
若将来再出现缺测区间，规矩不变（返回/离开/被阻断/未知分开记）。

顺带两条，都在记录里摆着、不是找出来的：

1. **细化给出了更好的鞍点种子。** 种子从 1.0 Å 那一帧换成 0.8 Å 那一帧，
   切向重叠 **0.699 → 0.986**、steps 17 → 14，同一个鞍点（−535.65 → −535.82 cm⁻¹）。
   离开盆地的幅度越窄，帧越贴着脊。此前没有观测到过。
2. **位移不再自己制造碎片**：probe2 有 4 次交付把分子撕成三块（全在被拒的 9 次里，
   所以从未进入网络），probe3 的 47 次交付全是 (2→2)。
   此前"这类化学会碎裂"的读法是错的 —— 是位移实现自己撕的。

逐条数字：`benchmarks/b04_sn2_chloride/PROBE_2026-09-04.md` 末节。

## 已完成：REVIEW §7.4（`free_bonds` 并集误合并）· 2026-09-09

`8418b82b1f6197d6` → **`093bf311a094b3dd`** · `238 passed, 1 skipped` · digest 零差异。

评审把它列为"待验证的风险"。**验证结果：机制是真的，实例是零的。** 两半分开答：

### 机制：会误合并，而且不是勉强越线

乙醇 anti 与 gauche（同一物质、能量差 2 meV、`chemical_key` 相同）：

```
plain RMSD                     0.4383 Å     ← basin_rmsd_A = 0.15 的 2.9 倍，两个真构象
把 C–O 扭转商掉之后             0.0614 Å     ← 判为同一微态
只要有一侧报了这根扭转 free      → existing_microstate（实测）
```

`free_aligned_symmetric_rmsd` 的 docstring 自己写着这件事不许做：
*"Passing a stiff or merely soft coordinate would wrongly merge two real conformers,
which is why nothing infers this from the graph alone."*
**而并集就是一次推断** —— 从"一侧测到平"推到"两侧都平"。

### 实例：语料里从未发生，而且是在当前代码下从未发生

扫了盘上全部 128 个带 quench 报告的微态结构（104 个 (run, node) 组、474 次扭转测量）：

```
tier 计数        soft 473 · free 1
同一节点的微态之间改过档的扭转    0 根
```

那唯一一次 `free`（`runs/p2_round3/chemical/c0001_m0000.extxyz`，`[4,3,8,10]`）
`|梯度| = 0.014 eV/rad`，而 free 的界是 `0.01 × 0.0175 = 1.75e-4` ——
**它是 `flat_biased` 那一档补上之前的产物**（§8.18.1 里冻住 9.97 meV 的那个甲基转子）。
按当前代码它会被判 `flat_biased`，不是 `free`。

**所以当前代码下语料里 `free` 是 0 次，`free_bonds` 恒为空，这条并集分支一次都没执行过。**
（probe3 的两个微态也都是 `[]`。）这不是"没查出问题"，是**这个分支从未运行**。

### 处置：把假设改成 fail-closed，而不是留着等它第一次执行

```
并集照旧，但不与测量对抗：
  另一侧测过同一根扭转、且判的不是 free  →  不商掉（两个几何的测量不一致）
  没人测过                              →  不反对（保住原意：真自由转子不许生出伪构象）
```

新增 `measured_bonds_of`（"测过哪些"，对 `free_bonds_of` 的"哪些是平的"）——
与身份层那次"未知不等于空"是同一个区分，所以注释直接引它。

**对已记录结果零影响**（digest 零差异、语料里没有分歧对），
所以这是预防不是修复，写清楚免得后来人当成一次真实误合并的修复来读。

变异检验：改回朴素并集 → `existing_microstate` 报出。测试两半都钉：
分歧时不合并、无人测量时照旧合并。

**顺带一条值得记的观测**：`free_aligned_symmetric_rmsd` 这套机器
（连同 `free_alignment_samples = 24` 这个旋钮）**在任何已记录的运行里都没有真正商掉过东西**。
它的理由（近自由转子会生出无穷多伪构象）是论证，不是这份语料上的实测 ——
语料里每一根转子要么 `soft`，要么曲率在地板下但**有偏**（`flat_biased`）。
这不构成删掉它的理由，但**它是未被实测覆盖的一块**，比"待验证的风险"这个说法更准确。

---

## 已完成：第 1 项（digest 修正与回归网补齐）

树 `cf29f66b` → **`686a5c4ffb50c9d5`**，`227 passed, 1 skipped`，digest 零差异。

| | |
| --- | --- |
| **分类调用修正** | `classify` 字段**退休**：不重算、不删除，留在 v1 基线里作为缺陷记录，比对时跳过（`RETIRED_EVENT_FIELDS`）。新增正确计算的 `transition_class`，由 v2 守 |
| **v1 写保护** | `--write` 遇到已存在的 v1 基线**拒绝覆盖**并说明理由 |
| **v2 必须存在** | 缺失时**硬退出**，不再静默跳过。理由写在错误信息里：v1 只记 parity 的**个数** |
| **完整立体内容** | v2 记 `locked_bond_parity` / `tetrahedral_parity` 的**内容**；新分类 `parity_content` 与 `parity_only`（个数）分开，`--allow parity_only` 从此只放行个数变化 |
| **投影比对** | v1 只比它自己有的字段，所以加字段不会把冻结基线弄红 |
| **身份关系断言 ×4** | 不查基线，能活过基线重写 |

### 两项针对实现的变异检验

| 变异 | 结果 |
| --- | --- |
| `_locked_bond_parities` 恒返回 `[]` | **抓住** —— digest 2009 个 `parity_content`，pytest 1 failed |
| `record_channel` 不发布自环通道 | **初次存活**（225 passed 全绿）→ 补断言后**抓住** |

第二项证实了一条真实缺口：既有守卫全在高一层 —— `classify_transition` 测分类器、
"P1 两端同 key" 测状态身份，**没有一条说事件被发布过**。
而退化反应恰恰是状态身份扛不动证据的情形：两端同一物质，通道是唯一的记录。
新断言 `test_a_degenerate_reaction_is_published_as_a_self_loop_channel` 用
`SharedProton` 合成体系端到端跑一次搜索。

### 附带修掉的一个真实缺陷：负 seed

补上面那条断言时撞出来的，与 VALENCE / collateral / IRC 都无关。

```
seek_saddle → register_saddle(seed = cfg.seed + saddle_searches)      小整数，≤ 50
            → descend_saddle(seed = seed + 7919 * sign)               sign=-1 必为负
            → confirm_minimum(seed + hop) → _probe_minimum
            → np.random.default_rng(负数)  ValueError
```

**为什么没打到任何已记录的结果**：`confirm_minimum` 在读 seed 之前就返回了，
除非 `minimum_check == "probe"`；而 p1_round8 / p2_round7 / e2e_ethanol 全用 `"hessian"`。
但 `"probe"` 是文档承认的配置（cheap fallback），在默认 `cfg.seed = 42` 下
`sign = -1` 侧**必然**为负（43 − 7919 = −7876）。

修法是 `runner.side_seed(seed, sign)`，用 `SeedSequence` 而不是 `abs()`：
**abs() 在 seed = 0 处碰撞**，两侧都映到 7919，两条下降会共享同一个流 ——
那正是这个偏移要防的事。回归测试
`test_the_two_sides_of_a_saddle_get_distinct_non_negative_streams`。

---

## 已完成：第 2 项（拆开准入权、立体解析、诊断）

树 `686a5c4f` → **`5d905aaa84236c65`**，`228 passed, 1 skipped`，digest 零差异。
**E/Z 功能原样保留**，身份定义未改。

### 2a 准入权与解析状态分离

`Registry.admit` 的那条分支现在**按名声明它属于哪一层**：`refused_by =
"model_domain_declaration"`、`declaration = "closed_shell_only"`。
docstring 写明三者的分工（身份内容 / 身份解析状态 / 模型适用域证据），
以及那条禁止的反向读法 —— **未解析不是域外的证据，已解析也不是域内的证据**；
域属于 `model.json` 与 preflight（checkpoint 元素覆盖、实测电荷响应、可用时的不确定度）。
`_ADMIT_MEANING` 改写成同时说出两句话：注释器提供证据，模型的声明才使它成为拒绝。

### 2b 预算类诊断不再是身份 —— 实测的缺陷

```
bond_order_limit = 1000  ->  key 8710064d12ce6a5a
bond_order_limit =  100  ->  key 38584b6ee5c4eca2
```

同一分子、同一几何、同一代码。**化学身份成了算力的函数。**
事件键有同一缺陷的安静版本：`automorphisms` 放弃时只返回恒等置换，
键差在不完整的群上正则化，键就变了，而哈希里没有任何字段记录这件事。

修法是把两类理由分开：

| 理由 | 归类 | 处理 |
| --- | --- | --- |
| `unsupported_elements` · `no_assignment_satisfies_valences` · `degree_exceeds_valence` | **分子的属性** | 照旧记进 key 的 `stereo_unresolved` |
| `enumeration_limit` · `candidate_limit` | **运行的属性** | **抛 `IdentityUnavailable`** |

按 §2 的 fail-closed：判不出来就拒绝并说明理由，绝不给一个降级的答案。

**不影响任何已记录结果**：语料 2257 帧全部在默认预算内解析
（只有 `None` 2130 与 `no_assignment_satisfies_valences` 127），
四对事件全部 `symmetry_enumerated=True`。所以这是在**未发生过的情形**上钉住行为。

回归测试 `test_a_budget_overflow_is_refused_rather_than_named`，
同时钉住"充足预算下逐字节不变"与"化学类解析失败仍然只记录不抛出"。

---

## 已完成：第 3 项第一步（未知记为未知，v2 身份迁移）

**这是本轮唯一真正改动身份定义的一步，也是唯一需要迁移审查的一步。**

### 改了什么

`chemical_key` 里 `locked or frozenset()` —— `locked` 为 `None`（注释器判不出哪些键承载构型）
时传入空集，算出 `locked_bond_parity = []`，**而那正是"确实没有 locked 键"的分子给出的值**。

身份本身没错：`stereo_unresolved` 在哈希里，两种情形的 key 早就不同。
**错的是字段的形状** —— 而 v2 记的就是字段内容，用一个混淆两种状态的字段没法审查迁移。

现在三种情形可分辨：

| 情形 | `locked_bond_parity` | `stereo_unresolved` |
| --- | --- | --- |
| 解析成功，确实没有 locked 键 | `[]` | `None` |
| 解析成功，有 locked 键 | `[[[hashA, hashB], ±1]]` | `None` |
| **判不出** | **`None`** | 理由 |

### 迁移审查（验收要求逐条对照）

| 要求 | 结果 |
| --- | --- |
| 文件集合、帧数、图、碎片、事件字段一致 | ✅ `graph_invariant` 0 · `shape` 0 · `event` 0 |
| 每项差异逐帧解释 | ✅ **127 帧，恰好 1 个 (旧 key → 新 key) 组合**，全部 `no_assignment_satisfies_valences` |
| 旧状态如何合并或拆分 | ✅ **双射，无合并无拆分**；新 key 与 6 个未受影响的 key 不撞车 |
| 自由转子仍等价 / 应区分的构型仍区分 / 退化事件仍保留 | ✅ 四条身份关系断言全部通过 |

```
65673127b9550ac7… → 4787e68e5cb1cc71…    127 帧（p3_paths 120 · p3_three_model 7）
```

### 登记方式：`MIGRATIONS`，不是 `--allow`

**旧基线一个字节没动，而且仍在守另外 2009 帧。** 迁移写成脚本里一条枚举记录，
含旧/新 key、前置状态、允许移动的字段与其确切取值、帧数、审查结论。
匹配得**精确**：一帧若移动超出声明范围就不算迁移。

**变异检验**：在一个受迁移影响的帧上额外改动 `fragments` →
没有被吞掉，以 `graph_invariant NOT ALLOWED` 报出。

报告也修了：迁移过的帧**不再被说成 "unchanged"**，
而是 `plus N comparison(s) on reviewed migration <id>`，并注明 v1/v2 各数一次。

### 顺带纠正一条把缺陷写进断言的测试

`test_an_unsupported_valence_is_stated_rather_than_guessed` 原先断言
`locked_bond_parity == []` —— **它的名字说"stated rather than guessed"，
而它断言的恰恰是那个 guess。** 改成 `is None`，并补一条对照
（能读懂的分子、没有 locked 键，仍然报 `[]`），那条对照才让前一句有意义。

---

## 已完成：审计 P0 + P1 的四条 · 2026-09-13

`277 passed, 1 skipped` · ruff 绿 · **50/50 已记录 `event_key` 逐字节复现**（见下）。

| | 修了什么 | 检验 |
| --- | --- | --- |
| **P0** | `reaction_event_key` 原来把两端自同构群 **拼成一个 list** 取 min。两个子群的并集不是群，对并集取 min 连并集自身的元素都不是不变式。改成在 `_generated_group()` 的**生成群闭包**上正则化；闭包用同一个 `automorphism_limit` 兜底，溢出照旧 `raise IdentityUnavailable`（部分闭包＝更大的并集，还是同一个病） | `test_event_key_is_canonical_under_the_GENERATED_group_not_the_union`：K1,3 / P4 四顶点对，`|union|=7`、`|⟨union⟩|=24`；用 m ∈ Aut(P4) 重标号后 union 版给出**两个** form，生成群版给出**一个**。`test_generated_group_refuses_rather_than_closing_partially` 钉住预算溢出 |
| **P1-6** | `internal.rotate_about_vertex` 补 `(0, π)` 目的地守卫，与 `displace()` 逐字一致。原来只有 `norm < 1e-12` 的退化轴检查 —— 178° 过得去，0.35 rad bend 折过 π，arccos 回读 `2π−(θ+δ)`，`target_error ≈ 0.63` 只记录不设门，幅度阶梯仍按**请求**幅度标记。这是全部审计发现里唯一会把错误数字写进命运-幅度记录的路径 | `test_rigid_bend_refuses_to_fold_a_near_linear_angle_through_pi`：178° + 0.35 rad 抛 `ValueError`；反向 −0.35 rad 仍精确交付到 1e-12 |
| **P1-1** | `b04/check.py` clean-ledger 分支原来打印 "this IS a statement about the surface" —— `CRITERIA.md` §5 明文写着"**它仍然不是关于势能面的陈述**"。改成逐字引 §5 | selftest 的 "clean" 例仍报 `searched_no_hit`，措辞改为"基内、逐微态" |
| **P1-2** | `CRITERIA.md` §5 有**三条**可采纳条件，`check.py` 只实现了两条。补第三条（gate refusal）：`network.json` 不逐微态归属 refusal，只有 run 级 `outcome_counts`，所以保守读 —— 有 refusal 就censor 全部微态；**`outcome_counts` 缺失 ⇒ UNKNOWN ⇒ 同样 censor**（"没有记录"不是"没有发生"，与 no-ledger 分支同一口径） | selftest 从 3 例扩到 5 例：新增 "gate refused" 与 "refusals unknown"，两者都必须拿掉 `searched_no_hit` |

**P0 不改动任何既有判定**，这是测的不是推的：对 `runs/` 下 50 个带 `event_key` 的通道重算，
key 逐字节相同，且每一对的 `|union| == |group|`（b04 与 P2 是 6==6，P1 是 1==1）——
审计说的"群相等或嵌套的安全区"现在有数了。**已记录的 run 不需要重新裁决。**

### 规范文本（审计 P1-5）

2026-09-09 裁决当时只落到 `b04/CRITERIA.md`、`check.py`、ledger，没回头改**引用源头**：
`benchmarks/README.md` 规矩 5 与 `TAXONOMY.md` 仍写着"不存在 ⇒ 记为协同，两者都是 pass"。
照旧文起草 b06 就会把这个无效推理形式重新进口。两份都改了，走 `benchmarks/README.md`
文末**修订记录**。现在这条轴是**三个**答案：分步 / 协同 / `unanswered`，
协同只能由 `observed_direct` 的正证据认定。

### 没动的（故意）

- **P1-4**（`network.py:304-315` 水库路径自同构溢出返回 `[identity]`，`info` 存进 `components` 但全仓无下游读取，对称匹配静默退化为朴素 RMSD）——
  改成 raise 会改变已有 run 的行为，需要先决定溢出时是拒绝建节点还是降级并标记。
- **P1-3**（四色互异当立体中心的充分条件；WL 色相等 ⇏ 存在自同构）、**P1-7**（方向去重用 WL 色）、
  **max_depth 不进账本**、以及 P2 的全部条目。

## 已完成：商层置换不变性 fuzzing（审计脑洞 3c）· 2026-09-13

`tests/test_quotient_invariance.py` · `409 passed, 1 skipped, 9 xfailed` · ruff 绿 · 9.8 s。

六个商点原来各有各的重标号测试，是各自建的时候顺手写的。**它们之间的缝就是 P0 那个
并集 bug 活了一辈子的地方**：每个商点单独都有测试，六个商点**共有的那条性质**没有。
所以这个文件测的是性质不是站点 —— 置换原子编号、重算、要求逐位回来。

覆盖：`canonical_labels`（色跟着原子走）、`chemical_key`（含喂给它的每个 component，
防摘要碰撞掩盖）、`automorphisms`（群阶 + 色类大小）、`symmetric_rmsd` 与
`free_aligned_symmetric_rmsd`（水库聚类用的两个数）、`reaction_event_key`、
`_bond_directions`（严格：候选集逐位）与 `_angle/_dihedral_directions`（商不变：类的多重集）。
另加一个纯群论 fuzz：300 组随机子群对，`⟨A∪B⟩` 必须与共轭交换。

**故意不断言的三条**（写进文件头，免得后人当 bug"修"弱了断言）：色去重挑**哪个**代表、
`propose` 选**哪些**方向（seeded shuffle 跑在枚举序上）、水库 leader 聚类的到达序 ——
三条都是构造上就依赖编号的，要不变的是它们聚类**所依据的那个数**，那个查了。

### 第一次跑就抓到一条：`reaction_event_key` 是**编号内**身份，不是结构身份

**与 P0 的并集 bug 无关，是同一个函数里更深的第二条。** `canonical` 把每个原子渲染成
`(Z, 原始原子下标)`，而在自同构群上取 min **只能把原子在自己的轨道内移动** ——
不同轨道的原子各自保留输入给它的下标。实测 C2H6，同一根 C–H 断裂：

```
编号 A: canonical = ((), ((1, 2, 6, 0),))
编号 B: canonical = ((), ((1, 0, 6, 2),))     ← 不同 event key
endpoint_graphs 两边相同（那一半是 WL 色建的，不变）
```

**影响范围精确**：一次 run 里所有结构共享同一套编号，所以**任何已记录的 run 都没被影响**
（50/50 复现那次已经证到这一点），它买的那件事也还在 —— 乙烷六根等价 C–H 仍然合成一个 key，
已钉成 `test_event_key_still_collapses_equivalent_atoms_within_one_numbering`。
坏的是**跨 run / 跨 build 比 `event_key`**：`build.py` 换个原子顺序，同一个反应换个 key。
而 `README.zh-CN.md` 与 `_REACTION_EVENT_KEY_MEANING` 都写得像它与编号无关。

**没有顺手修**，因为修法是个有代价的设计选择，不是 bug fix：要按色诱导的正则序渲染
（先按 `canonical_labels` 排原子、再在群上取 min 消掉色类内的任意打破），**这会改掉
历史上每一个 `event_key` 摘要**，包括 b04 冻结判据引用的那些。轨道 ⊊ 色类时（WL 不完备，
即审计 P1-3/P1-7 同一个洞）仍不正则，也要一并写进限定。记成 `xfail(strict=True)`
而不是删掉，理由写在 reason 里 —— **strict 是关键：哪天真修好了，这条会因为"意外通过"而报错，
逼人回来改文档。**

## 已完成：b04 冻结 + 脑洞 A（谱提案）+ 脑洞 B（碎片刚体转动）· 2026-09-13

`414 passed, 1 skipped, 9 xfailed` · `make lint` 绿。**`src/` 动过，digest 基线仍待重跑。**

### ① b04 冻结（用户 2026-09-13 落笔）

判据条目**一字未改**，只加两处：文首冻结戳、§5 末的实现说明。`MANIFEST` 状态词
`preflight` → `frozen`。实现说明记的是：§5 第三条要求**逐微态**归属 gate refusal，
而 `network.json` 只有 run 级 `outcome_counts`，所以 `check.py` 保守实现 ——
run 内有任何 refusal 就 censor 全部微态，**字段缺失同样 censor**。
**实现比文本更严，不是更松**，只会让可采纳的负记录更少。
**冻结不追认既有 run**：probe1–4 都在冻结之前跑的，是探针不是验收。

### ② 脑洞 A：谱提案（转子投影版）—— 先量后建，量到了

`docs/experiments/spectral_proposal_probe.py`。A 原本**没有失败案例支撑**，只有可信度论证，
按 `HANDOFF.md` §2（旋钮先验证再设计实验）必须先量。量的问题很窄：
**在一个已知答案的 run 上，实际反应的那根方向在谱排序里排第几？**

b04 探针结果（源微态，解析 Hessian，`asymmetry 3.0e-11`，平凡模地板残差 1.87e-05）：

```
1  stretch:0-1:+   overlap 0.9651   run order 1  <== 实际反应的那根
4  compress:0-5:-  overlap 0.6893
7  compress:2-5:-  overlap 0.2558
...
28 kick:0-4:+      overlap 0.0031
```

**28 选 1，第一名，0.965 对第二个不同方向的 0.689。** 而运行里它排第一是**运气**：
Cl⁻+CH₃Cl 没有氢键给体，所有候选都在 band 1，族内顺序纯粹是 seeded shuffle。
**A 把彩票换成机制。**

探针自己先踩了一个坑，值得记：`curvature_spectrum` 把六个平凡模**投影掉而不是删掉**，
它们作为零本征值留在谱的**软端**——"最软的四个模"字面取出来就是平动转动，
与所有内坐标正交。第一版探针因此把 28 根方向全部打成 `0.0000`。
这个坑已钉成 `test_soft_internal_modes_drops_the_rigid_body_block`。

落地：`runner.soft_internal_modes`（投掉平凡 + 扭转子空间，返回**质量加权**单位向量）·
`perturbations.spectral_overlap` + `propose(..., soft_modes=)`（改写族内 rank，
**band 机制不动**——氢键给体照旧凭自己的证据拿前排）· `search` 在提案处按微态**记忆化**计算。
新旋钮 `spectral_rank_modes = 4`，进 `RATIONALE`，理由是实测的谱隙
（0.0186 / 0.0187 / 0.0392 然后跳到 1.5697）。**拿不到谱就回落到原顺序** —— 排序不值得让一次搜索失败。

### ③ 脑洞 B：碎片刚体转动 —— 缺口是真的，**但没解决 b05**

`internal.rotate_rigid`（第三种精确实现：轴过碎片自身质心，不需要轴上有原子）·
`perturbations.fragment_rotations`（多碎片时每个碎片绕自身惯量主轴三根；主轴而非实验室系，
因为实验室系不是结构的性质）· `polish_soft_modes` 统一走 `soft_mode_step` / `soft_mode_response`。
精确性已测：碎片内距离与碎片外距离变化 **0.00e+00**，质心漂移 **0.00e+00**。

**定义域缺口确实存在且已关闭**，但 `docs/experiments/b05_fragment_rotation_polish.py`
在三个卡住的端点上实测：

```
                     打磨前 -> 打磨后        能量        实际转动
t000002   order 1, -76.8 -> order 1, -75.1   -0.43 meV   1 步 0.0193 rad
t000058   order 1, -77.3 -> order 1, -75.6   -0.44 meV   1 步 0.0196 rad
t000080   order 1, -75.4 -> order 1, -73.6   -0.42 meV   1 步 0.0192 rad
```

**鞍点没动。** 原因清楚：`polish_soft_modes` **故意不碰负曲率**——"沿这根坐标不是极小"
是极小检查的裁决，不是打磨的。所以同一个脚本问了第二个问题：刚体转动是**精确且周期**的，
全周期扫描在局部曲率不可用时依然良定义（`flat_biased` 那一档已经是这个论证）：

```
全周期扫描：+29.1 / +28.7 / +28.5 meV，转动 1.571 rad（四分之一圈），三个样本一致
            结果仍然 order 1，-54 cm^-1
```

**所以 b05 的阻塞不是打磨定义域问题。** 找到了一个低 29 meV 的水取向（三个对称等价样本
一致到 0.6 meV，是实测不是巧合），残余的 −54 cm⁻¹ 是**耦合的碎片间运动**，
一次动一根坐标去不掉 —— 而碎片**平动**正是这次故意留在外面的那一半
（单位是 Å，需要自己的曲率地板、自己的容差、以及支撑这两者的自己的测量；
一个三档机器喂两套单位穿同一个阈值是最静默的那类错误）。

**`DECISIONS_PENDING` 里原本的诊断依然成立**：先放开 `quench_steps` 与 `irc_relay_hops`
重跑六条下降。B 不是它的替代品。

**未决的语义选择（不是我能定的）**：打磨要不要对**精确周期**坐标扫过负曲率？
要，就等于改变 `confirm_minimum` 被交到手里的是什么；29 meV 那个数已经量好放在这里了。

### 两个下游消费者的防护

`free_bonds_of` 与 `measured_bonds_of` 都从打磨报告里取 `indices[1:3]` 当键。
碎片转动的 `indices` 是碎片的原子表，`[1:3]` 是两个**未必成键**的原子。
两处都加了 `coordinate == "torsion"` 的门（缺字段默认 torsion，老记录都是扭转）。

## 已完成：digest 基线重跑 + b04 冻结后首次运行（kasuga02）· 2026-09-13

### ① 跨机协议照走，两个坑各踩一次

`env/BASELINE.md` 记着的两条，这次原样又遇到：**`ssh 不 cd`**（落在 `/home/korakuen`，
`src/prrs/*.py` 匹配不到）与 locale 排序差异。核完两机 **17 个源文件 + 2 个新测试逐字节一致**，
NFS 缓存没有滞后，才开始跑。

### ② digest：七个文件改动，零摘要移动

kasuga02 上 `python -B`：

```
332 files, 2136 frames, 4 event pairs —— digests unchanged
  + 254 comparisons on reviewed migration 2026-09-04-unknown-is-not-empty
四条命名检查全过
implementation_sha256 = c393da782a1350a1   (probe4: fa15a829869cdc29)
pytest 两机一致：416 passed, 1 skipped, 9 xfailed
```

改的是 `chemistry / internal / config / network / perturbations / runner / search`
（生成群闭包、`(0,π)` 守卫、谱提案、碎片刚体转动），**一个持久化摘要都没动**。

### ③ b04 的 config 也必须冻，而且缺一个键

冻结 CRITERIA 那一步漏了 `config.json` —— 它文首自己写着"冻结与 CRITERIA.md **同时**发生"。
而且按 `benchmarks/README.md` 它必须是**完整**的 SearchConfig dict，实际 100/101，
缺的正是新加的 `spectral_rank_modes`。**不补进去就等于冻结一份静默省略了
"预算花在哪些方向上"这个旋钮的记录。** 已补（默认值 4），文首写明这是冻结时**新增**
而非修改，理由指向实测。现在 101/101。

### ④ `runs/b04_frozen1`：冻结后第一次运行，`check.py` 退出码 0

```
取代通道 : observed_direct   D1-D4 全 PASS
           ts0000 / ts0001 两个 order-1 鞍点，都支持 rc0000（自环）
           event_key = 847306d3f8c57bd4  与 probe2/3/4 逐字节相同
分步侧   : unanswered   candidate_not_selected, proposed_not_executed
           逐微态 uncovered=14 / not_executed=5，refusals=no
```

`event_key` 一致，是**生成群闭包改动在一次活运行里的确认**。

### ⑤ A 的端到端验收：能读什么，不能读什么

运行自己记着排序来源（每条 scan 的 `proposal.budget.ranked_by`）：
`spectral overlap with 4 soft internal mode(s) (analytic; 0 torsional removed)`。
族内头名 4 个族里 2 个换成了 overlap 显著更高的方向（compress 0.2558→0.6893、
kick 0.2558→0.9651），另 2 个原本就是。**这是 A 在一次运行里唯一能证明的事。**

`ts` 1→2、`reactive` 2→4、`executed` 17→18 **不得归因给 A**：两次运行
`implementation_sha256` 不同，且这一轮同时动了七个文件 —— 这不是受控对照。
要隔离 A 需要一次只切 `spectral_rank_modes` 的 A/B，**没做**。

**没变的**：基覆盖仍 14/28，负侧照旧 `unanswered`。A 不改覆盖率只改顺序，
读成"覆盖变好了"是错的。

`MANIFEST` b04 状态词 `frozen` → **`done`**。

## 已完成：PLAN 七决策落笔 + 第 4 件（schema 显式化）+ P1-4 · 2026-09-20

`c393da782a1350a1` → **`b4e92642c66c9bab`** · `420 passed, 7 skipped, 9 xfailed`（yayoi；
7 个 skip 全是环境性——本机无 MACE checkpoint / 无 CUDA，kasuga02 上为 1）· **digest 零差异**
（332 文件 / 2136 帧 / 4 事件对）。两次变异检验（P1-4 守卫删除、schema 校验器灌假问题）都被
测试杀死后恢复，恢复后指纹即 `b4e92642`。

### 前置：七决策落笔（开写的门）

用户授权按 PLAN 推荐项落笔（出处与七条决定见 `DECISIONS_PENDING.md` 2026-09-20 收口节）。
要点：① `structural_key` 加字段不换 `key`；② 1.3 差异非零不动 `chemical_key`；③ committor
走 (a) Maxwell–Boltzmann T=300 K 独立探针，N=16 分辨率开跑前写死；④ 证书命题草案待批、
批准前不发；⑤ 碎片平动方向基 = 碎片间连线 + 正交补；⑥ B′ 暂不启用（先跑 5.0）；⑦ P1-4
溢出拒建节点。

### 第 4 件：`network.json` schema 显式化

**问题实测**：`schema_version: 2` 自 2026-09-02 起写在每个网络上，但全仓无一处定义 schema 2
是什么。对 `runs/` 70 个 `network.json` 逐个校验（见下）：2 个是 schema 1（`nodes`/`edges`
时代）、3 个 schema-2 文件没有 `reaction_channels`、只有 5 个带方向账本——**版本戳盖在一个
未定义的对象上**。

**做了什么**（`src/prrs/schema.py`，无 jsonschema 无 pydantic，一张字面量表 + 一个校验器）：

| 部件 | 内容 |
| --- | --- |
| `SCHEMA_2` | 按**当前写入器**逐字段描述七个集合 + 微态 summary + 方向账本；`required` 只收写入器**每次都写**的字段，代码里有分支的可选字段（`ts_support`、`connects`、`transition_class`）在 schema 里同样可选——optionality 必须能追到 search.py 里的一处分支，不是乐观 |
| `validate(network)` | 返回问题列表不抛；lifecycle 规则：`completed` ⇒ termination/termination_detail/三张计数表/budgets/两级方向账本，`running`/`partial`/`failed` 快照不被要求这些 |
| `publish()` 末尾挂钩 | 生产路径把问题记进 `network["schema_warnings"]`（干净时不写该键——文件与无此检查时逐字节相同）；`tests/test_schema.py` 与 CLI `--strict-schema` 才把非空当失败。五小时的 run 不因描述滞后一个字段而死 |
| `read_network(path)` | 给新消费者共用；**现有 13 个消费者一个没改**——改它们没有收益 |
| `tests/test_schema.py`（9 条） | **锚测试**：`runs/` 全部 70 个 network.json 过 validate；不在漂移清单里的文件必须零问题（前向门），清单里的文件必须**逐字复现**其记录的偏差模式 |

**漂移清单**（锚测试钉死，索引起成 `[*]` 后去重）：65 个 schema-2 文件的偏差**全部可用写入
年代解释**——老 run 缺的是后来才引入的字段（`withheld_candidates`←5dc83ed4、两级账本与
`termination_detail`←fa15a829/89f41c80、`mechanisms`/`continuations`←d7d53219、
`found_by`/`curvature_source` 等更新），最老的三个连 `reaction_channels` 都没有。
**没有一处未解释的不一致**。当前写入器家族的 5 个 run（`b04_frozen1` + 四个 b05 probe）全绿。

**变异检验**：校验器改成永远返回一个假问题 → 9 条里 7 条红（新鲜 run 测试抓住 publish 挂钩
把假警告写进文件、锚测试的 clean 断言全炸）；校验器改成永远返回空 → 锚测试在 legacy 文件上
红（漂移清单断言不匹配）。两个方向的绿都死得掉。

### P1-4：水库路径自同构溢出拒建节点（决策 7）

`automorphisms()` 溢出返回 `[identity]` + `reason="candidate_limit"`，事件路径有
`BUDGET_REASONS` 把关，但 `network.py` 建节点这条路径没人读 `reason`——群从"未认证"静默
变成 `[identity]`，构象匹配随之退化为朴素 RMSD（甲基转 120° 永远像新构象），与
`chemistry.py` "An incomplete group is never used" 的自我声明直接矛盾。

修法（决策 7）：`admit` 的新节点路径上 `info["reason"] ∈ BUDGET_REASONS` ⇒ 返回
`Admission("rejected", reason="automorphism_budget", detail=info)`。**不 raise**——三个
消费点已有完整的纸面追踪（trial 端点 → `network_exclusion` + withhold；两侧下降与边界延续
→ `admission_reason` + withhold），run 不死，几何不丢，节点不建。源结构若溢出，
`run()` 原有的 `new_chemical_node` 检查把它变成显式 RuntimeError——同样 fail closed。

语料里该分支从未执行（2136 帧最坏 57 search nodes，b05 反应物 3889，都在 20000 内）——
这是预防不是修复，digest 零差异是它的证据。

**测试**（`tests/test_automorphisms.py` 第 18 条）：乙醇 + `automorphism_limit=2` ⇒
`rejected` / `automorphism_budget` / detail 带 `candidate_limit` / `registry.nodes == []`；
同一结构在默认预算下照常 `new_chemical_node`（拒绝是预算的不是一刀切的）。
**变异检验**：删守卫 → 测试红；恢复 → 绿，指纹回到 `b4e92642`。

### 没做的（照 PLAN 顺序排队）

第 1 件正则标号（最大的一块，唯一挡队列的）→ 第 2 件 committor → 5.0 零代码诊断 → 第 5 件
碎片平动 → 第 6 件续跑（先写 6.1）→ 第 3 件证书。随手件剩三：logging 静默、manifest 去
hostname、`max_depth` 记账词。

## 已完成：第 1 件（正则标号）+ 第 2/3/5/6 件的机器与协议 + 三随手件 · 2026-09-20 下午

`b4e92642` → **`b1c4d44009f8785d`** · `462 passed, 7 skipped（环境性）, 9 xfailed` · **digest 零差异** · lint 绿

### 第 1 件：正则标号取代 WL 色（B / P1-3 / P1-7 一次修三条）

**1.1 引擎**（`src/prrs/canonical.py`，个体化-细化，不写 nauty 剪枝）。实现由 GLM-5.3-Flash
子代理按 PLAN §1.1 写出（工作流 dwfrun-dd1e3eef 起、dwfrun-5a1b17dc 因 flash 配额中断后由
主会话接管验收——配额 20:05 重置，事件如实记录）；**验收全在主会话**：

| 验收 | 结果 |
| --- | --- |
| 独立 fuzz（编排自写、实现者不可见；修正过一处**门自己的逆置换 bug**——元素的映射约定与边不一致会构造非同构图冤枉引擎） | certificate 重标号不变、无边色 colours 与旧 `canonical_labels` 逐字节同、群元皆保边双射 ✓ |
| 群逐位 oracle（`tests/test_automorphisms.py` ORACLE：b05 1296 / 乙醇 12 / b04 6 的完整 permutation 集 sha256） | 逐位相同 ✓ |
| 性能 | 引擎 IR 树与旧回溯内核同量级（语料最坏仍远低于 20000 预算） |
| 变异检验 | **种子随机的目标类选择数学上惰性**（整类展开则群仍完备——已如实记录）；改用无种子随机 → 6 条测试红（同图两次调用不再逐字节相同）→ 恢复绿。另补两条更锐的确定性测试：canonical labelling 是"图性质至多差一个群元"（轨道隶属断言）、同图可复现 |
| 引擎自行堵的洞 | 仅边集渲染的 certificate 区分不了水与次氟酸——证书带色序列与边色，docstring 有论证 |

**1.2 `structural_key`**（决策 1：加字段不换 key）。事件编码为一张边着色叠合图
（kept/broken/formed），`canonical_form` 出正则形式，方向商掉（broken↔formed 对换取小证书）。
`reaction_event_key` 增 `structural_key` / `structural_form` 两个字段，`key` 一个字节未动——
b04 冻结判据引用的历史摘要全部复现（digest 为证）。channel 落盘增 `structural_key`，
去重仍按 `event_key`（一次 run 内共享编号，语义正确）。
**`tests/test_quotient_invariance.py` 的 `xfail(strict=True)` 在同一次改动里摘掉**，
换成 structural_key 的正面断言 + "旧字段允许随编号变"的限定钉死。方向无关（反参数同 key，
两侧无序对相等）与区分性各有测试。

**1.3 轨道取代色类**（先量后改）。`docs/experiments/orbit_vs_colour_probe.py` 对 digest 语料
2136 帧实测：立体中心 / E/Z 取代基 / 方向去重三个判据，**色与图轨道零差异**——改动免费，
决策 2 不触发，`chemical_key` 一字未动（digest 零差异兑现了预言）。改动：

- P1-3：`_tetrahedral_parities` 色四互异（充分非必要）→ **Stab(中心) 轨道四互异**
  （`_stabiliser_orbits`）；`_locked_bond_parities` 色平手 break → 轨道判等、轨道内按正则序
  取最小。`chemical_key` 一次 `canonical_form` 同时供 colours（逐字节同旧）与认证群；
  群溢出时保持色判据（= 历史行为，不 fail-open 也不拒绝）。
- P1-7：新增 `geometric_automorphisms`（群里近似保持当前几何的子群，一次 Kabsch/元素，
  拒镜像）。**propose 的去重键换几何轨道**（固定构象可破图对称，色/图轨道并掉的端桥键是
  看不见的覆盖损失）；群溢出 → identity 子群 → 每原子自成一轨道 → **宁可多列不猜合并**
  （与事件 key 同一失败方向）。测量侧（`torsion_for` / `torsion_profile` / runner 三处）换
  **图轨道**——坐标集是物质的性质不是某个构象的性质。**转子周期（symmetry_order）刻意留在
  图对称**：几何子群会把甲基 C3 也拆掉，搜索域凭空 ×3 全是片段内部对称的副本（首次实现
  踩中此坑，`test_rotor_amplitudes` 抓住，已回退并记录理由）。

### 第 2 件：committor 机器 + 协议 + 接线（实跑排队 MACE 机）

- `runner.shoot`（同几何 N 发只换 thermal seed；**T=0 拒绝**，测试钉死）+ `wilson_interval`
  （16/16 → 下界 0.8064 与项目引用值逐位一致）。
- `config.committor_shots`（默认 0=关，RATIONALE 写明 N=16 分辨率限定与"只答接近 0/1"）。
- search 接线：已 bracket 窗口从 ridge 帧（复用 `_dividing_frame`）打 N 发，A/B/**其他**
  三分类（其他绝不进分母），`bracket["committor"]` 记录 + 协议引用。射击 trial 目录用
  bracket 前缀防撞。
- 协议冻结于跑前：`docs/P1_COMMITTOR_PROTOCOL.md`（系综/T/N/Wilson/判归/范围守门全先写）。
- e2e 测试（解析双阱，T=300）验证记录形状与算术。

### 第 3 件：完备性证书机器（解析势全链条，命题待批故不签发）

`src/prrs/completeness.py`：`cover`（覆盖半径 ≤ δ 的格点，**含球外 δ 环**——首版把边界点
的最近格点滤掉了，覆盖性测试当场抓住）、`membership`（quench+回盆 RMSD，未收敛/门拒
记自己的结局，绝不记 True）、`gradient_bound`（样本对上确界，声明为真 Lipschitz 常数的
**下界**；近重合点对跳过——浮点噪声会把比值放大成假曲率，实测抓过 124 eV/Å² 的噪声）、
`audit_basin`（三条件：membership 全回盆 / L·δ < g_min / 边界壳能量 > X；四量缺一即
ValueError）、`certificate_text`（**永远盖 DRAFT 章**——决策 4：命题批准前不发证书）。
解析三阱势验收：证书对已知答案成立；**δ 放大到覆盖论证撑不住的值时必须拒绝**（即使每个
采样点都回盆——PLAN 的变异检验，落在 (ii)）；球含第二个盆时被 membership 诚实拒绝。

### 第 5 件：平动的算子与枚举（阈值等测量）

- 5.1 `internal.translate_rigid`（第四种精确刚体实现，验收同 B：碎片内距离 0.00e+00、
  COM 还原、组间相对位移精确等于步长）。
- 5.4 `_tier` 抽出，**常数作参数**（弧度/Å 两套各传各的——脚注 [1] 点名的单位混喂错误
  从结构上排除；单测钉住两套常数同逻辑）。
- 5.2 `fragment_translations`（决策 5 的 COM 连线基 + 确定性正交补）**惰性交付**：只枚举
  不接线——三条阈值必须来自 b05 端点实测（脚本已备：
  `docs/experiments/b05_fragment_translation_measurement.py`，需 MACE），不是从弧度换算。
- 5.0 协议：`docs/experiments/b05_budget_rerun_protocol.md`（quench_steps 1500 /
  irc_relay_hops 8，kasuga02 执行，判读先写）。

### 第 6 件：6.1 + 6.2（resume 的地基）

- `Registry.from_network`（除 `permutations` 外全部从盘重建；permutations 由首结构重枚举——
  枚举确定性使这合法；群溢出时拒绝重建，与 P1-4 同一 fail-closed）。
- 6.1 等价性测试：活库逐字段比对（含重枚举群逐位相同、evicted 名单、轨道）+ 行为等价
  （重接纳 existing、编号续到 m0003 不与历史相撞）+ 真实数据验证（b04_frozen1 重建：
  1 节点 2 微态 |Aut|=6，重接纳源结构 existing_microstate）。
- **6.3–6.5（search.py 的 resume 接线）未做**，原因如实：除 Registry 外还要重建
  records/ts_candidates/saddle_searches/continuations/schedules/queue 的运行态，
  这是 PLAN 里"最难也不能省"的部分，不在长会话尾部赶工。6.1 的纪律是"写不出来就放弃"
  ——地基验证过了，接线另起一次会话做。

### 三随手件

- **logging**：全包原本零日志，五小时 run 全程静默。加 `prrs.search`/`prrs.runner` 两个
  logger（源接纳/扫描/bracket/鞍点/通道/收工/失败 + 下降落定/极小模式收敛），CLI 入口
  `basicConfig`（库导入保持静默），caplog 测试钉住。
- **hostname**：**前提被实测推翻**——src 里没有任何 `platform.node()`/`gethostname()`，
  `versions()` 只写 `platform.platform()`（不含主机名）。按规矩记"测量杀了它"，不硬做。
- **max_depth**：第五种预算截断有了词——`budgets["admissions_capped_by_depth"]`，
  三个队列门槛各计数；不进方向账本（它截的是接纳不是方向）。测试：max_depth=1 的双阱
  run 计数 ≥1，默认协议 =0。

### 配额事件（如实）

GLM-5.3-Flash 的 5 小时配额在引擎写完、验收门未跑时用尽（今日 20:05 重置）。引擎代码
已在盘上且独立 fuzz 预跑通过，剩余验收与全部接线由主会话模型完成。**分工模式的结论**：
flash 写 + 主会话出 spec 和门是可行的——引擎一次通过全部验收；但编排脚本里的独立门
本身也要被检查（这次门里的逆置换 bug 差点冤枉引擎）。


## 已完成：第 6 件收尾（6.3–6.5 resume 接线）· 2026-09-20 晚

`b1c4d440` → **`1a6a297143deb645`** · `464 passed, 7 skipped（环境性）, 9 xfailed` · **digest 零差异** · lint 绿

**验收就是 PLAN 6.3 的原话,唯一接受形式**:子进程在 `run_trial` 内 `os._exit(137)`(真 kill 语义——
异常处理器不跑,network 停在 `status: "running"`,与断电/SIGKILL 同形)→ 主进程 `run(resume=True)` 续完 →
**network.json 与 attempts.jsonl 与不中断孪生逐字节相同**(`tests/test_resume.py`)。

### 三条硬前提(任一不满足即拒)

network.json 存在且 `status == "running"`(异常快照 partial/failed 不可续——它们带着死因);
config 逐字段相同(json 往返的元组→列表已规范化处理);`implementation_sha256` 相同——跨实现续跑是混批,
本项目作废过一轮。另测:对已完成 run、改配置、目录不存在三种拒绝。

### 重建清单(全部从盘,零猜测)

| 状态 | 来源 |
| --- | --- |
| registry | `from_network`(6.2,已验证) |
| records | attempts.jsonl **裁剪到 scan 边界**(见回滚) |
| ts_candidates / reactions / channels / withheld | network.json 原 dict 续用 |
| queue | publish 新落 `network["queue"]`(6.4 的唯一新字段) |
| saddle_searches / continuations | 数 `paths/climb*` / `paths/continuation*`(PathRecorder 即刻建目录,数与计数器同源) |
| schedules | 各 scan 的 `proposal.budget`(首个 propose 的 report)+ `proposal.schedule` 的 offered/executed/halted 逐向恢复 |
| spectra | 惰性重算(确定性:同结构同 config 同结果;位精确性由 sidecar 保证) |
| depth_capped | 推导:cap 深度节点的 microstates+evicted(每次接纳都发生在节点定死的深度上) |

### 回滚(6.5)

在飞 trial = attempts.jsonl 里不被任何已发布 scan/ts/withheld/reactions 引用的行 → 从 jsonl 删除、
trial 目录删除;committor 射击目录按 scan 序号判定。**删除是唯一动已落盘数据的地方,逐条记进
`resume.json`——独立文件,不进 network.json,否则逐字节验收连陈述都不可能**(对 PLAN 6.5 原文
"记进 network['resume_log']"的一处显式偏离,理由即此)。

### 过程里被验收测试抓住并修掉的三个洞

1. **extxyz 文本精度**:续跑从盘读构象,位置丢末位 → 下游谱/Hessian 全漂 1e-9,逐字节验收不可达。
   解法:`write_state` 同时落 **npz sidecar**(positions+forces+energy,float64 位精确),`from_network`
   优先读 sidecar;旧 run 无 sidecar 照旧读文本(resume 不适用于它们)。sidecar 只为 resume 存在。
2. **重建构象无计算器**:run_trial 要 reference energy、谱要 forces——sidecar 里都有,挂
   SinglePointCalculator。
3. config 比较:json 往返把元组变列表,先同序列化再比。

### 状态

第 6 件全件完成。PLAN 六件里**只剩**需要 MACE 机的实跑(kasuga02:5.0 重跑 / 平动阈值测量 /
committor 实跑)与决策 4 的命题批准(证书签发的门)。本机可做的代码与机器,全部落地。


## 已完成：kasuga02 实跑队列 + committor 判归修复 · 2026-09-22

`ec32bc3df5d2bf14` → **`63be7064f9e09f82`** · `471 passed, 1 skipped` · **digest 零差异**
（332 files / 2136 frames / 4 event pairs）· lint 绿 · 变异检验过

PLAN 六件里「只剩 MACE 机实跑」的那三件全部执行完，外加 b04 在当前指纹上的重跑与
b04 的 T[E_m]。串行队列脚本 `.scratch/queue_183.sh`，fail-closed 在指纹上，逐步日志在
`.scratch/queue/STATUS`。

| 步 | 结果 | 落在哪 |
| --- | --- | --- |
| 1 b05 preflight | **7/7 pass**，gate 9 **1149 ms** | 见下「过期文本」 |
| 2 b04 重跑 | `observed_direct` / 分步侧 `unanswered`，退出码 0 | `runs/b04_frozen2` |
| 3 完备性审计 | 首跑作废（我的 bug），修正后 (i)(iii) 过、(ii) 差 11.7 倍 | `runs/b04_completeness_audit` |
| 4 b05 预算重跑 | 六条下降 **0/6 → 4/6 到确认极小** | `runs/b05_saddle_descent_budget` |
| 5 平动阈值测量 | `d²E/dΔ² = 1.737 eV/Å²`、`|dE/dΔ| = 0.0930 eV/Å` | `runs/b05_fragment_translation_measurement` |
| 6 b04 的 T[E_m] | 三次才收敛，结果与 P1 同形 | `docs/experiments/b04_qm/T_Em/` |
| 7 P1 committor | 判归有 bug，杀掉；落点重判后有数 | `runs/p1_committor1/rejudge_2026-09-22.json` |

### b04 重跑：四次指纹移位对编排层输出零影响

`runs/b04_frozen1`（`c393da78`）与 `runs/b04_frozen2`（`ec32bc3d`）逐字段相同 ——
`ts_candidates` 2 个、`outcome_counts {'completed': 64}`、`classification {'elastic': 60,
'reactive': 4}`、`admission {'existing': 63, 'new_microstate': 1}`、`termination
configured_scope_exhausted` 64/96、覆盖 `uncovered=14 / not_executed=5`、
`event_key = 847306d3f8c57bd4`。**唯一的差是多了 `budgets.admissions_capped_by_depth: 0`**
（09-20 那件随手件的新记账字段，值为 0 说明这一格没撞深度上限）。

所以「判定绑不绑指纹」那个口径不必落笔了：**不是"digest 零差异所以推定不变"，是实测不变。**

### b05 预算重跑（5.0 协议，零代码）

`quench_steps 500→1500` / `irc_relay_hops 3→8`，其余逐字节不动。三个 order-1 驻点、六条下降：

```
t000002: 两侧都到确认极小 (relays 1 / 4)
t000058: +1 到，−1 停在 fmax_not_reached_after_polish
t000080: −1 到，+1 停在 fmax_not_reached_after_polish
四条到达的落点全同: C4H8 + Cl + H2O, key a4cc3716d8bb73fe, order 0, dE ≈ −1.322 eV
```

2026-09-09 裁决「先放开预算再谈新机制」押对了。**剩下两条的停止理由是新的**
（`fmax_not_reached_after_polish`，不是旧的两道预算），按协议记为新事实，不沿用旧结论。

### 残余负曲率在转动上 —— 第 5 件（碎片平动）的前提被实测否掉

```
lowest eigenvalue −0.0217 eV/(A²·amu)
  99.4% 落在水碎片上
  水内部：转动 99.9% · 平动 0.1% · 内部形变 0.0%
```

平动阈值同日测出（`1.737 eV/Å²` / `0.0930 eV/Å`，三个端点一致），**但它量的自由度只占
残余模的 0.1%**。「平动可能是 −54 cm⁻¹ 的出路」这个假设到此为止。

**决策 6 的触发条件因此成立**：5.0 跑完，残余负曲率仍在精确周期坐标上（转动，`rotate_rigid`
正是精确周期实现）。

### b04 的 T[E_m]：ORCA 的优化器在 POLAR-1 面上自己找鞍点

P1 那条线建在 MACE-OFF24 上，而 P4+ 的默认模型是 POLAR-1 —— 换了模型，
`A_θ[E_m] vs T[E_m]` 这一行要重新建。三次才成，每次的失败都是记录：

| 次 | 输入 | 结果 |
| --- | --- | --- |
| 1 | `! ExtOpt OptTS`，无 Hessian | 跑满 50 cycle 无 HURRAY，几何漂到 C–Cl 1.765/2.548 |
| 2 | `+ %geom Calc_Hess true` | `ERROR (SHARK)` —— **ExtOpt 下任何由 ORCA 自产 Hessian 的任务都失败**，不止 `Freq`/`NumFreq`（P1 记录的范围要据此放宽） |
| 3 | `+ InHess Read`，Hessian 由 `scripts/prrs_to_orca_hess.py` 从 PRRS 解析二阶导写出 | **HURRAY after 4 CYCLES** |

```
                        r(C-Cl)            对称感知 RMSD (|Aut| = 12)
T[E_m]   ORCA RS-P-RFO   2.2805 / 2.2827   T[E_m] vs A_θ[E_m] : 0.000146 A
A_θ[E_m] PRRS min-mode   2.2802 / 2.2825   T[E_m] vs T[E_QM]  : 0.047980 A
T[E_QM]  r2SCAN-3c       2.3645 / 2.3645   A_θ    vs T[E_QM]  : 0.048124 A
虚频     −535.66 / −535.82 cm⁻¹（差 0.03%）
```

两个完全不同的算子、独立起点，差 **0.000146 Å**，而各自离 DFT 鞍点 0.048 Å ——
**是彼此距离的 330 倍**（P1 那次是 12 倍）。方法这一行干净，分歧全在模型那一行。

`guess.hess` 在起点上就有一个真负模 **−5.661 eV/Å²**，印证第一次失败是优化器没信息，
不是面上没鞍点。**独立性未被削弱**：独立的是优化器与起点，Hessian 是势能面的性质。

**服务端第一次上生产**：`External energy and gradient 92.3%`、约 1.3 s/次，对 P1 那次的
`100.0%`、约 5 s/次。剩下的 1.3 s 是 ORCA 每 cycle 重 exec 一次 python，ExtOpt 的模型决定的。

### 完备性证书：第一次上真实 MLP 面，(ii) 是付不起的那一条

`docs/experiments/b04_completeness_audit.py`（新），R=0.5 Å / δ=0.1 Å / d=3、619 个采样点。

**首跑作废**：我取「最低三个本征值」取到了刚体零空间 —— `hessian_spectrum` 把刚体模
**投影掉但不降维**，它们以 ~1e-15 留在谱里。那一版 484 个 `other_basin` 是基造出来的假象。
留档 `runs/b04_completeness_audit_INVALID_nullspace/`，脚本改为按 `|λ| > 1e-8` 滤除并把
滤掉的个数记进 record（六原子非线性分子应当正好 6 个，实测正是 6）。

修正后：

```
三个最软真实模: 0.0186 / 0.0187 / 0.0392 eV/(A²·amu)
membership: basin 618 / 619            (i)   met
                                       (iii) met
L_used·δ = 0.4121  vs  g_min = 0.0353  (ii)  NOT met —— 差 11.7 倍
verdict: refused
```

要过 (ii) 需 δ < 0.0086 Å ⇒ 球内约 53 万点、每点一次 quench ⇒ 数周。覆盖论证的代价是 δ³。
**而且那还是乐观版本**：实测 L 是抽样上确界，即真 Lipschitz 的**下界**，模块明令不许冒充
上界；换成真上界只会更差。

对决策 4 的意义：(i)(iii) 是干净的正证据，卡住的是把它们接成结论的那一步。
**草案那条路在真实面上不是精度问题，是代价问题** —— 这是第一次有数说这句话。

### committor：判归有两个 bug，A/B 在代码里到不了

首次实跑（`runs/p1_committor1`）四个 bracket 全报 `A 0 / B 0 / other 16`。不是统计不够：

1. **判归读了没人写的键**。`record["target_microstate"]` 只有扫描路径在 `registry.admit`
   之后才写（`search.py:1813`），射击不走那条路 ⇒ A/B 两个分支不可达。
2. 修好①之后仍是 0：**`find_node` 拿到的是 dict**。`key_of()` 返回分解后的身份字典，
   而 `node.key` 是其中的哈希字符串，两者比较静默不匹配 —— 一个明明坐在盆里的几何被报成
   「没有盆」。

修法：`Registry.identify()`，只读的 `key_of()["key"] → find_node → _match_microstate`。
**不用 `admit`** —— 那会把射击当新微态写进网络，而协议写着 committor「只进模型报告，
不进判据」；探针不能写进它要观测的那张网。

顺带补 `bracket["committor"]["landings"]`：记下每发落在哪个微态。原来 `other 16` 分不出
「第三个盆」「淬火没终止」「判归有 bug」，这次正是被这个分不出来坑的。

**变异检验**：把 `identify` 换回旧键 ⇒ 新测试 `test_a_shot_that_terminates_is_judged_a_basin_not_other`
如期红，其余 5 条照绿；还原 ⇒ 6 条全绿。旧那条 shape 测试为什么挡不住：它的 `p_B` 断言在
`if terminated:` 里，而 `terminated` 从来没真过。

### committor 的科学结论：ridge 代理有偏，不是有噪声

> **2026-09-23 被推翻，下面原文保留。** 这批射击每发都被二次扰动（+0.4 eV 量级），p_B 量的不是 ridge 帧。修复后重跑见本文件末节。

被杀那次的 64 个落点都在盘上，重判**只读**、不重跑动力学
（`docs/experiments/p1_committor_rejudge.py` → `runs/p1_committor1/rejudge_2026-09-22.json`）：

```
bracket [0.3000, 0.3250]  B=m0001  {m0001:14, m0000:2}  p_B 0.875  95% [0.640, 0.965]
bracket [0.7375, 0.8250]  B=m0001  {m0001:16}           p_B 1.000  95% [0.806, 1.000]
bracket [0.7375, 0.8250]  B=m0001  {m0001:14, m0000:2}  p_B 0.875  95% [0.640, 0.965]
bracket [0.2750, 0.3000]  B=m0000  {m0000:14, m0001:2}  p_B 0.875  95% [0.640, 0.965]
```

64 发全部识别成功，零 unidentified。协议规定 N=16 只能答「p_B 是否接近 0 或 1」，
**答案是接近 1，四个窗口无一例外**。

`seek_saddle` 的 ridge 帧本是个假设 p_B≈1/2 的 1-sample 代理；实测说它已经过了分水岭。
**这是代理有偏，不是代理有噪声** —— committor 要给的那根误差棒给出来了。

射击是在 `ec32bc3df5d2bf14` 上打的，重判用 `63be7064f9e09f82` 的只读判归，
两个指纹都记在 json 里。

### 过期文本（本轮实测推翻）

| 写在哪 | 写的是 | 实测 |
| --- | --- | --- |
| `MANIFEST.md` b05 行 · `b05/README.md` · `b05/preflight.md` | gate 9 FAIL，自同构枚举器不 scale，19.8 s | **1149 ms，7/7 pass**（修它的是 2026-09-09 的逐点回溯剪枝） |
| `DECISIONS_PENDING.md:453` · `ROUNDS.md:162` · `REFERENCES.md:159` · `SESSION_2026-09-13.md:68` | `T[E_m]` 从未算过 | P1 的 T[E_m] **2026-09-04 就跑过**（`p1_qm/orca/p1_T_Em/`）；09-13 那次审计判错了 |


## 已完成：committor 射击的二次扰动 + P1 重跑 · 2026-09-23

`63be7064f9e09f82` → **`4014c4b267b25142`** · kasuga02 `make check` 绿（466 passed / 7 skipped）
· **digest 零差异**（332 / 2136 / 4）· 变异检验过 · 协议执行记录在 `P1_COMMITTOR_PROTOCOL.md` 末两节

### 结论

1. **09-22 的「ridge 代理有偏，不是有噪声」不成立。** 那批射击每发都被多推了一次，p_B 量的是
   「ridge 帧 + 再推一次」。修复后，四个质子转移窗口**不再同向**：

   ```
   stretch (4,5) [0.300, 0.325]   m0000→m0001   5 / 10 / 1   p_B 0.667  [0.417, 0.848]   修复前 0.875
   kick    (4,5) [0.7375, 0.825]  m0000→m0001  11 /  5 / 0   p_B 0.313  [0.142, 0.556]   修复前 1.000
   kick    (0,4) [0.7375, 0.825]  m0000→m0001   1 / 15 / 0   p_B 0.938  [0.717, 0.989]   修复前 0.875
   stretch (0,5) [0.275, 0.300]   m0001→m0000  15 /  1 / 0   p_B 0.063  [0.011, 0.283]   修复前 0.875
   ```

   按协议 N = 16 只答「是否接近 0 或 1」：两个窗口区间含 0.5（**未排除在分水岭上**），一个接近 1
   （产物侧），一个接近 0（源侧）。**ridge 帧作为 p_B≈1/2 的代理，在这四个窗口上有散布、没有看得出的单向偏移**；
   N = 16 分辨不了更细，这句话也只到这一层。
2. **「射击点往低幅度移」那条待落笔撤回** —— 它是为一个不存在的偏移开的药。
3. 扭转窗口 torsion (5,0,1,2) [1.57, 2.09] 16/16 落到构象微态 m0002，那是构象分水岭，不参与上面的判断。

### 缺陷

`shoot()` 把 ridge 帧交给 `run_trial`，而 `run_trial` 照常施加 probe —— 传进去的是扫描的 `probes[0]`，
即最低一档粗幅度。实测每发首帧：

| bracket | 再施加的 | 首帧多出 |
| --- | --- | --- |
| stretch (4,5) / (0,5) | 0.15 Å | 势能 +0.405 / +0.491 eV |
| kick (4,5) / (0,4) | 0.3 eV | 动能 0.44–0.72 eV（纯 300 K 热动量 0.15–0.34） |

势垒 419 meV。09-22 修的两处判归 bug 是真的，但它们修好之后露出来的数仍然是假的 —— 这是同一段代码的第三个缺陷。
旧测试挡不住它：`test_shots_differ_by_seed_and_reproduce` 只查 seed 和可复现，从没问过射击从哪出发。

**修法**：`run_trial(..., deliver=False)` 跳过位移 / kick / pulse，只留 `thermal_momenta`；`shoot` 传 False；
射击记录带 `probe_delivery: "skipped:committor_shot"`（普通 trial 不写这个键，记录逐字节不变）。
**变异检验**：`shoot` 改回 `deliver=True` ⇒ `test_a_shot_starts_from_the_frame_it_is_given` 红；还原 ⇒ 绿。

### 运行

`.scratch/queue_committor2.sh`，kasuga02 串行，fail-closed 在指纹上：`make check`（1 h 25 min）→
`runs/p1_committor2`（同一份 `.scratch/p1_committor_config.json`，1 h 16 min，`status: completed`，
outcome `completed 52 / rejected 8`）。**80 发射击首帧势能与 ridge 帧逐位相同，0 例 dV ≠ 0。**
「其他」1 发是 `unconfirmed`（拓扑没撑过自由段尾部），按协议不进分母。

### 还开着的

- 分辨 0.5 与 0.25 需 N ≳ 50 —— 协议「分辨率先行」已写明是另一批，要开须先写协议并跑前冻。
- `runs/p1_committor1/` 与其 `rejudge_2026-09-22.json` 保留作记录，**不可引用为 p_B**。

## 已完成：外部审查核对与续跑回滚重做 · 2026-09-24

`4014c4b2` → **`ae47c4d9142203dc`** · `476 passed, 1 skipped` · digest 零差异 332 / 2136 / 4 · **本机 kasuga01**（kasuga02 ssh 不可达，用户当场授权）。

| | 修了什么 | 检验 |
| --- | --- | --- |
| 续跑回滚 | 旧回滚以 network.json 为边界保留被引用的 attempt；`publish()` 每 trial 都跑，scan 中途的 reaction / ts 先上盘 ⇒ 记录有洞、trial id 撞目录（kill@4 → `FileExistsError: trials/t000002`），在飞微态已出队被丢。改为微态出队前写 `resume_point.json`，按计数回滚 trial / committor / climb / continuation / descent / ts / withheld | 逐字节验收参数化 kill@4 / kill@6 + 新测试「climb 里被杀」；回滚改回保留全部行 ⇒ kill@4 红、kill@6 仍绿（旧验收测不到的证明）；保留越界 climb ⇒ 红 |
| resume 输入绑定 | 核对 `input_structure_sha256`，缺字段也拒，拒绝在回滚前 | 新测试（含「盘上一个文件没少」）；删检查 ⇒ 红 |
| 覆盖网格 | `cover` 拒非正交归一基底 | 新测试（缩放 / 斜交）；删检查 ⇒ 红 |
| committor 记录 | 加 `p_B_conditioned_on`、`other_fraction`、`max_nve_drift_eV_atom`；`p_B` 不改名 | 断言 +3；`other_fraction` 分母写错 ⇒ 红 |
| 曲率文档 | `PRRS_CURVATURE_ROUTES.md` B.7 / C 节加勘误引用块，原文不删 | — |

核对全文、未修项与计划在 [`REVIEW_2026-09-24.md`](REVIEW_2026-09-24.md)。

### 追加 · 2026-09-25：p1_committor2 的结论对 t_free 不稳健

`runs/p1_tfree_sensitivity`：配对重打，只改自由段长度。scan 8 的 p_B 在 40 / 80 / 160 fs 下为 0.06 / 0.57 / 0.67，
本文件上文 09-23 节「scan 8 接近 0」**不可再引**（原文保留）；scan 2 / 6 / 15 稳定。表与读法见
`REVIEW_2026-09-24.md` §4½「实测」与 `P1_COMMITTOR_PROTOCOL.md` 末节。

## 已完成：首次进入 committor 第 0 段 · 2026-09-25 · kasuga02

`docs/experiments/p1_first_entry_stage0.py` → `runs/p1_first_entry_stage0/summary.json`；日志 `.scratch/first_entry/stage0.log`。
零 MD：从 `runs/p1_tfree_sensitivity/L0640`（64 发 × 640 步，scan 0 / 2 / 6 / 8）逐步重建滞回键图，**精确重建门全过**；
quench 核对用进入后最近的已存帧。k ∈ {10, 20, 40, 80, 120, 160} × L_max ∈ {160, 320, 640}。

| 结论 | 数 |
| --- | --- |
| k ≤ 20 是假进入 | scan 0 / 8：16/16 进入源侧，其后 15 / 12 发又离开 |
| k = 40 不稳 | 相对 k = 20 翻 12 / 8 发，再离开 6/16 |
| **k ≥ 80 平台** | k = 80 / 120 / 160 逐发判归相同；再离开 ≤ 1/16（k = 160 全 0）；L_max = 640 时未定 0 |
| quench 核对 | k ≥ 40 起 0 例不一致 |
| `q_PT` 交叉检查 | k ≥ 80 时进入 A / B 的 `q_PT` 符号相反、离 0 ≥ 0.25 |
| 与固定 t_free | 平台 ≈ 固定 160 fs（scan 0：0.81 / 0.79；scan 8：0.63 / 0.67；scan 2 / 6 相同） |

**更正**：「孤立分子 NVE 会一直回越」不成立（`REVIEW_2026-09-24.md` §4½ 读法 2、`P1_COMMITTOR_PROTOCOL.md` 末节已就地更正）。
**这些数只用于选参数，不可引为 p_B** —— 参数在同一批 64 发上选，报数要等第 1 段新 seed。

## 已完成：首次进入 committor 确认批（A3 / A4）· 2026-09-26

`runs/p1_first_entry_confirm/verdict.json`（`docs/experiments/p1_first_entry_confirm.py --analyse`，用户 2026-09-26 16:23 运行）。射击 2026-09-25 19:54 → 09-26 10:08（14 h 14 min，kasuga02）。协议 sha256 `b685120f…` 与冻结时一致，实现指纹 `ae47c4d9`。每窗口 N = 50，新 seed（`seed0_published + 100000 + i`），1280 步不提前停、每步存帧；k = 160 步。

| scan | 坐标 | A / B / 未定 | **p_B_given_AB** | 95% Wilson（本窗口自己的） | V2 翻转 80↔160 / 160↔320 | V3 再离开 | V4 quench 不一致 | V5 未定 | 可引用 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | stretch (4,5) + | 12 / 38 / 0 | **0.76** | [0.63, 0.86] | 1 / 1 | 6% | 0% | 0% | ✅ |
| 2 | kick (4,5) + | 44 / 6 / 0 | **0.12** | [0.06, 0.24] | 1 / 0 | 6% | 0% | 0% | ✅ |
| 6 | kick (0,4) − | 0 / 50 / 0 | **1.00** | [0.93, 1.00] | 0 / 0 | 0% | 0% | 0% | ✅ |
| 8 | stretch (0,5) + | 14 / 36 / 0 | **0.72** | [0.58, 0.83] | 1 / 2 | 8% | 0% | 0% | ✅ |

- **V1（整批）通过**；四个窗口 V2–V5 **全过**，按冻结协议 §4 四个窗口的 `p_B_given_AB` 都可作为确认结果引用，
  同时报原始计数与各自的区间（不声称四个区间同时 95% 覆盖）。
- 区间重叠（预定检查）：k = 80 / 160 / 320 三档两两重叠，四个窗口都是。
- 最紧的一项是 V3：scan 8 再离开 4 / 50 = 8%，阈值 10%。
- quench 核对 200 发里 0 例不一致；`q_PT` 交叉检查：进入 A 与进入 B 的 `q_PT` 符号相反、离 0 ≥ 0.24。
- 进入时刻中位数：scan 0 / 8 约 312 / 315 步（78 fs），最晚 661 / 744 步；scan 2 为 175 步；scan 6 每发恰在第 160 步（一出发就在 B 图上）。
- 200 发里 199 发 `completed`、1 发 `unconfirmed`（第 1280 步末尾拓扑没撑住）—— 首次进入按逐步键图判，不受它影响。
  最大 NVE 漂移 9.1e-4 eV/atom（×9 原子 ≈ 8 meV，势垒 419 meV）。
- **读法**：四个窗口的 ridge 帧都不在分水岭上 —— scan 0 / 8 偏 B 侧（区间下界 0.63 / 0.58），scan 2 偏 A 侧（上界 0.24），scan 6 已在 B。
  已撤回的 40 fs 值里，scan 8 的「0.06 · 接近 0」与确认值 0.72 方向相反。

## 已完成：首次进入 committor 生产实现（A5）· 2026-09-26

指纹 `ae47c4d9` → **`b601b124e14538a2`**。`run_trial(first_entry=…)` 逐步计数、连续 k 步在 A / B 图上即停，之后原有 quench 就是入口核对；`shoot` 透传；配置 `committor_entry_steps`（k，0 = 关，默认）/ `committor_max_steps`（1280）；search 的 committor 分支按进入侧判，核对不一致只计数不改判，没有 B 或两端同图时退回固定规则并在 `termination` 写明；R1：`p_B` → `p_B_given_AB`。

kasuga02：`make check` 绿、digest 零差异 332 / 2136 / 4；变异检验两条都红（`>=`→`>` 进入步差一；删 `break` 跑满上限）；**真实数据核对** `runs/p1_first_entry_production_check`：生产路径重打 A3 的 200 发（同帧同 seed、提前停），**进入侧与进入步 200/200 相同、入口 quench 落点 200/200 相同**，约 23 min。

测试（`tests/test_committor.py`）：双阱上提前停止 = 同 seed 全长离线判（k = 5 / 20）、永不进入跑满上限记未定、
端到端记录（k = 0 用原夹具；k = 10 用幅度 0.1 / 0.9 —— 原夹具细化后唯一窗口高端是 0.600 的瞬态越界，没有 B，
已由诊断实测确认），并断言至少一个 bracket 真正走了首次进入。

## 已完成：支线 B0 · C2 普查 · C3 探索批 · 2026-09-26 · kasuga02（用户运行）

### B0 · 搜索的 bracket 随 t_free 移动 —— `runs/p1_bracket_tfree`

`docs/experiments/p1_bracket_tfree.py`：只重放各 scan 已采样的幅度，三档长度统一按 quench 落点判。
160 步基线**逐个复现已发布端点**（门通过）。窗口（相邻 returned → departed）：

| scan | 40 fs（发布，160 步） | 20 fs（80 步） | 80 fs（320 步） |
| --- | --- | --- | --- |
| 0 | [0.300, 0.325] | **[0.35, 0.50]** | **[0.15, 0.25] + [0.300, 0.325]**（0.25 离开、0.30 又回来） |
| 2 | [0.7375, 0.825] | **[0.65, 0.7375]** | [0.7375, 0.825] |
| 6 | [0.7375, 0.825] | **[0.825, 1.0]** | [0.7375, 0.825] |
| 8 | [0.275, 0.300] | **[0.30, 0.35]** | **[0.30, 0.35]** |

- **20 fs 时四个 scan 的窗口全挪**；80 fs 时 scan 2 / 6 与 40 fs 相同，**scan 0 / 8 不同** —— 与 committor 的格局一致（2 / 6 早定，0 / 8 慢定）。
- scan 0 在 80 fs 出现 `0.25 离开 → 0.30 返回`：一个现行规则不建窗的 reverse 窗口。
- 结论：搜索判「回到 / 离开」同样依赖 t_free；**按计划不在此改**，另立决策 R8（`DECISIONS_PENDING.md` 末节）。

### C2 · 窗口普查 —— `runs/window_census/census.json`

`docs/experiments/window_census.py`，零 MD，`runs/` 全部有 scan 的网络：**688 个 scan**，相邻不同结果中
现行规则建窗的 **151**；它忽略的：reverse **10**、switch（两个不同产物）**1**、other（如 departed ↔ transient）**13**，
以及**一侧 unresolved 的 154**；细化被第三种结果打断的 bracket 14 个。reverse 主要在 ethanol / p0 系列（`ridge → returned`）
与 `p1_malonaldehyde`（`departed → returned`）；unresolved 最多的是 b05 probe3 / probe4 与 p2 T300 系列。
**读法**：窗口泛化（reverse / switch / other）多出约 24 个，比现行 151 多 16%；**更大的缺口是 unresolved（154）** ——
先做 C1（把未判定细分记录）比先做 C2 的窗口泛化更有信息量。

### C3 · 沿鞍点不稳定模的 committor（探索批，**不可引用**）—— `runs/p1_path_committor_explore`

`docs/experiments/p1_path_committor_explore.py`，首次进入 k = 160 / 上限 1280，每帧 16 发。A = m0000、B = m0001。
plus 分支几何因 D1 丢失，用 minus 首帧对鞍点镜像代替；镜像帧 q_PT −0.304 与记录的 plus 首步 −0.304 一致。

| 沿不稳定模 | q_PT | ts0000 A / B | ts0001 A / B |
| --- | --- | --- | --- |
| saddle − δ | −0.30 | 0 / 16 | 0 / 16 |
| − 2/3 δ | −0.20 | 0 / 16 | 0 / 16 |
| − 1/3 δ | −0.10 | 0 / 16 | 0 / 16 |
| **saddle** | 0.00 | **10 / 6（0.38）** | **12 / 4（0.25）** |
| + 1/3 δ | +0.10 | 16 / 0 | 15 / 1 |
| + 2/3 δ | +0.20 | 16 / 0 | 16 / 0 |
| + δ | +0.30 | 16 / 0 | 16 / 0 |

- p_B 的全部变化发生在 **|q_PT| < 0.1** 之内（q = −0.1 已 16/16 进 B，+0.1 已 31/32 进 A）。
- 鞍点本身 0.38 / 0.25（合并 10/32 = 0.31）。N = 16 不足以判断是否偏离对称预期的 0.5；两个鞍点可能是同一个（q 值逐位相同）。
- 对固定批（N = 50）的含义：帧要加密在 q_PT ∈ [−0.1, +0.1]；现在 ±δ/3 的间距太粗。固定批须先立协议、跑前冻（R5）。

## 已完成：C1 样本分类记录 + D1 鞍点下降几何 · 2026-09-26

指纹 `b601b124` → **`99edddde34294b69`**。D1：鞍点双侧下降改用**一个** `PathRecorder` 管两支，帧计数连续，`descent.extxyz` 不再被第二支覆盖（此前历来只剩 minus 支的几何）。 C1：scan 每个样本加 `outcome`（returned / departed / ridge / transient / unresolved）与 `outcome_reason`（未定时的失败码或 `admission:…`），只记录不路由。

测试 `tests/test_scan_records.py`（双阱真实搜索里两支帧都在且与日志一一对应；各分类单测；端到端字段一致）。变异检验两条都红（D1：记录器挪回循环 → 帧数 ≠ 行数；C1：returned / departed 对调 → 两条测试红），脚本还原并核指纹 `99edddde`；变异检验在 **yayoi**（PBS 作业）上跑。未变异全量：**用户报告通过**（未贴 digest 行）。

C1 的来由：窗口普查（688 scan）中一侧 unresolved 的相邻对 154 个，与建窗的 151 个一样多 —— 现在每个未定样本都带着原因落盘，
下次普查可以把 154 拆开。

## 已完成：b04 的 t_free 敏感性（B2，R8 依据）· 2026-09-26 · kasuga02（用户运行）

`docs/experiments/bracket_tfree.py --run runs/b04_frozen2 --calculator benchmarks.b04_sn2_chloride.calculator:polar1 --lengths 100,400`
→ `runs/b04_bracket_tfree/summary.json`。18 scan、64 个已采样幅度 × 3 档长度；T = 0 K，逐发确定。
**基线门通过**：50 fs 重放逐个复现 64 个已发布端点。

| scan | 50 fs（冻结 run，200 步） | 25 fs | 100 fs |
| --- | --- | --- | --- |
| 0 stretch (0,1) + / 9 stretch (0,5) +（SN2 正反向） | [0.75, 0.80] | 相同 | **[0.3, 0.6]**；0.6 / 0.7 / 0.75 离开，**0.8 / 1.0（50 fs 时离开）回来了** |
| 6 kick (0,1) + / 15 kick (0,5) + | 无窗口 | 相同 | **新窗口 [0.5, 1.0]**（1.0 eV 离开） |
| 其余 14 个 scan | 无窗口 | 相同 | 相同 |

- **缩短到 25 fs 无一改变**；**加长到 100 fs 四个 scan 变了**，而且变在 b04 唯一的两个 bracket 上。
- scan 0 / 9 在 100 fs 同时出现 returned → departed（[0.3, 0.6]）和 departed → returned（[0.75, 0.8]，现行规则不建窗）：
  恒等 SN2（Cl⁻ + CH₃Cl）两端对称，T = 0 且无耗散，越过去的复合物在 50–100 fs 间又回越 —— quench 读到的是「100 fs 那一刻在哪边」。
- kick 1.0 eV 要到 100 fs 才定下来离开，50 fs 时判为回到。
- **只比了窗口**：`check.py` 要整张网络，重放给不出，所以 b04 冻结判据本身是否翻转**没有测**。定性上 100 fs 仍观察到离开事件，但出现在不同幅度。

## 已完成：电荷 / 多重度整数校验 · 2026-09-27 · 本机（用户授权）

指纹 `99edddde` → **`b2c056524a704ecd`**。`SearchConfig.__post_init__`：`total_charge` 设了值必须是 `int`，`multiplicity` 设了值必须是 >= 1 的 `int`；与其他整数字段一样用 `type(v) is int`，所以 `bool` 也拒。

来由（外部审查 2026-09-27 提出，核实为真）：身份层经 `network.Registry` 把 `total_charge` 原样交给 `chemical_key`，势能侧 `_stamped_with_charge` 取 `int()`。`1.5` 会以 1.5 记身份、以 1 算能量；`True` 静默成 +1。`_stamped_with_charge` 防的是「结构与配置两种说法」，防不了「配置自己一种说法、两处读法」。

测试 `tests/test_validation.py::test_charge_and_multiplicity_must_be_integers_when_set`，修复前先跑红（`DID NOT RAISE`，卡在 `total_charge=1.5`）。全量 484 passed / 1 skipped，fmt-check / lint 绿，digest 零差异 332 / 2136 / 4。仓内配置（b04、b05 的 `config.json` 与 `preflight.py`）都是整数，不受影响。

## 已完成:路径成熟层 M0(清点)+ M1(协议与落盘)· 2026-10-03 · 本机

独立模块 `mature/`(骨架同日落地,见施工计划 §9 修订:平铺模块改独立模块,单向 `import prrs`)。
**本块不碰 `src/prrs/`,指纹停在 `ebff4ec0`,digest 零差异(332 / 2136 / 4,见下)。**

### M0:preflight 清点(只读,四个指定 run,不递归 runs/)

`docs/experiments/maturation_preflight.py` → `docs/experiments/maturation_preflight/preflight.json`
(+ `baseline_2026-10-03.json`)。逐 run 判定(支持/不支持 + 原因):

- **runs/p1_committor2 = supported(源 run)**:config/network 完整,模型 checkpoint sha 一致,
  140/140 trial 带全原子 response.traj(抽查 n_atoms = 输入 n_atoms),2 个 ts_candidates,
  1 个通道 `rc0000`(退化自环,broken/formed = (4,5)/(0,5),event_key `a04410b76e16b117`)。
  注意事项:run 指纹 `4014c4b2` 是历史来源;M3 射击跑在当前树,指纹在批次 manifest 记录。
- **runs/p1_first_entry_confirm = supported(确认批参照)**:batch.json 协议哈希 `b685120f…`、
  N=50 × 4 窗口、k=160 / cap=1280、seed_offset=100000;200 个 trial 目录 = N × 窗口数。
  **区域是无序对 {m0000, m0001}**:窗口 8 的 A/B 标签随源态对调(记录为
  `windows_with_swapped_labels: [8]`,不是第二个区域)。
- **runs/p1_path_committor_explore = supported_reference_only(探索批,不可引)**:
  summary 自带 not-citable 标记;2 个鞍点 × N=16 逐帧计数。
- **runs/pars_ground_truth = supported(精化基线)**:27/27 记录 `ci_neb_converged`,
  barrier 0.419128–0.419208 eV(跨度 0.08 meV < 1 meV,"同一鞍点基线"成立);
  来源 run = p1_committor1(`ec32bc3d`)/ p1_round7(`8eb7d82f`),provenance 链接完整。
- **区域资格(映射微态 + 键图,不用 chemical_key)**:m0000/m0001 的 `chemical_key` 相同
  (P1 退化自环的预期),但 encode 后图哈希不同(各 8 条边)→ **首次进入射击 eligible**。
- **PARS selection:不存在** → M2 前须用 10-02 修复后实现重新 `pars select` 到新目录并记录哈希。
- **基线**:`make test` = 537 passed / 7 skipped(373 s);`make digest` 零差异 332 / 2136 / 4。

### M1:schema 与落盘(`mature/src/mature/schema.py` + `io.py`)

- **schema.py**:MaturationProtocol(batch_role 四类)、SourceRef(含逐文件 sha256)、
  RegionDefinition(恰好 A/B 两侧,**同图区域定义期即拒** = unsupported_same_graph)、
  DynamicsSpec(NVE / MB / VV,T>0,**pulse_steps≠0 拒绝**)、CandidateRule(端点弧长位置拒绝,
  全 H 保留钉死)、ShootingRule(确认/抽查批必须带跑前冻结的 selection_rule;**漂移容差只有
  零号测量批允许留空**)、Budget(至少一项限额)、PathProposal(proposal_kind 四类,
  observed_frame 必须带 frame_index,插值候选不得伪装真实 MD 帧)、ShotRecord(八种 outcome
  分列,A/B 必须带证据;排除项必须带原因)、ShootingEstimate(p_B_given_AB 在 A+B=0 时为 null,
  Wilson 区间随行;带分类只用区间不用点估计)、TransitionDiagnosis(词表里没有"已验证 TS")。
  单位钉死 eV/Å/fs/K;canonical JSON 哈希全程 allow_nan=False。
- **io.py**:输出目录布局照计划 §3.2;原子写(tmp+fsync+replace);completion 只认
  commits/ 提交标记;**同配置任务表可重建**(种子 = seed_base+序号,按 (proposal_id, shot)
  排序,与执行顺序无关);换协议哈希/源文件内容/几何内容/实现指纹**拒绝续跑**;
  半写 shot 目录隔离到 interrupted/,不删除;重复 proposal 同内容幂等、异内容拒绝。
- **测试**:mature/ 自有 **53 passed**(边界 6 + M1 47),fmt/lint 绿(openmm_dev + .venv-tools)。
  统计测试引用已发布数字:Wilson(38,50) = 确认批 verdict 的 [0.6259, 0.8570],
  Wilson(8,16) ≈ [0.28, 0.72](可行性评估修正第 4 条)。
  边界变异检验:`import mature` 塞入 pars/src/pars → 边界测试红,移除后恢复绿。
- **不做**(属 M2–M4):候选提取、射击接线、CLI 实装(现 fail-closed 退 2)、诊断写盘编排。

## 已完成：PARS 10-05 复审三项修复 + 真实 b04（POLAR-1）核验 · 2026-10-05（2026-10-08 补记）· 本机 / 183

问题清单见 [REREVIEW_2026-10-05.md](REREVIEW_2026-10-05.md)（两项 P1、一项 P2）。修复于 10-05 12:0x 落地，
**当天没有写文档**，本节是 10-08 根据代码、测试和结果目录补的记录。
PARS 源码指纹：修复前 `57fcad4d…`（REREVIEW 记录）→ 修复后 **`d7493683c7683cd3…`**
（10-08 现算；`pars/src` 最新文件时间 10-05 12:09，晚间真实 b04 跑在这份代码上）。不碰 `src/prrs/`。

### 修了什么

- **P1 电子态交接**：新增 `pars/src/pars/electronic.py`；`adapters/prrs_run.py:200` `_declared_electronic_state`
  从 manifest 读 total_charge / multiplicity，与逐帧、端点的 `atoms.info` 交叉核对，冲突、部分打标、
  缺失或非法都在计算前拒绝；状态随 candidate → 代表文件 → seed → 后端所有构型 → 输出一路传递，
  `refinement/ase_neb.py:298` 把 `electronic_state` 写进精修记录（没有声明时写 `unspecified`）。
- **P1 Hessian 采样不确定度**：`ase_neb.py:167-194` `_curvature` 每次正负位移的 force evaluation 当场检查
  `force_uncertainty_eV_A`：非有限或负值直接拒绝；设了拒绝阈值时，缺测或超阈值不交付曲率成功；
  采样的不确定度计入报告的最大值。
- **P2 先核验再过滤**：`paths/ingest.py:32` 改为遍历 `run.checked_attempts()`，账本与 result.json 冲突在过滤前就报错；
  一致的 elastic 试验照常跳过。
- 回归测试：`pars/tests/test_rereview_regressions.py`（电子态 9 条、Hessian 门 3 条、过滤顺序 2 条）。

### 真实分子核验：b04（Cl⁻ + CH₃Cl，charge −1，单重态），MACE-POLAR-1，本机 GPU，10-05 晚

脚本 `benchmarks/pars_integrity/real_b04/`，结果 `benchmarks/pars_integrity/results/2026-10-05/real_b04/`。
`logged_polar1.py` 记录 POLAR-1 实际收到的每个结构的 (charge, spin)。

| 条件 | 代表 | 结果 | POLAR-1 收到的电子态 |
| --- | --- | --- | --- |
| **修复前做法**（`counterfactual.py`：剥掉 electronic_state，probe4_ledger） | t3 | CI 阶段被拒，退回 plain NEB，势垒 0.618 eV（CI 0.715） | `(None, None)` |
| 同上 | t4 | `failed`：连接检查下降不能唯一匹配两端盆地，不报势垒 | `(None, None)` |
| **修复后**，probe4_ledger（`b04_real/`） | t3、t4 | 都是 `ci_neb_converged`，0.59670 / 0.59681 eV，1 个负模，连接 `matched` | `(-1, 1)` × 8969 |
| 修复后，b04_probe | t3 | `ci_neb_converged` 0.59666 eV，1 负模，`matched` | `(-1, 1)` × 5222 |
| 修复后，b04_probe2_openceiling | t3 | `ci_neb_converged` 0.59666 eV，1 负模，`matched` | `(-1, 1)` × 5222 |
| 修复后，b04_probe3_rigidbond | t3 | `ci_neb_converged` 0.59670 eV，1 负模，`matched` | 未落 seen.json（批次中断） |

读法：
- 电子态丢失**在真实带电分子上确实改变结果**（同一组代表：修复前一个退化成 plain NEB、一个连接失败；
  修复后全部 CI 收敛并连上两端）。修复前的结果不是"数值略差"，而是错误 PES 上的计算。
- 修复后 5 条精修势垒跨度 0.15 meV；10-08 补齐到 14 条后跨度 2.1 meV（见下）。b04 仍是**单鞍点**数据，
  和 P1 一样**不能**验证"选择保留多机理多样性"。
- 势垒、负模数值不进 b04 判据（CRITERIA.md §2 规矩 2）；与 QM 的比较见 `docs/experiments/b04_qm/`。
- `ingest_audit.json`：用修复后的严格导入把 `runs/` 下 97 个 run 全部导入一遍，没有一个被拒；
  新格式 run 的电子态都读出来了（b04/b05 为 (−1, 1)，其余 (0, 1)），旧格式与无账本的照常跳过。

**10-08 补齐（183，2080 Ti）**：`gpu_batch.sh … b04_probe3_rigidbond b04_frozen1 b04_frozen2`，16:40–17:47 一批跑完
（log `benchmarks/pars_integrity/results/2026-10-08/gpu_batch_183.log`，末尾 `BATCH_DONE`）。9 个代表约 7 min/个。

| run | 代表 | 结果 |
| --- | --- | --- |
| b04_probe3_rigidbond | t4 | `ci_neb_converged` 0.59681 eV，1 负模，`matched` |
| b04_frozen1 | t3 / t4 / t35 / t36 | 全部 `ci_neb_converged`，0.59670 / 0.59681 / 0.59866 / 0.59878 eV，各 1 负模，`matched` |
| b04_frozen2 | t3 / t4 / t35 / t36 | 与 frozen1 逐个相同（到 0.01 meV） |

全部 9 条记录的 `electronic_state` 为 (−1, 1)。合计 b04 上修复后 14 条精修全部 CI 收敛、单负模、连接 matched。
t35/t36 比 t3/t4 高约 2 meV，可能是端点极小不同，**未核实**；不当作第二个鞍点。
POLAR-1 收到的电子态的**独立计数**（`seen.json`）对这 9 条**丢失**：后面两次续跑（全部"already done"、没有计算）
在退出时把 `{}` 覆盖了上去。`logged_polar1.py` 已改为累加写入。这 9 条的电子态证据只剩精修记录里的 `electronic_state`。

**更正（10-08 当天）**：我一度在本文记录"183 上补跑无声死亡、原因未查明"，**这是错的**。
两台机器的 NFS 都挂载为 `acregmin=acdirmin=21600`（属性缓存 6 h）+ `nocto`，181 上看到的 log 长度和目录内容是
过期缓存（log 显示 618 B，183 上实为 2493 B 且有 `BATCH_DONE`）。之后两次"重跑"都只是发现结果已在、几十秒内正常退出。
10-05 那一批在 probe3 处停下的原因仍不清楚（当时全在本机），但与此无关。
教训：**跨主机读结果要在写结果的那台机器上读**（ssh 过去读），不要从另一台的 NFS 视图下结论。

### 验证

- 2026-10-08 · 183（kasuga02）：`make -C pars fmt-check lint` 通过（47 文件）；
  `pytest tests/test_rereview_regressions.py tests/test_review_regressions.py tests/test_ingest.py tests/test_cli.py`
  **81 passed**，18 s。日志 `benchmarks/pars_integrity/results/2026-10-08/smoke_183.log`。
  本次没有跑 PARS 全量测试，也没有跑根目录 `make check`（不涉及 `src/prrs/`）。
- `benchmarks/pars_integrity/rereview.py` 的断言**固定的是修复前的缺陷**，修复后它理应失败，不能再当回归用；
  另外它要执行 `git rev-parse HEAD`，在 NFS 工作区上跑不起来（10-08 在 183 上就停在这一步）。
  修复后的行为以 `test_rereview_regressions.py` 和上面的真实 b04 对照为准。
- `real_b04/gpu_batch.sh` 原本写死旧会话 scratchpad 路径；10-08 改成仓库内路径 + `<out dir>` 参数，可续跑。
  10-08 用它在 183 补齐了剩下 9 个代表（见上）；之后又加了 `RUN_EXIT` 退出码记录、`TRACE_EVERY` 栈转储（`traced_pars.py`），
  以及 `logged_polar1.py` 的累加写入。

### 仍未处理

- REREVIEW 末尾的两项 provenance 缺口（PRRS 绑定 factory 指纹相同但实例状态不同；PARS 模型非权重 scale 不同但指纹相同），
  见 [DEFERRED_FIXES_2026-10-04.md](DEFERRED_FIXES_2026-10-04.md)。
- 多机理选择有效性仍缺多鞍点参照集（b05 方向），见项目总览 §9 第 3 条。

## 已完成：PRRS 谱排序改为按实际执行的运动打分（2026-10-09）

起因：b09 探针 1（`runs/b09_probe1`）0 反应，卡宾 C → 环烯碳的方向从没被选中；外部审查（用户转来）指出排序打分用端点梯度，
而位移类 probe 由 `apply` 以整片段刚体运动执行，两者不是同一个向量。

- `src/prrs/perturbations.py`：新增 `delivered_direction`（位移类按 `apply` 的刚体平移 / 转动求闭式切向量，环内回退仍为 g/m；
  kick / pulse 沿 g/√m 不变；投影掉刚体分量）；`spectral_overlap` 改签名为 `(atoms, family, indices, modes, bond_scale, active, graph)`；
  `propose` 传入已建好的图。先试过对 `apply` 做数值差分，b09 上每个方向 7.5 s（`collateral_metrics` 占大头），弃用。
- `tests/test_spectral_rank.py`：闭式切向量与 `apply` 数值差分逐族对齐（compress / stretch / bend / 两个 torsion），
  分子间 compress 的分数不再等于端点梯度，kick 保持梯度方向。相关 8 个测试文件 261 passed（181，CPU）。
- 实测（`benchmarks/carbene_ambimodal/results/2026-10-09/rank_b0{4,9}.log`）：b04 的 C–Cl 拉伸仍排第一（0.990）；
  b09 的 `compress:2-15` 0.186 → 0.515，名次 21 → 36/84，**仍选不上**。分子间 compress 都是两分子沿轴线的刚体平移，
  软模先验只看轴向，分不出成键原子对。详见 `benchmarks/carbene_ambimodal/probe_b09/README.md`。

### 同一份审查的其余各条（10-09 下午由 Codex 修，用户选方案 B）

见 [SEARCH_DESIGN_REPAIR_2026-10-09.md](SEARCH_DESIGN_REPAIR_2026-10-09.md)：
第 2 条受限软模重新对角化、第 3/4 条方案 B（源微态内先全部粗扫再细化）、第 5 条环内保键交付（门限不变）。
我在 181 上复跑相关 11 个测试文件 277 passed、ruff 通过；Codex 的 `make check` 592 passed / 2 skipped 日志与三组变异检验
原放在 `/tmp`，10-09 拷进 `benchmarks/carbene_ambimodal/results/2026-10-09/search_repair/`。

仍未处理：双分子体系的成键排序先验（没有替代先验，先做全候选粗扫基线）；
`runner.py:2279` 的 `log.info("... %d relay(s)", relays)` 把列表传给 `%d`（日志报 Logging error，不影响计算）。

## 已完成：b09 compress 族全候选粗扫（方案 B，2026-10-09 15:57 → 10-10 06:51，183）

`runs/b09_compress_coverage`，配置 `benchmarks/carbene_ambimodal/probe_b09/config_compress_coverage.json`。
336/336 trial，14 h 53 min（2.66 min/trial）；84 个 compress 方向全部执行，账本无截断。
0 个反应、0 个门拒；37 个 `ts_candidate` 全是复合物软模（离源微态 ≤ 0.23 Å，+0.4–23 meV）。
卡宾 C → Cb（`compress:2-15`）在最大振幅 1.0 Å 时推到 C:–Cb 2.12 Å，POLAR-1 TS-2OMe 是 1.98 Å，随后弹回：
阴性结论只覆盖"振幅 ≤ 1.0 Å 的刚体压缩"。分析脚本 `benchmarks/carbene_ambimodal/analyze_coverage_b09.py`，
结果 `benchmarks/carbene_ambimodal/results/2026-10-10/`，详见 `benchmarks/carbene_ambimodal/probe_b09/README.md`。
- 10-10 加大力度冒烟（`benchmarks/carbene_ambimodal/push_b09.py`，8 个 trial）：compress 2-15 推 1.5 Å（推完 +6.16 eV）
  撞出 (2+1) 环丙烷（比源低 1.17 eV），但淬火 500 步未收敛被记为 `quench_failed`，search 不会记成反应；1.2 Å 与
  0.5 / 1.0 eV 的 kick 都弹回；compress 2-13 1.5 Å 被碎片数门拒。覆盖扫描的 26 个淬火失败里没有卡宾–环成键，
  但 `t000327` 有一个甲氧基 H 被撞飞（C5–H6 末帧 10.22 Å），同样未收录——"没有漏掉反应"不成立（10-10 外部复核指出，已核实）。
  p001 单纯力优化：同一起点 FIRE 834 步、BFGS 139 步收敛到同一个环丙烷（`quench_p001/continue_quench.log`；
  导出文件的能量/力写错过，已修脚本重跑）；未经 runner 的 polish / 极小确认，完整验收见 replay_quench。
  更正：此前"TS-2OMe 能垒 ≈ 0.16 eV"错误，POLAR-1 TS-2OMe 相对本次源微态 +0.477 eV。
  淬火重放冒烟（`benchmarks/carbene_ambimodal/replay_quench_b09.py`，quench_steps 1500）：p001 经完整 runner 链判为 reactive、极小确认通过，与 POLAR-1 (2+1) 产物 RMSD 0.029 Å，134 s。
- 10-10 覆盖扫描 26 个 `quench_failed` 重放（`replay_quench_b09.py`，quench_steps 1500，62 min，183）：25 个收敛（527–1102 步），
  全是无键变化的构象重排、0 阶极小，22 个落到比源低 0.21–0.22 eV 的构象；t000327（H 被撞飞，+7.07 eV）2129 步后软模 polish
  仍未定住。已收敛的 25 个重放都不是反应，即这 25 个没有因淬火上限漏收反应；t000327 仍未决，不能当作已证实无反应（2026-10-11 按审计 §5 收窄）。结果 `benchmarks/carbene_ambimodal/results/2026-10-10/replay_quench_coverage/`。
- 10-10 交接：`docs/HANDOFF_B09_GAP_2026-10-10.md`（b09 定位为终极验收，只测差距；用户在新对话中分析差距）。
- 10-10 差距分析：`docs/B09_GAP_ANALYSIS_2026-10-10.md`（读完交接必读清单后写；结论：盲搜第 1 级"发现 (2+1) 连接"未过、
  验收判据未写；差距分参数级 / 设计稿 §3.4 协同扰动等未实现 / 同键图多鞍点超出发现层三类；用户 10-10 表示不着急做）。

## 已完成：项目第一轮审计（2026-10-11；未实施修复）

报告 `docs/AUDIT_2026-10-11.md`，反例脚本 `benchmarks/project_audit/reproduce.py`，
有效证据 `benchmarks/project_audit/results/2026-10-11/counterexamples-v2/report.json`。
确认四项数据契约缺口：PARS 代表路径逐帧元素序列未验、重复 attempt ID 遮蔽账本冲突、
mature 提交不绑定本批 protocol/task seed/实现、commit marker 内容 hash 在恢复与汇总时未重验。
另指出 b09 分析中“TS 来源起点最有利／差距下界”的推论缺依据；已知 factory/model 指纹缺口仍按既有裁决延后。
四项新缺口均为合成 CPU 反例，没有证明历史运行已污染。算法源码、生产配置和历史结果未改。
PRRS/PARS 回归 805 passed / 2 skipped；mature 53 passed；共 858 passed / 2 skipped。
三包源码/测试 122 文件格式与 lint 通过，审计脚本格式/lint 通过；历史 digest 332 文件 / 2136 帧 / 4 事件对不变，254 次既有迁移比较通过。
首次三包联合 pytest 因同名 conftest/test 模块出现收集冲突；有效结果来自 PRRS/PARS 联合与 mature 独立调用，详见报告。

同日算法补充：审计报告 §6 汇总 ALG-01…10，区分当前必须修正、一般双分子／多鞍点／分支目标所需能力，
以及需测出瓶颈才实施的性能候选。覆盖严重度跨族可比性、几何反应管、接近姿态与协同提案、覆盖保底与响应反馈、
振幅域／时窗、完整淬火链、鞍点种子分配、全链成本账本、长扫描续跑与性能、射击成熟层及模型资格。
新增 `benchmarks/project_audit/algorithm_probes.py` 与 `results/2026-10-11/algorithm-probes/report.json`，
直接复现混合单位严重度的任意排序和同键序列不同几何路径的距离为零／归同管。
证据是描述符级 CPU 探针，没有运行势能／生产候选导入或证明真实鞍点丢失；脚本格式与 lint 通过。
源码、既有设计、生产配置和历史记录未改；未启动 GPU／新科学实验，未重复跑前述同源码回归。

随后补充甲基／普通扭转／碎片软转动的节流讨论稿：建议分账、限额、去重复升级，
保留真实构象与接近姿态入口，不按甲基名一律忽略、冻结或降低端点接受门。
源码检查注明谱分数会覆盖原 genuine/rotor 排序字段；尚无真实谱反例或历史浪费因果结论。
记录现有周期过滤、free 判别与构象容量已实现的范围，以及需要成本分解和同预算对照才冻结的规则。
用户尚未确定方案，本次未实施代码／配置改动或新科学实验。

同日 b09 差距续查：按 `CLAUDE.md` 在结果写入主机 183（kasuga02）经 SSH 只读核验，
新证据 `benchmarks/project_audit/results/2026-10-11/b09-writer-snapshot.json`，
读取脚本 `benchmarks/project_audit/b09_snapshot.py`，报告 §7。
07:27 JST 的快照保存 41 个文件哈希和单文件读取稳定性；未加载势能模型或重算轨迹。
compress 主批 manifest 与当前 PRRS hash 一致，84/84 方向、336/336 coarse、无重复 attempt，
反应边／反应通道／已连接 TS 通道均为 0；26 个重放的 summary 与逐 trial JSON 相同，
25 个已确认极小且无键变化，t000327 仍未决；p001 的有效产物记录仍属于人工强推诊断。
旧 probe1 是旧源码结果；coarse_smoke 留着 running 状态但其历史 PID 不存在，不能当作活跃进展。
差距清单分开参照资格、盲搜发现、同连接多入口保留、精修连接、分岔与 QM 门；
明确当前经典 shoot 与论文准经典／溶剂条件尚不匹配，四个参照重优化点不等于完成盲搜或分支认证。
仅补审计与依赖顺序，没有实施修复、更改协议或运行新分子实验；新脚本格式／lint 通过。

随后对“拉伸像盲试、探索又不足、每个已执行候选需要 MD”的问题补审计 §8 讨论稿。
回读 stretch/compress 的共用 Å 振幅、几何交付、MD 前门和完整 run_trial 链，
区分可执行性筛选与搜索价值判断。提出有账本／探索保底的分级投入对照，
及另立证据合同的 PES 激活／鞍点发现支路，分别说明短时／静态筛选不能推出无反应或动力学分支。
dimer／activation-relaxation 能力说明核对官方文档／原始论文；没有安装库或实施／运行新方法。

## 已完成：审计 2026-10-11 的 §1/§2/§3/§5 与 ALG-01 修复（2026-10-11）

用户 10-11 决定修这五项；§4（marker 内容 hash 重验）不修，用户判为过度工程（marker 内嵌整条 record，只有手改 marker 才会触发）。
- §1 `pars/src/pars/cli.py::_candidate_from_seed`：新增参数 candidate_id；逐帧核对元素序列、candidate/channel/provenance 与 image 序号，
  并核对文件属于 selection 点名的候选。之前链条按第 0 帧的元素重建。
- §2 `pars/src/pars/adapters/prrs_run.py`：attempts.jsonl 出现重复 attempt ID 时拒绝。PRRS 续跑会先把账本截到快照点，正常运行不会产生重复。
- §3 `mature/src/mature/io.py::OutputStore.commit`：提交前核对任务在本批 tasks.jsonl 中、seed/batch_role 与分配一致、
  protocol hash 与 manifest 一致、实现 hash 与 manifest 一致（manifest 记录了的话）。只在提交侧核对，读取 marker 侧不加（理由同 §4）。
- ALG-01 `pars/src/pars/frontier/objectives.py`：触发严重度按族/单位分键 `trigger_severity[<family>/<unit>]`，
  不同族的候选可用性掩码不同，支配/cover 中互不可比（保守保留）；缺族或单位时报 None。选择记录里的键名随之改变。
- §5：`docs/B09_GAP_ANALYSIS_2026-10-10.md` §1 去掉"最有利/下界"；本文件 10-10 重放条目改为"已收敛 25 个无漏收，t000327 未决"。
测试：mature 58 passed（新增 5）；PARS 219 passed（新增反例：换元素帧、乱序帧、错配代表文件、重复 ID、跨族严重度）；ruff format/check 通过。
审计的 `reproduce.py` 断言"缺陷仍在"，且调用 `_candidate_from_seed` 旧签名，修后不能再用来复现，以新增回归测试为准。

## 已完成：PRRS 两处小修（2026-10-11；审计 §6 / 已知待修项）

- `src/prrs/runner.py:2279`：`relays` 是列表，原来传给 `%d`，日志打到这里时会在 stderr 输出 "Logging error"（不抛异常、不影响结果）；改为 `len(relays)`。
- `src/prrs/perturbations.py::propose` 注 [4]：有软模谱时，原来用 `-spectral_overlap` 整列替换结构档位（genuine 扭转 0 / 对称 rotor 1），
  注 [1] 里"rotor 排后"的意图就失效了。现在排序键改为 `(结构档位, -overlap)`，谱分数只在同一档位内部排序。
  乙醇反例：软模正好是甲基转子的运动时，转子重叠 1.0、genuine 扭转 0.017；修改前转子排第一，修改后 genuine 排前面。
  这是排序策略的改动，只影响带谱、并且同时有 genuine 和 rotor 的扭转族。没有测过它对历史 run 的预算实际有多大影响。
- **PRRS 实现指纹因此改变**：改动前开始的 PRRS run 不能跨这次改动续跑，这是原有续跑门的设计。
- 测试：新增 `tests/test_spectral_rank.py::test_the_spectrum_orders_inside_the_rotor_rank_not_over_it`；
  与提案/谱/下降相关的 12 个测试文件共 210 passed（54.9 s，本机）。ruff 通过。
  **补验（同日 10:00–10:11，本机，含 mature 显式候选清单之后的源码）**：根目录 `make check` 全部四步通过——fmt-check、lint、`594 passed, 1 skipped in 643.18s`、digest `332 files, 2136 frames, 4 event pairs` 不变（另 254 次迁移比较通过，1 分 45 秒）。日志在 `benchmarks/prrs_check/results/2026-10-11/`。

## 已完成：mature M2–M4（2026-10-11；用户指示施工）

依据：`docs/PLAN_GEN_COMPAS_MATURATION_2026-10-03.md` §3、§4 M2–M4、§7（已全文读完）。代码在 `mature/src/mature/maturation.py`、`cli.py`；schema 小改。
- **M2 prepare**：PARS selection 只用于选路径（`refinement_order`，每通道最多 `max_paths_per_channel` 条，不同反应管优先）；
  候选几何取自源 PRRS `response.traj` 的自由段（`free_start`/`free`），按弧长 1/6…5/6。
  没有使用 PARS 代表文件，因为那是重采样、刚体对齐后的链，从里面取点等于把插值冒充成 MD 帧。
  插值候选记为 `path_interpolation`，带来源帧和权重；`allow_interpolation=false` 时取最近的真实帧（`observed_frame`）。
  区域两侧的结构文件 hash 和编码出的键图必须与协议一致，否则拒绝。账本列出每个通道（`not_this_region` / `sampled`）、未采样的反应管、跳过的位置。
- **M3 shoot**：每个任务调用一次 `prrs.runner.shoot(shots=1)`，`deliver=False`，种子取自任务表，首次进入用区域两侧的键图；
  判归复用只读的 `Registry.identify`。A/B = 首次进入侧；入口与淬火落点不一致只计数，不改判；
  未进入时淬火落到第三盆记 `other`，否则记 `unabsorbed`；门拒、失败分开计数。力评估用计数 wrapper 统计（不改变 ASE 缓存）。
- **M4**：`mature prepare|shoot|report`，prepare/shoot 支持 `--dry-run`。续跑沿用 M1 的拒绝规则；半写的 shot 目录隔离后重做；report 每次都从 marker 重建。
  探索批的序贯停止写进协议（`shooting.interim_look`）：前 n 发的 Wilson 区间整段落在诊断带外就停；
  排除项已经多到全样本必然超过 `other_fraction_max` 也停。确认批、抽查批固定发数。
  预算在每发之前检查，力评估和墙钟按已见到的单发最大值预留。
- **schema 改动**：`ShootingRule.interim_look`（只允许探索批）、`ShotRecord.cost_free_steps`、插值候选可记录 `frames`；
  **修复 `band_class`**：区间从带内跨到带上方（例如 4A/12B，约 [0.46, 0.91]）原来判 `unresolved`，按计划 M3 应判 `transition_region_compatible`。
- 测试：mature `82 passed`（23 s），其中新增 `tests/test_maturation.py`，在 ProtonRelay 解析势上真跑 PRRS（16 trial）+ `pars select` 后端到端测试。
  覆盖：自由段取点且保留全部原子、插值权重可复核、区域键图篡改被拒、dry-run 不落盘、同几何同种子与直接调用 `prrs.runner.shoot` 结果一致（首次进入和端点都一致）、
  probe 没有被投递、9 种判归情形、停止规则、预算在发射前拦截、kill 后重做且重报不重复计数、实现指纹变化后拒绝续跑、CLI 三个动词。ruff 通过。
- **没做**：M5（P1 真实 MACE 实跑），需要先冻结协议并报计划；`runs/` 和 GPU 都没有动。

## 已完成：审计 B/C 后续清单核对与测量建议（2026-10-11）

新增 `docs/AUDIT_NEXT_STEPS_2026-10-11.md`，审计报告顶部加入链接。
按当前源码 diff 与交付记录区分已修转子排序、已实现 M2–M4 和仍待实跑的 M5；
将 ALG-08 归为工程测量前置，区分已有解析势双鞍点工程台与真实多分子机制验证。
记录零号 8 发→正式 410 发的协议冻结顺序、GPU 计划建议和完整成本要求；
补低振幅返回不能直接剪除高振幅的漏检对照、转子分账及性能候选处理。
两个延后身份问题仍按先提出后修；PERF-01/09 的闭项来自用户清单，本次未找到原测量工件。
没有实施算法改动、执行新科学实验或重跑修复回归；旧审计测试不作为当前源码全量验收。

## 已完成：mature 显式候选清单（2026-10-11；用户同意施工）

依据：`docs/PLAN_GEN_COMPAS_MATURATION_2026-10-03.md` §4 M5 表与"确认候选规则"段（已全文读完）。
缺口：零号批要固定 1 个候选，确认批、抽查批要打探索批诊断过的候选，原来每批都从 PARS 路径重新抽候选，表达不了。
- `CandidateRule.proposal_ids`：只准备列出的候选；列了产不出来的 id 在写盘前拒绝；未列出的位置在账本记 `not in CandidateRule.proposal_ids`。
- `CandidateRule.from_batch`（必须同时给 `proposal_ids`）：从那一批的 `proposals.jsonl` 取记录、按字节拷几何，`batch_role` 改成本批，
  新增 `PathProposal.origin_batch`（目录、协议 hash、角色），精化交接里也带上。拒绝：id 不在那一批、候选的区域键图与本协议不同（A/B 含义会变）、几何文件被改过。
- 两个字段都只在设置时写进协议 JSON，未设置的协议 hash 不变。
- 测试：mature `90 passed`（34 s），新增 8 个：列出一个候选只准备它、未知 id 写盘前拒绝、确认批几何与探索批逐字节相同且溯源完整并能射击和报告、
  三种拒绝（id 不在、区域换边、几何被改）、schema 三种非法清单与往返 hash。ruff 通过。

## 已完成：修复后审计清单状态同步（2026-10-11）

用户告知一批修复完成后，核对当前源码入口、显式候选清单交付记录和 PRRS 补验日志。
更新 `docs/AUDIT_NEXT_STEPS_2026-10-11.md`：显式候选／跨批几何复用标为已有，
验证记录更新为 mature 90（交付记录）、PRRS 594 passed / 1 skipped（已读日志）及历史 digest 不变。
仍保留 ALG-08 全链测量、M5 实跑、节流收益、一般协同提案和真实多通道研究门。
本次未修改算法、重跑测试或启动实验。

## 已完成：整体算法深挖 ALG-11…22（2026-10-11）

用户要求继续开挖整体算法。新增 `docs/AUDIT_ALGORITHM_DEEP_DIVE_2026-10-11.md`，
原审计与后续清单加链接。读取当前提案／响应／网络水库与 PARS 选择／种子／精修链，
新增四个独立 CPU 探针，12 个算法阶段／几何／配置合成检查成功。
新证据在 `benchmarks/project_audit/results/2026-10-11/` 的 `deep-dig`、`seed-stages-v2`、
`audit-population`、`scale-domain`；报告均绑定相关源码和脚本哈希，检查运行前后未变。
覆盖不变背景键稀释、随机流碰撞、响应奖励、IQR 池依赖与负备用尺度、leader 命名敏感性、
交付跳跃吞掉自由段内部图像、施力段事件／阶段名的种子定位、外场曲率条件、
源图返回计数与抽审人口的范围；候补水库／实际交付覆盖与结果反馈仅提出测量方向。
种子脚本初次漏填电子态，被现有导出门正确拒绝；保留失败说明与部分工件，另开 v2 验证。
四脚本 format/lint 通过。核对 AutoNEB 与 geomeTRIC 官方资料作为精修对照候选，未安装或运行。
未改生产算法、历史数据或冻结协议，未加载势能模型、运行 MD／GPU 或复跑全量回归；
未证明历史通道丢失／污染或新策略的真实科学收益。

## 已完成：算法深挖第一批（ALG-12 / ALG-17 阶段名 / ALG-14 负尺度，2026-10-11；用户指示）

依据：`docs/AUDIT_ALGORITHM_DEEP_DIVE_2026-10-11.md`，三处都先在源码核实。
- **ALG-12**：`search.py` 的方向 seed 原来只取节点 ID 倒数第二个字节、不含微态 ID，100 个节点只有 10 种 seed，同节点微态共用一种。
  新增 `_direction_seed(seed, microstate_id, n_microstates, trials_spent)`：完整微态 ID（含节点 ID）经 sha256 取 8 字节进 `SeedSequence`，其余输入不变。
  这个 seed 只决定 `propose` 档位内的打乱顺序，不决定热动量。
- **ALG-17（只修阶段名）**：`_dividing_frame` 无切换帧时的回退查 `("perturbed", "pulse")`，`run_trial` 写脉冲用的是 `"forced"`，pulse 族永远拿不到末步切向；改为 `"forced"`。
  "最早事件抢占"属于策略问题，没有改。
- **ALG-14（只修取值域）**：`ParsConfig` 拒绝非有限或 ≤0 的 `objective_scales`（IQR 为零时拿来做除数，负值会翻转好坏）。IQR 随池漂移的问题没有动。
- **PRRS 实现指纹因此改变**，改动前开始的 run 不能跨改动续跑。
- 测试：新增 `tests/test_search.py` 两个（真 pulse trial 的回退切向、800 个微态 seed 互不相同）、`pars/tests/test_config.py` 五个（四种非法尺度、正尺度照常）。
  变异检查：阶段名改回 `"pulse"` 时新测试失败；旧 seed 公式在 800 个微态上只有 10 种 seed。PARS `224 passed`，mature `90 passed`。
  根目录 `make check` 四步全过（本机）：fmt-check、lint、`596 passed, 1 skipped in 487.60s`、digest `332 files, 2136 frames, 4 event pairs` 不变（另 254 次迁移比较通过）。日志 `benchmarks/prrs_check/results/2026-10-11/make_*_alg.log`。

## 已交付文档：评估协议与发现层设计 v0.1（2026-10-13，待用户审阅）

用户授权接续 `PROJECT_DIRECTION_REVIEW_2026-10-13.md` §8 的暂停点补写设计。
新增 `EVALUATION_PROTOCOL_DESIGN_2026-10-13.md`：两项比较、全链成本/恢复、目标 PES 参照、收益与重复设计、身份/来源隔离；
新增 `DISCOVERY_LAYER_DESIGN_2026-10-13.md`：姿态/意图、有限约束松弛、疑似脊/释放、分级升级、抽查、自适应和静态工件契约。
文档索引加入口，方向记录末尾追加交付和新暂停点。预算、数值门、样本量、N3 分子等明确待审或待测。
核对两份新稿的本地链接与尾部空白，0 错误；另核对当前 D0 的 `relax_source` 与现有 pulse 表达范围，
查询官方约束优化/Dimer 文档和公开数据集原始来源，仅作为设计依据。
本次仅文档交付，未改生产源码、历史结果或冻结协议，未加载模型、执行 MD/GPU 或重跑科学回归。
两份稿尚未冻结；用户先审评估再审发现，审过后才写施工计划。


## 已完成：发现层 v1 施工 M1–M2（计数账本 + 约束松弛探针）与 b04 冒烟（2026-10-13；冒烟，未做验收）

按 [施工计划](PLAN_EXPLORE_V1_2026-10-13.md) Task 1–6 以子 agent 施工、逐任务审查。新独立包 `explore/`（单向 import prrs；PRRS 源码未动，指纹/digest 不变）。
- `explore/src/explore/ledger.py`：forward 账本（预算许可、预留、崩溃按预留上界计费；重放按 key 用栈配对 open/close，防同 key 重试“退款”）；`instrument`/`counted_factory`（`calculate` 前计费，缓存命中只计 request，后端 `get_hessian` 计 1）。
- `explore/src/explore/probe.py`：`Intent`、`constrained_scan`（单标量约束 `FixBondCombo` + ASE BFGS，maxstep 0.05；**偏离施工计划**：ASE `FixInternals` 在线性几何上报 IndexError，双原子台和 b04 共线背面进攻都会碰到，故自写投影约束。审查者另用 Ar5 三对组合临时核查：残差 ≤2.2e-16 Å、dE*/dq 与有限差分差 2e-6——该核查不在测试套件内）、每点保存无偏 PES 力；`read_scan`（`ridge_suspected` / `jump`）、`bracket_jump`、`reverse_check`；`release`（解约束→淬火→极小确认→分类，返回释放坐标、释放点无偏力、seed、移除约束清单）。
- 测试：`make -C explore check` → fmt-check、lint、**31 passed**（账本 15、探针 16）。“pytest 通过”只说明这些合成用例成立。

**b04 真实模型冒烟**（本机 181 GPU；183 的 GPU 被他人 PID 1314669 占用，未抢卡）。结果：`benchmarks/explore_smoke/results/2026-10-13/run/`，PID 3400640 已退出，约 3 分钟。
- 意图：亲核 Cl5 → C0 距离压缩，3.6 → 1.9 Å，步长 0.1 Å，18 点全部收敛（BFGS 4–11 步）。
- 成本：扫描 157 forward / 18 点 = **8.7 forward/点**，段内墙钟 19.3 s ≈ **1.07 s/点**（2080 Ti，float64）；噪声重复 18–23 forward/次；释放：淬火 56、极小确认 13；全程 288 forward（含 1 个 unattributed，即来源快照）。模型每次 factory() 重新加载，加载时间不在扫描段内。
- 能量噪声：中点（q = 2.7 Å）每个笛卡尔分量加 σ = 0.02 Å 高斯扰动后拉回同一 q，3 次收敛重复的极差 **4.7e-5 eV**，直接作 `read_scan` 的 noise_eV。只测了一个点、3 次，靠近脊处的噪声可能更大。
- 判读：`read_scan` 报 `ridge_suspected`，峰在 q = 2.3 Å（相对扫描起点 +0.543 eV；相对扫描中能量最低点 q = 3.2 Å 为 +0.589 eV）。**但能量与反力积分在 2.3→2.2 Å 一步不一致**：ΔE = −0.085 eV，梯形积分 −0.039 eV，差 0.045 eV（其余各步 ≤0.007 eV），反力斜率 −21.9 eV/Å²（其余 ±3 以内），指向这一步发生换支；扫描未保存几何，无法确证。按设计 §4.3，跳支旁的峰不满足连续脊触发，所以**不能说已显示连续的脊**。`read_scan` 目前没有连续性检查。另记的 `jump`（2.2→2.1 Å）只是键图阈值越过，所在这一步本身是平滑的；真正的不连续（2.3→2.2）没有键图变化，`jump` 没抓到。`reverse_check` 没有运行。
- 释放（q = 2.2 Å 点）：`c1_candidate`，label `reactive`，断 C0–Cl1、成 C0–Cl5，即 SN2 产物，极小已确认。
- 事后对照（仅诊断，未进扫描）：b04 已认证 POLAR-1 势垒 0.5967–0.5988 eV（PROJECT_OVERVIEW §6，以精修的反应物端点为参照），本扫描峰相对扫描中最低约束点为 0.589 eV——两者参照点不同，且峰旁疑似换支，只能算量级相符。施工计划 Step 6 要求的“与已知 TS 的距离对照”没有做：脚本未保存几何。未做 C2。

结论：**冒烟跑通，未做验收**。约束松弛扫描在 b04 真实 PES 上全程收敛，越过高能区后释放到了正确的 SN2 产物（极小已确认）；单点约 9 forward、约 1 s。**“出现连续的脊”没有被显示**：峰旁疑似换支，判读缺连续性检查。最终审查列出 M4/M5 前要修的代码项：I1 探针函数透传 attempt/reserve/artifact；I2 缓存命中真正计数；I3 失败码和被拒几何写进账本；I4 多对意图的可行性检查（在 M3 计划中定域）。另需在 M5 前给 `read_scan` 加连续性检查（涉及新阈值，按全局规则 §8 先报用户）。
