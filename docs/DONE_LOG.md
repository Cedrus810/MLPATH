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
  Wilson(8,16) ≈ [0.28, 0.72]（可行性评估修正第 4 条）。
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

## 已完成：b04 约束扫描“疑似换支”复查（2026-10-18；冒烟，183 GPU）

起因：10-13 冒烟中 2.3→2.2 Å 一步能量与反力梯形积分差 0.045 eV，曾怀疑换支（上节）。
脚本 `benchmarks/explore_smoke/b04_branch_check.py`（新脚本，10-13 脚本与结果未动）；结果 `benchmarks/explore_smoke/results/2026-10-18/run/`
（`branch_check.json` + 四段扫描的逐点几何 `*.extxyz`），183 PID 994445，已退出；全程 698 forward。
同设置重做正向 3.6→1.9 Å（0.1 Å），再从 1.9 Å 反扫回 3.6 Å，再在 2.4↔2.1 Å 以 0.01 Å 正反各扫一遍。
- **不是换支**：0.01 Å 细扫每一步能量与梯形积分差 ≤ 0.0014 eV，几何连续（相邻点对齐 RMSD ≤ 0.048 Å）；
  正反细扫同 q 能量差 ≤ 5e-5 eV，粗扫正反也重合，没有滞回。0.045 eV 是 0.1 Å 步长跨过 dE/dq 急剧变号（−0.70 → +1.49 eV/Å）的离散误差：
  峰附近 q 2.29→2.27 Å 间离去基 C0–Cl1 从 2.21 拉到 2.44 Å，变化陡但连续。
- 峰位细化到 q ≈ 2.28 Å，相对扫描最低点（q = 3.2 Å）+0.597 eV。与已认证 POLAR-1 势垒 0.5967–0.5988 eV 数值一致，
  但参照点不同（约束扫描最低点 vs 精修反应物端点），也没有拿峰几何与已认证 TS 对照，**不能据此说峰就是该 TS**。
- 对 `read_scan` 连续性检查的含义：单看“能量–反力积分差”会在陡峭但连续的脊上误报；本例中缩小步长后差值从 0.045 降到 ~1e-3 eV。
  怎么改判读（例如差值大时加密复扫）属于新判据，按全局规则 §8 待用户决定，本次未改代码。
- **10-22 补：峰几何与已认证 TS 对照**（只读，未加载模型；`benchmarks/explore_smoke/compare_ts.py`，在 183 上读，
  结果 `benchmarks/explore_smoke/results/2026-10-22/ts_compare.json`）。参照 `runs/b04_probe2_openceiling/ts_candidates/ts0000.extxyz`
  （POLAR-1 面上 Hessian order 1，经 QM 精化与 IRC 确认，见 `docs/experiments/b04_qm/README.md`）。
  细扫 q = 2.28 Å 点（正反两遍）：与 ts0000 对齐 RMSD 0.0038 / 0.0045 Å（两 Cl 标号取小者），r(C–Cl) 2.280 / 2.293 Å 对 TS 的 2.2805 / 2.2827 Å，
  能量比 TS 低 0.10 meV（约束极小不高于同 q 上任何结构，方向相符），dE/dq = +0.08 eV/Å；相邻 0.01 Å 的点 RMSD 已到 0.023–0.052 Å。
  即 b04 上约束扫描的峰点与已认证 TS 几何重合。**这不是 C2 认证**：扫描点没有做 Hessian、两侧下降。

## 已完成：发现层 v1 M3 姿态与意图（2026-10-22；纯几何，未加载模型）

计划 [PLAN_EXPLORE_M3_2026-10-22](PLAN_EXPLORE_M3_2026-10-22.md)（用户未逐条审，指示“先推主线”后按 Subagent-driven 执行；设计稿 §1.1/§2/§3/§9/§10 已全文读完）。每个 Task 实现→审查→修复→复审，最后整体终审（opus）+ 一轮修复 + 复审。
- `explore/src/explore/streams.py`：按键随机流（追加样本不挪动已有样本）。
- `explore/src/explore/poses.py`：片段（含 ordinal）、扭转网格构象（genuine 扭转 60° 格，转子与环不枚举并上报）、构象资格（淬火 + 极小确认 + PRRS 盆规则去重且要求键图相同；片段电荷/多重度必须由调用方显式给，缺失则 `electronic_state_missing`、0 forward；失败保留几何与诊断）、`conformers_from_qualification`（姿态用已确认极小）、双片段刚体摆放（最近接触比值 = 壳层，Lipschitz 安全步长、循环有界、失败为 `placement_failed`）、均匀方向与 Haar 朝向采样、单分子姿态、L0 JSONL（设计 §2.3 字段 + 片段原子表）。
- `explore/src/explore/intents.py`：距离原语（键两向；非键按 `pair_cutoff_A` 分近/远层；二三跳内非键对按 PRRS 规则排除并计数；远层按颜色类整类抽样）、二元组合（系数 (1,1)/(1,−1)，三种局部规则 shared_atom ⊂ adjacent ⊂ all）、`IntentSpec.deltas()` 是取扫描符号的唯一入口。
- `explore/src/explore/probe.py`：I4 可行域——同号组合的目标低于重叠下限之和时，求值前判 `infeasible`（下限即现有 `min_distance_ratio·ΣR`，不新设阈值）。
- 测试：`make -C explore check` → fmt-check、lint、**80 passed**（本机，1.55 s，我亲自跑过）。“pytest 通过”只说明合成用例成立。
- 施工中的裁决（R1–R11，详见 `.scratch/sdd/explore-m3/progress.md`）：R1 原碰撞测试不可能通过，改为 monkeypatch；R2 去重加键图相同；R4 相同片段的构象 ID 加 ordinal；R5 片段电子态显式传入（原实现会把体系总电荷 −1 套到 b04 的中性 CH3Cl 上）；R6 计划里 `place()` 的外界写错（重原子时姿态不在壳层却标 placed），已改并加反例测试；R7 I4 让 M2 的重叠门测试改用 monkeypatch 覆盖；R8 二元意图 ID 改为显式传 pose_id（原按 `|` 切分会撞 ID）；R9 姿态从已确认极小摆；R10 符号统一由 `IntentSpec.deltas` 给；R11 终审小项一并修。

**纯几何计数表**（`benchmarks/explore_counts/count_l0.py`，本机 CPU 12 s；结果 `benchmarks/explore_counts/results/2026-10-22/COUNTS.md`）：

| 体系 | 原子 | 实测接触比 f | 构象配对 | 单原语/姿态 | +二元(adjacent)/姿态 | 32 个方向样本下 L1 外推 GPU 小时（单原语 / +二元） |
| --- | --- | --- | --- | --- | --- | --- |
| b04 | 6 | 2.02 | 1 | 11 | 66 | 1.9 / 11 |
| b05 | 16 | 2.40 | 1（转子未枚举，下界） | 91 | 2 024 | 15 / 346 |
| N1 | 19 | （3.67 是脚本 8 Å 摆放的假值） | 6 | 43 | 369 | 44 / 379 |
| b09 | 29 | 2.06 | 210 | 185 | 6 321 | 6 660 / 227 000 |

外推按 b04（6 原子）实测：每意图 18 点 × 8.7 forward、1.07 s/点，大体系每点更慢，只作量级；二元意图按单坐标扫描计价（可能偏低）；每姿态意图数只在网格 0 构象配对上算，乘以全部配对。b09 的 210 个构象配对里有对称副本（两条等价 OMe），既有高估也有低估。

**结论：在 b04/b05/N1 上 v1 规模可行；b09 全覆盖不可行**，量级来自“构象配对 × 方向样本 × 意图”的乘积，L1 必须抽样而不是全覆盖。待用户定：U1 构象网格（终审提醒：PRRS 把 symmetry_order>1 的扭转都算转子，所以二阶转子的非等价构象不枚举，链中双键若成桥会被网格扭转成 E/Z；开发集已核实不受影响）、U3 壳层（还须大于 `bond_scale` 1.2，否则片段间会凭空出现“键”）、U4 局部规则、U5 `n_far`，以及 b09 的抽样方式。

M4/M5 前的注意事项（终审）：`_feasible` 只要任一 delta 越界就拒整条意图，大压缩幅度会让短键（C–H 下限约 0.48 Å）的压缩意图整批 infeasible，是否按意图截断扫描域先报用户；`qualify_conformers` 须传计数 factory；Pose 不带电子态，跑探针前要写入 `info`；意图未做对称商（甲基三根 C–H 是三条）。
- **用户 2026-10-22 决定**：同意“b09 全覆盖不可行、L1 必须在预算内抽样”（原话“的确 09 太慢了……同意吧”），指示继续往下推（M4）。U1/U3/U4/U5 的具体取值尚未给出，M5 编排前再定。

## 已完成：发现层 v1 M4 C2 认证（2026-10-22；b04 真实模型冒烟通过，未做验收）

计划 [PLAN_EXPLORE_M4_2026-10-22](PLAN_EXPLORE_M4_2026-10-22.md)（用户 10-22 指示“继续往下”；评估协议 §1.1/§2/§3.2、发现层设计 §1/§5/§9 已全文读完），Subagent-driven，每 Task 审查 + 修复 + 复审，裁决 M4-R1…R6 记在 `.scratch/sdd/explore-m4/progress.md`。
- `explore/src/explore/certify.py`（新）：`certify()` 依次 `saddle_refine`（`follow_min_mode`）→ `hessian`（独立 `confirm_minimum`，恰一个低于 `minimum_check_eigenvalue_tol` 的负模，取它的不稳定模）→ `descent`（两侧 `descend_saddle`，PRRS `side_seed`）→ 身份（`same_basin`、`classify`、`chemical_key`）→ 判重（端点 key 无序对相同、能量差 ≤ `basin_energy_eV`、自同构下鞍点 RMSD ≤ `basin_rmsd_A`，只和已有 c2 比）。状态分开：`c2 / duplicate_c2 / unresolved:{saddle_not_found, hessian_<原因>, no_hessian, descent_failed, relayed, both_sides_same_basin, identity} / rejected:order_<n> / gate_rejected:<code> / budget_exhausted`。接力下降不算直接连接；自同构群溢出不拿来判重（记 identity 未决）。不新设阈值。
- 已知 v1 限制（写在 docstring）：判重没比负模方向（协议 §3.2 第 5 步列了它，缺一项只会让合并更容易）；两端极小确认在 PRRS `descend_saddle` 内部做，费用记在 `descent` 段；`hessian_source="analytic"` 时 PRRS runner.py ~684 会吞掉预算耗尽，可能报成 `gate_rejected:analytic_hessian_unavailable`；判重要等两侧下降做完才能判，重复项要付全额认证费。
- 账本/探针遗留项：I1 `constrained_scan/release/reverse_check` 透传 attempt/reserve/artifact；I2 段内记 `cache_hits`——**注意它只是“物理 calculator 的属性重读”，在 GuardedCalculator 下约等于 forward 数，量不到结构层面的缓存复用**，所以协议 §2.1 的缓存命中字段 v1 实际没有测到（PRRS 原语内部自己建 GuardedCalculator，不改 PRRS 就看不到）；I3 门拒/预算耗尽的失败码写进账本段 outcome，被拒几何随扫描结果返回。`reserve` 的语义是“预先扣住、先从这里扣”的下限，不是每段上限。
- 测试势：`explore/tests/twosaddle.py`，与冻结双鞍点台同参数，但改成旋转平移不变（H 的坐标取在三个锚原子的分子坐标系里）。**与上位计划不同**：冻结台是外场势，PRRS 的 Hessian 总投影掉刚体模、没有开关，不能直接用；冻结台原件未动。两个极小在对称置换下是同一结构，两个鞍点不同（自同构下 RMSD 0.33 Å > 0.15 Å）。
- 测试：`make -C explore check` → fmt、lint、**111 passed**（本机，我亲自跑，含终审修复）。双鞍点势上两个鞍点都 c2、互不判重（**更正**：起初两鞍点能量差 0.0745 eV 已超过 `basin_energy_eV` 0.05，是能量而不是几何把它们分开，我在裁决 M4-R6 里写的“几何余量 0.33 vs 0.15 Å”当时并没有被测到；终审后补了把能量门放宽到 0.2 eV 的测试，并做过去掉 RMSD 条件就失败的变异检查），同一鞍点第二次为 duplicate_c2；ProtonRelay 质子跃迁 c2 且有断键成键；从极小出发得 `unresolved:saddle_wrong_index`（49 forward）。
- 终审后补修：精修失败按原因分状态 `unresolved:saddle_<原因>`；`release()` 的极小确认遇到门拒记 `gate_rejected`；认证前去掉初猜上的约束；记录初猜到鞍点的位移；从 JSON 重载的已认证项仍能判重；`b04_certify.py` 的 certify.json 加入有效 SearchConfig、seed、预算、模型路径与类型、设备/精度、主机（这些字段是冒烟跑完后才加的，GPU 上没有跑过新版脚本）。
- **待用户裁决（M5 计分前）**：判“两侧同一盆”目前按固定原子标号（PRRS `same_basin`），而 PRRS 网络登记按对称商后的微态；于是“两端是置换镜像、无键变化”的鞍点会记成 c2（label rearrangement），PRRS 自己则记为构象转换。b04 这种有键变化的退化反应不受影响。

**b04 真实模型冒烟**（183 GPU，PID 1032893 已退出，墙钟约 7 分钟；脚本 `benchmarks/explore_smoke/b04_certify.py`，结果 `benchmarks/explore_smoke/results/2026-10-22/certify_b04/`，在 183 上读）：初猜取 10-18 细扫 q = 2.28 Å 点。
- 状态 **c2**：Hessian 恰一负模，虚频 −535.63 cm⁻¹（ts0000 记录 −535.65）；鞍点与 `runs/b04_probe2_openceiling/ts_candidates/ts0000.extxyz` 对齐 RMSD **5.0e-5 Å**，能量差 −6.5e-8 eV；两侧都是确认极小、无接力；事件 reactive，断 C0–Cl5、成 C0–Cl1（另一侧反之）；两侧势垒 0.5985 / 0.5964 eV，落在已认证 0.5967–0.5988 eV 范围内。
- 成本：200 forward（prepare 1、saddle_refine 27、hessian 30、descent 142）。墙钟 prepare 段 **296 s 只做了 1 个 forward**（冷启动加载模型 + CUDA 初始化；当时 183 有别人的 4 个作业占约 30 核），其余三段合计 99 s。**加载开销远大于认证本身**，跟“主线慢”的感受直接相关；是否每次 `factory()` 都重新加载、能否复用，需要先测再说，本次没测。
- 这是单通道、已知答案的控制，不说明多通道体系上的表现；也不是 QM 核验。

## 已完成：模型常驻（复用 calculator）的收益实测（2026-10-22；183 GPU，冒烟）

脚本 `benchmarks/explore_smoke/factory_reuse.py`，结果 `benchmarks/explore_smoke/results/2026-10-22/factory_reuse/factory_reuse.json`（在 183 上读；PID 1033222 已退出）。同一进程里先连续建 3 次 POLAR-1 calculator 计时，再把同一个 b04 C2 认证跑两遍：每次 `factory()` 都新建（现状） vs 只建一次、每次返回同一实例。
- 加载：第一次 **30.5 s**（冷启动），之后 1.4 s、0.8 s。
- 一次认证调用 `factory()` **13 次**；墙钟 现状 39.5 s → 复用 31.3 s，**省 8.2 s（约 21%）**，与 12 次热加载 × 0.6–0.8 s 相符；两遍结果完全一致（c2，鞍点能量相同），forward 数都是 199。复用后约 0.16 s/forward，与 10-13 扫描的约 0.12 s/forward 同量级。
- **更正**：我此前根据“认证约 0.5 s/forward vs 扫描 0.12 s/forward”推测重复加载让认证慢了约 8 倍，这个推测不对。同一份认证今天只用 39.5 s，而下午那次用了 99 s（另加 296 s 冷启动），差别主要来自运行环境（冷启动、当时 183 上别的作业），不是重复加载。
- 结论：进程内复用 calculator 值得做（省约 20%，风险小）；更大的一项是**冷启动**（这次 30 s，下午 296 s），所以 M5 应让一个长寿命进程在一张 GPU 上连续处理多个候选，而不是每个候选起一个新进程。复用后时间主要花在 forward 本身；再往下要靠别的手段（例如模型本是 float32、现在转成 float64 在跑），那需要先评估精度对 Hessian 的影响，本次没做。
- 用户 2026-10-22 同意：`explore/src/explore/ledger.py` 新增 `reused_counted_factory`（calculator 只建一次，仍逐次计 forward；每次调用清掉缓存结果，因为 ASE 判断是否重算时不看 `atoms.info` 里的电荷/自旋），2 个新测试，`make -C explore check` 113 passed（本机）。M5 起“一个进程守一张 GPU”写进 `PLAN_EXPLORE_V1` 的执行约束。现有脚本与 certify 默认仍用 `counted_factory`，没有改调用方。

## 已完成：POLAR-1 float32 vs float64 实测（2026-10-22；183 GPU，冒烟，各跑一次）

脚本 `benchmarks/explore_smoke/fp32_vs_fp64.py`（每种精度一个进程，因 MACE 会改 torch 全局默认 dtype），结果 `benchmarks/explore_smoke/results/2026-10-22/fp32_vs_fp64/compare.json`（在 183 上读；链 PID 1033436 已退出）。
- **单次 forward 中位数**（各 20 次）：b04（6 原子）0.278 → 0.111 s（**2.5×**）；b09（29 原子）0.266 → 0.191 s（**1.4×**）。只测了一轮，当时 183 上有别人的作业占约 30 核。
- **能量**：fp32 的能量就是 float32 精度，−26132 eV 处最小可分辨约 0.002 eV。相对能量（同一组几何间的差）最大误差 b04 **2.3e-3 eV**（50 个几何）、b09 **2.8e-3 eV**（11 个），是 fp64 扫描实测噪声 4.7e-5 eV 的约 50 倍。
- **力**：与 fp64 的差中位数 ~1e-4 eV/Å，最大 b04 5.4e-3、b09 1.7e-3 eV/Å（现行 `quench_fmax` 0.01 eV/Å 的约一半）；fp32 同一几何两次求值的力差 1.9e-4 eV/Å（GPU 非确定性），fp64 为 2e-13。
- **b04 C2 认证**：fp32 仍是 **c2**；虚频 −534.81 vs −535.63 cm⁻¹，最低本征值 −1.0518 vs −1.0551；鞍点与 ts0000 RMSD 2.5e-4 vs 5.0e-5 Å；势垒（鞍点 − 一侧极小）0.5977 vs 0.5984 eV。forward 数 **235 vs 199（+18%）**，墙钟 40.9 vs 73.4 s（同一份 fp64 认证上一小时测的是 39.5 s，墙钟受机器负载影响大，以 forward 数为准）。
- 加载 45.6 s（fp64，先跑）vs 5.6 s（fp32，后跑）：先后顺序造成的冷热缓存差，不是精度的效果。
- 结论（就 b04/b09 这两组）：fp32 不影响 b04 上 C2 的结论与鞍点几何，单次 forward 快 1.4–2.5 倍，但收敛要多花约 18% forward；能量分辨率降到约 2–3 meV、力的偏差可到 fmax 的一半。凡是依赖细能量差的判断（`read_scan` 的噪声门、0.01 Å 细扫）在 fp32 下要重新测噪声。混合精度（搜索段 fp32、认证段 fp64）是否值得、以及在哪些阶段切换，待用户定；本次未改任何默认。
- 用户 2026-10-22 同意：M5 起搜索段 float32、认证段 float64，已写进 `PLAN_EXPLORE_V1` 的执行约束，并要求 M5 计划先在 float32 下重测扫描噪声。

## 已完成：软模结构上 float32 的 Hessian 判读（2026-10-22；183 GPU，冒烟）

脚本 `benchmarks/explore_smoke/fp32_hessian_softmodes.py`，结果 `benchmarks/explore_smoke/results/2026-10-22/fp32_hessian/compare.json`（在 183 上读；链 PID 1033657 已退出，墙钟约 10 分钟，比我事先估的“几分钟”长）。7 个结构（10-08 在 POLAR-1 上 fp64 重优化的 b09 TS-2OMe、TS-2OMe-cis、二甲氧基卡宾、环己烯甲醛，b09 复合物极小，b08 TS-CCl₂-C4、TS-CCl₂），每个在 fp64 和 fp32 下各算一次 PRRS `hessian_spectrum` 和 `confirm_minimum`。
- **符号没有翻转**：7 个结构里，|本征值| 超过 `minimum_check_eigenvalue_tol`（1e-3）的模式，fp32 和 fp64 符号全部一致，负模个数全部相同；最低 12 个本征值的最大差 4.3e-4（最大出现在 TS-2OMe，已接近容差的一半）。虚频差 ≤1.4 cm⁻¹。
- **但有 2/7 在 fp32 下过不了 PRRS 的数值自检**：TS-2OMe-cis 报 `trivial_mode_floor_above_tolerance`（刚体模的数值下限超标），TS-CCl₂-C4 报 `internal_mode_recheck_disagrees`（第二种方法复核不一致）；fp64 下两者都正常判定。也就是说 fp32 没有给出错的阶数，而是让判定变成“无法确定”，在 C2 里会落到 `unresolved:hessian_*`。
- 附带发现（与精度无关，两种精度一致）：10-08 重优化的“二甲氧基卡宾”几何有一个 −81 cm⁻¹ 虚频，不是极小；b08 的 TS-CCl₂-C4 在这个几何上 fp64 判为极小（10-08 表里记的是 −4 cm⁻¹ 的近平坦鞍点）。这两条只是记录，没有追查。
- 结论：认证段继续用 fp64（与 10-22 的决定一致）；若想认证也用 fp32，需要按 fp32 实测噪声重新推导 PRRS 那两个自检的容差，属于改门，按全局规则 §8 要先报用户。

## 已交付：约束扫描性能挖掘与有界冒烟（2026-10-23；比较未完整完成）

用户明确关注当前扫描速度，误发的“节省 VRAM”回复撤销。
报告 [SCAN_PERFORMANCE_2026-10-23](SCAN_PERFORMANCE_2026-10-23.md)，测量脚本
`benchmarks/explore_smoke/scan_performance.py`；结果在
`benchmarks/explore_smoke/results/2026-10-23/scan_performance/`，直接在结果生成节点 183 读取。

- 已有结果核对：calculator 复用省约 21% 测的是认证，不能套给单条扫描；
  搜索 fp32/认证 fp64 与常驻进程已由 10-22 决定，M5 尚待接入。
- 新测：benchmark 内关闭 MACE 默认计算的 stress，fp32 固定几何中位求值耗时
  b04 0.07905→0.07521 s（1.051×），b09 0.07670→0.07435 s（1.032×）。
  小样本有慢调用，未证明稳定收益；b09 能量/力差最高 1.953e-3 eV / 1.609e-3 eV/Å，
  未做同后端重复噪声控制，未宣布等价，也未接入生产。
- b04 stress-on 基线 18/18 点收敛、158 forward、15.024 s（probe 段 13.451 s），
  每点 BFGS 平均 7.78 步；源码每点重建 BFGS 丢失曲率信息，是下一项待测候选，未改优化器。
- 两次冒烟均限 180 s：首次停在加载，末次进程状态为 NFS RPC 等待，未捕获退出码；
  同协议重试加载 105.961 s，完成固定几何与一条扫描后超时，退出码 124。
  第二条扫描已启动未闭合（forward 数未知、预算上界 400），另两条未启动，ABBA 不完整。
  全部保留，未第三次重试；现场确认本轮所有 PID 已退出。
- 实测脚本归档并核对来源指纹；之后改善脚本的启动日志/周期栈回溯和段 reservation，
  只做静态检查，未重新 GPU 验证。ruff 格式/lint、语法检查通过；未改生产源码、环境或阈值。

## 已完成：相邻点曲率复用与常驻模型四路扫描（2026-10-23；b04 工程冒烟）

用户指示修复每点曲率丢失，并允许模型常驻、多个小体系同时推进，检查控制器耗时。
实现及边界见 [SCAN_THROUGHPUT_2026-10-23](SCAN_THROUGHPUT_2026-10-23.md)。

- `explore/probe.py`：默认传递前一点已收敛的曲率并投影到新切空间，重置 secant 历史；
  `ScanSettings.reuse_curvature=False` 可保留旧基线，反扫透传该设置。
- `explore/parallel.py`：独立扫描控制器、一个常驻模型主线程、显式四路合批；
  各任务独立缓存/约束/预算，在物理 dispatch 前逐结构计费，异常通过 Future 完成回调回收。
- 真实冒烟首轮在新增剖析路径失败：Python 的多个 profiler enable 冲突，worker 在 try 前
  退出却未回收 active，造成空闲 GPU 等待。最小 CPU 复现确认；失败日志和初版源码保留。
  按确切 PID 停止后，最终退出码 137；没有完成扫描组，不算性能结果。
- 修复后 183 2080 Ti、fp32、4 CPU 线程，单进程一次模型加载，全部完成且退出码 0。
  同一起点每点重置 **158 forward / 17.482 s** → 曲率复用 **84 / 9.527 s**（probe 段）。
  同一组四条扫描串行 **37.161 s / 328 forward** → 合批 **16.600 s / 323 forward**，
  吞吐 **2.239×**，按实际结构数归一化 **2.205×**；批次实际模型调用 83 次，仍付 323 单位。
- 9 条扫描 **162/162 点收敛**，逐点键图一致；合批/串行最大能量差 3.906 meV、
  原子位移 0.00861 Å、反力差 0.00942 eV/Å，未宣布逐位相同或一般科学等价。
- 控制器 evaluator 外墙钟从基线 1.365/17.537 s（约 8%）到合批 3.323/16.600 s（约 20%）；
  不是全部耗时都来自控制器。导入/加载分别 30.160/69.078 s，常驻模型摊销这些一次性成本。
- 结果 `benchmarks/explore_smoke/results/2026-10-23/scan_throughput/SUMMARY.json`，
  原始结果在 `fp32_fixed/`；实测源码归档并逐个匹配 manifest，现场确认所有 PID 已退出、GPU 释放。
  相关 26 项回归通过；最终反扫透传后的完整包 **118 passed in 15.35 s**，格式/lint 通过，
  日志 `scan_throughput_check_final.log`。
- 不改 PRRS 源码、环境、搜索/认证精度决策或科学阈值。执行原语可用，M5 编排、全局预算恢复、
  多通道收益和长期负载尚未验证；本次不是发现层科学验收。

## 已完成：发现层 v1 M5 D2 编排（2026-10-23；解析台 kill→resume 通过 + b04 真实模型冒烟，未做验收）

计划 [PLAN_EXPLORE_M5_2026-10-23](PLAN_EXPLORE_M5_2026-10-23.md)（用户 10-23“可以 你执行完了”：按默认 M5-D1…D6 施工）；发现层设计全文、评估协议 §2/§4.1 已读完。Subagent-driven，每 Task 实现→审查（含修复轮）→终审（opus）→一轮修复→复审；裁决 M5-R1…R11 记在 `.scratch/sdd/explore-m5/progress.md`。
- `explore/src/explore/ledger.py`：`Ledger(replay=False)`（多个按任务账本共用一个文件）、写文件加模块锁、`tally()`（按 task/attempt 汇总已关闭 forward 与未关闭段，忽略被 kill 截断的末行）、`CountedView` / `view_factory`（一个常驻 calculator、每个任务各自计数，先计费后求值，每个视图独立缓存）。
- `explore/src/explore/journal.py`（新）：运行日志（每条带 host、t）、搜索额 / 认证额两账户、`repair_tail`。
- `explore/src/explore/probe.py`：`feasible_prefix`（M5-D2：幅度格越过重叠下限时截断定义域，剩不足 2 点才 infeasible）、`reverse_plan`（`reverse_check` 改用它，行为不变）。
- `explore/src/explore/d2.py`（新）：`init_run` 冻结人口表（按键随机顺序、每姿态配额轮转、冻结时抽查签、定义域截断；缺电子态在任何求值前拒绝）；`run` 严格分相 L1 → L2（跳支反扫）→ 升级（脊 > 夹住的跳支 > 抽查；L3 静态释放记搜索额，L4 认证记认证额，`audit_release` 标记）→ `report.json`；dispatch 先扣上限、settle 按实数结算；续跑时未结算的 attempt 能按账本核对就记实数，有未关闭段就按上限计入并列 unreconciled，用新 attempt 号重做；首次命中成本 = C2 结算时的累计支出。
- 测试：`make -C explore check` → fmt、lint、**149 passed**（本机，我亲自跑，66.7 s）。其中 `test_d2_resume.py` 用子进程 `os._exit(137)` 真 kill（L1 中途、认证中途各一次）后用同一条命令续跑：C2 集合、释放/认证状态与不中断跑完全相同；总支出 = 不中断跑支出 + 被打断 attempt 的计费（精确相等）；外部计数器 ≤ 总支出。“pytest 通过”只说明这些合成用例成立。
- 与设计/上位计划不同（计划里写明）：L3 不做 MD（R8 时窗未定）；§6.1 第 3 条（重复回同盆减额）不做；L2 不做峰附近幅度细化；出口在旋转不变双鞍点势上验。

**fp32 扫描噪声与同进程混合精度**（Task 7；183 GPU，PID 1037632 已退出，约 1 分钟；脚本 `benchmarks/explore_d2/fp32_noise.py`，结果 `benchmarks/explore_d2/results/2026-10-23/fp32_noise/`，在 183 上读）：
- b04 C0–Cl5 压缩扫描 3.6→1.9 Å，fp32 + 曲率复用：18/18 收敛，80 forward，9.6 s；峰 +0.543 eV（q = 2.3 Å），与 10-13 fp64 扫描的 +0.543 eV 一致。
- 噪声（q = 3.2/2.7/2.4/2.3 Å 各 5 次扰动重复，σ 0.02 Å，20/20 收敛）：极差 1.953e-3 / 0 / 0 / 0 eV。1.953e-3 正是 −2.6e4 eV 处 float32 的能量量子（2⁻⁹ eV）：弛豫噪声低于 fp32 的能量分辨率。**按 M5-D5 默认规则 noise_eV(b04, fp32) = 1.953125e-3 eV**，read_scan 要求升降大于它，即至少 2 个量子（3.9 meV）。换体系要重测。
- 同进程共存：先建 fp32 再建 fp64，fp64 与单独 fp64 进程对 3 个几何能量差 0、力差 ≤2.5e-12；fp32 在建 fp64 前后能量差 0 / 1.953e-3 / 0（1 个量子）、力差 ≤8.1e-4 eV/Å。即 10-22 担心的“MACE 改全局 dtype”在当前 MACE 版本里不成立（构图用上下文管理器）；一进程同时常驻两种精度可用。

**b04 D2 端到端冒烟**（Task 8；183 GPU；脚本 `benchmarks/explore_d2/run_b04_d2.py`，run 目录 `benchmarks/explore_d2/results/2026-10-23/b04_d2/`，在 183 上读）：2 个姿态（壳层 f = 2.0）× 17 条距离单原语 = 34 条 L1；第一个 L4 认证开始后按 PID（1037959，ps 核对）杀掉，同一条命令续跑（PID 1038323）到完成，退出码 0。
- L1：19 completed / 15 censored_energy，4786 forward，197.6 s（4 路合批 fp32 约 0.041 s/forward）；定义域被截断 8 条；无跳支，L2 0。
- 升级：脊 3 + 抽查 3；暂缓 no_trigger 28。L3 释放 6 个 c1_candidate（568 forward，fp64 约 0.12 s/forward）。L4：**c2 1**、duplicate_c2 1（另一姿态的同一鞍点）、descent_failed 1（就是被 kill 后重做的那个，重做花 785 forward）、seed_above_ceiling 2（抽查项）、not_converged 1（1692 forward，fp64 约 0.20 s/forward）。
- C2 鞍点与 `ts0000` 对齐 RMSD **3.2e-4 Å**（重复项 2.8e-4 Å），事件 reactive（C0–Cl 一断一成，两端同 chemical_key，退化自环）。首次命中成本 6750 forward（累计臂内支出）。
- 账目：搜索 5354 + 认证 2492 = 7846 = 账本实计 7046 + 被 kill attempt 按上限计入 800（该 attempt 账本可核对的 forward 为 0，按协议 §2.3 列 unreconciled）。模型加载 43 s/进程。
- **终审 Important #2 在这里实测到**：冒烟取 `certify_max_rise_eV` = 1.0（相对姿态起点），低于扫描上限 3.0，于是 2 个抽查项的初猜在第 0 步就被认证天花板拒绝（各花 8 forward）。之后的运行取 certify 上限 ≥ 扫描上限（裁决 M5-R10）。
- 这是单通道已知答案的工程冒烟，不是发现层验收；多通道收益、正式预算、U1–U5 都没定。

终审另两点已按裁决处理：同主机续跑（NFS 属性缓存 6 h + nocto，换主机续跑可能读到旧 journal）——写进 d2 与脚本文档，每条 journal 事件与 manifest 记 host，报告列 `hosts`，没有加拒跑门；journal 带完整扫描点，b04 每条 L1 settle 约 15 KB，b05 以上规模要再看（M5-R11，暂缓）。

## 已完成：b05 fp32 扫描噪声冒烟（2026-10-23；180 GPU，D2 长任务的前置）

用户 10-23 说“跑”，并说 181、183 给别人用，改在 180（kasuga，2080 Ti，同一 NFS，属性缓存同样 6 h，结果在 180 上读）。脚本 `benchmarks/explore_d2/run_b05_d2.py --measure-noise`（L0 生成抽到 `benchmarks/explore_d2/l0.py`；OH⁻ 只有 2 个原子，不做构象资格确认）。
- 第一次（PID 139762，结果 `results/2026-10-23/b05_noise/`）在噪声重复阶段崩溃：噪声模式把未设置的 `--max-active`（None）传给了 `scan_many`。4 条扫描已经做完并记进账本，但 noise.json 没写出。改为与 `fp32_noise.py` 相同的固定 4 路并发（只改这一行），重跑到 `results/2026-10-23/b05_noise_r2/`（PID 139927，退出码 0）。
- L0：壳层 f = 2.4（M3 计数表里 b05 实测接触比）、1 个姿态；准备费用 52 forward；加载 38 s。
- 4 条意图（按键随机抽取）各扫到 ±0.5 Å，每条 150–160 forward，4 条合计 14.9 s；每条取最后一个收敛点做 5 次扰动重复（σ 0.02 Å），20/20 收敛，15.4 s。fp32 4 路合批实测约 **0.024 s/forward**。
- 极差 3.906e-3 / 1.953e-3 / 3.906e-3 / 0 eV，标准差 ≤ 1.6e-3 eV。b05 能量处 float32 的能量量子也是 1.953e-3 eV，所以 b05 的弛豫噪声约 2 个量子，比 b04（1 个量子）大。**按 M5-D5 规则 noise_eV(b05, fp32) = 3.90625e-3 eV**，read_scan 要求升降至少 3 个量子（5.9 meV）。

## 已完成：b05 D2 v1 两个配对块（2026-10-23；PBS yayoi42，RTX 4080 Super，作业 18345，14:02–19:05）

用户给出空节点 yayoi42 并说“投下这个任务”；脚本 `benchmarks/explore_d2/b05_d2_yayoi42.zsh`（照抄 PBS 模板，只改 ← 行），两个块（seed 20261023 / 20261024）同卡并行，各 8 线程。参数：壳层 2.4、8 个方向样本、距离单原语全取、noise_eV 3.90625e-3（b05 实测）、搜索额 300k / 认证额 40k forward、上限 scan 800 / reverse 200 / release 400 / certify 1500、抽查 0.1、认证天花板 = 扫描天花板 3.0 eV。结果 `benchmarks/explore_d2/results/2026-10-23/b05_d2_s{A,B}/`（181 上直接读，目录此前未被 181 访问过）。两块 unreconciled 都是 0。

| | sA | sB |
| --- | --- | --- |
| 人口（ready） | 932（113 条定义域被截断） | 894 |
| L1 扫描 / 未轮到 | 863 / 69（搜索额用完） | 861 / 33 |
| L1 forward、墙钟 | 299 404，2.25 h（平均 347、中位 335 forward/条） | 299 295，2.23 h |
| 升级：脊 / 夹住的跳支 / 抽查 | 223 / 13 / 57 | 237 / 23 / 56 |
| L3 释放 | 0（293 条全因搜索额用完跳过） | 0（316 条跳过） |
| L4 认证 / 因认证额跳过 | 48 / 245（平均 803 forward/个，2.72 h，约 0.25 s/forward fp64） | 50 / 266 |
| L4 状态 | c2 2、saddle_not_converged 18、relayed 11、descent_failed 10、wrong_index 5、budget_exhausted 2 | c2 2、not_converged 25、relayed 14、descent_failed 6、wrong_index 1、budget_exhausted 2 |
| 首个 C2 时的累计支出 | 300 826 | 315 581 |

- **C2 事件**（端点键变化取自认证记录的 fingerprint）：sA 两个、sB 一个都是“一侧有 C0–Cl1 与 C(甲基)–H 键，另一侧该 H 在 O14 上、片段数 3↔2”，即**协同 E2 消除**鞍点（t-BuCl + OH⁻ ↔ 异丁烯 + H₂O + Cl⁻）。sB 另一个是无键变化的 rearrangement 鞍点，能量高 0.22 eV。
- 对称商后 RMSD：sA 第二个 E2 鞍点与 sB 的 E2 鞍点 **0.002 Å**（ΔE 0.5 meV），即两个独立块找到同一个鞍点；sA 内两个 E2 鞍点 ΔE 0.9 meV 但 RMSD 1.67 Å，按判重规则是不同鞍点，是否同一 TS 的不同构象（例如水/Cl⁻ 的相对位置）未核实。
- **取代通道没有得到 C2。** C0–O14 压缩意图两块共 6 条全部扫完：5 条无触发；sB 1 条出脊但认证未收敛；另有 3 条 C0–O14 意图被升级但因认证额用完未认证。
- **预算是瓶颈**：搜索额全部被 L1 基础覆盖用完（严格分相 M5-D6），所以一个释放都没做（C1 = 0），每块还有 33–69 条意图没扫到；脊的升级率约 26%（b04 冒烟是 9%），认证额只够约 50 / 290 个升级项，抽查项和夹住的跳支一个也没轮到认证（升级顺序排在脊之后）。
- 这不是验收：b05 没有冻结参照集 S，算不出召回率；QM 核验没做；POLAR-1 在氢氧根/消除产物上的适用性本身未验证（`benchmarks/b05_tertiary_halide/model.json`）。

## 已交付：HF TS 核心问题、缺陷与接手路线（2026-10-24；文档整理）

用户明确当前核心是面向 MLP 的 Hessian-free TS 搜索，完整路径与连接在其后，并要求把缺陷、修法、整体方向和未来工作写给接手者。

- 新增 [TS_HESSIAN_FREE_DESIGN_RECORD_2026-10-24.md](TS_HESSIAN_FREE_DESIGN_RECORD_2026-10-24.md)：方向来源、术语与边界、当前实现地图、F1–F8、连续分支 E/F 反力零点内核建议、HF-1 至 HF-6、论文验证与未来方向。
- 源码复核发现 `certify → follow_min_mode` 不仅验收，而且周期性调用完整谱并按本征值精修；准备／端点确认也存在 Hessian 依赖。记录区分算法目标与当前实现，不宣称已做到端到端 HF。
- 本轮 review 的四个 CPU 行为复现写入附录：同图几何跳支仍触发脊；未收敛尾点可提供下降证据；低于高能起点的局部峰被漏掉；精确解析 TS 在完整谱不可用时未返回几何。合成反例不证明历史分子结果已受污染。
- 更新项目总览、文档索引、旧 explore 计划的入口提示与项目上下文；修正索引中 mature M2–M4 的过期未实施描述。原计划、冻结协议、C1/C2 与旧结果保留。
- 本次只写设计／文档，不实施算法、数值阈值、预算或状态格式修改，不运行 GPU／分子实验；新 HF 路线的验收条件是待实施要求。文档附录四个 CPU 行为已用生产解释器重新复现并核对，127 个本地链接有效、代码围栏闭合；提交说明在 `.commit_drafts/2026-10-24_hf-ts-design-record.md`。顺带修正旧 Wilson 区间注释误形成链接的括号，数值未变。

## 已完成：b05 D2 两块重跑，全程一个 fp32 模型（2026-10-23→24；PBS yayoi32，作业 18465，19:46–05:35）

起因：用户 10-23 指出 POLAR-1 权重本身是 float32（已核实，见 `PLAN_EXPLORE_V1` 精度一条），fp64 只是升格。`run_b05_d2.py` / `l0.py` 改为只建一个 `polar1("cuda", "float32")`，L1 合批与释放/认证共用它（释放和认证仍串行，`run_many` 合批尚未接入 d2）。先在 181 冒烟一次（PID 4047131，`results/2026-10-23/b05_fp32_smoke181/`，1 姿态 12 意图：L1 3442 forward、1 个脊、释放 c1、认证 not_converged，端到端跑通）；PBS 脚本 `benchmarks/explore_d2/b05_fp32_pbs.zsh`（照抄模板）。参数同第一轮，只把搜索额 / 认证额改为 340k / 240k。PRRS 容差未改（Task 4 的 fp32 Hessian 噪声测量已在 180 跑完，结果 `results/2026-10-23/fp32_hessian_noise/`，尚未读出、未采用）。

| | sA | sB |
| --- | --- | --- |
| L1 扫描 | 932 条全部扫完（683 completed / 243 censored / 6 budget），2.0 h | 894 条全部扫完，2.0 h |
| 升级：脊 / 跳支 / 抽查 | 234 / 30 / 58 | 251 / 26 / 58 |
| L3 释放 | 121 做了（c1 82、not_minimum 33），201 因搜索额跳过 | 206 做了（c1 141），129 跳过 |
| L4 认证 | **322 个全部做完**，233 125 forward，7.3 h | 318 做完、17 因认证额跳过，239 403 forward，7.7 h |
| C2 | 3：E2 1（−18880.0195 eV，与第一轮 E2 差 1.5 meV）、无键变化 2 | 2（均无键变化）+ duplicate 1 |
| 主要失败 | not_converged 112、descent_failed 94、relayed 40、seed_above_ceiling 26、wrong_index 17、same_basin 8、fp32 Hessian 自检（trivial floor / recheck）10 | not_converged 114、descent_failed 78、relayed 35、seed_above_ceiling 16、wrong_index 18、same_basin 8、fp32 自检 21 |

- fp32 认证速度：saddle_refine 0.16–0.17 s/forward、descent 0.06、hessian 0.10（第一轮升格计算 0.40 / 0.14 / 0.23，约快 2.4 倍）；释放 quench 0.053、minimum_check 0.143。
- **结论**：认证额加到 6 倍、几乎所有升级项都认证了，C2 仍只有 2–3 个/块，认证成功率约 1%；E2 这次只在 sA 出现，sB 没有；取代通道仍无 C2。失败主要是精修不收敛（约 35%）和下降失败 / 接力（约 35%），不是预算或精度问题。fp32 下 PRRS 的 Hessian 数值自检拦下 9–21 个（约 3–6%），这部分按 fp32 噪声重定容差可以收回一部分，但不改变整体结论。

## 已完成：b05 旧批次中 F3 假脊的离线统计（2026-10-24；只读 journal，未加载模型）

用户 10-24 同意。依据 [HF 设计记录](TS_HESSIAN_FREE_DESIGN_RECORD_2026-10-24.md) §4 F3：`read_scan` 的峰后下降取“所有后续点”的最低能量，未收敛尾点也能提供下降证据。脚本 `benchmarks/explore_d2/f3_false_ridges.py` 对每条已结算 L1 扫描分别用现行规则和“下降只取峰后连续收敛段”的修正读法判一次，其余条件不变；修正读法先在记录附录 A 的 failed_tail 构造上核对过（旧规则触发、修正读法不触发；正常峰两者都触发）。结果 `benchmarks/explore_d2/results/2026-10-24/f3_false_ridges.json`。
- 四批（第一轮 sA/sB、fp32 重跑 sA/sB）共 945 个 `ridge_suspected`：**F3 假脊 0 个**，峰位也全部不变。
- 原因：这几批扫描几乎没有未收敛点——第一轮 sA 13 452 个扫描点里 6 个未收敛，fp32 sB 13 905 个里 10 个，各分布在 6 / 10 条扫描上。
- 含义：F3 是确实存在的实现漏洞（附录 A 可复现），但**解释不了 b05 认证的低成功率**；b05 的失败要从 F1（精修依赖完整谱）、F6（方向未传入）、F4/F5（连续性、局部峰）和跳支去找。F4 同图突变的比例本次没有统计。

## 已完成：HF-1 漏洞修复 F3 + F6，D2 规则版本 v1.1（2026-10-24；181 冒烟跑通）

用户 10-24：「先修漏洞 然后准备重跑」。依据 [HF 设计记录](TS_HESSIAN_FREE_DESIGN_RECORD_2026-10-24.md) §4、§7 HF-1。
- `explore/src/explore/probe.py` `read_scan`（F3）：峰后下降只取峰后连续收敛段的能量；未收敛点保留在扫描里，但不作下降证据。
- `explore/src/explore/d2.py`（F6）：`_promotion_points` 同时给出初猜处的扫描切向（脊：峰后点 − 峰前点；夹住的跳支：产物侧 − 反应物侧；抽查：峰两侧点之差），`_certify_task` 以 `direction=` 传给 `certify` → `follow_min_mode`，使精修追踪扫描坐标而非最软模。L4 payload 记 `"direction": "scan_tangent"`。`RULE_VERSION` 改为 `D2-v1.1`（旧批次 manifest 仍是 D2-v1）。
- 测试：F3 两条（假下降不触发、真下降仍触发）；HF-1 夹具两条，固定 F5（峰低于扫描起点不触发）与 F4（同图横向突变仍触发）的**现行**行为，不改语义；F6 一条（`certify` 收到非零 (N,3) 方向）。`make -C explore check` → fmt、lint、**158 passed**（本机 181，72 s）。F1（精修依赖完整谱）未动，属 HF-2/3。
- 181 冒烟（PID 64230，`benchmarks/explore_d2/results/2026-10-24/b05_v11_smoke181/`，1 姿态 12 意图）：端到端跑通，rule D2-v1.1；1 个脊 → 释放 c1 → 认证收到切向（`target_supplied` True，初始重叠 0.51），结果 `unresolved:saddle_lost_target_mode`（452 forward）；账目闭合，无 unreconciled。只说明代码路径通，不说明 F6 的效果。
- PBS 脚本 `benchmarks/explore_d2/b05_v11_pbs.zsh`（未提交）：与 `b05_fp32_s{A,B}` 同种子同参数，F3 在 b05 上影响为 0（见上一节），所以这一轮是 F6 单变量的 A/B。

## 已完成：b05 D2-v1.1 两块（F6 单变量 A/B）（2026-10-24→25；PBS yayoi10，RTX 3090 ×1，16 线程）

用户给出 yayoi10（“一张卡”）后提交 `benchmarks/explore_d2/b05_v11_pbs.zsh`；与 `b05_fp32_s{A,B}` 同种子同参数，唯一变化是 F6（扫描切向传入精修器；F3 在 b05 上影响为 0）。结果 `benchmarks/explore_d2/results/2026-10-24/b05_v11_s{A,B}/`，两块 unreconciled 都是 0。对比脚本 `benchmarks/explore_d2/f6_ab.py`（只读 journal，按意图配对 L4），输出 `results/2026-10-24/f6_ab.json`。

| | fp32 sA → v1.1 sA | fp32 sB → v1.1 sB |
| --- | --- | --- |
| L1 / 升级 | 基本相同（脊 234→235，跳支 30→28，抽查 58→61） | 同左量级 |
| L4 认证个数 / forward | 322 / 233 125 → 324 / **140 126** | 318 / 239 403 → 329 / **152 547** |
| 精修到一阶鞍点（c2 + duplicate_c2） | 3 + 0 → **4 + 17** | 2 + 1 → **4 + 18** |
| C2 标签 | rearrangement 2、**reactive（E2）1** → rearrangement 4 | rearrangement 2 → rearrangement 4 |
| `saddle_lost_target_mode` | 0 → **122**（中位 35 步；初始重叠中位 0.73） | 0 → **111**（中位 25 步；0.70） |
| `saddle_ceiling_blocked` | 0 → 17 | 0 → 13 |

- **能精修到鞍点的多了，但几乎都是同一个无键变化的鞍点**：sA 21 个鞍点里 15 个能量 −18880.5 eV、sB 22 个里 16 个在 −18880.50 eV（虚频约 290 cm⁻¹）；它是什么运动（例如甲基转动）没有核实。
- **反应性鞍点丢了**：fp32 sA 的 E2 鞍点（意图 `distance:10-13:+1:+`）在 v1.1 里是 `lost_target_mode`；fp32 两块的 6 个 C2/重复意图，v1.1 里全部没再得到（4 个 lost_target_mode、1 个 not_converged、1 个 relayed）。
- 失败方式换了：原来的 not_converged / descent_failed / relayed 有约 120 个/块变成 `lost_target_mode`，即跟踪扫描切向 25–35 步后没有任何本征模与之重叠 ≥ `min_mode_overlap` 0.5（`src/prrs/config.py:191`），精修中止。约 30 个/块的 lost_target_mode 记录里 `target_overlap` 为空，原因未查。
- 墙钟（L4 13.0 h / 14.7 h，fp32 时 7.3 / 7.7 h）不可比：硬件是 3090 + 两块共用 16 线程，forward 反而更少。
- **结论**：F6 单独接入，换来的是“沿距离坐标切向追踪最软模”的旧精修器（F1）更容易收敛到同一个低势垒构象鞍点，而对 E2 这种多键协同的反应坐标，单一距离切向追踪会在途中丢模而中止。C2 的去重数（4/块）略多于 fp32（2–3/块），但没有新的反应通道，反应性 C2 从 1 变 0。不支持“F6 本身能修好 b05 认证”；问题在精修器本身（F1，HF-2/3 的 E/F 定位内核）。

## 已完成：b05 丢模诊断修复与 HF-2/3 独立 E/F 定位原语（2026-10-25；解析验证，未接 D2）

实现及下一关见 [HF-2/3 记录](HF_EF_LOCATOR_2026-10-25.md)。

- 只读审计 `benchmarks/explore_d2/hf_f1_audit.py`：233 次丢模中，59 次在第 0 步退出
  （sA 30、sB 29），其初始重叠全部为空；另 174 次发生于成功初选之后。
  空值已复现为返回字典的旧 `**bounds` 覆盖新诊断值，`prrs/runner.py` 只调整展开顺序，
  不改变搜索行为或阈值；历史 journal 不修改、不补造。审计 JSON 放 `results/2026-10-25/`。
- 对上一节解释作限定：旧精修第一次对照输入扫描切向，之后对照上次选中的本征向量，
  并非始终追踪原始切向，也不一定始终是最软模。E2 原反应性命中对应的 v1.1 意图在第 80 步
  停止（初始重叠 0.7030、末次 0.01818）。尚未区分物理漂移、模混合及谱数值问题，
  也未确认重复鞍点的运动类型；F6 单独未解决问题的结论保持。
- 新增 `explore/ts_ef.py`：独立 `refine_ts_ef` / `EFResult`，重新松弛两侧括区、反力正→负
  判别、双侧延续到中点、有界二分、分支几何／反力一致性诊断，同时检查完整无偏力、横向力、
  反力与约束残差；保存括区、E/F、方向与漂移、历史、失败结构、设置和逐结构成本。
  搜索不读取完整谱或模型 Hessian，BFGS 只从 E/F 历史学习近似；没有旧精修器回退。
- `converged_candidate` 与驻点性质／连接分开记录，不改变 C1/C2 定义。
  旧 D2-v1.1 和 `certify` 保持原链。当前原语不含 D2 调度、自动落盘与 kill→resume。
  本轮未运行 GPU、真实模型或 QM；不能写成 b05 E2 已找回或端到端 HF 搜索已验收。
- 验证：解析耦合距离势覆盖完整谱禁用、完整残余力、刚体变换、反扫／负坐标、
  极小／平坦、同图双分支、噪声底、预算／步数／能量截断、可靠性门与坐标退化。
  `make -C explore check` 格式／lint／全包 **177 passed**（含新增 19 项），
  PRRS 相关回归 **96 passed**；修改前后摘要检查均为 332 文件、2136 帧、4 事件对不变。
  日志和源码指纹见 `benchmarks/explore_d2/results/2026-10-25/hf_ef_checks.json`。
  分支诊断不是连续性／横向稳定的证明，下一步以已保存 b04／b05 E2 括区做有界真实模型对照。

## 已完成：HF-2/3 E/F 定位原语的有界真实模型对照（2026-10-25；本机 2080 Ti，冒烟）

依据 [HF_EF_LOCATOR_2026-10-25.md](HF_EF_LOCATOR_2026-10-25.md) 末节“下一关”。脚本 `benchmarks/explore_d2/hf_ef_bounded.py`：只读冻结 journal（不调 `open_run`，不追加历史 run），对历史上 5 个反应性 C2 意图取扫描峰两侧点 (peak−1, peak+1) 交给 `refine_ts_ef`（一个常驻 fp32 POLAR-1，搜索上限 6000 forward；仅 `converged_candidate` 才做独立 Hessian 与两侧下降）。对照量：候选与同一意图旧 C2 鞍点的 Kabsch RMSD 与 ΔE。结果 `benchmarks/explore_d2/results/2026-10-25/hf_ef_bounded/`（PID 513716 跑完 b04 后因我脚本按后缀取意图取错姿态崩溃，修正为按 report 的 c2 任务号取后以 PID 515863 续跑；日志 `hf_ef_bounded_run1.log`、`hf_ef_bounded.log`）。

| 括区 | 状态 | forward / 墙钟 | 中点（左延续 / 右延续） |
| --- | --- | --- | --- |
| b04 SN2 `0-1` | `branch_force_disagreement` | 34 / 6 s | q 2.2814：E 都是 −26132.201（旧 TS −26132.201），反力 −0.057 / −0.065 eV/Å，RMSD 5.6e-4 Å |
| b05 第一轮 sA `6-9` | `branch_unresolved` | 187 / 21 s | q 1.389：E −18880.027 / −18881.305，反力 +0.46 / −3.69，RMSD 1.04 Å |
| b05 第一轮 sA `10-13` | `branch_unresolved` | 211 / 22 s | 同上量级（−18880.027 / −18881.301，RMSD 1.04） |
| b05 第一轮 sB `10-13` | `branch_unresolved` | 158 / 16 s | 同上 |
| b05 fp32 sA `10-13` | `branch_unresolved` | 170 / 17 s | 同上 |

- **b04 是门误拦**：两次延续几何一致（5.6e-4 Å）、反力同号（峰在左半区），但反力差 0.008 eV/Å 超过 `ts_ef.py` 的 `2 × reaction_force_tol` = 0.006，第一个中点就停。中点能量已等于旧 TS，最后括区与旧 TS 的 RMSD 0.059 Å。我的估计（未测）：横向只收敛到 0.003 eV/Å，残余几何误差经横向–纵向耦合即可造成几 meV/Å 的反力差，这个门比它所比较的量的噪声更紧（全局 CLAUDE.md §8.2）。`ts_ef.py` 是另一 AI 所写，未改。
- **b05 E2 是真实跳支**：4 个括区（不同意图、不同批次，实际是同一几何区间 q 1.289–1.489）表现一致：括区右端已在产物侧（比左端低 1.49 eV），中点处左右延续落在相距 1.04 Å、能量差 1.28 eV 的两支上；左支在中点能量 −18880.027 已接近旧 E2 鞍点（−18880.018～−18880.020），反力仍 +0.46 eV/Å。即：沿单一距离坐标的约束极小路径不连续，E2 鞍点不在任一连续支上，原语按设计返回未决、不回退旧精修。要找回 E2 需要 HF 记录 §5.2 的跳支定位（或多坐标），不是调本原语的阈值。
- 5 个括区均未出候选，独立性质/连接检查都没触发。新内核真实尺度：每次约束松弛中点 10–100 步，单括区 34–211 forward、6–22 s（含 BFGS 开销，2080 Ti）。

### 补测：b04 中点的反力抖动（2026-10-25；本机 2080 Ti，PID 519762，1182 forward，145 s）

用户 10-25 “可以”。脚本 `benchmarks/explore_d2/hf_ef_rf_jitter.py`，结果 `results/2026-10-25/hf_ef_rf_jitter.json`。量：b04 括区第一个中点 q = 2.2814 处约束极小的反力 dE/dq（与 `branch_force_disagreement` 门比较的量相同），fp32 POLAR-1，约束与 BFGS 同 `ts_ef`。

| 测法 | 结果（eV/Å） |
| --- | --- |
| A. 同一几何 10 次随机刚体变换（物理上不变） | 极差 0.048，标准差 0.016 |
| B. 同一坐标系，8 次 σ 0.01 Å 扰动后重收敛，τ = 3e-3 | 极差 0.011，标准差 0.0034（几何两两 RMSD ≤ 6e-4 Å） |
| 同上，τ = 1e-3 | 极差 0.0086，标准差 0.0031 |
| C. 左右两次延续继续收紧：τ = 1e-3 / 3e-4 / 1e-4 | 反力差 0.0032 / 0.0036 / 0.0043；1e-4 两侧 400 步都收不到（fp32 底） |

- **更正上一节我的估计**：反力差不是横向收敛残差造成的——把 τ 收紧 30 倍反力差不降，是 fp32 力本身的噪声底（同一坐标系 σ ≈ 0.003 eV/Å；换坐标系还有约 0.016 的系统偏移，MACE 在 fp32 下旋转等变只近似成立）。
- 对门的含义：同支两次延续的反力差标准差约 √2 × 0.0034 ≈ 0.005，实测同支极差 0.011，而门是 2 × 0.003 = 0.006，所以同一支也会被判分歧；b04 的 0.008 落在噪声内。真实跳支（b05）的反力差约 4.2，比噪声大三个量级。`reaction_force_tol` = 0.003 本身约等于 1σ，近零处的符号判断由噪声决定。
- A 项还意味着 `fmax_eV_A` = 0.01 在某些取向下可能被 fp32 力误差占满，未进一步测。以上只在 b04 一个点上测，b05 体系（能量处 2 个量子）需另测。`ts_ef.py` 未改。

### `branch_force_disagreement` 门改为按实测抖动缩放，b04 重跑（2026-10-25；本机，冒烟）

用户 10-25 “可以”。`explore/src/explore/ts_ef.py`：`EFSettings` 新增 `reaction_force_sigma`（调用方按体系实测的同支反力 1σ）；给定时门为 `3·√2·σ`（两次独立延续之差的 3 倍标准差），未给定时保持原 `2 × reaction_force_tol`，阈值写进 `provenance["branch_force_limit"]`。`explore/tests/test_ts_ef.py` 加一项（交替 ±2e-3 eV/Å 反力抖动：旧门判分歧，σ = 2e-3 时不判）；`test_ts_ef.py` 20 passed，ruff 通过；全包 `make check` 未跑。`hf_ef_bounded.py` 加可选参数 `SIGMA [CASE ...]`。

b04 重跑（σ = 0.0034 → 门 0.0144 eV/Å；结果 `results/2026-10-25/hf_ef_bounded_sigma/b04.json`，166 forward、17 s）：
- 二分走了 7 个中点（原来第 1 个就停），括区缩到 q ∈ [2.2782, 2.2798]，**最后括区与旧 b04 TS 的 RMSD 5.1e-4 Å**，能量相同（−26132.201 eV）；每个中点两次延续几何差 2–7e-4 Å。
- 仍以 `branch_force_disagreement` 停在第 8 个中点 q = 2.27902：左右延续反力 +0.038 / −0.0047（差 0.043，约 13σ），几何差 2.3e-3 Å；两次延续都只走了 1–2 步 BFGS 就满足横向 0.003。原因我不清楚（猜测是近根处软横向模的残余位移经耦合进入反力，未测）。
- 前一个中点 q = 2.2798 的右延续已是完整力 fmax 0.0058 eV/Å（≤ `fmax_eV_A` 0.01），但 |RF| = 0.0053 > `reaction_force_tol` 0.003，且要求左右两次都满足，所以不出候选。`reaction_force_tol` ≈ 1σ，这一判据在 fp32 下基本达不到——这是设计层面的问题（候选判据），没有再叠补丁。

### 已修复：HF E/F 候选判据分离分支一致性与驻点收敛（2026-10-25）

- `ts_ef.py` 改为 `HF-EF-bracket-v2`：保留双侧几何、步间几何和反力差检查；通过后，
  任一延续结构满足完整无偏力、横向力和约束残差即可返回候选，选择达标者中完整力最小的。
  移除候选的独立 `|RF| <= reaction_force_tol` 门与双侧同时驻点收敛要求。
  `reaction_force_tol` 仍用于括区符号分辨；σ 仍只决定分支反力差门，不放宽完整力。
  provenance 记录新算法版本与 `candidate_rule`。
- 对历史 b04 保存帧做只读判据复核：第 7 个中点 q = 2.2798038977 的右延续合格
  （完整力 0.00580621、横向力 0.00257815 eV/Å；反力差 0.00835765 < 门 0.01442498）。
  这不是新模型重跑，不修改旧结果。第 8 个中点 0.04318588 的分歧仍被拒绝，原因仍未测清。
- 新增 4 项解析回归：单侧达标的两种延续顺序、单侧达标不能绕过分支反力分歧、σ 不能放宽完整力。
  `test_ts_ef.py` **24 passed**；`make -C explore check` 格式、lint、全包 **182 passed**（73.33 s）。
  源码／输入 SHA256、逐中点离线结果与检查记录：
  `benchmarks/explore_d2/results/2026-10-25/hf_ef_candidate_v2_checks.json`。
- 当前执行环境 NVIDIA 驱动不可用（`nvidia-smi` 退出 9），未做 GPU 重跑或新噪声测量。
  驻点性质与两侧连接未运行；下一步是在可用 GPU 上重跑 b04，并独立核验候选及第 8 点差异。

### b04 用 HF-EF-bracket-v2 重跑：E/F 搜索得到候选，独立核验为一阶鞍点、连通 SN2（2026-10-25；本机 2080 Ti，冒烟）

对面改为 v2 候选判据（上一节）后，按其“下一步”在本机 GPU 重跑：`hf_ef_bounded.py results/2026-10-25/hf_ef_bounded_v2 0.0034 b04`。第一次在核验阶段因我的脚本没给候选挂 calculator 而崩（搜索已出候选；ledger/log 改名 `*_crashed` 保留），补一行后重跑通过。结果 `results/2026-10-25/hf_ef_bounded_v2/b04.json`。
- 搜索（只用 E/F）：`converged_candidate`，159 forward、19 s；q 2.27980，完整力 fmax 0.0031、横向 0.0030 eV/Å，反力 −0.0029；相对初始括区切向漂移 21°。
- 与旧 b04 C2 鞍点（fp64 升格精修）比：**Kabsch RMSD 2.1e-4 Å，ΔE −0.4 meV**（小于一个 fp32 能量量子 2e-3 eV）。
- 独立性质（完整 Hessian，只在搜索之后）：**一阶**，负本征值 −1.058，虚频 536i cm⁻¹。
- 独立连接（两侧下降）：两侧都到极小、无接力，`reactive`，键变化 C0–Cl1 ↔ C0–Cl5（SN2 氯交换），与旧 C2 指纹一致（方向相反）；两端能量 −26132.799 / −26132.797 eV。核验 173 forward；全程 45 s。
- 只是 b04 一个括区；b05 E2 的 4 个括区是几何跳支（`branch_unresolved`，第一个中点即两支相距 1.04 Å），v2 判据不影响该门，未重跑。

### J0：b05 E2 左支的反力抖动（2026-10-25；本机 2080 Ti，671 forward，100 s）

[跳支定位施工计划](PLAN_HF_JUMP_2026-10-25.md) J0，用户 10-25「是的」。`hf_ef_rf_jitter.py` 加 CASE 参数；点为 `b05_fp32_sA_10-13` 第一个中点 q = 1.3893 的左支（左右延续相距 1.04 Å，是两支，C 项跳过）。结果 `results/2026-10-25/hf_ef_rf_jitter_b05.json`。

| 测法 | 反力（eV/Å） | 同支几何两两最大 RMSD |
| --- | --- | --- |
| A. 10 次随机刚体变换 | 极差 0.027，σ 0.0082 | — |
| B. 8 次扰动重收敛，τ = 3e-3 | 极差 0.0087，**σ 0.0025** | **7.0e-3 Å** |
| 同上，τ = 1e-3 | 极差 0.0036，σ 0.0012 | 4.4e-3 Å |

- b05 的 σ_RF 与 b04（0.0034）同量级；但同支几何噪声 δ_geo 比 b04（6e-4 Å）大约 10 倍（b05 有软的片段 / 甲基运动），且收紧 τ 时反力抖动随之下降——与 b04 不同，b05 的抖动有一部分来自横向未松弛完。
- 供后续使用（默认 τ = 3e-3）：b05 `reaction_force_sigma` = 0.0025；J1 同支几何噪声尺度 δ_geo = 7e-3 Å。只在这一点测得。

### J1：双侧分支延续原语与真实模型对照（2026-10-25；本机 2080 Ti，与另一窗口的 make_data.py 共用 GPU）

[跳支定位施工计划](PLAN_HF_JUMP_2026-10-25.md) J1 / J3。
- 代码：新增 `explore/src/explore/ts_branch.py`（`continue_branches_ef`、`BranchSettings`、`BranchResult`，算法 `HF-EF-branches-v1`）：两端各自向对侧做约束松弛延续；本步对齐位移超过 `jump_factor × max(近几步位移率 × 步长, 实测同支几何噪声)` 时步长折半，折到 `min_step_q` 仍跳才判该支终止；某支反力由上坡转下坡即把该支上的这对点交给 `refine_ts_ef`，否则返回 `unresolved:hysteresis`。`ts_ef.py` 只把 `relax` 闭包提成模块函数 `_relax`（行为不变）。测试 `explore/tests/test_ts_branch.py`（解析折叠势报 hysteresis、终止点落在理论折叠点附近；光滑势转交得到候选；设置校验），`test_ts_branch.py` + `test_ts_ef.py` **30 passed**；全包 `make check` 未跑。
- 真实模型（`benchmarks/explore_d2/hf_jump_bounded.py`，PID 565630，结果 `results/2026-10-25/hf_jump_bounded/`）：σ / 几何噪声用 J0 实测值（b04 0.0034 / 6.2e-4 Å；b05 0.0025 / 7e-3 Å）。

| 括区 | 延续结果 | 交给 `refine_ts_ef` 后 | forward / 墙钟 |
| --- | --- | --- | --- |
| b04 | 低端支在 q 2.2695 判终止（误判，见下）；高端支在 q 2.2764–2.2864 反力变号 | `branch_force_disagreement` | 370 / 93 s |
| b05 第一轮 sA `6-9` | **低端（反应物）支连续，q 1.4293→1.4393 反力 +0.033→−0.028** | 同上 | 753 / 158 s |
| b05 第一轮 sA `10-13` | 同上 | 同上 | 780 / 143 s |
| b05 第一轮 sB `10-13` | 同上 | 同上 | 723 / 64 s |
| b05 fp32 sA `10-13` | 同上 | 同上 | 707 / 63 s |

- **E2 鞍点在反应物侧连续支上**：四个 b05 括区，该支变号括区两端与旧 E2 鞍点 Kabsch RMSD 0.006–0.020 Å、ΔE ±2 meV（≤ 1 个 fp32 量子），完整力 0.029–0.046 eV/Å。此前 `refine_ts_ef` 报 `branch_unresolved` 是因为粗括区（peak±1）的右端在产物支上。计划“未确认”一节中“E2 不在 q 约束极小支上”的推断不成立。
- b04 低端支终止是误判：步长折到 6.25e-4 时位移仍 6.7e-3 Å（门 3.1e-3），该处反力 0.87→0.77 变化陡，而高端支随后连续越过同一区域。跳支判据在陡而连续的区段会误报；本次不影响结果（另一侧找到变号）。
- **剩下的阻碍都是近根处的 `branch_force_disagreement`**：最后一个中点左右延续反力差 0.013–0.032（门 b05 0.0106、b04 0.0144）；其中 `6-9` 的右延续完整力已 0.0088（< `fmax_eV_A` 0.01），因门在先未出候选。
- 近根点补测（`hf_ef_rf_jitter.py` 加 SOURCE_JSON 参数；`results/2026-10-25/hf_ef_rf_jitter_root_{b05_fp32_sA_10-13,b04}.json`，2691 / 1272 forward）：

| 点 | B 重复（τ 3e-3）σ / 极差 | B 重复（τ 1e-3）σ / 极差 | 保存的那对延续 ΔRF | C：那对再收紧到 τ 1e-3 后 ΔRF |
| --- | --- | --- | --- | --- |
| b05 q 1.4343 | 0.0034 / 0.0097 | 0.0022 / 0.0066 | 0.027 | **0.0018** |
| b04 q 2.27949 | 0.0037 / 0.0136 | 0.0042 / 0.0148 | 0.032 | 0.010 |

  近根处重复收敛的 σ 与远处相近（b05 0.0025→0.0034，b04 0.0034→0.0037），**“近根噪声变大”的假设不成立**。保存的那对延续之差（0.027 / 0.032）比同点 8 次重复的极差大 2–3 倍；b04 那次左延续 0 步就满足 τ（从相邻 q 的结构直接投影过来）。b05 上把同一对继续收紧到 τ = 1e-3，差从 0.027 降到 0.0018。即差异来自“松弛一满足 τ 就停”时残留的起点依赖（软横向方向的残余位移），不是力噪声；b04 上收紧后 0.010 仍接近其重复极差，只部分支持。

### 二分中点收紧到 τ = 1e-3：5 个括区全部得到 E/F 候选，b05 `6-9` 找回 E2（2026-10-25→26；本机 2080 Ti，与另一窗口共用）

用户 10-25「可以试试」。`ts_ef.py`：`EFSettings.midpoint_transverse_fmax_eV_A`（可选，只用于二分中点的两次延续；端点与延续不变），`_relax` 加 `fmax` 参数；`test_ts_ef.py` 加一项，`test_ts_ef.py` + `test_ts_branch.py` **31 passed**。`hf_jump_bounded.py` 加 `--mid-tau`，σ 取近根处 τ = 1e-3 的实测（b04 0.0042、b05 0.0022）。PID 574627，结果 `results/2026-10-25/hf_jump_bounded_tau1e-3/`。

| 括区 | 搜索（只用 E/F） | 与旧鞍点 RMSD / ΔE | 独立 Hessian | 两侧下降 | forward：搜索 / 核验 |
| --- | --- | --- | --- | --- | --- |
| b04 | 候选，fmax 0.0012 | 2.5e-4 Å / −0.4 meV | 未决：`trivial_mode_floor_above_tolerance`（fp32 Hessian 刚体底门；虚频 536i） | 未执行 | 389 / 44 |
| b05 第一轮 sA `6-9` | 候选，fmax 0.0024 | 2.0e-3 Å / −1.2 meV | 一阶，885i | **连通，reactive：断 9–14，成 0–1、6–9**（与旧 C2 指纹相同） | 854 / 775 |
| b05 第一轮 sA `10-13` | 候选，fmax 0.0009 | 9e-4 Å / −0.1 meV | 一阶，891i | 一侧失败（`soft_mode_not_converged`） | 856 / 2120 |
| b05 第一轮 sB `10-13` | 候选，fmax 0.0084 | 2.3e-3 Å / −1.5 meV | 一阶，885i | 一侧 4 次接力后失败 | 861 / 1144 |
| b05 fp32 sA `10-13` | 候选，fmax 0.0010 | 7.3e-3 Å / 0 | 一阶，894i | 一侧失败（`soft_mode_moved_on_the_last_round`） | 840 / 1593 |

- **搜索**：5/5 由 E/F 延续 + 二分得到收敛候选，全部与旧精修器的鞍点重合（≤ 7e-3 Å，≤ 1.5 meV）；搜索不调用完整谱。
- **E2 找回（`6-9`）**：按计划 §4 的三条：E/F 搜索得到、独立 Hessian 一阶、两侧下降连通 E2 事件——满足。其余三个 b05 的鞍点同样一阶，只是连接未闭合；下降失败原因是 `descend_saddle` 内的软模跟随（PRRS 原有下降工具，HF 记录 F7 所说的仍依赖 Hessian 的连接环节），不是搜索失败。b04 的核验未决是 fp32 Hessian 自检门（fp32 重跑中 3–6% 出现过），v2 那次同一鞍点核验为一阶。
- **松弛步数（成本主项）**：b05 延续每点 ≈ 43 步（相邻点相距 0.01 Å，τ 3e-3，16 点 ≈ 680 forward）；τ 1e-3 中点每次 56–76 步；b04 延续每点 ≈ 12 步、中点 6–27 步。每次松弛都是新建 BFGS（初始 Hessian = 单位阵 × α），`constrained_scan` 已有的相邻点曲率复用这里没用上。核验（Hessian + 两侧下降）775–2120 forward，比搜索还贵。

### 已完成：分支曲率复用与切向预测，补齐松弛成本优化（2026-10-26；本机 2080 Ti）

- `ts_ef.py` / `ts_branch.py`：每条分支保存独立运行期状态，复用 E/F 学得的 BFGS 曲率，
  投影到新约束切空间并重置法向；用同支最近两个已接受点作有界切向预测。
  每次换 q 都新建优化器，不携带跨约束 secant 历史；零步收敛保留继承的曲率。
  试步先复制状态，几何／反力分支检查接受后才提交；拒绝、折半和失败不污染原状态。
  同支括区转交二分时传递两端各自的状态。完整力、横向力、噪声门与分支判据未放宽。
- `EFSettings.reuse_curvature`、`predict_tangent` 默认开启，均关闭为冷启动对照。
  算法版本 `HF-EF-bracket-v3` / `HF-EF-branches-v2`；设置入 provenance，逐帧保存
  `curvature_reused`、`predictor_step_A`。运行期提示不作恢复证书，未接 D2。
- 沙箱内 `nvidia-smi` 报驱动不可用；沙箱外只读检查确认本机 2080 Ti 可用且空闲，
  本轮模型与 GPU 实测均在沙箱外执行。没有因此改用另一 SSH 入口或复制工作区。
- 新脚本 `benchmarks/explore_d2/hf_warm_start_bounded.py`：同一常驻 fp32 POLAR-1，
  只读同批冻结括区，端点 τ 3e-3、中点 τ 1e-3；σ 用此前对应体系／阈值实测值。
  每个搜索独立账本、上限 3000 forward；性质／连接核验另开账本，上限 4000。

同轮 A/B（单次／模式；表中为 **搜索 forward / 秒**）：

| 括区 | 冷启动 | 仅曲率 | 仅预测 | 两者合用 |
| --- | --- | --- | --- | --- |
| b04 | 413 / 37.9 | 242 / 21.5 | 341 / 30.1 | **215 / 18.9** |
| b05 第一轮 sA `6-9` | 859 / 76.1 | 303 / 27.2 | 445 / 39.9 | **284 / 25.5** |

- 两例搜索 forward 分别减少 **47.9% / 66.9%**，墙钟约 **2.00× / 2.99×**。
  b05 `6-9` 延续每点中位步数 43 → 6；两次中点松弛 64／74 → 22／27。
  这是一轮有界性能实测，不是逐位一致性或普适加速声明。
- 合用方案 **5/5** 得到达标候选：b04 215、b05 `6-9` 284、第一轮 sA `10-13` 259、
  第一轮 sB `10-13` 244、fp32 sA `10-13` 256 forward。
  所有新候选完整力 ≤ 0.01、中点横向力 ≤ 0.001 eV/Å；约束残差合格。
- 两例独立核验：b04 **一阶、536i、连通 SN2**（断 0–1、成 0–5；173 forward）；
  b05 `6-9` **一阶、886i、连通 E2**（断 9–14、成 0–1／6–9；715 forward）。
  新 b05 候选相对旧鞍点 RMSD 0.00181 Å、ΔE −3.11 meV；b04 0.000703 Å、−0.415 meV。
  另外三个新候选只复核搜索与旧鞍点对照，本轮没有重复其独立核验。
- 新增 8 项解析回归：实际计费降低并保持解析驻点、冷／热启动均保留真实滞回、
  拒绝试步回滚、预算／能量失败、曲率法向投影与零步保留、坐标系变更重置、开关类型校验。
  相关 **39 passed**；`make -C explore check` 格式、lint、全包 **197 passed**（89.05 s）。
- 结果、逐阶段账本、两次运行日志、检查日志、前后源代码快照与汇总均在
  `benchmarks/explore_d2/results/2026-10-26/hf_warm_start/`，入口 `SUMMARY.json`。
  初版性能驱动错误地从不存在的 `property.status` 读取摘要标签，已按实际 `saddle_order=1`
  修正；保留初版源码和修正说明，原始搜索／性质／连接数据未变，历史结果未覆盖。
  账本费用、候选阈值及源码 SHA256（含初版归档）逐项核对通过。

搜索成本缺口已补齐；b05 `6-9` 当前搜索 284、核验 715 forward，后者仍为主要成本。
F7 的 Hessian 依赖连接追踪仍是后续工作，本轮没有改 `descend_saddle` 或添加 J2 回退。

## 已完成：F7 只用 E/F 的连接追踪 `trace_sides_ef`（K1–K3；2026-10-26；本机 2080 Ti，冒烟）

[连接追踪施工计划](PLAN_HF_CONNECT_2026-10-26.md)，用户 10-26「好的」。
- 代码：新增 `explore/src/explore/ts_connect.py`（`trace_sides_ef`、`branch_starts`、`ConnectSettings`、`ConnectionResult`，算法 `HF-EF-connect-v1`）。起点取候选自身分支上最近的已收敛帧，要求 |dE/dq| ≥ 3σ 且完整力 > 下降判据（`quench_fmax_eV_A` 0.03）；某侧没有满足条件的帧时，沿搜索切向取最小的、能量下降超过 3 个 fp32 量子的位移。撤约束后用 BFGS 只靠 E/F 下降，不做扭转打磨。端点资格：3 次 σ 0.03 Å 扰动后重下降，要求同键图，且没有超过“3 量子 + ‖F_end‖₂·‖Δx‖₂”的能量下降（凸盆一阶界）；`same_basin` 只记录，不作判据。不调用 `confirm_minimum` / `curvature_spectrum` / `descend_saddle`，无接力。
- 测试 `explore/tests/test_ts_connect.py`（二维双盆解析势：分支帧起步两侧到两个极小；无帧时退回切向；鞍点作端点被扰动复核识别；预算截断计费；起点完整力条件；设置校验）。TS 相关 `test_ts_connect/ts_ef/ts_branch/ts_warm_start` **48 passed**；全包 `make check` 未跑。
- 第一轮（`results/2026-10-26/hf_connect_v1/`）暴露两个问题，已改：① 起点帧完整力已低于 0.03，下降 0 步停在 TS 旁（4 个括区中出现），扰动后能量掉 0.4–2.4 eV、被判 `ef_unstable`，没有误判，但也没有结论 → 加“起点完整力 > 下降判据”；② 松散的 b05 产物复合物在 0.03 Å 扰动下回到相邻构象（0.21 Å，+0.002 eV），被 `same_basin` 判为不同盆 → 改为同键图 + 无超界能量下降。
- 第二轮（`benchmarks/explore_d2/hf_connect_bounded.py`，结果 `results/2026-10-26/hf_connect/`，输入为对面热启动的 5 个候选，只读）：

| 候选 | 连接 | 与旧 C2 键变化 | 连接 forward（新） | 旧核验 forward（Hessian + `descend_saddle`） |
| --- | --- | --- | --- | --- |
| b04 | connected，SN2 | 相同 | **90** | 173（热启动核验）/ 44（τ 1e-3 那次：Hessian 自检未决，未下降） |
| b05 `6-9` | connected，E2 | 相同 | **330**（高侧用切向起步） | 715 / 775 |
| b05 第一轮 sA `10-13` | connected，E2 | 相同 | **410** | 2120（下降失败） |
| b05 第一轮 sB `10-13` | connected，E2 | 相同 | **362** | 1144（下降失败） |
| b05 fp32 sA `10-13` | connected，E2 | 相同 | **400** | 1593（下降失败） |

- 5/5 连通，键变化全部与旧 C2 一致；旧链在后三个上下降失败。标签显示 b05 为 `dissociation`、旧为 `reactive`：`classify` 按“终点片段数多于起点”判 dissociation，新链按低侧→高侧（2→3 个片段），旧链方向相反，事件相同。
- 等级说明：这是 E/F 等级的连接（端点通过有限扰动复核），不是 Hessian 认证的极小，也不是 IRC；候选的一阶性质仍需独立检验（计划 K4，未做）。σ 0.03 Å 扰动的能量升高 0.05–0.9 eV，全部可分辨。

## 已完成：HF-4 H1–H4，E/F 链接入 D2（规则 D2-HF-v1；2026-10-26；本机 2080 Ti 冒烟）

[HF-4 施工计划](PLAN_HF4_D2_2026-10-26.md)，用户 10-26「开始」。
- `explore/src/explore/d2.py`：`HFSettings`（σ、几何噪声、中点 τ、`hs_cap`、`hc_cap`，全部必填）；`init_run(..., hf=)` 冻结 `D2-HF-v1`，旧 v1.1 按 manifest 分流不变；HF 模式第二账户名 `refine`（额度 = `certify_budget`）；升级项 L3 释放照旧，L4 换成 `HS|`（`continue_branches_ef` → `refine_ts_ef`，括区：脊 peak±1、夹住的跳支 反应物侧 + L2 产物侧、抽查 peak±1）与 `HC|`（`branch_starts` + `trace_sides_ef`，只在 HS 出候选时派发）；HC 连通后按 化学对 + 鞍点能量 + 对称 RMSD 去重（`hf_duplicate`）；HS payload 去掉逐帧力（只留位置）。报告新增 `hs_status`、`hc_status`、`hf_forwards`、`hf_events`、`property_status: not_run`。
- `explore/src/explore/ts_branch.py`：修“走到对侧端点”的判断——原用 1e-12，而约束只收敛到 `constraint_tol` 1e-8，b04 冒烟里同一端点被重复松弛到 400 点（白花约 380 forward）；改为 `constraint_tol`，`test_ts_branch.py` 加回归项。
- 测试：`explore/tests/test_d2_hf.py`（TwoSaddle 解析台 HF 端到端：两个解析鞍点中定位并连通、外部计数 == 账本 == 两账户之和；HS 中 / HS 后第一个 HC forward / HC 中三个杀点续跑：每个任务只结算一次、已结算的 HS 几何不变、费用 = 不中断 + 浪费、外部计数不超付）；`d2_bench.py` 的 HF 输入按真实体系做法在解析台近根点实测（σ 6.7e-4、几何噪声 3.3e-4 Å，我起初拍的 1e-4 让 HS 被误拦）；`d2_driver.py` 加 `hf` 参数。`make -C explore check`：格式、lint、**213 passed**（H3 后）；ts_branch 修复后相关 `test_ts_branch` + `test_d2_hf` 14 passed。
- b04 冒烟（`run_b04_d2.py --hf`：一个 fp32 模型；HF 时 refine 额度 20000，因为每次派发预留整个 3000 上限，第一次用旧 3000 额度时 6 项里 5 项被跳过，结果改名 `results/2026-10-26/b04_hf_smoke_v0*` 保留）：`results/2026-10-26/b04_hf_smoke/`：6 个升级项 → HS converged 2 / hysteresis 2 / above_ceiling 2；HC connected 1（**SN2，鞍点 −26132.2012 eV，断 0–1、成 0–5**）+ hf_duplicate 1（另一意图找到同一鞍点）；refine 共 791 forward（HS 610、HC 181）；HF 段约 0.089 s/forward；unreconciled 空。
- b05 冒烟（`run_b05_d2.py --hf`，与 10-24 v1.1 冒烟同参数：1 姿态 12 意图，refine 12000）：`results/2026-10-26/b05_hf_smoke/`：端到端跑通，1 个脊 → HS `hysteresis`（959 forward；低侧在 q 4.2432 终止，下一步 6e-4 处位移 0.60 Å、ΔE −0.28 eV、反力 +0.86 → −0.22，是真实跳支；高侧到端点无变号），HF 段 0.089 s/forward，账目闭合。真实滞回需要 J2（Dimer / 有限幅度响应），本轮不做。
- H5：PBS 脚本 `benchmarks/explore_d2/b05_hf_pbs.zsh`（照模板，只改 ← 行；节点待填 `NODE_TBD`），两块与 `b05_fp32_s{A,B}`、`b05_v11_s{A,B}` 同种子同参数，`--hf`，结果写 `results/2026-10-26/b05_hf_s{A,B}/`。

## 已完成：重复响应离线统计（2026-10-26；只读 journal，不调模型）

起因：用户 10-26 提出“在昂贵精修前识别同一种响应、每组精修少量代表”，要求先离线统计。脚本 `benchmarks/explore_d2/repeat_response_census.py`，结果 `benchmarks/explore_d2/results/2026-10-26/repeat_response_census/`（每个运行一个 json + log）。
- 数据：已完成的 b05 PBS 运行 `b05_fp32_s{A,B}`（D2-v1）、`b05_v11_s{A,B}`（D2-v1.1），均为 Hessian 路径 L4；b04_d2 只当程序检查。HF 的 `b05_hf_s{A,B}` 还在跑，没纳入。
- 分组代理：同一姿态 + L3 释放终点落在同一盆（现成盆地规则：同化学键 key、|ΔE| ≤ 0.05 eV、对称 RMSD ≤ 0.15 Å），按 journal 顺序贪心建组。没有 L3 终点的项（L3 因 search 额度耗尽被跳过 95–201 项，或 not_minimum）不分组、全部保留。
- 结果（4 个运行合计）：L4 共 1293 项、765 201 forward。可分组 487 项落进 72 组，最大组 15–23 项；72 组里 58 组终点与姿态同化学键 key（构象响应）。“每组只精修第一项”跳过 251 276 forward（各运行 21–43%），但 13 个已认证事件只保留 11 个：丢掉的 2 个都是无键变化的 c2（fp32_sA ΔE −0.96 eV、RMSD 2.63 Å；v11_sB ΔE +0.008 eV），每组保留 2 项也找不回。事后去重 `duplicate_c2` 的费用：v1.1 合计 18 736 forward（certify 的 6.4%），v1 合计 895。
- 代理质量：36 个 `duplicate_c2` 里只有 8 个和原 c2 在同一终点组。同一终点组里 c2、duplicate、各类 unresolved 混在一起，说明 L3 终点既不决定 L4 结果，也不能代表“同一个鞍点”。
- 局限：只有 b05 一个体系，sA/sB 共用同一个种群；Hessian 路径的 L4 费用不能直接换算成 HF 的 HS/HC 费用；代理是 L3 终点，不是扫描过程中的几何响应。
- **续：按 L4 起始猜测分组**（`benchmarks/explore_d2/guess_separability.py`、`attempts_per_start.py`；结果同目录 `guess_b05_v11_s{A,B}.{json,npz}`、`attempts_per_start.json`；只算 v1.1 两块，36 个 duplicate_c2 里 35 个在这两块；每块两两对称 RMSD 约 200 s，20 进程）：
  - 起始猜测 = 冻结规则 M5-D1 `_promotion_points` 的猜测点，重算的分类与 payload 记录 0 处不一致。跨姿态比较（所有姿态同一化学 key，|Aut| = 1296）。
  - 同一事件内的猜测两两 RMSD 中位 0.23 Å、最大 0.58 Å；不同事件之间最小 0.35 Å。猜测离自己的姿态 0.2–3 Å，但彼此常只差 0.003–0.05 Å（< 0.05 Å 的对：sA 701、sB 874，大部分跨姿态）：很多扫描被推到同一个约束几何。
  - 起点几乎相同、结果不同：< 0.05 Å 的对里，“一个找到事件、一个 unresolved”的有 67（sA）/ 111（sB）；这些对的扫描切向 |cos| 多在 0.1 左右，也有 0.95–0.98 仍然不同。含事件的簇里，单次尝试的成功率从 1/18 到 6/9 不等。重复精修实际上起了“重试”的作用。
  - 分组半径取实测几何抖动 0.0057 Å（同一扫描点重复约束弛豫的最大两两 RMSD，`results/2026-10-25/hf_ef_rf_jitter_root_b05_fp32_sA_10-13.json`，τ 3e-3）的倍数。每簇最多精修 m 次（按 journal 顺序）：

| 半径 | m | sA 省 forward / 事件 | sB 省 forward / 事件 |
| --- | --- | --- | --- |
| 2× 抖动（0.0114 Å） | 1 | 38% / 3 of 4 | 41% / 3 of 4 |
| 2× 抖动 | 2 | 26% / 4 of 4 | 25% / 4 of 4 |
| 4× 抖动（0.0228 Å） | 2 | 37% / 4 of 4 | 41% / 4 of 4 |
| 10× 抖动（0.057 Å） | 3 | 36% / 4 of 4 | 35% / 4 of 4 |

  - 局限：两块共 8 个事件，全部无键变化；结论是 Hessian 路径 certify 上的事后重放，不等于 HF 路径；v1.1 certify 额度未用满，省下的 forward 当时换不来更多结果。

## 已完成：HF 运行上的重复响应统计 + 去重规则漏判（2026-10-28；只读 journal，不调模型）

TODO 第 1 条（[TODO.md](TODO.md)）。输入：`benchmarks/explore_d2/results/2026-10-26/b05_hf_s{A,B}/`（10-26 在 yayoi28 跑完，`report.json` 已生成）。脚本：`guess_separability.py` 加 HF 模式（项 = HS 任务；结果取 HC 状态，没有 HC 时取 HS 状态；费用 = HS + HC），新增 `hf_event_symmetry.py`，`attempts_per_start.py` 加 `--merge-symmetric`。结果在 `benchmarks/explore_d2/results/2026-10-26/repeat_response_census/`：`guess_b05_hf_s{A,B}.*`、`hf_event_symmetry.json`、`attempts_per_start_hf{,_merged}.json`。
- 运行概况：refine 额度用尽（237k / 238k of 240k），HS 因额度被跳过 89 / 96 项。HS 占 refine 的 99%。HS 结果 hysteresis 170 / 174、transverse_not_converged 56 / 55、converged_candidate 8 / 6。HC connected 4 / 6，全部是 E2（断 C–Cl，H 转移到 O 14），化学对全部相同。
- **去重漏判（已坐实）**：“不同”事件之间，用姿态图自同构群（|Aut| 1296）算对称 RMSD，同类只差 0.002–0.006 Å，能量差 ≤ 0.004 eV，远在规则自己的 0.15 Å / 0.05 eV 之内；每块实际只有 2 个不同鞍点（两类相距约 0.57 Å）。原因：`d2.py:724` `_hf_duplicate` 在鞍点几何上求自同构群，转移中的 H 打破了甲基对称（|Aut| 144）。同一对鞍点：鞍点图群下 1.673 Å，姿态图群下 0.0021 Å。同种原子的置换是势能面的严格对称，所以这些是同一个反应换了编号。`certify.py:224`（Hessian 路径）写法相同，v1/v1.1 事件有没有同样的问题【未核实】。未改代码。
- 按起始猜测分组、每簇限 m 次（半径为实测几何抖动 0.0057 Å 的倍数），事件按对称合并后计数：

| 半径 | m | sA 省 refine forward / 事件 | sB 省 refine forward / 事件 |
| --- | --- | --- | --- |
| 2× 抖动 | 1 | 35% / 2 of 2 | 34% / 2 of 2 |
| 4× 抖动 | 1 | 48% / 2 of 2 | 53% / 2 of 2 |
| 10× 抖动 | 1 | 62% / 2 of 2 | 61% / 2 of 2 |

  不合并时 m = 1 会“丢” 1–4 个事件，丢的全是对称重复。含事件的簇里，第一项就连通（first found at 1）。这和 Hessian 路径（单次成功率 1/18–6/9）不同。
- 局限：每块只有 2 个不同事件，全是同一化学对的 E2；这是事后重放，没有验证省下的额度换成那 89 / 96 个被跳过的项之后能多找到什么。

## 已完成：去重规则改用 鞍点 + 端点 的对称群并集（2026-10-28；用户「可以 开修」）

TODO 第 6 条。问题与证据见上一节。
- `explore/src/explore/certify.py`：新增 `pair_symmetries(saddle, endpoints, config)`，返回鞍点图与各端点极小图自同构群的并集（去重后）和第一次出现的预算溢出信息。certify 第 5 步改用它（端点为两侧下降得到的极小），溢出仍判 `unresolved:identity`，与原来一致。
- `explore/src/explore/d2.py`：`_hf_duplicate` 增加 `ends` 参数（HC 两侧下降终点），改用 `pair_symmetries`；调用处传 `ends.values()`。
- 并集包含原来的鞍点图群，所以原规则判出的重复一个不少，只会多认出对称等价的鞍点。同种原子的置换都是势能面的严格对称，判据没有放宽阈值（仍是 basin_rmsd_A 0.15 Å、basin_energy_eV 0.05 eV）。
- 测试：`tests/test_certify.py::test_pair_symmetries_keep_the_minimum_symmetry_the_saddle_graph_lost`（CH3…O 桥连 H：鞍点图群下 > 0.15 Å，并集下 < 1e-6，并集 ⊇ 鞍点图群）；`tests/test_d2_hf.py::test_hf_duplicate_merges_a_symmetry_equivalent_saddle`。变异检验：把 d2 换回只用鞍点图群，新测试失败。ruff format/check 通过；`test_certify.py test_d2.py test_d2_hf.py test_d2_resume.py` **52 passed**；全包 `make check` 未跑。
- 已有 journal 重放（`benchmarks/explore_d2/replay_duplicates.py`，结果 `results/2026-10-26/repeat_response_census/replay_duplicates.json`；HF 直接调用修改后的 `_hf_duplicate`）：b05_hf_sA 4 → 2、b05_hf_sB 6 → 2；b05 v1.1 s{A,B}、v1 s{A,B}、b04_d2、b04_hf_smoke 均不变（Hessian 路径没有这种漏判）。
- 未做：没有用新代码跑一次真实 D2；没有改 rule_version（`D2-HF-v1` / `D2-v1.1` 的 manifest 不变，已有运行的 report 不改写）。
- 冒烟（本机 181，2080 Ti；`run_b04_d2.py --hf`，noise_eV 0.001953125，与 10-26 b04 HF 冒烟同参数；`benchmarks/explore_d2/results/2026-10-28/b04_hf_smoke_dedupfix/`）：退出码 0。升级 3 脊 + 3 抽查；HS converged 2 / hysteresis 2 / above_ceiling 2；HC connected 1（SN2，鞍点 −26132.2012 eV，断 0–1、成 0–5）+ hf_duplicate 1，与 10-26 相同，去重新路径执行到了。forward：HS 599、HC 189（10-26 是 610 / 181，fp32 GPU 不确定性引起的小差异）；unreconciled 空。b05 上的漏判只在 journal 重放里验证过，没有用新代码重跑 b05。

## 已完成：D2-HF-v1.1，同一起点的 HS 两轮派发（2026-10-28；用户「可以」）

[施工计划](PLAN_START_CLUSTER_2026-10-28.md)。用户 10-28 原话：“先减少相似起点上的重复失败，再研究单项 HS 的早停。”
- 代码：`explore/src/explore/d2.py` 新增 `StartCluster`、`HF_SC_RULE_VERSION = "D2-HF-v1.1"`、`_start_clusters`、`_rmsd_lower_bound`、`_hf_tasks`、`_tags`；`init_run(..., start_cluster=)`；`_promote` 两轮派发（第二轮跳过原因 `HS_deferred:same_start`）；HS payload 记 `start_cluster`；报告增 `start_clusters`。旧 manifest（无 `start_cluster`）照 v1 走。`benchmarks/explore_d2/run_b0{4,5}_d2.py` 加 `--start-cluster`（`StartCluster(2.0, 1)`，半径 = 2 × 冻结的 `geometry_noise_A`）。
- 测试：`explore/tests/test_d2_hf.py` 新增 4 项（设置与冻结、下界 ≤ 对称 RMSD、额度充足时 HS 状态与 v1 完全相同、额度只够一项时其余记 deferred 且外部计数 == 账本）；两个变异均使测试失败。ruff 通过；`test_certify.py test_d2.py test_d2_hf.py test_d2_resume.py` **56 passed**；全包 `make check` 未跑。
- 重放（`benchmarks/explore_d2/replay_start_clusters.py` 调用施工后的 `_start_clusters`，结果 `results/2026-10-26/repeat_response_census/replay_start_clusters.json`；半径 0.014 Å、每簇 1 项）：sA 325 个升级项 → 203 簇，sB 336 → 202 簇；10-26 实际跑过 HS 的项里属于第一轮的 141/236、140/240；原 refine 花费有 40% / 42% 落在第二轮；对称合并后的 2 个不同事件都在第一轮。分簇本身 97 s / 112 s CPU。
- b04 冒烟（本机 181，`run_b04_d2.py --hf --start-cluster`，`results/2026-10-28/b04_hf_v11_smoke/`）：退出码 0，规则 `D2-HF-v1.1`；半径 0.00124 Å 下 6 个升级项各成一簇、第二轮为空，所以结果与 v1 相同（HS converged 2 / hysteresis 2 / above_ceiling 2；HC connected 1 SN2 + hf_duplicate 1；refine 740 forward；unreconciled 空）。**第二轮推迟在真实模型上没有被触发**，只有解析台测试覆盖；b05 运行脚本的 `--start-cluster` 只过了参数解析。
- 计划 §5 ④：PBS 脚本 `benchmarks/explore_d2/b05_hf_v11_pbs.zsh`（照 `b05_hf_pbs.zsh`，只改 ← 行；节点 `NODE_TBD`；与 10-26 `b05_hf_s{A,B}` 同种子同参数，只差 `--start-cluster`），结果写 `results/2026-10-28/b05_hf_v11_s{A,B}/`。

## 已核查：b05 D2-HF-v1.1 正式结果（2026-10-28；离线读盘，无模型调用）

- 用户告知新任务已算完。新两块 `results/2026-10-28/b05_hf_v11_s{A,B}/` 最终 report 齐全；四块新旧运行的 journal 结算费用、ledger close 与 report 账户逐项一致，held=0、unreconciled 为空。写入节点为 yayoi04；本轮直连时 DNS 解析失败，只核对 kasuga01 可见的 NFS 文件，未作写入节点实时/退出码核验。
- 同样 refine 额度下，HS 结算 236→290 / 240→275（排除 0-forward 异常后 283/271）；实际计费起点簇 141→196 / 140→201。第一轮派发 203/205 簇，第二轮执行 87/70 项，暂缓 42/55 项。旧覆盖数取已有相同半径规则的重放；实际新旧 L1/L2 有差异，不能声称完全相同轨迹的纯顺序对照。
- refine 237127→237635 / 238176→237193 forward，几乎用满，没有实际节省总额度。两块各保住 2/2 类 E/F 连通事件，没有新增事件；独立对称几何匹配最大 RMSD 0.00585 Å、能量差 0.001953125 eV，chemical_pair 全匹配。property_status 仍 not_run。
- **新反例：sB 的第二类事件由第二轮补回。** 簇首 HS 收敛、HC `failed:tangent_no_resolved_drop`，rank=1 成员 connected；不能将每簇首轮一项改成永久删除后续成员。
- **新缺陷：HS 7/4 项 `error:ValueError`，均 0 forward、branch_bracketed。** 原始报错 `bracket must have positive width`；从 L1/L2 重建，10 项 q 宽度为 0，1 项为 2.22e-16 Å。同 q 两支共存证据被 `_hf_bracket` 交成正宽搜索区间，接口语义需要修复；未改代码或历史结果。
- 详细结论与限制：[HF_START_CLUSTER_RESULTS_2026-10-28.md](HF_START_CLUSTER_RESULTS_2026-10-28.md)。原始指针、逐项状态/费用、几何匹配、异常括区与分析脚本快照：`benchmarks/explore_d2/results/2026-10-28/hf_v11_comparison/summary.json`。下一步先修括区交接、离线回归，不启动新大批运行；尚不能宣布完整验收或提速。
- 补充（Claude）：上述数字已用 report 独立复核。括区缺陷在 10-26 v1 两块中被预算挡住、从未触发（`branch_bracketed` 项一个都没结算 HS），v1.1 第一轮才第一次覆盖到，不是两轮派发引入的回归。新两块 `hf_duplicate` 7 / 3、`connected` 各 2：去重修复首次在真实运行中生效。见结果记录“补充”节。

## 已完成：跳支括区交给 HS 改用正向扫描的跳变区间（TODO 第 8 条；2026-10-28；用户「是的继续」）

- 原因：`_hf_bracket` 对 `branch_bracketed` 用（反应物点，L2 反扫点）作 HS 括区。L2 反扫点就在反应物点的 q 上，它只是“两支在该 q 共存”的证据，不是区间端点；交给 `continue_branches_ef` 后宽度为 0 或机器精度。离线重建（`benchmarks/explore_d2/replay_hf_brackets.py`，结果 `results/2026-10-28/hf_bracket_fix/replay_brackets.json`）：b05 HF 四块共 123 个跳支项，**旧括区全部是 0 或 ≤ 9e-16 Å**；v1.1 里 7/4 项报 `positive width` 错误，另外 4/3 项 HS 跑在机器精度宽度的区间上，结果没有意义。10-26 两块 v1 的跳支项被额度挡住，从未进 HS。
- 修法（`explore/src/explore/d2.py`）：`_classify` 在 `branch_bracketed` 的 info 里加 `product_side_scan`（正向扫描跳变后一点的下标）；`_hf_bracket` 取（`points[reactant_side]`，`points[product_side_scan]`），即正向扫描自己测到的跳变区间，两端在两条支上，各自延续（HF 记录 §5.2）。不扰动 q；L2 证据仍保留在 journal 与分类里。`_hf_bracket` 去掉不再使用的 `intent` 参数。Hessian 路径（`_promotion_points`）不变。重建后 123 项新括区宽度全部为一个扫描步长 0.1 Å，两端键图都不同。
- 规则版本：原地修复，没有新版本号。改动只影响 HF 规则下跳支项的 HS 输入；10-28 及以前的 HF 运行，跳支项的 HS 都是用旧括区跑的，引用时注明。
- 测试：`explore/tests/test_d2_hf.py::test_bracketed_jump_hands_the_forward_jump_interval_to_hs`（经 `_classify` 的真实跳支路径；L2 点与反应物点 q 完全相等 / 差 2.22e-16；正向与反向扫描，共 4 组）。变异检验：换回旧括区，4 组全部失败。ruff 通过；`test_certify test_d2 test_d2_hf test_d2_resume test_probe test_ts_branch` **92 passed**；全包 `make check` 未跑。
- 冒烟（本机 181，`benchmarks/explore_d2/hf_bracket_smoke.py`，只读 v1.1 sA journal，取前 3 个报错项，按 `_hs_task` 的设置跑 HS，每项上限 3000；`results/2026-10-28/hf_bracket_fix/smoke/`）：退出码 0，3 项都不再报错，实际计算 201 / 218 / 241 forward（20–24 s）。三项都是 `unresolved:transverse_not_converged`：产物侧端点在端点松弛的步数上限内横向力没收敛（第 3 项反应物端先走了一段）。这与普通 HS 的 55 个同类失败属同一类，不在本次修复范围内；没有找到新的 TS。

## 已完成：HS `transverse_not_converged` 离线诊断（2026-10-28；只读 journal，不调模型；用户「可以」）

脚本 `benchmarks/explore_d2/hs_transverse_diag.py`，结果 `results/2026-10-28/hs_transverse_diag/diag.json`；输入 b05 HF 四块（10-26 v1 s{A,B}、10-28 v1.1 s{A,B}），共 221 个此状态的 HS。失败的松弛按所在阶段分两类，四块一致：

| 阶段（每块项数） | 阈值 | 步数 | 最小横向力 | 能量变化 | 首末对齐 RMSD |
| --- | --- | --- | --- | --- | --- |
| refine 二分中点（10 / 11 / 11 / 11，共 43） | 1e-3 | 全部 200 | 0.0010–0.0017 | −0.01 ～ +0.002 eV（≤ 5 个 fp32 量子） | 0.004–0.085 Å |
| branch 端点 / 延续（46 / 44 / 44 / 44，共 178；其中端点 2 / 4 / 3 / 3） | 0.003 | 全部 200 | 0.003–0.03 | 中位 −0.27 ～ −0.32 eV，最多 −2.4 eV | 0.16–2.1 Å（中位 0.6–1.1） |

- refine 类：两侧延续已经找到连续最大值括区，第一个中点的松弛没收敛。上表“最小横向力”是过程最低值；**末帧横向力中位 0.0075、最高 0.0475 eV/Å（末/最低中位 5.7 倍），只有 6 项末帧 ≤ 0.003**，即“曾接近阈值，随后反弹”，不是末态平台。结构与能量基本不变。是 fp32 噪声还是优化器振荡，尚未区分（10-28 更正：先前写“停成平台、贴着噪声底”不成立）。同期每块 `converged_candidate` 6–9 个；这 43 项是待挽救的搜索，不能预计为新增候选。
- branch 类：约束松弛中结构持续大幅移动并下滑，**疑似失支，未确认换盆**：大位移与能量下降也可能来自同一分支上的缓慢松弛，需缩小延续步长看响应才能区分（10-28 更正：先前写“换盆”越过了证据）。这时 `_relax` 不返回，延续的步长折半（J1）没有机会处理，整个 HS 直接结束。10-28 跳支括区冒烟的 3 项也属于这一类。
- 成功的松弛（末帧 ≤ 0.003，每块约 9.3–9.8k 次）：中位 4–5 步，p99 161–170 步，≥150 步的约 2%。**注意：脚本统一用 0.003 判成功，中点实际阈值是 1e-3，失败中点段里有 6 段末帧 ≤ 0.003 被误计为成功；且成功样本的分布不能排除失败尾部需要更多步。** 此条统计待按阶段修正。
- 未改代码。两类的处理都涉及判据或语义，见 TODO 第 9、10 条。

## 已完成：延续试步松弛失败时有界折半重试（TODO 第 10 条；2026-10-28；用户「是的」）

[施工计划](PLAN_HF_WALK_RETRY_2026-10-28.md)。
- 代码（`explore/src/explore/ts_branch.py`）：`_walk` 只接试步里的 `_Stop("transverse_not_converged")`；当前点、state、rates 不变（失败试步写过的曲率在 `trial_state` 里，被丢弃），步长折半；降到 `min_step_q` 仍失败才终止该侧，`termination = "relax_unconverged"`，诊断（`_unconverged`：试步 q、步长、末帧横向力 / 完整力、ΔE、对齐 RMSD、位置）写进 `relax_unconverged` 与 `unconverged_trials`。疑似跳变终止的一侧加 `termination = "jump"`，原 `jump` 字段不变。失败试步的帧照常计费、留在 history，不进 points / rates / 反力判断。都没找到最大值且有一侧 `relax_unconverged` 时，总状态为 `unresolved:walk_relax_unconverged`。其他 `_Stop`、`BudgetExhausted`、`GateRejected` 照旧上抛；端点松弛不改。provenance `HF-EF-branches-v3`；D2 规则版本不变。无失败试步时，计算路径和 forward 不变，输出多一个 `termination` 字段。
- 测试（`explore/tests/test_ts_branch.py` 新增 5 个测试函数，其中一个参数化 3 组）：大步失败折半后走完、失败帧不进 points、已接受 state 不被失败试步污染；每步都失败时按折半序列降到 `min_step_q` 后终止、诊断齐全；失败帧反力反号也不能括出最大值；其他停止原因照旧上抛；解析台把 `max_relax_steps` 设为 1 时端到端为 `walk_relax_unconverged`，外部计数 == 账本。变异检验（变异代码放在临时目录，用 `-o pythonpath=` 覆盖 pyproject 的 pythonpath，同参数的真实源码对照 4/4 通过）：M1“失败帧当已接受点”4 项失败；M2“失败即终止、不折半”2 项失败。ruff 通过；`test_ts_branch test_ts_ef test_ts_connect test_ts_warm_start test_d2_hf test_d2 test_d2_resume` **94 passed**；全包 `make check` 未跑。
- 冒烟（本机 181 2080 Ti，`benchmarks/explore_d2/hs_walk_retry_smoke.py`，取 v1.1 sA 中延续阶段失败的前 3 项，按 `_hs_task` 设置重跑；`results/2026-10-28/hf_walk_retry/smoke/`）：退出码 0。① 原 1108 forward 失败 → 新 `budget_exhausted` 用满 3000（一侧跳变终止，另一侧走到额度用完；本次无失败试步）；② 原 245 → 新 `walk_relax_unconverged` 1972（低侧失败试步折半 9 次至 0.000625 仍失败而终止，高侧到端点）；③ 原 1020 → 新 `hysteresis` 1732（本次无失败试步）。没有新的最大值或候选，单项花费上升。**3 项中 2 项的原失败在本机没有复现**（原运行在 yayoi04 3090，fp32 数值随硬件不同），所以旧记录与新运行的差别不能归因于代码；有界对照须在同一台机器上把新旧代码各跑一遍。

## 已核查：b05 D2-HF-v1.1 用当前代码重跑（yayoi08，2026-10-28→29；离线读盘）

当前代码 = 10-28 去重修复 + 跳支括区修复 + 延续试步折半重试（HF-EF-branches-v3）。PBS `b05_hf_v11r.o19794`（项目根目录）末行 `done 2026-10-29 01:26:42`；`benchmarks/explore_d2/results/2026-10-28/b05_hf_v11r_s{A,B}/` 两份 report 齐全，`ledger_total` == search + refine，unreconciled 空。对照：同种子、同型号节点的 10-28 yayoi04 `b05_hf_v11_s{A,B}`。两次的 L1/L2 与分簇已有差别（sB 簇数 205 → 215），所以不是只差代码的对照。

| 旧 → 新 | sA | sB |
| --- | --- | --- |
| HS 结算 / 暂缓 | 290 / 42 → 239 / 81 | 275 / 55 → 237 / 97 |
| `error:ValueError` | 7 → 0 | 4 → 0 |
| `transverse_not_converged` / `walk_relax_unconverged` | 55 / — → 10 / 25 | 55 / — → 18 / 19 |
| `converged_candidate` | 9 → 6 | 6 → 3 |
| 对称合并后的连通事件类 | 2 → 2 | **2 → 1** |
| 第一轮脊类 HS 每项 forward | 1000 → 1280 | 1028 → 1198 |
| 跳支类 HS 总 forward（每项） | 412（37）→ 12 640（1264） | 357（51）→ 15 345（1395） |
| 第二轮总 forward（项数） | 82 936（87）→ 30 998（34） | 65 721（70）→ 20 926（22） |

- sB 少掉的一类事件（旧 `k5 / 10-11` 那一类）：该类组首这次 HS 仍出候选，但 HC `failed:tangent_no_resolved_drop`（旧运行里另一项 `k3 / 10-13` 也出过同样失败，这次反而连通）；旧运行中该类的三个备份项（`2-3`、`6-7`、`2-5`，当时都出候选并判 `hf_duplicate`）这次全部被排进第二轮并以 `HS_deferred:same_start` 跳过。即：组首 HC 失败 + 备份在第二轮 + 第二轮额度被第一轮挤掉。
- 折半重试的同次运行内证据（不受换机影响）：出现过失败试步的 HS 两块 37 / 44 项，**没有一项成为候选**；其中 12 / 24 条边在折半后继续前进，最后都以 hysteresis 结束；这些 HS 共花 65k / 77k forward（refine 的 27–32%）。`walk_relax_unconverged` 每项约 1.7–1.8k forward。（10-29 更正：65k / 77k 是这些 HS 的全部费用，不是重试的增量；同次运行的证据只说明“重试未观察到候选收益”，不能据此量化净损失或承诺节省比例。）
- 跳支括区修复后，跳支类 HS 真正计算了（11 / 11 项），但全部未出候选：`walk_relax_unconverged` 6 / 7、`above_ceiling` 3 / 3、`transverse_not_converged` 1 / 1。
- 小结：在 b05、refine 额度用满的条件下，两项修复后第一轮每项更贵、没有观察到新候选，第二轮被挤小；sB 少的一类事件来自“组首 HC 失败 + 备份在第二轮被暂缓”，其中备份被暂缓与第二轮额度变小有关，但取消重试能否找回它未验证。这只是 b05 一个体系、每块两类事件的结果。

## 已完成：延续试步未收敛时首次失败即终止该侧（HF-EF-branches-v4；2026-10-29；用户指示）

用户 10-29 原话：“改成‘未收敛试步首次失败，就结束该侧’。我之前建议的折半重试，在这批数据里没有显示出候选收益，不宜继续作为默认策略。”
- 代码（`explore/src/explore/ts_branch.py` `_walk`）：试步 `transverse_not_converged` 后不再折半，直接 `terminated` / `termination = "relax_unconverged"`，保留诊断（`unconverged_trials`、`relax_unconverged`），当前已接受点与曲率不变，另一侧照常延续；无最大值时总状态仍为 `unresolved:walk_relax_unconverged`。已收敛试步的几何跳变仍按原规则折半定位（`termination = "jump"`）。跳支括区修复保留。provenance `HF-EF-branches-v4`。中点阈值与调度不变。
- 测试：`test_ts_branch.py` 中原两个折半测试改为“首次失败即终止、只记一次失败试步、已接受 state 不被污染”；折叠势测试加断言“跳变一侧 `termination == "jump"`”。变异检验 M3（放回 v3 折半，临时目录 + `-o pythonpath=`）使 2 个新测试失败，真实源码对照 14/14 通过。ruff 通过；`test_ts_branch test_ts_ef test_ts_connect test_ts_warm_start test_d2_hf test_d2 test_d2_resume` **94 passed**；全包 `make check` 未跑。
- 对照工具：`benchmarks/explore_d2/legacy_walk_v3.py`（v3 `_walk` 冻结副本，只差折半那一处，已 diff 核对，仅供对照）；`benchmarks/explore_d2/hs_policy_compare.py`（冻结 HS 输入、同机 v3 / v4 各跑一次、每块按种子抽“出现过失败试步”与“对照”两组、逐项逐策略落盘可续跑）。

## 已完成：v3 / v4 冻结输入同机对照（183，2026-10-29；在 183 上读结果）

脚本 `benchmarks/explore_d2/hs_policy_compare.py`（运行）与 `hs_policy_compare_report.py`（汇总，只读）；结果 `benchmarks/explore_d2/results/2026-10-29/hs_policy_compare/full/`、`report.json`。输入：yayoi08 `b05_hf_v11r_s{A,B}` 的 HS 括区（冻结、只读）；每块按种子 20261029 抽 15 个“记录中出现过失败试步”的项与 3 个对照项，每项 v3 / v4 各跑一次，单 2080 Ti，16 线程（02:34 起，进程 54624）。72 个结果文件齐全，日志无报错。

| | 失败组（30 项） | 对照组（6 项） |
| --- | --- | --- |
| v3 / v4 总 forward | 52 628 / 26 772（省 25 856） | 3 968 / 3 740 |
| 本机 v3 实际遇到失败试步的项 | 20（这些项共省 19 408，平均约 970） | 0 |
| v3 未遇到失败试步的项（两策略同一代码路径） | 10（差 6 448，含机会与噪声） | 6 |
| 候选 `converged_candidate` v3 / v4 | 0 / 0 | 1 / 1（同一项） |
| `continuous_maximum` v3 / v4 | 1 / 2 | 1 / 1 |

- 同机复现噪声：对照组同一代码路径，v3 − v4 为 +9、−72、+49、+416、−171、−3 forward，状态全部相同。即单项 forward 差在约 ±400 以内不能归因于策略。
- 状态变化（失败组）：hysteresis → `walk_relax_unconverged` 4 项、`budget_exhausted` → `walk_relax_unconverged` 4 项、不变 20 项、hysteresis → `transverse_not_converged` 2 项（`sA_05` 像端点失败，未逐项查；`sB_12` 是 v4 找到连续最大值后停在二分中点。10-29 更正：先前写“2 项都是端点失败”不准确）。
- 连续最大值：`sB_14` 两策略都有，二分中点都没收敛（TODO 第 9 条）；`sB_12` 只有 v4 找到，同样停在中点。没有观察到 v4 丢失 v3 的最大值或候选。
- 结论：在这 30 个冻结输入上，v4 少花约一半 forward（25.9k / 52.6k），候选两边都是 0，只能说“未观察到候选丢失”，不能证明 v4 不会丢候选。省下的额度能否在整跑中换来更多覆盖或找回 sB 的事件，未验证。样本外推不做。

- 整理：这条线（诊断 → v3 → yayoi08 整跑 → v4 → 同机对照）的结果记录写成 [HF_WALK_POLICY_RESULTS_2026-10-29.md](HF_WALK_POLICY_RESULTS_2026-10-29.md)。

## 已核查：b05 v4 整跑（yayoi32，2026-10-29→30；离线读盘）

PBS `b05_hf_v4.o20601`：18:26:51 → 23:52:55（5.4 h，单张 4080 SUPER 两块共用）；`benchmarks/explore_d2/results/2026-10-29/b05_hf_v4_s{A,B}/` 报告齐全，账目闭合。HS 结算 267 / 272（yayoi08 v3：239 / 237），第二轮暂缓 60 / 58，候选 8 / 6，连通事件类 3 / 2。sB 两类已知 E2 都找到。sA 多出一个新鞍点（`k1 / distance:10-12`，−18880.41 eV，比已知两类低 0.39 eV，同一化学对 E2，与已知类对称 RMSD 约 1.40 Å）；它在三次运行里延续都已夹住最大值、0 次失败试步，前两次卡在二分中点，这次收敛，**与 v4 策略无关**，属 TODO 第 9 条那类。新鞍点无独立阶数与 QM 检验。三次运行硬件不同，v3 / v4 的整跑差别不能全归因于策略。详见 [HF_WALK_POLICY_RESULTS_2026-10-29.md](HF_WALK_POLICY_RESULTS_2026-10-29.md) §6b。

## 已完成：新 E2 鞍点的独立阶数检验（181，2026-10-30；只读 run，不改 journal）

脚本 `benchmarks/explore_d2/hf_saddle_validate.py`（复用 `hf_ef_bounded.py` 的 `check`：`confirm_minimum` 完整 Hessian 判阶数，一阶再用 `descend_saddle` 两侧下降；模型为运行时同一 fp32 POLAR-1；每个鞍点上限 4000 forward），结果 `benchmarks/explore_d2/results/2026-10-30/hf_saddle_validate/`。退出码 0。
- 对照（已知类 sA `k6 / distance:6-9`，−18880.0195 eV）：**一阶**，唯一负本征值 −2.877 eV/(Å²·amu)（884i cm⁻¹），内部模复核（步长 0.005 / 0.01 / 0.02 Å）一致；刚体模残差 1.8e-5。旧下降法一侧到极小（−18880.394 eV）、另一侧失败 → `descent_failed`（b05 上旧下降法以前就有此情况）。751 forward、89 s。
- 新鞍点（sA `k1 / distance:10-12`，−18880.4102 eV）：**阶数未判定**，`confirm_minimum` 返回 `trivial_mode_floor_above_tolerance`：刚体模残差最大 1.14e-3，超过 `trivial_floor_fraction` 0.5 × `minimum_check_eigenvalue_tol` 1e-3 = 5e-4。观察到唯一一个明显负本征值 −1.887（716i cm⁻¹），但第一个正内部本征值离 0 多远未输出，内部模复核未执行，所以“只有一个负模”未被证实。未做下降。104 forward、14 s。
- 这道门在 `src/prrs/runner.py:1588`（既有，非本轮新增），未改动。残差偏大的原因（推测与 fp32 Hessian 噪声有关）未验证。

## 已完成：新 E2 鞍点的内部模复核（181，2026-10-30；用户「跑」）

脚本 `benchmarks/explore_d2/hf_saddle_recheck.py`（直接调用 `prrs.runner.curvature_spectrum` 与 `internal_mode_recheck`，即 `confirm_minimum` 在刚体模门之后本会执行的复核；输出完整内部本征值谱；不推翻门的判定），结果 `benchmarks/explore_d2/results/2026-10-30/hf_saddle_recheck/`。退出码 0，每个鞍点 20 forward。
- **两次结果不同的原因**：`hessian_source = "auto"` 先算解析 Hessian，沿 3 个随机方向与有限差分交叉核验，相对偏差 > 0.01（`analytic_hessian_check`）就改用 fd Hessian（步长 0.01 Å）。上一轮 `confirm_minimum` 中新鞍点的交叉核验偏差为 1.086e-2 → 改用 fd → 刚体模残差 1.14e-3 → 被门拦下；本轮偏差 6.2e-3 → 用解析 Hessian → 刚体模残差 2.6e-5。即该交叉核验门限在此 fp32 体系上贴着噪声，同一几何可能走不同来源。门未改动，只报告。
- 新鞍点（sA `k1 / distance:10-12`，−18880.41 eV，当前几何 fmax 0.006 eV/Å）：解析 Hessian 下**只有一个负内部本征值 −1.952**；最低正内部本征值 0.0112、0.0271、0.0361、0.209（容差 1e-3）；内部模复核 `agree = true`、`resolution_stable = true`：负模沿模曲率 −1.62 / −2.03 / −2.00（步长 0.005 / 0.01 / 0.02 Å），第一个正模 0.025 / 0.013 / 0.014，三步长符号一致。
- 对照（已知类 sA `k6 / distance:6-9`）：负内部本征值 −2.862，最低正 0.0189、0.0220、0.179；复核一致。
- 结论：按 `confirm_minimum` 自身判据（解析 Hessian 交叉核验通过、刚体残差远低于门限、内部模复核一致），新鞍点为**一阶鞍点**。限定：这是单独调用得到的，不是一次完整通过的 `confirm_minimum`（正式那次因改用 fd 被门挡下）；第一个正模 0.011 是软模；旧 Hessian 下降未做（E/F 连接两侧 `ef_stable`）；无 QM 核验。

## 已完成：一阶性检验为什么把新鞍点拦下——解析 Hessian 交叉核验的噪声测量（181，2026-10-30）

脚本 `benchmarks/explore_d2/hessian_check_noise.py`（只读 run；fp32 POLAR-1），结果 `benchmarks/explore_d2/results/2026-10-30/hessian_check_noise/noise.json`。对象：新鞍点 sA `10-12` 与已知对照 sA `6-9`。
- 路径：`hessian_source = "auto"` → `_analytic_hessian`（`src/prrs/runner.py:649`）用 `config.seed` 播种的 3 个随机方向，比较解析 Hv 与中心差分 HVP（步长 `minimum_check_step_A` 0.01 Å），最坏相对偏差 > `analytic_hessian_check` 0.01 就改用 6N 逐坐标 fd Hessian（同一步长），再过刚体模门（`runner.py:1588` / `2655`）。
- 测量（每种步长重复 10 次，种子 42 与 43）：同一几何重复算力的最大抖动 0.0015 / 0.0027 eV/Å；解析 Hessian 两次调用相对差 1.4e-4 / 1.0e-4。交叉核验偏差随步长：0.0025 Å 中位约 0.033（100% 超 0.01）；0.005 Å 0.011–0.026（70–100% 超）；**0.01 Å（现行）中位 0.0065–0.012、最大 0.019，20–60% 超**；0.02 Å 0.004–0.0065（0–10% 超）；0.04 Å 0.003–0.005、几乎不抖（0% 超）。两个鞍点（含已判一阶的对照）表现相同。
- 解读：解析 Hessian 稳定；偏差在小步长由 fp32 力噪声主导（约按 1/h），0.04 Å 时只剩截断误差。现行步长下门限 0.01 落在噪声中，同一几何同一组方向约 20–60% 概率被判“解析不可信”，而回退得到的同步长 fd Hessian 在 fp32 下更差（新鞍点刚体残差 1.1e-3 对解析 2.6e-5），继而被刚体模门拦下。这道门比噪声还紧，挡住了正常检验（全局规矩 §8）。门在 PRRS 共享代码中，未改，记为 TODO 第 12 条，待用户决定。
- QM 准备：183 上 ORCA 6.1.1（`/home/korakuen/orca_6_1_1_linux_x86-64_shared_openmpi418_avx2/`，PATH 与 LD_LIBRARY_PATH 加 ORCA 目录，系统 Open MPI 5.0.11）2 核 H₂ 测试正常结束（20 s；`docs/experiments/b05_qm/orca_check/`，183 写入）。两个鞍点的 r2SCAN-3c OptTS 输入已备（`docs/experiments/b05_qm/ts_{new_10-12,ctrl_6-9}.inp`，几何取自 v4 sA journal），未运行。
- QM 启动（用户 10-30「可以」）：183 上两个 r2SCAN-3c `OptTS + Freq` 并行（各 16 核），18:57 启动，PID 67009（`ts_new_10-12`）/ 67010（`ts_ctrl_6-9`），记在 `docs/experiments/b05_qm/pids.txt`（首行混入一行 uptime 输出）；启动后 20 s 内两份 .out 均在增长。耗时未知（b04 6 原子 89 s；b05 16 原子未测）。结果须在 183 上读（NFS 缓存）。

## 已完成：新 E2 鞍点的 QM TS 优化（183，r2SCAN-3c，2026-10-30→31）

`docs/experiments/b05_qm/`（183 写入；协议同 `docs/experiments/b04_qm/`：r2SCAN-3c / TightSCF / `OptTS + Freq`，`Calc_Hess true`，总电荷 −1，各 16 核）。起点为 v4 sA journal 中的 POLAR-1 鞍点几何。
- 对照 `ctrl_6-9`：15 个优化循环收敛，34 min 36 s；**一个虚频 −529.0 cm⁻¹**（下一模 +86.8）；E = −693.783382 Eh；QM 优化位移（Kabsch，无置换）0.08 Å。
- 新鞍点 `new_10-12`：28 个循环收敛，47 min 57 s；**一个虚频 −127.1 cm⁻¹**（下一模 +75.2）；E = −693.799081 Eh，**比对照低 0.01570 Eh（约 0.43 eV / 9.9 kcal/mol，r2SCAN-3c，未补 ωB97M-V 单点）**；QM 优化位移 0.224 Å。
- 几何（QM TS）：新 C–Cl 2.015、Cβ10–H12 1.176、H12–O14 1.672 Å；对照 C–Cl 2.005、Cβ6–H9 1.282、H9–O14 1.385 Å。新 TS 更“早”。
- 虚频模组成（从 ORCA NORMAL MODES 读取）：对照以转移 H9 为主（振幅 0.89），沿模 Cβ–H / H–O / C–Cl = −0.81 / +1.01 / −0.13，近乎纯质子转移；新鞍点转移 H12（0.48）之外混有 OH 的 H15（0.37）与甲基 H（0.30–0.35），沿模 +0.16 / −0.56 / +0.19。两者都是 Cβ–H 与 H–O 反向、C–Cl 与 Cβ–H 同向的协同 E2 方向。
- 结论：新鞍点在 r2SCAN-3c 面上仍是一阶、E2 型、与已知类不同的 TS。未做：IRC（连接确认）、端点优化、ωB97M-V/def2-TZVPD 单点（能量按协议须以此为准）。
- QM 链启动（用户 10-31「可以 跑起来 183跑」）：`docs/experiments/b05_qm/run_chain.sh NAME 16`（IRC → 两端 `OPT Freq` → ωB97M-V/def2-TZVPD DEFGRID2 单点于 TS 与两端；每步输出正常结束即跳过，可续跑；任一步失败即停）。与 b04 协议唯一差别：IRC 初始 Hessian 读 TS `Freq` 的 `.hess`（同级别同几何），不再 `calc_anfreq` 重算。183 两链并行各 16 核，20:50:47 启动，PID 92004（`new_10-12`）/ 92005（`ctrl_6-9`），日志 `chain_*.log`，结果在 183 上读。

## QM 链改走第二种形式（183，2026-11-01；用户「可以」）

- 原链（`run_chain.sh`）跑了 9 h 13 min 仍在第一步 IRC：每步约 9.5 min；产物一侧是很平的长尾（新鞍点正向 57 步时比 TS 低 28.5 kcal/mol、梯度约 0.0016；对照反向 49 步时低 42.5 kcal/mol、梯度约 0.0028）。最后一帧几何：新鞍点正向 C–Cl 3.26 / Cβ–H 2.17 / H–O 0.98 Å，对照反向 3.29 / 2.31 / 0.97 Å，即 Cl⁻ 离开、水生成的 **E2 产物一侧**；对照正向 6 步结束于 1.94 / 1.18 / 1.60 Å（反应物一侧）。新鞍点的反应物一侧（反向）尚未走。
- 按 PID 核对后终止两条原链：链脚本 92004 / 92005、`orca` 92011 / 92012、`prterun` 235320 / 235662；之后两个会话无残留进程。
- 新链 `docs/experiments/b05_qm/run_chain2.sh NAME 16`（同级别同规则，可续跑）：新鞍点 = 只走反向 IRC（`Direction backward`，读 TS `.hess`）→ 产物端 `OPT Freq`（起点 = 正向 IRC 最后一帧 `prod_new_10-12_start.xyz`）→ 反应物端 `OPT Freq`（起点 = 反向 IRC 终点）；对照 = 反应物端 `OPT Freq`（`react_ctrl_6-9_start.xyz`）→ 产物端 `OPT Freq`（`prod_ctrl_6-9_start.xyz`）；两者最后对 TS 与两端做 ωB97M-V/def2-TZVPD DEFGRID2 单点。06:05:19 启动，PID 235826 / 235827，日志 `chain2_*.log`，结果在 183 上读。
- 改到 PBS yayoi14（用户 11-01：“不用给183 给这个节点”，“module load orca/6.1.1”）：新增 `docs/experiments/b05_qm/run_chain3.sh`（= `run_chain2.sh`，ORCA 路径改从 `PATH` 取；未改正在运行的 `run_chain2.sh`）与 PBS 脚本 `docs/experiments/b05_qm/b05_qm_pbs.zsh`（照模板只改 ← 行：yayoi14 纯 CPU 48 核、150 GB、48 h；`module load orca/6.1.1`，缺 mpirun 时补 `openmpi5.0.10`；两链并行各 24 进程，`OMP_NUM_THREADS=1`）。183 上 `run_chain2` 两链按 PID 停止（链脚本 235826 / 235827、`orca` 235832 / 235834、残留 `prterun` 238222），两个会话无残留。183 那 17.5 min 未完成任何一步（反向 IRC 与对照反应物端优化），PBS 上从头重跑这两步。
- yayoi14 第一次提交（作业 22002）13 s 即结束：两链第一步都 `FAILED`。原因是 Open MPI（模块 `openmpi/5.0.7_gcc13.3.0`）报 “not enough slots … 24 slots requested”：select 只写 `ncpus=48`，PBS 默认 `mpiprocs` 为 1，只给 1 个 MPI 槽位。修正：select 加 `mpiprocs=48`。（我先前说的“模块定义的是 zsh 别名、bash 找不到 orca”是错的：ORCA 实际已经启动，`/home/apps/orca_6_1_1_avx2/orca`。）`libtorque.so.2` 的 dlopen 警告日志中注明为 ignored。失败的输出会在重提交时覆盖。

## QM 链结果（yayoi14 作业 22006，2026-11-01 06:40→09:19；在 183 上读，183 看到的文件为新）

`docs/experiments/b05_qm/summarize.py`（标准库解析）→ `summary.json`（183 写入）。两条链每一步 ORCA 均正常结束，但：
- **三个端点几何优化未收敛**：新鞍点反应物端、新鞍点产物端、对照产物端都停在 ORCA 默认 50 个循环上限（RMS 梯度已达标，最大梯度 / 步长未达标，能量每循环仍降 1e-5～5e-5 Eh，松散复合物在平坦势能面上漂移），因此没有算频率；只有对照反应物端收敛（0 个虚频，最低 64.7 cm⁻¹）。这是链脚本的漏洞：`done_ok` 只看 “ORCA TERMINATED NORMALLY”，优化未收敛时 ORCA 也正常结束，链照常往下走并对未收敛几何做了单点。
- 新鞍点反向 IRC 正常收敛（7 点）。
- 连接（按几何）：新鞍点反应物端 C–Cl 1.95 / Cβ–H 1.12 / H–O 2.05 Å（OH⁻···H–C 复合物），产物端 C–Cl 3.97 / Cβ–H 2.26 / H–O 0.97 / C0–Cβ 1.342 Å（Cl⁻ 离开、水生成、C=C 形成）；对照两端同型（反应物 1.91 / 1.12 / 1.96；产物 4.01 / 2.64 / 0.96 / 1.338）。即两个 TS 都连 E2 反应物复合物与 E2 产物。
- **可靠的能量**（两个 TS 都已收敛且一阶）：TS_new − TS_ctrl = **−11.75 kcal/mol（ωB97M-V/def2-TZVPD）**，−9.85（r2SCAN-3c）。
- **暂定（依赖未收敛端点）**：势垒（TS − 各自反应物复合物）新 +4.07 / 对照 +7.59 kcal/mol（ωB97M-V；r2SCAN +0.57 / +2.35）；反应能新 −30.88 / 对照 −47.61（ωB97M-V）；反应物复合物 new − ctrl = −8.22（ωB97M-V）。即新 TS 低的 11.75 kcal/mol 中约 8 kcal/mol 来自它的反应物复合物本身更低（不同的 OH⁻ 位置），按暂定值势垒差约 3.5 kcal/mol。
- 补跑三个未收敛端点（用户 11-01「可以」，节点 yayoi41）：`docs/experiments/b05_qm/run_endpoint2.sh TAG NPROCS`，从上次最后几何 `ep_<TAG>.xyz` 接着 `OPT Freq`（`%geom MaxIter 300`），**优化步须见 “OPTIMIZATION RUN DONE” 才算完成**，未收敛即 `FAILED`；随后 ωB97M-V 单点；新文件名 `ep2_*` / `wb_ep2_*`，不覆盖旧输出。PBS 脚本 `docs/experiments/b05_qm/b05_qm_ep2_pbs.zsh`（照 `b05_qm_pbs.zsh` 只改 ← 行；按用户建议改到纯 CPU 节点 yayoi13：ncpus 48 / mpiprocs 48 / 150 GB / 24 h；三个端点并行各 16 进程）。

## 已完成：b05 新 E2 TS 的 QM 核验收尾（yayoi13 作业 22053，2026-11-01 12:04→12:40）

三个端点用 `run_endpoint2.sh` 重新优化（须真正收敛才算完成）：新鞍点反应物端 17.5 min、产物端 35 min、对照产物端 12.5 min，全部 “OPTIMIZATION RUN DONE”，频率 0 个虚频（最低 96.0 / 28.3 / 25.5 cm⁻¹）；对照反应物端原已收敛（最低 64.7 cm⁻¹）。`docs/experiments/b05_qm/summarize.py` 改为优先取 `ep2_*`，输出 `summary.json`。

| 量（电子能，无 ZPE / 热校正） | 新 TS `10-12` | 对照 TS `6-9` |
| --- | --- | --- |
| 虚频（r2SCAN-3c） | −127.1 cm⁻¹ | −529.0 cm⁻¹ |
| 势垒 = TS − 各自 IRC 反应物复合物，ωB97M-V/def2-TZVPD（r2SCAN-3c） | **+4.08**（+0.57）kcal/mol | **+7.59**（+2.35）kcal/mol |
| 反应能 = 各自 IRC 产物复合物 − 反应物复合物，ωB97M-V（r2SCAN-3c） | −31.39（−34.65） | −47.61（−51.16） |

- TS_new − TS_ctrl = **−11.75 kcal/mol（ωB97M-V）**、−9.85（r2SCAN-3c）；POLAR-1 面上两者差 0.39 eV ≈ 9.0 kcal/mol。反应物复合物 new − ctrl = −8.23（ωB97M-V）。即新 TS 更低的 11.75 kcal/mol 中约 8.2 来自其反应物复合物更低（OH⁻ 位置不同），按各自复合物计的势垒差约 3.5 kcal/mol。
- 连接：两个 TS 的 IRC 两端都是 E2 反应物复合物（OH⁻···H–Cβ）与 E2 产物复合物（Cl⁻ 离开、H₂O、C=C 1.34 Å），端点均为收敛极小。
- 限定：势垒以各自 IRC 连到的局部复合物为参照，两者参照态不同；没有统一参照（分离反应物或全局最低复合物），没有 ZPE / 热校正。两个产物复合物的排布不同，反应能之差不代表反应热力学差异。
- 结论：sA `k1 / distance:10-12` 在 r2SCAN-3c 面上是独立的一阶 E2 TS，IRC 连通 E2 反应物与产物，且比此前已知的 E2 TS 低 11.75 kcal/mol（ωB97M-V 单点）。

- 交接：本工作线（10-26 至 11-01）的交接文档 [HANDOFF_HF_D2_2026-11-01.md](HANDOFF_HF_D2_2026-11-01.md)；提交草稿补 `2026-11-01_1400_b05-qm-and-saddle-checks`。

## 已完成：二分中点失败的按阶段诊断（TODO 第 9 条；2026-11-01；只读 journal，不调模型；用户「可以」）

- 脚本修正 `benchmarks/explore_d2/hs_transverse_diag.py`：成功判据按 `_relax` 实际用的阈值分阶段：二分中点（refine 的 `from_left` / `from_right`）取 manifest `hf.midpoint_transverse_fmax_eV_A`（b05 为 1e-3），refine 两端重松弛（`left` / `right`）与所有延续 / 端点松弛取 `EFSettings.transverse_fmax_eV_A`（0.003）。修正了 10-28 版统一用 0.003 的错误。新增每段末 50 步的统计：横向力分位、落在 1.5× / 2× 阈值内的帧数、每步最大原子位移、能量极差（fp32 量子）、反力极差。
- 输入 b05 HF 八块：10-26 v1 s{A,B}、10-28 v1.1 s{A,B}、10-28 v1.1r（yayoi08）s{A,B}、10-29 v4（yayoi32）s{A,B}。结果 `benchmarks/explore_d2/results/2026-11-01/hs_transverse_diag/diag.json`、`diag.log`（10-28 的旧 `diag.json` 保留不动）。
- 二分中点松弛（八块合计）：**197 / 268 次收敛到 ≤ 1e-3**；收敛的步数中位 12–26，p90 多为 24–34，≥ 150 步的只有 5 次。失败 71 次（前四块 43 次，与 10-28 一致），全部用满 200 步。refine 两端重松弛 0 次失败。
- 失败中点的末 50 步（71 段）：
  - 过程最低横向力为阈值的 1.00–1.75 倍（中位 1.23）；64 / 71 段的最低点出现在末 50 步之前，末段反而回升。
  - 末段横向力中位数在各段间为 0.0033–0.032 eV/Å（中位 0.0096，即阈值的约 10 倍）；3550 帧里只有 16 帧在 1.5 倍阈值内、79 帧在 2 倍内。
  - 末段每步最大原子位移中位 0.011 Å，p90 中位 0.0495 Å，贴着 BFGS `max_step_A` 0.05 Å 上限。
  - 同一 q 上末段反力极差中位 0.065 eV/Å（0.011–0.18），是冻结 σ（0.0025）的约 26 倍；整段首末对齐 RMSD 0.0014–0.085 Å；能量变化 −5 至 +1 个 fp32 量子（0.00195 eV）。
- 对照的噪声量（不同几何，只作量级参考）：10-30 在两个 b05 鞍点上同一几何重复算力，最大分量偏差 0.0015 / 0.0027 eV/Å，逐分量标准差的最大值 0.0007 / 0.0014（`results/2026-10-30/hessian_check_noise/noise.json`）。
- 读法：失败中点的末段不是停在噪声底的平台。横向力是阈值的约 10 倍、高于测得的力抖动；优化器还在走接近上限的大步，反力随之大幅变化，即在固定 q 上来回游走。阈值 1e-3 本身是可达的（74% 的中点达到，多数 35 步以内），把步数上限加大也看不出能挽救（收敛的几乎都在 150 步以内）。所以"阈值贴着噪声"这个说法不被支持，放宽阈值没有依据。游走的原因本轮没有区分（见 TODO 第 9 条未确认部分）。
- 补充（同日，只读）：71 次中点失败只来自 22 个意图，同一意图在不同硬件和代码版本上反复卡在中点（如 `k0 / distance:2-4` 八块全部失败）；这 22 个意图走到中点共 72 次，失败 71 次，唯一例外是 v4 sA `k1 / distance:10-12`（新 E2 鞍点）。下一步的对照计划：[PLAN_HF_MIDPOINT_COLD_2026-11-01.md](PLAN_HF_MIDPOINT_COLD_2026-11-01.md)（待用户说"跑"）。

## 已完成：二分中点"带入曲率 vs 同帧冷启动"同机对照（183，2026-11-01 19:26→20:53；用户「183 跑起来」）

[计划](PLAN_HF_MIDPOINT_COLD_2026-11-01.md)。脚本 `benchmarks/explore_d2/midpoint_cold_replay.py`（A/B/C/D 四组，逐项落盘可续跑），启动脚本 `benchmarks/explore_d2/midpoint_cold_run.sh AB|CD`。
- 加载核对（本机，无模型）：八块读出失败中点 71、已收敛中点 197、涉及 22 个意图，与 11-01 诊断一致。
- 冒烟（183，每组 1 项、重复 1 次，`results/2026-11-01/midpoint_cold/smoke/`）：两种模式端到端跑通。A：v4 sB `k0 / distance:10-12` 整条 HS 399 forward（延续 196 + 搜索 203，账本闭合），183 上重现了中点失败；B：同帧冷启动 200 步仍未收敛（201 forward）；C：1 项未收敛；D：1 项 16 步收敛。B 的第 0 帧与 A 的第 0 帧几何相同，能量差 1 个 fp32 量子，横向力 0.0780 对 0.0785（GPU fp32 重复求值差）。实测 0.090–0.106 s/forward；183 读入八块 journal 约 7.5 min（首次读，之后有缓存）。
- 正式运行：19:26 启动 AB（python PID 241912）、19:29 启动 CD（PID 242006），各 16 线程，记在 `results/2026-11-01/midpoint_cold/full/pids.txt`，结果在 183 上读。
- 操作注意：`ssh ... "bash -c 'cd X && setsid nohup CMD > log 2>&1 < /dev/null & echo \$!'"` 这种写法里，`&` 作用于整个 `cd && ...` 列表，派生出的子 shell 继承了 ssh 的输出，ssh 不返回，`$!` 拿到的也是子 shell 的 PID。改成 `cd X && { setsid nohup CMD > log 2>&1 < /dev/null & echo \$!; }` 就正常返回，拿到的是 python 本体的 PID。
- 结束：AB 20:24、CD 20:53 完成（11-02 09:16 在 183 上核对：两进程已退出，A 22/22、B 84/84、C 213/213、D 90/90 文件齐全，日志无 error / traceback）。共 73.9k forward（A 13 615、B 14 883、C 39 455、D 5 981）。汇总为只读读取，在 183 上做。
- A（现行代码，183 整条 HS 重跑）：22 个意图中 21 个仍以 `transverse_not_converged` 结束，1 个为 `walk_relax_unconverged`（没走到中点）；183 上共 28 次中点松弛，8 次收敛、20 次失败。原失败在同机重现。
- B（同机、同一起始帧配对）：带入曲率失败的 20 次，冷启动 3 次中收敛 0/3 的 16 次、1/3 的 2 次、2/3 的 1 次、3/3 的 1 次；带入曲率收敛的 8 次，冷启动 1/3 的 6 次、0/3 的 2 次。
- C（记录的 71 个失败帧，冷启动 ×3）：213 次中收敛 24 次（11%）；53 帧 0/3、13 帧 1/3、4 帧 2/3、1 帧 3/3。收敛的步数中位 42.5、p90 154。失败段末 50 步横向力中位 0.014 eV/Å、每步位移中位 0.012 Å，过程最低横向力 0.00101–0.0024，与原失败同样是在固定 q 上游走。按意图：k3 的几个意图较常收敛（`k3 / distance:8-14` 8/12、`4-14` 4/9），k0、k2 几乎全为 0；新鞍点所在 `k1 / distance:10-12` 3/12。
- D（30 个原本收敛的帧，冷启动 ×3）：90 次中收敛 70 次（78%）；23 帧 3/3、1 帧 1/3、6 帧 0/3；收敛步数中位 27.5、最多 82，与原记录步数相近。
- 读法：冷启动在原本收敛的帧上收敛率 78%，在原本失败的帧上 11%，差别跟着起始帧走，不跟着曲率来源走；带入曲率在 B、D 里反而略有帮助。所以原因 ①（带入曲率不适合中点）不成立。71 帧中 18 帧三次重复结局不一致，说明这部分帧的结局受 GPU fp32 不确定性影响；多数失败帧稳定失败。原因 ②（噪声污染曲率更新）与 ③（中点处横向很软或不是极小）本实验未区分。


## 已完成：失败中点与收敛中点的横向曲率（TODO 第 9 条；2026-11-02，PBS yayoi35 作业 22656，10:31→10:55；用户「冒烟结束后给这个节点」= yayoi41；yayoi41 排队，用户改到 yayoi35）

[计划](PLAN_HF_MIDPOINT_TRANSVERSE_2026-11-02.md)。脚本 `benchmarks/explore_d2/midpoint_transverse_hessian.py`，PBS 脚本 `benchmarks/explore_d2/midpoint_transverse_hessian_pbs.zsh`（照模板只改 ← 行：yayoi35（原 yayoi41 排队），ncpus 32 / ngpus 1 / 56 GB / 24 h），结果写 `results/2026-11-02/midpoint_transverse_hessian/full/`。
- 本机核对：本机 torch 2.12.0 / CUDA 12.9 的架构列表含 sm_120（RTX 5080）。
- 冒烟（183，F、S 各 1 帧，`results/2026-11-02/midpoint_transverse_hessian/smoke/`，09:26→09:37，其中读入 journal 约 11 min）：读出 F 71、S 197，与诊断一致；每帧 16–34 s（解析 Hessian 15–28 s）；横向 41 个本征值（16 原子，3N−7）。有限差分抽查（0.04 Å）偏差 0.0018–0.011，Hessian 不对称度 4e-5–7e-5。F_000 最低横向本征值 0.016、S_000 0.061 eV/Å²（各 1 帧，不作比较）。
- 冒烟后修了一处输出：游走位移原来没去掉刚体漂移（F_000 未对齐 1.38 Å，几乎全是整体移动，横向占比 1.4%）；改成先 Kabsch 对齐到该帧再分解，离线在同一帧核对对齐后为 0.099 Å。Hessian 部分未改，没有重新冒烟。
- 耗时依据：183 上每帧约 16 s × 268 ≈ 72 min，加读入约 11 min；5080 上未测。
- 运行：PBS `b05_midT.o22656`（项目根目录）：yayoi35，RTX 5080，10:31:26 → 10:55:01，268 帧全部完成，日志无报错；首帧 40 s，之后每帧约 3.2 s。在 183 上读结果（183 上 268 个文件齐全）。汇总脚本 `benchmarks/explore_d2/midpoint_transverse_summary.py`（只读，并按 run / task / q 与 11-01 冷启动 C 组对上 71/71），输出 `results/2026-11-02/midpoint_transverse_hessian/summary.json`。
- 解析 Hessian 抽查：0.04 Å 有限差分偏差最大 0.014（中位 0.004–0.005），不对称度 ~1e-4 以下。
- 失败中点 F（71 帧，取横向力最小帧，横向力 0.0010–0.0018）：横向最低本征值全部为正，0.013–0.024 eV/Å²（中位 0.0187），次低 0.113–0.153，二者之比 4.9–9.0，即有一个与其余明显分开的软横向模。对齐后的游走位移 0.001–0.34 Å（中位 0.068），100% 在横向子空间内，**落在最低一个横向模上的比例中位 96%**（p10 68%），最低 3 个模合计中位 98%。
- 收敛中点 S（197 帧，末帧）：横向最低本征值 0.005–0.123（中位 0.060），没有负值。按意图分两类：与 F 同一批意图（如 `k0/k1/k3/k7 · distance:10-12`、`k3 · 2-4 / 6-8 / 8-14`）的收敛中点最低值也是 0.016–0.021，与 F 相同；其余意图约 0.06；另有几个意图（`k4 · 3-7 / 7-11`、`k6 · 3-13 / 9-13`）低到 0.005–0.007 仍然收敛。
- F 与冷启动结局的关系（同一帧）：冷启动 0/3 的 53 帧最低值中位 0.0179，1/3 的 13 帧 0.0195，2/3 的 4 帧 0.0206，3/3 的 1 帧 0.0218。样本小，只是方向一致。
- 读法：失败中点都不是横向鞍点（无负曲率），而是有一个单独的软横向模，游走几乎全部沿着这个模。但同样软、甚至更软的收敛中点也存在，所以"软"本身不足以决定失败。原因 ③ 只得到部分支持，② 没有被排除；两者可能是同一件事（见 TODO 第 9 条未确认部分）。软模的原子组成没有输出（脚本没存本征向量）。


## 已完成：中点阈值放宽的离线重放与改法选项（TODO 第 9 条；2026-11-02；只读 journal，不调模型；用户「可以 这个部分可以仔细研究」）

- 脚本 `benchmarks/explore_d2/midpoint_threshold_replay.py`，结果 `results/2026-11-02/midpoint_threshold_replay/replay.json`（八块 b05 HF 合并）。依据：`_relax` 遇到第一帧横向力 ≤ τ 就返回，所以放宽到 τ' 时停的帧就是记录轨迹中第一帧 ≤ τ' 的帧。
- 原本收敛的 197 段、92 对：放宽到 2e-3，反力与 1e-3 结果之差中位 0.0013、p90 0.0038；成对反力差超门（0.0106）由 1 对增至 2 对，3e-3 时 6 对；成对符号类别改变 13 / 92（2e-3）、22 / 92（3e-3），几乎都在近根处。1 对在 2e-3 时反力 −0.025 / −0.036，继续松弛 138 / 156 步后变为 +0.42（长松弛中换支）。
- 原本失败的 71 段：2e-3 下 71 段都有帧达到（中位第 6 步）；**65 段在所有横向力 ≤ 2e-3（或 ≤ 3e-3）的帧上反力同号且 |RF| > 0.003**，括区方向其实已定（多为第一个中点，|RF| 约 0.1）；6 段在近根处跨过 ±0.003。新 E2 鞍点所在意图 `k1·10-12` 的 4 次失败都在那 65 段里。
- 代码阅读：中点阈值只决定中点松弛何时算完成；候选判据另用横向 0.003、完整力 0.01（`ts_ef.py:350`）。
- 选项文档：[HF_MIDPOINT_THRESHOLD_OPTIONS_2026-11-02.md](HF_MIDPOINT_THRESHOLD_OPTIONS_2026-11-02.md)（O1 不改 / O2 统一放宽 / O3 分级终止 / O4 软模感知 / O5 换优化器），等用户定。

## 已完成：二分中点分级终止（O3，HF-EF-bracket-v4）施工、界线测量与同机验收（TODO 第 9 条；2026-11-02；用户决定见计划）

[施工计划](PLAN_HF_MIDPOINT_COARSE_2026-11-02.md)（含用户 11-02 原文）。
- 代码 `explore/src/explore/ts_ef.py`：`EFSettings` 新增 `midpoint_coarse_transverse_fmax_eV_A`、`midpoint_sign_bound_eV_A`（要么都给、要么都不给；粗阈值须大于中点阈值；界线不小于 `reaction_force_tol`）。启用时每个中点两侧先松弛到粗阈值（origin `from_left_coarse` / `from_right_coarse`），`coarse_sign` 依次检查：括区宽度、几何（`branch_rmsd_A`、`branch_step_A`）、双侧反力差门、两侧都不满足候选判据、两侧反力同号且都超过界线。全部通过才按符号收括区；任一不过就从粗阶段的结构和曲率继续到中点阈值，之后完全走 v3 逻辑。约束残差仍在 `_relax` 内检查。启用时 provenance 为 `HF-EF-bracket-v4`，不启用仍为 v3。D2 / `HFSettings` 未接。
- 测试 `explore/tests/test_ts_ef.py` 新增 9 项：设置校验 4 组；界线不可达时结果与 v3 相同、每个中点仍有精阶段；可分辨时在粗阶段收括区、候选与 v3 相同、精阶段次数减少（解析面上松弛只要几步，粗阶段自身的求值可能抵消节省，费用不在解析台断言）；粗阶段反力差超门时转精阶段而不报分歧；粗阶段一侧已满足候选判据时不提前收；精阶段失败仍为 `transverse_not_converged`。变异检验（临时目录 + `-o pythonpath=`）：去掉粗阶段反力差门、去掉候选让位，各使 1 个新测试失败。`test_ts_ef test_ts_branch test_ts_connect test_ts_warm_start test_d2_hf test_d2 test_d2_resume` **103 passed**；全包 `make check` 未跑。
- 粗、精反力偏差的离线估计（journal 中 197 个收敛中点；同一次松弛里第一帧 ≤ 2e-3 与最终 ≤ 1e-3 的反力之差；去掉 2 个换支离群 |差| ≈ 0.45）：带符号均值 3e-5、中位 0、标准差 0.0023；按 |RF| 分层均值 +0.0004 / −0.0002 / −0.0005（标准差 0.0018 / 0.0025 / 0.0030）；朝 |RF| 增大方向的均值 −1e-4。未见一致偏差，差的分散与 σ 同量级。限定：只含最终收敛到 1e-3 的段（有选择偏差），带入曲率的松弛。
- 新脚本：`benchmarks/explore_d2/midpoint_sign_resolution.py`（界线测量：同点 8 次扰动重复松弛，一次松弛同时给出粗、精两阶段的反力；收敛中点按 |RF| 三层各 8 个 + 失败中点每意图 1 个）；`benchmarks/explore_d2/midpoint_coarse_validate.py`（冻结括区新旧同机对照：22 个失败意图 + 按种子抽的原本出候选的对照项，每项 v3 / v4 各跑完整 HS）。`midpoint_cold_replay.py` 的读入加了 `recorded_rf` 字段。
- 冒烟（183，12:15→12:23，`results/2026-11-02/midpoint_coarse_smoke/`，PID 244774，脚本 `smoke.sh`；183 CPU 负载 15–25 来自别人，GPU 空闲）：两个脚本端到端跑通。验收脚本（临时界线 0.0105，仅为跑通，非定值）：失败项 `k0·10-12`（v4 sB 块）v3 405 forward 仍卡第一个中点；v4 442 forward，前三个中点在粗阶段收括区（|RF| 0.11 → 0.046 → 0.014），第四个中点 q 1.3232 处粗阶段反力 +0.0003 / +0.0041 未达界线，转精阶段后 200 步未收敛（精阶段反力 −0.013）。即该项卡点后移到近根。对照项：v3 260 / v4 252 forward，候选能量差 0、RMSD 3e-4 Å。界线测量脚本：S 点 1 次重复 51 forward（粗 RF 0.0036，第 37 步；精 0.0050，第 50 步）；F 点 1 次重复 201 forward（粗 RF −0.115，第 54 步；精未收敛）。
- 界线测量正式运行（用户 11-02「183丢过去 让他自己跑去」）：183，12:55 启动 S（PID 245207）与 F（PID 245215），各 16 线程，启动脚本 `benchmarks/explore_d2/midpoint_sign_resolution_run.sh S|F`，结果 `results/2026-11-02/midpoint_sign_resolution/`（`pids.txt` 同目录；在 183 上读）。启动前 12:53 实测 183 负载 0.00、GPU 空闲。
- 界线测量结果（S 13:42、F 13:59 结束；46 点 × 8 次 = 368 个结果，日志无报错，共 55.5k forward；在 183 上读）。汇总脚本 `benchmarks/explore_d2/midpoint_sign_resolution_summary.py` → `results/2026-11-02/midpoint_sign_resolution/summary.json`。注意：183 上列这个目录时每个文件名重复 8 次（NFS 读目录异常，文件只有一份），汇总按文件名去重；在 183 上 glob 读结果都要去重。
  - σ_c（同点 8 次重复的粗阶段反力标准差，各点中位 / p90 / 最大，eV/Å）：收敛中点 |RF| < 0.01：0.0032 / 0.0043 / 0.0062；0.01–0.05：0.0041 / 0.0052 / 0.0065；> 0.05：0.0060 / 0.021 / 0.033（最大值来自换支点 `k4·7-11`：journal 中 +0.42，粗阶段均值 −0.038；去掉后全体最大 0.010）；失败中点：0.0040 / 0.0075 / 0.0100。全部 46 点中位 0.0040、p90 0.0073、合并 0.0068。24 次重复 200 步内未达到 2e-3（S_18、S_20 各 7 次，F_06 5 次，F_08 4 次，S_23 1 次）。
  - 粗、精偏差（精阶段收敛的 140 次重复）：均值 −0.0006、中位 −0.0002、标准差 0.0047、5–95% 为 −0.0088 至 +0.0042，无 |差| > 0.02。未见一致偏差，分散小于 3σ_c。限定：精阶段收敛的样本偏向小 |RF| 的点（|RF| > 0.05 层只有 12 次）；冷启动 BFGS 加 0.01 Å 扰动，与实际运行带入曲率不同。
  - 代入 B = 0.003 + 3σ_c：σ_c 取中位 → 0.015、合并 → 0.023、p90 → 0.025。三者下，粗阶段超过 B 的精收敛重复（62 / 40 / 37 次）中，精阶段符号相反的都是 0；没有点出现两次重复以相反符号超过 B；失败中点中 8 次重复全部达到 2e-3 且都超过 B、同号的有 16 / 22（另有 2 个点的未达到 2e-3 的重复除外后也全部超过；此前口头报的 18 / 22 是这个较宽的数法）。
- **用户决定（11-02「可以」，回答"B 取 0.025，可以吗"）：B = 0.025 eV/Å**（σ_c 取各点 p90 0.0073）。
- 验收正式运行（用户 11-02「183跑起来」）：183，14:33 启动，PID 247165（`results/2026-11-02/midpoint_coarse_validate/pids.txt`），16 线程，`benchmarks/explore_d2/midpoint_coarse_validate_run.sh`（B 0.025，22 失败意图 + 22 对照，各 v3 / v4 完整 HS）。启动前负载 0.18、GPU 空闲；`ts_ef.py` md5 76e3fc89…、`midpoint_coarse_validate.py` dc92ce6c…（两边相同）。结果在 183 上读。
- **验收结果**（183，14:33→15:39，74 次运行全部完成，日志无报错；18:51 核对进程已退出；在 183 上读）。汇总脚本 `benchmarks/explore_d2/midpoint_coarse_validate_summary.py` → `results/2026-11-02/midpoint_coarse_validate/summary.json`。
  - 与计划不同：对照只有 15 项（原本出候选且不在失败名单里的不同意图共 15 个，脚本全部取用），不是 22 项。
  - 对照 15 项：v3、v4 全部 `converged_candidate`，候选能量差 0 或 −0.002 eV（1 个 fp32 量子），Kabsch RMSD ≤ 0.0017 Å；forward 4 225 / 4 368。对照项的粗阶段都未提前收括区。
  - 失败 22 项：v4 救回 2 项为 `converged_candidate`（`k1·8-14` −18880.4219 eV、`k3·4-14` −18880.4082 eV，与新 E2 鞍点 −18880.41 同一能量段，比已知两类低约 0.4 eV；**是否为同一鞍点、两侧连到哪里都未查，本次只跑 HS 未跑 HC**）；16 项两版都是 `transverse_not_converged`；另 4 项的状态变化发生在两版相同的端点阶段（`no_resolved_maximum_bracket` ↔ `transverse_not_converged`，以及 v3 的一次 `branch_force_disagreement`），属同机 fp32 不可复现。forward 13 977 / 14 167，两组合计 18.2k / 18.5k，成本持平。
  - 失败项的 v4 粗阶段：提前收括区 27 次；未收 22 次，其中 20 次两侧 |RF| 为 0.0002–0.026、低于 B（多数 0.01–0.02，两侧同号且互差 ≤ 0.004），2 次两侧都超过 B 且同号，但反力差 0.0117、0.0178 超过门 0.0106 被挡下。未收之后转精阶段，在精阶段失败，即卡点后移到 |RF| < B 的近根中点。
  - 门的观察（只报告，未改）：反力差门用冻结 σ 0.0025 得 0.0106；粗阶段实测 σ_c 中位 0.004、p90 0.0073，两侧之差的标准差约 √2·σ_c ≈ 0.006–0.010，按此门在粗阶段经常触发。门按用户决定保留。
  - 小结：v4 不改变对照组候选，成本持平，救回 2 / 22；剩下的失败集中在近根中点的精阶段松弛。按用户 11-02"若剩下的失败仍集中在近根松弛，再做 O5 的优化器对照"，下一步是 O5（未开始）。

## 已完成：近根中点精松弛的优化器对照（O5，TODO 第 9 条；2026-11-02 19:27→21:19；用户「慢慢跑 你看看183 空吗」）

[计划](PLAN_HF_MIDPOINT_O5_2026-11-02.md)。脚本 `benchmarks/explore_d2/midpoint_optimizer_compare.py`（A0 BFGS 0.05 Å / A1 BFGS 0.01 Å / A2 FIRE；FIRE 通过在 `ts_ef` 内临时替换 `BFGS` 符号接入，`_relax` 其余语义不变），启动脚本 `midpoint_optimizer_compare_run.sh OUT LIMIT REPS_N`；`midpoint_cold_replay.py` 的 `cold`、`rerun` 加了可选 `settings` 参数（默认不变）。
- 183 负载：19:19 有 THERMO 项目约 10 个单线程作业（负载 10–14），GPU 空闲；按"慢慢跑"只开 1 个进程、8 线程，N → F → S 依次。
- 冒烟（每组 1 项、每臂 1 次，`results/2026-11-02/midpoint_optimizer_compare/smoke/`，19:2x 跑完）：N 组 HS 重跑 423 forward，精阶段卡点 1 个；该帧 A0、A1 200 步未收敛，A2（FIRE）8 步收敛。F 组 1 帧：A0、A1 未收敛，A2 79 步收敛。S 组 1 帧：三臂都收敛（31 / 9 / 16 步）。各 1 帧，不作结论。
- 正式运行：19:27 启动，PID 267531（`results/2026-11-02/midpoint_optimizer_compare/pids.txt`），输出 `full/`，在 183 上读。
- **结果**（11-02 21:19 结束；11-04 12:59 在 183 上核对：进程已退出，439 条日志、417 次优化器运行 + 22 次 HS 重跑齐全，无报错）。汇总脚本 `benchmarks/explore_d2/midpoint_optimizer_compare_summary.py` → `results/2026-11-02/midpoint_optimizer_compare/summary.json`（按文件名去重）。

| 组（收敛 / 次数，forward） | A0 BFGS 0.05 Å（现行） | A1 BFGS 0.01 Å | A2 FIRE（默认参数） |
| --- | --- | --- | --- |
| N 近根卡点（19 帧 × 2 次） | 3 / 38，7 071 | 5 / 38，6 789 | **33 / 38，2 244** |
| F 远根失败中点（71 帧） | 3 / 71，13 799 | 4 / 71，13 751 | **54 / 71，8 127** |
| S 原本收敛的对照（30 帧） | 23 / 30，2 050 | 24 / 30，1 889 | **28 / 30，1 755** |

  - FIRE 收敛步数中位 19（N）/ 73（F）/ 35（S），最多 198；缩小 BFGS 步长（A1）只多收敛 1–2 次，问题不只是步长过冲。
  - 收敛时的反力：同一帧两次重复都收敛的 17 对（各臂合计）相差中位 0.0010、最大 0.0035，在 σ_c（中位 0.004）以内；S 组 A0 与 A2 都收敛的 23 帧相差中位 0.0010、最大 0.0039，与 A0 对 A1（中位 0.0012、最大 0.0040）同量级。无独立真值核对。
  - FIRE 仍未收敛：N 组 `k0·8-14`、`k2·6-8` 两帧两次都未收敛，`k1·2-4` 一次未收敛；F 组 17 / 71。
  - N 组为取得卡点几何重跑的 22 次 v4 HS：2 次 `converged_candidate`，19 次 `transverse_not_converged`，1 次 `walk_relax_unconverged`；与 11-02 验收的 2 个救回项是否同为那两个意图，未核对。
  - 限定：各臂冷启动（实际运行带入曲率）；FIRE 默认参数未调；只回答"能否收敛"，把 FIRE 接进中点松弛后整条 HS 的救回数、候选一致性和成本未测。

- 交接：本工作线（11-01 至 11-04）的现状与待办写成 [HANDOFF_HF_MIDPOINT_2026-11-04.md](HANDOFF_HF_MIDPOINT_2026-11-04.md)（用户 11-04「你把当前东西落盘一下」）；提交草稿 `2026-11-02_2000_hf-midpoint-coarse` 补入 O5 的脚本、结果与交接文档。

## 已完成：中点精阶段改用 FIRE（HF-EF-bracket-v5，TODO 第 9 条 9a）与同机验收（2026-11-04 14:34→17:48；用户决定见计划）

[施工计划](PLAN_HF_MIDPOINT_FINE_FIRE_2026-11-04.md)（含用户 11-04 原文）。
- 代码 `explore/src/explore/ts_ef.py`：`EFSettings.midpoint_fine_optimizer`（`bfgs` 默认 / `fire`，`fire` 须启用两阶段规则）；`_relax(..., optimizer=)`：FIRE 每次新建、速度从零、不用带入曲率，收敛后清掉 state 的 BFGS 曲率、保留分支几何，失败不改 state；每帧记 `optimizer`。只有中点精阶段用它；粗阶段、括区端点、延续仍为 BFGS。provenance `HF-EF-bracket-v5`。
- 测试 `test_ts_ef.py` 新增 5 项（设置校验 2 组；FIRE 只在精阶段、候选与 v3 相同；接受后保留几何并清曲率；失败 state 不变）。变异检验：粗阶段也换 FIRE、接受后保留旧曲率，各使 1 个测试失败；"去掉显式清除"的变异与原代码等价（FIRE 无 H，else 分支同样置空），未被测出，属预期。`test_ts_ef test_ts_branch test_ts_connect test_ts_warm_start test_d2_hf test_d2 test_d2_resume` **108 passed**；全包 `make check` 未跑。
- 验收脚本 `benchmarks/explore_d2/midpoint_fire_validate.py`（v4 / v5 同机配对完整 HS；候选有无或差异超过 `basin_energy_eV` / `basin_rmsd_A` 的项两版各重复 REPEATS 次）、`midpoint_fire_validate_run.sh`；`midpoint_coarse_validate.py` 加 `settings_for`（v3 / v4 / v5），midpoint 记录加 `optimizer` 字段，原 v3 / v4 行为不变。
- 冒烟：183，14:26 启动，PID 6415（启动前负载 0.00、GPU 空闲），`results/2026-11-04/midpoint_fire_validate/smoke/`。
- 冒烟结果（14:26→14:34，失败 1 项 + 对照 1 项，REPEATS 1）：失败项 `k0·10-12`：v4 两次都 `transverse_not_converged`（467 / 496 forward）；v5 两次都 `converged_candidate`（293 / 302 forward），精阶段帧均为 FIRE，粗阶段与其他为 BFGS。对照项：v4 / v5 都出候选（255 / 237）。候选有无不同的那一项触发了复核，复核路径跑通。
- 正式运行：183，14:34 启动，PID 6622（`results/2026-11-04/midpoint_fire_validate/pids.txt`），8 线程，`midpoint_fire_validate_run.sh full 0 2`（复核每项两版再各 2 次）。启动前负载 3.2、GPU 空闲。结果在 183 上读。
- **验收结果**（17:48 结束，126 次运行 = 首轮 74 + 复核 13 项 × 4，日志无报错；18:15 在 183 上核对进程已退出；在 183 上读）。汇总脚本 `benchmarks/explore_d2/midpoint_fire_validate_summary.py` → `results/2026-11-04/midpoint_fire_validate/summary.json`。
  - 失败 22 项（首轮 v4 → v5）：卡住 → 候选 13（其中 12 项原为 `transverse_not_converged`，1 项原为延续阶段 `walk_relax_unconverged`）；两版都出候选 1（`k3·10-12`；11-02 验收中两版都没出，本次 v4 也出了，属 fp32 波动）；两版都没出 8。
  - 复核（13 项，每项两版共 3 次）：v5 有 12 项 3 / 3 次出候选，`k0·6-8` 为 2 / 3（失败那次在延续阶段 `walk_relax_unconverged`，还没进中点）；v4 只有 `k1·6-8`、`k3·8-14` 各 1 / 3 次出候选，其余 0 / 3。救回主要是算法影响，不是 fp32 波动。
  - v5 在失败组出的 14 个候选，能量 −18880.406 至 −18880.414 eV，与新 E2 鞍点（−18880.41）同段，比已知两类低约 0.4 eV。**几何身份未比**；很可能多个意图找到同一个 TS，新增候选（14）与新增事件要分别计数。
  - 对照 15 项：两版全部出候选，能量差 ≤ 0.002 eV（1 个 fp32 量子），Kabsch RMSD ≤ 0.0046 Å，无一项触发复核。
  - 剩余 8 项失败：中点精阶段 FIRE 跑满 200 步 4 项（`k0·8-14` |RF| 0.0105、`k1·8-14` −0.030、`k2·6-8` −0.018、`k2·8-14` −0.0017；末帧横向力 0.0014–0.0027；`k0·8-14`、`k2·6-8` 即 O5 中 FIRE 未救回的两帧）；延续阶段 `walk_relax_unconverged` 3 项（还没进中点，v5 改动不涉及；v4 / v5 结局不同来自 fp32 波动）；两版都是 `no_resolved_maximum_bracket` 1 项。
  - 成本（首轮）：失败组 forward 15 736 → 14 452（−8%），墙钟 1 918 → 2 609 s（+36%，即每次 forward 0.12 → 0.18 s）；对照组 4 495 → 4 555 forward，416 → 423 s。两版在同一台机器上先后跑，183 上别人的负载在变（负载 3–5），墙钟差不能归因于 FIRE。
  - 按用户 11-04 的分流：剩余失败没有转到粗阶段，与中点相关的 4 项仍在精阶段，下一步是集中诊断剩余帧；14 个候选先比几何身份，再做 HC。


## 已完成：救回候选的几何身份（TODO 第 9 条 9b-1；2026-11-06；只读，不调模型；用户「继续做一做」）

用户 11-04 原文："已有两个救回候选也应先比较几何身份，再做 HC——多个意图可能找到同一个 TS，新增候选和新增事件要分别计数。"
- 脚本 `benchmarks/explore_d2/candidate_identity.py`，结果 `results/2026-11-06/candidate_identity/identity.json`（本机 181 读盘；源文件 11-04 写入，已过 NFS 6 h 缓存期，并按文件名去重）。判据用运行自己的 HF 判重规则：位姿图自同构群（|Aut| 1296）下的对称 RMSD ≤ `basin_rmsd_A` 0.15 Å 且能量差 ≤ `basin_energy_eV` 0.05 eV。
- 参照：八块 b05 HF 运行中全部 22 个 HC `connected` 事件，按同一规则分成 3 类：
  - 类 0：已知 E2，断 0–1、6–9，成 9–14，−18880.0195 eV，11 个事件；
  - 类 1：已知 E2，断 0–1、2–3，成 3–14，−18880.0195 eV，10 个事件；
  - 类 2：新 E2 鞍点，断 0–1、10–12，成 12–14，−18880.4102 eV，1 个事件（v4 sA `k1·10-12`，QM 已确认）。
  - 类间对称 RMSD：0.574（0–1）、1.406（0–2）、1.400（1–2）Å。
- 候选：11-02 v4 验收与 11-04 v4 / v5 验收中的全部 HS 候选，共 104 个（含重复运行）。**全部落在这 3 类里，没有新类。**
  - 对照组候选（60 个）分属类 0、类 1，对称 RMSD ≤ 0.007 Å。
  - 失败组候选（44 个，来自 15 个意图：v5 39 个、11-04 v4 3 个、11-02 v4 2 个）**全部是类 2，即新 E2 鞍点**：对称 RMSD 中位 0.014、最大 0.038 Å，能量差 −0.012 至 +0.004 eV。11-02 v4 救回的 `k1·8-14`、`k3·4-14` 也在其中。
- 计数：**新增候选 15 个意图，新增事件 0 个**——都是此前已由 HC 连通、QM 确认的同一个新 E2 TS。变化在于：这个 TS 此前在四次整跑里只被找到 1 次（二分中点碰巧收敛），现在从 15 个意图都能稳定找到。
- 限定：身份基于几何与能量，这些候选本身没有跑 HC；类 2 只有 1 个参照事件；类 2 内部的 RMSD 分散（≤ 0.038 Å）比类 0 / 1（≤ 0.007 Å）大，与其软横向模一致（推断）。

## 已完成：v5 救回候选的 HC（TODO 第 9 条 9b-2；2026-11-06 14:05→14:36；用户「可以」）

- 脚本 `benchmarks/explore_d2/midpoint_fire_hc.py`（22 个失败意图各用 v5 重跑一次完整 HS——验收记录没存 HC 需要的分支帧；出候选就照 `d2._hc_task` 做 HC：`branch_starts` + `trace_sides_ef`，各自账本 hs_cap / hc_cap；记录 HC 状态、两端、化学对，以及化学对是否与参照新 E2 事件（v4 sA `distance:10-12`）相同），启动脚本 `midpoint_fire_hc_run.sh OUT LIMIT`。
- 冒烟（183，13:55 时负载 0、GPU 空闲；1 项，`results/2026-11-06/midpoint_fire_hc/smoke/`）：`k0·10-12` HS 出候选（291 forward，−18880.4121 eV），HC `connected`（349 forward），化学对与参照相同。
- 正式运行：183，14:05 启动（PID 见 `results/2026-11-06/midpoint_fire_hc/pids.txt`），8 线程，输出 `full/`，在 183 上读。
- **结果**（14:36 结束；15:57 在 183 上核对：进程已退出，22 条结果齐全，日志无报错）。汇总脚本 `benchmarks/explore_d2/midpoint_fire_hc_summary.py` → `results/2026-11-06/midpoint_fire_hc/summary.json`。
  - HS（v5 重跑）：15 个意图出候选（能量 −18880.406 至 −18880.414 eV）；4 个停在延续阶段 `walk_relax_unconverged`；3 个停在中点精阶段 `transverse_not_converged`。
  - HC：13 个 `connected`，两端化学对**全部**与参照新 E2 事件（v4 sA `distance:10-12`）相同；2 个 `unresolved:endpoint_unstable`（`k2·2-4`、`k2·6-8`）：低能侧都下降到稳定端点（约 −18880.74 eV），高能侧终点能量仍为 −18880.410 eV，即没有离开鞍点，扰动复核判为不稳定；原因未查。
  - 与 11-04 验收的差别：`k2·8-14` 这次出候选（11-04 v5 三次都卡在精阶段）；`k1·8-14`、`k0·6-8` 这次没出（11-04 分别出过 / 2 / 3）；另有几项的失败位置在延续阶段与中点精阶段之间变化。延续阶段代码 v4 / v5 相同，这些差别归为 fp32 波动（推断）。
  - 计数：新增候选 15 个意图，新增 HC 连通 13 个，**新增事件 0**（都连通到同一个新 E2 TS）。
  - 成本：HS 14 373 forward，HC 4 393 forward（每个连通 216–555）。
  - 未做：这次重跑得到的候选没有再用对称 RMSD 比几何（`candidate_identity.py` 只认验收结果的文件名格式）；身份依据是 HC 化学对与能量。

- 交接（11-06，用户「你写一个 handoff 我要开新的窗」）：[HANDOFF_HF_MIDPOINT_2026-11-06.md](HANDOFF_HF_MIDPOINT_2026-11-06.md)，接替 11-04 版；提交草稿 `2026-11-06_1752_hf-midpoint-identity-hc` 补入该文件与 README。

## 已完成：v5 剩余精阶段卡点的诊断（TODO 第 9 条 9c；2026-11-06；用户「顺着handoff往下做」，全量用户「跑」）

用户 11-04 原文："若仍卡在精阶段，就集中诊断剩余帧"。用户 11-06 另贴了 6 条性能建议要求审视，其中第 2 条（FIRE 速度投影到约束切空间、记录被约束修正掉的位移）的诊断量并入本步测量。
- 离线统计（只读，11-04 v5 验收结果，本机读，源文件已过 6 h 缓存期）：
  - HS 费用构成（出候选的运行）：失败组 39 次共 19 117 forward，其中二分（refine）5 688（30%），延续 13 429（70%）；对照组 15 次 4 555 forward，二分 876（19%）。二分内部，失败组精阶段 3 732、粗阶段 1 853 forward；二分轮数中位：失败组 4、对照组 1。
  - **粗阶段已有候选、仍被要求精松弛**：某轮粗阶段两侧都收敛、未提前收括区、且至少一侧已满足候选判据（完整力 ≤ 0.01、横向力 ≤ 0.003；粗阶段横向力 ≤ 2e-3，故等价于完整力 ≤ 0.01），并且两侧反力差 ≤ 门 0.0106。这种轮次：v5 失败组 40 次（精阶段随后收敛 38、失败 2：`k0·8-14`、`k2·8-14`）、对照组 16 次；其后的精阶段共花 2 764 / 447 forward（占两组 HS 总费用的 14% / 10%）。v4 臂中这种轮次的精阶段失败 13 次，即 v4 的剩余失败多数属此类。原因是 v4 规则里"一侧已满足候选判据就转精阶段"（COARSE 计划 §1 第 2 步，我加的，列在该计划"未确认"）。几何一致性在离线数据里没有存，未核对。
  - 11-04 v5 剩余 4 个精阶段失败的末帧：`k0·8-14` RF 0.0105、横向 0.0027；`k1·8-14` −0.030、0.0014；`k2·6-8` −0.018、0.0017；`k2·8-14` −0.0017、0.0026（另一侧精阶段已收敛到 0.0008）。
- 新脚本：`benchmarks/explore_d2/midpoint_fine_diag.py`（按 11-04 `items.json` 取项，只读所需源 run；v5 重跑完整 HS，逐帧记录全部中点松弛；精阶段 FIRE 换成数值不变、只记录内部状态的子类：dt、α、速度重置、|v|、速度沿当前约束法向的比例、提议步长、是否截断、提议位移被约束修正掉的比例；精阶段未收敛的段在最小横向力帧做同几何 16 次重复求力、解析 Hessian 横向谱及力的模分解、16 个 0.002 Å 横向随机位移并扣除 Hessian 线性项的"粗糙度"，在末帧做 16 次重复求力），`midpoint_fine_diag_run.sh OUT REPS ITEMS`，`midpoint_fine_diag_summary.py`（只读，按文件名去重）。Hessian 只作离线诊断。
- 冒烟（183，18:2x，PID 508857，启动前负载 0.00、GPU 空闲；`results/2026-11-06/midpoint_fine_diag/smoke/`，汇总 `smoke_summary.json`，在 183 上读）：项 3（`k0·8-14`，v4 sA）1 次。HS 775 forward、76 s（0.10 s/forward），诊断 12 s。结局重现：第一个中点（q 1.3344）粗阶段两侧 RF 0.0082 / 0.0117，左侧完整力 0.0098（**已满足候选判据**），两侧对齐 RMSD 0.011 Å、反力差 0.0035，按现行规则转精阶段；左侧 FIRE 200 步未收敛（末帧横向 0.0014）。
  - 第 2 条建议的诊断量：速度沿约束法向的比例中位 2.7e-6、p90 1.3e-5；提议位移被约束修正掉的比例中位 1e-11、p90 2e-5；无截断。**这一帧上速度带入法向分量可以忽略。**
  - FIRE 状态：200 步内速度重置 18 次（前后各 9），dt 在 0.02–0.17 间反复，实际每步最大原子位移中位 4e-5 Å；末 100 步 log(横向力) 斜率 −3e-4 / 步，即停滞，不是慢收敛。
  - 最小横向力帧（第 199 步，记录 0.00115）：同几何 16 次重复的横向力（逐原子最大）中位 0.00155、最大 0.0026；只取重复间的差（纯噪声）中位 0.00087、最大 0.0017；逐分量标准差最大 0.0010、中位 0.00014（噪声集中在少数原子）；反力标准差 0.0008。16 次平均力的横向力 0.0014，其中最软横向模（0.0199 eV/Å²，次低 0.145）上的分量逐原子最大 0.0009（对应沿该模位置偏 0.053 Å），其余模 0.0011；噪声在软模上的标准差 0.00022（折合位置 0.011 Å）。0.002 Å 位移、扣除线性项后的横向力中位 0.00145、最大 0.0030。末帧重复求力：横向力中位 0.0016，纯噪声中位 0.00096。
  - 读法（只此 1 帧）：这一帧单次求值的纯噪声对逐原子最大横向力的贡献（中位 0.0009）已接近 1e-3 阈值，FIRE 在此反复重置、步长塌到 1e-4 Å 以下。与 11-01 BFGS 失败段（末段横向力约 10 倍阈值、大步游走）不同。是否所有剩余帧都如此，要看全量诊断。
- 全量（用户「跑」）：183，22:55 启动，PID 509622（`results/2026-11-06/midpoint_fine_diag/pids.txt`），启动前负载 0.00、GPU 空闲；`midpoint_fine_diag_run.sh full 2 3,6,8,11,12`，23:14 结束，10 个结果齐全，日志无 error / traceback。汇总 `results/2026-11-06/midpoint_fine_diag/full_summary.json`（在 183 上读）。HS 合计 8 247 forward；另有诊断求值（每个卡住段 48 次求力 + 1 次解析 Hessian，5 段共 240 次，不入账本）。
  - 结局（5 项 × 2 次）：出候选 3（`k2·6-8` r1、`k2·8-14` r0 / r1）；精阶段 FIRE 200 步未收敛 5（`k0·8-14` ×2、`k1·2-4` r0、`k1·8-14` r0、`k2·6-8` r0）；粗阶段 BFGS 200 步未收敛 1（`k1·8-14` r1，右侧末帧横向 0.0041）；延续阶段 `walk_relax_unconverged` 1（`k1·2-4` r1）。同项两次结局不同的有 3 项。
  - **第 2 条建议（速度法向分量）**：测到的 9 段 FIRE（卡住 5、收敛 4）速度沿约束法向比例 p90 ≤ 2.6e-5，提议位移被约束修正掉的比例 p90 ≤ 5e-5，无截断。可以排除。
  - **卡住的 5 段**：200 步内速度重置 16–18 次，实际每步最大原子位移中位 5e-6–3e-5 Å，末 100 步 log(横向力) 斜率 −0.002 至 +0.002 / 步（不下降）。最小横向力帧（记录 0.00101–0.00106）上同几何 16 次重复：单次横向力中位 0.0013–0.0019（**5 段全部 > 1e-3**），16 次平均力的横向力 0.00115–0.00168（也全部 > 1e-3），纯噪声部分中位 0.0007–0.0015；反力标准差 0.0009–0.0022。记录的 ≈ 0.00103 是 200 次带噪求值里取最小的选择效应。0.002 Å 横向位移扣除线性项后横向力中位 0.0011–0.0022，与噪声同量级。
  - 最软横向模 0.019–0.021 eV/Å²（次低 0.12–0.15）；平均横向力中软模所占平方比例 0.04–0.35，其余模上逐原子最大 0.0009–0.0017。残余力不主要在软模上。
  - **收敛的 FIRE 段也一样**：末 50 帧横向力中位 0.0014–0.0023，重置 2–16 次；收敛发生在某一帧带噪读数碰巧 ≤ 1e-3 时。
  - 其他停止：`k1·2-4` r0 第二个中点（q 1.3251）两侧粗阶段 RF −0.075 / −0.059 都超过 B，但反力差 0.016 超门 0.0106（两侧对齐 RMSD 0.014 Å）而转精阶段；精阶段左侧收敛（−0.0755）、右侧卡住（−0.0588），即使收敛也会因反力差超门停止（9e 的门）。`k1·8-14` / `k2·6-8` r0 卡在 |RF| 0.019–0.033 的近根中点，粗阶段两侧都不满足候选判据。
  - 粗阶段已有候选：`k0·8-14` r1（两侧完整力 0.0099 / 0.0095，RMSD 0.0008 Å）属此类；`k2·8-14` 两次也属此类，这次精阶段 19–24 步就收敛了。
  - 读法：在这些帧上，1e-3 低于实测的单次求值横向力水平，FIRE 已经停止移动，不是慢收敛；精阶段能否收敛取决于噪声读数是否碰巧低于阈值，这与同项重复运行结局不同一致。所以预条件（第 3 条）对这些帧帮不上。剩余停止分三类：①粗阶段已有候选（`k0·8-14`）；②近根、无候选、精阶段阈值在噪声内（`k1·8-14`、`k2·6-8`、`k1·2-4` 一次）；③反力差门（`k1·2-4` 一次）与粗阶段 BFGS 未收敛（`k1·8-14` 一次）。改法（判据变更）需用户决定，未动代码。
- 改法选项文档（用户 11-06「可以」）：[HF_MIDPOINT_FINE_STAGE_OPTIONS_2026-11-06.md](HF_MIDPOINT_FINE_STAGE_OPTIONS_2026-11-06.md)。补测（只读，在 183 上读）：FIRE 段末 50 帧首末对齐 RMSD 卡住段 1.4e-4–5.7e-4 Å、收敛段 1.9e-5–8.4e-4 Å，窗口内反力标准差 0.0013–0.0026；卡住段窗口平均反力与同中点另一侧 / 粗阶段读数差 0.001–0.0025。

## 已完成：粗阶段出候选（P1，HF-EF-bracket-v6，TODO 第 9 条 9c）施工、冒烟与同机验收（2026-11-06→07；用户决定见计划）

[施工计划](PLAN_HF_MIDPOINT_COARSE_CANDIDATE_2026-11-06.md)。用户 11-06：问"P1 可以定下来吗？"，回答「可以」。
- 代码 `explore/src/explore/ts_ef.py`：`EFSettings.midpoint_coarse_candidate`（bool，默认 False；只允许在 v5 上开：两阶段规则 + 精阶段 FIRE）。`coarse_sign` 在括区宽度、几何、双侧反力差门都通过后，若任一侧满足 `stationary`，开关开时给两帧记 `coarse_candidate = True`（`coarse_sign` 仍为 0）；二分循环随即取满足一侧中完整力最小者出候选，HS 结束。开关关时 v3 / v4 / v5 不变。开启时 provenance `HF-EF-bracket-v6`。D2 / `HFSettings` 未接。
- 测试 `test_ts_ef.py` 新增 6 项：设置校验 3 组；粗阶段已有候选且门通过时 v6 在粗阶段出候选、没有精阶段帧、费用低于 v5，同输入 v5 仍转精阶段；一侧满足候选判据但反力差超门时不在粗阶段出候选。`test_ts_ef test_ts_branch test_ts_connect test_ts_warm_start test_d2_hf test_d2 test_d2_resume` **113 passed**。变异检验（临时目录 + `-o pythonpath=`）：候选检查挪到反力差门之前 → 1 项失败；忽略开关 → 4 项失败。全包 `make check` 未跑。
- 验收脚本：`benchmarks/explore_d2/midpoint_coarse_candidate_validate.py`（v5 / v6 同机配对，项与复核规则同 11-04）、`_run.sh`、`_summary.py`（含对称 RMSD、v6 候选出自粗 / 精阶段、v5 中"粗阶段已有候选之后"的精阶段 forward）。`midpoint_coarse_validate.settings_for` 加 v6、midpoint 记录加 `coarse_candidate`；`candidate_identity.py` 文件名规则放宽到 v3–v6。原 v3–v5 行为不变。
- 冒烟（183，00:03 启动，PID 510606，启动前负载 0.00、GPU 空闲；`results/2026-11-07/midpoint_coarse_candidate_validate/smoke/`，汇总 `smoke_summary.json`，在 183 上读；日志无报错）：失败项 `k0·10-12`：v5 288 / v6 293 forward，v6 末轮在粗阶段出候选，与 v5 候选 ΔE 0.002 eV（1 个 fp32 量子）、对称 RMSD 0.0018 Å。对照项 `k2·10-13`：v5 230 / v6 207 forward，v6 在粗阶段出候选，ΔE 0、RMSD 0.0008 Å。0.09–0.10 s/forward；读入八块 journal 约 7 min。
- 正式验收改交 PBS（用户 11-07："yayoi40 … 跑这个上面"）：脚本 `benchmarks/explore_d2/midpoint_coarse_candidate_validate_pbs.zsh`（照模板只改 ← 行：yayoi40，RTX 5080，ncpus 32 / ngpus 1 / 56 GB / 24 h；参数与 `_run.sh full 0 2` 相同），输出 `results/2026-11-07/midpoint_coarse_candidate_validate/full/`。`qsub` 由用户执行。
- **验收结果**（PBS 作业 24581，yayoi40 RTX 5080；日志 `b05_v6val.o24581`（项目根目录）有 done 行、无 error / traceback；节点时钟记 18:07:38 → 18:39:15，与本机时钟不一致，墙钟按差值约 32 min；94 个结果齐全；用户「跑完了」后在 183 上读）。汇总 `results/2026-11-07/midpoint_coarse_candidate_validate/summary.json`，身份 `identity.json`。0.036 s/forward。
  - 失败 22 项首轮：v5 15 项出候选，v6 16 项。v5 卡在精阶段 FIRE → v6 粗阶段出候选 3 项（`k0·2-4`、`k0·8-14`、`k2·6-8`）；v5 出候选 → v6 `no_resolved_maximum_bracket` 2 项（`k5·11-14`、`k7·6-8`），这两次停在第一个中点之前，v5 / v6 在那里代码相同，归为 fp32 波动（推断）。
  - 复核 5 项（每项每版共 3 次，v5 → v6 出候选次数）：`k0·2-4` 1 → 3，`k0·8-14` 0 → 2，`k2·6-8` 0 → 3，`k5·11-14` 1 → 2，`k7·6-8` 1 → 1。
  - 对照 15 项：两版全部出候选，v6 全部在粗阶段出；两版候选 ΔE ≤ 0.0039 eV（2 个 fp32 量子）、对称 RMSD ≤ 0.0049 Å。失败组两版都出候选的 13 项：ΔE ≤ 0.002 eV、对称 RMSD ≤ 0.0040 Å。
  - v6 候选来源：失败组 24 个中 23 个在粗阶段、1 个在精阶段；对照组 15 / 15 在粗阶段。
  - 几何身份（`candidate_identity.py`，运行判重规则）：70 个候选（含重复）全部落在八块 b05 HF 的 3 个已连通类里；失败组 v5 16 个、v6 24 个全部是类 2（新 E2 鞍点，对称 RMSD ≤ 0.0375 / 0.0368 Å），对照组分属类 0、1。新增事件 0。
  - 剩余失败停止阶段（失败组全部运行，各 32 次）：v5——精阶段 FIRE 8、粗阶段 BFGS 1、反力差门 1、第一个中点之前 6（`walk_relax_unconverged` 6）；v6——精阶段 FIRE 1（`k0·8-14` 一次，RF 0.009、横向 0.0022）、粗阶段 BFGS 200 步未收敛 2（`k1·2-4` 横向 0.0021、`k1·8-14` 0.0155）、第一个中点之前 5（`no_resolved_maximum_bracket` 3、`walk_relax_unconverged` 1、`hysteresis` 1）。
  - 成本（首轮）：失败组 13 699 → 11 240 forward（−18%），墙钟 493 → 402 s；对照组 4 459 → 4 101（−8%），159 → 151 s。v5 首轮在"粗阶段已有候选"之后的精阶段花了 1 551 / 398 forward。两组结局不完全相同，费用差不全是 P1 的作用。
  - 小结：v6 不改变对照组候选，精阶段失败由 8 次降到 1 次，救回的候选都是已知的新 E2 鞍点，forward 减少。剩余失败多在第一个中点之前（5 / 8），中点相关的 3 次中 2 次在粗阶段。

## 已完成：第一个中点之前的 HS 停止普查（TODO 第 9 条 9c'' / 第 10 条；2026-11-07；只读，不调模型；用户「可以」）

- 脚本 `benchmarks/explore_d2/hs_prewalk_census.py`（只读 journal），输入 10-29 yayoi32 v4 两块（现行延续代码 HF-EF-branches-v4），结果 `results/2026-11-07/hs_prewalk_census/census.json`（本机读，源文件 10-29 写入）。
- HS 结局（sA / sB，项数与 forward）：`hysteresis` 174 / 179 项，177.8k / 179.1k；`walk_relax_unconverged` 48 / 43 项，47.2k / 40.4k；`transverse_not_converged` 13 / 17 项（其中端点失败 4 / 9，其余是 v4 以前的中点失败），6.7k / 6.6k；`above_ceiling` 24 / 22；`no_resolved_maximum_bracket` 0 / 3；`budget_exhausted` 0 / 2；出候选 8 / 6 项，2.6k / 1.7k。两块 HS 共约 470k forward，`hysteresis` 与 `walk_relax_unconverged` 合计约 95%。收敛的延续试步共 19k 次，步数中位 4、p90 86–88、p99 163–169。
- `hysteresis`（353 项）：两侧结束方式为一侧跳变一侧走到头 226、两侧跳变 64、两侧走到头 63。反力沿延续点的符号：低 q 一侧 328 / 353 全为正、高 q 一侧 318 / 353 全为负，两侧都在各自分支上没有翻转；只有 2 项某一侧同时出现正负，其中 1 项是经过 |RF| ≤ 0.003 的点（"相邻两点检查漏掉翻转"只此 1 例）。每侧最小 |RF|：中位 0.35、p10 0.096，只有 5% 低于 0.05。每侧 20 点左右，单项中位约 1 000 forward。跳变：对齐 RMSD 中位 0.52 Å，能量变化中位 −0.09 eV。
- `walk_relax_unconverged`（91 项，95 次失败试步）：失败试步全部跑满 200 步；从试步起点算，对齐 RMSD 中位 1.0 Å（p10 0.24、最大 2.0），单原子最大位移中位 2.7 Å；能量变化中位 −0.33 eV（p10 −2.36）；57 次末 50 步能量仍下降 > 0.01 eV；44 次反力在试步内变号；末 50 步每步位移中位 0.022 Å（上限 0.05）；过程最低横向力中位 0.0038。终止一侧的反力全程同号；另一侧走到头 66、跳变 21、也是试步失败 8，**没有一侧出现反力翻转（0 / 95）**。每侧最小 |RF| 中位 0.42。读法：结构在固定 q 上滑向另一个盆地，与 `hysteresis` 同属"分支在此终止、两侧都没有最大值"，不是数值收敛问题（推断）。
- `no_resolved_maximum_bracket`（sB 3 项：`k7·6-8`、`k1·6-8`、`k5·11-14`）：延续括区的反力为 +0.53 / −0.005、+0.50 / −0.0046、−0.0058 / +0.074，有一端只比 tol 0.003 略大；refine 重新松弛这一端得到 −0.0022 / +0.0011 / −0.0020，落在 ±tol 内，于是停止。**这三个端点帧都满足候选判据**：完整力 0.0028–0.0030、横向力 0.0024–0.0029，能量 −18880.41 至 −18880.42 eV（新 E2 鞍点所在能量段，几何身份未比）。11-07 v6 验收中的 3 次 `no_resolved_maximum_bracket`（`k0·6-8`、`k5·11-14`、`k7·6-8`）记录里没存端点帧，是否同一情况未核对。
- 端点 `transverse_not_converged`（13 项）：过程最低横向力 0.0031–0.0043（阈值 0.003），对齐 RMSD 中位 0.38 Å，能量变化中位 −0.008 eV，每步位移中位 0.008 Å；与 11-01 BFGS 中点游走相似（推断）。

## 已完成：refine 端点出候选（R1，HF-EF-bracket-v7，TODO 第 9 条 9c''）施工与定向冒烟（2026-11-07；用户决定见计划）

[施工计划](PLAN_HF_ENDPOINT_CANDIDATE_2026-11-07.md)。用户 11-07：问"R1 可以定吗？"，回答「可以」。
- 代码 `explore/src/explore/ts_ef.py`：`EFSettings.endpoint_candidate`（bool，默认 False，要求 v6）。两端重松弛后括区未分辨时：右端满足 `stationary` 且左端反力 > tol，或左端满足且右端 < −tol，就由该端出候选（都成立取完整力小者），记 `endpoint_candidate = True`；否则照旧 `no_resolved_maximum_bracket`。已分辨的括区照旧二分。provenance `HF-EF-bracket-v7`。"另一端须指向最大值"是我加的限制（计划"未确认"）。D2 未接。
- 测试 `test_ts_ef.py` 新增 5 项（设置校验 3 组；端点恰在 TS 时 v6 停止、v7 由右端出候选且无中点帧；q 方向极小处的端点不出候选）。相关测试 **118 passed**。变异检验：去掉"另一端须指向最大值"→ 1 项失败；忽略开关 → 1 项失败。全包 `make check` 未跑。
- `midpoint_coarse_validate.settings_for` 加 v7；`candidate_identity.py` 文件名规则放宽到 v3–v7。
- 定向冒烟（183，09:0x，启动前负载 0.00、GPU 空闲；脚本 `benchmarks/explore_d2/endpoint_candidate_smoke.py`，输出 `results/2026-11-07/endpoint_candidate_smoke/`，在 183 上读）：第一次启动（PID 513840）在第一次求力前因坐标是 list 而崩溃，修正后重启（PID 513949），日志无报错。10-29 v4 sB 的 3 个 `no_resolved_maximum_bracket` 项，取延续括区冷启动 refine：
  - `k1·6-8`：v6 右端 RF −0.0017 → `no_resolved_maximum_bracket`（3 forward）；v7 右端 RF −0.0001 → **右端出候选**（3 forward），对称 RMSD 0.024 Å 匹配类 2（新 E2 鞍点）。
  - `k7·6-8`：冷启动后右端 RF −0.0046 / −0.0053，括区已分辨，R1 不起作用；v6 粗阶段出候选（类 2），v7 在二分中 `branch_force_disagreement`（括区已分辨后 v6 / v7 代码相同，归为 fp32 波动，推断）。
  - `k5·11-14`：括区已分辨，两版都在粗阶段出候选（类 2）。
  - 即：原记录中的 3 项冷启动重跑后只有 1 项仍是未分辨括区，v7 在这一项上由端点出候选。

## 已完成：HS 两侧无最大值（R3）的测量与改法选项（2026-11-07；只读，不调模型；用户「可」）

- 选项文档：[HF_HS_NO_MAXIMUM_OPTIONS_2026-11-07.md](HF_HS_NO_MAXIMUM_OPTIONS_2026-11-07.md)。新脚本 `benchmarks/explore_d2/hs_walk_profile.py`（每侧延续轨迹、跳变与失败试步前后的成键变化，共价半径 × 1.2）、`hs_ridge_continuity.py`（ridge 类 HS 的 L1 扫描峰前后几何与扫描自带 `edges`），结果在 `results/2026-11-07/hs_walk_profile/`、`hs_ridge_continuity/`（本机读，源文件均早于 6 h）。八块 b05 HF 共 2 056 个 HS，其中 ridge 类 1 702 个。
- v4 两块：ridge 类 HS 441k forward，其中 `hysteresis` 301 项 337k、`walk_relax_unconverged` 75 项 80k。延续中的 384 次跳变全部无成键变化；失败试步有成键变化的 ridge 类 15 / 79、`branch_bracketed` 类 15 / 15。
- L1 扫描峰前后成键图不变的比例：找到连续最大值的 139 项中 0.04（5 项，都以 `branch_unresolved` 结束、无候选），`hysteresis` 1 255 项中 0.98，`walk_relax_unconverged` 106 项中 0.75，`transverse_not_converged` 194 项中 0.75。52 个 HS 候选全部来自成键改变的脊；22 个新 E2 意图在各块中峰前后成键都改变。成键不变的脊占 ridge 类 HS forward 的 92%（八块 1.64M / 1.78M；v4 两块 411k / 441k，362 / 420 项）。峰前后相邻扫描点的对齐位移分不开两组（最大值组 0.82–1.13 Å，是其他步的 15–28 倍；`hysteresis` 组 0.57–1.39 Å）。

## 已完成：同键图的脊最后派发（R3a-2，D2-HF-v1.2）施工与离线重放（2026-11-07；用户决定见计划）

[施工计划](PLAN_HF_SAME_GRAPH_RIDGE_2026-11-07.md)。用户 11-07：问"R3a 用 R3a-2（降优先级、不拒跑），可以吗？"，回答「可以」。
- 代码 `explore/src/explore/d2.py`：`init_run(..., same_graph_ridge_last=True)`（需要 `hf` 与 `start_cluster`）冻结为 `D2-HF-v1.2`，manifest 加 `same_graph_ridge_last`（旧规则不写此键）。`_same_graph_ridge(classified, points)`：`ridge_suspected` 且扫描峰点与下一点 `edges` 相同。`_promote`：同键图的脊在起点分簇中排最后参与；派发在第一、二轮之后，额度不够记 `HS_deferred:same_graph_ridge`；HS payload 加 `same_graph_ridge` 标签；`_hf_report` 加 `same_graph_ridge: {settled, deferred}`。`tests/d2_bench.py` 的 `init` 透传新参数。
- 测试 `test_d2_hf.py` 新增 4 项（判定函数、冻结、额度充足时最后派发且不当簇首、额度不足时暂缓）。变异检验：不重排起点分簇 / 不放最后 / 暂缓原因写错，分别使 2 / 2 / 1 项失败。`test_ts_ef test_ts_branch test_ts_connect test_ts_warm_start test_d2_hf test_d2 test_d2_resume test_probe` **147 passed**；全包 `make check` 未跑。
- 离线重放 `benchmarks/explore_d2/same_graph_ridge_replay.py`（八块 b05 HF，用新代码的 `_classify` 与 `_same_graph_ridge`；结果 `results/2026-11-07/same_graph_ridge_replay/replay.json`，本机读）：2 634 个升级项中 1 663 个被排到最后（ridge 类 1 932 个），原 HS 费用占 87.3%（1.64M / 1.88M）；其中候选 0 个（52 个候选都不受影响）；22 个新 E2 意图与 15 个对照意图不受影响（`k5·3-14` 在 v4 sA 命中的是同一缩写下另一个起始姿态，以 `hysteresis` 结束）。HC 未读入（`load` 只保留 L1 / L2 / HS）。v4 两块各有 32–34 个非同键图的项原来因额度用完没跑到 HS（`branch_bracketed` 18–21、`ridge_suspected` 3–4、`audit` 8–12），v1.2 下排在同键图的脊之前。b05 额度下同键图的脊仍会花掉剩余额度，v1.2 本身省不了多少 forward（推断）。整跑验收未做。

## 已完成：HS 括区规则 v7 接入 D2（D2-HF-v1.3，TODO 第 9 条 9d）施工、冒烟与整跑准备（2026-11-07；用户「准备」）

[施工计划](PLAN_HF_D2_V13_2026-11-07.md)。
- 代码 `explore/src/explore/d2.py`：`HSBracket`（粗阈值、符号界线，无默认值）；`init_run(..., hs_bracket=)` 需要 `same_graph_ridge_last`，冻结时用 `EFSettings` 校验，manifest 加 `hs_bracket`，规则版本 `D2-HF-v1.3`；`_ef_settings(hf, hs_bracket=None)`（不给时与原来相同）；`_hs_task` 用冻结的括区规则。`run_b05_d2.py` 加 `--same-graph-ridge`、`--hs-bracket`（`HS_BRACKET_B05` 2e-3 / 0.025）。`tests/d2_bench.py` 透传。
- 测试 `test_d2_hf.py` 新增 2 项（冻结与校验；解析台 v1.3 整跑结局与参照相同、HS provenance 为 v7、计数闭合）。变异检验：HS 不用冻结规则、漏掉端点出候选，各使 1 项失败。相关测试 **149 passed**；全包 `make check` 未跑。
- 冒烟（183，14:53 启动，PID 516074，启动前负载 5.1（别的用户）、GPU 空闲；与 10-26 b05 HF 冒烟同参数加 v1.3 开关；`results/2026-11-07/b05_hf_v13_smoke/`，日志 `b05_hf_v13_smoke.log` 无 traceback，在 183 上读）：D2-HF-v1.3；1 个升级项（同键图的脊），最后派发，HS `hysteresis` 977 forward；search 3 354 + refine 977 = ledger 4 331，unreconciled 空。HS 未进入二分，v7 括区规则在 D2 里未被真实模型走到。
- 整跑 PBS 脚本 `benchmarks/explore_d2/b05_hf_v13_pbs.zsh`（照 10-29 v4 脚本只改 ← 行；host 暂填 yayoi40；两块同种子同参数），未提交，等用户给节点。

## 已完成：b05 D2-HF-v1.3 整跑（PBS yayoi32 作业 24794，2026-11-07 15:54:47 → 21:10:01；用户「跑完了」后在 183 上读）

[计划](PLAN_HF_D2_V13_2026-11-07.md)。日志 `b05_hf_v13.o24794`（项目根目录），两块 `results/2026-11-07/b05_hf_v13_s{A,B}/report.json` 齐全，日志无 error / traceback，unreconciled 空。与 10-29 yayoi32 v4 整跑同机、同种子、同额度。
- 账目：sA search 339 647 / refine 237 076（ledger 576 723）；sB 339 623 / 238 220（ledger 577 843）。两块的 refine 额度都基本用完。
- HS 出候选：sA 8 → 21，sB 6 → 18。候选阶段：sA 粗阶段 19、精阶段 2，sB 粗阶段 18；端点出候选 0 次（v4 全部在精阶段）。按派发组：第一轮成键改变的脊 11 / 8，第二轮 10 / 10（v4 第二轮 3 / 2）；同键图的脊 0（sA 180 个 HS：`hysteresis` 147、`walk_relax_unconverged` 32、`transverse_not_converged` 1；sB 183 个：150 / 27 / 6）。
- HC：sA 连通 3、重复 13、`endpoint_unstable` 3、`tangent_no_resolved_drop` 2；sB 连通 3、重复 10、`endpoint_unstable` 1、`tangent_no_resolved_drop` 4（v4：sA 3 / 4 / 0 / 1，sB 2 / 2 / 1 / 1）。
- 不重复事件的身份（按八块参照类，对称 RMSD + 能量，`results/2026-11-07/b05_hf_v13_events_identity.json`）：sA 类 0、1、2（与 v4 相同）；sB 类 1、0、2，**新 E2 鞍点（类 2）第一次在 sB 种子出现**（v4 sB 只有类 0、1）。无新类。
- 派发：第二轮因额度用完暂缓（`HS_deferred:same_start`）60 / 58 → 0 / 0；同键图的脊暂缓（`HS_deferred:same_graph_ridge`）32 / 26。第二轮的 `branch_bracketed` 这次跑到（8.1k / 8.4k forward），以 `walk_relax_unconverged` 为主，无候选。
- forward：HS 235.1k → 230.6k（sA）、235.5k → 233.5k（sB），其中同键图的脊 193.4k / 188.9k；HC 2.7k → 6.5k、2.5k → 4.7k。墙钟 5.4 h → 5.3 h。
- 读法：v1.3 在同额度下候选与不重复事件增加，没有新类；同键图的脊仍花掉约 190k forward / 块，产出 0。候选增加中 v7 括区规则与派发顺序各占多少，这次分不开。
- 选项页（用户 11-07「可」）：[HF_SAME_GRAPH_RIDGE_SPEND_OPTIONS_2026-11-07.md](HF_SAME_GRAPH_RIDGE_SPEND_OPTIONS_2026-11-07.md)。只整理数字：非同键图项花费 43.7k / 49.3k（HS + HC）/ 块，同键图的脊 193.4k / 188.9k、0 候选；选项 S0–S3，未定。

## 计划：同键图的脊抽查（S3）+ HS / HC 并行（D2-HF-v1.4）（2026-11-07；用户「可以」；未施工）

[施工计划](PLAN_HF_D2_V14_2026-11-07.md)（含用户原话）。只写了计划：S3 按种子、以 `audit_fraction` 抽同键图的脊，其余记 `HS_skipped:same_graph_ridge`；HS / HC 照 `_run_scans` 按组分波次、`run_many` 共用批量求值器并行，按输入顺序结算，HC 判重在结算时按顺序做；先同构（HS 波 → HC 波），测量后再定是否流水线（异构）。先量后建：基线 GPU 利用率与每 forward 墙钟未测。

## 进行中：D2-HF-v1.4（S3 + HS / HC 并行）施工、测试与测量冒烟（2026-11-07；用户「可以」）

[施工计划](PLAN_HF_D2_V14_2026-11-07.md)。
- 代码 `explore/src/explore/d2.py`：`init_run(..., same_graph_sample=True, hf_parallel=True)`（两者同给，需要 v1.3）冻结为 `D2-HF-v1.4`，manifest 加两个键（旧规则不写）。S3：同键图的脊按 `_key(seed, "same_graph_sample", intent_id) < audit_fraction` 抽样，抽中的仍在最后一组，其余记 `HS_skipped:same_graph_ridge`（HC 不跑），HS payload 加 `same_graph_sample`，报告 `same_graph_ridge` 加 `skipped`、`sample_candidates`。并行：HS / HC 拆成执行（`_hs_run`、`_hc_run`，接受任意计数工厂）与结算（HC 的化学对与判重在 `_hc_finish`，结算时按输入顺序做）；串行路径（v1.1–v1.3）改为调用同一组函数，逻辑不变；`_hf_parallel` 按组（第一轮 → 第二轮 → 同键图的脊样本）分波次，每波预留上限、`run_many` 共用 `search_batch` 批量求值器（`max_active` / `wave_size` 沿用冻结值）、按输入顺序结算，HS 波后接该组的 HC 波；一项都派不出去的波才按现行原因跳过。`tests/d2_bench.py`、`tests/d2_driver.py`（新 `v14` 模式）透传。
- 测试 `test_d2_hf.py` 新增 8 项，现有 kill-resume 测试主体抽成 `_kill_and_resume` 复用：冻结与校验；v1.4（波次 1 与 4）每个 HS / HC 结局与 v1.3 串行逐项相同、按输入顺序结算、波次 1 时 journal 顺序与串行相同；S3 比例 0 / 1；额度只够一个 HS 时组间优先级；v1.4 在 HS 波中、HC 波中 kill → resume 只付一次。变异检验：结算顺序倒置（1 项失败）、三组顺序倒置（3 项失败）、S3 不跳过（1 项失败）；"判重放进并行阶段"测试台上无重复事件，测不出（缺口）。`test_ts_ef test_ts_branch test_ts_connect test_ts_warm_start test_d2_hf test_d2 test_d2_resume test_probe test_run_many test_scan_performance` **166 passed**；全包 `make check` 未跑。
- 测量脚本 `benchmarks/explore_d2/hs_parallel_bench.py`（11-04 的 22 + 15 项冻结括区，v1.3 HS 设置；串行臂走单计算器计数视图，pK 臂走 `run_many` + `MACEBatchEvaluator`；每秒采 `nvidia-smi`）。
- 冒烟（181，23:26 启动，PID 1781266，启动前负载 5.1、GPU 上有他人进程；每组 2 项，`results/2026-11-07/hs_parallel_bench/smoke/`，日志无 traceback）：串行 144 s、1 162 forward、0.124 s / forward；p4 91 s、1 510 forward、0.060 s / forward；两臂结局相同（3 候选、1 精阶段卡住）。p4 共 672 轮批量求值、平均每轮 2.2 个结构、求值共 81 s（每轮约 0.12 s，与串行单结构相同），求值外 9.8 s。GPU 利用率不可靠：串行臂期间他人进程占用约 10 GB 显存。
- 测量（用户 11-08 给 yayoi43；PBS 作业 24948，04:16:17 → 04:59:22，日志无报错；`benchmarks/explore_d2/hs_parallel_bench_pbs.zsh` 照模板只改 ← 行；结果 `results/2026-11-08/hs_parallel_bench/full/arm_*.json`，在 183 上读）：22 + 15 项冻结括区，v1.3 HS 设置。与计划不同：加了 p16（看饱和点）。
  - 串行 1 345 s、18.2k forward、0.074 s / forward、GPU 利用率均值 17%（p90 20%）；p4 475 s（2.8×）、16.8k、0.028 s、28%（36%）、平均批 3.9、求值 356 s / 求值外 118 s；p8 310 s（4.3×）、16.3k、0.019 s、37%（55%）、批 7.1、208 / 102 s；p16 337 s（4.0×）、19.4k、0.017 s、40%（73%）、批 9.5、236 / 101 s。
  - 出候选 31 / 32 / 33 / 28；与串行结局不同的项 p4 4、p8 5、p16 10，全部在失败组；两边都出候选的候选能量差 ≤ 0.0039 eV（2 个 fp32 量子）。p16 的差异偏多，未重复运行，不能区分 fp32 波动与批量求值的噪声差别。
  - 读法：串行时 GPU 很空；p8 墙钟最短；求值外时间在 p4–p16 停在约 100–118 s（线程下 Python 侧串行化，推断），GPU 最高约 40%，未饱和。
- v1.4 b05 冒烟（183，04:05 启动，PID 519580，启动前负载 0.00、GPU 空闲；`run_b05_d2.py` 新增 `--sample-parallel`；参数同 10-26 冒烟但 `--audit-fraction 1.0`（否则唯一的同键图脊多半被跳过、并行路径走不到）；`results/2026-11-08/b05_hf_v14_smoke/`，04:17 结束，日志无 traceback，在 183 上读）：D2-HF-v1.4；search 5 217 + refine 1 692 = ledger 6 909，unreconciled 空；升级 12 项（audit 11、ridge 1），HS 按波次派发 4 / 3 / 3 / 1 / 1（refine 12 000 / hs_cap 3 000 限制同时预留数），同键图的脊在最后一组；HS 结局 `hysteresis` 7、`above_ceiling` 5，无候选。
- 整跑 PBS 脚本 `benchmarks/explore_d2/b05_hf_v14_pbs.zsh`（用户 11-08：并发「8」、yayoi05「装得下」；照 v1.3 脚本只改 ← 行：yayoi05 2×3090、ncpus 32 / ngpus 2 / 120 GB，两块各占一张卡、各 16 线程，`--sample-parallel --max-active 8`，种子与额度同 v1.3）。
- 批量求值噪声（用户 11-08「A 跑」）：脚本 `benchmarks/explore_d2/batch_noise.py`。第一次在 183（08:13）因他人作业（PID 520460，`experiments/2026-11-08_ppi_r4/`，08:15 起占 8.3 GB 显存）CUDA OOM，未出结果；改在 181（08:20，GPU 空闲）跑完，结果 `results/2026-11-08/batch_noise/noise.json`。11-07 验收的 4 个对照候选结构，各单独算 8 次、作为批大小 1 / 8 / 16 的首个结构各算 8 次：力分量重复标准差最大值单独 0.0004–0.0011、批量 0.0004–0.0010，中位 6–8e-5 对 6.6–7.9e-5；批量与单独平均力之差逐原子最大 0.00014–0.0011 eV/Å（与 8 次重复平均值的统计误差同量级）；能量差 ≤ 0.0005 eV。未见批量求值噪声更大或有系统偏差；并发 16 的 10 项结局差异不能由此解释。限定：4 个结构，都是对照组候选，失败组卡住的帧没测。
- B / C 的 PBS 脚本（用户 11-08「BC 我全给你3090 跑」，yayoi03 / yayoi08）：`benchmarks/explore_d2/hs_parallel_bench_rep_pbs.zsh`（并行测量重复两次，各占一张 3090，输出 `results/2026-11-08/hs_parallel_bench/full_r{1,2}/`）；`benchmarks/explore_d2/b05_hf_v13r_pbs.zsh`（v1.3 在 3090 上重跑，布局与参数照搬 yayoi05 v1.4：两块各占一张卡、16 线程、`--max-active 8`、同种子同额度；输出 `results/2026-11-08/b05_hf_v13r_s{A,B}/`）。与先前所说不同：C 原为 yayoi32 上的 v1.4，节点是 3090，改为 v1.3 同硬件重跑，使 v1.3 / v1.4 只差规则。是否已提交未核实。
