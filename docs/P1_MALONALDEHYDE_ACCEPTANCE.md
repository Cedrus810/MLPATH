# P1 丙二醛验收判据（跑前冻结）

冻结时间：2026-09-01，**在第一次搜索运行之前**。照 §5 的规矩：判据先写死，不边跑边改。
体系：丙二醛烯醇（3-羟基丙烯醛），C3H4O2，中性闭壳，MACE-OFF24_medium。

## 这个基准要验证什么

P0 乙醇只验证了构象发现与架构，`reactions` 一直是 0 —— **"扰动-响应能找到化学反应"这个核心
主张至今没有证据**。P1 补的正是这一条，并且它同时检验四件事：

1. PRRS 真能找到 **bond-changing** 事件；
2. 反应检测**不依赖已知产物**；
3. 自同构商**不制造假产物**（不因退化重排开新节点）；
4. 自同构商**也不把真实反应事件吃掉**（§8.11 的第三层）。

第 4 条是这个体系被选中的原因：丙二醛的质子转移是**退化重排**，两个互变异构体的
`chemical_key` 相同（preflight 实测，§8.11.4），所以它只能被 `reaction_channels` 的
self-loop 捕获，捕不到就说明第三层没起作用。

## 必须通过（跑前冻结）

```text
A. 结构与网络
   A1  chemical_nodes == 1
   A2  reactions == 0                     普通跨节点边应当没有（退化反应不跨节点）
   A3  没有任何 microstate 的 chemical_key 与源不同

B. 反应事件（这个基准的核心）
   B1  len(reaction_channels) >= 1
   B2  存在一条通道满足全部：
         class     == "degenerate_reaction"
         self_loop == True
         ends      == [c0000, c0000]
         broken 恰含一条 O-H 键
         formed 恰含一条 O-H 键，且与 broken 的那条不是同一条
         broken 与 formed 涉及同一个 H
   B3  该通道的 event key 正反一致（由构造保证，作为回归项记录）

C. 不得发生（负判据，同等重要）
   C1  任何 bond_delta 为空的事件都不得进入 reaction_channels
   C2  reaction_channels 里不得出现 class == "key_change_without_bond_change"
   C3  不得出现自由基/开壳层产物：任何 microstate 的 fragments 必须是闭壳中性式
       （出现 C3H3O2 之类奇电子碎片即判失败，说明跑出了 MACE-OFF 的域）

D. 曲率层（沿用第一、二批的门）
   D1  所有被采纳的 TS 满足 saddle_order == 1
   D2  所有 curvature_source == "analytic"
   D3  所有 trivial_mode_floor_worst_residual < 5e-4
   D4  若有 TS 支持该反应通道，其 response.sign_stable == True
```

## 观察项（不作判据，只记录）

- 该通道被哪个家族、哪个幅度首次命中（`best_by_family` 语义：**最小观测成功幅度，不是势垒**）；
- TS 的 $\hat\kappa$、$\sigma$、$\hat c$ —— 这是第一次在**成键/断键区**测这些量；
- D1 门的残差在质子转移 TS 上的大小（此前只在乙醇的构象势垒上测过）；
- 是否同时出现构象通道（烯醇的 C–C 转动等），以及它们是否正确留在 `conformer_transitions`。

## 预期的失败模式（先写下来，免得事后合理化）

1. **搜不到**：幅度不够或家族覆盖不到 O–H 方向 → 表现为 `reaction_channels == []`，
   A/C/D 全过。这是**判据失败**，不是"体系不对"。
2. **跑过头**：大幅度 kick 把分子拆了 → C3 触发。
3. **TS 找不到但产物找得到**：B 过、D 空。记录为部分通过，说明 TS pool 的命中率问题。
4. **质子转移被判成构象**：C1 触发 —— 说明第三层接线有误。

## 配置

见 `runs/p1_malonaldehyde/config.json`（随运行写出）。相对 P0 的改动只有：
`families` 增加 `compress`/`bend`、`geometry_amplitudes_A` 增加 0.5、`max_trials` 提到 60。


---

# 运行结果（2026-09-01，判据冻结后第一次运行）

`runs/p1_malonaldehyde/`（seed 17，60 次试验）。**判据一条没有改过**；下面记录的是判据
原样评估的结果，包括我写得不够精确的那一条。

## 7 / 12 通过

```
FAIL  A1  chemical_nodes == 1                             got 2
FAIL  A2  reactions == 0                                  got 1
FAIL  A3  no microstate under a different key             ['005c151a', 'ad42b1d9']
PASS  B1  reaction_channels >= 1                          got 2
PASS  B2  degenerate self-loop channel                    rc0001
PASS  B3  event key recorded                              a04410b76e16b117
FAIL  C1  no empty-bond-delta event in reaction_channels  ['rc0000']
FAIL  C2  no key_change_without_bond_change               ['rc0000']
PASS  C3  every fragment closed-shell
PASS  D1  every TS saddle_order == 1                      0 TS
PASS  D2  every curvature_source analytic                 0 TS
PASS  D3  trivial floor residual < 5e-4                   0 TS
```

## B 通过：PRRS 第一次找到了化学反应

```json
{"id": "rc0001", "class": "degenerate_reaction", "ends": ["c0000", "c0000"],
 "self_loop": true, "broken": [[4, 5]], "formed": [[0, 5]],
 "event_key": "a04410b76e16b117", "observed_by": ["t000022", "t000023"]}
```

`broken` 是 O2–H、`formed` 是 O1–H，event key 与 preflight 预测的**逐字符相同**。
命中它的两次试验都是 **`stretch` 家族作用在 O2–H（indices [4,5]）**，幅度 **0.5 和 0.8 Å**
—— 0.8 这一档是量到"质子要走 0.76 Å"之后才加的，0.5 也命中了。

60 次试验的分类分布：`conformational_transition` 36、
`key_change_without_bond_change` 6、`degenerate_reaction` 2；
结局 44 completed / 15 rejected / 1 unconfirmed。

**这是"扰动-响应能找到化学反应"第一次有证据**，而且它只能被第三层的 self-loop 捕获 ——
改动之前它会被塞进 `conformer_transitions`。

## C2 失败：负判据抓到 `chemical_key` 的一个真实缺陷

那个"原则上不可能"的第四格触发了 6 次，并且诊断正确。

`chemical_key` 的 `locked_bond_parity` 用**一个硬几何阈值**决定哪些键算锁定
（`chemistry.py:146`）：

$$\text{locked} \iff d_{ab} < 0.93\,(r_a + r_b)$$

C–O 的阈值是 $0.93 \times (0.76+0.66) = 1.3206$ Å。而烯醇的 C3–O2 就坐在上面：

| 结构 | C1–O1 | C3–O2 | C1–C2 | C2–C3 | locked 集合 | key |
| --- | --- | --- | --- | --- | --- | --- |
| c0000/m0000 | 1.2283 | **1.320812** | 1.4430 | 1.3545 | C1-O1, C2-C3 | 005c151a |
| c0001/m0000 | 1.2286 | **1.320495** | 1.4426 | 1.3548 | C1-O1, **C3-O2**, C2-C3 | ad42b1d9 |

两个结构**每根键都相差不到 0.0004 Å**、能量差 0.07 meV，却因为 C3–O2 分别在阈值
**+0.00021 Å** 与 **−0.00011 Å** 两侧，落进了不同的 chemical node。
A1/A2/A3 的失败全部是这一个缺陷的后果（多出的节点 c0001，以及它带来的一条假跨节点边）。

**这不是调阈值能修的。**同一族结构里 C3–O2 从 1.2280 连续变到 1.3540，任何静态阈值都会被
它扫过。**"这根键锁不锁"和 §8.8.2 的"这个坐标自不自由"是同一个问题的两端 ——
后者已经确立必须靠测量而不是靠图猜，前者现在还在靠一个键长猜。**

## D 空过：产物找到了，TS 没找到

60 次试验产出 **0 个 TS candidate**，所以 D1–D3 是空过，D4 无法评估。这正是冻结时写下的
预期失败模式第 3 条（"TS 找不到但产物找得到"），记为**部分通过**。

## 我写得不够精确的一条判据（不追改冻结文件）

C1 写的是"bond_delta 为空的事件不得进入 `reaction_channels`"，但第四格
`key_change_without_bond_change` 是**故意记录下来供诊断的**，所以它进 `reaction_channels`
是设计如此。C1 把"不得算作反应"和"不得被记录"混为一谈了。
**冻结的文件保持原样**；下一次冻结时这条应写成"不得出现 `class == conformational_transition`
的条目"，并把第四格单列为诊断项。


---

# 追查：60 次试验为什么一个 TS 都没抓到

## 先确认鞍点存在

约束 $d(\mathrm{O_1H}) - d(\mathrm{O_2H}) = 0$ 后松弛（`ase.constraints.FixInternals` 的
`bondcombos`），再去掉约束算解析 Hessian
（`docs/experiments/`，结构存在 `runs/p1_malonaldehyde/p1_ts.extxyz`）：

```
O1-H = O2-H = 1.1910 A     O1...O2 = 2.3311 A   （从极小的 2.6099 收缩）
fmax（去掉约束后）= 0.0019 eV/A            -> 真驻点
E = -7274.130236 eV        势垒 = 419.14 meV
negatives = 1              lambda_1 = -37.323      icm = -3186 cm^-1
trivial floor residual = 5.89e-16
```

**所以是机制失败，不是鞍点不存在。**

朴素的笛卡尔中点**不是** TS（高 1.03 eV、两个负模，最低 −49.5 = 3670 cm⁻¹ 的 O–H 伸缩）：
中点把 H 放在两个 O 中间却没让 O···O 收缩，而 **O···O 从 2.61 收到 2.33 是反应坐标的必要部分**。

顺带一条模型精度观察：MACE-OFF24 给的势垒是 **419 meV**，而本文档冻结时凭印象引的文献值是
约 170 meV。这是 2.5 倍的差距，也正是 §10.3 第 1 条说的"平滑偏差只有外部参照能查"的实例
——**这个数字本身没有被独立验证过**，不要当作已知量使用。

## 三条独立的原因

### 1. 反应坐标是刚性的，脊比乙醇陡 120 倍

命中反应的两次是 `stretch` 作用在 O2–H，注入 **1.59 / 1.41 eV** 去做一个 **0.42 eV** 的反应。
更要命的是曲率：

| | $\lambda_1$ | icm |
| --- | --- | --- |
| 乙醇 syn 势垒（P0 命中的那类） | $-0.311$ | $-291$ |
| 丙二醛质子转移 | $-37.3$ | $-3186$ |

从势垒顶上方 1.5 eV 出发的淬火，落在陡 120 倍的脊上的概率实际为零。
P0 之所以能抓到 TS，是因为命中的家族是 `torsion` —— **软坐标**，注入能量与 47 meV 势垒同量级。

$$\boxed{\text{"淬火停在鞍点"对软反应坐标有效，对刚性反应坐标失效}}$$

### 2. 最接近分界面的那一档被整个丢弃

O2–H 上的四档：

```
0.15  -> conformational_transition          dE = 0.30 eV
0.35  -> unconfirmed (transient_topology)   dE = 1.28 eV   <- 质子确实转移了
0.50  -> degenerate_reaction                dE = 1.59 eV
0.80  -> degenerate_reaction                dE = 1.41 eV
```

0.35 那次越过了势垒，但新图没撑过自由 MD 尾部三步，于是 `run_trial` 返回**不带 endpoint** 的
`Outcome`，`confirm_minimum` 根本没跑到它。**最可能停在脊上的幅度，连 TS 检验的门都没进。**
而且它不算"确认成功"，所以幅度括号没打开、细化没发生。

"新图有没有撑过自由 MD 尾部"是怀疑**反应主张**的合理理由，但不是丢弃**淬火端点**的理由 ——
那个端点是个完好的驻点。

### 3. `refinement_steps = 1`

整轮唯一打开的幅度括号是 `[0.425, 0.5]`，还在另一个坐标上。

## 外加一条独立的预算泄漏

15/60 被 `collateral_bond` 挡掉，其中 **12 次是 `compress` 作用在非键对**上
（`[7,8]` = H···H，`[4,7]` = O···H）。continuation 交付把共价键拖动 0.08–0.41 Å，
而预算是 0.05 Å。**门是对的，交付满足不了它** —— 20% 的预算花在对这个分子无法交付的家族上。

## 由此产生的四项待办

1. **`transient_topology` 不应连 endpoint 一起丢弃。**仍标 `unconfirmed`，但把端点交给
   `confirm_minimum`。直接命中本次失败的那个幅度。
2. **`unconfirmed` 的穿越也应开幅度括号**，让细化朝分界面收敛；`refinement_steps` 调大。
3. **在已知端点之间搜 TS**，而不是指望淬火停下来。反应通道一旦有两个端点，
   `lowest_response_direction`（§8.8.3，每步 2 次力评估）+ 解析 Hessian 就够做 dimer 式攀爬。
   这让 TS 发现**不再依赖 probe 幅度**，是通用解法。
4. **非键接近的交付要尊重 collateral 预算**：把共价键伸缩分量投影掉，或在交付不了时跳过该
   方向而不是烧掉试验。


---

# 第三轮运行（identity 层修复 + 保留 transient crossing 之后）

`runs/p1_round3/`，同一输入、同一 seed、`refinement_steps` 由 1 提到 3。
**判据仍然一字未改。**

## 12 / 12，但 D 组是空过

```
PASS  A1 chemical_nodes == 1        got 1     <- 第一轮 got 2
PASS  A2 reactions == 0             got 0     <- 第一轮 got 1
PASS  A3 no microstate under a different key   ['420aa409']
PASS  B1/B2/B3   rc0000 degenerate_reaction self_loop=True
                 broken=[[4,5]] formed=[[0,5]] event_key=a04410b76e16b117
PASS  C1/C2      []  []                       <- 第一轮第四格触发 6 次
PASS  C3
PASS  D1/D2/D3   0 TS —— 空过；D4 无法评估
```

A 组从全灭到全过，C 组的第四格从 6 次降到 0 次：**bond-order 本体论（§8.13）修掉了它**。
B 组保持，同一通道被观测 **4 次**（第一轮 2 次），event key 与第一轮、与 preflight 逐字符相同。

## transient crossing 确实打开了 bracket

```
stretch      [4,5]      0.250 -> 0.275   transient_crossing
torsion      [0,1,2,7]  1.700 -> 1.830   confirmed
torsion      [2,3,4,5]  1.440 -> 1.570   confirmed
torsion_kick [2,3,4,5]  0.300 -> 0.387   confirmed
```

O–H 的完整阶梯（含细化）：

| 幅度 | 阶段 | 结果 | $\Delta E_{\rm geom}$ | `saddle_order` |
| --- | --- | --- | --- | --- |
| 0.150 | coarse | 弹性 | 0.302 | 0 |
| **0.250** | refine | 弹性 | 0.749 | 0 |
| **0.275** | refine | 反应 | 0.887 | 0 |
| 0.300 | refine | 反应 | 1.027 | 0 |
| 0.350 | coarse | transient | 1.278 | 0 |
| 0.500 | coarse | 反应 | 1.592 | 0 |
| 0.800 | coarse | 反应 | 1.409 | 0 |

第一轮里 0.15–0.35 之间一次都没探过；现在分界面被夹到 **0.025 Å 宽**，
注入能量从 1.59 eV 降到 0.89 eV。注意 **0.275/0.300 已确认反应而 0.350 只是 transient**
—— 不单调，所以 bracket 的 `meaning` 写的"no monotonicity assumed"是必要的。

**但每一档的 `saddle_order` 仍然是 0。**幅度已经收到最优，quench-to-saddle 依然不响 ——
这是"机制问题不是幅度问题"的决定性证据。

## min-mode following 从 transient 帧爬到了同一个鞍点

用三个 `raw_endpoint`（自由 MD 之后、淬火之前）作种子跑 `follow_min_mode`：

| 种子 | $E - E_{\min}$ | 结果 |
| --- | --- | --- |
| 0.250 弹性 | 69.3 meV | `wrong_index`，7 步收敛到 order 0 |
| 0.275 越过 | 493.5 meV | `wrong_index`，9 步收敛到 order 0 |
| **0.350 transient** | **796.3 meV** | **8 步到达鞍点** |

```
fmax = 0.00080     order = 1     lambda_1 = -37.32855     icm = -3186.0
barrier = 419.14 meV             floor residual = 7.1e-16
O1-H = 1.1910   O2-H = 1.1910   O1-O2 = 2.3311
```

**与独立构造的 ground truth（约束 $d(\mathrm{O_1H})-d(\mathrm{O_2H})=0$ 松弛）逐位相同。**
两个完全不同的方法给出同一个鞍点。

### 由此得到一条比"把 bracket 夹得更紧"更准确的结论

Newton 步收敛到**最近**的驻点。0.25 的帧只高出 69 meV，最近的驻点是极小；
0.275 的帧高出 494 meV，最近的仍是极小；**只有 0.35 那个 transient 帧（796 meV、拓扑正在
切换）最近的驻点才是鞍点**。

$$\boxed{\text{有用的种子是 transient-crossing 的原始帧，不是 bracket 的端点}}$$

**而第一轮把这一帧连同端点整个丢掉了。**所以"保留 transient crossing"不只是省下一次试验，
它保住的是唯一能爬到鞍点的那个几何。

## 剩下的就是接线

`follow_min_mode` 已实现并有测试，种子来源已经确认。尚未做的：
在观测到 reaction event 之后，用该试验的 transient/近分界面帧作种子调用它，
把结果送进 TS pool 并接上现有的端点下降。**这一轮 TS 仍是 0 是预期之内的。**


---

# 第五、六轮：TS 发现接线完成，13/13 且 D 组不再空过

判据仍然一字未改。`runs/p1_round5/`、`runs/p1_round6/`，同输入同 seed。

## 结果

两轮都是 **13/13**，且 **D 组第一次不是空过**：

```
ts0000  found_by=min_mode_following  saddle_order=1
        icm = -3189.6            barrier vs m0000 = 419.14 meV
        curvature_source=analytic  trivial floor residual = 2.8e-16 / 6.2e-16
        response  kappa = -36.75 +- 0.45   sign_stable=True  significant=True
        connects  (+1, c0000/m0002, existing_microstate)
                  (-1, c0000/m0000, existing_microstate)
channel rc0000  degenerate_reaction  self_loop=True  ts_support=['ts0000']
```

**这是本文档最关键的一条：质子转移通道由两侧下降独立支持。**两侧落在 m0002 与 m0000
（两个互变异构体），所以 `ts_support` 是**下降建立的事实**，不是发起爬山的意图。
icm $-3189.6$ 与独立构造的 ground truth $-3186.0$ 差 0.1%，势垒 419.14 meV 逐位相同。

## 三个分母不同的指标

单看"有没有 TS"会掩盖种子选错通道的问题，所以三项分别报
（`docs/experiments/saddle_metrics.py`，可对任意 run 目录运行）：

| | round 5 | round 6 |
| --- | --- | --- |
| `ceiling_blocked` | 0（当时是 `above_ceiling`，一步越界即终止） | **1** |
| `seek_success_rate`（一阶鞍点 / seek 次数） | 1/2 = 0.50 | 1/2 = 0.50 |
| `attribution_precision`（命中 / 可完成归属） | 1/1 = 1.00（unknown 0） | 1/1 = 1.00（unknown 0） |
| `requested_channel_coverage` | 1/1 = 1.00 | 1/1 = 1.00 |

```
transient  attempts=1 found=0 hits=0 unknown=0     (0.35)
confirmed  attempts=1 found=1 hits=1 unknown=0     (0.275)
```

**成功率相同而归因相反**：round 5 的失败是我的实现缺陷（一步过冲即终止），
round 6 的失败是一个起点过高的种子被界限正确拒绝。**单看 `seek_success_rate`
会把这两轮当成一样** —— 这正是要三个不同分母的理由。

## ceiling 的机制修复，以及它把结论改了

原实现在一步越过 ceiling 后**下一轮才检测并直接终止**，把一次可回退的步长过冲变成整个 seek
失败 —— 正是它砍掉了那个唯一到过质子转移鞍点的 transient 种子。改成四条机制
（`saddle_search_max_rise_eV` 保持 **1.5 eV** 不变）：

1. 种子起点即越界 → `seed_above_ceiling`，**跳过换下一个种子**（判决关于种子，不该被继承）；
2. 爬山步越界 → **步长减半回溯**（`min_mode_backtracks = 8`），找不到可接受步长才
   `ceiling_blocked`；
3. **收敛后的鞍点仍必须低于 ceiling**（`above_ceiling_at_saddle`）—— 否则界限可以用更小的步
   绕过去，而界限是关于**区域**的；
4. 记录 `initial_rise_eV` / `max_rise_eV` / `ceiling_hit_step` / `backtracks`。

第 4 条立刻还了债。round 6 的 transient 种子：

```
t000022  initial_rise_eV = 1.2782   max_rise_eV = 1.4998
         backtracks = 15   ceiling_hit_step = 2   -> ceiling_blocked
```

质子转移鞍点在 0.419 eV，而这个种子起点 **1.278 eV** —— 已在目标鞍点上方 0.86 eV，
距界只剩 0.22 eV。回溯 15 次（步长降到 $2^{-15}$）仍找不到不越界的方向。
**所以这一次的 `ceiling_blocked` 是界限正确工作，不是上一轮那个 bug。**
两者的区别完全由这三个字段分开；没有它们，两种情况显示成同一个失败。

顺带解释了为什么命中的是 0.275 那个 **confirmed** 种子：它的 `initial_rise` 约 0.887 eV，
余量 0.6 eV，够走。**"transient 才是好种子"这个结论下得太早** ——
帧的选择（`topology_switch:perturbed:1` vs `raw_endpoint`）比 transient 标签更起决定作用。

## 三个语义边界（已实现，对后续运行生效）

- `requested_channel` 是**搜索意图**；
- `supports` 是**两侧下降得到的事实**；
- `supports_requested` 是**种子策略的命中统计**，不参与事实归属；
- 任一侧下降失败时 `supports_requested` 记 **null 而非 false** —— false 是一个主张
  （端点证据完整、但支持了别的通道），在证据缺失时断言它等于把"没测"变成"负结果"。
  另记 `attribution_complete`。

## 中间交付帧：按判据不动

判据是"seek 成功率高但 `supports_requested` 明显低 → 补帧"。实测相反：
`attribution_precision` 与 `coverage` 都是 1.00，唯一失败是种子起点过高。
补中间帧会给出 `initial_rise` 低得多的种子（连续化步骤在 0.15–0.35 之间），
所以它能提高 `seek_success_rate` —— 但那是**效率**问题不是正确性问题。
真要做，只加只读轨迹钩子（最后一个未穿越帧 + 第一个穿越帧），不改 stepper。


---

# 第七轮：P2 round 4 之前的回归（同 seed、同输入、四项干预）

跑这一轮的唯一目的是**先确认 P2 round 4 的四项干预没有动坏已经过的东西**，再去跑 P2。
判据一个字未改，`docs/experiments/p1_check.py` 也未改。

干预：I1 `free` 档改为"平且驻" + `flat_biased` 扫描；I2 步长上限 $\pi/n_{\rm sym}$、
`soft_polish_rounds` 3→6；I3 `quench` 必须"打磨既满意又没动"才算收敛；
I4 boundary continuation。详见 PRRS_STATUS §8.18–8.19。

## 13/13，而且关键量逐字段不变

```
PASS A1 A2 A3   PASS B1 B2 B3   PASS C1 C2 C3   PASS D1 D2 D3 D4      13/13
```

| | round 6 | round 7 |
| --- | --- | --- |
| `event_key` | `a04410b76e16b117` | **同** |
| ts0000 `icm` | $-3189.6$ | **同** |
| ts0000 $\kappa$ / $\sigma$ / $c$ | $-36.751254$ / $0.4470498$ / $1433.0669$ | **逐字段同** |
| TS 能量 | $-7274.13024$ | **同** |
| 势垒（自 m0000） | 419.14 meV | **419.14 meV** |
| `chemical_key` | `420aa409…` | **同** |

**TS 数 1 → 4**，四个都是 index-one、`curvature_source=analytic`、刚体地板残差
$\le 4.2\times10^{-16}$，虚频 $-3189.6/-3188.2/-3187.3/-3186.7$ cm⁻¹
（DFT/实验区间的参考值 $-3186$）。
**通道的 `observed_by` 从 4 个试验涨到 17 个。**

`reactions` 仍为 0 —— 这是**对的**：P1 的质子转移是退化的，两端是同一个化学态，
所以它是自环通道而不是有向边。A2 判的就是这件事。

## boundary continuation 在 P1 上的表现（真值已知，所以是好的检验）

两次续行，两次都到达退化产物，而且**逐档分辨**：

```
seed t000001 (transient, amp 0.35, tangent=topology_switch_difference)
  rung 1 @ 0.050 A  returned_to_source
  rung 2 @ 0.100 A  reached_degenerate_product          <- 到这里才进入
  q_PT: +0.063 -> -0.749

seed t000040 (transient, amp 1.0)
  rung 1..3 @ 0.05/0.10/0.15 A  returned_to_source
  rung 4 @ 0.200 A  reached_degenerate_product
  q_PT: -0.192 -> +0.741
```

两个种子都是**产物声明被扣下的那两次试验**（`unconfirmed` / `transient_topology`,
`target=None`）。也就是说：续行把那两次被扣下的试验里真正包含的信息取回来了，
而且没有把"扣下"改成"接纳" —— 那两次试验的记录仍然是 `target=None`。

`q_PT` 的轨迹是这件事第一次可以直接读出来：$+0.063 \to -0.749$，质子确实换了氧。

## 与 round 6 的差异，逐条说清楚

- **microstates 4 → 2。**不是塌缩。round 6 的 m0001/m0003 是**真实的**开链（非螯合）
  烯醇构象，在新代码下重新淬火仍收敛、仍是 order 0，位于 $+359.06$ / $+359.84$ meV
  （O0–H5 = 3.79 / 5.02 Å）。round 7 换成了 m0001 = **镜像螯合烯醇**（$+0.08$ meV，
  O0–H5 = 0.987、O4–H5 = 1.752），也就是退化产物本身。两者都合法，判据都不涉及；
  这是端点一变、试验序列就分叉的正常后果。
- `rejected` 14 → 8，`completed` 45 → 50，`unconfirmed` 1 → 2。
- 分类分布从 `{elastic 33, unphysical 14, rearrangement 8, reactive 5}` 变为
  `{elastic 35, reactive 17, unphysical 8}`。

**这一轮不声称是改进**，它只声称四项干预没有破坏 P1 的任何一条冻结判据，
并且在真值已知的体系上验证了 I4 的机制。


## 第八轮：对最终代码的重跑，与第七轮逐字段相同

第七轮之后代码还改了两处：两侧下降的中继逻辑从 `connect_saddle` 提取成
`descend_saddle`（两处共用），以及续行判决回调改成纯函数。提取本应行为等价，
但"本应"不是测量，所以同 seed 同输入重跑一次（`runs/p1_round8`）。

**13/13**，且与第七轮逐字段相同：

```
attempts  SAME     outcome_counts  SAME     reactions  SAME
ts        SAME  (4 个，icm -3189.6/-3188.2/-3187.3/-3186.7，E=-7274.13024)
microstates SAME  (-7274.549375 / -7274.549293，free_bonds 均为空)
channels  SAME  (rc0000 degenerate_reaction, observed_by=17, 4 个 TS 支持)
event_key a04410b76e16b117
```

所以那两处改动在 P1 上没有可观测影响，提取是等价的。
（刚体地板残差从 $6.2\times10^{-16}$ 变到 $3.4\times10^{-16}$ —— 那是
`confirm_minimum` 的探针 seed 依赖项，量级远低于判据的 $5\times10^{-4}$。）
