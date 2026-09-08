# P2 3-氧代丁醛烯醇验收判据（跑前冻结）

冻结时间：2026-09-01，**在第一次搜索运行之前**，且在 preflight 结果之前定稿。
体系：3-氧代丁醛的螯合烯醇 CH₃–C(OH)=CH–CHO，C4H6O2，中性闭壳，MACE-OFF24_medium。

## 这个基准要验证什么

P0 只验证构象发现与架构；P1 验证了 **reaction-event 发现**，但它的质子转移是**退化重排**，
两端 `chemical_key` 相同，所以只能证明 self-loop 通道可被捕获，**证明不了反应网络会扩张**。

P2 补的就是这一条：

$$C_A \rightarrow C_B,\qquad k_A \neq k_B$$

把丙二醛一端的 H 换成甲基后，两个氧不再能由分子自同构交换，质子从一端转到另一端**必须**
产生一个新的 chemical node。它仍然是 C/H/O、中性、闭壳、β-二羰基质子转移化学 ——
比为了测试反应发现硬上带电 SN2 干净得多（SN2 被 §2.5 的电荷阻塞挡住）。

**这个基准的验收对象是 `chemical_key` 本身**，所以它排在 identity 层稳定之后：
bond-order 本体论（§8.13）构象无关、P1 连续两轮 `chemical_nodes == 1`、域门有真实自由基
测试钉住。

## 必须通过（跑前冻结）

```text
A. 网络扩张（这个基准的核心）
   A1  chemical_nodes >= 2
   A2  存在两个节点，chemical_key 不同而 fragments 相同（都是 ['C4H6O2']）
       —— 变的是图不是碎片，这正是 P1 无法产生的情形
   A3  len(reactions) >= 1，且该有向边 source != target
   A4  存在一条 reaction_channels 条目满足全部：
         class     == "reaction"
         self_loop == False
         ends      是 A2 那两个节点
   A5  该通道 broken 恰含一条 O-H，formed 恰含另一条 O-H，且涉及同一个 H

B. 不得发生（负判据，同等重要）
   B1  reaction_channels 里不得出现 class == "key_change_without_bond_change"
   B2  所有 microstate 的 fragments 必须闭壳；不得出现自由基碎片
       （域门应当在入网前就拦住，出现即说明门失效）
   B3  bond_delta 为空的事件（甲基旋转、羟基转子等）不得出现在 reaction_channels 里

C. 曲率层（沿用 P1）
   C1  所有 TS saddle_order == 1
   C2  所有 curvature_source == "analytic"
   C3  所有 trivial_mode_floor_worst_residual < 5e-4
   C4  若有 TS 支持 A4 那条通道，其 response.sign_stable == True

D. TS 与归因
   D1  至少一个 TS 由两侧下降支持 A4 那条跨节点通道
   D2  该 TS 的两侧下降落在两个**不同**的 chemical node
   D3  归因三态语义生效：任一侧下降失败的 saddle 记 supports_requested = null，
       不得记成 false
```

## 观察项（不作判据，只记录）

- 三个分母不同的指标：`seek_success_rate`、`attribution_precision`、
  `requested_channel_coverage`，按种子类型与幅度分开（`docs/experiments/saddle_metrics.py`）；
- 首次命中该反应的家族与幅度（`best_by_family` 语义：**最小观测成功幅度，不是势垒**）；
- TS 的 $\hat\kappa$、$\sigma$、$\hat c$、势垒、虚频；
- 是否同时出现构象通道（甲基转子、羟基转子），以及它们是否正确留在
  `conformer_transitions`；
- `ceiling_blocked` / `seed_above_ceiling` 的次数与 `initial_rise_eV` 分布。

## 预期的失败模式（先写下来，免得事后合理化）

1. **搜不到**：幅度或家族覆盖不到 O–H 方向 → A3/A4 失败，A1 可能仍过（构象节点不算）。
   这是**判据失败**，不是"体系不对"。
2. **只找到构象**：产物是同一节点的另一个 microstate（羟基转子、甲基转子），
   `reactions` 仍为 0 → A3 失败。
3. **跑出域**：大幅度探针拆出自由基 → B2 触发，或域门把端点拒了导致 A1 达不到 2。
4. **TS 归因到别的通道**：D1/D2 失败但 A/B 通过 —— 说明种子策略选错通道，
   记 `supports_requested = false`（证据完整时）而不是掩盖。
5. **新节点其实是同一物种**：A2 通过但两个 key 的差别来自立体化学而非图 ——
   §8.12 那一类缺陷的复发。A2 明确要求 `graph_hash` 也不同。

## preflight 门（零搜索预算，必须先过）

```text
P1'  两个互变异构体各自 optimize 后 n_minus == 0
P2'  chemical_key 不同，且 graph_hash 也不同
P3'  两者 fragments 都是 ['C4H6O2']、stereo_unresolved is None（域内、价态可解）
P4'  reaction_event_key 正反一致，broken/formed 各一条 O-H
```

**这四条不过就不跑**，验收对象本身还在动的时候烧搜索预算没有意义。

## 配置

见 `runs/p2_oxobutanal/config.json`（随运行写出）。相对 P1 的改动只有分子本身；
`families`、`geometry_amplitudes_A`、`max_trials`、`refinement_steps`、seed 全部保持一致，
以便与 P1 的指标直接比较。


---

# 源结构准备与 preflight（判据冻结之后，正式搜索之前）

留档在 `runs/p2_oxobutanal/`（`p2_input.extxyz` = A、`p2_product_reference.extxyz` = B、
`p2_A_saddle.extxyz` = 最初那个鞍点、`p2_source_provenance.json`、`p2_preflight2.json`），
脚本在 `docs/experiments/p2_preflight.py`、`p2_fix_a.py`、`p2_preflight2.py`。

## 第一次 preflight：P1' 失败

```
A (质子在 O2)  fmax=1.31e-03  n_minus = 1  lambda_1 = -0.067119   <- 不是极小
B (质子在 O1)  fmax=9.27e-04  n_minus = 0
```

**措辞要准确**：`fmax` 小只说明接近**某个**驻点，不说明它是极小。A 很可能就在一个真实的
一阶鞍点附近，而不是"还没走到驻点"；$\lambda_1 = -0.067$ 比零模噪声（$10^{-15}$）大 13 个
数量级，不是数值残渣。这正是 §3.5 那条原则的实例，也正是 P1' 这道门存在的理由 ——
它拦住了一次"用鞍点当反应物"的搜索。

## 源结构准备：沿负模两侧下降

同一个 calculator、同一个 Hessian 来源（analytic）、同一套松弛配置（LBFGS，`fmax = 2e-3`）；
位移按**最大单原子位移**归一化，幅度从小到大 0.05 / 0.10 / 0.20，两侧各一次；
**每个端点都跑 `confirm_minimum`**，不按优化器的 `fmax` 判定 —— 那正是产生错误判断的那个量。

```
delta   E (eV)          confirmed  order  key           proton   O1-H   O2-H
+0.05   -8345.351864    True       0      afd590f8e107  [[4,5]]  1.693  0.991
-0.05   -8345.351864    True       0      afd590f8e107  [[4,5]]  1.693  0.991
+0.10   -8345.351865    True       0      afd590f8e107  [[4,5]]  1.694  0.991
-0.10   -8345.351865    True       0      afd590f8e107  [[4,5]]  1.694  0.991
+0.20   -8345.351865    True       0      afd590f8e107  [[4,5]]  1.694  0.991
-0.20   -8345.351865    True       0      afd590f8e107  [[4,5]]  1.694  0.991
```

**六个端点全部收敛到同一个极小**，$\pm$ 两侧一致说明那个负模是对称的（甲基转子或面内对称模）。
鞍点比该极小高 **42.2 meV**。

取用之前验了四项身份：`chemical_key`、`graph_hash`、`fragments` / `stereo_unresolved`、
以及**质子挂在哪个氧**（`proton_bonds == [[4,5]]`，仍在 O2）。四项都没变，所以低能端点确实
还是 A，不是 B 冒充的。

**若两侧都落进 B 或别的化学态、域内不存在 A 极小，就应判定这个基准不成立**，
而不是继续加工出一个 A。脚本里是这个逻辑；这一次不需要它。

## 第二次 preflight：完整重跑 P1'–P4'

只重算 `n_minus` 是不够的 —— **构象变化把能量差改了 4.6 倍**：

$$\Delta E(B-A):\quad -53.96\ \mathrm{meV}\ \longrightarrow\ -11.72\ \mathrm{meV}$$

而且 B 这一次用同样的严格度复验（第一次只看了 `n_minus`，没走 `confirm_minimum`），
以保证两端是用同一标准准备的。B 一次通过，无需修复。

```
A  E=-8345.351865  confirmed=True order=0 floor=6.6e-16
   key=afd590f8e107  graph_hash=6b9cbfe364d83a0b  frags=['C4H6O2']  stereo=None
   proton=[[4,5]]  O1-H=1.694  O2-H=0.991           螯合，烯醇偏酮端
B  E=-8345.363587  confirmed=True order=0 floor=5.3e-16
   key=573a10fe35da  graph_hash=b097b7d49519f89f  frags=['C4H6O2']  stereo=None
   proton=[[0,5]]  O1-H=0.991  O2-H=1.697           螯合，烯醇偏醛端

PASS  P1'  both are confirmed minima
PASS  P2'  keys differ and graph hashes differ
PASS  P3'  both C4H6O2 and valence-resolved
PASS  P4'  event key direction-free, one O-H each way

broken=[[4,5]]  formed=[[0,5]]  event_key=b4c3ba7d108c52ad
classification: reaction      ALL GATES PASS
```

`fragments` 相同而 `graph_hash` 不同 —— 正是 A2 判据要的情形：**变的是图，不是碎片，
也不是立体化学**。

结构本身也核对过确实是 β-二羰基烯醇，双键位置随互变异构正确互换，甲基只挂在 C3 一端：

```
A: O1=C1 1.234 / C1-C2 1.432 / C2=C3 1.365 / C3-O2 1.320 / O2-H 0.991 / C3-CH3 1.496
B: O1-H 0.991 / O1-C1 1.322 / C1=C2 1.352 / C2-C3 1.452 / C3=O2 1.233 / C3-CH3 1.505
```

## 搜索配置

与 P1 **完全一致**（`families` 六个、`geometry_amplitudes_A` = [0.15, 0.35, 0.5, 0.8]、
`max_trials` = 60、`refinement_steps` = 3、seed = 17），只换分子，以便指标与 P1 直接比较。


---

# 第一轮运行结果：A 组全灭，原因是反应坐标从未被采样

`runs/p2_oxobutanal/`（seed 17，60 次试验，配置与 P1 完全一致）。**判据一字未改。**

## 逐条评估

```
FAIL  A1  chemical_nodes >= 2                  got 1
FAIL  A2  两节点 key/graph 不同                 只有一个节点
FAIL  A3  reactions >= 1                       got 0
FAIL  A4  跨节点 reaction channel               got 0
FAIL  A5  broken/formed 各一条 O-H              无通道
PASS  B1  无 key_change_without_bond_change     无通道
PASS  B2  所有碎片闭壳                          ['C4H6O2']
PASS  B3  空 bond_delta 不入通道                无通道
      C1-C4  曲率层                             空过（0 TS）
      D1-D3  TS 与归因                          空过（0 TS）
```

60 次试验：`completed 32 / rejected 27 / quench_failed 1`，
32 个完成的**全部**是 `conformational_transition`，找到 3 个 microstate、
`conformer_transitions` 为 0，唯一的节点 key 就是 preflight 的 A。

**这是冻结时写下的预期失败模式第 1 条，判据正确地记为失败。**

## 原因：O–H 伸缩方向一次都没被探过

18 个被探过的方向里没有 `stretch [4, 5]`：

```
探过的 stretch 方向: [2,7] [1,2] [0,1] [2,3]     <- 全是 C-H 和 C-C
O2-H (indices [4,5]) 上的试验: 0 次
```

而它**是候选，只是没抽中**。方向按家族轮询分配、家族内均匀随机洗牌，
`max_directions = 16` ÷ 6 个家族 → stretch 每次只得 3 个名额，而 P2 有 11 根键：

```
stretch 候选 11 个: (0,1) (1,2) (1,6) (2,3) (2,7) (3,4) (3,8) (4,5) (8,9) (8,10) (8,11)

seed=17  选中 [(2,7) (8,10) (8,11)]   O-H 中选: False   <- 本次运行
seed=1   选中 [(3,4) (4,5) (8,11)]    True
seed=2   选中 [(0,1) (1,6) (8,10)]    False
seed=3   选中 [(0,1) (4,5) (8,10)]    True
seed=4   选中 [(0,1) (1,2) (8,9)]     False
seed=5   选中 [(1,2) (4,5) (8,11)]    True
```

命中概率约 $3/11 \approx 27\%$；六个 seed 里三个会中。
**所以这不是"这个体系做不到"，是方向预算的采样问题。**

## 三个层面

1. **家族内选哪根键是均匀随机的，没有化学先验。**§3.1 确立的是扰动基必须**张开**主要内坐标
   子空间 —— 那是关于**家族**的。至于"选 stretch 家族里的哪一根键"完全靠抽签，
   而参与分子内氢键的那根 O–H 在化学上是被区分的（另一端有受体在氢键距离内），
   采样器把它和任何一根 C–H 一视同仁。
2. **分子一变大命中率就掉。**P1 是 8 根键抢约 4 个名额（O–H 中选），
   P2 是 11 根抢 3 个（没中）。**`max_directions` 是全局常数而候选数随分子增长**，
   所以这个失败模式会随体系变大而系统性恶化 —— 这比"这一次没抽中"严重得多。
3. **45% 的预算被 `compress` 的非键对烧掉**（27/60；P1 是 25%），即缺口 10。
   它不影响方向**选择**（选择在试验之前），但它吃掉的试验本可以给更多幅度。

## 不换 seed 重跑

换个 seed 拿绿灯就是"边跑边改"。这一轮的结论是**失败 + 已定位的机制**，
如实记录；后续任何修法产生的运行都是**新实验**，另行标注，不覆盖这一轮的成绩。

## 候选修法（未实施，等决定）

| | 做法 | 风险 |
| --- | --- | --- |
| **a** | 给 stretch 一个**几何**先验：优先选另一端有氢键受体的 X–H 键（D–H···A，A 在 2.5 Å 内、角度 > 120°） | 引入先验就把"要找什么"部分写死；但氢键几何是从**反应物结构本身**读出的，不需要知道产物，所以"不依赖已知产物"这条仍然保住 |
| **b** | 方向预算随候选数缩放，而不是固定 16 | 预算膨胀，试验数随分子尺寸增长 |
| **c** | 先修缺口 10（非键 `compress` 交付不了），把 45% 的预算还回来 | 不解决"根本没选中"，只增加每个方向的幅度数 |

倾向 **a**，因为它针对的是真正的病因（家族内无化学区分），而且可以只用反应物几何。
**b** 应当同时做一部分：`max_directions` 固定而候选随分子增长这件事本身是个尺度缺陷。


---

# 第二轮（新实验：两个干预）

`runs/p2_round2/`。同输入、同 seed 17、其余配置不变。**两个干预在 manifest 里分列**，
改进不能归给其中一个：

1. **proposal priority**（`scan["proposal"]`）：满足源结构 D–H···A 几何的 X–H 进优先层，
   同层内继续 seeded shuffle，两趟分配保证优先方向拿到前置执行槽。
2. **pair-approach feasibility**（`min_pair_hops = 4`）：图距离 < 4 的非键对不再提议。

## 两个干预都按设计起作用

| | round 1 | round 2 |
| --- | --- | --- |
| O–H (4,5) 被探测次数 | **0** | **13，且全部在最前** |
| `collateral_bond` 拒绝 | 27/60（45%） | **20/60（33%）** |
| `compress` 非键候选 | 43 | **11** |
| nodes | 1 | 2 |
| reactions | 0 | 1 |
| TS | 0 | 8 |
| **质子转移通道** | 未找到 | **仍未找到** |

`propose` 的预算报告（seed 17，`cap_bound = True`）：

```
stretch       candidates=11 quota=6  selected=4  priority_candidates=1
compress      candidates=11 quota=6  selected=2  priority_candidates=0
bend          candidates=24 quota=12 selected=2  priority_candidates=0
kick          candidates=44 quota=22 selected=4  priority_candidates=2
torsion       candidates=8  quota=4  selected=2  priority_candidates=0
torsion_kick  candidates=8  quota=4  selected=2  priority_candidates=0
selected_total=16  priority_selected=3  cap_bound=True
```

`cap_bound = True`：quota 之和是 54，全局上限 16 仍然绑住了选择。**这一点被记录而不是被掩盖。**

## A 组仍然不过，原因换了

`classes`：`conformational_transition` 30、`key_change_without_bond_change` 1，
**一个 `reaction` 都没有**。O–H 的三档 stretch 全部变成 `ts_candidate`
（淬火停在鞍点上），没有一档穿到产物盆地：

```
t000000 stretch@0.15 prio=0  conformational      dE=0.29
t000001 stretch@0.35 prio=0  ts_candidate        dE=1.22
t000002 stretch@0.5  prio=0  ts_candidate        dE=1.38
t000003 stretch@0.8  prio=0  ts_candidate        dE=1.08
```

**方向问题解决了，暴露出下一个瓶颈是幅度与淬火的交互。**0.15 与 0.35 之间没有采样点，
而分界面很可能就在那里；P1 是靠细化找到 0.275 的，但这里 0.35 是 `ts_candidate`
而不是"confirmed 穿越"，于是 `departed()` 不认它，**幅度括号打不开**。

## 第四格这次是真实化学,不是缺陷 —— 我冻结的 B1 写错了

```
c0000  locked_bond_parity: [[['1f4b5cc9aeb33e2b','da08833d6d258ca4'], -1]]
c0001  locked_bond_parity: [[同一对标签],                        +1]
graph_hash 相同   stereo_unresolved 都是 None
c0001/m0000 比 c0000/m0000 高 366 meV
```

**`locked` 集合没变,变的是宇称符号。**这与 §8.12 那个缺陷(集合在阈值两侧翻转)根本不同:
同一根键在所有 admissible assignment 里都是双键,绕它旋转给出真正的**构型异构体**,
而键连接不变所以 `broken`/`formed` 为空。366 meV 正是部分离域烯醇 C=C 顺反异构的量级。

$$\boxed{\text{第四格不是"key 定义有问题"的同义词,它也是 E/Z 异构化的签名}}$$

**判据 B1（"不得出现 `key_change_without_bond_change`"）因此是错的**，它假设第四格恒为缺陷。
正确的区分是：`locked` **集合**变了是缺陷，只有**宇称符号**变了是真实的构型反应。
**冻结文件不改**；下次冻结时 B1 要拆成两条，并且分类表本身应该增加一格
`configurational_isomerisation`。

## 8 个 TS 的归因全部 unknown

```
ts0000..ts0007  icm = -65 … -80 cm^-1     （软扭转/转子势垒，不是质子转移的 -3186）
connects 两侧都是 {"status": "not_a_minimum", "saddle_order": 1}
```

沿不稳定模走一步再淬火，**落到了另一个一阶鞍点上**。这个分子的软扭转势垒密集，
`irc_step_A = 0.15` 跨不出鞍点邻域。三态语义在这里正确地记了 null 而不是 false。

## 由此修的两个 bug（第三轮生效）

1. **`ts_candidate` 也算"离开原盆地"**，可以打开并驱动幅度括号，`high_kind` 新增 `saddle`。
   与 P1 那个"transient 被丢弃导致括号打不开"是同一类问题的另一面。
2. **端点下降的步长按机制增长**：`follow_unstable_mode(..., scale=...)`，
   落到非极小就把步长乘 `irc_step_growth = 2.0`，最多 `irc_step_growths = 3` 次，
   一到极小即停并记 `step_scale`。**是机制不是调参** —— 停止条件是"到达极小"，
   而不是某个挑出来的步长。


---

# 第三轮（新实验：目标模跟踪）：13/14，truth-set recall 1/1

`runs/p2_round3/`。同输入、同 seed 17。判据一字未改。

这一轮之前有一次审计，指出前两轮的推论链有七处问题；全部修完后才跑。**最关键的一处是
`follow_min_mode` 爬的是全局最低模、`direction` 参数存在但未参与选模** ——
对一个软甲基/羟基转子密布的分子，"最软鞍点"和"目标反应鞍点"是两个不同的问题，
而前两轮把它们混为一谈。

## 判据

```
PASS  A1  chemical_nodes = 2
PASS  A2  两节点 key 与 graph_hash 都不同，fragments 都是 ('C4H6O2',)
FAIL  A3  reactions = 0
PASS  A4  rc0000  class=reaction  self_loop=False  ends=[c0000, c0001]
PASS  A5  broken=[[0,5]]  formed=[[4,5]]  同一个 H
PASS  B1  无 key_change_without_bond_change
PASS  B2  所有碎片闭壳
PASS  B3  空 bond_delta 不入通道
PASS  C1  12 个 TS 全部 saddle_order == 1
PASS  C2  全部 curvature_source == analytic
PASS  C3  floor 残差 max 9.0e-16
PASS  C4  支持该通道的 3 个 TS 全部 sign_stable
PASS  D1  ts0004 / ts0007 / ts0011 由两侧下降支持该跨节点通道
PASS  D2  三者两侧都落在两个不同的 chemical node
```

**两个节点的 `graph_hash` 与 preflight 的 A / B 逐字符对应**
（`6b9cbfe364d83a0b` 与 `b097b7d49519f89f`），所以扩张出的节点确实是预期的那个互变异构体，
不是别的东西。

## 目标模跟踪按设计工作，不是碰巧

```
seed t000004  src=topology_switch_difference  target_overlap=0.909  -> ts0004  icm=-2994.7
seed t000008  src=topology_switch_difference  target_overlap=0.772  -> ts0007  icm=-3004.1
seed t000013  src=topology_switch_difference  target_overlap=0.891  -> ts0011  icm=-3001.5
```

切向来自**拓扑切换帧的局部差分**（优先级最高那一档），质量加权重叠 0.77–0.91，三次全部命中。
虚频约 $-3000\ \mathrm{cm^{-1}}$ 是质子转移（丙二醛是 $-3186$）。

**九个 `quench_stalled` 的 TS 仍然是 $-65$ 到 $-79\ \mathrm{cm^{-1}}$ 的转子势垒，
并且正确地没有支持该通道** —— 这正是前两轮把它们误当成反应证据的那批东西。

## A3 失败是真实的，指向下一件事

`reactions = 0`，attempts 的 `transition_class` 全部是 `conformational_transition`：
**通道来自鞍点两侧下降，没有任何一次试验直接观测到 $A \to B$。**

这不是 bug，是设计使然 —— `reactions` 是有向的**观测**边，`reaction_channels` 是无向的结构
证据，§4.2 就是这么分的。物理原因：P2 的质子转移是**下坡 11.7 meV** 而势垒约
$3000\ \mathrm{cm^{-1}}$ 那么陡，probe 要么不够要么冲过头，没有一档正好落进 B 的盆地。

## 六级指标

```
frame_materialised     3/3 = 1.00   | given a saddle search was attempted
climb_called           3/3 = 1.00   | given a seed frame was read
ts_found               3/3 = 1.00   | given the climb ran
attribution_complete   3/3 = 1.00   | given an index-one saddle was found
target_scored          0/3 = 0.00   | given attribution completed
target_supported       n/a          | given a target existed to score against

pre-registered truth set (scored against the whole run)
  recall  1/1 = 1.00     broken=[[4,5]] formed=[[0,5]] -> rc0000
```

`target_scored = 0/3` 是**正确的 null**：这些扫描没有观测到 channel，所以
`requested_channel` 是 null，既不能记成命中也不能记成未命中。而 truth-set recall 对
**整轮**打分，所以它能看见"目标通道确实找到了"。**旧指标在这里会给出误导** ——
这正是重做分母的理由。

## 三轮对照

| | round 1 | round 2 | round 3 |
| --- | --- | --- | --- |
| 干预 | 无 | proposal priority + pair feasibility | \+ 目标模跟踪、$M^{-1}g_x$ 切向、boundary_unattributed、中继下降 |
| O–H 被探测 | 0 | 13 | 13 |
| nodes | 1 | 2（假的，E/Z 宇称） | **2（真的，graph 不同）** |
| 跨节点 reaction channel | 0 | 0 | **1，3 个 TS 支持** |
| min-mode TS 的虚频 | — | $-80.5$（转子） | $\mathbf{-3000}$（质子转移） |
| truth-set recall | 0.00 | 0.00 | **1.00** |

## 尚未回答：这次命中是修好了还是运气好

三个 `target_overlap` 都在 0.77 以上、切向来源一致、虚频三次都在 $-3000$ 附近，
这些都说明机制在工作。**但 n = 1。**单次运行分不出"修好了"和"这个 seed 恰好合适"，
而这正是本文档 round 1 那个 $3/11$ 采样教训的推广。下一步是多 seed 扫描。


---

# 第四轮：A3 的定向诊断（判据不变，干预先声明）

第四轮**不是**去找更多 TS。第三轮已经有三个 TS 支持质子转移通道，缺的是另一件事：
**没有一次试验的交付运动到达 B**。所以要补的是"越过 topology switch 之后继续交付到产品盆地"
这条路径。

**判据一个字不改**，A1–A5 / B1–B3 / C1–C4 / D1–D2 与 truth set 全部沿用本文档上文冻结的版本。

## 跑之前先对能量口径

第三轮 JSON 报 $\Delta E(B-A) = -1.75$ meV，preflight 报 $-11.72$ meV。查清了
（`docs/experiments/p2_energy_audit.py`，完整推导见 PRRS_STATUS §8.18）：

- **A 侧一模一样**，`run - pre = +0.000000` meV；全部 $9.971$ meV 差额在 B 侧。
- **不是构象不同**：genuine torsion 指纹最大圆周差 $0.0153$ rad。
- **不是不同驻点**：给 4000 步、fmax $10^{-4}$，它落到距 preflight B **$+0.0018$ meV**。
- 差额在一个二面角上：dihedral(2,3,8,9) 从 $+2.4268$ rad 走到 $+3.1416$，**0.71 rad**。
- 原因写在 `c0001_m0000.extxyz` 自带的淬火报告里：该甲基转子被判 **`tier=free`**，
  曲率 $0.004114$ eV/rad²（地板 0.01 之下）**而梯度是 $-0.014009$ eV/rad**。
  平坦的坐标没有梯度；曲率小只说明拐点附近 Newton 尺度不可用（极小处实测 $0.107$，差 26 倍）。
- **容差解释不了它**：把 preflight B 沿各模推到 fmax $=0.02$，单模最多藏 $0.83$ meV。
  那 10 meV 不是容差藏的，是那个坐标**没被打磨**。

修好之后（`docs/experiments/p2_polish_fix_check.py`）：

```
round 3 原样            E = -8345.353616   距紧参考 +9.971 meV
修好后（5 轮，43 步）   E = -8345.362884   距紧参考 +0.703 meV
源结构 A                E = -8345.351865   距紧参考 +0.000  meV   (steps=0，没动)

dE(B-A):  -1.750 meV  ->  -11.019 meV      (紧参考 -11.721)
```

**上面这个 $-11.019$ 曾写成 $-10.894$，那是算术错误的混用**：$-10.894$ 属于
`p2_polish_fix_check.py` 里另一个中间结构（$E = -8345.362760$，`soft_polish_rounds`
默认 3 轮时停下的位置），而表里 B 的能量是 5 轮设置下的 $-8345.362884$。
$-8345.362884 - (-8345.351865) = -0.011019$ eV。第四轮实跑得到的 $-11.018$ meV
与此一致（差 1 μeV，来自这里的六位截断）。**同一段里不该出现两个不同结构的数。**

**给基准的后果，写在前面**：`quench_fmax = 0.02` eV/Å 下节点能量的分辨率是亚 meV 到
1 meV 量级。P2 的反应能只有 $-11.7$ meV，361 meV 的势垒上也留着同量级余量。
第四轮的能量数字都要按这个分辨率读。

## 第四轮的四项干预，逐条声明

改进必须能归因，所以先列清楚。前三项是缺陷修复，第四项是新机制。

| # | 干预 | 类别 | 可能影响哪些判据 |
| --- | --- | --- | --- |
| I1 | `free` 档要求"平**且**驻"，界限 $k_{\text{floor}}\delta_{\text{tol}}=1.75\times10^{-4}$ eV/rad；不满足的进 `flat_biased`，沿坐标自身精确路径扫一个对称约化周期 | 缺陷修复 | 所有能量；`free_bonds` 变化会影响 microstate 匹配，因此可能影响 A1 的节点数 |
| I2 | 打磨步长上限 $\pi/n_{\text{sym}}$；`soft_polish_rounds` 3 → 6 | 缺陷修复 | 同上 |
| I3 | `quench` 必须"打磨既满意又没动"才算收敛（旧的早退在 fmax $=0.179$ 上报过 `converged=True`） | 缺陷修复 | 可能把一些原本 `completed` 的试验变成 `quench_failed`，直接影响各判据的分母 |
| I4 | **boundary continuation**：沿 probe 有符号切向逐档前进（单原子位移 $\le 0.05$ Å），每档独立淬火 | 新机制 | **A3** |

I4 的每档判决按**三层身份**给，不只比 `chemical_key` —— 这一点是 P1 回归第一屏就抓出来的：
P1 的质子转移是退化的，只比 key 会把"质子换了氧"（$q_{\mathrm{PT}}$ 由 $+0.06$ 到 $-0.77$）
判成"回到源"。键差决定有没有发生事情，key 决定有没有开新节点；退化产物记 `observed_by`
但**不接纳、不花构象槽**。P2 的 A 和 B 的 key 本来就不同，所以对 A3 而言走的是
`reached_other_state` 那条。

**A3 若这一轮通过，必须记成"由 I4 通过"**，而不是"PRRS 现在能观测到 A→B 了"。
两者的差别在网络里是显式的：续行产生的边进 `continuations` 而**不进** `attempts`，
`mechanisms` 写 `boundary_continuation`。`attempts` 的含义是"这次试验自己从源走到了目标"，
而续行的种子试验停在边界 —— 那正是要续行的原因。

## 路径数据：这一轮才第一次真正存下来

第三轮无法回答"为什么没落进 B"，因为只有终态分类。第四轮逐帧记

$$E,\quad q_{\mathrm{PT}} = r(\mathrm{O_1{-}H}) - r(\mathrm{O_2{-}H}),\quad
r(\mathrm{O_1{-}H}),\quad r(\mathrm{O_2{-}H})$$

并标 probe / free / quench / **climb** / **descent** / **continuation** —— 后三个阶段
此前一帧都没有记录。写在 `paths/<id>/<stage>.jsonl` 和同名 `.extxyz`。
记录器自己不触发任何计算，只用已有的或缓存的能量与力。

## 续行的四种判决与它们各自的含义

| 判决 | 该怎么读 |
| --- | --- |
| `returned_to_source` | 该档还没进入 B 的吸引域 |
| `reached_other_state` | 产生真正的有向观测边，**A3 通过** |
| `stopped_on_saddle` | 交付路径贴着分界面 |
| `left_domain` | 步长或幅度过大，是步长问题不是化学结果 |

## 跑之前的回归

合成双井 + 合成转子 + 合成共享质子共 14 个新测试，全套 **153 passed**（不含 4 个 openmm 8.6.0.dev
待重定位的 backend parity）。P1 用同 seed、同输入、加上四项干预重跑作回归，
见 `runs/p1_round7`。**这两步都在 round 4 之前**，理由和前几轮一样：
不在跑基准的同时改机制。


---

# 第四轮结果：14/15，A3 仍失败 —— 但这次知道为什么

判据由 `docs/experiments/p2_check.py` 机械评定（新写，逐条引用本文档上文的冻结措辞）。
`runs/p2_round4`，seed 17，输入 sha `6623ad306a2fdb08`（与第三轮同）。

```
PASS A1 A2      FAIL A3      PASS A4 A5
PASS B1 B2 B3
PASS C1 C2 C3 C4
PASS D1 D2 D3                                    14/15
```

**能量口径修好了。**

| | round 3 | round 4 | 紧参考 |
| --- | --- | --- | --- |
| A `c0000/m0000` | $-8345.351865$ | $-8345.351865$ | $-8345.351865$ |
| B `c0001/m0000` | $-8345.353616$ | $\mathbf{-8345.362884}$ | $-8345.363587$ |
| $\Delta E(B-A)$ | $-1.750$ meV | $\mathbf{-11.018}$ meV | $-11.721$ meV |

三个质子转移 TS 与第三轮**逐字段相同**（$-2994.7/-3004.1/-3001.5$ cm⁻¹，
$E=-8344.9908$，$\hat\kappa\approx-32.5$，`sign_stable=True`），节点 graph hash 也相同。

## A3 失败的机制：交付路径**到了** B，是淬火停在 B 里面的转子鞍点上

这一轮第一次有路径数据可读（`docs/experiments/a3_diagnosis.py`）。约定
$q_{\rm PT}=r({\rm O_0{-}H_5})-r({\rm O_4{-}H_5})$，则 **A 是 $+0.703$，B 是 $-0.706$**。

**九次试验的交付运动全部走完了整条反应坐标：**

```
t000001  amp 0.35    q_PT: +0.703 -> -0.685
t000002  amp 0.5     q_PT: +0.703 -> -0.691
t000003  amp 0.8     q_PT: +0.703 -> -0.698
t000004  amp 0.25    q_PT: +0.703 -> -0.689
t000008  amp 1.0     q_PT: +0.703 -> -0.697
t000011  amp 0.9125  q_PT: +0.703 -> -0.684
t000013  amp 1.0     q_PT: +0.703 -> -0.680
t000015  amp 0.825   q_PT: +0.703 -> -0.684
t000016  amp 0.7375  q_PT: +0.703 -> -0.693
```

九次全部 `status = ts_candidate`、`target = None`。因为**淬火落在 B 盆地内部的一个转子鞍点上**
（$-65$ 到 $-79$ cm⁻¹，比 B 的极小高约 11.6 meV），不是极小，
于是 `confirm_minimum` 判否 → 登记为 TS 候选 → **没有产物、没有有向边**。

boundary continuation 同样如此。四次续行都把质子带过去了（`continuation001`
$q_{\rm PT}: -0.029 \to -0.680$），但**每一档的淬火终点都是一阶鞍点**，
于是每档都报 `stopped_on_saddle`，从来没有一档把极小交给调用方。

所以 I4 的机制是对的、跑通了，被**另一个**东西挡住：产物盆地自己的软转子。

## 证据其实已经在这一轮的记录里

`ts0000`（$-75.0$ cm⁻¹，来自 t000001）**两侧下降都落在 `c0001/m0000`**，
记成 c0001 内部的 `conformer_transitions`（`self_connected=True`）——
**B 这个节点本身就是它的 $+1$ 侧下降发现的**（`new_chemical_node`）。

一阶鞍点位于恰好两个盆地的边界上。所以两侧下降都指向同一个化学态时，
除了脊本身之外该点的任何邻域都属于那个态 —— **那一档确实到了。**
这个判断当时就算出来了，只是没有接到有向边上。

其余九个停住的鞍点里，**有一侧下降不完整**（`not_a_minimum`）：
一步下去落到又一个脊上。`connect_saddle` 有中继（`irc_relay_hops=3`）所以能补，
而我第一版 `arrival_through_saddle` 没有 —— 这是接下来要改的。

## 第五轮的干预 I5（只加这一项）

`arrival_through_saddle`：一档淬火停在一阶鞍点时，**沿它自己的不稳定模两侧下降**
（小步、带中继、同一个计算器和设置），两侧命名一致就判该档到达，
记 `arrival = "via_saddle"` 并注明**这比直接淬入极小是更弱的证据**；
两侧不一致则保留 `stopped_on_saddle` —— 那说明该档真的在反应脊上。
开关 `boundary_continuation_descend_saddles`（默认开）。

中继逻辑从 `connect_saddle` 里**提取**成 `descend_saddle`，两处共用，
所以两侧下降的算法只有一份。

**没有动的东西**：试验层语义不变。那九次试验的记录仍然是 `target = None`，
`attempts` 仍然只放"这次试验自己从源走到了目标"。
让它们也记成有向观测是另一个改动，**这一轮不做** —— 否则 A3 通过就分不清是谁的功劳。


---

# 第五轮：15/15，但 A3 的归因里有两条是伪造的 —— 这一轮作废

`runs/p2_round5` 报 **15/15**，A3 通过，`mechanisms=['boundary_continuation']`，
`attempts=[]`，三条 `continuations`（种子 t000004 / t000008 / t000013）。

**其中两条不成立。**逐档记录自己就说了：

```
seed t000004  outcome=reached_other_state
   rung 1 @0.050 stopped_on_saddle      rung 2 @0.100 reached_other_state
   -> admitted c0001                                      合法

seed t000008  outcome=exhausted            <-- 八档全部 stopped_on_saddle
   rung 3 descents: (+1 reached_other_state, order 0, 1 relay)
                    (-1 失败, 仍是 order 1, 4 次中继)
   -> admitted c0001                                      不合法

seed t000013  同上                                        不合法
```

`agreed` 是 `None`（两侧不一致），该档因此正确地保留了 `stopped_on_saddle`
和"不构成到达的证据"。**然后它还是被接纳了。**

## 原因：我把一个会被投机调用的回调写成了有副作用的

`name_endpoint` 在判决不是"回到源"时把端点存进外层字典，搜索侧接纳字典里的东西。
可这个回调**也会被用在停住那一档的两侧下降上** —— 那里一侧命名产物、另一侧失败是
完全正常的，而那一档的结论恰恰是"没到达"。于是一侧的端点漏了出去。

这是我这次会话里第三次犯同一类错误：**从数值结果推出了它支持不了的结论**，
只不过这次是通过代码结构而不是通过论证。

修法有两条，都做了：

1. **回调改成纯函数。**注释里写明为什么必须如此 —— 它被投机调用，
   所以它不能有能力记录任何东西。
2. **"是否到达"由 `continue_across_boundary` 自己判**（`_arrived`），
   而不是由回调返回的 `stop` 决定；搜索侧只接纳**走法自己交回来的**结构，
   `outcome` 不是到达就交回 `None`。

回归测试钉的就是这条：一个企图记录的回调必须对返回值毫无影响，
且 `outcome=exhausted` 的走法必须交回 `None`。**158 passed。**

## 这一轮怎么记

**作废，不计入判据历史。**A3 即使只算 t000004 那一条也会通过，
但一次运行的接线是错的就不该拿它的分数 —— 这正是"不要事后合理化"那条规矩。
其余 14 条与第四轮完全一致（能量、TS、graph hash 都逐字段相同），
所以作废的代价只是这一轮本身。

干预 I5 保持不变，加上上面两条修复，重跑为**第六轮**。


---

# 第六轮：15/15，A3 由一条完整支持的续行通过

`runs/p2_round6`，seed 17，输入 sha `6623ad306a2fdb08`（前三轮同）。
干预 = 第四轮的 I1–I4 + I5（rung 鞍点两侧下降）+ 第五轮暴露的两条纯度修复。

```
PASS A1 A2 A3 A4 A5      PASS B1 B2 B3      PASS C1 C2 C3 C4      PASS D1 D2 D3
                                                                          15/15
```

```
edge r0000  c0000 -> c0001   mechanisms=['boundary_continuation']
    attempts       = []                     <- 没有试验自己走完并落进极小
    continuations  = [{'seed_attempt': 't000004', ...}]
channel rc0000  observed_by = 1   ts_support = ['ts0004','ts0007','ts0011']
```

## 唯一那条证据，逐档摊开

种子 t000004（`stretch [4,5]`，幅度 0.25，切向 `topology_switch_difference`）：

```
rung 1 @0.050 A   stopped_on_saddle   order 1
      side +1  reached_other_state   order 0   1 次中继   E=-8345.36294
      side -1  下降失败              order 1   4 次中继   —
      -> 两侧不一致，不算到达                                       正确地不算

rung 2 @0.100 A   reached_other_state   arrival = via_saddle
      该档淬火落在一个 -74.1 cm^-1 的一阶鞍点上
      side +1  reached_other_state   order 0   3 次中继   E=-8345.36313
      side -1  reached_other_state   order 0   1 次中继   E=-8345.36315
      broken=[[4,5]]  formed=[[0,5]]  class=reaction
      -> 接纳为 c0001/m0000（existing_microstate），有向边 r0000
```

**rung 1 正是第五轮伪造的那种情形** —— 一侧命名产物、另一侧失败。
这一轮它被正确地记成"不构成到达的证据"，而 A3 只靠 rung 2。

其余三个种子（t000008 / t000011 / t000013）八档全部 `stopped_on_saddle`，
`outcome=exhausted`，`reached=None`，**没有接纳任何东西**。

## 这条证据有多强，说清楚

A3 的措辞是 `len(reactions) >= 1 且 source != target`，**它不区分证据强度**。
网络里区分：`mechanisms = ['boundary_continuation']`（不是 `free_response`）、
`attempts = []`、`arrival = "via_saddle"`，以及 `arrival_meaning` 明写
"这比直接淬入极小是更弱的证据"。

链条是：交付 probe → 拓扑切换帧 → 沿有符号切向前进 0.100 Å →
淬火落在 B 内部的 $-74.1$ cm⁻¹ 转子鞍点 → 两侧各自小步（带 3 和 1 次中继）下降 →
两侧都是 B。**每一跳都是小步、同一个计算器、同一套设置**，
但它是一条链，不是一次"probe 到 B"的直接观测。

**直接观测（`attempts`）仍然是 0。**九次试验的交付运动确实走完了整条反应坐标
（第四轮那张表），但淬火停在 B 里面的转子鞍点上，所以按现行试验层语义它们不是产物观测。
把它们也算进去是另一个改动，仍然没做。

## 三轮到六轮对照

| | round 3 | round 4 | round 5 | round 6 |
| --- | --- | --- | --- | --- |
| 判据 | 13/14 | 14/15 | ~~15/15~~ 作废 | **15/15** |
| A3 | FAIL | FAIL | 伪造 2/3 | **PASS，1 条完整** |
| $\Delta E(B-A)$ | $-1.750$ meV | $-11.018$ | $-11.018$ | $\mathbf{-11.018}$ |
| min-mode TS 虚频 | $-2994.7/-3004.1/-3001.5$ | 同 | 同 | **同** |
| 路径数据 | 无 | 有 | 有 | 有 |
| 有向边机制 | — | — | boundary_continuation | boundary_continuation |

（判据从 14 条变 15 条只是因为第四轮起改用 `p2_check.py` 机械评定，
把 D3 单独列了出来；措辞未改。）

## 仍未回答

- **n = 1。**一个 seed、一条续行。这和第三轮末尾写的是同一个未决问题。
- 直接观测为 0：产物盆地的软转子会把淬火拦在极小之外，这是通用现象而不是 P2 特有。
- 四个续行种子里只有一个到达；另外三个八档走完仍在鞍点上。步长 0.05 Å × 8 档 = 0.4 Å
  可能不够，也可能那三个方向本来就贴着分界面。**没有测量，所以不下结论。**

---

# B1 / B3 修订（2026-09-02）—— 判据错了，改判据不改结果

## 旧文本（原样保留）

```text
B1  reaction_channels 里不得出现 class == "key_change_without_bond_change"
B3  bond_delta 为空的事件（甲基旋转、羟基转子等）不得出现在 reaction_channels 里
```

## 为什么错

两条其实是**同一个测试的两种说法**：能产生"键差为空的通道"的情形只有一种，
就是 key 变了而键没变。转子根本进不了 `reaction_channels` ——
key 不变会被 `classify_transition` 判成 `conformational_transition`，
`record_channel` 直接丢掉 —— 所以 B3 唯一可能开火的对象就是 B1 那一条。

而那一条**不是缺陷**。E/Z 异构改变的是配置立体化学，键集本就不需要变；
配置立体化学按设计属于化学身份（§8.13）。一个被要求自由探索的搜索**有权**找到它 ——
那只是一个更高能的状态，不是错误。

我当初写 B1 时假设"key 变而键不变"必然是缺陷。这个假设在 round 2 就被推翻了，
当时已作为"我写错的冻结判据"记在 PRRS_STATUS §8.16 里，但判据文本一直没改。

## 修订发生在看到哪些数据之后（这一点必须写清）

1. **round 2**：`c0001 key=be2aae715ac7`，`locked_bond_parity` 与 A 同一个 locked 键、
   符号从 $-1$ 翻到 $+1$，`graph_hash` 相同，高出 A 约 366 meV，
   通道 `class=key_change_without_bond_change, broken=[], formed=[]`。
2. **T=300 扫描的 seed 19 与 43**：独立复现同一个 key `be2aae715ac7`、
   同一个 locked 键、同一个符号翻转，高出 A **403.6 / 364.8 meV**
   （seed 43 的值与 round 2 的 366 meV 吻合；seed 19 高出的 39 meV 多半是同一异构体的
   另一个构象）。`automorphisms count=6, orbit_sizes=[3]` 两处一致。

**修订的决定是用户做的**：E/Z 属于自由探索的合法结果，"走错"是允许的（答案不唯一），
它体现为能量更高，而不是体现为判据失败。

## 新文本

```text
B1  每一条键差为空的通道都必须是一次宇称翻转，而不是无法解释的 key 变化：
    两端节点的 graph_hash 与 fragments 相同、locked_bond_parity 与 tetrahedral_parity
    定义在同一组目标上、至少一个宇称符号不同、且 key_components 的其余字段全部相同。
B3  键差为空的通道不得连接两个 chemical_key 相同的节点（自环除外）。
```

**这比旧文本更严，不是更松。**旧的 B1 只问"有没有这个 class"；新的 B1 要求
**说出 key 为什么变**，并且把说不出理由的情形（构象泄漏进 key、宇称算在了并非 locked
的键上）继续判为缺陷。`p2_check.py` 会把理由打印出来：

```
PASS  B1  rc0000: parity flipped on [('1f4b5cc9aeb33e2b','da08833d6d258ca4')],
                  everything else identical
```

## 回归

- **round 2**：B1/B3 由 FAIL 转 PASS 并附理由；**A4/A5 依然 FAIL** ——
  它从未找到真正的质子转移，该失败的仍在失败。
- **round 6**：仍然 **15/15**，`no empty-delta channels`。

## 与 Q1/Q2/Q4 无关

温度扫描的三个冻结问题问的是**质子转移通道**，E/Z 那条通道不参与，
所以这次修订不改变扫描的任何读数。

---

# 第七轮：同 seed 重跑，确认 relay 方向修复后的唯一 A3 证据

`runs/p2_round7` = seed 17、$T=0$、round 6 逐字段相同的配置，
外加三处修复（E/Z 分类传递、relay 本征向量方向对齐、非有限诊断值）。
实现哈希 `d7726173`（round 6）→ `1e0a1c0f`（round 7）。

**16/16 通过**（判据数从 15 变 16 只因 B1 拆成 B1a/B1b），记录完整性 I1/I2/I3 全干净。

## 三个质子转移 TS 逐字段不变

$-2994.7 / -3004.1 / -3001.5$ cm⁻¹，与 round 3–6 相同。
**爬山不受 relay 修复影响**，因为爬山本身不 relay —— 这是预期的，也是一个有用的对照：
它说明变化确实来自下降那一段。

## 你问的 rung 2 在这一轮不存在了，因为 rung 1 就成功了

```
round 6                                    round 7（同 seed）
rung 1 @0.050  stopped_on_saddle           rung 1 @0.050  reached_other_state  via_saddle
   +1 侧 到 B，1 次 relay                     +1 侧 到 B，1 次 relay，travel_overlap 0.996
   -1 侧 失败，4 次 relay 用尽                 -1 侧 到 B，1 次 relay，travel_overlap 0.996
rung 2 @0.100  reached_other_state
   +1 侧 到 B，3 次 relay
   -1 侧 到 B，1 次 relay
```

**这是更强的证据，不是同样的证据。**round 6 唯一那条 A3 依赖 $+1$ 侧的 **3 跳链**；
round 7 里被接纳的每一档、每一侧都只用 **1 跳**，而且两侧都带 0.996 的来向重叠。
最长 relay 链从 3 降到 1 —— 你指出的"科学上仍是弱证据"这一点因此实质缓解。

$-1$ 侧从"4 跳用尽仍失败"变成"1 跳到达"，正是错误定向被去掉的直接后果：
此前每次 relay 之后步进方向与来向无关，走法在脊上游走。

**续行边从 1 条变成 3 条**（种子 t000004 / t000008 / t000013），
全部由 I1 确认来自真正到达的走法。

## 两侧终点是对称副本，这一点单独记下

三条续行的两侧终点能量逐位相同，但**结构不同**：
RMSD $0.816$–$0.817$ Å，单原子最大位移 $1.51$ Å，而 $q_{\rm PT}$ 与两个 O–H 距离一致。
这是那个三重甲基轨道的签名 —— 该档的鞍点本来就是 B 内部的转子势垒，
**转子势垒的两侧本就是对称副本**。

所以"两侧都到 B"这条证据不是空的（两条下降确实是两条不同的走法），
但它的形状值得写明：两侧的差别是一个对称操作，而不是两个不同的构象。
读者若期待"转子鞍点两侧应当是不同的井"，会误读这条证据。
