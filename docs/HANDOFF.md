# PRRS — 交接文档

**最后更新** 2026-09-05（三项正确性修复之后；接手请先读 §12 与 `DECISIONS_PENDING.md`）· 环境 `openmm_dev` · 工作目录 `/home/ruigengji/MLPATH`
（不是 git 仓库；`implementation_sha256` 是对 `src/prrs/*.py` 的哈希，是唯一的代码指纹）

**当前实现 `5b69dfed41c4bdb2`**。别手抄这个数 —— 用 `make hash`
（`scripts/implementation_hash.py`）现算，它和 manifest 用的是同一个算法。
这一行以前是手抄的，抄完就过期了：它写着 `9c13fdb5` 时树已经是 `dba26be3`，中间隔了两轮修复。

沿革（每一步都是一次 `src/prrs/*.py` 的改动）：

| 哈希 | 日期 | 内容 |
| --- | --- | --- |
| `0058158d` | 2026-09-02 | **P2 冻结点**，`runs/p2_round7` 与 16 个 `p2_T300_seed*` 都在它上面 |
| `158f2f3b` | 2026-09-02 | P3 全部探针记录在它下面产生 |
| `9c13fdb5` | 2026-09-03 | 四项确定性缺陷修复 |
| `2f02fef9` | 2026-09-03 | `internal_mode_recheck` 的 fail-open 修复（`runs/p3_descent_recheck`） |
| `dba26be3` | 2026-09-03 | `quench` 两处收敛缺陷 + `follow_unstable_mode` 软模接纳（`runs/p3_descent_recheck_fixed`） |
| `23179dea` | 2026-09-04 | **工程化重构，零行为改动**，见 §10 |
| `cf29f66b` | 2026-09-04 | `_stamped_with_charge`：声明的总电荷必须到达势能，冲突则拒绝开跑 |
| `686a5c4f` | 2026-09-04 | `side_seed`：修负 seed（只在 `minimum_check="probe"` 下可达，未打到任何冻结结果） |
| `5d905aaa` | 2026-09-04 | `IdentityUnavailable`：枚举预算溢出不再成为身份；`admit` 的三层分离 |
| `6b372359` | 2026-09-04 | 未知记为未知（`locked_bond_parity` 为 `None` 而非 `[]`）—— **本轮唯一的身份迁移** |
| `69d53501` | 2026-09-04 | 失败原因移出哈希；`withheld_candidates` 检索路径 |
| `5b69dfed` | 2026-09-05 | **三项正确性修复**（见 §12）：验证不完整不再当通过 · 候选出口统一 · 失败落盘独立于网络快照 |

改 `src/prrs/` 前先跑 `make digest`（= `python docs/experiments/digest_regression.py`）
（332 文件 / 2136 帧的 digest 黄金样本，不匹配即非零退出）。
它把差异**分类**：`graph_invariant`（图变了，P1/P2 冻结在上面）/
`parity_only`（只有立体化学 parity 动了）/ `key_only`（组件没变而 key 变了）。
有意改变 digest 的工作用 `--allow parity_only` 之类放行 ——
**但放行只确认差异的形状，不确认新值是对的**，那需要化学判断。

这份文档给接手的人。它不重复细节，只说**现在在哪、什么已冻结、什么还开着、以及这个项目的规矩**。
细节全部在 §7 列出的文档里，每一条结论都能追到一次运行。

---

## 1. 这是什么

**PRRS = Perturbation-Response Reaction Search.** 一个 MLP-native 的反应发现方法：
给分子一个受控扰动 → 无偏 NVE 自由响应 → 淬火 → 看它落到哪个化学态。
它**不**预设反应坐标，**不**预设产物，**不**从 TS 猜起。

核心设计取向（这些是设计决定，不是实现细节，改它们要重新验证一批结论）：

1. **Hessian 是验证极限，不是主观测量。** 有限幅度响应是 PES 上的尺度滤波器，
   只有 $a\to0$ 那个极限才有谱定理。有限幅度的本征值**不得**称为"频率"、
   **不得**报 cm⁻¹。见 `PRRS_CURVATURE_ROUTES.md` §E。
2. **三层身份**，不可互相折叠：
   - **化学态身份** `chemical_key` —— 图 + 碎片 + 配置立体化学，**构象被商掉**
   - **反应事件身份** `reaction_event_key` —— 键差 $B^-=E_0\setminus E_1$、$B^+=E_1\setminus E_0$，
     在 $\mathrm{Aut}(G_0)\cup\mathrm{Aut}(G_1)$ 下正则化、方向商掉
   - **通道** `reaction_channels` —— 允许自环（退化反应就是自环，不是构象变化）
   这一层被折叠回第一层过一次（P1 的退化转移），后果见 `PRRS_STATUS.md` §8.19.5。
3. **失败必须关闭（fail closed）。** 判不出来就拒绝并说明理由，绝不用默认值替代。
4. **有向观测与结构证据分开记。** `attempts` = 这次试验自己从源走到了目标；
   `continuations` = probe 被继续交付后到达；`ts_support` = 鞍点两侧下降（**无向**）。
   三者强度不同，网络里分三个字段。

---

## 2. 项目的规矩（这些是用户定的，不是建议）

**接手之前先读这一节。**违反它们产出的结果会被作废，本项目已经作废过一轮。

- **判据跑前冻结，逐字不改。** 要改必须作为**显式修订**记录：旧文本原样保留、
  改的理由、以及**修订发生在看到哪些数据之后**。已有两次修订（B1/B3、A4），
  格式照抄 `P2_OXOBUTANAL_ACCEPTANCE.md` 和 `P2_TEMPERATURE_SCAN_PROTOCOL.md`。
- **不许"边跑边改"。** 绝不为了得到绿色结果而换 seed 或调参。
- **干预必须分别命名**，这样改进能归因到具体某一项。
- **未完成 ≠ 失败。** 分母是**发起数**；未完成记未完成（0 分因为没回答），
  崩溃记崩溃并附原因。`scripts/scan_status.sh` 机械地区分这三态。
- **测量优先于推断。** 本项目有大量"我推断 X，检验杀了它"的记录，
  那些记录是资产不是耻辱 —— 见 §6。
- **旋钮先验证再设计实验。** 在确认被改的变量能推动被测机制之前设计方差实验，
  等于用算力买一个已经注定的答案。这条是 $T=0$ 扫描 8/8 那次教训的产物。
- **门（gate）不能因为它挡住了想要的结果就调。** 要动必须先量清楚它挡掉的是什么。

---

## 3. 已冻结的结论

### P0 乙醇 — 逐字段复现
频率经独立的 `ase.vibrations` 交叉验证（`src/prrs/validation.py`，**搜索侧从不 import 它**）。

### P1 丙二醛（分子内质子转移，退化）
`runs/p1_round8`，**13/13 判据**，与 `runs/p1_round7` 逐字段相同。

| | |
| --- | --- |
| 通道 | `degenerate_reaction`，自环，`event_key = a04410b76e16b117` |
| TS | 4 个，$-3189.6 / -3188.2 / -3187.3 / -3186.7$ cm⁻¹（参考区间 $-3186$） |
| 势垒 | **419.14 meV** |
| `reactions` | **0 —— 这是对的**：退化转移两端是同一化学态，所以是自环不是有向边 |

> **2026-09-04 追加，冻结文字未动**：上面这个 419.14 meV 第一次有了 QM 参照 ——
> `T[E_QM] = 134.4 meV`（ωB97M-V/def2-TZVPD//r2SCAN-3c）。
> **13/13 没有被推翻**：两张面的鞍点位置几乎重合（对称感知 RMSD **0.0151 Å**），
> 所以那 285 meV 的差是**纯能量的、落在模型上**，不是几何、不是方法。详见 §11 与
> `docs/experiments/p1_qm/README.md`。

### P2 3-氧代丁醛烯醇（不对称，跨节点 A→B）
判据由 `docs/experiments/p2_check.py` 机械评定（16 条 + 3 条记录完整性检查）。

| | |
| --- | --- |
| 最终 $T=0$ 轮 | `runs/p2_round7`，**16/16**，实现 `0058158d` |
| 最终 $T=300$ 批 | `runs/p2_T300_seed{16 个}`，单一实现 `0058158d`，16/16 完成，**0 失败** |
| $Q_1$ 有向边落在转移通道上 | **16/16**，95% Wilson $[0.81,1.00]$ |
| $Q_2$ 转移通道存在 | **16/16**，$[0.81,1.00]$ |
| $Q_4$ 有 min-mode TS 支持它 | **15/16**，$[0.72,0.99]$ |
| 全量评分 | 12 个 seed 16/16；3 个 15/16（仅 B1b）；1 个 12/16（seed 71） |
| 虚频 | $-3010.9$ 到 $-2998.8$ cm⁻¹（n = 21） |
| $\Delta E(B-A)$ | $-11.567$ 到 $-10.758$ meV（紧参考 $-11.721$） |
| 势垒 | **+361 meV** |

### 冻结时**明确不成立**的七条（不要引用成结论）

1. $n=16$ 的 16/16 **不是"接近 1"** —— 95% 下界 0.81。区间必须与点估计同行。
2. **B1b 记"不评"**（4 个含 E/Z 通道的 seed）。标签早于分类器修复；只重跑那四个
   会重新制造混合实现批次。
3. $Q_1$ 从 12/16 升到 16/16 **归于 relay 定向是旁证，不是隔离归因** ——
   三处改动同批进入，没做关掉定向的对照运行。
4. **0.81 meV 是本协议下的端点散布，不是能量分辨率** ——
   要成为分辨率必须把所有端点用同一紧容差重新松弛。没做。
5. 直接观测与续行是**两类强度不同的证据**，不可合并。
6. **seed 71 的 C4/D1/D2 是真实失败**（它没有任何 TS），已测明与 ceiling 无关。
7. **openmm 8.6.0.dev 下 backend parity 为何变红仍未查明** ——
   现在绿是因为环境回滚到 8.5.2。

---

## 4. 关键机制，以及它们为什么长这样

每一条都是从一次失败里长出来的。改之前先读对应章节。

| 机制 | 解决什么 | 在哪 |
| --- | --- | --- |
| 解析 Hessian（`hessian_source="auto"`） | MACE 的力本身是 autograd 梯度，所以二阶导已经在图里；比 FD 快 2 倍且到机器精度 | `PRRS_STATUS.md` §2.6b |
| 刚体地板 + 梯度扣除 | $\dot x^\top H\dot x=\omega^2(g\cdot x_\perp)$，转动曲率**正比于残余梯度**，扣掉后残差 $\sim10^{-16}$ | §2.6c / §2.6g |
| 目标模跟踪 | 不给目标方向就爬**全局最软**模：3-氧代丁醛上返回 8 个 $-65\sim-80$ cm⁻¹ 的转子鞍点，而质子转移在 $-3186$ | §8.16 第 1 条 |
| 切向的质量度量 $M^{-1}g_x$ | `internal.gradient` 返回协向量；直接用会给 O 比 H 大 4 倍的权重（实测 H/O = 15.87 = 质量比） | §8.16.1 |
| `flat_biased` 档 | 曲率在地板下**而梯度不为零**的坐标不是平的。一个甲基转子因此冻住 9.97 meV | §8.18 |
| boundary continuation | 鞍点两侧下降证明"两极小相连"，**不**证明"交付的运动到得了产物" | §8.19.2 |
| `arrival_through_saddle` | 一阶鞍点在恰好两个盆地的边界上；两侧同名 ⇒ 该点邻域属于那个态 | §8.20.1 |
| relay 定向对齐 | 本征向量符号任意；不对齐会让下降在 relay 后倒回来 | `P2_TEMPERATURE_SCAN_PROTOCOL.md` |
| `descend_saddle` | 两侧下降的**唯一**实现，`connect_saddle` 与续行共用 | `runner.py` |

---

## 5. 现在的状态与下一步

### P3：苯甲酰丙酮 + 离子化学（**判据仍未冻结**）

**探针已跑过一轮，记录在 `docs/P3_BENZOYLACETONE_PROBES.md`（378 行）。**
那份文档里没有一条是验收结果 —— 它测的是"判据能被写出来之前必须先知道的东西"。
下面只列结论，数字与出处全在那里。

已测掉的：身份层能处理芳香环（Kekulé 对返回恰好 2 个解，环上键正确地不被 lock）；
22 原子成本 42.5 ms/点、解析 Hessian 4.2–4.5 s、~1350 MiB/worker；
**闭合螯合盆地只有 1 个构象**（否掉了"更多软构象"这个预期，`conformer_max_per_node = 8` 够用）；
TS 存在且 order 1，势垒 **+316 meV**、虚频 **−3086 cm⁻¹**，与 P1（419 / −3186）、
P2（361 / −3010）同一带；TS 处 r(O···O) 从 2.5289 收缩到 2.3234 Å（**0.2055 Å**），**重原子门控有实测值**。

**2D 投影已跑通**（`docs/figures/p3_projection.png`，48×48，relaxed 2304/2304 全收敛，
曲面容差与路径同为 2e-3）：双阱 + 中间鞍谷，阱在 r_OO ≈ 2.52、通道底部 2.32–2.36。
**投影间隙最大 +32.9 meV、平均 +1.4，仅 41/490 帧为负**（容差修正前是 375/491）——
按 §0 这就是 (q_PT, r_OO) 的成绩单：位置层面够，能量层面在 TS 附近缺约 33 meV。
**注意该剖面不含势垒顶点**：`descend_saddle` 先走一整步 `irc_step_A` 再淬火，
首帧已比 TS 低 184.9 meV。刚性减松弛 451–4187 meV 处处为正。

**职责边界（用户 2026-09-02 界定）：PRRS 只负责 TS 在哪里、产物是什么。**
势垒和虚频的绝对值是势能面的属性，不是搜索方法的属性，交下游 QM。
三模型对照支持这条界定：**scope 内的量三模型稳健，不稳健的恰好是声明不负责的部分** ——
`chemical_key` / `event_key` / `broken` / `formed` / saddle order **逐字段一致**，
r(O···O)@TS 差 0.036 Å；而**势垒差 3.7 倍**（+316.3 / +84.8 / +96 meV）、
**虚频差 3.3 倍**（−3086 / −926 / −1436 cm⁻¹）。与 SN2 探针同一模式：分歧在过渡态本身。

**两个必须先解决的**（2026-09-04：两条都已处理，见 `benchmarks/b03_benzoylacetone/README.md`；
下面的原文保留不动。**注意 §5 开头那句「判据仍未冻结」指的是本节写作当时** ——
判据已于 2026-09-03 冻结并判定通过，见 `P3_BENZOYLACETONE_ACCEPTANCE.md`；
现在待办的是**标准化收编**：尺度指标要写成带分母/阈值/pass-fail 的逐条判据，
以及补一份机械 `check.py`，P3 的判定至今是人读出来的）：
- **双侧连接只在 POLAR-1 上拿到。** `connects_two_nodes = True`，鞍点两侧分别落回
  A（`8710064d…`）和 B（`b1a4338a…`）—— P3 第一次拿到完整的产物身份证据。
  但在 **OFF24（P0–P2 的参照模型）上 `side +1` 仍然到不了 A**，
  即使 `quench_steps = 5000`，所以不是预算问题。
  跨模型对照指向那个近自由甲基转子：**OFF24 在 A 处的最软内坐标模 0.00262，
  POLAR-1 是 0.00470**，软将近一倍。所以这个失败是**那个势能面的属性**，不是方法的 ——
  但两个模型差的不只是这一个数，仍非隔离归因。
- **必须是 Z 异构体。** E 异构体在 300 个构象下 r(O···O) 全在 4.3–4.7 Å，**根本不能螯合**，
  而它会返回一整套自洽全绿的数字 —— 身份层看图，E 和 Z 的图一样。
  **冻结的 preflight 必须含一条几何门 `r(O···O)`。**

  > **2026-09-04 已建**：`benchmarks/b03_benzoylacetone/preflight.py`。
  > 实测 E **0/300** closed（r 4.396–4.559 Å）、Z 169/300（2.328–4.073 Å）、
  > 两个源结构 2.5193 / 2.5287 Å。数值与上面原记录的差别来自 ETKDG 种子与版本，结论一致。
  > 那条 `side +1` 的失败**记为 MACE-OFF24 这张势能面的已知限制**，
  > 保留"两个模型差的不只是这一个数，仍非隔离归因"这句限定；
  > 要重开需要在同一模型上单独改那个转子软度的隔离实验。

用户定的三类指标（跑前要冻结）：
- **正确性**：仍看节点 / 通道 / TS / 有向边
- **尺度**：单独看候选数、方向 inclusion、Hessian 成本
- **不再要求所有 seed 16/16**

**Track A（中性闭壳，现有管道，MACE-OFF24）**
1. 苯甲酰丙酮质子转移 —— 不对称 β-二酮 + 苯环 + 更多软构象 + 芳香图 + 竞争异构体
2. Diels–Alder（丁二烯 + 乙烯）—— **碎片数 2→1**，PRRS 至今只做过单分子内重排

**Track B（离子化学，需要 charge-aware 模型）**
SN2 → SN1 → E2 → EA。**MACE-OFF24 做不了**：它对总电荷完全不响应
（实测 charge 0/−1/+1 逐位相同，模型里没有任何 charge/spin 模块），
而域门会按构造拒绝每一个离子端点。

### 2D 势能面 + 路径叠加（用户要的可视化方向）

坐标对建议：苯甲酰丙酮用 $(q_{\rm PT},\ r_{\rm O\cdots O})$（显示重原子门控）、
Diels–Alder 用两条成键 C–C 长 $(r_1,r_2)$（对角线上是同步 TS）。

**成本**（实测 12 原子 52.2 ms/点，OFF24）：刚性 40×40 = 1.4 min；
松弛（~100 calls/点）40×40 ≈ 2.3 h、60×60 ≈ 5.2 h。

**必须处理的方法问题**：PRRS 路径是 $3N$ 维的。投影到 2 个坐标后，
**画位置没问题，画能量不能当成曲面的能量** —— 松弛曲面在同一点的能量是其余
$3N-2$ 个坐标都松弛过的最低值，路径一般高于它。正确画法是两联图：
等高线上画路径位置，旁边画路径自身能量 vs 弧长**并叠加曲面沿同一投影线的能量**。
**两条线之间的间隙本身是诊断量** —— 它回答"选的两个坐标有没有张开这个反应"。

### SN2 结果（`runs/sn2_probe.json`，`docs/experiments/sn2_probe.py`）

气相同一性 SN2，Cl⁻ + CH₃Cl。参考态是**一次计算、总电荷 −1、离子放 20 Å**，
所以全局电荷嵌入的常数偏移在每个数里抵消。D3h 鞍点由**在二参数对称子空间内取能量极小**
得到 —— 在 D3h 几何上所有非对称方向的梯度按对称性自动为零，所以这是完整势能面上的
真驻点，不是近似。**内部判据是 Hessian 必须给出恰好一个负本征值，两个模型都通过。**

| | OMOL-0-4M | POLAR-1-M |
| --- | --- | --- |
| CH₃Cl $r$(C–Cl) | 1.7816 Å | 1.7816 Å |
| 复合物 $r$(C···Cl⁻) | 3.1355 Å | 3.1767 Å |
| **复合物深度** | $-9.57$ kcal/mol | $-9.23$ kcal/mol |
| `confirm_minimum` | True, order 0 | True, order 0 |
| **D3h 鞍点** $r$(C–Cl)×2 | 2.4236 Å | **2.2815 Å** |
| **势垒 vs 参考** | $+13.24$ | $+4.53$ |
| **势垒 vs 复合物** | $+22.81$ | $+13.76$ |
| 虚频 | $-149.5$ cm⁻¹ | $-535.6$ cm⁻¹ |
| Hessian 判定 | first-order saddle | first-order saddle |

**不依赖任何文献值的结论**：两个模型在**复合物上只差 0.34 kcal/mol**，
在**鞍点上差 8.71 kcal/mol**，鞍点键长差 0.14 Å。所以它们的分歧**不在长程静电**
（那在复合物距离上才起作用，而那里它们几乎一致），而在**五配位碳的 TS 本身**。
虚频差 3.6 倍，是第二个独立的判别量。

**依赖文献值的结论（暂不成立）**：我记忆里的气相值是复合物 $-10.5$、
势垒相对参考 $+3.0$（相对复合物约 $13.5$）。按这组数，POLAR-1-M 的
$+13.76$ 只差 0.26 kcal/mol，OMOL-0 的 $+22.81$ 差 9.3。
**但本次会话没有核实这两个数** —— 一次网络检索没有找到 Cl⁻+CH₃Cl 的具体值
（只找到 F⁻+CH₃Cl 的 15.8 / 9.6 / 3.3±0.3 kcal/mol，那是**另一个体系**）。

所以"POLAR-1-M 复现了 SN2 势垒"目前是**条件命题**。要坐实它，读这两篇里的数字：

- Gas-Phase Identity SN2 Reactions of Halide Anions and Methyl Halides,
  JACS 1996 — <https://pubs.acs.org/doi/abs/10.1021/ja9620191>
- Model Identity SN2 Reactions CH3X + X⁻ (X = F, Cl, CN, OH, SH, NH2, PH2),
  J. Phys. Chem. A — <https://pubs.acs.org/doi/10.1021/jp054734f>

两篇都在付费墙后（抓取工具进不去），有机构访问权的人一分钟就能读到
复合物深度、势垒、D3h 键长和虚频。**在那之前不要把 POLAR-1 的势垒当作已验证。**

### 模型现状（`/home/ruigengji/MLP/mace/`，1.1 GB）

| 模型 | 电荷 | 原子电荷输出 | Cl/Br/F EA 偏差 (eV) | 尺寸一致性 | 备注 |
| --- | --- | --- | --- | --- | --- |
| MACE-OFF24_medium | **无响应** | 无 | — | — | P0–P2 全部用它，**不要换** |
| mace-omol-0-xl-4M | 有 | **无** | $-0.021$ / $+0.013$ / $\mathbf{-0.858}$ | $-0.403$ kcal/mol 逐位恒定 | 78.6 ms/点（OFF24 的 1.51×） |
| MACE-POLAR-1-M | 有 | **本轮取不到，见下** | $+0.085$ / $+0.078$ / $+0.119$ | $-0.80$ + 距离依赖项 | 需 `graph_longrange 0.4.0`（已装） |

**`Qs` 在 P3 这一轮取不到**（本地 CPU 与节点 GPU 两侧一致）：只拿到 `dipole`，
`'Qs' in atoms.arrays` 是 `False`。键名 / 触发条件 / 尺寸依赖未查明。
所以本表"有，且正确定域"那一格来自更早的一次测量，**P3 没能复现它**，
而"SN1 只有 POLAR-1 能检验电荷去了哪"这条推论目前**不可用**。
偶极不能替代原子电荷 —— 偶极不变不排除电荷在内部重排。
详见 `P3_BENZOYLACETONE_PROBES.md` §3.14 与 §5.7。

**POLAR-1 的推理**（`model_type` 必传，否则走不到极化场路径）：

```python
calc = MACECalculator(model_paths=".../MACE-POLAR-1-M.model",
                      model_type="PolarMACE", device="cuda", default_dtype="float64")
atoms.info["charge"] = 0        # 总电荷
atoms.info["spin"] = 1          # 总自旋，约定未验，闭壳一律用 1
atoms.info["external_field"] = np.zeros(3)   # V/Å，可省
# 输出：calc.results["dipole"]（Debye）、atoms.arrays["Qs"]（Fukui 均衡后的原子电荷）
```

**已验证的三件（都有数）**
- 电荷完整穿过 PRRS 管道：`snapshot` 用 `atoms.copy()` 保留 `.info`，quench 正常
- POLAR-1 电荷正确定域：20 Å、总电荷 −1 → 氯上 $-0.987$、CH₃Cl 上 $-0.013$
  （**SN1 是异裂，只有 POLAR-1 能被检验电荷去了哪；OMOL-0 不输出原子电荷**）
- 离子–偶极复合物 Cl⁻···CH₃Cl（OMOL-0，完整松弛，`confirm_minimum` 确认极小）：
  **$-9.57$ kcal/mol**，文献约 $-10.5$

**未解决的三件**
- **`info["spin"]` 的约定不明**：OMOL-0 上 H₂ 在 `spin=3` 时比 `spin=1` **还低** 0.9 meV，
  物理上该高约 10 eV。**开壳体系在搞清之前不要碰。**
- **POLAR-1 的长程项机制不明**：`×r²` 单调增长（103→320），与用模型自己偶极算的
  离子–偶极理想值之比从 3.5 涨到 10.8 —— **不是 $1/r^2$**。减掉常数后大致 $1/r^{2.6\text{–}3.3}$。
  它在复合物距离（3–8 Å）上贡献 1–3 kcal/mol，任何结合能都带这个系统量。
- **SN2 已跑通，但文献值未核实。**见下节。

---

## 6. 本轮被检验杀掉的推断（保留，因为它们是方法记录）

不要重复这些。每一条都是"我有一个自洽的解释，然后一个测量把它否了"。

| 我推断的 | 杀掉它的测量 |
| --- | --- |
| 更快的 GPU 会抬高并发天花板 | RTX 5080 单点快一倍，天花板从 3.2× **降到** 1.7× |
| 争抢的是 CPU 线程 | 每进程一线程 + 绑一核：1.70× → 1.72×，**没变** |
| CUDA MPS 能解除序列化 | 交错 A/B：**0.15×**，慢 6.8 倍 |
| 那批 seed 被主机 OOM killer 杀了 | 125 GB 空闲、内核日志无 kill；日志里一直写着 `ValueError` |
| 八个 seed 被外部中断了 | **它们还在跑**。`pgrep` 跑在了共享文件系统但不共享进程表的机器上 |
| $Q_1$ 与 $Q_4$ 会互相交换 | 边界帧 vs 直接观测 $r=-0.891\,[-0.96,-0.71]$ 真的竞争；传导到 TS 的是 $-0.062\,[-0.54,+0.45]$ |
| POLAR-1 的距离漂移是真实离子–偶极物理 | 与理想 $-q\mu/r^2$ 之比从 3.5 涨到 10.8，**不是那个律** |
| 苯甲酰丙酮"更多软构象"，reservoir 容量要加 | 闭合盆地只有 **1** 个构象，80 个确认极小同一 key、跨度 0.003 meV |
| 我那条 SMILES 是螯合烯醇 | rdkit 判 **STEREOE**，300 个构象 r(O···O) 全在 4.3–4.7 Å，**0** 个闭合 |
| A 与 conf265 是相距 1.255 meV 的两个极小 | 紧容差下**收敛到同一点**，ΔE = 0.0000 meV |
| A 的负模是苯环扭转（22% 在环上） | 读错了本征向量（`vec` 是**列**索引且质量加权）。改对后 **99.83% 在甲基** |
| `climb_seed` 是个能推动走法的旋钮 | `seed` 在 `follow_min_mode` 函数体里从未被使用（`inspect` 证实） |
| 那 8 对重复会逐字段相同 | **6/8 相同**；amp 0.50 差 1.2e-7 eV，**amp 0.10 结论完全不同**（GPU 非确定性） |
| `TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD` 是 shell 继承的，pop 掉即可 | **mace 自己设**（`mace/__init__.py:5`），`import mace` 那一刻就设回来 |

**最严重的一次不是判断偏差**：我在一台看不到那些进程的机器上读 `pgrep`，
把"本机不可见"读成"没在跑"，据此建议 `rm -rf` 七个**正在跑**的作业目录。
教训写死在 `scripts/scan_status.sh` 的注释里：**跨主机时 liveness 不可知就不要声称。**

---

## 7. 文件地图

```
src/prrs/                 17 个模块，6008 行
  runner.py               单次试验、淬火、Hessian、min-mode、续行、路径记录
  search.py               ReactionSearch 主循环、通道/边登记、鞍点搜索
  chemistry.py            三层身份、classify_transition、stereo_only_difference
  network.py              Registry：接纳、微观态水库、域门
  perturbations.py        方向提案、优先级带、家族配额
  config.py               全部参数 + 校验（跑前拒绝不合法配置）
  validation.py           独立交叉校验（ase.vibrations）—— 搜索侧从不 import
tests/                    13 个文件，172 passed / 1 skipped / 0 failed
docs/
  HANDOFF.md              本文件
  PRRS_STATUS.md          1978 行，主叙事，§8.x 是每一批落地的干预与其证据
  PRRS_CURVATURE_ROUTES.md 1270 行，路线 A/B/C 与 §E 架构结论
  P1_MALONALDEHYDE_ACCEPTANCE.md   507 行，冻结判据 + 八轮记录
  P2_OXOBUTANAL_ACCEPTANCE.md      882 行，冻结判据 + 七轮 + 两次判据修订
  P2_MULTISEED_PROTOCOL.md         241 行，T=0 扫描（**已作废**，理由：旋钮不动机制）
  P2_TEMPERATURE_SCAN_PROTOCOL.md  593 行，T=300 协议 + 结果 + ceiling 诊断 + P2 冻结
  P3_BENZOYLACETONE_PROBES.md      378 行，P3 探针记录（**判据未冻结，无一条是验收结果**）
  experiments/            40 个脚本，全部只读运行产物或做单点测量
  figures/                5 张图 × SVG/PDF/PNG + make_figures.py（读盘，不算新东西）
scripts/
  p3_node_bundle/              P3 的重活，全部分片 + PER_GPU=8 工作池
    run_all.sh                 ts / quench-ab / models / pes / landscape / merge / status
    _env_guard.py              float64 断言、TF32 关闭、模型按 hash 认、torch.load 静音、hostname
    p3_ts_barrier.py           TS / 势垒 / 虚频 / 两侧下降
    p3_three_model.py          OFF24 / omol-0 / POLAR-1 对照，只报模型内部的差
    p3_pes2d.py                (q_PT, r_OO) 刚性 + 松弛曲面，含路径投影
    p3_landscape_full.py       300 构象松弛 + Hessian 确认
    p3_merge.py                合并分片；**缺片会报并以非零退出**；拒绝混容差
  run_p2_temperature_scan.sh   check / probe / run / provenance / report
  scan_status.sh               status + 进程存活 + trial 数三项交叉判断
  _saturation_probe.py         并发天花板测量（含两块卡的实测值）
  _provenance.py               混合实现批次会 exit 1
runs/                     442 MB；哪些计入结论见 runs/README.md
models/MACE-OFF24_medium.model   随仓库走，sha e5ccf5837f685899
```

## 8. 环境与硬件

```
python 3.12.13   Linux-7.2.0-1-cachyos   /home/ruigengji/miniforge3/envs/openmm_dev
mace-torch 0.3.16   graph_longrange 0.4.0   torch 2.12.0   ase 3.29.0
numpy 2.4.3   openmm 8.5.2   e3nn 0.4.4
本机 GPU  RTX 2080 Ti 11 GB    单点 52.2 ms（12 原子，OFF24）
节点 GPU  RTX 5080   16 GB     单点 23.8 ms
```

**并发实测（`scripts/run_p2_temperature_scan.sh probe`）**

```
RTX 2080 Ti   天花板 3.2×   （7 并发 3.7×、20 并发 3.2× —— 过 7 就不涨）
RTX 5080      天花板 1.7×   （1→24 并发：单点 23.8→328 ms，聚合 1.72×）
                            32 并发在 16 GB 上 OOM
```

**更快的卡天花板更低**，这排除了 GPU 计算是争抢源；绑核 + 单线程毫无变化，
排除了 CPU 线程争抢。剩下的嫌疑是多进程共享一卡的 CUDA 上下文切换，**未测**。
**不要开 MPS**（实测 0.15×）。杠杆是**多卡**，不是单卡堆并发。

## 8.1 P2 的 relay 方向判定：19 条不可判定，350 条已证明稳健

全仓 369 条 `travel_overlap` **全部来自 P2**（P3 一条没有）。
2026-09-03 由另一个 session 提出、本 session 独立复核。

**这件事不用 GPU 就判完了**，方法是换一个提法：不问"重叠有多小"，而问
**"投影误差是否大于重叠本身"**。两种投影的差异完全由质量加权质心位移决定，
那是纯几何量，逐帧存在 `runs/p2_*/paths/ts*/descent.extxyz` 里，不需要任何势能评估。

### 数据是双峰的，而且分界不是人挑的

| | 条数 | |
| --- | --- | --- |
| 下簇 | **19（5.1%）** | max \|overlap\| = **0.027051** |
| 上簇 | **350（94.9%）** | min \|overlap\| = **0.935910** |
| 间隙 | **0 条** | [0.0271, 0.9359]，宽 **0.9089** |

本 session 用了一个不预设阈值的检验：把 369 个值排序后找最大相邻跳变。
**它恰好落在 19/350 的分界，宽 0.9089，而第二大跳变只有 0.0148（小 61 倍）。**
所以这个划分是数据自身的结构，不是选出来的。

### 投影误差完全落在间隙内部

348 条下降轨迹上独立复现新旧两种投影并比较单位方向：

| | 本 session（24676 个 travel 向量） | 对方（49982 个） |
| --- | --- | --- |
| p0 | 0.0043 | 0.0043 |
| p50 | 0.1157 | 0.1029 |
| p100 | **0.3220** | **0.4420** |

**向量计数差约 2 倍、p100 差 37%，但结论对此不敏感**（两组数字都满足下面两条）：

- **350 条：裕度 0.9359 远高于误差上界（0.32 / 0.44）** → 符号不可能被这个缺陷翻转。
  **这一条现在是被证明的，不再是被假设的。**
- **19 条：裕度 ≤ 0.0271 低于误差的 5 分位（0.029 / 0.051）** → 记录下来的符号是由一个
  小于伪影的量决定的。

### "不可判定" ≠ "符号错了"

要真正定下这 19 条的符号需要 relay 的本征向量，**那个没有存盘**，所以那一步才需要势能。
复核脚本：`docs/experiments/p2_relay_overlap_margin.py`（只报告，不 gate）。

**对复算的意义：需要复算的是 19 条而不是 369 条**，`p2_round7` 的 3 条优先
（round7 是 P2 最终 T=0 轮，16/16 判据就是在它上面评的）。P2 冻结在 `0058158d`，
要动得走显式修订流程。

与 §8.16.1（切向的质量度量要用 $M^{-1}g$ 而非 $M^{1/2}g$，实测 H/O = 15.87 恰是质量比）
并排看：**这是同一条线上的第二个质量加权缺陷。**

## 9. 技术债

- **manifest 不记 hostname**，只有内核串。跨主机判断运行状态靠它碰巧能用。
  改它会改变 `implementation_sha256` —— P2 已冻结，**现在可以做了**。
- **扫描不可续跑。** `run` 拒绝非空目录，一次外部中断丢掉在飞的全部工作。
  要做需重建 registry 状态，而**接纳顺序会影响结果**，做错比不做更糟。
- **这个目录不是 git 仓库。** 所以"两份实现只差记录层"这类判断**无法独立核验** ——
  P2 的混合实现批次就是因此只能作废重跑。**建议接手后第一件事是 `git init`。**

  2026-09-04 又收到一次证据，而且比 P2 那次更难处理：`make hash --check` 报出
  `runs/` 下 **16 个实现指纹 / 107 个 run**，分布极偏 —— `d7726173` 占 26 个、
  `0058158d` 占 19 个，**其余 14 个实现分摊 62 个 run**，也就是一大半 run 的实现
  只出现过一两次。`scripts/_provenance.py` 只能在**一批之内**发现混批，
  跨批的实现漂移没有任何东西在看。

  同一天的逆 DA 探针把这件事推到了极限：扫描**开始时**树是 `0a7aa55dd2aa556e`，
  **跑完时**是 `23179deadf657452` —— 另一个 session 的注释轮次落在了运行期间。
  这次无害（那一轮只改注释；探针只 import `prrs.calculators` 与 `prrs.state`，
  两个文件在整轮重构中 **AST 逐字段相同**，在该运行窗口内**逐字节相同**），
  但它说明**"批"这个单位本身已经不成立**：漂移可以发生在一次运行的中间。
  一个 sha256 只能说"树变了"，说不出变的是注释还是判据 —— 那正是 git 提供的东西。
- **openmm 8.6.0.dev 下 4 个 backend parity 测试为何变红**，未查明。
- **`torch.load` 的包裹已在 `src/prrs/` 内收敛成一处**（`src/prrs/torch_guard.py`，
  2026-09-03）：单一 wrapper、安装一次、**没有 uninstall 路径**，device-move 改成同一个
  wrapper 上的作用域深度计数，所以退出作用域时无物可还原 ——
  原先那个"finally 静默吃掉另一处"的失效模式**结构上不再存在**。

  **两份实现现在通过一个共享信任表协作**（`_prrs_trusted` 属性），三条规则，
  每一条都是先量到失败才成为规则：

  1. 已有 guard 且暴露共享表 ⇒ **把自己的信任路径加进去，不再重复包裹**。
     `_env_guard` 原先看到 `_prrs_wrapped` 就 early-return 而不贡献任何东西，
     所以 `torch_guard` 先装时 bundle 的信任列表**完全不生效** ——
     而三模型探针要加载 omol-0 与 POLAR-1，那两个 `torch_guard` 没有理由信任过。
     用假 torch 的反事实实测：旧逻辑下两者都拿到 `weights_only=True` 而加载失败，
     **是静默失败不是告警**。
  2. **下面还有 guard 时，绝不强加限制性默认** —— 它可能知道那个文件可信。
     反向的同一条曾是真实的 `UnpicklingError`：外层设了 `weights_only=True`，
     内层按设计看到参数已存在就不再干预，于是一个内层信任表允许的 checkpoint 加载失败。
  3. 暴露 `_prrs_trusted`，让**后装**的 guard 能对自己做规则 1。

  两种安装顺序实测（三个 checkpoint：OFF24 / omol-0-xl-4M / POLAR-1-M，
  `-W error::UserWarning` armed）：

  | 顺序 | 共享表 | 三个 checkpoint | warning |
  | --- | --- | --- | --- |
  | bundle 先装，`torch_guard` 后 | 3 条 | 全部 OK | 零 |
  | **`torch_guard` 先装（信任表空），bundle 后** | 3 条 | 全部 OK | 零 |

  第二种顺序在规则 1 之前会静默失败。

  **`intact()` 一类的自检只能测到失效的一半。** 两个方向：
  *REMOVED*（第三方 finally 还原掉 guard）与 *NESTED*（第三方包在外层、guard 仍在链上工作）。
  比较 `torch.load is <自己的 wrapper>` 对两者都返回 False，后者是假警。
  `_env_guard` 的那个函数已按此更名为 `torch_load_outermost()`，
  env 字段改为 `torch_load_wrapper_outermost`，**False 是"去看一眼"的提示，不是判决**。

- **`tests/test_backend_parity.py` 有一个"静默不运行"的失效模式，比变红危险。**
  该文件顶部有四个 `pytest.importorskip`（`openmm` / `openmmml` / `mace` / `torch`）
  加两道 gate（`models/MACE-OFF24_medium.model` 存在、CUDA 可用）。
  **在只装了 numpy/ase/pytest 的环境里整个模块收集到 0 条**，pytest 照样报"全绿"。

  实测对照（2026-09-03，两个 session 各跑一次同一份代码）：

  | 环境 | 结果 | backend parity |
  | --- | --- | --- |
  | 只有 numpy/ase/pytest | 174 passed, 1 skipped | **0 条被收集** |
  | 全栈 + CUDA + 模型文件 | **180 passed**, 1 skipped | 6 passed, 1 skipped |

  差的恰好是这 6 条。**所以"测试通过"在不同环境下含义不同**，而变红至少看得见。
  交接或 CI 时 pytest 应带 `-rs`（列出 skip 原因），或显式断言收集到的用例数。
  唯一被跳过的第 7 条是 `test_report_backend_throughput`，需 `PRRS_RUN_BENCHMARKS=1`，
  那是有意的。

---

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
   没有一个被允许优化出自己的鞍点。要 ORCA `ExtOpt` 驱动 MLP，**尚未做**。
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

## 12. 2026-09-05 交接：三项正确性修复，与接手要先知道的事

树 `69d53501` → **`5b69dfed41c4bdb2`** · `233 passed, 1 skipped` · digest 绿（一条已审查迁移）
· 复跑基线在 `env/BASELINE.md`

### 12.1 接手先读这三份

| | |
| --- | --- |
| `docs/REVIEW_2026-09-05.md` | 用户的独立代码审查，逐项验收条件。**本节是对它前三项的执行。** |
| `docs/DECISIONS_PENDING.md` | 待决事项与迁移计划。**没做完的都在这里，别去别处找。** |
| `env/BASELINE.md` | 指纹、`versions()` 九项、复跑命令 |

### 12.2 修了什么（每项都有变异检验）

**① 验证不完整不再算通过**（`runner._probe_minimum`）

`_curvature_along` 被 gate 拒绝的方向原先被 `continue` 跳过，而
`confirmed = bool(measured) and not negative` —— **四个方向被拒、一个非负，就返回通过**。
反向也错：调用方读 `not confirmed` 后写 `not_a_minimum`，
于是"没测出来"被记成"不是极小"，几何被丢。

现在三态：`complete`（所有方向都读到）/ `confirmed`（且无负值）/
`verdict ∈ {minimum, not_a_minimum, unverified}`。
`unverified` 走 `status="unverified"`、`admission="withheld"`，几何进候选清单。

**② 候选出口统一**（`search.ReactionSearch.run`）

三处 `admit` 里，`execute` 会列候选，而**两侧下降与边界延续直接丢几何**。
现在一个 `withhold(structure, source, attempt, reason, detail, confirmed)`，四处调用，
结构写 `withheld/<序号>_<来源>.extxyz`，`source ∈
{trial_endpoint, saddle_descent, boundary_continuation}`。
**几何与出处先落盘，再决定进不进网络。**

**③ 失败落盘独立于网络快照**（异常处理）

旧顺序是先 `publish()` 再写 manifest。当异常**来自** `publish()`（网络里有非有限值，
`io.dumps_strict` 按设计拒绝）时，处理器再抛一次，manifest 根本写不下去，
`status` 停在 `running` —— 看起来像被信号杀掉。
**T=300 的三个 seed 就是这么死的**，`io.py` 记着现场，但当时只改了错误消息、没改顺序。

现在最小的先写：`failure.json`（纯字符串，含 traceback）→ manifest →
`publish()` 尽力而为，失败则在 `failure.json` 里记 `network_snapshot_failed`。

### 12.3 一条被推翻的分析，接手不要再引

**"b04 没有 bracket 是因为只有最高一级反应、没有下一级可比"——错的。**

```
0.1 → m0000 返回 · 0.3 → m0000 返回 · 0.6 rejected(collateral) · 1.0 → m0001 离开
```

细化只比相邻样本（`search.py`），(0.3,0.6) 与 (0.6,1.0) 都因 0.6 是拒绝点而不成对。
**是 collateral 拒绝隔断了唯一有效的 (0.3, 1.0) 对。**
若 0.6 有效返回，现有实现就能成 bracket —— 不用加档、不用改阶梯。

**所以 bracket 不独立于 collateral**，次序是先修扰动交付再评估扫描策略。

### 12.4 QM 复核这条线已经通了两例

| | |
| --- | --- |
| **P1**（`docs/experiments/p1_qm/`） | 双向 IRC 收敛、两端 0 虚频、正向端点对极小 RMSD 0.0004 Å；`event_key a04410b76e16b117` **复现冻结值** |
| **b04**（`docs/experiments/b04_qm/`） | 从 PRRS 的 TS 出发精化到 2.3645/2.3645，与独立求得的 QM 鞍点小数点后四位相同；双向 IRC 收敛、两端 0 虚频互为镜像；`event_key 847306d3f8c57bd4` 与 PRRS 一致；中心势垒 +13.02 kcal 两组几何相同 |

**IRC 在本项目里第一次真正执行**（§3.7 那条"由执行证实"此前只有文档）。

协议：**r2SCAN-3c 出几何与频率 → ωB97M-V/def2-TZVPD 在同一几何上补单点**。
r2SCAN-3c 单独用把 SN2 中心势垒错 11 kcal，**几何好、能量必须补**。
校准：SN2 上复现文献 0.4–0.9 kcal。

**两例同一形状：几何与事件身份属于方法且正确，能量与曲率属于模型。**

### 12.5 没做完的，按用户定的顺序

```
1 可复跑环境与版本记录、原版回归        ✅ env/BASELINE.md，本轮实跑
2 三项正确性问题                       ✅ 本节
3 单独修扰动交付，比较有效扰动比例与 bracket   ← 下一步
4 扫描策略 · b04 正式验收 · 扩展 B 组
5 效率用「每个经 QM 确认的独立通道成本」衡量
```

**第 3 步的具体建议**（用户）：对适用的桥键按分子图切分片段、整体移动，保持片段内部几何；
环内键与跨碎片接近分别处理。之后把扫描显式分为
**返回 / 离开 / 被阻断 / 未知** —— 不得把拒绝当返回，也不得无记录跳过。

诊断已完成：拉 C–Cl 0.600 Å 时只有 C（0.4482 Å）与 Cl（0.1518 Å）动，
三个 H 移动 **0.0000 Å**，每根 C–H 缩短 0.0628 Å；CH₃ 若随 C 刚体平移则变化 **0.000000 Å**。
**是位移实现，不是必要耦合。** torsion 已用 `rotate_fragment` 避开，bond 没有等价物。
契约已钉成测试：被驱动的键豁免、其余键的损伤仍拒绝。

### 12.6 接手要知道的三条规矩之外的事

- **这个目录仍然不是 git 仓库。** 2026-09-04 实测：一次扫描的**中途**树哈希变了
  （`0a7aa55d` → `23179dea`），连"批"这个单位都不成立。
  `runs/` 下 16 个实现 / 107 个 run，14 个实现只分摊 62 个。
- **`benchmarks/` 是新建的反应基准套件**，六条规矩在它的 `README.md`。
  b00–b03 已转 `fixture`（判据冻结过、判定通过过，验收职责移交 P4+ 的 QM 协议）。
  **b04 仍缺正式判据与机械评定**（`CRITERIA.md` / `check.py` 是空壳）。
- **默认模型从 P4 起是 `MACE-POLAR-1-M`**（用户 2026-09-04）。
  P0–P3 冻结在 MACE-OFF24 上，不重跑、不改标签。
  依据是两组 QM 对照，不是偏好：POLAR-1 在中性有机质子转移上偏 +12 meV、
  在离子 SN2 中心势垒上偏 +2.1 kcal，是唯一两边都成立的候选。
