# 外部参考

一律标注**本项目内是否核验过**。未核验的数字只能用来让比较可证伪，不能当结论引用。

---

## ACS Catal. 2026, 16, 8053–8066

Quirion, Kong, Stanley, Joy, Ess ·
*Density Functional Theory Surrogate Enables Fast and Broad Computational Evaluation of
Homogeneous Transition Metal Catalytic Energy Landscapes* ·
`cs-2025-090935.pdf`（工作区根目录）· **本项目未复现其中任何数字**

### 它做了什么

用 UMA（Universal Models for Atoms，Meta/FAIR；eSEN 架构，训练于 OMol25）
当 DFT 代理，评估五个均相过渡金属催化循环的能量面：
Ir 烷烃脱氢、Rh 氢甲酰化、Ru 烯烃复分解、Pd Buchwald–Hartwig、Cu 二氟卡宾插入。

对 ωB97M-V/def2-TZVPD 的平均绝对误差（论文自报，未在此核验）：

| 体系 | UMA-S | UMA-M |
| --- | --- | --- |
| Ir 烷烃脱氢 | 2.8 kcal/mol | 2.3 |
| Rh 氢甲酰化 | 1.9 | 1.6 |
| Ru 烯烃复分解 | 1.7 | 0.8 |
| Pd Buchwald–Hartwig | ~1.7 | — |

### 对本项目有用的四条

1. **模型层的存在性证明。** OMol25 覆盖电荷、自旋多重度、>80 元素，
   且论文明确指出它含 "systems with partially ruptured or forming bonds,
   which enables the calculation of transition structures"。
   这条论证的载体是**数据集不是架构**，所以它直接传递给 `MACE-omol-0`
   （同一数据集，MACE 架构）。**不自动传递给 POLAR-1** —— 那是另一个模型
   （`PolarMACE`，83 元素，`Default` head），训练集需另行查证。

2. **一个记录在案的失效模式。** UMA 无法优化 Cu 循环里的 `TS1e` ——
   论文判断原因是它"架在闭壳（2e）与开壳单重双自由基（3e）之间"，
   而把 3e 当单重态优化会给出非物理的构象剧变。
   **对本项目的意义**：任何可能走双自由基的通道，
   模型给出的"找不到"和"这里没有"是两回事，必须分开记（`README.md` 规矩 4）。

3. **验收协议可以抄。** 每个结构都做频率分析确认驻点性质
   （极小 = 无虚频，过渡态 = 恰一个虚频），并且把
   `UMA//DFT 单点` 与 `UMA//UMA 优化` 分开报 —— 几何误差与能量误差因此可分离。
   PRRS 的 `confirm_minimum` 已经在做前半条。

4. **溶剂化是加性修正。** 用 GFN2-xTB 的 ALPB 作为 UMA 能量与梯度上的加性项。
   本项目第一版是气相，这条现在用不上，但 B 组（离子化学）真做起来时
   气相与溶液的差别不是小量，届时这是一条已被验证过的做法。

### 它没有做、也因此没有触及 PRRS 新颖性的一条

**论文里每一个中间体和过渡态都来自前人已发表的机理。**
图 1–6 的结构取自 Goldman、van Leeuwen、Jensen、Norrby、Koh/Liu 等人的工作，
UMA 的任务是把它们重新优化并给出能量。**机理是输入。**

PRRS 不给产物、不给反应坐标、不从 TS 猜起（`docs/HANDOFF.md` §1）。**机理是输出。**

所以这篇论文与本项目的关系是：**它消掉了模型层的风险，同时没有碰搜索层。**
引用它时应当只引第 1–4 条，不要把它读成"这件事已经有人做过了"。

### 一条不成立的差异化

UMA **也是等变的** —— 论文原文 "UMA are equivariant Smooth Energy Network (eSEN)
based models"。所以"我们是 E3 等变"**不是**相对 UMA 的优势，不要这样写。

---

## 气相 SN2，Cl⁻ + CH₃Cl

文献值（高水平气相计算，**本项目未核验**）：
离子-偶极复合物 −10.5 kcal/mol，势垒 +3.0 kcal/mol。（本项目未核验，且下表两个数与它不可直接比较——见下。）

本项目实测（`runs/sn2_probe.json`，`docs/experiments/sn2_probe.py`，2026-09-02）：

| | MACE-omol-0-4M | MACE-POLAR-1-M |
| --- | --- | --- |
| 复合物 vs 参考 | −9.57 kcal/mol | −9.23 |
| `confirm_minimum` | True，order 0 | True，order 0 |
| D3h 点能量 vs 参考 | **+13.24** | **+4.53** |
| 最低 Hessian 本征值 | **−0.0822** | **−1.0549** eV/Å²/amu |

**这两行都不能按字面读，前提没建立。**

那个 D3h 几何是在**两参数对称子空间**里用 Nelder-Mead 最小化能量得到的
（`xatol=1e-4, fatol=1e-7`），它的**残余梯度从未被测量** ——
驻点性是从对称性**论证**出来的，不是测出来的（`sn2_probe.py` 的鞍点段没有任何 `get_forces`）。

- ω² = λ/μ 要求驻点，所以那两个本征值**不得称为频率、不得报 cm⁻¹**。
  Hessian 是解析的救不了：**非驻点上的解析 Hessian 精确且无意义。**
  `docs/HANDOFF.md` §4 的刚体地板量化过同一件事：`ẋᵀHẋ = ω²(g·x⊥)`，曲率正比于残余梯度。
- 一个未证明是驻点的几何上读到的能量，也不是势垒高度。

**与 P1–P3 的 −3186 / −3010 / −3086 cm⁻¹ 不是同一类数。**
那三个经过 `quench` 到 `fmax` 再 `confirm_minimum`，梯度界是**测出来的**；这两个不是。

站得住的只剩一条：**两个模型在身份层一致，在过渡区分歧。**
这是 `docs/HANDOFF.md` §5 划定"PRRS 只负责 TS 在哪里、产物是什么"的实测依据，
也是 `README.md` 规矩 2 的来源。**不要据此说哪个模型更好** ——
没有参考计算，n = 1，而且这两个数本身前提未建立。

---

## T[E_QM]：本套件的 QM 参照协议

2026-09-04 建立，ORCA 6.1.1，`scripts/orca_remote.sh`。
两份实测：`docs/experiments/sn2_qm/README.md`、`docs/experiments/p1_qm/README.md`。

### 协议（两段，第二段不是精修而是必需）

| | 方法 | 用途 |
| --- | --- | --- |
| 几何 + 频率 | r2SCAN-3c / TightSCF | 结构、驻点性、虚频个数 |
| 能量 | ωB97M-V / def2-TZVPD / DEFGRID2，**在上一段的几何上** | 势垒 |

**校准**：SN2 上复现气相文献到 0.4–0.9 kcal/mol。这是这把尺子可信的依据。

**r2SCAN-3c 单独用不行**：SN2 中心势垒偏低 **11 kcal**，丙二醛偏低 **72 meV**，
同一方向 —— meta-GGA 的自相互作用误差。它的**几何**是好的，能量必须补单点。

### 驻点性必须测量

两个体系的鞍点都满足 `README.md` 规矩 6：优化到收敛、频率确认恰好一个虚模。

| | RMS grad | 虚频 | 次低模 |
| --- | --- | --- | --- |
| SN2 D3h | 1.14e-05 | **1**，−314.42 cm⁻¹ | +195.87 |
| 丙二醛 TS | 7.74e-05 | **1**，−909.03 cm⁻¹ | +363.49 |

对比 `runs/sn2_probe.json`：那里的 D3h 点驻点性是从对称性**论证**的，
残余梯度从未测量，所以那两个本征值不得称频率。**这一节的可以。**

### 两组结果：模型对 T[E_QM] 的偏差（同一组几何，零优化器混淆）

| | 丙二醛质子转移 | 离子 SN2 中心势垒 |
| --- | --- | --- |
| MACE-OFF24_medium | **+276 meV** | 做不了（对电荷逐位无响应） |
| MACE-omol-0-4M | −1.8 meV | **+10.46 kcal/mol** |
| **MACE-POLAR-1-M** | **+12.1 meV** | **+2.09 kcal/mol** |

> **两个数不要混。** 上表的 `+276 meV` 是**固定几何**下的模型误差（同一模型在 DFT 几何上的单点 410.4 对 QM 的 134.4）。
> **P1 的势垒是 419 meV，相对 QM 的偏差约 285 meV** —— 那是 PRRS 在**自己找到的几何**上的值。
> 涉及 P1 势垒时用后者。


**用户 2026-09-04 定：从 P4 起统一用 POLAR-1。** 它是唯一两边都成立的候选。
P0–P3 冻结在 MACE-OFF24 上，不重跑、不改标签。

### 三条要记住的

1. **PRRS 被这组数洗清了，不是被质疑。**
   `A_θ[E_OFF24] = 419.1 meV`（P1 冻结）、同模型在 DFT 几何上 `410.4 meV`、
   `T[E_QM] = 134.4 meV`。**搜索找到的就是它那个模型自己的答案**（差 8.7 meV），
   285 meV 全部落在模型那一行。这是那个拆分第一次真正产出归因。
2. **"专用模型在自己域内更好"在这里不成立。**
   MACE-OFF 专训中性有机，却在中性有机质子转移上比两个通用模型差一个数量级。
3. **n = 2，且只做了 DFT 几何上的单点。** 这**不是** T[E_m] ——
   没有任何模型被允许优化出自己的鞍点。那要 ORCA `ExtOpt` 驱动 MLP，尚未做。
