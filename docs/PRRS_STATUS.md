# PRRS 实现状态

日期：2026-08-31  
基线：`84 passed, 1 skipped`；乙醇 P0 基准 `13/13`（`/home/ruigengji/miniforge3/envs/openmm_dev/bin/python -m pytest`）  
配套文档：[PRRS_ENGINEERING_PLAN.md](PRRS_ENGINEERING_PLAN.md)（计划稿，部分结论已被实测推翻，见第 6 节）  
曲率判定的两条路线（路线 A 无矩阵负曲率检测、路线 B 有限幅度响应算子 + Morse 指标；除两节实测外为理论组织）：[PRRS_CURVATURE_ROUTES.md](PRRS_CURVATURE_ROUTES.md)

本文只记录**已实测**的事实与**已落地**的代码。推断和待办分别标注。

---

## 1. 模块与职责

| 模块 | 职责 | ASE 依赖 |
| --- | --- | --- |
| `config.py` | 冻结的协议参数与校验；全部进 manifest | 无 |
| `internal.py` | primitive internal coordinate（bond/angle/dihedral）、精确刚体旋转、continuation 位移、精确 ΔK 冲量 | 无 |
| `chemistry.py` | `chemical_key`（自同构不变图哈希 + 碎片 + 构型立体化学）、自同构群枚举、对称感知 RMSD | 数据表 |
| `perturbations.py` | 六家族提案、对称商、分层方向预算、二维交付验收 | 无 |
| `network.py` | 两层网络（chemical node / conformer microstate）、reservoir 评分、双预算 | 无 |
| `state.py` | 活性子集图编码、Kabsch RMSD、响应分类 | `covalent_radii` |
| `reliability.py` | 每次力评估的失败关闭门控 | Calculator 基类 |
| `runner.py` | 单次试验编排、mode-resolved 淬火、Hessian/极小检验、鞍点下降 | 积分器/优化器 |
| `search.py` | 预算调度、幅度细化、TS pool、manifest | I/O |
| `openmm_backend.py` | 经实测验证的 openmm-ml/MACE 部署原语（交叉校验后端） | 无 |
| `validation.py` | 曲率层的独立交叉校验（ASE 振动机制），**搜索从不导入它** | `ase.vibrations` |
| `io.py` / `cli.py` | 持久化与命令行 | I/O |

搜索核心（前六个模块）**不依赖积分器**。ASE 保留在结构/分析/执行层，与计划书第 2 节的分工一致。

---

## 2. 已实测的后端事实

环境：`openmm_dev`（ase 3.29.0、openmm 8.5.2、openmm-ml 1.6、openmm-torch 1.5.1、mace-torch 0.3.16、torch 2.12 cuda129），GPU 为 RTX 2080 Ti。模型 `MACE-OFF24_medium`。

### 2.1 两个后端是同一个物理 Hamiltonian

四个构型（平衡、C–C 拉伸 0.8 Å 且 fmax=133 eV/Å、断键、rattle 0.15 Å）上：

```
最大 |dE|      = 2.610e-06 eV
最大 max|dF|   = 3.436e-05 eV/Å
有限差分梯度   最差 1.709e-05 eV/Å
```

在 numpy 2.2.6 与 2.4.3 上**逐位相同**，所以结论与 numpy 版本无关。

### 2.2 Context 可以跨键图变化复用

用过的 Context 与全新 Context 在断键构型上 `dE = 0.000e+00`、`max|dF| = 0.000e+00`。

原因是 `openmmml/models/macepotential.py:246` 每次评估都用当前坐标重建 MACE 邻域，没有冻结邻接。**计划书第 3 节第 9 条的疑问由此解决：不需要因键图变化重建 Context。**

### 2.3 OpenMM 路线没有吞吐优势

`macepotential.py:226` 用的是 `openmm.PythonForce` 而非 `TorchForce`（同文件 `:52` 的 docstring 已过时）。因此每步都回 Python 并做 GPU→CPU 同步，与 ASE 付同一笔开销。

三次测量 ase/openmm = 0.95 / 1.04 / 0.79 —— 单分子计时在共用 GPU 上有约 25% 批次间抖动，**诚实结论是"测不出差别"，不要引用任何具体比值**。

### 2.4 三个会静默改变物理的默认值（已在 `openmm_backend.py` 强制处理）

1. `createSystem` 默认加入 `CMMotionRemover`，会移除质心运动并破坏动量守恒验收 → 必须 `removeCMMotion=False`，且建完 System 后**清点 force 列表**确认。
2. 本地 `modelPath` 分支只做 `torch.load(map_location=...)`、缺 `.to(device)`（具名模型分支有），CUDA 上直接 device mismatch → 已用窄范围可逆补丁，并有回归测试钉住上游行为。
3. `except: platform = CPU` 式静默降级 → 平台不可用必须失败，不得回退。

### 2.5 MACE-OFF24 完全不响应总电荷

```
charge = 0   E = -4221.594512325 eV
charge = -1  E = -4221.594512325 eV
charge = +1  E = -4221.594512325 eV
```

逐位相同。**任何阴离子化学（如 SN2 水解）在这个模型下没有意义**，所以 `charge_sensitive = False`，且配置层拒绝非零 `total_charge`。元素覆盖为 `H C N O F P S Cl Br I`，`r_max = 6.0 Å`。

### 2.6c MACE-OFF24 的二阶导在乙醇上看不出噪声（四项无参考诊断）

脚本 `docs/experiments/second_deriv_reliability.py`（极小点）与 `docs/experiments/sdr_ts.py`
（鞍点），输出为同名 `.out`。"MLP 二阶导不可靠"是一个必须实测而非假设的问题；这四项诊断
**不需要参考数据**。

**鞍点的独立复现**：固定二面角 C–C–O–H = 0° 的约束极小化 → `fmax = 4.0e-05 eV/Å`、
$E = -4221.556394$ eV、相对极小点 $-4221.603394$ eV 即**势垒 47.0 meV，恰好一个负模
$\lambda_1 = -0.295964$**。§5 记的是同一个 syn 势垒 $+47.1$ meV / $-290.1\ \mathrm{cm^{-1}}$
—— 完全不同的路径（约束极小化 + 解析 Hessian）给出同一结果。

```
D1 不变性误差地板（严格零本征值的实测大小，eV/(A^2 amu)）
   极小点  平动 ~1e-15   转动 -3.7e-6 / -4.4e-6 / -2.6e-5   |max| = 2.5e-05
   鞍点    平动 ~1e-15                                       |max| = 2.9e-06
   对照    最低物理 |lambda| = 0.216 / 0.296    eigenvalue_tol = 1e-3
   ！这两个转动数字后来被证明不是模型误差，是残余梯度。见 2.6g。

D2 kappa_a 的偏差（鞍点，力形式，相对 lambda）
   a =            1e-3    1e-2    2e-2    5e-2    1e-1    2e-1    3e-1
   不稳定 -0.2960 0.01%   0.98%   3.9%    24%     95%*    354%    705%     (*此处变号)
   最低正 +0.2550 0.00%   0.13%   0.54%   3.4%    13%     53%     115%
   次低正 +0.6540 0.00%   0.00%   0.02%   0.12%   0.48%   1.9%    4.0%
   全部平滑单调，无小 a 抖动。

D3 Hessian 的粗糙度  ||H(R+dw)v - H(R)v||/d，v = 不稳定模，||Hv|| = 0.5252
   d =        1e-5      1e-4      1e-3      1e-2      1e-1
              22.6530   22.6510   22.6304   22.4229   20.2684
   跨两个数量级稳定到 4 位 -> 三阶导存在且平滑，没有噪声。
   但相对量是 43 / A：H 每 0.01 A 变化约 43%。

D4 float32 vs float64
   鞍点 lambda_1  -0.295964 (f64) vs -0.295966 (f32)   全谱 max|dlambda| = 1.5e-05
   float32 不改变负模的符号或大小。

D5 源码层排除  PolynomialCutoff (p=6) 在 r_max 处 f = f' = f'' = 0
   -> 邻域截断不产生 Hessian 跳变。f'''(1) = -336 != 0，只影响有限幅度的高阶项。
   加载的模型中与光滑性相关的子模块只有 PolynomialCutoff。
```

**结论**：在乙醇的极小点与鞍点上，"二阶导粗糙/有噪声"**没有得到证实**。误差地板
$10^{-6}$–$10^{-5}$，比 `eigenvalue_tol` 低 40–350 倍；三阶导平滑；float32 都够。
**所以 2.6b 第 2 条改用解析 Hessian 的前提成立。**

**但无参考诊断只能查不一致，不能查不准。**一个平滑而整体偏移的 Hessian 会完美通过
D1–D5。剩余风险是平滑偏差，唯一手段是外部参照（乙醇三个 TS 上一个 DFT 单点 Hessian）；
MACE-OFF 单模型，`Committee` 不可用（§6）。**在那之前只能声称曲率自洽，不能声称它准。**
且这只是**一个体系的两个几何**，成键/断键区（缺口 8）未测，不可外推。

### 2.6e 有限幅度下不存在对称算子，所以 $\kappa_a$ 不是本征值

脚本 `docs/experiments/finite_operator.py`（MACE 鞍点）、
`docs/experiments/finite_operator_analytic.py`（解析四次势）、
`docs/experiments/critical_points.py`（临界方向计数）。

代数：$D_a(v) = [g(x{+}av) - g(x{-}av)]/(2a) = Hv + \tfrac{a^2}{6}F[v,v,v] + O(a^4)$，
$F = \nabla^4 E$。于是

$$\langle v_1, D_a v_2\rangle - \langle v_2, D_a v_1\rangle
= \tfrac{a^2}{6}\big(F[v_1,v_2,v_2,v_2] - F[v_2,v_1,v_1,v_1]\big) \ne 0$$

而能量形式 $\kappa_a^E(v) = \langle v,Hv\rangle + \tfrac{a^2}{12}F[v,v,v,v]$ 是**四次型**。
三条式子在解析四次势上精确验证（非对称量测量/预测比值 1.00000）。

**MACE 鞍点实测**（$v_1,v_2$ 为 Hessian 特征向量，故 $v_1{\cdot}Hv_2 = 0$ 精确，测到的全是
$a^2$ 项）：

```
a        <v1,D_a v2>   <v2,D_a v1>   非对称量     /a^2
1e-2     1.982e-05     1.278e-04     -1.08e-04   -1.0796
2e-2     7.93e-05      5.11e-04      -4.31e-04   -1.0786
5e-2     4.95e-04      3.17e-03      -2.68e-03   -1.0716
1e-1     1.97e-03      1.24e-02      -1.05e-02   -1.0469
```

最小测试幅度上两者已差 **6.4 倍**；反推 $F[v_1,v_2,v_2,v_2] = 1.19$、
$F[v_2,v_1,v_1,v_1] = 7.67$，$(1.19{-}7.67)/6 = -1.0796$ 与实测逐位吻合。

**谱结构分两阶段崩掉**（合成 $N=8$，Newton 解 $\nabla_v\kappa_a^E = \mu v$，3000 随机起点）：
$a=0$ 时 8 个正交临界方向（= 特征向量）；$a=0.1/0.3/1.0$ 仍 8 个但
$|v_1{\cdot}v_2| = 0.0046/0.041/0.376$（$\propto a^2$，**正交性先丢**）；
$a=2/3$ 时临界方向数变 **28/92**（**计数也崩**）。

**后果**：
1. `lambda_min` 应改名 `kappa_min` 并必须带 $a$；`lambda` 与
   `imaginary_frequency_cm1` 只在 `method == "hessian"` 时允许出现。
2. $\nu_a$ 是 index-like 嵌套计数，依赖嵌套顺序，**无交错定理** —— 所以路线 A 里唯一严格的
   $\nu$ 下界不传到有限 $a$。
3. **底方向随 $a$ 转动，实测**：$\theta^* = -0.013°/-0.054°/-0.379°/-2.59°$
   对 $a = 10^{-2}/2{\times}10^{-2}/5{\times}10^{-2}/10^{-1}$（按 $a^2$ 增长），
   $a \ge 0.2$ 已跑出 $\pm28.6°$ 扫描窗。**2.6d 的热幅度 $a_{\rm rms} \approx 0.30$ 上，
   底方向已不是 Hessian 特征向量。**
4. **白拿的便利**：$\partial_v \kappa_a^E(v) = 2D_a(v)$（解析验证残差 $2\times10^{-10}$）。
   一对 $\pm av$ 评估同时给目标与梯度，所以球面最小化每步 **2 次力评估，与一次 $Hv$ 同价**，
   不建矩阵、不需要谱；$a \to 0$ 退化为经典 dimer。

**架构结论**：Hessian 是 validation limit 而不是 primary observable，理由不是"小体系便宜"，
而是 **$a \to 0$ 是唯一有谱定理的地方**。完整论证见
[PRRS_CURVATURE_ROUTES.md](PRRS_CURVATURE_ROUTES.md) §E。

### 2.6g 转动的"误差地板"几乎全是残余梯度 —— 2.6c 的两个数字要重新解读

$E(R(s)x) = E(x)$ 对一切 $s$ 成立，两次求导给出

$$\boxed{\dot x^\top H \dot x = \omega^2\,(\mathbf g \cdot x^\perp)}$$

所以**转动的 Rayleigh 商正比于梯度，只在驻点才为零**。平动没有这一项：能量只依赖坐标差，
恒等式在任何几何上都精确。

**解析路径实测**（乙醇，`docs/experiments/floor_split.py`）：

| 几何 | $f_{\max}$ (eV/Å) | 转动 measured | 梯度预言 | **残差** |
| --- | --- | --- | --- | --- |
| P0 输入原样 | $6.59\times10^{-1}$ | $7.75\times10^{-2}$ | $7.75\times10^{-2}$ | $8.0\times10^{-16}$ |
| 同上 | | $1.535\times10^{-2}$ | $1.535\times10^{-2}$ | $3.3\times10^{-17}$ |
| 同上 | | $1.773\times10^{-2}$ | $1.773\times10^{-2}$ | $-3.8\times10^{-16}$ |
| 松弛到 1e-4 | $7.26\times10^{-5}$ | $-2.516\times10^{-5}$ | $-2.516\times10^{-5}$ | $1.4\times10^{-16}$ |
| 同上 | | $-4.519\times10^{-6}$ | $-4.519\times10^{-6}$ | $-9.2\times10^{-17}$ |
| 同上 | | $-4.354\times10^{-6}$ | $-4.354\times10^{-6}$ | $-2.9\times10^{-16}$ |

**梯度项把转动商解释到机器精度。**所以：

1. **2.6c 记的 $2.5\times10^{-5}$ / $2.9\times10^{-6}$ 不是 MACE 的转动不变性误差**，而是那两个
   几何上残余梯度的表现。解析 Hessian 的真实刚体误差是 $\sim10^{-16}$ —— e3nn 等变性给出的
   严格不变性，与之一致。**结论方向不变**（解析路径干净），但数字的归因错了，量级差 11 位。
2. **原样的 D1 门是错的**，而且被 P0 跑抓了个正着：源结构松弛到 `quench_fmax_eV_A` 后
   $f_{\max} \approx 3.6\times10^{-3}$，转动商 $1.15\times10^{-3}$ 超过阈值 $5\times10^{-4}$ →
   在**源结构上**就失败关闭。原因不是曲率有噪声，是判据把梯度污染当成了误差。
3. 修法：**减掉梯度预言，只对残差设门**（§8.1 已改）。

**FD 路径上的交叉验证**（双井 C–H 对，`docs/experiments/floor_split_cpu.py`）：

| 几何 | $f_{\max}$ | measured | 梯度预言 | 残差 |
| --- | --- | --- | --- | --- |
| $r = 1.5$（非驻点） | 1.0000 | $7.168\times10^{-1}$ | $7.169\times10^{-1}$ | $-4.25\times10^{-5}$ |
| $r = 2.0$（非驻点） | 0.7901 | $-4.249\times10^{-1}$ | $-4.248\times10^{-1}$ | $-3.45\times10^{-5}$ |
| $r = 1.2$（极小） | 0.0000 | $3.318\times10^{-4}$ | $-1.1\times10^{-15}$ | $3.318\times10^{-4}$ |
| $r = 1.8$（势垒顶） | 0.0000 | $-7.375\times10^{-5}$ | $0$ | $-7.375\times10^{-5}$ |

FD 路径上残差不为零，且在驻点上就是 2.6f 的弦代弧误差。**所以减掉梯度项之后，残差恰好
分离出两件事：解析路径 $\to 10^{-16}$（模型严格不变），FD 路径 $\to$ 2.6f 的闭式。**
这才是 D1 本来要测的量。

### 2.6f FD Hessian 的误差地板有闭式，且乙醇尺度上超过本征值容差

把原子横向位移 $\delta$ 只在**二阶**改变键长（弦代弧），所以有限差分的横向元不是零而是
$\approx E''\delta^2/(2r^2)$。这不是模型误差 —— 双井势严格不变，这个地板**完全是微分方法的**。

实测（双井 $r=1.2$ 极小点，`docs/experiments/floor_mechanism.py`）：

| step (Å) | fd 地板 | $/\mathrm{step}^2$ | 闭式 $E''\,\mathrm{step}^2/(r^2 m)$ |
| --- | --- | --- | --- |
| 0.005 | $1.285\times10^{-5}$ | 0.5139 | $1.285\times10^{-5}$ |
| 0.010 | $5.139\times10^{-5}$ | 0.5139 | $5.139\times10^{-5}$ |
| 0.020 | $2.054\times10^{-4}$ | 0.5136 | $2.056\times10^{-4}$ |
| 0.040 | $8.202\times10^{-4}$ | 0.5126 | $8.223\times10^{-4}$ |

闭式吻合到 3 位有效数字，按 $\mathrm{step}^2$ 精确标度；解析路径上同一个地板 $< 10^{-12}$。

**乙醇尺度代入**（$E'' \approx 50$ eV/Å²、$r \approx 1$ Å、H 原子 $m = 1$、step $= 0.01$）：

$$\text{地板} \approx 2.5\times10^{-3}\ \mathrm{eV/(Å^2\,amu)} \;>\; \texttt{minimum\_check\_eigenvalue\_tol} = 10^{-3}$$

**这就是缺口 9 的机制**：不是"步长调得不好"，而是笛卡尔有限差分在这个步长上必然有的弦代弧
误差。两条出路只有一条便宜：走解析 Hessian（2.6b/2.6c）。调小 step 会按 $\mathrm{step}^2$
降低这一项，但按 $1/\mathrm{step}$ 放大相消误差，所以不能无限调。

**注意这条对旧结果的意义**：§5 的乙醇基准用的是 FD 路径，所以其中**接近零的模**
（§3.5 三级判据里的 `free` 边界）此前处在自身噪声之下。已实现的 D1 门现在会把这种情况
失败关闭，而不是静默通过。

### 2.6d 有限幅度响应就是对 PES 的 sinc² 低通滤波（已定量确认）

沿方向 $v$ 取 PES 的傅里叶分量 $e^{ikx}$，对称差分相对精确二阶导的传递函数是

$$T_a^{E}(k) = \operatorname{sinc}^2\!\left(\tfrac{ka}{2}\right) \quad\text{(能量形式)},
\qquad T_a^{F}(k) = \operatorname{sinc}(ka) \quad\text{(力形式)}$$

即幅度 $a$ 就是分辨率（截断波数约 $2/a$）。**只有能量形式是非负滤波**；力形式
$\operatorname{sinc}(ka) < 0$ 当 $ka \in (\pi, 2\pi)$，会把波长在 $a$ 与 $2a$ 之间的分量反号。

领头阶给出可检验的预言 $(\kappa_a^F - \kappa_0)/(\kappa_a^E - \kappa_0) \to 2$。
**实测在三个模式、跨一个半数量级幅度上都是 2.00**（$a = 10^{-2}$ 到 $3\times10^{-1}$：
不稳定模 2.004→1.857，最低正模 1.988→1.967，次低正模 1.923→1.956）。滤波图像因此是
被数据定量确认的，不是解释性说法。

三个直接后果：

1. **多尺度估计量应用能量形式**，不是力形式（非负性压倒条件数，因为 D1–D5 表明数值噪声
   不是本体系的问题）。
2. **滤波扣不掉平滑偏差**（作用在频率上，偏差在幅度上），也扣不掉真实的高频物理
   —— 它是**带宽选择**，不是去噪器。
3. **热幅度 $a_{\rm rms} = \sqrt{k_BT/\kappa}$ 落在偏差极大的区间**：300 K 下软模
   $a_{\rm rms} \approx 0.30$，该处最低正模 $\kappa_a$ 偏 115%，**鞍点不稳定模的
   $\kappa_a$ 连符号都不对**。对速率有意义的曲率是热可及幅度上的，$a \to 0$ 不是正确目标
   —— 这个结论不依赖 MLP 可靠性。

估计量与判据（带内回归、$\sigma$ 显著性检验、`sign_stable`、D1 门）见
[PRRS_CURVATURE_ROUTES.md](PRRS_CURVATURE_ROUTES.md) 路线 C。

### 2.6b MACE 的力是 autograd，所以解析二阶导本来就在

`mace/modules/utils.py:23` 的 `compute_forces` 是
`torch.autograd.grad(E, positions, create_graph=training)`，且 mace 自带
`compute_hessians_vmap`（`utils.py:113`）与 `MACECalculator.get_hessian`
（`calculators/mace.py:770`）。ASE Calculator 路径以 `training=self.use_compile`（默认
`False`）调模型，所以**图被丢掉的是接口，不是模型**。

乙醇 9 原子、float64、CUDA、一个 `rattle(0.05)` 后的非驻点（$f_{\max} = 6.61$ eV/Å）实测
（脚本 `docs/experiments/hvp_autograd_probe.py`，输出同名 `.out`）：

```
双反向 Hv 与 get_hessian 一致       max|H@v - Hv_exact| = 4.4e-15 .. 7.1e-15
对称性                              v0.H v1 - v1.H v0    = 2.665e-15
与 ASE 路径同一物理量               dE = 0.000e+00,  max|dF| = 1.110e-15 eV/A
FD 截断误差（|Hv| ~ 20-40 eV/A^2）  eps=1e-3: 1.4e-5   eps=1e-2: 1.4e-3   eps=5e-2: 3.5e-2
一次力评估                          35.4 ms
一次解析 Hv                         69.5 ms  = 1.96 次力评估
一次 FD Hv（2 次力）                70.8 ms  = 2.00 次力评估
完整解析 Hessian（get_hessian）     946 ms   = 26.7 次力评估
完整 FD Hessian（6N = 54）          1912 ms  = 54 次力评估
峰值显存                            117.9 MiB
```

三条结论：

1. **解析 $Hv$ 与 FD $Hv$ 成本相同（1.96 vs 2.00），精度差 12 个数量级。**
2. **完整解析 Hessian 比 `hessian_spectrum` 的 FD 版快 2 倍**，且是机器精度。
   `hessian_spectrum`（`runner.py:105`）在解析导数已可得的情况下做数值微分。
3. **当前 `minimum_check_step_A = 0.01` 下的 FD 误差（$1.4\text{–}7.8 \times 10^{-3}$ eV/Å²）
   大于 `minimum_check_eigenvalue_tol = 1e-3`。** 质量加权压小重原子分量但不压 H。
   §3.5 三级判据里 `free`（$|k| < 0.01$）与 `negative curvature` 的分界因此处在 FD 噪声里。
   这是非驻点上的单点测量，$\epsilon^2$ 标度普适、系数不普适，**必须在乙醇三个 TS 上复测**。

4. **`get_hessian` 不做单位换算**。`calculate` 在 `mace.py:657` 起对
   energy/forces/stress 逐项乘 `energy_units_to_eV` / `length_units_to_A`；
   `get_hessian`（`mace.py:770`）**一项都不乘**。MACE-OFF 两个因子都是默认 1.0，所以现在
   静默正确 —— 这正是 §2.4 那类"默认值恰好对，所以没人发现"的陷阱。切换到解析 Hessian
   时必须自己乘，并按 §3.7 用执行证实（拿一个已知谱交叉校验），不能信参数。

路线与优先级见 [PRRS_CURVATURE_ROUTES.md](PRRS_CURVATURE_ROUTES.md) §A.13。

### 2.6 势能面刚度实测（乙醇，收敛到 fmax 0.00359 eV/Å）

| 坐标 | 刚度 |
| --- | --- |
| OH 扭转 | 0.212 eV/rad² |
| 甲基转子 | 0.625 eV/rad² |
| C–H 伸缩 | 33.2 eV/Å² |
| O–H 伸缩 | 52.7 eV/Å² |
| C–C 伸缩 | 35.3 eV/Å² |

扭转比键伸缩软约 150 倍。这个量级差是第 3.5 节判据存在的原因。

---

## 3. 已确立的设计原则

每条都由实测或缺陷驱动，不是先验设计。

### 3.1 扰动基必须张开主要内坐标子空间

$$\mathrm{span}(\mathcal P) \supseteq \{q_{\rm stretch},\, q_{\rm bend},\, q_{\rm torsion}\}$$

**证据**：乙醇首轮 24 次试验 6/6 全弹性，而 anti / gauche⁺ / gauche⁻ 三个真实盆地就在 4 meV 处。原三家族（键拉伸、非键接近、pair kick）**全是成对径向扰动**，在 100 fs 内耦合不到扭转。

按家族统计"移动到另一 microstate"的比率：`torsion 4/10、torsion_kick 0/6、stretch 0/4、kick 0/4`。

**推论（尚未实现）**：位移版对浅势垒优于动量版 —— 对 47 meV 的势垒，kick 注入的能量让体系在自由段游走，淬火时回到原处或落在鞍点上。

### 3.2 target-coordinate fidelity ≠ geometric fidelity

$$Q_{\rm perturb} = (Q_{\rm target},\, Q_{\rm collateral})$$

**证据**：2.09 rad 的扭转位移把 O–H 从 0.957 Å 拉到 **2.346 Å**，而请求二面角精确达标。当时所有检查都通过。

现在 `apply()` 返回 `target_error` 与 `collateral` 两组量。门控取舍：

| 量 | 处理 | 理由 |
| --- | --- | --- |
| 共价键长（非目标） | 硬门 0.05 Å | 键是刚的，任何探针都不该拉伸 |
| 原子重叠 | 硬门 | |
| 碎片数（排除目标坐标自身造成的） | 硬门 | stretch 拉断目标键是目的，不是损伤 |
| 非目标键角、非目标扭转 | 只测量 | 合法耦合：同顶点角度有求和约束 |

### 3.3 大扰动 = 一串局部有效的小扰动

$$R_0 \to R_1 \to \cdots \to R_n$$

沿起点切线一口气走到终点等于用圆弧在起点的切线近似整个圆。位移增量上限 `MAX_STEP = 0.05`（Å 或 rad），每步重算梯度。

**并且优先用精确实现**：桥键扭转走刚体 fragment 旋转 —— j、k 在轴上不动，故**所有键长和所有键角精确保持**，只有该扭转改变。实测在 ±2.09 rad 下 `target_error ≤ 4e-16`、键长变化 ≤ 2.2e-16 Å、键角变化 ≤ 2.2e-16 rad。

rigid 与 continuation 的分界**恰好是 bridge 判据**：桥键 ⟺ 存在刚体分解。环上的键不是桥，走 continuation（苯环实测键长变化 < 0.05 Å）。

### 3.4 在分子对称性的商空间上搜索

$$\phi \in S^1/C_n,\qquad n = \mathrm{lcm}(n_{\rm left},\, n_{\rm right})$$

一端有 $m$ 个等价取代基则贡献 $C_m$；两端各 $m$、$n$ 时由 $2\pi/m$ 与 $2\pi/n$ 生成的加法子群是 $2\pi/\mathrm{lcm}(m,n)\cdot\mathbb Z$。乙醇：OH 扭转 $n=1$，甲基 $n=3$；乙烷两端甲基 $n=\mathrm{lcm}(3,3)=3$。

**证据**：甲基扭转吃掉 15/24 次试验且不可能有发现（转出来全是对称副本）。

判据不是"幅度小于基本域宽度"，而是**到最近对称副本的距离**：$\min_m |a - m\cdot 2\pi/n| > $ margin。2.09 rad 是 119.75°，距 $2\pi/3$ 只有 0.25°，形式上在域内、物理上是恒等操作。

**对称转子降级但保留**。它对盆地发现无用，但实测贡献了一个真实一阶鞍点（甲基旋转 TS，+130.5 meV，虚频 −247.5 cm⁻¹）。所以：

$$U_{\rm basin} \approx 0 \quad\text{而}\quad U_{\rm TS} > 0$$

### 3.5 mode-resolved 收敛判据

$$\mathrm{converged} = \Big(f_{\max} < f_{\rm base}\Big) \wedge \Big[\forall v \in \mathcal S_{\rm soft}:\ |g_v|/k_v < q_{\rm tol}\Big]$$

固定的力容差只控制"力有多小"，不控制"离驻点多远"：$|\delta q| \approx |g_q|/k_q$。

三级：

```
k > 50 eV/rad²        stiff               全局 fmax 已管住
0.01 ≤ k ≤ 50         soft                要求 |g/k| < 1°，沿自身路径 Newton 步
|k| < 0.01            free                记录，不要求收敛（值不携带信息）
k < -0.01             negative curvature  交给极小检验，polish 不碰
```

torsion 直接用内坐标单位的 $k_\phi$（刚体旋转给出精确路径，能量中心差分即可，无链式法则歧义，2 次能量评估）。

**实测对照**（解析链，$k = 0.225$ eV/rad²，两者都满足 fmax ≤ 0.02）：

| | fmax (eV/Å) | 预测 $\|g/k\|$ | 实际 φ 偏移 |
| --- | --- | --- | --- |
| 仅 fmax | 0.01994 | 3.930° | −3.877° |
| + polish | 0.00692 | 0.053° | 0.053° |

预测与实际吻合，直接验证了线性响应关系。代价 6 次额外能量评估；硬坐标原地不动（键长偏差 0.34 mÅ、键角 0.007°）。

**真实体系效果**：乙醇 gauche 角度的批次间散布从 ±4° 收到 ±0.15°，且 300.15 + 59.85 = 360.00 精确（对映体必须满足的镜像关系）。polish 开销约每试验 2%。

### 3.6 商掉冗余自由度

这是全套设计中唯一反复出现的抽象，已在六处独立救场：

1. `chemical_key` 的图部分用自同构不变哈希 —— 否则转移三个等价甲基 H 中的哪一个会产生三个"产物"
2. 构象去重用对称感知 RMSD —— 否则甲基转 120° 会被判为新构象（plain RMSD 0.914 Å vs symmetric 0.006 Å）
3. 角度方向枚举按正则色去重 —— 三个等价 H–C–H 角不占三份方向额度
4. 对称转子降级 + 基本域限制
5. reservoir 多样性度量排除转子扭转 —— 否则对称副本会显得"多样"
6. 立体化学签名只取**构型**手性（四面体奇偶性 + 锁定键 E/Z），绝不含可旋转单键的扭转符号 —— 否则 gauche⁺/gauche⁻ 会变成两个化学节点

### 3.7 存在/生效必须由执行证实，不能由记录证实

`importlib.metadata` 对 conda 装的 openmm-ml/openmm-torch/nnpops 假阴性（但 openmm/numpy 有 metadata，所以这不是 conda 的普遍规律 —— 正因如此 metadata 不能当存在性判据）。已改为 `find_spec` 探模块存在性，结果诚实标为 `"installed; no version metadata"` 而非 `None`。

同理 `create_mace_context` 不信 `removeCMMotion=False` 这个**参数**，而是建完 System 后**清点 force 列表**。

---

## 4. 两层网络与双预算

```
Chemical network  G_chem = (V_chem, E_rxn)
    │
    └── chemical node C_i        由 chemical_key 判定
            │
            └── conformer reservoir  {R_i1, ..., R_iK}，K = conformer_max_per_node
```

**准入裁决权归 `chemical_key`，不归分类标签**。等价 H 转移会得到 `reactive` 标签但相同的 key，此时绝不开新节点。

预算分离：

```yaml
chemical:   max_nodes: 16   trials_per_node: 96
conformer:  max_per_node: 8  trials_per_microstate: 32
```

同一节点内新增构象**不扣 chemical 预算**。构象仍是合法扰动起点（反应可及性依赖构象：$P(B|A,\text{conformer})$ 而非 $P(B|A)$）。

Reservoir 评分 $S = w_E S_E + w_D S_{\rm diversity} + w_R S_{\rm response}$：

- $S_E = \exp(-\Delta E/kT_{\rm eff})$，相对节点内最低能构象
- $S_{\rm diversity}$：可旋转扭转上的 RMS 圆周距离（不是笛卡尔 RMSD；角度是周期的）
- $S_{\rm response}$：未搜索者给乐观先验 1.0，已搜索者取"独有响应占比"。这是**在线选择**而非一次性聚类，因为响应必须先花预算才知道

### 4.1 三种响应

| 响应 | 归属 |
| --- | --- |
| 回到原极小 | 弹性证据 |
| 松弛到另一个极小 | 盆地发现（同 key → microstate；异 key → chemical node） |
| 落在 $n_- = 1$ 的驻点 | TS pool |
| $n_- > 1$ | 拒绝，可诊断 |

### 4.2 鞍点下降定端点

对每个 $n_-=1$ 候选，沿虚频特征向量 $\pm\epsilon$ 各淬火一次。**TS 连接是无向结构性证据，不自动生成有向边**（计划书 §6.2：已观测 A→B 不自动产生 B→A）。它进 `ts_connected_channels`，并在 publish 时交叉引用为有向边的 `ts_support`（在发现当时标注会漏掉后出现的边）。

端点同属一个 chemical node 时记为该节点的 `conformer_transitions`，不进反应网络。

**副产品**：polish 的分级数据直接指出每个 TS 的反应坐标是哪个内坐标（负曲率所在），无需投影特征向量，且是化学可读的名字。

---

## 5. 乙醇 P0 基准

验收判据在跑前冻结，不边跑边改。最近一次完整运行：**13/13 通过**，24 次试验 512 s。

```
1 chemical node  key=ec63499308a3  frags=['C2H6O']  reservoir 3/8
    m0000  C-C-O-H = 180.00°   E = -4221.60312 eV   anti
    m0001  C-C-O-H = 300.15°   E = -4221.59941 eV   gauche−
    m0002  C-C-O-H =  59.60°   E = -4221.59952 eV   gauche+
reactions: []                      ← 构象之间没有虚假反应边
ts_connected_channels: []          ← 乙醇只有一种化学，故不应有化学通道
TS pool: 3，全部两端解析成功（6/6 端点定名）
    φ≈0°/359°   +47.1 meV   -290.1 cm⁻¹   OH 扭转负曲率    → m0001 ↔ m0002
    φ=180.03°  +130.5 meV   -247.5 cm⁻¹   甲基扭转负曲率  → m0000 ↔ m0000（self_connected）
    ！这两个虚频由 FD Hessian 得到，各带约 1 cm⁻¹ 截断误差。解析值 -291.03 / -246.86，见 8.7。
budgets: chemical_nodes_used = 1,  microstates_total = 3
按家族的"移动到另一 microstate"率： torsion 4/10, torsion_kick 0/6, stretch 0/4, kick 0/4
```

两条核心不变式：

$$\text{不同极小} \neq \text{不同 chemical state}$$
$$\text{conformer discovery} \not\Rightarrow \text{chemical budget consumption}$$

加上 torsion coverage、finite-amplitude geometry fidelity、symmetry reduction、mode-resolved
convergence，以及鞍点端点解析。

### 5.1 两条事前写下并命中的预测

- **syn 势垒连接 gauche⁺ ↔ gauche⁻**（`self_connected=False`）
- **甲基势垒自连接**（`self_connected=True`）—— 它连接一个结构与它自己的对称副本，
  对称感知去重正确判定两端为同一 microstate。这是对第 3.4 / 3.6 节对称商逻辑的独立确认。

### 5.2 数值精度的诚实边界

gauche⁺/gauche⁻ 是对映体，理论上应精确镜像。实测 59.60° 与 300.15° 的镜像关系差 **0.25°**，
能量差 0.11 meV。这在声明的 1° 容差之内，是 `torsion_tolerance_rad` 设定的边界而非缺陷；
上一轮同一配置恰好落到 360.00° 精确。收紧该容差会收紧此关系。

---

## 6. 与计划书的分歧（以实测为准）

| 计划书 | 实测结论 |
| --- | --- |
| §2「优先 OpenMM 原生轨迹……避免每步经 Python 往返搬运 GPU 力数据」 | 对 MACE 不成立。openmm-ml 的 MACE 力本身就是 `PythonForce`。**ASE 为主执行后端，OpenMM 为独立数值交叉校验**，其长期价值在 `createMixedSystem` 的 ML/MM 分区（QM/MM 接口），不在吞吐 |
| §3 第 9 条「不能因键图变化就一律重建 Context；需 P0/P1 验证」 | 已验证：不需要重建，dE = dF = 0 |
| §4「集体扰动」列为后续增加 | **提为 MVP 必需项**。见 3.1 |
| §10「`src/prrs/` 为未验收草稿」 | 已过时。当前 82 条测试通过，含两个已修实质 bug（势垒顶被当产物盆地、manifest 漏记 conda 包） |

尚未采纳但已识别的计划书要点：force pulse 的外功验收、ensemble 不确定度（`Committee` 已在 `calculators.py`，但 MACE-OFF 单模型下不可用）、周期体系。

---

## 7. 已知缺口

1. ~~**自由坐标未商掉**~~ —— **已关闭（§8.8.2）**。`free_aligned_symmetric_rmsd` 在比较前沿自由坐标对齐，自由度由淬火报告的 `tier == "free"` 提供并随几何走；零力评估。测试同时钉住"不过度归并"（仅软的坐标不商）。原文：若某扭转真自由（$k \to 0$），只在该坐标上不同的结构由零势垒路径相连、应属同一盆地，但 `symmetric_rmsd` 会判为不同 microstate，于是自由转子能生成无穷多假构象。乙醇不触发（甲基 $k = 0.625$）。 理论出处见 [PRRS_CURVATURE_ROUTES.md](PRRS_CURVATURE_ROUTES.md) §B.6：自由坐标就是 Morse–Bott 分解里的中性子空间 $\mathcal N$，TS 判据只要求 $\dim\mathcal U = 1$，不要求 $\dim\mathcal N = 0$。
2. **`free` 与 `stiff` 两级只有解析测试覆盖** —— `free` 现在有一条端到端走通的解析测试（§8.8.2 的四原子链，势垒设 0），但**真实体系运行仍未触发**：乙醇的所有扭转都落在 soft 带（0.19–0.69 eV/rad²）。`stiff` 一级依旧只有解析覆盖。
3. **排序仍是启发式**。应改为 $U(P) = (U_{\rm basin}, U_{\rm reaction}, U_{\rm TS})$ 三分量。$U_{\rm TS}$ 有物理代理 $\sim k_\phi \cdot (\text{基本域})^2$；$U_{\rm basin}$ 可用方向刚度（`direction_stiffness` 已实现，2 次力评估）。
4. **大体系没有 $n_-$ 判定路径**。完整 Hessian 是 $6N$ 次力评估；projected Lanczos 只需数十个 $Hv$，crossover 约在 15 原子。`hessian_vector_product` 原语已抽出，Lanczos 未实现。 完整路线（三档部署、认证条件的残差界、诊断字段改名）见 [PRRS_CURVATURE_ROUTES.md](PRRS_CURVATURE_ROUTES.md) 路线 A；结论是核心层只需 $\nu = \#\{\lambda_i < 0\}$，虚频降为可选 reporting quantity。 **但 2.6b 的实测把优先级倒过来了**：解析完整 Hessian 在 $N \lesssim 30$ 上比 Lanczos 认证便宜 3–4 倍且严格，Lanczos 的适用尺寸目前没有基准体系。第一件该做的事是让 `hessian_spectrum` 改用 `calc.get_hessian()`，不是实现 Lanczos。
5. **mode-resolved 判据只覆盖 torsion**。其他软模需要 projected-Hessian 路线。任意方向的曲率现在可由 `_curvature_along`（2 次力评估）或 `response_curvature` 取得（§8.3、§8.8.1），但收敛判据本身尚未推广到非扭转坐标。
6. **`S_response` 只有乐观先验那一半在起作用**。乙醇三个构象里两个 `trials=0`，在线更新未被真实验证。
7. **`pulse` 家族未验收**。外功与撤力切换的连续性未测。
8. **真实化学基准：P1 已把核心主张证实了一半，P2 仍未跑**。乙醇 P0 只验证构象发现与架构；**P1 丙二醛在 2026-09-01 第一次给出了 bond-changing 反应的证据** —— `reaction_channels` 里的 `rc0001`：`class=degenerate_reaction`、`self_loop=true`、`broken=[[4,5]]`（O2–H）、`formed=[[0,5]]`（O1–H），由 `stretch` 在 0.5 与 0.8 Å 命中，event key 与 preflight 预测逐字符相同。见 [P1_MALONALDEHYDE_ACCEPTANCE.md](P1_MALONALDEHYDE_ACCEPTANCE.md)。

   **剩余部分**：P1 第一轮 7/12（A/C 的失败查出了真实的 identity 缺陷，见 §8.12/§8.13）；TS 一个都没找到（缺口 9）；**P2（不对称 β-二羰基，验证 $k_A \neq k_B$ 的网络扩张）尚未运行** —— 它的验收对象是 `chemical_key` 本身，所以要等 identity 层冻结之后再跑。SN2 仍被电荷阻塞（§2.5）。

   **势垒数值的两个来源必须分开**：MACE-OFF24 自己测得的丙二醛质子转移势垒是 **419.14 meV**（`runs/p1_malonaldehyde/p1_ts.extxyz`，$\lambda_1 = -37.323$、$-3186\ \mathrm{cm^{-1}}$）。验收文档里出现过的"约 0.17 eV"是**外部印象值、未经独立验证**，不是本项目的已知量，两者相差 2.5 倍这件事本身正是 §10.3 第 1 条的实例。

9. **TS 发现依赖"淬火停在鞍点"，对刚性反应坐标失效**。P1 实测：60 次试验 0 个 TS candidate，而鞍点确实存在（$\lambda_1 = -37.3$、$-3186\ \mathrm{cm^{-1}}$、势垒 419 meV）。三条原因与四项待办见 [P1_MALONALDEHYDE_ACCEPTANCE.md](P1_MALONALDEHYDE_ACCEPTANCE.md)"追查"一节。通用解法是在已知端点之间搜 TS（`lowest_response_direction` + 解析 Hessian），让 TS 发现不再依赖 probe 幅度。

10. **非键接近家族在密堆分子上交付不了**。P1 实测 12/60 试验被 `collateral_bond` 挡掉，全部是 `compress` 作用在非键对上：continuation 交付把共价键拖动 0.08–0.41 Å，预算是 0.05 Å。门是对的，交付满足不了它。

11. **FD 步长与本征值容差互相矛盾** —— 机制已定（2.6f 的闭式），D1 门已落地（§8.1），解析路径已默认（§8.2）。剩余：`minimum_check_step_A` 对非 MACE 后端仍是同一个问题。见 2.6b 第 3 条：`minimum_check_step_A = 0.01` 下 FD 的 $Hv$ 误差大于 `minimum_check_eigenvalue_tol = 1e-3`，于是接近零的模式（§3.5 的 `free` 边界、以及弱负曲率）处在数值噪声里。改用解析 Hessian 直接消除这一项（2.6c 表明该解析 Hessian 在 $10^{-5}$ 水平可信，前提成立）。另加 D1 门（6 个平凡模的 Rayleigh 商，零额外力评估）作为永久的无参考质量监控。

---

## 8. 曲率层改动（已落地）

基线：`105 passed, 1 skipped`（原 `84 passed, 1 skipped`），新增 `tests/test_curvature.py` 21 条。
第一批是 8.1–8.7（`96 passed`），第二批是 8.8。
路线依据见 [PRRS_CURVATURE_ROUTES.md](PRRS_CURVATURE_ROUTES.md) §A.13.6 / §C.10。

### 8.1 D1 门：平凡模误差地板（`runner.py: trivial_mode_floor`）

6 个（线性分子 5 个）刚体模的 Rayleigh 商。对任何不变势严格为零，所以它们的实测大小就是
整个曲率计算的误差地板，**零额外力评估**（刚体基本来就要构造来投影掉）。
`confirm_minimum` 现在在

$$\max_i \big|\lambda_i^{\rm trivial} - \omega_i^2(\mathbf g\cdot x^\perp_i)\big|
> \texttt{trivial\_floor\_fraction} \times \texttt{minimum\_check\_eigenvalue\_tol}$$

时**失败关闭**，`reason = "trivial_mode_floor_above_tolerance"`、`saddle_order = None`。
减掉的那一项是 2.6g 的梯度预言：不减它，任何按普通力容差收敛的几何都会被误判
（P0 跑实测到了这一点）。每个模式都记 `kind` / `measured` / `gradient_expectation` /
`residual`，门只看 `residual`。额外代价是一次梯度，而调用点刚淬火完通常已缓存。

### 8.2 解析 Hessian（`hessian_source = "auto" | "analytic" | "fd"`）

后端自带二阶导时直接用（MACE 的力本身就是 autograd，见 2.6b）。**按 §3.7 用执行证实**：
沿一个固定方向的 $Hv$ 必须与有限差分一致到 `analytic_hessian_check`（默认 1e-2 相对），
否则 `"analytic"` 抛 `GateRejected`、`"auto"` 退回 FD 并把原因记进 provenance。
这一层专门防 2.6b 第 4 条的单位陷阱（`get_hessian` 不做换算，而 MACE-OFF 两个因子恰好是 1.0
所以静默正确）。测试里用 96.485 倍的假 Hessian 钉住了这条。

### 8.3 能量形式带内回归（`runner.py: response_curvature`）

$\kappa_a = \kappa_0 + c\,a^2$ 对 $a^2$ 的最小二乘，返回
`kappa` / `sigma` / `anharmonicity` / `sign_stable` / `significant` / `residual_rms` / `samples`。
用能量形式而非力形式，因为只有 $\operatorname{sinc}^2$ 非负（2.6d）。
$2K$ 次评估（只用能量，但每次都带回力，因为可靠性门需要力来判断探针是否物理）。`response_gradient` 实现 $\partial_v\kappa_a^E = 2D_a(v)$（2.6e 第 4 条）。

### 8.4 命名规则（2.6e）

TS 记录新增 `curvature_source`、`trivial_mode_floor_worst`、`response`；
`imaginary_wavenumbers_icm` 改为 `.get()`，只有 `method == "hessian"` 的路径产生它。

### 8.5 修掉一个已有 bug：质量加权方向反了

`direction_stiffness` 原来对笛卡尔方向**除以** $\sqrt m$。质量加权坐标是 $q = \sqrt m\,x$，
所以应当**乘**。同质量双原子上两者只差一个被归一化吃掉的标度，所以原测试（C₂）看不出来；
质量不等时它们是**不同的方向**。已抽出 `_mass_weighted_direction` 统一，并用 C–H 对
（质量比 > 10）钉住。

发现路径值得记：这个 bug 是在真实后端上验证 8.3 时暴露的 —— 带内回归给出
$\hat\kappa = 0.3995$ 而 $\lambda = 0.2164$，偏 85%，而 $\sigma = 3.7\times10^{-6}$ 毫无察觉。
**$\sigma$ 只度量"二次模型解释样本的程度"，不度量"方向对不对"**，这是它的边界。

### 8.6 真实后端实测（乙醇极小点，`fmax = 7.3\times10^{-5}`）

```
source=analytic   1.62 s   finite_difference_deviation = 5.19e-05
    raw trivial floor = 2.516e-05   residual after subtracting g-term ~ 1e-16   PASS
source=fd         2.92 s   force_evaluations = 54
    raw trivial floor = 1.211e-03   residual ~ 2.6f 的弦代弧误差         PASS/FAIL 由 step 决定

两条路径的 6 个最低物理本征值
    analytic  0.21642 0.29237 0.64861 2.44093 3.04340 4.05355
    fd        0.21621 0.29472 0.64778 2.43816 3.04307 4.05096

OH 扭转模上的带内回归（band = [0.01, 0.1]，5 点）
    lambda = +0.216425      kappa = +0.216429 +- 3.21e-06      c = +1.742
    kappa - lambda = +4.4e-06        residual_rms = 4.36e-06
```

三条：

1. **解析路径生效且更快**（1.62 vs 2.92 s），交叉校验偏差 $5.2\times10^{-5}$。
2. **上表的 `PASS/FAIL` 是这一批第一版的结果，已被 2.6g 修正。**第一版对**原始**地板设门，
   于是 P0 跑在源结构上就失败关闭（$1.15\times10^{-3} > 5\times10^{-4}$）——
   而那个数字是残余梯度，不是误差。现在门只看减掉梯度项后的残差：解析路径 $\sim10^{-16}$，
   FD 路径留下 2.6f 的闭式。**修正过程本身是这一批最有价值的产出**：一个真实跑把一个
   概念错误的判据顶了出来。
3. **回归的 $\sigma$ 诚实**：实际误差 $4.4\times10^{-6}$ 与 $\sigma = 3.2\times10^{-6}$ 同量级。

### 8.8 第二批：曲率探针、自由坐标商、球面下降

基线 `105 passed, 1 skipped`（第一批后是 96）。三项对应路线 A.4 / B.6 / A.5+E.4。

#### 8.8.1 A.4：`_probe_minimum` 从"位移后重淬火"改成测曲率

原来每个方向要一次完整淬火，返回值只有"回来了/没回来"。乙醇势垒顶上两个等价鞍点只抓到
一个 —— 滑走是比负数弱的信号。现在每个方向 **2 次力评估**给出 $v^\top H v$，符号直接判定。

方向来源不再是纯随机：先取分子自己的可旋转扭转的切向（软模与不稳定模都在那里，切向由
`step_torsion` 的精确刚体旋转差分得到，零能量评估），再用种子随机方向补足
`minimum_check_probes` 个。每个方向都投影掉刚体分量（平动曲率恒为零，转动曲率正比于残余
梯度，留着只会稀释信号）。

诊断里每个方向记 `source` / `indices` / `curvature_eV_A2_amu` / `negative`，外加
`negative_directions`、`lowest_curvature`。**判定的不对称性写进了 `meaning`**：
一个负值是鞍点的决定性证据，全正**不是**极小的证明。

实测钉住：双井三个几何上 `lowest_curvature` 与闭式 $E''/\mu$ 相对偏差 $<10^{-3}$；
四原子链在 $\phi = 60°$（$(V/2)(1-\cos3\phi)$ 的极大，曲率 $-4.5V$ 已知）上**由扭转切向那个
方向**测出负曲率，并报出 `indices = [0,1,2,3]` —— 反应坐标直接带化学可读的名字，不需要
投影特征向量。

#### 8.8.2 B.6：自由坐标商掉，缺口 1 关闭

刚度落到零的扭转把它两侧的结构用无势垒路径连起来，所以它们是**同一个 microstate**，
不管坐标差多远。原来 `symmetric_rmsd` 判它们不同，于是真自由转子会生成**无穷多假构象**。

落地：

- `chemistry.free_aligned_symmetric_rmsd(...)`：在 `symmetric_rmsd` 之上再对自由坐标的
  旋转取极小（坐标轮转下降，两轮，`free_alignment_samples = 24` 个角度）。**零力评估** ——
  桥键旋转是精确几何。
- `perturbations.torsion_on_bond(...)`：坐标由**中心键**识别而不是四个原子。同一键上所有
  二面角只差常数，`rotatable_torsions` 选哪个代表只依赖正则着色，所以两个同分子结构对它一致
  —— 这才使得"在一个结构上量到的刚度"能用到另一个结构的几何上。
- `network.free_bonds_of(structure)`：从 `structure.info["quench"]` 最后一轮的
  `tier == "free"` 条目读出。**这个判断是能量性的，所以只在量过它的几何上成立**，
  因此随几何走，不从图上猜。没有淬火报告的结构贡献空集，比较退回原来的对称感知 RMSD。
- `Microstate.free_bonds` 存下来；`_match_microstate` 用**两者的并集**（哪一边发现坐标是
  平的都算，平的对双方都平）。
- `runner.run_trial` 与 `follow_unstable_mode` 现在给 endpoint 盖上
  `info["quench"]`，与 `relax_source` 早就在做的一致。

**为什么不能从图上猜**：把仅仅"软"的坐标也商掉会摧毁真实构象 —— 乙醇的 anti 与 gauche
就差一个 $0.21$ eV/rad² 的扭转。测试同时钉住了这一面：链的势垒设为 0 时
`free_bonds_of` 给 `((1,2),)`、四个起始构象归并成 **1** 个 microstate；设为 0.05 eV 时给
`()`、四个起始构象落到三个真实极小上，仍是 **3** 个 microstate。

**乙醇不触发**（所有扭转 0.19–0.69 eV/rad²，全在 soft 带），所以这条对 P0 是恒等操作 ——
§8.9 的回归确认了这一点。缺口 2（`free` 一级只有解析测试）也因此部分改善：现在有一条
真实走通 `free` 分支的端到端测试。

#### 8.8.3 A.5 + E.4：球面下降，以及只报告不裁决的 persistence

`_response_and_gradient`：一对 $\pm av$ 的评估**同时**给出 $\kappa_a^E(v)$ 与它的梯度
$2D_a(v)$（门本来就随能量返回力）。所以

$$\boxed{\text{每步 2 次力评估 —— 与一次 } Hv \text{ 同价，不建矩阵、不需要谱}}$$

`lowest_response_direction`：球面投影梯度下降 + 回溯。$a \to 0$ 时退化为经典 dimer /
Rayleigh 商最小化。**收敛判据是"梯度中还有多少比例是切向的"，不是切向范数的绝对值** ——
后者带单位，一个固定容差不可能同时服务软模和硬模（与 §3.5 同一个理由）。

`response_persistence`：跨幅度追踪底方向，报每个 $a$ 上的 $\kappa$、与上一个 $a$ 的转角、
以及 `sign_stable`。**不给裁决** —— 允许转多少度这个阈值在乙醇上定不了（所有扭转都软，
没有东西去检验它），所以数字报出来，阈值留给有真实基准的人。

实测钉住：四原子链在 $\phi = 0°$（极小）和 $60°$（鞍点）上，从随机方向出发的下降都收敛到
投影 Hessian 的**最低本征值**（相对偏差 $< 2\%$），初始转角 $> 5°$（确实动了），
而全程没有构造或对角化任何矩阵。成本核对：实际物理评估数 $=$ 报告的
`force_evaluations` $+ 1$（参考能量）。

#### 8.8.4 这一批仍未做

- **$\nu_{\rm band}$ 的嵌套计数**：只有单方向，还没在正交补里逐层往下。E.1.4 已说明有限
  $a$ 下这个计数依赖嵌套顺序，所以做它之前要先决定报什么。
- **`lowest_response_direction` 还没接进搜索**：`confirm_minimum` 仍用 $a\to0$ 的解析
  Hessian 取不稳定模。接线要等 persistence 阈值有依据。
- **$\operatorname{PT}$（B.3 的梯度输运）**：仍是投影。
- **外部 DFT 参照**：平滑偏差依旧不可见。整条线最大空缺，位置没变。

### 8.9 乙醇 P0 回归：结构量全部复现，虚频被修正

同一输入、同一配置（`runs/e2e_ethanol/config.json` 的 `config` 段）、同一种子重跑 24 次试验。

| 量 | 基准（FD） | 本次（解析 + D1 门 + 带内回归） |
| --- | --- | --- |
| status / trials | completed / 24 | **完全相同** |
| chemical key / fragments | `ec63499308a3` / `['C2H6O']` | **完全相同** |
| microstates / reservoir | 3 / 3-of-8 | **完全相同** |
| microstate 能量 | $-4221.60312$ / $-4221.59941$ / $-4221.59952$ | **逐位相同** |
| reactions / ts_connected_channels | 0 / 0 | **完全相同** |
| TS 数 / 能量 | 3 / $-4221.55634$ ×2, $-4221.47261$ | **逐位相同** |
| TS 连接 | gauche⁺↔gauche⁻ ×2；甲基 self_connected | **完全相同** |
| conformer_transitions | 3 | **完全相同** |

**两条核心不变式与 §5.1 的两条事前预测全部保持。**唯一变的是曲率数值：

| TS | FD $\mathrm{cm^{-1}}$ | 解析 $\mathrm{cm^{-1}}$ | $\lambda_{\rm fd}$ | $\lambda_{\rm analytic}$ | $\Delta\lambda$ |
| --- | --- | --- | --- | --- | --- |
| syn（OH 扭转） | $-290.08$ | $-291.03$ | $-0.309439$ | $-0.311469$ | $-2.03\times10^{-3}$ |
| 甲基转子 | $-247.55$ | $-246.86$ | $-0.225354$ | $-0.224100$ | $+1.25\times10^{-3}$ |

两个偏移**方向相反**，所以是截断误差而非系统偏置；量级 $1.3\text{–}2.0\times10^{-3}$，与 2.6f 闭式
给的 $\sim2.5\times10^{-3}$ 一致。**结论：§5 记录的 $-290.1$ / $-247.5\ \mathrm{cm^{-1}}$ 各带约
1 $\mathrm{cm^{-1}}$ 的有限差分误差；解析值是 $-291.03$ / $-246.86$。**

新字段实测：

```
              floor 残差   kappa            sigma     c        sign_stable  significant
ts0000/1      1.2e-15      -0.311287    1.4e-04    +14.28    True         True
              2.5e-16
ts0002        6.1e-16      -0.224094    3.6e-06     +1.77    True         True
```

三条读法：

1. **D1 门在解析路径上以 12 个数量级的余量通过**（残差 $10^{-16}$–$10^{-15}$ vs 阈值
   $5\times10^{-4}$）。它不是噪声计，是防 bug 与防非不变模型的门。
2. **带内回归的 $\sigma$ 在两个 TS 上都诚实**。ts0000：$\hat\kappa - \lambda = 1.9\times10^{-4}$
   对 $\sigma = 1.4\times10^{-4}$。ts0002：两者到 6 位小数相同，因为 $c = 1.77$ 小，
   带内几乎不需要外推。
3. **$c$ 把两个 TS 区分开了**：OH 扭转 $c = 14.28$，甲基转子 $c = 1.77$，差 8 倍。
   这是简谐数值完全丢弃的信息，也是缺口 3（排序启发式）可用的新信号。

**第二批后重跑（`p0_run3`）逐字段与第一批后（`p0_run2`）相同，也就是与基准相同。**
逐字段对照：`status` / `trials` / `chemical_key` / `fragments` / 三个 microstate 能量 /
`reservoir` / `reactions` / `ts_connected_channels` / `conformer_transitions` /
三个 TS 的能量与连接关系 —— **全部 identical**；唯一不同的仍然只有虚频，且偏移量与第一批
完全一致。`free_bonds` 三个 microstate 均为 `()`，与预测一致（乙醇所有扭转都在 soft 带），
所以自由坐标商在这个体系上是恒等操作 —— **它没有被这个基准检验过，只被解析测试检验过**。

运行环境注意：本次 24 次试验的墙钟约 10 分钟，但机器 load average 一度到 60 且
`/home/ruigengji` 是 NFS（§9），所以不可与 §5 的 512 s 直接比较。

### 8.10 第一批没有做、第二批做掉的

- **球面最小化**（找 persistent bottom direction）。原语齐了（`response_gradient`，2 次力评估/步），
  但 persistence 的容差在乙醇上定不了，需要真实化学基准（缺口 8）。
- **$\nu_{\rm band}$ 的嵌套计数**。只对单个方向做回归，还没有在正交补里逐层往下。
- **Lanczos**。优先级已下调（2.6b 第 2 条、§A.13.5）。
- **外部参照**。平滑偏差仍未测，见 2.6c 末段。这是整条线最大的空缺。

## 8.11 第三层：reaction-event identity（已落地）

丙二醛暴露的**不是 `chemical_key` 的 bug，而是 reaction-channel 层缺失**。
两个问题从来就不等价：

$$\text{“产物是不是新的 chemical state?”} \qquad\text{与}\qquad \text{“轨迹有没有发生 chemical reaction?”}$$

`chemical_key` 商掉状态的一切冗余 —— 包括让退化重排回到出发点的那个自同构 —— 这是**对的**，
否则转移三个等价 H 中的哪一个会产生三个产物（§4 早就写了这条）。但同一个商也让它对**路径**
失明：质子从一个氧移到等价的另一个氧，回来的 key 完全相同，而它握着的键**变了**。

所以 key 决定节点，**键差决定事件**。§3.6 一行不用退。

### 8.11.1 判据

对固定原子身份的图 $G_0=(V,E_0)$、$G_1=(V,E_1)$，令 $B^-=E_0\setminus E_1$、$B^+=E_1\setminus E_0$：

| `chemical_key` | $B^-\cup B^+$ | 分类 |
| --- | --- | --- |
| 相同 | 空 | `conformational_transition` |
| 相同 | 非空 | **`degenerate_reaction`**（automerization） |
| 不同 | 非空 | `reaction` |
| 不同 | 空 | `key_change_without_bond_change` —— key 是键图的函数时不可能到达，**命名出来是为了诊断** |

`chemistry.classify_transition`。乙醇甲基旋转 $B^-=B^+=\varnothing$，仍然是
`conformational_transition`，**已验证过的 §3.6 结论一条不变**。

### 8.11.2 事件身份也要再商一次，商的是 fingerprint 不是状态

直接存原子索引会让等价原子重造重复通道（转移三个等价 H 中的哪一个 = 三条通道）。
`chemistry.reaction_event_key` 做两重商：

$$k_{\rm rxn}=\operatorname{Canon}_{\operatorname{Aut}(G_0)\cup\operatorname{Aut}(G_1)}(B^-,B^+),
\qquad k_{\rm channel}=\min\big[k_{\rm rxn}(B^-,B^+),\,k_{\rm rxn}(B^+,B^-)\big]$$

方向也商掉：正反两次观测是**同一条通道的两个 event**，与 TS 通道"无向结构性证据"的既有
逻辑一致。自同构枚举超限时只用恒等，于是**多算通道而不是错误合并** —— 与其余对称处理同一个
失败方向。

**一个实测到的碰撞及其修法**：丙二醛与 3-氧代丁醛的事件在同样的索引上断 O–H、成 O–H，
所以只编码键模式时**拿到了同一个 key**（`1d2c84d4794adb85`）。通道身份必须同时指明是哪个物种，
所以摘要里加了**端点图哈希的无序对**（无序 → 不给 key 方向）。已写成回归测试。

### 8.11.3 网络层：允许 self-loop 的 multigraph

$$G_{\rm chem}=(V,E_{\rm rxn}),\qquad e: C_i\to C_i \text{ 合法}$$

新增 `network["reaction_channels"]`，**退化与普通反应统一放在一起**，不另造一个
`degenerate_reactions` 世界：

```text
reaction_channels:
  rc0000:
    event_key: ...        class: degenerate_reaction | reaction
    ends: [c0000, c0000]  self_loop: true
    broken: [[4,5]]       formed: [[0,5]]
    undirected: true      observed_by: [t000012]   ts_support: [ts0003]
```

**self-loop 不开新节点、不消耗 chemical-node 预算**，所以两层网络 + 双预算的全部逻辑保持不变。

行为上的实质修正只有一处：原来"端点同属一个 chemical node → `conformer_transitions`"
现在改成先分类 —— 只有 $B^-\cup B^+$ 为空才进 `conformer_transitions`，非空进
`reaction_channels`（`self_loop = true`）。TS 记录与 attempt 记录都新增 `transition_class`
与 `reaction_channel`。

### 8.11.4 preflight（零搜索预算，已实测）

```
P1 丙二醛（对称）        A/B 松弛后 E 都是 -7274.549375 eV（逐位相同）
                        graph_hash 同  chemical_key 同  broken [[4,5]] formed [[0,5]]
                        >>> degenerate_reaction
P2 3-氧代丁醛烯醇（不对称） graph_hash 异  chemical_key 异  broken [[4,5]] formed [[0,5]]
                        >>> reaction              dE = -53.96 meV
```

两个体系都是 C/H/O、中性、闭壳层，在 MACE-OFF 域内。P1 两个互变异构体的能量**逐位相同**，
正是退化互变异构应有的样子。

（早先用手搭的未松弛产物几何做同一判断时，`locked_bond_parity` 曾给出"key 不同"的假象 ——
松弛后消失。**这类判断必须在松弛结构上做。**）

## 8.12 locked-bond 判据改成测量（P1 抓到的缺陷，已修）

### 8.12.1 缺陷

`chemical_key` 的 `locked_bond_parity` 用**一个硬几何阈值**决定哪些键能承载 E/Z：

$$\text{locked} \iff d_{ab} < 0.93\,(r_a + r_b)$$

C–O 的阈值是 $0.93 \times (0.76+0.66) = 1.3206$ Å，而丙二醛烯醇的 C3–O2 就坐在上面。
P1 运行里两个**每根键相差不到 0.0004 Å、能量差 0.07 meV** 的结构，因为该键分别在阈值
**+0.00021 Å** 与 **−0.00011 Å** 两侧，拿到了不同的 key，于是开出了一个假 chemical node。

**而且翻来翻去的那根键承载的是构象自由度**：C3–O2 上挂的是羟基转子，同一次运行里的
m0001/m0003 就是转出去的构象（高 0.36 eV）。把它的取向记成构型宇称，**直接违反 §3.6 自己
写下的规则**——"立体化学签名只取构型手性，绝不含可旋转单键的扭转符号"。

调阈值修不了：同族结构里 C3–O2 从 1.2280 连续变到 1.3540，任何静态阈值都会被扫过。

### 8.12.2 测量把设计定死了

| 键 | 长度 | 0.93 阈值 | $k_\phi$ (eV/rad²)，c0000 → c0001 |
| --- | --- | --- | --- |
| C1–C2 | 1.4430 | 1.4136 | 2.129 → 2.145 |
| C2–C3 | 1.3545 | 1.4136 | 4.566 → 4.591 |
| C3–O2 | **1.3208** | **1.3206** | 2.441 → 2.469 |

长度判据在 C3–O2 上的余量是 **0.015%**；刚度的结构间散布约 **1%**，对阈值 4.0 的距离是
**40%**。阈值 4.0 的依据：乙烯 C=C 的二重扭转势垒 2.7 eV 给 $k = 2V = 5.4$ eV/rad²，
所以 4.0 收得下真双键、排得掉共轭单键。

**注意 `soft_mode_curvature_max_eV_rad2 = 50` 不能拿来当这个阈值** —— 它回答的是"fmax 是否
已经管住了这个坐标"，是另一个问题；任何真实扭转都到不了 50。

### 8.12.3 实现：与 §8.8.2 完全对称

自由坐标的判据是"测量出来的、随几何走、绝不从图上猜"。**锁定键是同一个问题的另一端**，
现在用同一个模式：

- `network.locked_bonds_of(structure, threshold)` 从 `structure.info["quench"]` 最后一轮的
  `curvature_eV_rad2` 读出，取超过阈值的中心键。
- **没有淬火报告时返回 `None` 而不是空集** ——"没测过"和"没有锁定键"是两个不同的主张，
  返回空集会把前者悄悄变成后者。
- `chemical_key(..., locked_bonds=...)`：给了就用，没给就退回长度判据，
  并在 key 里记 `locked_bond_source: "measured" | "geometric"`。
  **两种来源的 key 因此不可比** —— 混用会立刻表现为不匹配，而不是悄悄合并两个从未在同一
  标准下比较过的状态。
- `Registry.key_of` 传入测量结果。

### 8.12.4 实测：修复生效，但刀口只是搬了家

两个 P1 结构携带的淬火报告给出
`{(1,2): 2.129/2.145, (2,3): 4.566/4.591, (3,4): 2.441/2.469}`，
`locked@4.0` 都是 `{(2,3)}`，**key 相同**。（`info["quench"]` 挺过了 extxyz 往返。）
基线 `122 passed, 1 skipped`，新增 3 条测试直接钉在暴露缺陷的那两个结构上。

**但整轮重跑仍然开出两个节点。**把全部四个 microstate 的刚度量出来：

| microstate | $E$ (eV) | C1–C2 | C2–C3 | C3–O2 | locked@4.0 | O1···H |
| --- | --- | --- | --- | --- | --- | --- |
| c0000/m0000 | $-7274.54937$ | 2.129 | **4.566** | 2.441 | {(2,3)} | 1.744（螯合） |
| c0000/m0001 | $-7274.54929$ | **4.548** | 2.119 | — | {(1,2)} | 0.987（螯合，另一互变体） |
| c0001/m0000 | $-7274.19020$ | 0.729 | 2.107 | 0.543 | {} | 3.788（羟基转出） |
| c0001/m0001 | $-7274.23618$ | 2.733 | 0.658 | — | {} | 0.962（羟基转出） |

两件事：

1. **在螯合家族内部修复是对的**：哪根键刚随互变体互换（m0000 是 C2–C3、m0001 是 C1–C2），
   而自同构不变的标签让两者都给出**一条**条目 → 同一个 key。
2. **但断开分子内氢键会摧毁共轭**：同一根 C=C 的刚度从 **4.57 掉到 2.11，甚至 0.66** ——
   跨越 7 倍。于是非螯合构象一条锁定键都没有，key 又不同了。

$$\boxed{\text{按结构测出的“锁不锁”不是物种的属性，因此不能定义化学身份}}$$

同一个分子的同一根键，刚度在构象之间横跨 $0.66 \to 4.59$ eV/rad²。**任何阈值都会被扫过。**

**这也更正了 8.12.2 里的一个说法**：那里写"刚度对阈值 4.0 有 40% 余量"，那是只量了两个
**同族**结构得出的；把全部构象算进来，同一根键的跨度是 7 倍，余量根本不存在。

### 8.12.5 结论：几何与刚度都不能定义身份

| 判据 | 在 P1 上的表现 |
| --- | --- |
| 键长 $d < 0.93(r_a+r_b)$ | C3–O2 相差 $3\times10^{-4}$ Å 跨过阈值 → 两个假节点 |
| 扭转刚度 $k_\phi > 4.0$ | 同一根 C=C 在螯合/非螯合构象间 $4.59 \to 0.66$，跨 7 倍 → 仍然两个节点 |

$$\boxed{\text{“锁不锁”不是构象的物理测量属性，而是电子结构/键级语义}}$$

**成对判定（把 E/Z 交给节点内逐对比较）也被否掉了，理由不是代码不漂亮而是传递性**：
设 $A = \mathrm E$、$B$ 是高度扭曲的中间构象、$C = \mathrm Z$。若 $B$ 上那根键一度被测成可旋转，
就会出现 $A \sim B$、$B \sim C$、$A \not\sim C$。节点身份必须是等价关系，
而**搜索算法恰恰会主动生成这种 $B$**，所以这是现实危险不是数学洁癖。

$$\boxed{\text{chemical identity 必须由构象无关的离散电子/拓扑语义定义}}$$

几何、刚度、Hessian 可以**验证**身份，不能**定义**身份。

## 8.13 立体化学本体论：admissible bond-order assignments（已落地）

### 8.13.1 不推断唯一 Lewis 结构，推断一个集合

$$\mathcal B(G) = \{\text{所有满足中性闭壳价态约束的 bond-order assignment}\}$$

边 $e$ 只有在

$$\forall B \in \mathcal B(G):\quad b_e \ge 2$$

时才允许贡献 E/Z parity。

读单个 assignment 会让结果变成"挑了哪个 Kekulé 结构"的产物；问"所有 assignment 上都成立
什么"才是良定义的，而且**构象无关** —— 这正是身份需要的性质。

### 8.13.2 各种坑自然消失（实测）

```
c0000_m0000（螯合）   solutions=1  locked=[O1-C1, C2-C3]   assignment: O1=C1 2, C1-C2 1, C2=C3 2, C3-O2 1
c0001_m0000（非螯合） solutions=1  locked=[O1-C1, C2-C3]   assignment 与上完全相同
benzene              solutions=2  locked=[]               两个 Kekulé，无一条边恒为双键
ethylene             solutions=1  locked=[C-C]
ethane / ethanol     solutions=1  locked=[]
```

**两个 P1 结构给出逐项相同的 assignment**，所以它们的 key 相同 —— 缺陷消除，
且消除的方式与构象无关。

### 8.13.3 无法解析时显式声明，绝不偷猜

`VALENCE = {H:1, C:4, O:2, F:1, Cl:1, Br:1, I:1}` —— C/H/O 覆盖 P0/P1/P2；单价卤素是终端
原子、在中性有机物里无歧义，一并纳入。**N、P、S 故意留在外面**（它们在这个设定下的价态不是
单值的）。

超出价表、价态无解、或枚举超限时返回 `None`，`chemical_key` 记
`stereo_unresolved: "<reason>"` 并**不贡献任何 parity**。这个标记**在 key 里面**，
所以未解析的物种绝不会和已解析的物种悄悄比较。代价要说清：**在未解析类内部，
两个 E/Z 异构体会被合并** —— 这是"显式声明"换来的，不是免费的。

### 8.13.4 删掉的东西

`locked_bond_ratio` 几何判据、`network.locked_bonds_of` / `_stiffness_by_bond` 刚度判据、
`config.locked_bond_curvature_eV_rad2` 全部移除。**不再有任何按结构测的刀口判据参与身份。**

基线 `123 passed, 1 skipped`。回归测试直接钉在暴露缺陷的两个结构上（长度仍然把它们分开、
键级把它们合并），另加苯的双 Kekulé、乙烯、乙烷，以及 NH₃ 的 `unsupported_elements`。

## 8.14 保留 transient crossing（P1 缺口 9 的第一优先级修复）

`transient_topology` 原来把**端点一起丢掉**：P1 里 0.35 Å 那次质子确实转移了，
却因为新图没撑过自由 MD 尾部三帧而返回不带 endpoint 的 `Outcome`，
`confirm_minimum` 根本没跑到它 —— **最接近分界面的那个幅度连 TS 检验的门都没进**。

现在：

- 端点照常交给 `confirm_minimum`（鞍点就是鞍点，与产物是否稳定无关），
- 但**产物主张不被信任**：`admission = "withheld_untrusted_topology"`，不进网络，
- 记录 `crossing: {observed: true, tail_persistence_steps: ...}`。

幅度细化随之改为把

$$\text{elastic} \to \text{transient crossing} \to \text{confirmed product}$$

当作**有序信息**而不是"失败、失败、成功"：transient crossing 算作"离开原盆地"，
可以**打开并驱动 bracket**，`high_kind` 记 `confirmed | transient_crossing`。
三者之中 transient 那一档离分界面最近，信息量最大。

---

## 8.15 TS 发现接线：min-mode following 接入搜索（已落地）

### 8.15.1 五步接线

1. `register_saddle` 从 `execute()` 里抽出，供两条路径共用；`found_by` 记
   `quench_stalled | min_mode_following`。
2. 新增 `seek_saddle`：一次 amplitude scan 观测到 reaction event 之后，从该扫描的帧取种子
   调 `follow_min_mode`。
3. 成功则走同一套 TS 登记 + `connect_saddle` 两侧下降。
4. 搜索级集成测试（`SharpWell`，见 8.15.3）。
5. P1 同输入同 seed 重跑（§8.15.5）。

预算独立于 `max_trials`：`saddle_search_enabled` / `saddle_search_per_run = 8` /
`saddle_search_seeds = 2`。它花的是力评估不是试验，种子来自已经做过的工作。

### 8.15.2 我把 min-mode following 换成纯 Newton，测试把它打回来了

第一版用"解析 Hessian 本征基里的纯 Newton 步"，理由是定步长的
$F_{TS}=F-2(F\cdot v)v$ 发散。集成测试直接否掉了它：

```
SharpWell  seed r=1.9146  E=0.6122   climb -> wrong_index, order=0
```

诊断链：探针交付**不保存中间帧**，轨迹从 r=1.2 直接跳到 1.9；而该势的拐点在 r=1.885，
所以所有可用帧都落在**正曲率**区。**Newton 收敛到最近驻点，正曲率处就下山** ——
它是个从高处出发的最小化器，不是爬鞍算法。

改回真正的 min-mode following：最软模上 $\lambda < 0$ 用 Newton（本身即上爬）、
$\lambda > 0$ **逆着力走固定一步**；其余模各走各的 Newton。同一个种子 **8 步到 r=1.8000,
order 1**。

$$\boxed{\text{这个区分就是 min-mode following 的全部；没有它就只是个最小化器}}$$

顺带查清：`t000002` 的拓扑切换被记录在**淬火相**，因为**键图迟滞**（0.08）把切换点从
1.824 推到 1.97。所以"取拓扑切换帧"这个种子策略在迟滞 + 单步交付下并不可靠；
现在是切换帧优先、退回 `raw_endpoint`，而 stepper 的鲁棒性让种子选择不再是成败关键。

### 8.15.3 集成测试

`SharpWell`：同样两个极小（1.2 / 2.4），脊上加窄高斯，$\lambda_1$ 从 $-0.74$ 变
$-14.52$、势垒 1.0 eV，且剖面无多余驻点（已逐点核对）。三条搜索级测试 ——
闭环（reaction → 帧 → 爬 → 鞍点 → 两侧下降 → channel 引用）、
关掉新路径时 `found_by` 全是 `quench_stalled`、transient 产物"检验但不入网"。

注意这个势**关掉新路径仍会产生 TS**（一维体系只有一个内坐标，停在脊上远比九原子分子容易），
所以对照写的是"鞍点从哪来"而不是"有没有鞍点"。

### 8.15.4 种子决定爬到哪个鞍点（实测）

同样三个丙二醛种子，换成真正的 min-mode following 之后**全部到达一阶鞍点，但不是同一个**：

| 种子 | $E-E_{\min}$ | 步数 | 到达的鞍点 | $\lambda_1$ | icm | 势垒 |
| --- | --- | --- | --- | --- | --- | --- |
| 0.250 弹性 | 69.3 meV | 97 | O1–H 3.04 / O2–H 0.96 / O···O 3.64 | $-0.134$ | $-191$ | 697 meV |
| 0.275 越过 | 493.5 meV | 97 | 同上（镜像） | $-0.135$ | $-191$ | 697 meV |
| **0.350 transient** | 796.3 meV | **8** | **O1–H = O2–H = 1.191 / O···O 2.331** | $\mathbf{-37.33}$ | $\mathbf{-3186}$ | **419 meV** |

前两个爬到的是**羟基转子势垒**（六元环打开），不是质子转移。

$$\boxed{\text{爬山总能找到“一个”鞍点；只有两侧下降才能说它属于哪条通道}}$$

### 8.15.5 P1 第四轮：9/13，失败项全部是真实信息

`runs/p1_round4/`，同输入同 seed。

```
FAIL  A1  chemical_nodes == 1     got 2
PASS  A2  reactions == 0          got 0
FAIL  A3  no microstate under a different key   ['420aa409', 'af73cca8']
PASS  B1/B2/B3   rc0000 degenerate_reaction self_loop=True  event_key=a04410b76e16b117
PASS  C1/C2      []  []
FAIL  C3  every fragment closed-shell           ['C3H3O2', 'H']
PASS  D1  every TS saddle_order == 1            1 TS      <- 不再空过
PASS  D2  every curvature_source analytic       analytic
PASS  D3  trivial floor residual < 5e-4         4.46e-16
FAIL  D4  supporting TS sign_stable             response 缺失
```

**D 组第一次不是空过** —— `ts0000`、`found_by=min_mode_following`、order 1、解析 Hessian、
地板残差 $4.5\times10^{-16}$。接线本身成立。

**但三个新问题，都由这一轮暴露：**

1. **爬山走出了模型域。**`ts0000` 的 icm 是 $-769.7$、$E$ 高出极小 **2.34 eV**，两侧下降给出
   `c0001 = ['C3H3O2', 'H']` —— **O–H 均裂的自由基对**。C3 判据正确触发（这是冻结时写下的
   预期失败模式第 2 条）。A1/A3 的失败是它的后果，**不是 identity 层退化**：第三轮同一份
   identity 代码给出 `nodes == 1`，这一轮多出来的节点是真正的新化学，只是那个化学在
   MACE-OFF 域外（§2.5）。
   **信号已经存在但没有被用作门**：bond-order solver 对该节点给出
   `stereo_unresolved: "no_assignment_satisfies_valences"` —— 开壳层没有中性闭壳价键解。
2. **爬来的鞍点没有带内回归。**`follow_min_mode` 的诊断里没有 `response` 字段，
   所以 `ts_entry["response"]` 是 `None`，D4 无法评估。
3. **我引入并已修掉的一个错误归属。**`seek_saddle` 曾把找到的 TS 直接挂到"发起这次爬山的
   那条通道"的 `ts_support` 上。**通道归属只能由两侧下降决定**（`register_saddle` 里的
   `record_channel` 本来就在做），按发起原因归属是编造。现在只记
   `requested_channel` / `supports` / `supports_requested` 三个字段，事实与意图分开。

### 8.15.6 三个 bug 已修（域门、爬来鞍点的带内回归、爬山能量界）

1. **域门**（`network.admit`）：bond-order solver 给出
   `no_assignment_satisfies_valences` 或 `degree_exceeds_valence` 时拒绝入网，
   `reason = "out_of_domain"`。信号本来就有，只是没被用作门。
   `unsupported_elements`（N/P/S）是**另一种判决，不拒绝** —— 那只意味着立体化学未解析。
   开关 `closed_shell_only`，因为电荷/自旋感知的势可以合法处理自由基。
2. **爬来的鞍点补带内回归**：`follow_min_mode` 确认 index 1 后沿不稳定模跑
   `response_curvature`，D4 因此可评估。
3. **爬山能量界**：见 8.15.7。

**开域门时 10 个测试全挂**，因为合成双井的源结构 `C2`（两个碳、度数 1、需四重键）
与链的 `C4` 都没有中性闭壳价键解。`DoubleWell` 的 docstring 自己写着
"software benchmark, **NOT carbon chemistry**"，所以正确处理是在这些 fixture 上显式关门
并注明理由，而不是放松门去迁就一个明说不是化学的势。CLI demo 同样处理。
新增真实域门测试：把丙二醛质子拽走 6 Å → `rejected / out_of_domain`、`registry.nodes == []`；
NH₃ 走 `unsupported_elements` 分支不被拒绝。

### 8.15.7 ceiling 修的是机制，不是把 1.5 eV 调高

原实现一步越界后**下一轮才检测并直接终止**，把可回退的步长过冲变成整个 seek 失败。
四条机制（上限保持 **1.5 eV**）：种子起点越界 → `seed_above_ceiling` 跳过换种子；
步越界 → 减半回溯（`min_mode_backtracks = 8`）才报 `ceiling_blocked`；
**收敛后的鞍点仍须低于 ceiling**（否则界限可用更小的步绕过，而界限是关于区域的）；
记录 `initial_rise_eV` / `max_rise_eV` / `ceiling_hit_step` / `backtracks`。

两类回归测试都是实测调出来的（`min_mode_step_A = 0.30`、ceiling $+0.06$ →
`backtracks = 2`、`max_rise = 0.0593` 仍收敛；ceiling $+0.02$ → 无论怎么回溯都拿不到结构）。
写测试时踩的坑值得记：第一版用 φ=48° 作种子时 `backtracks == 0`，因为**沿反应坐标越过顶点
后能量是下降的**，冲过鞍点并不触发能量上界；真正触发它的是过长的步在正交（键/角）方向抬高
能量。所以这个测试必须实测调参，不能推。

### 8.15.8 结果与三个语义边界

P1 第五、六轮 **13/13、D 组不再空过**，质子转移通道由两侧下降独立支持
（`icm = -3189.6` vs ground truth $-3186.0$，势垒 419.14 meV 逐位相同）。
完整指标与归因见 [P1_MALONALDEHYDE_ACCEPTANCE.md](P1_MALONALDEHYDE_ACCEPTANCE.md)。

**`requested_channel`（意图）/ `supports`（两侧下降的事实）/ `supports_requested`
（种子策略命中统计）三者分列**，且任一侧下降失败时后者记 **null 而非 false** ——
false 是一个主张，在证据缺失时断言它等于把"没测"变成"负结果"。

三个分母不同的指标（`docs/experiments/saddle_metrics.py`）：round 5 与 round 6 的
`seek_success_rate` 都是 0.50 而**归因相反**（前者是实现缺陷，后者是起点 1.278 eV 的种子被
界限正确拒绝）—— 单看成功率会把两轮当成一样，这正是要三个分母的理由。

---

---

## 9. 环境

```
解释器  /home/ruigengji/miniforge3/envs/openmm_dev/bin/python
模型    /home/kasuga/.cache/mace/MACE-OFF24_medium.model
参考    /home/ruigengji/ABFE_IBS/ABFE_IBS        （只读；行号 6379/6409/7940 已核对未漂）
GPU     RTX 2080 Ti，OpenMM 平台 Reference/CPU/CUDA/OpenCL
```

仓库根的 `.venv/` 是 Windows 布局的遗留物，Linux 下不可用。

**`/home/ruigengji` 是 NFS 挂载**（`192.168.1.100:/home/ruigengji`，vers=3、hard、timeo=600），
解释器、项目和输出都在上面；`/home/kasuga`（模型缓存）在本地 btrfs。后果：

- **进程启动可能极慢**，且卡在 `D (disk sleep)` / `rpc_wait_bit_killable`。实测在机器
  load average 60 左右时，一个只做单点计算的脚本从启动到第一行输出超过 10 分钟，
  全部耗在 NFS 读 conda 环境的动态库上，与 CUDA 无关。
- **任何"墙钟时间"测量必须区分启动与计算**。§8.6 的 1.62 s / 2.92 s 是模型载入之后的
  GPU 侧计算，不受影响；§5 的"24 次试验 512 s"包含启动，不可跨机器比较。
- 判断卡住与否看 `/proc/<pid>/wchan`：`rpc_wait_bit_killable` 是 NFS 等待，不是死锁。

单位契约：对外 Å、eV、eV/Å、fs、amu、K、rad。OpenMM 边界内部 nm、kJ/mol、kJ/mol/nm、ps，转换只在 `openmm_backend.py` 发生（`EV_PER_KJ_PER_MOL`、`NM_PER_ANGSTROM`）。

## 8.16 一次审计与七处修复（目标反应模 vs 全局最软模）

P2 前两轮的公式都是对的，错的是**"从数值结果能推出什么"的推论链**。审计指出七处，
全部修完才跑第三轮。按严重程度：

1. **`follow_min_mode` 爬的是全局最低模，`direction` 参数存在但未参与选模。**
   对软转子密布的分子，"最软鞍点"与"目标反应鞍点"是两个不同的问题。改成目标模跟踪：
   第一次 Hessian 更新按 $\arg\max_i|v_i^\top \hat d_q|$ 在**内部模**里选，之后按与上一轮所选
   模的重叠连续跟踪、符号翻转保持分支；重叠低于 `min_mode_overlap` 返回
   `lost_target_mode`，**不静默换主题**。不给 `direction` 时保留原语义（最软非平凡模），
   因为那是另一个合法问题，由调用方选。
2. **"发现任意 `ts_candidate` ⇒ 沿目标坐标已 departed"不成立。**转子鞍点只证明到了某个盆地
   边界。改为 `high_kind = "boundary_unattributed"`：仍可开幅度括号（"返回在此之间停止"成立），
   但不冒充目标方向的 departure。
3. **端点下降靠把位移放大到 1,2,4,8 倍也不能证明局部连通性** —— 大位移可能跨过多个分界面。
   改为**鞍点间中继下降**：步长始终小，落到另一个一阶鞍点就从**那个**鞍点沿它自己的不稳定模
   继续（`irc_relay_hops`），每跳记录，并注明"一串跳比一步是更弱的连通性证据"。
4. **接线**：`seek_saddle` 原本要求扫描中已有 `reaction_channel`，纯 `ts_candidate` 进不了
   种子过滤 —— 所以第 2 条的修改本来根本不会生效。现在接纳，但**无目标证据时
   `requested_channel` 保持 null**，不凭"探测了哪个方向"编造 intent。
5. **`sigma` 不是校准过的统计不确定度**，是确定性 OLS 拟合的残差尺度，对"方向选错"这类系统
   误差是盲的（实测它读 $3.7\times10^{-6}$ 而值错了 0.18）。`significant` 改名
   `exceeds_residual_scale`，并附 `sigma_meaning`。
6. **球面搜索线搜索一步都没接受时标了 `converged=True`** —— 改为 `stalled`。
7. **爬来的 TS 没有过极小检验那道刚体残差门**，数值噪声负模可能被当成 index-one saddle。
   已补。顺带发现 `follow_min_mode` 会把 `GateRejected` 抛进搜索循环，改为失败关闭返回
   `gate_*`。

### 8.16.1 切向的质量度量：$M^{-1}g_x$ 而不是 $M^{1/2}g_x$

`internal.gradient` 返回的是协向量 $g_x = \partial q/\partial x$。在质量度量下使该内坐标变化最快
的笛卡尔位移是

$$d_x \propto M^{-1}g_x, \qquad d_q = M^{1/2}d_x \propto M^{-1/2}g_x$$

原实现把 $g_x$ 直接交给 `_projected_direction`（乘 $M^{1/2}$），**质量因子方向反了**。
真实 P2 帧上的实测：修正后 O–H 切向的逐原子范数是 **O = 0.0625、H = 0.99206，H/O = 15.87**
—— 正好是质量比 $15.999/1.008$。改之前会给 O 比 H 大 4 倍的权重，也就是把爬山瞄向几乎不动
的那个原子。

**帧差分得到的是位移（逆变），不再乘 $M^{-1}$；只有梯度那条回退路径需要。**

### 8.16.2 内坐标奇点必须显式判据，`isfinite` 抓不住

共线二面角上 `internal.gradient` **不报错**，返回一个有限但荒谬的向量：范数 $3.3\times10^6$、
**在全部 12 个原子上都非零**（叉积是 $3.5\times10^{-18}$）。加了两道判据，都通向 fail closed：

- **几何奇点**：相邻两腿的 $|b_i\times b_{i+1}|/(|b_i||b_{i+1}|) < 10^{-3}$ 即拒绝 ——
  不是调出来的数，就是公式开始除零的地方；
- **支撑检查**：内坐标梯度只该落在它自己的原子上，泄漏到别处就是同一奇点的残渣。

### 8.16.3 指标分母重做

旧的三个率各自静默排除了会让它难看的情形。最要紧的一处：
`requested_channel_coverage` 只把**已观测到的** requested channel 放进分母，所以
"目标坐标完全没采到"报 `n/a` 而不是 0 —— **它本该抓的那个失败，在它里面是不可见的**。

改为六级链条，每级注明条件总体，外加**预注册 truth-set recall**（对整轮打分）：

$$\text{considered} \to \text{frame materialised} \to \text{climb called} \to
\text{TS found} \to \text{attribution complete} \to \text{target supported}$$

三轮对照（`docs/experiments/saddle_metrics.py`）：

| | round 1 | round 2 | round 3 |
| --- | --- | --- | --- |
| 旧 coverage | n/a | 1.00 | — |
| **truth-set recall** | **0.00** | **0.00** | **1.00** |
| 旧 attribution_precision | n/a | 1.00 | — |
| attribution_complete | n/a | 0/1 = 0.00 | 3/3 = 1.00 |

### 8.16.4 采样率的正确表述

$3/11$ 只是**一次、无优先级 stretch 抽样中"目标方向被纳入"的条件概率**，不是整轮的反应命中
率。多次 proposal 时方向纳入率才近似 $1-\prod_j(1-\pi_j)$，反应发现率还要再乘
$P(\text{reaction}\mid\text{direction executed})$。文档里凡出现"27% 命中率"的说法都按此收紧。

### 8.16.5 修的过程中我自己又犯了三次同类错误

值得单记，因为它们与被审计的那类错误同构 —— **从数值结果推出了它支持不了的结论**：

1. `selected_order = 0` 取的是全谱最低本征值，而极小附近那是刚体模、会被跳过 →
   没有任何模走上爬步，整个走法退化成纯 Newton 下降。原来的 `lowest is None` 取的是第一个
   **非平凡**模。
2. 第一版测试拿**初始几何**的模序号去断言**收敛后**的序号 —— 鞍点上负模会排到 0。
3. 第一版测试势用了绑定绝对坐标的锚定项，破坏平移不变性 → 刚体地板 0.389 →
   **新加的地板门把它拒了**。门在做该做的事，测试势才是错的。

另外两次 `bounds[...]` 赋值因缩进不匹配**静默没生效**（`str.replace` 不报错，字段全是
None）。已改用按行定位；这类静默失败在本次会话里出现三次，是个方法问题不是偶然。

### 8.16.6 结果

P2 第三轮 **13/14**、truth-set recall **1/1**，跨节点质子转移通道由三个 TS 的两侧下降支持，
虚频 $-2994.7/-3004.1/-3001.5\ \mathrm{cm^{-1}}$，`target_overlap` 0.77–0.91。
唯一失败的 A3（`reactions >= 1`）是真实的：通道来自鞍点下降，没有一次试验直接观测到
$A\to B$。详见 [P2_OXOBUTANAL_ACCEPTANCE.md](P2_OXOBUTANAL_ACCEPTANCE.md)。

**尚未回答：这次命中是修好了还是运气好。**三个种子的诊断一致说明机制在工作，但 $n=1$。

## 8.17 效率：瓶颈是工作形状，不是算力

实测（`docs/experiments/batch_probe.py`，12 原子、float64、RTX 2080 Ti）：

```
single point            47.4 ms
batch  1   28.0 ms/config   1.7x
batch  4   25.1             1.9x
batch  8   19.5             2.4x
batch 16   13.3             3.6x
```

批处理有效但**次线性**（16 个只有 3.6 倍），因为开销一部分是每批的而非每次的。
更要紧的是**可批的评估只占几个百分点**：

| 环节 | 每次试验的量级 | 能否批 |
| --- | --- | --- |
| 自由 MD | 160 | 不能，严格串行 |
| 淬火 + polish | ~54（实测 60 次试验共 3222） | 不能，串行 |
| `response_curvature` | $2K = 10$ | 能 |
| 曲率探针 | 2/方向 | 能 |
| 解析 Hessian | 1 次调用 | 已在内部 `vmap` 批过 |
| **不同试验之间** | 60 | 粗扫阶梯独立，细化依赖前序 |
| **不同 seed 之间** | — | 完全独立 |

$$\boxed{\text{杠杆在多 seed 并行，不在单点加速}}$$

**多 seed 并行**零共享状态，而且它回答一个单次运行回答不了的问题（8.16.6）。

> **实测更正（§8.21）：它不是近线性的。**这句话当时是从"零共享状态"推出来的，
> 不是量出来的。32-seed 扫描的实测是并发 20 个只有 **3.2×** 聚合吞吐，
> 和批处理 16 个的 3.6× 是同一个量级 —— **GPU 在这个工作负载上大约 3–4× 就饱和**，
> 再加进程只是把同样的吞吐摊得更薄。零共享状态保证的是**正确性**可以并行，
> 不保证**吞吐**可以线性。
**试验级并行有结构性风险**：并行会改变 registry 的接纳顺序，进而改变 microstate 编号与
reservoir 淘汰谁 —— 同一个 seed 跑两次结果不同比慢更糟。要做必须"并行计算、按确定顺序接纳"，
且**应等 P2 的科学结论稳定之后**，否则"改进来自哪"会再次变糊。

**不建议做的两件**：换 float32 或更小的模型（会动物理，而 §2.6c 那批可靠性诊断都在 float64
上做的）；绕开 ASE 直接跑 MD（收益是那 47 ms 里的一部分开销，代价是绕过 `GuardedCalculator`
的失败关闭门，用可靠性换速度）。

---

## 8.18 P2 的能量口径：−11.7 meV 与 −1.75 meV 差在哪

第三轮的 JSON 报 $\Delta E(B-A) = -1.75$ meV，preflight 报 $-11.72$ meV，差 6.7 倍。
先把差额定位到结构和阶段，再改任何东西（`docs/experiments/p2_energy_audit.py`）。

**A 侧一模一样，全部差额在 B 侧。**

```
run   A -8345.351865374   B -8345.353615824   dE  -1.750 meV
pre   A -8345.351865374   B -8345.363586553   dE -11.721 meV
A 侧 run - pre = +0.000000 meV      B 侧 run - pre = +9.970729 meV
```

**不是构象不同。**同一套 genuine torsion 指纹在两个 B 上最大圆周差 **0.0153 rad**。
**也不是不同的驻点。**给它 4000 步 FIRE、fmax $10^{-4}$，它落到
$-8345.363584742$，距 preflight B 的极小 **+0.0018 meV**。所以 round 3 的 B
就是同一个极小，只是**没走完**。

### 8.18.1 走不完的那个坐标：被判为"free"的甲基转子

残差几乎全在一个二面角上：

| | dihedral(2,3,8,9) |
| --- | --- |
| round 3 的 B | +2.4268 rad |
| 收紧后 | +3.1337 |
| preflight B | +3.1416 |

**0.71 rad（41°）。**而 `chemical/c0001_m0000.extxyz` 自带的淬火报告直接写着原因：

```
[4, 3, 8, 10]  tier=free  curvature=0.004114 eV/rad^2  initial_gradient=-0.014009 eV/rad
               note: "flat coordinate; its value carries no information"
```

**曲率在地板（0.01）之下，梯度却是 $-0.014$ eV/rad。**
`polish_soft_modes` 只用曲率判"free"，于是把一个还带着梯度的坐标标成平的、退出打磨，
并且——因为 `free_bonds_of` 读的就是这个 tier——**同时把它从 microstate 比较里商掉了**。

平坦的坐标没有梯度；曲率小只说明**在三重转子势的拐点附近 Newton 尺度不可用**。
实测两者差 26 倍：拐点 0.0041，极小处 0.107 eV/rad²。

这也解释了为什么"fmax 太松"不足以解释 10 meV。把 preflight B 沿各个投影 Hessian 模推到
fmax $=0.02$ 为止，能藏住的能量是：

```
mode  6  k= 0.016  位移 0.036 A   0.046 meV
mode  7  k= 0.055       0.052     0.825 meV
mode  8  k= 0.247       0.009     0.176 meV
mode  9  k= 0.486       0.012     0.137 meV
mode 35  k=39.082       0.001     0.008 meV
```

单模最多 0.83 meV。那 9.97 meV 不是容差**藏**住的，是那个坐标**根本没被打磨**。

### 8.18.2 判据改成"平且驻"，界限是推出来的不是挑的

`free` 现在要求曲率在地板之下**并且**梯度小到即使按地板曲率算也移不出容差：

$$|g| \le k_{\mathrm{floor}} \cdot \delta_{\mathrm{tol}} = 0.01 \times 0.0175
= 1.75\times10^{-4}\ \mathrm{eV/rad}$$

实测的 $0.014$ 是这个界限的 **80 倍**。不满足的进第四档 **`flat_biased`**：
沿该坐标**自己精确的路径**（绕键刚性旋转）扫一个对称约化周期取最低点，再 Newton 收尾。
扫描不依赖刚被证伪的局部模型。

步长上限也来自坐标本身而不是信任半径：**周期的一半 $\pi/n_{\mathrm{sym}}$**——
再长只是绕过一个对称副本。没有这个上限，地板夹住的 Newton 步和原始的一样离谱：
真实体系 $-0.014/0.0041$ 要 3.4 rad，测试转子 $-0.15/0.01$ 要 15 rad。

### 8.18.3 顺带暴露的第二个缺陷：打磨动了几何，淬火却说收敛了

装上扫描后第一次跑，`quench` 返回 **`converged=True` 而 fmax $=0.179$**，容差是 0.02。
原因在 `quench` 的早退：打磨报告自己满意就直接 return，**不再重新最小化**——
可它刚把结构挪了 0.35 rad。

这个洞**不是新的**：旧打磨的步长本来就在容差量级，所以洞一直在，只是很小。
现在的判据是打磨**既满意又没动**才算收敛，否则重新最小化再打磨。

实测 P2 的产物极小需要**五轮**这种交替才停（默认 `soft_polish_rounds` 因此从 3 提到 6；
不动的结构第一轮就退出，A 就是 `steps=0`）：

| | E | 距紧参考 |
| --- | --- | --- |
| round 3 原样 | −8345.353616 | **+9.971 meV** |
| 修好后（5 轮，43 步） | −8345.362884 | **+0.703 meV** |

$$\Delta E(B-A):\ -1.750\ \text{meV}\ \longrightarrow\ \mathbf{-10.894}\ \text{meV}
\quad(\text{紧参考} -11.721)$$

源结构 A **一动没动**（$+0.0000$ meV），所以第四轮和前三轮比的是同一个反应物。

**给基准的直接后果**：`quench_fmax = 0.02` eV/Å 下的节点能量分辨率是**亚 meV 到 1 meV**，
而 P2 的反应能只有 $-11.7$ meV。同一条容差也在 361 meV 的势垒上留着同量级的余量。

## 8.19 交付路径本身要被记下来

第三轮的 A3 失败（`reactions = 0`）和"三个 TS 支持质子转移通道"**同时成立且不矛盾**：
通道是靠鞍点两侧下降命名的，那是**无向的结构证据**；没有一次试验的交付运动到达 B。
而当时无法回答"为什么没落进 B"，因为**只有终态分类，没有路径**。

### 8.19.1 命名坐标与三个此前完全无记录的阶段

试验本来就写轨迹和逐帧观测，但只存了 `pair_distances_A` —— 上三角向量，没有名字。
现在 `SearchConfig.tracked_coordinates` 声明命名内坐标（`bond` / `angle` / `dihedral` /
`bond_difference`），逐帧随 $E$ 一起记。P1/P2 用的是

$$q_{\mathrm{PT}} = r(\mathrm{O_1{-}H}) - r(\mathrm{O_2{-}H}),\qquad
r(\mathrm{O_1{-}H}),\quad r(\mathrm{O_2{-}H})$$

哪个坐标是反应坐标是化学陈述，所以**声明**而不是猜。

更要紧的是**爬山、两侧下降、边界续行这三个阶段此前一帧都没记**。`PathRecorder`
给它们同样的字段、同样的命名坐标，每阶段一个文件（`paths/<id>/<stage>.jsonl` 和 `.extxyz`），
并且**自己绝不触发计算**——只用调用方已有的或计算器缓存里的能量与力，
否则观测会改变被观测的走法的代价。

### 8.19.2 boundary continuation：继续交付，而不是再找一个 TS

一阶鞍点加两侧下降证明的是"两个极小通过它相连"；它**不**证明"被交付的那个运动到得了产物"。
这是两个claim。`continue_across_boundary` 只回答第二个：

从拓扑切换帧出发，沿**probe 自己的有符号局部切向**（切换帧与前一帧之差；没留相邻帧时
退回坐标梯度那条路），每档前进量以**单原子最大位移**设上限（默认 0.05 Å），
**每档从位移后的几何独立淬火**——不是接着上一档的终点，所以一档的判决只关于那一档的位置。

四种判决就是这件事的全部内容：

| 判决 | 含义 |
| --- | --- |
| `returned_to_source` | 该处还没进入产物吸引域 |
| `reached_other_state` | 交付路径到达——**有向观测** |
| `stopped_on_saddle` | 交付路径贴着分界面走 |
| `left_domain` | 步长或幅度把结构推出模型域 |

它**明确不是**两侧下降换个说法：下降从鞍点出发沿本征向量下坡，续行从 probe 帧出发沿 probe
切向前进。把前者当成"交付路径的观测"正是这个函数存在要避免的伪造。

### 8.19.3 有向边分成两类证据，不合并

续行产生的边**不进 `attempts`**。`attempts` 的含义是"这次试验自己从源走到了目标"，
而续行的种子试验恰恰**没有**——它停在边界，这才是要续行的原因。
把它列进 `attempts` 会让种子试验声称一个它从未做出的观测，而且任何按试验记录复核边的消费者
都会发现两者矛盾（写完这条时测试立刻抓到了，见下）。

所以边上是：

- `attempts` —— 自由响应观测（probe → 无偏 NVE → 淬火）
- `continuations` —— `{seed_attempt, mechanism}`，续行观测
- `mechanisms` —— `free_response` / `boundary_continuation`

不变量随之改成"两类之中至少有一类非空，且各自与自己的试验记录相符"：
`attempts` 里的试验记录必须指向同一个 target，`continuations` 的种子试验必须**不**指向它。

### 8.19.4 测试

新增 10 个（合成双井 + 合成转子），总数 **149 passed**：

- 切向方向带符号：同一帧给相反切向必须给出不同判决（正向到达另一态，反向八档全部返回源）；
- 每档位移严格等于 `0.05 × 档数`，且记录里的 $r$ 逐档对得上（上限是单原子位移，不是范数——
  否则 12 原子分子每个原子都能走满一步）；
- 压到原子相撞时报 `left_domain` 并带 gate 码，不当成化学结果；
- 纯平移切向（无内部分量）直接拒绝，不退化成"每档都回到源"；
- 停在鞍点的档**不交给调用方命名**，所以无法从它接纳出一个状态；
- 拐点上曲率在地板下而梯度活着 → `flat_biased`，不是 `free`；
- 真正无势垒的转子仍然是 `free`（修的不能把所有软转子都变成工作量）；
- 扫描落到转子极小 0.1 rad 内；
- 步长不超过 $\pi/n_{\mathrm{sym}}$；
- `flat_biased` 的键**不**被 `free_bonds_of` 商掉。

---

### 8.19.5 第一次跑 P1 就抓到我自己犯的同一类错误

装好路径记录后先在 P1 上跑回归，`docs/experiments/a3_diagnosis.py` 读出来的第一屏就是：

```
continuation001/continuation  q_PT: +0.063 -> -0.770   [-0.792 .. +0.754]
    rung 1 at 0.050 A: returned_to_source order=0
    rung 2 at 0.100 A: returned_to_source order=0
    ... 八档全部 returned_to_source
```

**质子明明换了氧（$q_{\mathrm{PT}}$ 从 $+0.06$ 到 $-0.77$），八档却全被判"回到源"。**
因为我在续行判决里只比了 `chemical_key`，而 P1 的质子转移是**退化的** ——
两侧是对称副本，key 不变。这正是第三层存在的理由，而我在新代码里又把它折叠回了第一层。

改成：**键差决定有没有发生事情，key 决定有没有开新节点。**

| key | 键差 | 判决 | 后果 |
| --- | --- | --- | --- |
| 相同 | 无 | `returned_to_source` | 还没进别的盆地 |
| 相同 | 有 | `reached_degenerate_product` | 自环通道，记 `observed_by`，**不接纳、不花构象槽** |
| 不同 | 有/无 | `reached_other_state` | 接纳 + 通道 + 有向边 |

分类逻辑提到模块级 `rung_verdict`，这样它能对着**已知退化的体系**测，而不是只能通过整轮跑间接测。

### 8.19.6 两个测试势的力符号错了

为退化转移写的 `SharedProton` 一开始每档都淬火到**自己的势垒顶**（$E=0.5$，saddle order 3）。
原因不在续行，在测试势：`forces[]` 一半按 $+\partial V/\partial x$ 累加、一半按
$-\partial V/\partial x$，最后统一取负。**`Rotor` 也一样错**，只是它的键从一开始就在平衡位置，
键力为零，所以完全看不出来。

现在两个势都对数值梯度验过（最大偏差 $2.8\times10^{-10}$ 和 $4.6\times10^{-11}$），
并且这条验证本身成了参数化测试 —— 力和能量不一致的测试势是在拿优化器测优化器。

总数 **153 passed**。

## 8.20 A3 的答案：交付路径到了产物，卡在产物盆地自己的软转子上

P2 第四轮第一次有路径数据可读，A3 的失败机制随即清楚。约定
$q_{\rm PT}=r({\rm O_0{-}H_5})-r({\rm O_4{-}H_5})$，A 是 $+0.703$、B 是 $-0.706$。

**九次试验的交付运动全部走完整条反应坐标**（幅度 0.25 到 1.0，全部
$q_{\rm PT}: +0.703 \to -0.68\ldots-0.70$），九次全部 `status=ts_candidate`、
`target=None`。原因不是没到，是**淬火停在 B 盆地内部的转子鞍点上**
（$-65$ 到 $-79$ cm⁻¹，比 B 的极小高约 11.6 meV）：不是极小 →
`confirm_minimum` 判否 → 登记 TS 候选 → 没有产物、没有有向边。

boundary continuation 同样如此：四次续行都把质子带过去了，但每一档的淬火终点都是一阶鞍点。
**I4 的机制跑通了，被另一个东西挡住。**

### 8.20.1 一阶鞍点的邻域属于一个态，这可以判

证据当时就在第四轮的记录里：`ts0000`（$-75.0$ cm⁻¹）**两侧下降都落在 `c0001/m0000`**，
而且 **B 这个节点本身就是它 $+1$ 侧下降发现的**（`new_chemical_node`）。

一阶鞍点位于恰好两个盆地的边界上。所以两侧下降命名同一个化学态时，
除脊本身之外该点的任何邻域都属于那个态 —— 那一档确实到了。
两侧**不**一致则相反：该档真的在反应脊上，判决保持 `stopped_on_saddle`。

`arrival_through_saddle` 就是这条规则，`arrival = "via_saddle"`，
并在记录里明写**这比直接淬入极小是更弱的证据**。
中继逻辑从 `connect_saddle` 提取成 `descend_saddle`，两处共用一份 ——
第四轮十个停住的鞍点里有九个正是因为**没有**中继而一侧下降失败。

### 8.20.2 被投机调用的回调不能有副作用（第五轮作废）

第五轮报 15/15，A3 通过，三条 `continuations` 里**两条是伪造的**：
`outcome=exhausted`、八档全部 `stopped_on_saddle`，却被接纳了。

原因在结构上：`name_endpoint` 在判决非"回到源"时把端点存进外层字典，搜索侧接纳字典内容。
而这个回调**也会被用在停住那一档的两侧下降上**，那里一侧命名产物、另一侧失败完全正常,
该档的结论恰恰是"没到达"。一侧的端点就这样漏了出去。

这是本次会话第三次同一类错误 —— **从结果推出它支持不了的结论** ——
只是这次经由代码结构而非论证。两条修复：

1. 回调改成**纯函数**，注释写明它被投机调用所以不能有记录能力；
2. **是否到达由 `continue_across_boundary` 自己判**（`_arrived`），
   搜索侧只接纳走法自己交回来的结构；`outcome` 不是到达就交回 `None`。

回归测试：企图记录的回调必须对返回值毫无影响，且 `outcome=exhausted` 必须交回 `None`。

**第五轮作废，不计入判据历史。**A3 即使只算合法的那一条也会通过，
但接线错了就不该拿那一轮的分数。

### 8.20.3 第六轮：15/15，一条完整支持的证据

```
rung 1 @0.050 A  stopped_on_saddle   (+1 命名产物 / -1 下降失败) -> 正确地不算
rung 2 @0.100 A  reached_other_state  arrival=via_saddle
                 该档淬火落在 -74.1 cm^-1 的一阶鞍点
                 +1 侧 3 次中继 -> B (-8345.36313)
                 -1 侧 1 次中继 -> B (-8345.36315)
```

有向边 `r0000`：`mechanisms=['boundary_continuation']`、**`attempts=[]`**、
一条 `continuations`。其余三个种子八档走完仍在鞍点上，`reached=None`，不接纳任何东西。

$\Delta E(B-A) = -11.018$ meV（紧参考 $-11.721$；第三轮读 $-1.750$）。
三个质子转移 TS 与第三轮逐字段相同。**158 passed。**

**仍未回答**：$n=1$；直接观测（`attempts`）仍是 0；四个续行种子只有一个到达，
另外三个是步长不够还是本来贴着分界面，**没有测量所以不下结论**。

---

## 8.22 一个只有热噪声才能暴露的缺陷，以及它如何让运行无法记录自己失败了

P2 温度扫描的 seed 23/47/71 各跑到第 7 个 trial 后死亡，`status` 停在 `running`，
日志三份一致：`ValueError: Out of range float values are not JSON compliant: -inf`。

**根因**：`follow_min_mode` 的 `bounds["max_rise_eV"]` 初始化为 $-\infty$，
而 $-\infty \to$ `None` 的收尾只写在**成功返回**那一条路径上。该函数十一处返回里
**十处是提前返回且都带 `**bounds`** —— `no_internal_mode` 与 `lost_target_mode`
在 step 0 就返回，此时能量块尚未执行，`bounds` 还是初始值。

**为什么 $T=0$ 从来没出现**：那些提前返回在确定性路径上走不到。
热动量把爬山起点打散后才出现 —— **这个缺陷是温度扫描本身暴露的**，
而温度扫描的目的恰恰是"换一个能真正扰动系统的旋钮"。旋钮一旦真的能动，
它就顺手翻出了一处在旧旋钮下永远看不见的缺陷。

### 8.22.1 写记录的失败会吃掉失败本身

$-\infty$ 进入记录后，`atomic_json`（`allow_nan=False`）抛错。而**抛错的那次写正是
`publish()`** —— 失败处理器 `except BaseException` 里做的第一件事也是 `publish()`，
于是它同样失败，**终态状态永远写不出去**。

所以 `status == "running"` 是一个**歧义状态**：正在跑 / 被信号杀 / 写记录时崩了。

我据此推错了一次：从"status 停在 running"推断"连 except 处理器都没跑到 ⇒ 被信号杀"，
再猜是主机内存 OOM，并建议把并发从 8 降到 4。实测节点 125 GB 内存、单进程 RSS 1.9 GB、
`journalctl -k` 无任何 `killed process` —— **两层推断全错，而且那个"修法"对真因毫无作用**。
`prrs: ValueError:` 这个前缀本来就是 `cli.py` 的 `except` 打印的，
它一直在说"这是正常的异常退出"。**是日志推翻了推断，而我在读日志前已经准备好了一个自洽的错误故事。**

### 8.22.2 修法：守卫不放松，但必须说出位置

1. 去掉哨兵值。`runner.py` 中已无 `float("-inf")`。
2. `append_json` / `atomic_json` 失败时遍历记录、**点名每一个非有限字段**：

```
Out of range float values are not JSON compliant: -inf.
Non-finite values at: amplitude_scans[0].saddle_search.seeds[0].max_rise_eV = -inf.
A record may not carry NaN or infinity; fix the source rather than relaxing this check.
```

守卫本身一点没松 —— 记录里出现 NaN/inf 就是缺陷。**改的是它说不说清楚。**
一个不说位置的正确守卫，代价是一次错误的死因诊断。

3. 测试钉性质而不是那一条路径：零曲率势逼出 `no_internal_mode` 早退，
   断言返回的 diagnostics 可写且 `max_rise_eV is None`；另一个断言报错必须点名字段。

### 8.22.3 同一个混淆出现了第二次

`run` 子命令原来用"`network.json` 存在"判"这个 seed 已经跑过"。
但 `publish()` **全程持续重写** `network.json`（status 为 `running`），
所以死在半路的运行也有这个文件 —— 于是它被报成 `already run, skipping`，静默跳过。

这是"完成 vs 失败"混淆的**第二次**出现：上一轮我把它从"目录存在"挪到了"文件存在"，
只是换了个位置，没有解决。现在的判据是 **`status == "completed"`**，
外加 `scripts/scan_status.sh` 用 status + 进程存活 + trial 数三项交叉判断，
把 `DEAD: status running but no process` 单独列出来。

**教训是同一条**：判据要基于**能分辨的信号**。`network.json` 的存在分辨不了完成与死亡，
而它的 `status` 字段可以。我在设计判据时没有去看 `publish()` 什么时候被调用 ——
而那是我早就知道的事。

---

## 8.23 P2 温度扫描：证据阶梯，以及"找到"与"证成"的分离

`temperature_K = 0 → 300`，16 个 seed（素数序列前 16 个），其余配置逐字段等于 round 6。
16/16 完成。协议与逐项失败机制在 [P2_TEMPERATURE_SCAN_PROTOCOL.md](P2_TEMPERATURE_SCAN_PROTOCOL.md)。

$$Q_1 = \tfrac{12}{16}\ [0.51,\,0.90],\qquad
Q_2 = \tfrac{15}{16}\ [0.72,\,0.99],\qquad
Q_4 = \tfrac{14}{16}\ [0.64,\,0.97]$$

方括号是 95% Wilson 区间。**$n=16$ 很小，区间必须和点估计一起报** ——
$Q_1$ 的真值可能低到刚过一半、也可能高到九成，这个宽度不支撑
"0.75 是这个机制的命中率"这类阶段性结论。它支撑的是**与 $T=0$ 的 8/8 的对比**：
那里根本没有波动可测。

对照 $T=0$ 的 8/8/8 —— 那是**同一条确定性路径重复八次**：八个 seed 的 TS 能量到第 7 位
小数相同、续行种子全是 `t000004`、trial 0–16 的 probe 逐个相同。
$T=0$ 时 `thermal_momenta` 返回零动量，而 proposal priority 又把带氢键受体的那条 O–H
确定性地排在最前，所以**被测机制整条路径无随机性可采**。

### 8.23.1 证据阶梯

```
B 作为节点被发现        16/16 = 1.00
质子转移通道存在        15/16 = 0.94
有 min-mode TS 支持它   14/16 = 0.88
有向边落在它上面        12/16 = 0.75
```

**B 作为节点每次都被发现；掉下来的是证据链的层级。**这句话是这次扫描最该留下的一条：
$Q$ 测的不是"能不能找到反应"，而是"**能不能把找到的事情证成**"。
两者在 $T=0$ 下无法区分，因为那里两个都是 100%。

**措辞收紧**：此前写的是"物理每次都被找到"，太强。seed 71 发现了 B 这个节点，
却**没有建立任何 A↔B 通道**（四个鞍点全部只有一侧下降成功），
所以"反应被找到了"在那个 seed 上并不成立。可以断言的只有节点层。

TS 支持与有向边是通道之下**两条独立分支**，不是嵌套 —— seed 73 有有向边而无 TS。

### 8.23.2 四种失败，四种机制

- **seed 71**（唯一 Q2 失败）：四个 quench-stalled 鞍点**全部只有 $+1$ 侧下降成功**，
  $-1$ 侧中继用尽仍非极小 → `connect_saddle` 拒绝命名通道。
  **这是三态归因按设计工作**，单侧证据不足以命名通道。B 确实作为节点存在
  （由那个 $+1$ 侧下降发现），只是没有任何 A↔B 通道。
- **seed 59 / 61**（Q1）：通道存在、有 TS 支持，`observed_by = 0`，续行也没产出。
  6 个边界帧而仅 14–17 个 completed —— 试验大多停在脊上，没有一次干净落进 B。
- **seed 29**（Q1）：通道有 TS 但无直接观测；它**有**一条有向边，是 E/Z 异构那条。
  判据要求边落在**转移通道**上，所以正确地不算 —— 判据在做该做的事。
- **seed 73**（Q4）：直接证据最强（`observed_by = 7`），但唯一的 min-mode TS
  两侧下降没立起通道，且只有 2 个边界帧可作种子。

### 8.23.3 又一个被自己的数据否掉的预测

看到 seed 73（7 次直接观测、0 个 TS）与 59/61（0 次直接观测、有 TS）之后，
我预测 $Q_1$ 与 $Q_4$ **互相交换**：同一个 trial 不可能既干净落进 B 又停在脊上，
而后者才是爬山的种子。

| | $r$ | 95% CI (Fisher $z$, $n=16$) | |
| --- | --- | --- | --- |
| 边界帧 vs 直接观测 | $-0.891$ | $[-0.96,\,-0.71]$ | 不含 0 |
| 直接观测 vs TS 支持 | $+0.101$ | $[-0.42,\,+0.57]$ | **含 0** |
| 边界帧 vs TS 支持 | $-0.062$ | $[-0.54,\,+0.45]$ | **含 0** |

**前提对，结论过强。**竞争确实存在（区间不含 0，而且它本来就是机械必然）；
但后两条只能说**"本样本没看出关系"，不能说"没有关系"** —— 区间宽到 $\pm0.5$，
一个中等强度的真实相关完全落在里面。

机制上的理由仍值得记：`seek_saddle` 只取 2 个种子，且接受
`reacted` / `crossings` / `ts_candidate` 三类，所以一两个边界帧就够 ——
0–1 个边界帧的 seed 17/53 照样各拿 1–2 个 TS。但这是**解释**，不是这组数证实的结论。

### 8.23.4 `attempts` 从 0 到 59，以及它对 boundary continuation 的含义

| | $T=0$ (n=8) | $T=300$ (n=16) |
| --- | --- | --- |
| 直接观测总数 (`attempts`) | **0** | **59** |
| 续行贡献的有向边 | 每 seed 1 条 | 16 个里 3 条 |
| 结果随 seed 变化 | 否（第 7 位小数相同） | 是 |
| $\Delta E(B-A)$ | 恒 $-11.018$ meV | $-10.758 \sim -11.567$（跨度 0.81） |
| 转移 TS 虚频 | 恒 $-2994.7/-3001.5/-3004.1$ | $-3010.9 \sim -2998.8$（n=19） |
| 每 seed trial 数 | 恒 32 | 32 或 60 |

$T=0$ 下淬火总停在 B 内部那个 ~11.6 meV 的转子鞍点上，所以试验自己从不构成产物观测，
A3 只能靠 boundary continuation；$T=300$ 下热动量让淬火逃出去，
**A3 主要由自由响应通过，续行退为少数情形的补充。**

这确认了 §8.20 的判断：**续行解决的问题很大程度上是 $T=0$ 造出来的。**
但它不是没用 —— 283 档里 233 档仍是 `stopped_on_saddle`，淬火真卡住时它仍是唯一手段。

### 8.23.5 $\Delta E$ 的散布就是这套设置的能量分辨率

跨度 **0.81 meV**，紧参考 $-11.721$。$T=0$ 时每个 seed 恰好 $-11.018$（同一条路径），
**所以这个散布是第一次被量出来**。

**它应当称为"本协议下的端点能量散布"，不是能量分辨率。**它至少混了三样：
松弛终止误差、端点落在哪个构象上、以及不同 seed 到达 B 的路径不同。
要把它变成分辨率，必须**把所有端点用同一个紧容差重新松弛**再看剩下多少散布 —— 没做。
数值上它与 §8.18 的独立估计同量级，但"同量级"不是"就是它"。

### 8.23.6 仍然开着的

1. Q1 的 4 次失败里 3 次是"有 TS 无直接观测"，续行在这三个里都没救回来 —— 没查。
2. seed 71 的 4/4 单侧下降，`irc_relay_hops = 3` 在 $-1$ 侧全部用尽 ——
   深度不够，还是 $-1$ 侧本就通向域外，没查。
3. **`ceiling_blocked` 是最常见的爬山终止原因**（每个失败 seed 2–3 次）。
   1.5 eV 是否过紧没测，而且**不能因为它挡住了想要的结果就调高它**。
4. E/Z 异构（key `be2aae715ac7`）在 16 个 seed 里出现 4 次，位于 A 之上 365–404 meV，
   B1 修订后合法通过；它是否还有别的通道未被记录，没查。

---

---

## 9.1 环境变更：openmm 8.5.2 → 8.6.0.dev（2026-09-01，用户执行）

维护后实测：

```
openmm       8.6.0.dev-c6173db      （原 8.5.2）
openmmml     1.6                    （不变）
openmmtorch  1.5                    （不变）
NNPOps       import 成功
torch 2.12.0 / cuda True / ase 3.29.0 / numpy 2.4.3
```

**`tests/test_backend_parity.py` 的 4 项（1 failed + 3 errors）因此变红。搜索侧不受影响：
排除该文件后 `132 passed`。**

> **2026-09-02 复核：环境已回到 openmm 8.5.2，那 4 项现在是绿的
> （全套 `164 passed, 1 skipped, 0 failed`，不排除任何文件）。**
> 但**这不是"修好了"** —— 是环境回滚。8.6.0.dev 下那 4 项为什么红仍然没查，
> 环境再升上去就会重现。上面第 3 条那个"测试探测方式过时"与"上游真修了"的区分，
> 仍然是未完成的判断，**不能因为现在是绿的就当它已经解决**。

三条与既有结论的核对：

1. **§2.4 第 2 条那个窄范围补丁仍然需要** —— `openmmml/models/macepotential.py` 的本地
   `modelPath` 分支依旧没有 `.to(device)`，上游没修。
2. **§2.3 的结论不变** —— 源码里仍是 `PythonForce` 而非 `TorchForce`，
   所以"OpenMM 路线没有吞吐优势"照旧。
3. **`test_upstream_still_needs_the_explicit_device_move` 红了，但被钉住的上游行为并没有
   消失**（见第 1 条）。所以这条红是测试的探测方式碰到了 8.6 的变化，
   **不是"上游修了所以补丁可以删"** —— 后者才是这条测试红了本该表达的意思，
   两者必须分开，否则会据此误删一个仍然必要的补丁。

顺带一条方法论自纠：我最初用 `importlib.util.find_spec("nnpops")` 探测，报 MISSING；
`import NNPOps` 却成功 —— 包名大小写不同。**§3.7"存在性只能由执行证实"这条，
我自己那次探测就违反了。**

### 待办（与 P2 分开，不混在同一次改动里）

- 把 `test_backend_parity.py` 的 4 项在 openmm 8.6 上重新定位：区分"上游行为变了"
  与"测试的探测方式过期了"。**交叉校验后端与 P2 是两件独立的事，混着改会让归因变糊。**
- 在此期间 §2.1 的后端一致性交叉校验**不可用**，所以那一层目前只有单实现。

---

## 10. TODO：接入 geomeTRIC 作外部对照

<https://github.com/leeping/geomeTRIC>。**当前环境未安装**（`geometric` 无 spec，pip 里也没有），
所以下面全是计划，不是实测。

### 10.1 它落在架构的哪一格

$$\text{benchmark 档（validation limit）}$$

geomeTRIC 里的一切都是 $a \to 0$ 的：梯度、Hessian、本征向量跟随。所以它**不能**校验路线 C/E
的有限幅度估计量，它校验的正是 §E.5 那张图上方那一档 —— 那一档也正是唯一有谱定理的地方。
两者不冲突，是分工。

按 §1 的分工原则，它属于**分析/校验层，不得成为搜索核心的运行时依赖**，与 OpenMM 在 §6 里的
定位一致。

### 10.2 值得对照的四件事（按价值排序）

1. **独立的频率分析**，对照 `curvature_spectrum` 的投影 + 质量加权。
   这与 §2.1 的后端一致性检查同一类型：两个独立实现算同一个物理量，逐位比。
   **这是唯一能查出我们投影方案本身有没有错的手段。**
2. **真正的 IRC**，对照 §4.2 的 `follow_unstable_mode`（沿不稳定模 $\pm\epsilon$ 淬火）。
   我们那个是 IRC 的廉价近似；IRC 会告诉我们它命名的端点是否就是 IRC 找到的端点。
   §5 的 6/6 端点定名目前**只有自我一致性**，没有外部确认。
3. **鞍点精化**。TS 记录里的状态字符串至今写着
   `"unrefined first-order saddle candidate"` —— geomeTRIC 的内坐标本征向量跟随
   （P-RFO 一类）能精化它们，并量出"未精化"到底差多少。
4. **内坐标约束优化**。§C.5 那个固定二面角的鞍点构造给出的 $\lambda_1$ 与真驻点差 2%
   （$-283.69$ vs $-291.035\ \mathrm{cm^{-1}}$，因为约束把一个方向排除在收敛判据之外）。
   geomeTRIC 的内坐标约束优化会让这类构造可信。

另外 TRIC 坐标本身是**对 $SE(3)$ 商的一个独立、已发表的实现** —— 也就是路线 B.1 那个
$\mathcal X = \mathcal C/(SE(3)\times\Gamma)$ 的 $SE(3)$ 部分。读它的实现对我们那六处 patch
是有参考价值的。

### 10.3 它明确**不能**做的四件事（不要误读对照结果）

1. **不解决平滑偏差。**它是优化器，不是参照势。用同一个 MACE calculator 跑它，校验的是
   **算法**（坐标、鞍点搜索、IRC、频率投影），**不是 PES**。§2.6c 末段那个空缺原封不动。
   要查曲率准不准，仍然只有 DFT。
2. **没有置换商 $\Gamma$。**它不对自同构取商，所以 §3.6 那半边（chemical key、
   对称感知 RMSD、角度枚举去重）它帮不上。
3. **没有有限幅度响应。**见 10.1。
4. **它的 Hessian 默认是梯度的有限差分。**所以它报的频率会带 §2.6f 那个截断误差
   （闭式 $E''\,\mathrm{step}^2/(r^2m)$，乙醇尺度约 $2.5\times10^{-3}$）。
   **不做处理直接比，会在 $10^{-3}$ 量级"不一致"，而那是它的误差不是我们的。**
   要么把我们的解析 Hessian 喂给它，要么在对照时显式扣掉这一项。

### 10.4 零新依赖的那一步：已落地（`src/prrs/validation.py`）

`ase.vibrations` 已装，自带有限差分 Hessian 与频率，且**不投影刚体模**（我们投影），
所以它是一个真正不同的实现。10.2 第 1 条因此不需要 geomeTRIC 就能做第一轮。

新模块 `src/prrs/validation.py`，**搜索从不导入它**（与 `openmm_backend.py` 同样的定位）。
三个对照，彼此不可替代：

| 对照 | 隔离出什么 | 代价 |
| --- | --- | --- |
| 我们的 Hessian 走 ASE 的谱机制 | 质量加权、模式排序、单位换算 | 0 次力评估 |
| ASE 自己的 Hessian，**同一步长** | 我们把 $3N$ 个 $Hv$ 装配成矩阵的方式 | $6N$ 次 |
| ASE 自己的 Hessian，**不同步长** | 微分本身（截断行为） | $6N$ 次 |

第二条必须用同一步长才有意义：中心差分在同步长下就是同一个计算，那时"一致"证明的是装配
而不是极限。第三条才是对微分的独立检验。

**实测（四原子链，6 个内坐标模，`tests/test_validation.py`）**：

```
phi=0   我们的 Hessian 走 ASE   max|d| = 8.07e-09
        ASE 自己的, delta=0.01  max|d| = 8.07e-09   <- 与上一行相同，装配确认
        ASE 自己的, delta=0.002 max|d| = 1.400e-04   implied C = 1.458
        ASE 自己的, delta=0.05  max|d| = 3.500e-03   implied C = 1.458
phi=60  我们的 Hessian 走 ASE   max|d| = 5.44e-08
        ASE 自己的, delta=0.002 max|d| = 6.749e-05   implied C = 0.703
        ASE 自己的, delta=0.05  max|d| = 1.650e-03   implied C = 0.688
```

$C$ 由 $\max|\Delta\lambda| = C\,|h^2 - h_0^2|$ 反推（$h_0 = $ `minimum_check_step_A` $= 0.01$）。
**同一几何上两个相距 25 倍 $h^2$ 的步长给出同一个 $C$（1.458/1.458 与 0.703/0.688）** ——
所以步长依赖的差异就是截断，不是缺陷。这条已写成回归测试（若 $C$ 不再一致即为缺陷）。

顺带确认：双原子留 **5** 个刚体模而不是 6（绕自身轴的转动不动任何东西），两侧都算对；
虚频与 ASE 的 eV 能量之间的换算往返一致到 $10^{-9}$。

**这一步不覆盖的**：$3N-6$ 个模的**符号结构**两侧一致，但两侧用的是同一个 PES，
所以 10.3 第 1 条不变 —— 平滑偏差依旧不可见。第 2、3、4 条（真 IRC、鞍点精化、
内坐标约束优化）仍需要 geomeTRIC。

### 10.5 前置条件

- geomeTRIC 的 ASE engine 绑定是否存在、在哪个版本引入 —— **文档确认存在，执行未证实**
  （2026-09-04：`geometric 1.1.1` 已装入 `openmm_dev`；文档给出
  `--engine ase --ase-class=... --ase-kwargs=<json>`。按 §3.7 这还不算数，
  必须由执行证实，尚未做。装入后 `pytest 224 passed, 1 skipped` 与 digest 均无变化。）
- **10.3 第 4 条有了绕过办法**：`--hessian file:/path`（3N×3N，numpy 可读文本方阵）
  可以喂入我们的解析 Hessian，从而消掉它默认有限差分 Hessian 的 §2.6f 截断误差。
- **10.3 第 1 条已被独立证实为真且要紧**：2026-09-04 的 ORCA 参照
  （`docs/experiments/{sn2_qm,p1_qm}/`）表明模型误差可达 285 meV 量级，
  而 geomeTRIC 用同一个 calculator 跑**看不见**这一项。参照势只能是 QM。
- 按 §3.7：装上之后存在性与生效必须由执行证实，不能查 metadata。
- 单位契约（§9）在边界上转换，与 `openmm_backend.py` 同样处理。

### 10.6 真实体系上的对照结果（乙醇极小点 + 三个 TS）

脚本 `docs/experiments/validate_ethanol.py`，输出同名 `.out`。解析 Hessian 对照
ASE 的两条路径。

```
                        我们的 Hessian 走 ASE     ASE 自己的 FD Hessian (delta=0.01)
                        max|dlambda|   符号一致   max|dlambda|   负模数 (我们/ASE)
ethanol minimum         1.63e-09       是         5.92e-03       0 / 0
TS ts0000, ts0001       8.39e-05       是         3.33e-02       1 / 1
TS ts0002               2.09e-05       是         8.24e-03       1 / 1
```

**两个独立确认，都不是我预期的那一个：**

#### 确认一：ASE 的 FD 落回基准记的旧值，所以 §8.9 的归因是对的

| TS | 基准记的 FD | ASE 自己的 FD | 差 | 我们的解析 |
| --- | --- | --- | --- | --- |
| syn（ts0000/1） | $-290.08$ | $-290.1061$ | $-0.026$ | $-291.0350$ |
| 甲基（ts0002） | $-247.55$ | $-247.5518$ | $-0.002$ | $-246.8596$ |

**一个完全独立的有限差分实现，在同一步长上落到基准记的那两个数上**（差 0.026 和
0.002 $\mathrm{cm^{-1}}$）。所以 §8.9 说"虚频那 1 $\mathrm{cm^{-1}}$ 的变化是 FD 截断、
不是我们改代码改坏了"**已被独立确认**，不再只是我们自己的推断。

#### 确认二：不投影刚体模的频率分析在驻点上带约 $10^{-4}$ 的转动污染

同一个 Hessian 两条路径的差异**随残余梯度走**：极小点（$f_{\max} = 7.3\times10^{-5}$）
$1.6\times10^{-9}$，两个 TS（$f_{\max} \approx 2\times10^{-2}$）$8.4\times10^{-5}$ 与
$2.1\times10^{-5}$ —— 差了 4 个数量级。

原因就是 §2.6g：**ASE 不投影刚体模**，而转动曲率正比于残余梯度
（$\dot x^\top H\dot x = \omega^2(\mathbf g\cdot x^\perp)$），于是在未投影的对角化里泄漏进物理模。
所以这个差异不是约定错误（符号结构与整条谱两侧一致），而是**我们的投影正在移除一项别人默认
不移除的误差**。对按搜索自身容差收敛的驻点，那一项约 $10^{-4}$。

这条是对照送来的、事先没料到的结论：它给投影方案本身提供了正面证据，而这正是 10.2 第 1 条
存在的理由。

#### 仍然没被覆盖的

负模**计数**两侧全部一致（0/1/1/1），符号结构一致 —— 但两侧跑的是同一个 PES。
**10.3 第 1 条一字不改：平滑偏差依旧不可见。**第 2、3、4 条（真 IRC、鞍点精化、内坐标约束
优化）仍需要 geomeTRIC；查曲率准不准仍然只有 DFT。
