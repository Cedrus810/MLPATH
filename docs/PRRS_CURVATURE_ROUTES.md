# PRRS 曲率判定的三条路线与架构结论

日期：2026-08-31
配套文档：[PRRS_STATUS.md](PRRS_STATUS.md)（实测状态）、[PRRS_ENGINEERING_PLAN.md](PRRS_ENGINEERING_PLAN.md)（计划稿）

状态稿第 7 节缺口 4（"大体系没有 $n_-$ 判定路径"）和缺口 1（"自由坐标未商掉"）不是两个独立的
工程待办，它们是同一件事在两个层面的表现：**PRRS 需要曲率信息，但不需要 Hessian 矩阵**。

本文记录三条路线。**已实测的是 A.12、A.13.4 与整个路线 C**；这些测量反过来修订了
A.6、A.7、A.10、A.11、A.13.4 与 B.8，修订处都已标注。

- **路线 A（工程）**：无矩阵曲率判定。目标从"算虚频"改成"检测负曲率方向"。可立即实现。
- **路线 B（理论）**：有限幅度响应算子 + Morse 指标。把 Hessian 降级为响应算子的
  $a \to 0$ 极限，从而 torsion 的弧/弦问题、对称商、自由转子、软模、虚频、TS 端点
  被同一个数学对象覆盖。
- **路线 C（估计量，含本文主要实测）**：$\mathcal L_a$ 就是对 PES 的
  $\operatorname{sinc}^2$ 低通滤波，$a$ 是分辨率而不是待消除的误差。以带内回归取代点估计，
  自带误差棒。附 MACE-OFF24 二阶导可靠性的四项无参考诊断。
- **架构结论（§E）**：Hessian 是 validation limit 而不是 primary observable ——
  因为 $a \to 0$ 是**唯一有谱定理**的地方。有限 $a$ 下不存在对称算子，所以
  $\kappa_a$ 连"本征值"都不是。含 persistent bottom modes 的可测定义与实测。

两条路线的关系：**A 是 B 的一个数值实现**。先落 A，B 提供命名与接口纪律，并指出 A
在有限幅度下的自然推广。

一句话：

> 不要"避免 Hessian 物理"，而要"避免 Hessian 矩阵"。

---

# 路线 A：无矩阵负曲率检测

## A.1 观测量的重新定义

虚频不是判据，它是判据的一个报告形式。虚频来自质量加权 Hessian 的负本征值：

$$H_m v_i = \lambda_i v_i, \qquad \lambda_i < 0 \;\Rightarrow\; \omega_i \text{ imaginary}.$$

PRRS 不需要知道 $-290.1\ \mathrm{cm^{-1}}$ 这个数字。它需要知道的是

$$\boxed{\lambda_1 < 0, \qquad \lambda_2 > 0}$$

即"只有一个不稳定方向"。$-250$ 还是 $-300\ \mathrm{cm^{-1}}$ 对 TS pool 的裁决没有影响。

## A.2 原语：Hessian-vector product，不建矩阵

已实现于 `src/prrs/runner.py:64` (`hessian_vector_product`)。中心差分力：

$$\boxed{Hv \approx -\frac{F(R + \epsilon v) - F(R - \epsilon v)}{2\epsilon}}$$

每个 $Hv$ 恰好 **2 次力评估**，与 $3N$ 无关。质量加权版本是
$H_m v = M^{-1/2} H (M^{-1/2} v)$，`direction_stiffness` (`runner.py:84`) 已经这样做。

投影版本（去掉平动/转动）：

$$H_{\rm proj} v = P H (P v), \qquad P = I - B^\top B$$

其中 $B$ 是 `_trivial_modes` 给出的刚体模正交基。**投影必须作用在输入和输出两侧**，
否则求解器会收敛到刚体零模而不是化学模式。

## A.3 标量版本：directional curvature

对单个方向 $v$（$\|v\| = 1$）：

$$\boxed{\kappa(v) = -\frac{[F(R + \epsilon v) - F(R - \epsilon v)] \cdot v}{2\epsilon} \approx v^\top H v}$$

$\kappa(v) < 0$ 就说明这是一个不稳定方向。代价 2 次力评估，**与体系大小无关**。

也可以从能量写（3 次能量评估，参考能量已知时 2 次）：

$$\kappa(v) \approx \frac{E(R + \epsilon v) + E(R - \epsilon v) - 2E(R)}{\epsilon^2}$$

两种写法的取舍是数值的，不是概念的：能量式分母 $\epsilon^2$ 放大相消误差，力式分母
$\epsilon$ 只放大一次。`torsion_response` (`runner.py:143`) 用能量式，理由是刚体旋转给出
**精确路径**，无链式法则歧义，这个好处压倒了条件数的损失。笛卡尔方向没有这个好处，
应当用力式。

## A.4 第一层筛选可以不做 Lanczos

TS pool 的入池判据可以直接降为

$$\boxed{\text{stationary point} + \text{observed negative directional response}}$$

而不是"先算 Hessian、再找虚频"。这就是当前 `_probe_minimum` (`runner.py:378`) 的方向，
但目前它用的是**随机种子位移 + 重新淬火**（乙醇势垒顶实测 2 个等价鞍点只抓到 1 个）。
改成显式的 $\kappa(v)$ 测量后：

- 不需要淬火，只需 2 次力评估；
- 结果是有符号的数，不是"跑掉了/没跑掉"的二值；
- 探针方向可以由 PRRS 已有的六家族提案给出（软模优先），不必随机。

## A.5 主动寻找最负方向（dimer-like）

当前方向 $v$，算 $\kappa(v)$，然后旋转 $v$ 求

$$v_{\min} = \arg\min_{\|v\| = 1} \kappa(v)$$

这就是最低 Hessian 本征模，但全程没有：构造 Hessian、对角化 Hessian、频率分析。
只有 **force probes**。

工程注意：旋转必须在投影子空间内进行（$v \leftarrow Pv$ 并重新归一），否则 dimer
会滑向刚体模。

## A.6 但一个负方向不能证明只有一个

$$\kappa(v_1) < 0 \;\Rightarrow\; n_- \ge 1$$

**不能**给出 $n_- = 1$。要把候选正式认证成一阶鞍点，必须看第二低模。但这依然不需要
完整 Hessian，只需要最低两个：

$$\boxed{\lambda_1, \lambda_2}$$

用 Lanczos（带重正交化）或 LOBPCG（块大小 3）作用在 $H_{\rm proj}$ 上。

### 认证条件与它的诚实边界（已用合成谱实测，见 A.12）

Ritz 值满足 Cauchy 交错 $\theta_k \ge \lambda_k$。**这个方向只帮一半的忙**：

- **下界是严格的**。每个负的 Ritz 值都证明存在一个负本征值，所以
  $\nu \ge \#\{\theta_k < -\lambda_{\rm tol}\}$ **无条件成立**。"至少有一个逃逸方向"是可以
  形式化断言的。
- **上界（"没有第二个负模"）拿不到无矩阵证书**。$\theta_2 > 0$ 不证明 $\lambda_2 > 0$，
  而残差界 $\|r\|$ 只保证**存在**某个本征值落在 $[\theta - \|r\|, \theta + \|r\|]$ 内，
  **不保证没有别的本征值更低**。

A.12 的合成谱实测把这一点钉死：在 $\lambda_1 = -0.50$、$\lambda_2 = -0.004$（真实
$\nu = 2$）上，带全重正交化的 Lanczos 在 $m = 20, 30, 40$ 处给出 $\theta_2 > 0$ 且
$\theta_2 - \|r_2\| > 0$ —— **残差修正后的判据会错误地认证 $\nu = 1$**，直到 $m = 60$ 才翻。
"count 在 $m \to m+10$ 下稳定"同样失效（$m = 10$ 到 $40$ 一直稳定在 1）。

实测中唯一没有误判的是带**余量因子**的判据：

$$\boxed{\begin{aligned}
&\|\mathbf F\| < F_{\rm stat} \\
&\theta_1 + \|r_1\| < -\lambda_{\rm tol} \\
&\theta_2 - \gamma\|r_2\| > -\lambda_{\rm tol}, \qquad \gamma \approx 10
\end{aligned}}$$

$\gamma$ 是**启发式余量，不是定理**。它在三个测试谱上没有误判，但没有证明它不会误判。
定理级的选项只有两个，都要额外东西：

| 手段 | 需要什么 | 给出什么 |
| --- | --- | --- |
| Cauchy 交错（免费） | 无 | $\nu$ 的**严格下界** |
| Kato–Temple | 一个先验的隔离间隙估计 | $\lambda_2$ 的真下界 |
| 惯性计数（Sylvester，$H - \sigma I$ 的 $LDL^\top$） | **完整矩阵** | 严格的 $\nu$ |

**结论（比原来弱，但是真的）**：无矩阵路线上，$\nu = 1$ 只能作为
**index-1 candidate + 余量证据**，不是认证。必须记录 $\theta_1, \theta_2, \|r_1\|, \|r_2\|, m$；
只写一个 $\nu = 1$ 就是断言。小分子保留 full projected Hessian 作 gold-standard
不是冗余，正是这个区别的直接后果。

判定树（注意 "1 negative" 一支的标签）：

```text
stationary?
    |
    yes
    |
lowest curvature modes
    |
    +-- 0 negative                        -> minimum candidate
    |
    +-- 1 negative, margin criterion met  -> index-1 candidate / TS candidate
    |
    +-- 1 negative, margin not met        -> undecided; raise m or fall back to full Hessian
    |
    +-- >=2 negative                      -> higher-order saddle (rigorous, reject)
```

注意不对称：**">=2 negative" 是严格的**（下界），"恰好 1" 不是。全程不需要转成
$\mathrm{cm^{-1}}$。

## A.7 代价

**注意 $\lambda_1$ 便宜、$\nu = 1$ 判定贵，两者差 2–4 倍。** A.12 的合成谱实测：
$\theta_1$ 在 $m \approx 15$ 已收敛到 5 位，而余量判据要到 $m = 40$（清晰间隙）或
$m = 60$（紧间隙）才通过。

| 路线 | 力评估次数（FD 版） | 随 $3N$ 增长 |
| --- | --- | --- |
| full projected Hessian（现状，`hessian_spectrum`） | $6N$ | 线性 |
| Lanczos 只要 $\lambda_1$ | $2m$，$m \approx 15$ → 约 30 | 不显式增长 |
| Lanczos 到可判 $\nu$ | $2m$，$m \approx 40\text{–}60$ → 80–120 | 不显式增长 |
| 单方向 $\kappa(v)$ 探针 | 2 | 不增长 |

所以 crossover 分两个：找反应坐标（$\lambda_1$ 与 $v_1$）约在 $6N > 30$，即 $N \gtrsim 5$
原子就划算；而判 $\nu$ 要 $6N > 100$，即 $N \gtrsim 17$ 原子。状态稿写的"约 15 原子"是后者
的量级，成立。**乙醇（9 原子，$6N = 54$）上 full Hessian 更便宜**，这不是妥协而是正解。

$m$ 依赖谱间隙而非维数，这是"不显式增长"的准确含义，且已实测确认：紧间隙
（$\lambda_2 = +0.004$）比清晰间隙（$\lambda_2 = +0.05$）多要 20 次迭代。近退化时应改块方法。
若能拿到解析 HVP（A.13），$2m$ 这一列要按实测的单次 HVP 成本重算，不再是"力评估次数"。

## A.8 接口改动：把虚频从核心层移出

现在（`runner.py:363`、`search.py:284`）核心诊断里带

```text
imaginary_wavenumbers_icm: [-290.1]
```

建议核心层只存：

```text
stationary        = true
negative_modes    = 1            # nu, Morse 指标
lambda_min        = ...          # eV/(A^2 amu)
lambda_second     = ...
residual_min      = ...          # 认证界，见 A.6
residual_second   = ...
unstable_direction = v           # 已有 unstable_mode 字段
method            = "hessian" | "lanczos" | "probe"
```

`imaginary_frequency_cm1` 降为**可选 reporting quantity**，用于与量化化学软件对话。
真正的核心信息是：

$$\boxed{\text{这个驻点有没有且只有一个逃逸方向？}}$$

## A.9 与已有 endpoint test 合并

状态稿 §4.2 已经在做

$$R^\ddagger + \epsilon v_- \to A, \qquad R^\ddagger - \epsilon v_- \to B$$

这比"有一个虚频"给出更多化学信息。完整 TS certification：

$$\boxed{\begin{aligned}
&\|\mathbf F\| < F_{\rm tol} \\
&\lambda_1 < -\lambda_{\rm tol} \\
&\lambda_2 > -\lambda_{\rm tol} \\
&R^\ddagger + \epsilon v_1 \to A \\
&R^\ddagger - \epsilon v_1 \to B
\end{aligned}}$$

全程不必计算真正的 vibrational frequencies。

## A.10 三档部署（已按 A.13 实测修订）

| 档 | 方法 | 理由 |
| --- | --- | --- |
| $N \lesssim 30$ 原子（含现有全部基准） | **解析 full projected Hessian**（`calc.get_hessian()`） | 实测约 $3N$ 次力评估等价、机器精度、给严格惯性计数。乙醇上比 Lanczos 认证便宜 3–4 倍 |
| $N \gtrsim 30$，要 $\nu$ | projected HVP + Lanczos，余量判据（A.6） | $2m \approx 80\text{–}120$ 次力评估；得到的是 index-1 candidate，不是认证 |
| 只要反应坐标，或超大体系初筛 | 单方向 $\kappa(v)$ / dimer-like 探针 | 2 次力评估拿一个负方向；$\lambda_1$ 单独收敛只要 $m \approx 15$ |

**与原稿的差别**：原稿把"小分子 full Hessian"当成便宜的调试手段、把 Lanczos 当主路线。
实测把这个顺序倒过来了 —— 解析 Hessian 在当前所有基准尺寸上既更便宜又更严格，Lanczos
的适用区间被推到还没有基准体系的尺寸上（A.13.5）。

流水线（不变，只是入口档位换了）：

$$\text{cheap probe} \to \text{negative curvature candidate} \to \text{full analytic Hessian (small) or 2-mode Lanczos (large)} \to \text{endpoint descent}$$

## A.11 落地清单

> **A.4 已落地**（`_probe_minimum` 改为方向明确的 $\kappa(v)$ 测量，扭转切向优先）、
> **A.5 已落地**（球面下降 `lowest_response_direction`）。A.6 的 Lanczos 仍未实现且已降级。
> 见 [PRRS_STATUS.md](PRRS_STATUS.md) §8.8.1、§8.8.3。

1. `runner.py`：在 `hessian_vector_product` 之上加 `projected_hvp(atoms, guard, v)`
   —— 两侧投影 + 质量加权，作为唯一的曲率入口。
2. `runner.py`：`lowest_curvature_modes(candidate, k=2)`，Lanczos/LOBPCG，返回
   $(\theta_i, \|r_i\|, u_i)$。
3. `config.py`：`minimum_check` 增加 `"lanczos"`（现有 `"hessian" | "probe" | "none"`
   在 `config.py:139` 校验），新增 `lanczos_max_iterations`、`lanczos_residual_tol`。
4. `_probe_minimum` 的随机位移改为方向明确的 $\kappa(v)$ 测量（A.4）。
5. 诊断字段按 A.8 改名；虚频保留但移到 reporting 区。
6. 回归门：乙醇三个 TS 上 `lanczos` 与 `hessian` 必须给出相同的 $\nu$，且
   $\lambda_1$ 相对偏差在容差内。这是 A 路线唯一可信的验收方式。

> **优先级已被 A.13.6 覆盖。**上面第 1–3 条仍然对，但顺序错了：实测表明第一件该做的事是
> 让 `hessian_spectrum` 改用 mace 的解析 Hessian（2 倍加速 + 机器精度），而 Lanczos 的
> 适用尺寸还没有对应的基准体系。按 A.13.6 的顺序执行。

## A.12 合成谱实测：Lanczos 判据

脚本 [`experiments/ritz_check.py`](experiments/ritz_check.py) 与 [`experiments/ritz_check2.py`](experiments/ritz_check2.py)（$N = 300$ 随机正交相似变换构造已知谱，Lanczos 全重正交化）。
**这是合成矩阵的测量，不是 MACE 上的测量** —— 它检验的是判据的逻辑，不是本体系的谱。

$\lambda_{\rm tol} = 10^{-3}$，三个谱，各自跑 $m = 5 \ldots 120$：

| 谱 | $\lambda_1$ | $\lambda_2$ | 真 $\nu$ | $\theta_1$ 收敛于 | 判据 A 首次 PASS | 判据 B 首次 PASS |
| --- | --- | --- | --- | --- | --- | --- |
| 清晰间隙 | $-0.50$ | $+0.050$ | 1 | $m \approx 15$ | $m = 5$（侥幸）/ 20 | $m = 40$ |
| 紧间隙 | $-0.50$ | $+0.004$ | 1 | $m \approx 15$ | $m = 5$（侥幸）/ 15 | $m = 60$ |
| 二阶鞍点 | $-0.50$ | $-0.004$ | **2** | $m \approx 15$ | **$m = 5, 20, 30, 40$ 全部误判 $\nu = 1$** | 从不 PASS（正确） |

判据 A：$\theta_2 - \|r_2\| > -\lambda_{\rm tol}$。判据 B：$\theta_2 - 10\|r_2\| > -\lambda_{\rm tol}$。

第三行是关键。二阶鞍点上，Lanczos 在 $m = 20/30/40$ 处报 $\theta_2 = +0.057/+0.052/+0.050$，
残差 $\|r_2\| = 5.0/1.3/2.6 \times 10^{-2}$，判据 A 通过 —— **一个真正的二阶鞍点会被当成 TS
收进 pool**。$\theta_2$ 直到 $m = 60$ 才翻成 $-0.004$。

同时否掉了另一个直觉判据：**"count 在 $m$ 增大时稳定" 不安全**。二阶鞍点上负 Ritz 值个数
从 $m = 10$ 到 $m = 40$ 一直是 1，$m = 60$ 才变 2。所以若要用稳定性作证据，步长必须是
$m \to 2m$ 级别，不是 $m \to m + 10$。

结论已写回 A.6（余量因子 $\gamma \approx 10$，且 $\nu = 1$ 只能是 candidate）与 A.7（两个
crossover）。

## A.13 能不能直接用解析二阶自动微分拿 $Hv$？

一句 `torch.autograd.functional.vhp` / `jax.jvp(jax.grad(E), (R,), (v,))` 就能给出**没有
$\epsilon$、没有截断误差**的精确 $Hv$。原则上对，而且比中心差分好。但落到本体系有几条硬约束。

### A.13.1 模型侧：可以，因为 MACE 的力本身就是 autograd

`mace/modules/utils.py:23` 的 `compute_forces` 就是
`torch.autograd.grad(E, positions, create_graph=training)`。所以只要以 `training=True`
调模型，二阶图就在，`torch.autograd.grad((g \cdot v).sum(), R)` 直接给精确 $Hv$。
mace 自己也已经这么用：`compute_hessians_vmap`（`utils.py:113`）用 `vmap` 把 $3N$ 个
单位向量的 vjp 批起来算完整解析 Hessian，`MACECalculator.get_hessian`
（`calculators/mace.py:770`）是它的入口。

**所以本体系其实连 FD 完整 Hessian 都不必要**：`get_hessian` 已经是解析的。
当前 `hessian_spectrum`（`runner.py:105`）用 $6N$ 次力评估重建它，是在解析导数已经可得
的情况下做数值微分。这是路线 A 之外的一个独立发现，见 A.13.5。

### A.13.2 但 ASE Calculator 接口把图丢了

`MACECalculator.calculate`（`calculators/mace.py:636`）调模型时传
`training=self.use_compile`，默认 `False` → `create_graph=False`，且结果 `.detach()`
成 numpy。**所以不是"加一行"就行**：`runner.py` 只持有 ASE Calculator，拿不到 torch 图。
需要一个新的后端原语，绕过 ASE 直接和 torch 模型对话。

代价是三条设计约束都被触到：

1. **搜索核心目前不依赖积分器，也不依赖 torch**（状态稿 §1）。解析 HVP 会把 torch 引进
   曲率路径。可接受的形式是：`calculators.py` 暴露**可选**的 `exact_hvp(atoms, v)`
   能力，探测不到就退回 FD，与 §3.7"存在/生效必须由执行证实"一致。
2. **`GuardedCalculator` 的失败关闭门控按每次力评估计数**（`reliability.py`）。一次解析
   HVP 不是一次力评估，能量/力的合理性检查无处挂。需要给 HVP 定义自己的门控
   （至少：$Hv$ 有限、$\|Hv\|$ 上界、对称性抽检）。
3. **交叉校验必须保留**。解析 vs FD 的一致性正是状态稿 §2.1 那种验收；不做这个校验就
   引入解析路径，等于换掉一个已验证组件。

### A.13.3 JAX 这一条不适用

MACE-OFF24 是 torch checkpoint（`torch.load` 的 `nn.Module`，`get_hessian` 走
`torch.vmap`）。用 JAX 要重实现模型，不是写一行。**torch 侧的对应写法**：

| 写法 | 模式 | 对 MACE 的可用性 |
| --- | --- | --- |
| `torch.autograd.grad(g·v, R)`（手写二阶） | reverse-over-reverse | **推荐**，走的正是 mace 自己 `compute_hessians_vmap` 用的路径 |
| `torch.autograd.functional.vhp(E, R, v)` | reverse-over-reverse | 可用，但要把 batch/邻域构造包进一个纯函数 `E(R)`，比手写麻烦 |
| `torch.func.jvp(torch.func.grad(E), (R,), (v,))` | forward-over-reverse | 内存更省，但依赖 forward-mode 对 e3nn/scatter 的完整覆盖，**未验证** |

$E$ 是标量，Hessian 对称，所以 `vhp` 与 `hvp` 数学上同一个东西；差别只在实现模式与内存。

### A.13.4 已实测：MACE-OFF24 上双反向可行，精度到机器精度，成本与 FD 相同

脚本 [`experiments/hvp_autograd_probe.py`](experiments/hvp_autograd_probe.py)，
输出 [`experiments/hvp_autograd_probe.out`](experiments/hvp_autograd_probe.out)。
乙醇（9 原子，`rattle(0.05)` 后的非驻点，$f_{\max} = 6.61$ eV/Å），
MACE-OFF24_medium、float64、CUDA。

**双反向跑通，并且与 ASE 路径是同一个物理量**：

```
energy via direct model call: dE vs calculator = 0.000e+00
force  via direct model call: max|dF| vs calculator = 1.110e-15 eV/A
```

**与 mace 自己的解析完整 Hessian 一致到机器精度**（三个探针方向）：

```
asymmetry of get_hessian: max|H-H^T| = 1.066e-14
max|H@v - Hv_exact| = 4.4e-15 / 1.8e-15 / 7.1e-15
symmetry  v0.H v1 - v1.H v0 = 2.665e-15
```

**FD 的截断误差是 $O(\epsilon^2)$，实测确认**（最大绝对偏差，eV/Å²，$|Hv|_\infty \approx 20/13/40$）：

| 探针方向 | $\epsilon = 10^{-3}$ | $\epsilon = 10^{-2}$ | $\epsilon = 5\times10^{-2}$ |
| --- | --- | --- | --- |
| v0 | 1.44e-05 | 1.44e-03 | 3.54e-02 |
| v1 | 9.16e-06 | 9.16e-04 | 2.27e-02 |
| v2 | 7.75e-05 | 7.75e-03 | 1.92e-01 |

$\epsilon$ 每涨 10 倍误差涨 100 倍，干净的二阶。

> **这条要单独记住**：当前 `minimum_check_step_A = 0.01`，而
> `minimum_check_eigenvalue_tol = 1e-3`。在这个步长下 FD 的 $Hv$ 误差是
> $1.4 \times 10^{-3} \sim 7.8 \times 10^{-3}$ eV/Å² —— **比本征值容差还大**。
> 质量加权会把重原子分量按 $1/\sqrt{m_i m_j}$ 压小，但 H 原子（$m = 1$）不压。
> 所以现在判 $|\lambda| \lessapprox 10^{-3}$ 的模式（正是 §3.5 三级判据里 `free` 与
> `negative curvature` 的分界）**处在 FD 噪声里**。这在乙醇的三个 TS 上必须复测：
> 上面是一个非驻点上的单点测量，$\epsilon^2$ 的标度是普适的，系数不是。

**成本与显存**（20 次力评估 / 10 次 HVP 取平均；机器 load average 约 18，绝对值有抖动，
比值可信）：

| 量 | 墙钟 | 折合力评估 |
| --- | --- | --- |
| 一次力评估 | 35.4 ms | 1 |
| **一次解析 $Hv$（含重建图）** | **69.5 ms** | **1.96** |
| 一次 FD $Hv$（2 次力） | 70.8 ms | 2 |
| **完整解析 Hessian（`get_hessian`）** | **946 ms** | **26.7** |
| 完整 FD Hessian（$6N = 54$） | 1912 ms | 54 |

两个结论：

1. **解析 $Hv$ 和 FD $Hv$ 成本相同（1.96 vs 2.00 次力评估），精度差 12 个数量级。**
   没有取舍 —— 解析版是免费的精度。A.7 那张表里"力评估次数"这一列因此对两种实现都成立。
2. **完整解析 Hessian 比完整 FD Hessian 快 2 倍**（26.7 vs 54），因为 `vmap` 把 $3N$ 个 vjp
   批起来，单行只要约 1.0 次力评估而不是 1.96。

显存峰值 **117.9 MiB**（9 原子）。所以 A.13 早先写的"内存是真实风险"在这个尺寸上**不成立**，
11 GB 卡上离上限很远。但这只是 $N = 9$ 的一个点：解析 HVP 的图与 `get_hessian` 的
$O((3N)^2)$ 输出都随 $N$ 增长，`compute_hessians_vmap` 的 `chunk_size` 还是硬编码的
（`utils.py:132`）。**大体系上的显存标度未测**，不要外推。

### A.13.4b 一个必须一起处理的单位陷阱

`MACECalculator.calculate` 在 `mace.py:657` 起对 energy/forces/stress 逐项乘
`energy_units_to_eV` / `length_units_to_A`。**`get_hessian`（`mace.py:770`）一项都不乘。**
MACE-OFF 的两个因子都是默认 1.0，所以上面的机器精度一致性是真的，但它同时是
状态稿 §2.4 那一类"默认值恰好对"的陷阱。改用解析 Hessian 时必须自己补换算，
并按 §3.7 用执行证实（对一个已知谱交叉校验），不能信参数。

### A.13.5 由此产生的 crossover 修订

第 2 条结论直接改掉路线 A 的经济性论证：竞争对手不是 $6N$ 次力评估的 FD Hessian，而是
约 $1.0 \times 3N$ 次力评估的**解析** Hessian。

| 方法 | 折合力评估（外推自 $N=9$ 实测） |
| --- | --- |
| 解析完整 Hessian | $\approx 3N$ |
| Lanczos 只要 $\lambda_1$（$m \approx 15$） | $\approx 30$ |
| Lanczos 到可判 $\nu$（$m \approx 40\text{–}60$） | $\approx 80\text{–}120$ |

于是：

- **找反应坐标**（$\lambda_1, v_1$）：$3N > 30$，即 $N \gtrsim 10$ 原子 —— Lanczos 划算。
- **判 $\nu$**：$3N > 100$，即 $N \gtrsim 33$ 原子。

**比 A.7 用 FD 基线算出的 17 原子晚了一倍。**乙醇（9 原子）上解析完整 Hessian
只要约 27 次力评估，Lanczos 认证要 80–120 —— **full Hessian 便宜 3–4 倍**，而且是严格的
惯性计数而非余量证据（A.6）。

### A.13.6 落地清单修订（覆盖 A.11 的优先级）

1. **先做零成本的那一步，它和 Lanczos 无关**：`hessian_spectrum`（`runner.py:105`）在
   MACE 后端可用时改调 `calc.get_hessian()`。实测 2 倍加速、精度从 $O(\epsilon^2)$ 变成机器
   精度，且顺手解决上面那条"FD 噪声大于 `eigenvalue_tol`"的问题。FD 版保留为交叉校验与
   非 MACE 后端实现。**这一步的收益比整个 Lanczos 路线更确定。**
2. `calculators.py` 增加可选能力 `exact_hvp` / `exact_hessian`，探测方式是**试算一次并比对
   FD**，不是查版本号（§3.7）。
3. `lowest_curvature_modes` 的 HVP 入口做成 `exact | fd` 可替换；两者在乙醇三个 TS 上必须
   给出相同 $\nu$ 与相符的 $\lambda_1$。
4. Lanczos 的优先级**下调**：它的适用区间被推到 $N \gtrsim 33$，而当前唯一真实基准是 9 原子。
   在有那个尺寸的基准体系之前，实现它是在没有验收对象的情况下写代码。
5. 重测清单：乙醇三个 TS 上 `get_hessian` vs `hessian_spectrum` 的 $\nu$ 与 $\lambda_1$；
   以及 FD 误差系数在**驻点**上的实际大小（决定 `minimum_check_step_A` 该不该改）。

### A.13.7 仍未验证

- **forward-over-reverse 路线**（`torch.func.jvp(torch.func.grad(E), ...)`）没测。它内存更省，
  但依赖 forward-mode 对 e3nn/scatter 的完整覆盖。上面所有数字都来自
  reverse-over-reverse。
- **大体系显存与 `get_hessian` 的实际标度**。见 A.13.4 末段。
- **`GuardedCalculator` 对 HVP 的门控**（A.13.2 第 2 条）仍未设计。一次解析 HVP 不是一次力
  评估，现有失败关闭逻辑挂不上去。
- **openmm-ml 路径下的解析 HVP**。`openmmml` 的 `PythonForce` 只回传力，二阶图不出边界；
  上面走的是 `mace.calculators.MACECalculator`。

---

# 路线 B：有限幅度响应算子 + Morse 指标

路线 A 解决"怎么算"。路线 B 解决"算的是什么"，并且让状态稿里六处分别打的 patch
变成一个定义。

## B.1 坐标：对称约化的 configuration manifold

不要把体系看成裸的 $R \in \mathbb R^{3N}$。有意义的是商掉整体平移、转动与离散对称等价后的空间：

$$\mathcal X = \mathcal C / \bigl(SE(3) \times \Gamma\bigr)$$

$\Gamma$ 是原子置换/分子对称群。一个状态是 $x = [R] \in \mathcal X$。

状态稿 §3.6 列的六处"商掉冗余自由度"（chemical key 的自同构不变哈希、对称感知 RMSD、
角度枚举正则色去重、对称转子基本域、reservoir 多样性度量、立体化学签名）**都是这一个
定义的实例**。它们已经在实际做这件事，只是没有共同的名字。

## B.2 Perturbation：retraction，不是 $R + \delta R$

定义 retraction

$$\mathscr R_x(a v), \qquad v \in T_x \mathcal X$$

满足

$$\mathscr R_x(0) = x, \qquad \left.\frac{d}{da}\mathscr R_x(a v)\right|_{a=0} = v$$

直觉：从 $x$ 沿物理允许的方向 $v$ 走**有限**距离 $a$。

| 实现 | 地位 |
| --- | --- |
| 桥键 torsion 的刚体 fragment 旋转 | 几乎精确的 retraction（实测 `target_error` $\le 4\times10^{-16}$，键长变化 $\le 2.2\times10^{-16}$ Å） |
| continuation（$R_0 \to R_1 \to \cdots \to R_n$，`MAX_STEP = 0.05`） | 数值 retraction |
| 单纯 $R + a v$ | 只是欧氏空间的一阶近似 |

这直接给出状态稿 §3.3 "走弦不走弧" bug 的数学修复，也解释了 §3.2
"target fidelity $\ne$ geometric fidelity"：$R + a v$ 未必留在合理流形上（2.09 rad 扭转把
O–H 从 0.957 拉到 2.346 Å），而 $\mathscr R_x(a v)$ 按定义就是合法的 perturbation delivery。

## B.3 核心对象：finite-amplitude response operator

$$\boxed{\mathcal L_a(x) v = \frac{\operatorname{PT}_+ \nabla E(\mathscr R_x(av)) - \operatorname{PT}_- \nabla E(\mathscr R_x(-av))}{2a}}$$

$\operatorname{PT}_\pm$ 把两个扰动点的梯度输运回 $T_x \mathcal X$。在当前笛卡尔/投影实现里
$\operatorname{PT}$ 就是投影回来（$P$，见 A.2）。关键性质：

$$\boxed{\lim_{a \to 0} \mathcal L_a(x) v = \operatorname{Hess} E(x)\, v} \qquad\Longrightarrow\qquad \boxed{H = \mathcal L_0}$$

Hessian 不是算法从一开始就必须拥有的东西，它是响应算子的无穷小极限。
`hessian_vector_product` 的 docstring 已经写着这句话，B 路线只是把它提为定义。

## B.4 标量形式与尺度参数

$$\boxed{\kappa_a(x; v) = \langle v, \mathcal L_a(x) v\rangle}$$

能量形式：

$$\kappa_a(x; v) \approx \frac{E(\mathscr R_x(av)) + E(\mathscr R_x(-av)) - 2E(x)}{a^2}$$

$a \to 0$ 时 $\kappa_0 = v^\top H v$（即 A.3）。有限 $a$ 下它包含的是
**finite-amplitude nonlinear response**，不只是 harmonic curvature。

于是整个 PRRS 统一成一个尺度参数：

$$\boxed{a: 0 \longrightarrow \text{finite} \longrightarrow \text{basin crossing}}$$

$$\text{local curvature} \longrightarrow \text{anharmonic response} \longrightarrow \text{reaction}$$

同一个 $\kappa_a$ 在 $a$ 小的时候是 §3.5 的 mode-resolved 收敛判据里的 $k_v$，在 $a$ 大的
时候是搜索本身。状态稿里 `direction_stiffness`（$a$ 小）与幅度细化（$a$ 大）是同一个测量
的两端，不是两个模块。

## B.5 TS 不叫"有一个虚频"，叫 Morse 指标

驻点 $\nabla E(x^\ddagger) = 0$。传统条件是 Hessian 有一个负本征值。数学上的名字是
Morse index：

$$\nu(x) = \#\{\lambda_i < 0\}$$

$$\nu = 0 \Rightarrow \text{minimum}, \qquad \boxed{\nu = 1 \Rightarrow \text{first-order saddle}}, \qquad \nu > 1 \Rightarrow \text{higher-order saddle}$$

PRRS 核心只需要 $\nu(x^\ddagger) = 1$。虚频只是为了和量化化学软件交流才转换出来的
reporting quantity（同 A.8）。

## B.6 应该用 Morse–Bott，不是纯 Morse —— 这就是缺口 1 的修法

> **已落地**：`chemistry.free_aligned_symmetric_rmsd` + `network.free_bonds_of`，
> 缺口 1 关闭。$\mathcal N$ 由淬火报告的 `tier == "free"` 给出并随几何走，因为这个判断是
> 能量性的、不能从图上猜。见 [PRRS_STATUS.md](PRRS_STATUS.md) §8.8.2。

普通 Morse 理论假设临界点孤立，即 $\ker H = 0$。真实分子可能有近似自由转子/连续中性方向。
状态稿缺口 1 正是这个：若某扭转真自由（$k \to 0$），只在该坐标上不同的结构由零势垒路径
相连、应属同一盆地，但 `symmetric_rmsd` 会判为不同 microstate，于是自由转子生成无穷多假构象。

Morse–Bott 的分解：

$$T_x \mathcal X = \mathcal N \oplus \mathcal U \oplus \mathcal S$$

$$\mathcal N: \text{neutral/free}, \qquad \mathcal U: \text{unstable}, \qquad \mathcal S: \text{stable}$$

**不要求** $\dim \mathcal N = 0$。TS 只要求

$$\boxed{\nabla E(x^\ddagger) \approx 0, \qquad \dim \mathcal U = 1}, \qquad \dim \mathcal N \ge 0$$

这比人为规定 $|k| < 0.01 \Rightarrow$ free（状态稿 §3.5 的三级判据）更有理论出处：
free 一级不是容差 hack，它是 $\mathcal N$。

**并且给出缺口 1 的具体修法**：$\mathcal X$ 的商不止 $SE(3) \times \Gamma$，还要沿 $\mathcal N$
对齐 —— microstate 去重前先在 $\mathcal N$ 上取商（RMSD 比较前沿自由坐标对齐），与 §3.6
第 2 条同源。乙醇不触发（甲基 $k = 0.625$），所以这条目前只能靠解析测试覆盖。

## B.7 怎样不用 Hessian 找 $\mathcal U$

$$v_1 = \arg\min_{\|v\|=1} \kappa_a(x; v)$$

若 $\kappa_a(v_1) < 0$，找到一个 unstable response direction。然后在正交补里继续：

$$v_2 = \arg\min_{\|v\|=1,\ v \perp v_1} \kappa_a(x; v)$$

若 $\kappa_a(v_1) < 0$ 且 $\kappa_a(v_2) > 0$，则 $\boxed{\dim \mathcal U = 1}$
—— **在已收敛的响应子空间内**成立（这就是 A.6 那个残差界的几何说法）。

这正是 Lanczos / dimer / Rayleigh–Ritz 背后的数学。区别在于：

$$\text{Lanczos 是 } \mathcal L_0 \text{ 最低响应模式的 solver，不是理论本身。}$$

理论层不应该出现"Lanczos Hessian"这种命名。

## B.8 scale-dependent response index

$$\boxed{\nu_a(x) = \#\{\text{negative modes of } \mathcal L_a(x)\}}$$

$a \to 0$ 时 $\nu_a \to \nu_{\rm Morse}$。有限 $a$ 下 $\nu_a$ 描述的是**这个尺度下**体系有几个
不稳定响应方向。于是对一个候选可以观察 $\nu_a$ 对 $a$ 的谱：

| $a$ | $\nu_a$ |
| --- | --- |
| 0.01 | 1 |
| 0.05 | 1 |
| 0.10 | 1 |
| 0.30 | 2 |

读法：局部是干净的一阶鞍点，但更远处出现第二个 nonlinear escape channel。
这比一次传统 Hessian frequency analysis 信息更丰富，且**只用 PRRS 已有的原语**
——同一个 $\kappa_a$ 换 $a$。

> **上表是猜测，并且方向被实测证否。**C.5 的 D2 在乙醇 syn 鞍点上实测到 $\nu_a$ 随 $a$
> **减小**（$a \lesssim 0.05$ 时 $\nu_a = 1$，$a \gtrsim 0.1$ 时 $\nu_a = 0$），因为该不稳定
> 模是有界周期坐标，走远会翻过邻近极小再爬升。$\nu_a$ 对 $a$ **非单调**，两个方向都可能，
> 不能假设。$\kappa_a$ 的滤波解释见路线 C。

## B.9 TS certification 的最终形式

不再是

$$\text{optimize} \to \text{Hessian} \to \text{frequency} \to \text{count imaginary frequencies}$$

而是

$$\boxed{x^\ddagger \to \begin{cases}
\|\nabla E\| \approx 0 \\
\dim \mathcal U = 1 \\
+ v_u \to A \\
- v_u \to B
\end{cases}} \qquad v_u = \arg\min_v \kappa_0(x^\ddagger; v)$$

一个 TS 有三个数学属性：

$$\boxed{\text{stationarity} + \text{index} + \text{connectivity}}$$

$$\nabla E = 0, \qquad \dim \mathcal U = 1, \qquad \Phi(+v_u) = A,\ \Phi(-v_u) = B$$

$\Phi$ 是 quench/basin map。状态稿 §4.2 的鞍点下降就是 $\Phi$，且已经明确记录
"TS 连接是无向结构性证据，不自动生成有向边" —— 这条纪律在 B 的语言里是：
$\Phi$ 给出 $\mathcal X$ 上的无向关联，有向反应边需要另外的证据。

## B.10 统一对象

$$\boxed{\mathcal R_x(v, a, t)}$$

在状态 $x$、沿 perturbation $v$、幅度 $a$、经过时间 $t$ 后的响应。

| 极限 | 内容 |
| --- | --- |
| $a \to 0,\ t \to 0$ | linear response / Hessian |
| 有限 $a$，短 $t$ | anharmonic susceptibility |
| 更大 $a, t$ | basin escape |
| 最终 | reaction |

所以

$$\boxed{\text{Hessian} \subset \text{PRRS response theory}}$$

而不是 $\text{PRRS} + \text{Hessian module}$。

$t$ 这一维对应状态稿 §3.1 的位移版 vs 动量版之争（"kick 注入的能量让体系在自由段游走"）：
位移版是 $t \to 0$ 的 $\mathcal R$，动量版是 $t > 0$ 的 $\mathcal R$。该节的推论"位移版对浅势垒优于
动量版"在这个记号下是一个关于 $t$ 的陈述，仍未实现、未实测。

## B.11 骨架两句话

> **PRRS studies the finite-amplitude response operator on the symmetry-reduced
> molecular configuration manifold.**

> **The Hessian is merely its infinitesimal limit, and a transition state is identified
> by an index-one unstable response subspace plus two-sided basin connectivity.**

这样 torsion 的弧/弦问题、对称商、free rotor、soft mode、虚频、TS endpoint、finite
perturbation 落在同一个体系里。

---

# 路线 C：有限幅度响应 = 对 PES 的尺度滤波

路线 C 不是第三条并列路线，而是把路线 B 的 $\mathcal L_a$ 从"理论组织"变成**主估计量**，
并给出它为什么能替代 $a\to0$ 的精确理由。

起点是一句观察：**有限幅度响应算子实际上对 PES 做了一次低通滤波**。这不是类比，是可以
写出传递函数并实测验证的。

## C.1 传递函数：$\kappa_a$ 是 PES 曲率谱的低通滤波

沿方向 $v$ 取 PES 的一个空间傅里叶分量 $E \sim e^{ikx}$（$k$ 是沿 $v$ 的空间波数）。

**能量形式**（对称二阶差分）：

$$\frac{e^{ika} + e^{-ika} - 2}{a^2} = -\frac{4}{a^2}\sin^2\frac{ka}{2}
= -k^2 \left[\frac{\sin(ka/2)}{ka/2}\right]^2$$

精确二阶导给 $-k^2$，所以

$$\boxed{\kappa_a^{E}: \quad T_a^{E}(k) = \operatorname{sinc}^2\!\left(\frac{ka}{2}\right)}$$

**力形式**（对称一阶差分的投影）：

$$\frac{ik\,e^{ika} - ik\,e^{-ika}}{2a} = -k\,\frac{\sin(ka)}{a} = -k^2\,\frac{\sin(ka)}{ka}$$

$$\boxed{\kappa_a^{F}: \quad T_a^{F}(k) = \operatorname{sinc}(ka)}$$

两者都是低通：$T \to 1$ 当 $ka \ll 1$，$T \to 0$ 当 $ka \gg 1$。**截断波数约 $2/a$**，
即幅度 $a$ 直接就是分辨率。$T^E$ 的首个零点在 $ka = 2\pi$，即**波长恰为 $a$ 的分量被完全消掉**。

所以"给体系一个有限幅度的扰动，看它怎么响应"这件事，等价于**先把 PES 带限到波数
$\lesssim 2/a$，再取它的曲率**。$a$ 不是需要消掉的误差参数，$a$ 是所选的分辨率。

## C.2 两个形式的重要差别：只有能量形式是非负滤波

$$T_a^{E} = \operatorname{sinc}^2 \ge 0, \qquad T_a^{F} = \operatorname{sinc} \text{ 会变号}$$

$T^F_a(k) < 0$ 当 $ka \in (\pi, 2\pi)$ —— **力形式会把波长在 $a$ 与 $2a$ 之间的分量反号**，
能量形式永远不会。

原稿（A.3）以条件数为理由推荐力形式（分母 $\epsilon$ 而非 $\epsilon^2$）。**在滤波这个目的下这个
推荐反了**：C.5 的实测表明数值噪声根本不是本体系的问题，于是条件数的论证失效，而非负性
的论证保留。**多尺度估计量应当用能量形式**；力形式留给单点 $\kappa(v)$ 快筛（A.3）与
Lanczos 的 $Hv$（那里要的是 $H$ 本身，不是滤波）。

领头阶展开给出一个可直接检验的预言：

$$T^E_a = 1 - \frac{(ka)^2}{12} + \cdots, \qquad T^F_a = 1 - \frac{(ka)^2}{6} + \cdots$$

$$\Longrightarrow \boxed{\frac{\kappa_a^F - \kappa_0}{\kappa_a^E - \kappa_0} \to 2}$$

## C.3 实测确认这个因子就是 2

用 C.5 的鞍点数据（`experiments/sdr_ts.out`）直接算 $\text{dev}_F/\text{dev}_E$：

| 模式 | $a{=}10^{-2}$ | $2\times10^{-2}$ | $5\times10^{-2}$ | $10^{-1}$ | $2\times10^{-1}$ | $3\times10^{-1}$ |
| --- | --- | --- | --- | --- | --- | --- |
| 不稳定模 $\lambda_1 = -0.2960$ | 2.004 | 2.001 | 1.996 | 1.983 | 1.934 | 1.857 |
| 最低正模 $+0.2550$ | 1.988 | 2.012 | 2.002 | 1.996 | 1.985 | 1.967 |
| 次低正模 $+0.6540$ | 2.333 | 1.923 | 1.987 | 1.990 | 1.976 | 1.956 |

三个模式、跨一个半数量级的幅度，比值都是 **2.00**，大 $a$ 端随高阶项进入缓慢偏离。
（$a \le 3\times10^{-3}$ 的列因为原始表只印到 5 位小数，比值是打印精度伪影。）

**所以滤波图像不是解释性说法，它是被数据定量确认的。**

## C.4 这个滤波扣掉什么、扣不掉什么

**能扣掉**：PES 的高空间频率成分。如果 MLP 的误差是高频涟漪（这是通常的担心），
$\operatorname{sinc}^2$ 确实按 $1/(ka)^2$ 压制它。你的推理这一步成立。

**扣不掉的三件事，必须写清**：

1. **平滑偏差**。一个曲率整体偏硬 10% 的模型，通过任何低通滤波都还是偏硬 10%。滤波作用在
   频率上，偏差在幅度上。**这是本体系真正剩下的风险**，而且对所有无参考诊断都不可见。
2. **真实的高频物理**。滤波不区分"模型涟漪"与"真实非谐结构"，波长 $\approx a$ 的分量一律
   消掉。所以它不是去噪器，**它是一个带宽选择**。诚实的报告是"在分辨率 $a$ 下的曲率"，
   不是"曲率"。
3. **本体系的高频成分本来就不多**（C.5 的 D3）。所以这个滤波在这里是**零成本的保险**，
   不是对已证实问题的修复。

## C.5 实测：MACE-OFF24 的二阶导在乙醇上看不出噪声

四个**无需参考数据**的诊断，脚本
[`experiments/second_deriv_reliability.py`](experiments/second_deriv_reliability.py)（极小点）与
[`experiments/sdr_ts.py`](experiments/sdr_ts.py)（鞍点），输出为同名 `.out`。

**鞍点独立复现（限于能量）**：固定二面角 $\text{C–C–O–H} = 0°$ 的约束极小化给出
`fmax = 4.0e-05 eV/Å`（**约束下的** fmax），$E = -4221.556394$ eV，相对极小点
$-4221.603394$ eV 即**势垒 47.0 meV，恰好一个负模 $\lambda_1 = -0.295964$**。
状态稿 §5 记的是 $+47.1$ meV 的 syn 势垒，**势垒能量复现**。

**但曲率没有复现，也不该期待它复现**：$\lambda_1 = -0.295964$ 换算是
$-283.69\ \mathrm{cm^{-1}}$，而 §5 记的是 $-290.1\ \mathrm{cm^{-1}}$。原因是约束极小化的
驻点不是无约束驻点 —— 约束把一个方向排除在收敛判据之外，所以那个 `4.0e-05` 不覆盖它。
在真实搜索找到的同一个鞍点上（`fmax` 按搜索自己的容差），解析 Hessian 给
$\lambda_1 = -0.311480 \leftrightarrow -291.035\ \mathrm{cm^{-1}}$（状态稿 §8.8）。
本节 D2/D3/C.3/C.6 的所有相对结论都用同一个几何内部一致，因此不受影响；
**受影响的只有"复现了 $-290.1\ \mathrm{cm^{-1}}$"这个说法，它不成立。**

### D1 不变性误差地板（严格零本征值的实测大小）

对任何精确不变的 PES，3 个平动 + 3 个转动模的本征值**严格为零**。它们的实测大小是解析
Hessian 的硬误差下界，**不需要任何参照数据**。

| 几何 | 平动 | 转动 | $|\lambda_{\rm trivial}|_{\max}$ | 最低物理 $|\lambda|$ | `eigenvalue_tol` |
| --- | --- | --- | --- | --- | --- |
| 极小点 | $\sim10^{-15}$ | $-3.7$e$-6$, $-4.4$e$-6$, $-2.6$e$-5$ | $2.5\times10^{-5}$ | $0.216$ | $10^{-3}$ |
| 鞍点（约束几何） | $\sim10^{-15}$ | — | $2.9\times10^{-6}$ | $0.296$ | $10^{-3}$ |

平动到机器精度（能量只依赖坐标差，这是恒等式）；转动到 $10^{-6}\text{–}10^{-5}$。

> **这两个转动数字后来被证否为"模型误差"。**$E(R(s)x) = E(x)$ 两次求导给
> $\dot x^\top H\dot x = \omega^2(\mathbf g\cdot x^\perp)$，所以转动商正比于残余梯度，只在驻点为零。
> 减掉这一项后解析路径的残差是 $\sim10^{-16}$（e3nn 等变性给的严格不变性）。
> 详见状态稿 §2.6g；本节其余结论不依赖这两个数字的归因。

### D2 $\kappa_a$ 的尺度谱（鞍点，物理模）

$\kappa_a^F$，单位 eV/(Å² amu)，$a$ 单位 Å$\sqrt{\text{amu}}$；`dev%` 相对 $\lambda$：

| 模式 | $\lambda$ | $10^{-3}$ | $10^{-2}$ | $2\times10^{-2}$ | $5\times10^{-2}$ | $10^{-1}$ | $2\times10^{-1}$ | $3\times10^{-1}$ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 不稳定 | $-0.2960$ | 0.01% | 0.98% | 3.9% | 24% | **95%** | 354% | 705% |
| 最低正 | $+0.2550$ | 0.00% | 0.13% | 0.54% | 3.4% | 13% | 53% | 115% |
| 次低正 | $+0.6540$ | 0.00% | 0.00% | 0.02% | 0.12% | 0.48% | 1.9% | 4.0% |

**全部平滑、单调，没有小 $a$ 抖动。**不稳定模的 $\kappa_a$ 在 $a \approx 0.1$ 处**变号**
（$-0.0137 \to +0.319$）—— 这是真实物理：沿 OH 扭转从 syn 顶两侧走远会翻过
gauche$^\pm$ 极小再爬升，对称二阶差分因此转正。

**这是路线 B §B.8 的 $\nu_a$ 谱第一次被实测**，并且**修正了原稿的猜测方向**：B.8 那张示意表
猜 $\nu_a$ 随 $a$ 增大（远处出现新逃逸通道），实测这里是**减小**（$\nu_a = 1$ 当
$a \lesssim 0.05$，$\nu_a = 0$ 当 $a \gtrsim 0.1$），因为该模是有界周期坐标。两种情形都可能，
$\nu_a$ 对 $a$ **非单调**，不能假设方向。

### D3 Hessian 作为函数有多糙

$\|H(R + \delta w)v - H(R)v\| / \delta$，$v$ 为不稳定模，$\|Hv\| = 0.5252$：

| $\delta$ (Å) | $10^{-5}$ | $10^{-4}$ | $10^{-3}$ | $10^{-2}$ | $10^{-1}$ |
| --- | --- | --- | --- | --- | --- |
| $\|\Delta Hv\|/\delta$ | 22.6530 | 22.6510 | 22.6304 | 22.4229 | 20.2684 |

**跨两个数量级的 $\delta$ 稳定到 4 位** —— 三阶导存在且平滑，**没有噪声**。这是反对"MLP 二阶导
粗糙"的最直接证据。

同时它说了另一件事：相对量是 **43 / Å**，即 $H$ 沿这个方向每 0.01 Å 变化约 43%。
**"某一点的 Hessian"本身就是脆弱的描述符**，与它算得准不准无关。这是支持多尺度的
物理论证，独立于任何可靠性讨论。

### D4 float32 vs float64

| 几何 | $\lambda_1$ (f64) | $\lambda_1$ (f32) | 全谱 $\max\|\Delta\lambda\|$ | f32 误差地板 |
| --- | --- | --- | --- | --- |
| 极小点 | — | — | $1.9\times10^{-5}$ | $2.4\times10^{-5}$ |
| 鞍点 | $-0.295964$ | $-0.295966$ | $1.5\times10^{-5}$ | $3.4\times10^{-6}$ |

**float32 不改变负模的符号或大小。**部署精度不是风险来源。

### D5 源码层面排除的一个伪影

MACE 的 `PolynomialCutoff`（$p = 6$，`modules/radial.py:113`）在 $r_{\max}$ 处
$f = f' = f'' = 0$（$1-28+48-21$、$-168+336-168$、$-840+2016-1176$ 均为 0），
所以**邻域截断不破坏 Hessian 连续性**，邻居进出不产生二阶跳变。
但 $f'''(1) = -336 \ne 0$：三阶导在截断处跳变，只影响有限 $a$ 的高阶项。
加载的 MACE-OFF24_medium 中与光滑性相关的子模块只有 `PolynomialCutoff`
（无 `ZBLBasis` / `AgnesiTransform` / `SoftTransform`）。

### 结论

**"MLP 二阶导不可靠"在 MACE-OFF24 + 乙醇 + 极小点/鞍点上没有得到证实。**
误差地板 $10^{-6}$–$10^{-5}$，三阶导平滑到 4 位，float32 都够。

但**无参考诊断只能查不一致，不能查不准**。一个平滑但整体偏移的 Hessian 会完美通过
D1–D5。所以剩余风险是 C.4 第 1 条的平滑偏差，**唯一能查它的是外部参照**。

## C.6 估计量：带内回归，以及它自带的诚实误差棒

在幅度带 $[a_{\min}, a_{\max}]$ 上采样，用

$$\kappa_a = \kappa_0 + c\,a^2 + O(a^4)$$

对 $a^2$ 做线性最小二乘。实测（鞍点，力形式，`experiments/sdr_ts.out`）：

| 模式 | 带 | $\hat\kappa_0$ | $\sigma_{\hat\kappa_0}$（拟合给出） | $\hat\kappa_0 - \lambda$（真实误差） | $\hat c$ |
| --- | --- | --- | --- | --- | --- |
| 不稳定 | $[10^{-3}, 2\times10^{-2}]$ | $-0.295963$ | $8.2\times10^{-7}$ | $+7.4\times10^{-7}$ | $28.93$ |
| 不稳定 | $[10^{-2}, 10^{-1}]$ | $-0.295431$ | $4.9\times10^{-4}$ | $+5.3\times10^{-4}$ | $28.19$ |
| 不稳定 | $[5\times10^{-2}, 3\times10^{-1}]$ | $-0.228623$ | $3.7\times10^{-2}$ | $+6.7\times10^{-2}$ | $22.82$ |
| 最低正 | $[10^{-3}, 2\times10^{-2}]$ | $+0.254958$ | $2.2\times10^{-8}$ | $+2.0\times10^{-8}$ | $3.43$ |
| 最低正 | $[10^{-2}, 10^{-1}]$ | $+0.254973$ | $1.3\times10^{-5}$ | $+1.4\times10^{-5}$ | $3.41$ |
| 最低正 | $[5\times10^{-2}, 3\times10^{-1}]$ | $+0.257002$ | $1.2\times10^{-3}$ | $+2.0\times10^{-3}$ | $3.25$ |

两条要记住的：

1. **$\sigma$ 诚实**。三个带上 $\sigma$ 与真实误差 $|\hat\kappa_0 - \lambda|$ 同量级
   （$8.2$e$-7$ vs $7.4$e$-7$；$4.9$e$-4$ vs $5.3$e$-4$；$3.7$e$-2$ vs $6.7$e$-2$）。
   拟合自己就知道它有多不准 —— 这是当前流水线**完全没有**的东西。
2. **带太宽会把 $\hat\kappa_0$ 系统性拉向零**（$-0.2286$ vs $-0.2960$，偏 23%）。
   对一个弱 TS 这足以翻分类。$\sigma$ 是防线：第三行的 $\sigma = 3.7\times10^{-2}$
   已经在报警。

$\hat c$ 跨带稳定（28.9 / 28.2 / 22.8），是可提取的真实非谐系数，现在被整个丢掉。

## C.7 方向与数字分开：这是路线 C 自洽的关键

Rayleigh 商给出

$$\kappa(v) = \lambda_1 + (\lambda_2 - \lambda_1)\sin^2\theta$$

**本征值误差是方向误差的二次量。**所以即使不信 Hessian 的**本征值**，它的**特征向量**
仍然可用。于是：

$$\boxed{\text{方向}\ \leftarrow\ \text{Hessian（或 dimer）}, \qquad
\text{数字}\ \leftarrow\ \text{带内回归}}$$

这不是妥协，是各用其可靠的那一半。dimer 精化时应当最小化**回归后的** $\hat\kappa$，
而不是 $\kappa_0$。

## C.8 幅度带怎么选：热幅度是物理锚点

对曲率 $\kappa$ 的模式，温度 $T$ 下的均方幅度是
$a_{\rm rms} = \sqrt{k_B T / \kappa}$（质量加权坐标）。$T = 300$ K（$k_BT = 0.02585$ eV）：

| 模式 | $\kappa$ | $a_{\rm rms}$ | 该幅度上 $\kappa_a$ 的偏差 |
| --- | --- | --- | --- |
| 鞍点不稳定模 | $0.296$ | $0.296$ | **已变号** |
| 最低正模 | $0.255$ | $0.318$ | $+115\%$ |
| 次低正模 | $0.654$ | $0.199$ | $+1.9\%$ |
| 极小点 mid | $5.10$ | $0.071$ | $+0.1\%$ |
| 极小点 stiffest | $55.6$ | $0.022$ | $+0.3\%$ |

**这是路线 C 最重的一个数**：软模的热幅度落在 $\kappa_a$ 偏差 100% 以上的区间，
不稳定模的热幅度处 $\kappa_a$ **连符号都不对**。对速率有意义的曲率是热可及幅度上的，
**$a \to 0$ 从来就不是 PRRS 的正确目标** —— 这个结论不依赖 MLP 可靠不可靠。

代价是必须承认：**在热幅度上"简谐曲率"这个概念本身对软模不成立**。所以路线 C 的正确
输出不是"一个更好的 $\kappa_0$"，而是 $\kappa_a$ 连同 $a$、$\sigma$、$\hat c$ 一起报。

## C.9 核心接口（覆盖 A.8）

```text
stationary        = true
band              = [a_min, a_max]      # 分辨率，必须显式
kappa_hat         = [...]               # 每个探测方向；NOT an eigenvalue (E.2)
sigma_kappa       = [...]               # 拟合标准误
c_anharmonic      = [...]               # 领头非谐系数
nu_band           = 1                   # index-like nested count，依赖嵌套顺序 (E.2)
sign_stable       = true                # 带内是否有变号
bottom_rotation   = theta(a)            # persistence，见 E.3
trivial_floor     = 2.9e-06             # D1，无参考误差地板
unstable_direction = v                  # persistent bottom direction，非 eigenvector
method            = "regression" | "hessian" | "lanczos"

# 以下两项仅当 method == "hessian" 时允许出现 —— 只有那时它才真是本征值
lambda_min        = -0.295964
imaginary_frequency_cm1 = -290.1
```

判据加**显著性检验**，这是现在完全没有的：

$$\boxed{\begin{aligned}
&\|\nabla E\| < F_{\rm tol} \\
&\exists!\, i:\ \hat\kappa_i < -\lambda_{\rm tol}\ \wedge\ |\hat\kappa_i| > 3\sigma_i \\
&\forall j \ne i:\ \hat\kappa_j > -\lambda_{\rm tol}\ \wedge\ \hat\kappa_j > 3\sigma_j \\
&\text{band 内无变号（否则 } \texttt{scale\_unstable}\text{，失败关闭）} \\
&|\lambda_{\rm trivial}|_{\max} \ll \lambda_{\rm tol} \quad \text{(D1 门)} \\
&R^\ddagger \pm \epsilon v_i \to A,\ B
\end{aligned}}$$

现状是 $\lambda = -1.2\times10^{-3}$ 对 tol $10^{-3}$ 直接过关，没有任何不确定度概念。

## C.10 落地清单（覆盖 A.13.6 的优先级）

> **第 1、2、3、6 条已落地**（2026-09-01，`105 passed, 1 skipped`，乙醇 P0 回归通过：
> 结构量全部逐位复现，虚频被修正约 1 $\mathrm{cm^{-1}}$）。第 6 条的球面最小化见
> [PRRS_STATUS.md](PRRS_STATUS.md) §8.8.3：一对 $\pm av$ 同时给目标与梯度，每步 2 次力评估。实现细节、真实后端实测与
> 一个顺带修掉的质量加权 bug 见 [PRRS_STATUS.md](PRRS_STATUS.md) §8。第 4、5、7、8 条未做。
> 值得单记一条：D1 门在真实乙醇上**把 FD 路径关闭了**（地板 $1.21\times10^{-3}$ 超过阈值
> $5\times10^{-4}$），解析路径通过（$2.5\times10^{-5}$）—— 一个此前静默的数值问题。


1. **D1 门，立刻做，最便宜**。`hessian_spectrum` 已经在算刚体基（`_trivial_modes`,
   `runner.py:38`），加 6 个 Rayleigh 商是零额外力评估。它是永久开启的无参考质量门。
2. **`hessian_spectrum` 改用 `calc.get_hessian()`**（A.13.6 第 1 条，实测 2 倍加速 +
   机器精度，注意 A.13.4b 的单位陷阱）。D1–D5 表明这个解析 Hessian 在 $10^{-5}$ 水平可信，
   所以这一步的前提成立。
3. **多尺度回归**：`response_curvature(atoms, v, band, K)`，**能量形式**（C.2），
   返回 $(\hat\kappa, \sigma, \hat c, \text{sign\_stable})$。$2K$ 次能量评估。
4. **复用幅度细化的数据**。搜索已经沿同一方向采多个幅度，回归很大程度上是**对已有数据的
   重读**，不是新开销。这是路线 C 最 PRRS-native 的一点，也是它比 Lanczos 更该先做的原因。
5. **带的默认值由 $a_{\rm rms}$ 定，不由数值考虑定**（C.8），且必须进 manifest。
6. **球面最小化取代 Lanczos 作为"找底方向"的手段**（E.4）：一对 $\pm av$ 同时给目标与梯度，
   每步 2 次力评估，与一次 $Hv$ 同价，不建矩阵不需要谱。$a \to 0$ 时退化为经典 dimer。
7. **Lanczos 继续下调**。A.13.5 已把它推到 $N \gtrsim 33$；C 路线说明就算实现了它给的也是
   $\kappa_0$ 而非带内 $\hat\kappa$；E.1 进一步说明有限 $a$ 下没有谱可供 Lanczos 求。
   它的位置收缩到 benchmark 档的 $\lambda_1$ 加速。
8. **唯一能查平滑偏差的事：外部参照**。乙醇三个 TS 上一个 DFT 单点 Hessian
   （比较 $\lambda_1$ 与 $v_1$）。MACE-OFF 是单模型，`Committee` 那条不可用（状态稿 §6）。
   **在做这一步之前，不要声称曲率是准的，只能声称它是自洽的。**

## C.11 仍未验证

1. **只测了乙醇的一个极小点和一个鞍点。**MACE-OFF 的训练集以近平衡有机构象为主；
   真正的成键/断键区（丙二醛质子转移，状态稿缺口 8）没测。**D1–D5 的结论不能外推到那里。**
2. **平滑偏差完全未测**，因为没有外部参照。这是本节最大的空缺。
3. **回归的带选择只在三个人工带上试过**，自动选带（按 $\sigma$ 与残差自适应）未实现。
4. **能量形式的数值下限未测**。$\operatorname{sinc}^2$ 的非负性是好的，但 $a$ 很小时
   $s(a)/a^2$ 的相消误差会主导；带下限该由这个决定，未测。
5. **$\operatorname{PT}$ 仍按投影实现**（B.3）。滤波的傅里叶论证假设沿 $v$ 的直线路径；
   在真正的 retraction 上传递函数会带路径曲率的修正，未推导也未测。

---

# 架构结论：Hessian 是 validation limit，不是 primary observable

三条路线的收敛点。**这一节的每条数学断言都已实测或解析验证**，脚本
[`experiments/finite_operator.py`](experiments/finite_operator.py)（MACE 鞍点）、
[`experiments/finite_operator_analytic.py`](experiments/finite_operator_analytic.py)（解析四次势）、
[`experiments/critical_points.py`](experiments/critical_points.py)（临界方向计数），
输出 [`experiments/finite_operator.out`](experiments/finite_operator.out)。

$$\boxed{\text{multi-scale finite response} \to \text{persistent bottom modes} \to \text{endpoint test}}$$

小体系上仍算完整 Hessian，但它的角色是 **benchmark**：
"finite response 说这是稳定的一阶负响应；full Hessian 也给一个虚频。"

## E.1 为什么 Hessian 必须是 validation limit：$a \to 0$ 是唯一有谱定理的地方

这个理由比"小体系便宜"强得多。**在有限 $a$ 下根本不存在算子，所以也不存在本征值。**

### E.1.1 代数

$g(x+d) = g_0 + Hd + \tfrac12 T[d,d] + \tfrac16 F[d,d,d] + O(d^4)$（$T = \nabla^3E$、
$F = \nabla^4E$，都全对称）。于是力形式的有限响应

$$D_a(v) := \frac{g(x+av) - g(x-av)}{2a} = Hv + \frac{a^2}{6}F[v,v,v] + O(a^4)$$

它**既非线性**（$D_a(\alpha v) = \alpha D_{\alpha a}(v) \ne \alpha D_a(v)$），
**也非对称**：

$$\boxed{\langle v_1, D_a v_2\rangle - \langle v_2, D_a v_1\rangle
= \frac{a^2}{6}\Big(F[v_1,v_2,v_2,v_2] - F[v_2,v_1,v_1,v_1]\Big) \ne 0}$$

（$v_1 \cdot Hv_2 = v_2 \cdot Hv_1$ 抵消掉，剩下的是 $F$ 的两个**不同**收缩。）

能量形式则是从标量导出的，没有对称性问题，但

$$\boxed{\kappa_a^{E}(v) = \langle v, Hv\rangle + \frac{a^2}{12}F[v,v,v,v] + O(a^4)}$$

是球面上的**四次型**，不是二次型。

**解析验证**（随机 $H, T, F$，$N = 8$，`finite_operator_analytic.py`）：三条式子在
$a = 10^{-3} \ldots 1$ 上全部精确 —— 非对称量测量/预测比值 **1.00000**，
$\kappa^E_a$ 公式误差在机器精度，第三条见 E.4。

### E.1.2 MACE 鞍点上实测（T1）

$v_1, v_2$ 取 Hessian 特征向量，所以 $v_1 \cdot Hv_2 = 0$ **精确成立** —— 测到的整个数值
就是 $a^2$ 项。这是最干净的演示：一个**本该为零**的非对角元，按哪个方向探测给出不同的值。

| $a$ | $\langle v_1, D_a v_2\rangle$ | $\langle v_2, D_a v_1\rangle$ | 非对称量 | $/a^2$ |
| --- | --- | --- | --- | --- |
| $10^{-2}$ | $1.982\times10^{-5}$ | $1.278\times10^{-4}$ | $-1.08\times10^{-4}$ | $-1.0796$ |
| $2\times10^{-2}$ | $7.93\times10^{-5}$ | $5.11\times10^{-4}$ | $-4.31\times10^{-4}$ | $-1.0786$ |
| $5\times10^{-2}$ | $4.95\times10^{-4}$ | $3.17\times10^{-3}$ | $-2.68\times10^{-3}$ | $-1.0716$ |
| $10^{-1}$ | $1.97\times10^{-3}$ | $1.24\times10^{-2}$ | $-1.05\times10^{-2}$ | $-1.0469$ |

**在最小的测试幅度上两者已差 6.4 倍**，$/a^2$ 列跨半个数量级稳定在 $-1.08$，正是 E.1.1
的预言。反推 $F[v_1,v_2,v_2,v_2] = 1.19$、$F[v_2,v_1,v_1,v_1] = 7.67$，
$(1.19-7.67)/6 = -1.0796$ —— 与实测**逐位吻合**。

### E.1.3 能量形式不是二次型（T2）

在 $(v_1, v_2)$ 平面上扫 $\theta$，用最佳二次型拟合 $\kappa_a^E(\cos\theta\, v_1 + \sin\theta\, v_2)$：

| $a$ | 残差 rms | 残差 max | $\kappa$ 量程 | 相对 |
| --- | --- | --- | --- | --- |
| $10^{-2}$ | $1.67\times10^{-5}$ | $3.20\times10^{-5}$ | $0.175$ | $0.01\%$ |
| $10^{-1}$ | $1.60\times10^{-3}$ | $3.07\times10^{-3}$ | $0.111$ | $1.4\%$ |
| $3\times10^{-1}$ | $1.05\times10^{-2}$ | $1.98\times10^{-2}$ | $0.451$ | $2.3\%$ |

残差按 $a^2$ 增长（$1.67\times10^{-5} \to 1.60\times10^{-3}$，$a$ 涨 10 倍、残差涨 96 倍）。

### E.1.4 谱结构分两阶段崩掉（临界方向计数）

二次型在球面上有**恰好 $N$ 个**临界方向（其特征向量，差符号），且互相正交。对
$\kappa_a^E$ 用 Newton 解 $\nabla_v\kappa_a^E = \mu v$、$\|v\| = 1$，从 3000 个随机起点数：

| $a$（合成，$N=8$） | 临界方向数 | 最低两个的 $\|v_1\cdot v_2\|$ | 结论 |
| --- | --- | --- | --- |
| $0$ | **8** | $0.000000$ | 谱（就是特征向量） |
| $0.1$ | 8 | $0.004551$ | 计数对，**正交性已丢** |
| $0.3$ | 8 | $0.041095$ | 同上（$\propto a^2$：$\times9$ 对 $a\times3$） |
| $1.0$ | 8 | $0.376408$ | 同上，已严重非正交 |
| $2.0$ | **28** | $0.081$ | **计数本身崩掉** |
| $3.0$ | **92** | $0.571$ | 同上 |

**先丢正交性（$O(a^2)$），再丢计数。**所以"$\mathcal L_a$ 的负本征值个数"在小 $a$ 下顶多是
"$N$ 个非正交底方向里有几个负"，在大 $a$ 下连这个都不成立。

**一个必须承认的损失**：A.6 里唯一严格的东西是 Cauchy 交错给的 $\nu$ **下界**。
交错定理要谱，**所以它不传到有限 $a$**。有限 $a$ 下连下界都没有了。

## E.2 命名规则要比你说的更严一档

你说"finite-response eigenvalue 不能再叫 frequency、不能报 $\mathrm{cm^{-1}}$"。对，但还要加一条：

$$\boxed{\text{它也不是 eigenvalue。}}$$

正确的名字与对象：

| 禁止 | 正确 | 理由 |
| --- | --- | --- |
| frequency, $\mathrm{cm^{-1}}$ | response curvature $\kappa_a(v)$ | 非简谐、尺度依赖（C.5 D2） |
| eigenvalue $\lambda(a)$ | $\kappa_a$ 在球面上的**临界值** | 无对称算子（E.1.2）、四次型（E.1.3） |
| eigenvector | persistent bottom **direction** | 非正交（E.1.4） |
| $\nu_a$ = 负本征值个数 | index-like **nested count**，依赖嵌套顺序 | 无交错、无正交基 |

所以核心接口里的 `lambda_min` 应改名 `kappa_min`，并且**必须带 `a`**。
只有 `method = "hessian"` 时才允许出现 `lambda`、`imaginary_frequency_cm1`，
因为只有那时它真的是本征值。

## E.3 persistent bottom modes 需要一个可测定义 —— 已测出来

"persistent" 不能只是修辞。可操作定义：

$$v_1(a) = \arg\min_{\|v\|=1}\kappa_a^E(v), \qquad
\text{persistent} \iff \begin{cases}
\theta(a) := \angle\big(v_1(a), v_1(a_{\min})\big) \text{ 连续且小} \\
\operatorname{sign}\kappa_a^E(v_1(a)) \text{ 带内不变}
\end{cases}$$

**实测**（T3，在 $(v_1, v_2)$ 平面内扫 $\theta$ 求 argmin；$v_1, v_2$ 为 Hessian 的最低两模）：

| $a$ | $\theta^*$ | $\kappa$ at $\theta^*$ | $\kappa$ at $0$ | 增益 | $\theta^*/a^2$ |
| --- | --- | --- | --- | --- | --- |
| $10^{-2}$ | $-0.013°$ | $-0.294516$ | $-0.294516$ | $3.0\times10^{-8}$ | 130 |
| $2\times10^{-2}$ | $-0.054°$ | $-0.290177$ | $-0.290176$ | $4.8\times10^{-7}$ | 135 |
| $5\times10^{-2}$ | $-0.379°$ | $-0.259944$ | $-0.259923$ | $2.1\times10^{-5}$ | 152 |
| $10^{-1}$ | $-2.592°$ | $-0.154190$ | $-0.153625$ | $5.7\times10^{-4}$ | 259 |
| $2\times10^{-1}$ | $\le -28.6°$（出窗） | $+0.130502$ | $+0.245464$ | $1.2\times10^{-1}$ | — |
| $3\times10^{-1}$ | $\le -28.6°$（出窗） | $+0.466436$ | $+0.827313$ | $3.6\times10^{-1}$ | — |

读法：**底方向在 $a \lesssim 0.05$ 上稳（转动 $< 0.4°$，按 $a^2$ 增长），到 $a = 0.1$ 转 2.6°，
$a \ge 0.2$ 已跑出 $\pm28.6°$ 扫描窗。**

对照 $(v_1, v_{\rm stiffest})$ 平面：所有 $a$ 上 $\theta^* = 0$ 且增益恒为 0 —— 转动只发生在
软模子空间内，这与直觉一致，也说明上表不是数值噪声。

**关键对照**：C.8 的热幅度是 $a_{\rm rms} \approx 0.30$。**在热相关幅度上底方向已不是
Hessian 特征向量。**无量纲组合 $a^2 c/\kappa$（$c = 28.93$、$\kappa = 0.296$，见 C.6）：

| $a$ | $0.01$ | $0.05$ | $0.1$ | $0.2$ | $0.3$ |
| --- | --- | --- | --- | --- | --- |
| $a^2c/\kappa$ | $0.010$ | $0.244$ | $0.977$ | $3.91$ | $8.80$ |

E.1.4 的合成表在这个无量纲量 $\approx1$ 处正交性已丢到 $0.38$、$\approx4\text{–}9$ 处计数崩掉。
**启发式对应（合成 $F$ 是随机稠密张量，分子 $F$ 有结构，所以这只是量级提示）**：
热幅度带处于非谱区域。这正是"带内 persistence + endpoint test"是正确认证、而谱计数不是的原因。

## E.4 一个白拿的便利：目标函数与梯度同价

$$\frac{\partial}{\partial v}\kappa_a^{E}(v)
= \frac{a\,g(x+av) - a\,g(x-av)}{a^2} = \frac{g(x+av) - g(x-av)}{a}
= \boxed{2\,D_a(v)}$$

**解析验证**：`finite_operator_analytic.py` 上 $\|\nabla_v\kappa_a^E - 2D_a(v)\| \approx 2\times10^{-10}$
（即数值梯度自身的差分噪声），$a = 10^{-3} \ldots 1$ 全部如此。

后果很实在：**一对 $\pm av$ 的评估同时给出目标 $\kappa_a^E(v)$ 和它的梯度 $2D_a(v)$**
（力本来就随能量一起返回）。所以球面约束下最小化 $\kappa_a^E$ 的每步只要
**2 次力评估 —— 与一次 Hessian-vector product 完全同价**，而且不建矩阵、不需要谱。

球面投影梯度 $-(I - vv^\top)\,2D_a(v)$；$a \to 0$ 时 $D_a \to H$，退化成经典 dimer /
Rayleigh 商最小化。**于是路线 A 的成本结构原样搬到路线 C 的算子上**，这就是
"不用 Hessian 怎么找 persistent bottom modes"的具体答案。

注意此处力形式**不是**作为滤波器使用（C.2 的非负性论证只管目标函数），它是目标函数的
梯度。两个形式的分工因此清楚了：**能量形式定义要最小化的量，力形式提供它的梯度。**

## E.5 修订后的架构

```text
小体系 / benchmark
    full analytic Hessian  ->  lambda, nu, cm^-1        <- 唯一有谱定理的地方
                                   |
                                   | 交叉校验：符号与 persistence 是否一致
                                   v
大体系 / 正式跑
    multi-scale finite response     kappa_a^E(v), 带内回归 (C.6)
              |
              v
    persistent bottom modes         球面最小化，2 次力评估/步 (E.4)
              |                     persistence 判据 (E.3)
              v
    endpoint test                   Phi(+v) = A,  Phi(-v) = B
```

三个属性（B.9）保持不变，只是第二项换了对象：

$$\text{stationarity} + \underbrace{\text{persistent index}}_{\text{不再是 Morse index}} + \text{connectivity}$$

## E.6 这对 reaction discovery 是优点，但要说清是哪个意义上的优点

**是优点**：目标是"找到反应通道并定端点"，而端点由 $\Phi$（淬火/盆地映射）判定，
**不由曲率数值判定**。曲率只需要给出一个可靠的**方向**和一个稳定的**符号**。
$\kappa_a$ 在热相关幅度上给的正是这个，而且比 $\lambda_0$ 更贴近反应实际发生的尺度
（C.8：软模 $\kappa_a$ 偏 115%，不稳定模符号都变）。

**不是"更准"**：$\kappa_a$ 不是对 $\lambda$ 的更好估计，它是**另一个量**。谁需要
$\lambda$（比如报 $\mathrm{cm^{-1}}$、算简谐配分函数、做 TST 前因子），必须走 $a \to 0$，
那就回到 full Hessian，也就回到 validation limit 的角色。

**代价清单**：失去谱定理、失去交错下界（E.1.4）、失去正交模基、失去"虚频"这个与量化
软件对话的通用货币。前三项对 discovery 无碍，第四项要靠 benchmark 档补。

## E.6b 谁来充当 benchmark 档

$a\to0$ 那一档目前是我们自己的解析 Hessian，也就是**自我一致**。外部对照的计划见
[PRRS_STATUS.md](PRRS_STATUS.md) §10（geomeTRIC：独立频率分析、真 IRC、鞍点精化、
内坐标约束优化），以及 §10.4 那个零新依赖的第一步（`ase.vibrations`）。

注意它校验的是**算法**不是 PES：同一个 calculator 跑外部优化器，平滑偏差依然不可见。
那一项仍然只有 DFT 能查。

## E.7 仍未验证

1. **$v_1(a)$ 的真实最小化未做**。T3 只在两个二维平面内扫 $\theta$，是下界式的证据；
   全空间球面最小化（E.4 的算法）未实现，所以 $\theta(a)$ 的真实大小只会**更大**。
2. **$\nu$ 的嵌套计数对顺序的依赖未测**。E.1.4 只证明了正交性丢失，没有直接测出
   两种嵌套顺序给出不同的 $\nu$。
3. **合成 $\to$ 分子的无量纲对应是启发式**（E.3 末段）。
4. **`persistence` 的容差没有定**（$\theta$ 允许转多少度、需要多少个 $a$ 点）。这必须由
   丙二醛那类真实化学基准定，乙醇上定不了。
5. **E.4 的球面最小化在负曲率方向上的收敛性未测**。目标函数在 $\kappa < 0$ 区域的球面
   几何与经典 dimer 不同（四次项可能引入额外临界点，E.1.4 的 $a=2$ 行）。

---

# 三条路线与状态稿缺口的对应

| 缺口 | 路线 A | 路线 B | 路线 C |
| --- | --- | --- | --- |
| 1. 自由坐标未商掉 | — | $\mathcal N$（Morse–Bott），B.6 给出修法 | 带内 $\hat\kappa$ 与 $\sigma$ 给出"是否真自由"的显著性判据 |
| 2. `free`/`stiff` 两级只有解析测试 | — | $\mathcal N$ 的判据从容差变成定义 | `free` 变成 $|\hat\kappa| < 3\sigma$，可测 |
| 3. 排序仍是启发式 | $\kappa(v)$ 的 2 次力评估代理 | $\kappa_a$ 谱可给 $U_{\rm TS}$ | $\hat c$（非谐系数）是新的排序信号，实测跨带稳定 |
| 4. 大体系没有 $n_-$ 判定路径 | A.6 给判据，但 A.13.5 把 Lanczos 推到 $N \gtrsim 33$ | 命名与接口（$\nu$，非虚频） | $\nu_{\rm band}$ + `sign_stable` 取代点判定 |
| 5. mode-resolved 判据只覆盖 torsion | projected HVP 给任意软模 $k_v$ | $\mathcal L_a$ 的一般定义 | 回归给任意方向的 $(\hat\kappa, \sigma)$，不限 torsion |
| 9. FD 步长与容差矛盾 | — | — | **解析 Hessian 直接消除**（C.10 第 2 条）；D1 门永久监控 |

---

# 边界与未验证部分

**已实测的部分**（2026-08-31 至 09-01，脚本与原始输出都在 `docs/experiments/`）：

| 节 | 内容 | 性质 |
| --- | --- | --- |
| A.12 | Lanczos 判据的误判行为 | **合成矩阵**，$N = 300$ 已知谱，不是 MACE 谱 |
| A.13.4 | 解析双反向 $Hv$ 的可行性/精度/成本/显存 | 乙醇 9 原子，**一个非驻点** |
| C.3 | $\mathrm{dev}_F/\mathrm{dev}_E = 2.00$，滤波传递函数的定量确认 | 乙醇鞍点，3 个模式 × 6 个幅度 |
| C.5 | D1–D5 二阶导可靠性诊断 | 乙醇**一个极小点 + 一个鞍点** |
| C.6 | 带内回归的 $\hat\kappa$、$\sigma$ 与真实误差的一致性 | 同上，3 个人工带 |
| C.8 | 热幅度 $a_{\rm rms}$ 与 $\kappa_a$ 偏差的对照 | 由 C.5 的 $\lambda$ 与 D2 表算出 |

路线 A 的 Lanczos 未实现；路线 B 整体是理论组织，$\operatorname{PT}$ 未按输运实现、
$\mathcal N$ 对齐未实现；路线 C 的回归尚未进入代码。

具体未定的点：

1. **只有乙醇，且只有两个几何。**MACE-OFF 训练集以近平衡有机构象为主，成键/断键区
   （状态稿缺口 8 的丙二醛质子转移）没测。**D1–D5 的乐观结论不能外推到那里。**
2. **平滑偏差完全未测。**所有无参考诊断（D1–D5）都查不出它。唯一手段是外部参照：
   乙醇三个 TS 上一个 DFT 单点 Hessian。MACE-OFF 是单模型，`Committee` 不可用（状态稿 §6）。
   **在做这一步之前只能声称曲率自洽，不能声称它准。**
3. **$m$ 的实际值在真实体系上未测**（A.12 的 $40\text{–}60$ 来自合成谱）；块方法未实现。
4. **$\nu = 1$ 在无矩阵路线上不可认证**（A.6，已由合成谱确认）；余量因子
   $\gamma \approx 10$ 是启发式，没有证明。
5. **自动选带未实现**；能量形式在小 $a$ 端的相消误差下限未测（C.11 第 4 条）。
6. **大体系显存标度未测**（$N = 9$ 上 118 MiB，`get_hessian` 输出是 $O((3N)^2)$，
   `chunk_size` 硬编码）。
7. **forward-over-reverse**（`torch.func.jvp`）与 **openmm-ml 路径下的解析二阶导**未测。
8. **HVP / 回归的失败关闭门控**未设计（A.13.2 第 2 条）。
9. **滤波论证在真 retraction 上的修正**未推导（C.11 第 5 条）：C.1 的傅里叶论证假设沿
   $v$ 的直线路径。
10. **`get_hessian` 不做单位换算**（A.13.4b）。MACE-OFF 两个因子都是 1.0 所以现在静默正确。
