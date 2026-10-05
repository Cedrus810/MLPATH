# b04 的 T[E_m]：ORCA 的优化器在 POLAR-1 面上自己找鞍点

P1 那条线（`../../p1_qm/orca/p1_T_Em/`）建在 MACE-OFF24 上。P4+ 的默认模型是
POLAR-1，**换了模型，`A_θ[E_m] vs T[E_m]` 这一行要重新建** —— `../` 里的
b04 QM 复核做的是 `T[E_QM]` 与 PRRS TS 的 QM 精化，不是这一行。

- 算子：ORCA `ExtOpt OptTS`（RS-P-RFO，内坐标），经 `scripts/mace_engrad.sh`
- 势能面：MACE-POLAR-1-M，总电荷 −1（由 `* xyzfile -1 1` 传入，
  `mace_engrad.py` 写进 `atoms.info`，PolarMACE 从那里读）
- 起点 `g.xyz`：**独立初猜**，D3h、两个 C–Cl 各 2.30 Å ——
  既不是 POLAR-1 的 2.2805/2.2827，也不是 QM 的 2.3645。协议要求起点不取自任何一次结果
- 不加 `Freq` / `NumFreq` / `Calc_Hess`：ExtOpt 下频率任务必失败
  （外部程序没有波函数）。曲率用 PRRS 自己的解析 Hessian，形状见
  `../../p1_qm/compare_saddles.py`

判读**跑前写死**：`T[E_m]` 与 `A_θ[E_m]`（`runs/b04_frozen2/ts_candidates/ts0000`）
的对称感知 RMSD 应在 0.01 Å 量级以内，且两者各自到 `T[E_QM]`（2.3645 Å 的 D3h 鞍点）
明显更远。不出这个形状是新事实，不得沿用 P1 的结论。

## 第一次尝试：没收敛（2026-09-22，留档不覆盖）

`t_optts.inp` / `t_optts.out` —— `! ExtOpt OptTS`，无初始 Hessian。ORCA 正常退出，
但跑满 50 个 cycle 没有 HURRAY，几何漂到 C–Cl 1.765 / 2.548 Å，不是对称鞍点。
RS-P-RFO 没有 Hessian 就得自己猜跟哪个模式，而 P1 那次 13 步收敛是因为起点
本来就贴着鞍点。**这不是 POLAR-1 面上没有鞍点，是优化器没被给够信息。**

引擎侧的好消息在同一份 out 里：`External energy and gradient 66.719 sec (92.3%)`，
50 cycle 约 1.3 s/次，对 P1 那次的 293.8 s (100.0%)、约 5 s/次 ——
`scripts/mace_server.py` 第一次上生产，单次调用降到四分之一。剩下的 1.3 s 是
ORCA 每个 cycle 重新 exec 一次 python，ExtOpt 的模型决定的。

## 第二次：给初始 Hessian

`t_optts_hess.inp` —— 加 `%geom Calc_Hess true MaxIter 200 end`。
`Calc_Hess` 走梯度差分（6N = 36 次调用），ExtOpt 供得起；
P1 记录里挂掉的 `Freq` / `NumFreq` 是要波函数的属性积分，不是这个。
第一次的输出**不覆盖**，两次各自留档。

### `Calc_Hess` 也挂（2026-09-22，比 P1 记录的范围更宽）

`calchess_FAILED.inp` / `.out`：`%geom Calc_Hess true` 在第 1 个 cycle 就报
`ERROR (SHARK): Failed to read input file` —— 和 `Freq` / `NumFreq` 同一个错。
所以 P1 那条记录该改写成：**ExtOpt 下任何由 ORCA 自己产生 Hessian 的任务都失败**，
不止频率任务。它们都走属性积分，外部程序没有波函数。

### 第三次：外部 Hessian（`InHess Read`）

`scripts/prrs_to_orca_hess.py` 把 PRRS 的解析 Cartesian Hessian 写成 ORCA `.hess`
（`guess.hess`）。在 `g.xyz` 上它是：

```
最低六个本征值 (eV/A^2): [-5.661  -0.0128  -0.0112  -0.0109  -0.  0.]
负模 4 个 —— 一个真的（-5.661，SN2 反应坐标），三个是未投影的刚体近零
provenance: analytic，HVP 对有限差分的偏差 3.7e-4，非对称度 1.8e-9
```

**这印证了第一次失败的诊断**：起点确实带着正确的负模，是优化器没有 Hessian、
按模型 Hessian 猜了个模式去反转，才走开的。

**这不削弱 T[E_m] 的独立性**：独立的是**优化器**（RS-P-RFO / 内坐标）与**起点**
（不取自任何一次结果）。Hessian 是势能面的性质，两个算子本来就都在问同一张面要梯度。
