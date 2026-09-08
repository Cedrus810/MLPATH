# b04 — Cl⁻ + CH₃Cl

**这一格问什么**（按 `../README.md` 规矩 5，标签在输出侧）：

> 取代通道上，反应物与产物之间是否存在一个 `confirm_minimum` 通过的极小？
> 搜索必须给出答复。**两个答案都是 pass，答不出来才是 fail。**

教科书会说这是 SN2、是协同的。**那句话不进判据。** 它进 `reference.json`，
作为一个可被证伪的预期。

## 为什么它不能被 b01–b03 代替

| | b01–b03 | b04 |
| --- | --- | --- |
| 总电荷 | 0 | **−1** |
| 起始碎片数 | 1 | **2** |
| 反应坐标 | O–H 伸缩 | **C–Cl，重原子** |
| 模型 | MACE-OFF24 | **charge-aware** |

b01–b03 是同一类反应的三个装饰版本（见 `../TAXONOMY.md` A 组）。
b04 是这条轴上第一个换了化学的取样点，也是**离子块 b04–b06 的启用格**：
它要打通的三处代码缺口，后面两格全部搭车。

## 已知缺口（跑前必须解决，`preflight.py` 会逐条报）

实测于 2026-09-04，实现树 `23179deadf657452`。

### 缺口 1 —— 域门拒绝裸阴离子（**阻塞**）

```
Cl⁻ 单独        → stereo_unresolved = "no_assignment_satisfies_valences"
Cl⁻ + CH₃Cl     → 同上
```

`chemistry.VALENCE` 是**中性**闭壳价数表，Cl 记为 1 价。裸 Cl⁻ 度数为 0，
剩余价数放不下去，于是 `admissible_bond_orders` 返回 None，
而 `network.Registry.admit` 在 `closed_shell_only=True` 时把这个理由判为
`out_of_domain` 并拒绝。

**Cl⁻ 是完美的闭壳物种，只是不中性。** 这不是域外，是价键模型写死了中性。

`docs/experiments/sn2_probe.py` 当时设 `closed_shell_only=False` 绕过。
**基准不能这么做** —— 那样域门整个失效，就没有任何东西在守训练域了，
而守训练域正是这一格换模型之后最需要的东西。

修法：让价键模型接受总电荷，`total_charge = −1` 时 Cl⁻ 价数为 0。
这会改 `chemical_key`，因而**会改 digest** —— 属于有意改变，走 `--allow`。

### 缺口 2 —— `mace_factory` 装不了 charge-aware 模型（**阻塞**）

```python
MACECalculator(model_paths=..., device=..., default_dtype=...)   # 现状
```

没有 `model_type`。POLAR-1 需要 `model_type="PolarMACE"`（其 checkpoint 类是
`PolarMACE`），omol-0 需要 `"MACE"`。装错类型不会报错，会安静地算错。

### 缺口 3 —— 电荷传不到势能上（**阻塞**）

`SearchConfig(charge_sensitive=True, total_charge=-1, multiplicity=1)` **可以构造**
（已实测），且 `network.py:149` 把它送进 `chemical_key` ——
所以**身份层知道电荷**。

但没有任何东西把电荷送给 calculator。charge-aware MACE 从
`atoms.info["charge"]` / `["spin"]` 读，而 `src/prrs/` 里没有一处写这两个键
（`grep 'info\['` 只命中 `info["quench"]`）。探针是手工设的。

后果最坏的形式：**配置声明 −1、模型按 0 算，两边都不报错。**
`build.py` 因此把 charge/spin 写进结构的 `info`，但那只覆盖源结构 ——
搜索过程中每一个 `atoms.copy()` 出来的新结构是否带着它，未验证。

## 状态

`blocked`。前置是 **b03 收口**（用户 2026-09-04 定），
以及上面三个缺口。判据在缺口解决前不冻结 —— 现在写出来的每一条
都会被域门拒掉，报出来的只会是一列门拒绝。
