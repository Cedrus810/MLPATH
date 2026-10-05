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

## 已知缺口 —— **三条都已解决**（2026-09-04 / 09-09）

原文保留，处置写在每条末尾。实测于 2026-09-04，实现树 `23179deadf657452`；
当前树 `8418b82b1f6197d6`。

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

~~修法：让价键模型接受总电荷~~ —— **处置换了，而且更便宜（2026-09-04）**：
那道门的理由是 **MACE-OFF24 的中性训练域**，对 POLAR-1 不成立，
所以拒绝本来就不该由价键表发出。`5d905aaa` 把三层拆开（身份内容 / 身份解析状态 /
模型适用域证据），冻结配置声明 `closed_shell_only=False`，
**守训练域的职责移交 `model.json` 与 preflight**（checkpoint 元素覆盖 gate 5、
实测电荷响应 gate 4）—— 那才是域证据该待的地方。
身份层怎么说仍然逐条记录（`stereo_unresolved`），只是不再是否决。
**`chemical_key` 一字未改，digest 零差异，不需要 `--allow`。**
彻底重做 VALENCE 已另立为身份协议研究，不阻塞主线（`docs/DECISIONS_PENDING.md`）。

### 缺口 2 —— `mace_factory` 装不了 charge-aware 模型（**阻塞**）

```python
MACECalculator(model_paths=..., device=..., default_dtype=...)   # 现状
```

没有 `model_type`。POLAR-1 需要 `model_type="PolarMACE"`（其 checkpoint 类是
`PolarMACE`），omol-0 需要 `"MACE"`。装错类型不会报错，会安静地算错。

**已解决（2026-09-04）**：工厂放在本格的 `calculator.py`，经既有的
`--calculator module:function` 钩子进来，`mace_factory` **不动** ——
在有东西真正走过这条路之前不去拓宽它（理由写在 `RUN.md`）。三次探针都走的这条路。

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

**已解决（`cf29f66b`）**：`prrs.search._stamped_with_charge` 按 config 给结构盖章，
两者冲突则**拒绝开跑**。preflight gate 3 转绿就是这一条，三次探针都在它下面跑的。

## 状态

`preflight`（`MANIFEST.md` 同步）。**三个缺口都已解决，preflight 5/5，三次探针已跑。**
b03 收口也完成（`../b03_benzoylacetone/README.md`：判据已通过、几何门已建、
`side +1` 记为该势能面的已知限制）。

**判据文本已写定**（2026-09-09，`CRITERIA.md`，依用户裁决逐条给定），
`check.py` 已实现、含 3 例自检。**冻结＝逐字不改，生效由用户确认。**

在 `runs/b04_probe4_ledger` 上机械评定的结果（2026-09-09，181 上跑，带账本）：

```
取代通道   observed_direct      ts0000 order 1 · rc0000 自环 · 两端 m0000/m0001
                                D1 D2 D3 D4 全 PASS
分步侧     unanswered           理由 candidate_not_selected + proposed_not_executed
                                NOT searched_no_hit —— 基被截断
逐微态     m0000/m0001 各 uncovered 14 · not_executed 5/4 · censored_by
                                direction_quota_fraction  → 负结论不可采纳
exit 0                          记录可裁定
```

**这是 b04 第一次给出完整的、按裁决语义的答复。**
probe4 与 probe3 逐字段相同（唯一差别是虚频末几位，GPU 浮点非确定性），
所以它同时证明方向合同那批改动**没有动搜索行为**。
`runs/b04_probe3_rigidbond` 上跑同一个检查会得到 `NOT ADJUDICABLE` / exit 1 ——
那不是化学结论，是那次运行早于 `direction_ledger`。

摆在那里等的材料：

| | |
| --- | --- |
| preflight | 5/5，见 `preflight.md`；三条"还没测的"里两条已答，第三条据实记为仍未触发 |
| 探针 | `runs/b04_probe`（首次）· `runs/b04_probe2_openceiling`（上限不 binding）· `runs/b04_probe3_rigidbond`（交付修复后：0 拒绝、1 个 bracket `[0.75, 0.80]`、ts0000 独立复现） |
| 候选与 QM | `ts0000` order 1；`docs/experiments/b04_qm/` 精化到同一 QM 鞍点、双向 IRC 收敛、`event_key` 逐字一致 |
| 通道 | `rc0000` `degenerate_reaction`、自环、broken `[0,1]` / formed `[0,5]` |
| 要一起冻的 | `config.json`（probe3 逐字段用的就是它，只有探针期的 `saddle_search_max_rise_eV=10.0` 是覆盖值）与 `model.json` |
| 不进判据的 | 势垒、虚频数值（规矩 2）；机理标签只能是输出（规矩 5） |

**能量区域上限尚未冻结**，按 4a 它是 B 组跑完之后的输出：现有数据点 n = 2，
`max_rise_eV` 1.5479（重叠 0.699 的种子）与 1.6821（重叠 0.986 的种子）——
**种子质量与能量高低不同向**，所以从这两个数反推一个上限仍然无效。
