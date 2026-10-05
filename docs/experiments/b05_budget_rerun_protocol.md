# b05 两侧下降的预算重跑协议(PLAN 第 5 件 · 5.0,零代码)

**状态:协议已写,待在有 MACE checkpoint 的机器(kasuga02)执行。**
本机 yayoi 无 `MACE-OFF24_medium.model`,跑不了 POLAR-1。

## 为什么是这一步(出处,不改写)

`DECISIONS_PENDING.md` ② 的裁决原话:

> **quench 预算曾阻止三个样本得到完整结局;补算得到一致的消除产物型 order-1 驻点。**
> 瓶颈的正确说法:**对这三个已执行样本,首先暴露的是淬火收敛预算,随后是产物驻点的极小确认;
> 未选 155 个方向的覆盖限制仍然存在。**
> **③(给碎片软模加打磨机制)现在没有依据**:六条下降是被两道从未放开过的预算挡住的,
> 不是被证明处理不了。先放开预算重跑,再谈是否要新机制。

B(碎片刚体转动,2026-09-13)**不是**它的替代品:B 之后鞍点不动,残余 −54 cm⁻¹ 是耦合的
碎片间运动。B′(打磨扫负曲率)按 2026-09-20 决策 6 挂起,**触发条件就是本协议跑完之后**
残余负曲率仍在精确周期坐标上。

## 被挡住的实测(2026-09-09,`runs/b05_saddle_descent/`)

六条下降 = 三个 order-1 驻点 × 两侧,全部停在同一组成(`C4H8 + Cl + H2O`),无一到确认极小:

| 停止原因 | 条数 | 挡住它的预算 |
| --- | --- | --- |
| `fmax_not_reached` | 3 | `quench_steps = 500`(续算实测此处差约 59 步) |
| `order 1 / relays` 用尽 | 3 | `irc_relay_hops = 3`(淬火已收敛 fmax 0.0085–0.0096,落点仍 order 1) |

## 重跑参数(改动仅此两项,其余逐字节不动)

| 参数 | 原值 | 新值 | 依据 |
| --- | --- | --- | --- |
| `quench_steps` | 500 | **1500** | 续算实测差约 59 步(其 2×=1000 已收敛);下降起点更远,取 3× 留余量,不取无上限 |
| `irc_relay_hops` | 3 | **8** | 观测到中继用尽;取观测值之上仍有限的值,链跑飞仍会 fail closed 停下 |

## 执行步骤(kasuga02)

```
1. 先核实现:make hash 必须等于当时发布的指纹(跑之前记录在案的 implementation_sha256;
   跨机先对 src md5,见 env/BASELINE.md)
2. cp docs/experiments/b05_saddle_descent.py docs/experiments/b05_saddle_descent_budget.py
   只改两处:load_config() 里 replace 两项预算(1500 / 8),输出目录
   runs/b05_saddle_descent_budget(与原目录分开,三份结果各自保留,不覆盖)
3. python docs/experiments/b05_saddle_descent_budget.py
4. 结局逐条记:每侧 confirmed_minimum(附组成与能量)/ order-k 鞍点(k 与虚频)/
   预算停(附具体预算名——新预算下再停就是新事实,不是复现)
```

## 判读(先写下来,免得看完结果再定口径)

- 六条全部到确认极小 ⇒ 消除通道的极小确认补全;与 B 之后的 order-1 / −54 cm⁻¹ 对照,
  鞍点若仍不动,则「产物区停不下来」的剩余原因在耦合碎片间运动(B″/第 5 件的靶子),
  B′ 依然不启用。
- 仍有预算停 ⇒ 新预算下的停止是**新事实**,逐条记录,不得沿用旧结论。
- 判据层面:两条 b05 通道的正式判定**仍由 b05 的 check.py 说了算**,本协议只产出证据。
- 本协议是一次**独立探针**(预算放开),不回写 probe1–4 的任何历史结果。

---

## 执行结果 · 2026-09-22 · `runs/b05_saddle_descent_budget/`

kasuga02，`implementation_sha256 = ec32bc3df5d2bf14`（跑前两机核过）。
脚本 `docs/experiments/b05_saddle_descent_budget.py` —— 由本目录的 `b05_saddle_descent.py`
复制而来，**只改两处**（两项预算 + 输出目录），两处锚点都带 `assert`。原目录未覆盖。

```
                                              旧预算        新预算
到确认极小的下降条数                            0 / 6        4 / 6

t000002  lowest eigenvalue −0.02166   +1 confirmed order 0 relays 1
                                      −1 confirmed order 0 relays 4
t000058  lowest eigenvalue −0.02197   +1 confirmed order 0 relays 3
                                      −1 fmax_not_reached_after_polish
t000080  lowest eigenvalue −0.02088   +1 fmax_not_reached_after_polish
                                      −1 confirmed order 0 relays 1

四条到达的落点全同：C4H8 + Cl + H2O · key a4cc3716d8bb73fe · dE ≈ −1.322 eV
```

**按本协议的判读（跑前写下的那两条）逐条对：**

- 「六条全部到确认极小」**没有发生** —— 4/6。
- 「仍有预算停 ⇒ 新预算下的停止是新事实」**发生了**，而且停止理由是新的：
  `fmax_not_reached_after_polish`，**不是**旧的 `quench_steps` 或 `irc_relay_hops`。
  按协议逐条记录，**不得沿用旧结论**。

2026-09-09 裁决「六条下降是被两道从未放开过的预算挡住的，不是被证明处理不了；
先放开预算重跑，再谈是否要新机制」——**押对了**：放开预算之后三分之二的下降完成了。

### 残余负曲率在哪：转动，不是平动

```
lowest eigenvalue −0.0217 eV/(A²·amu)
  99.4% 落在水碎片上
  水内部：转动 99.9% · 平动 0.1% · 内部形变 0.0%
```

同日测出的平动阈值（`runs/b05_fragment_translation_measurement/`，
`d²E/dΔ² = 1.737 eV/Å²`、`|dE/dΔ| = 0.0930 eV/Å`，三个端点一致）**数是好的，但它量的自由度
只占残余模的 0.1%**。「碎片平动可能是那个 −54 cm⁻¹ 的出路」这个假设到此被实测否掉。

### 决策 6 的触发条件成立

本协议文首写着：B′（打磨扫过负曲率）按 2026-09-20 决策 6 挂起，
**触发条件就是本协议跑完之后残余负曲率仍在精确周期坐标上**。

跑完了，残余负曲率在**转动**上，而转动正是精确周期坐标（`internal.rotate_rigid`，
第三种精确刚体实现）。**条件成立，B′ 等用户落笔。**
