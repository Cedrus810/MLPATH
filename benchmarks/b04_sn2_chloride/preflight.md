# b04 preflight 结果

2026-09-04 · 环境 `openmm_dev` · RTX 2080 Ti · 实现树 `23179deadf657452`
脚本 `preflight.py`（gate 1–3 无需 GPU；gate 4 需 `--models`）

**这份文档没有一条是验收结果。** 它回答的只是"这一格能不能跑"。

## 5/5 通过（2026-09-04 复核后）

原报 3/5、两条阻塞。两条的性质在换 POLAR-1 与加 `_stamped_with_charge` 之后都变了，见下。

| gate | 结果 | 一句话 |
| --- | --- | --- |
| 1 裸氯离子是否被身份层否决 | **PASS** | 那道门的理由是 OFF24 的中性训练域，对 POLAR-1 不成立；冻结配置声明 `closed_shell_only=False`。身份层怎么说仍然记录，只是不再是否决 |
| 2 SearchConfig 接受非零总电荷 | PASS | 身份层本来就是电荷感知的 |
| 3 总电荷能否到达势能 | **PASS** | `prrs.search._stamped_with_charge`（`cf29f66b`）按 config 盖章，两者冲突则拒绝开跑 |
| 4 模型是否真的响应电荷 | PASS | 两个候选都强响应 |
| 5 结构里的元素模型是否都声明了 | PASS | 从 checkpoint 现读，不靠手抄清单 |

## gate 4 是这轮唯一的新测量

同一个 b04 输入几何，只改 `atoms.info["charge"]`：

| 模型 | E(0) | E(−1) | E(+1) | 跨度 |
| --- | --- | --- | --- | --- |
| MACE-omol-0-4M | −26128.656934 | −26132.720032 | −26124.438921 | **8.281 eV** |
| MACE-POLAR-1-M | −26128.290206 | −26132.741299 | −26116.074480 | **16.67 eV** |

对照 MACE-OFF24：同样三个电荷**逐位相同**（`PRRS_STATUS.md` §2.5）。

**这条只证明模型不瞎，不证明它算得对。** 跨度大小本身不是精度指标 ——
两个模型在 E(−1) 上差 0.021 eV，在 E(+1) 上差 8.4 eV，而没有参考计算，
所以哪个对不知道。**不要引用跨度去比较两个模型。**

## gate 1 与 gate 3 的形状不同，修法也不同

**gate 1 是设计缺口。** `chemistry.VALENCE` 是中性闭壳价数表，Cl 记 1 价；
裸 Cl⁻ 度数 0，剩余价数无处安放，`admissible_bond_orders` 返回 None，
`Registry.admit` 判 `out_of_domain` 并拒绝。Cl⁻ 是完美的闭壳物种，只是带电。
修它要动 `chemical_key`，**会改 digest，属于有意改变**，走 `--allow`。

`sn2_probe.py` 当时设 `closed_shell_only=False` 绕过。基准不能这么做：
那样域门整个失效，而换模型之后守训练域正是最需要的东西。

**gate 3 是管道缺口，而且是静默的那种。** 配置声明 −1、模型按 0 算，
两边都不报错 —— 这是最坏的失效形式，因为它产出的是一整套自洽的错数字，
和 P3 那个 E 异构体同一个形状（`docs/HANDOFF.md` §5）。
`build.py` 给源结构盖了 charge/spin，但搜索过程中每个 `atoms.copy()`
是否带着它**未验证**，而 gate 3 现在连第一步都没过。

## 还没测的

- 分子间 pair 提案在 3.6 Å 上真的会被提出来吗（`pair_cutoff_A = 4.0` 是上界，
  没有实测过两个碎片之间的提案路径）
- 两个自由碎片的相对平动/转动落进 conformer 储层是什么行为
  （`conformer_max_per_node = 8` 是在苯甲酰丙酮的单构象螯合盆地上定的）
- `collateral_gate` 对跨碎片探针的行为（`fragments_after != fragments_before` 会拒，
  但两碎片合成一个的方向没测过）

这三条要等 gate 1/3 修完才有意义 —— 现在测会全部撞在域门上。

## gate 5 顺带把一个零覆盖的检查补上了

`src/prrs/` 里**没有任何东西**检查结构的元素在不在模型的训练元素集里。
在一个模型永远没见过的元素上跑，会产出数字而不产出警告。

gate 5 从 checkpoint 现读（`scripts/model_registry.py` 用 `model.atomic_numbers`），
不依赖任何手抄清单：

| 模型 | 元素数 | `chemistry.VALENCE` 覆盖其中 |
| --- | --- | --- |
| MACE-OFF24_medium | 10 | 7 |
| MACE-POLAR-1-M | 83 | 7 |
| MACE-omol-0-4M | 82 | 7 |

**这两列回答的是不同的问题，不要合并。** 左列是"模型能不能算"（域），
右列是"能不能给出 Lewis 结构"（价键）。模型对后者无话可说 ——
Pd 没有单值价数，那是化学不是元数据。第二列不动也不是错：
`unsupported_elements` **不被拒绝**，key 退回图 + 碎片，
所以 76 个元素是"能搜、但没有立体化学"，不是"进不去"。

真正把 b04 挡住的是 gate 1，而它和元素覆盖无关 ——
Cl 就在那 7 个里面，问题是价数表写死了中性。
