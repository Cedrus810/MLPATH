# PRRS — 扰动–响应反应搜索（Perturbation–Response Reaction Search）

一个 MLP-native 的反应发现方法：**给分子一个受控扰动 → 让它在无偏 NVE 下自由响应 → 淬火 →
看它落到哪个化学态。**

不预设反应坐标，不预设产物，不从过渡态猜起。先发现连通性，之后才做 TS 与势垒的表征。

- 方法原稿：`mlp_native_perturbation_response_reaction_search.md`
- 实现状态（逐条实测）：`docs/PRRS_STATUS.md`
- 交接文档（什么已冻结、什么还开着）：`docs/HANDOFF.md`
- 英文版本：`README.md`

**这个工作目录不是 git 仓库。** 唯一的代码指纹是对 `src/prrs/*.py` 的
`implementation_sha256` —— 用 `make hash` 现算，**不要手抄**。

---

## 1. 机理

### 1.1 主循环

```
状态 R0  →  扰动 P_k  →  短 MLP 轨迹 τ_k(t)  →  淬火  →  分类
                                                         │
       构象 microstate ◄── chemical_key 相同 ── 盆地 ──► 新化学节点
                                                         │
                                      一阶鞍点（n_- = 1）──► TS pool
```

每次试验都是一个小数值实验，分两段：**扰动交付段**，然后是撤去扰动的**自由松弛段**。回到源
极小是弹性证据；松弛到另一个极小是盆地发现；停在恰好一个负本征值的驻点上是 TS 候选
（`n_- > 1` 记为可诊断的拒绝）。

### 1.2 扰动家族

六个家族在 primitive internal coordinate（键长/键角/二面角）上提案方向，而不是笛卡尔噪声。
两条实测结论定死了这一层：

- **扰动基必须张开主要内坐标子空间。** 乙醇首轮 24 次试验 6/6 全弹性，而 anti / gauche⁺ /
  gauche⁻ 三个真实盆地就在 4 meV 处 —— 原三家族（键拉伸、非键接近、pair kick）全是成对径向
  扰动，100 fs 内耦合不到扭转。集体/扭转扰动因此从"后续增加"提为 MVP 必需项。
- **target-coordinate fidelity ≠ geometric fidelity。** 2.09 rad 的扭转位移曾把请求的二面角
  精确达标，同时把 O–H 从 0.957 Å 拉到 2.346 Å。现在 `apply()` 分别返回 `target_error` 与
  `collateral`：共价键损伤与原子重叠是硬门（0.05 Å），非目标键角与扭转只测量（同顶点角度有
  求和约束，是合法耦合）。

交付是增量的（`MAX_STEP = 0.05` Å 或 rad，每步重算梯度），因为**大扰动 = 一串局部有效的小
扰动**。桥键扭转走精确刚体 fragment 旋转（所有键长键角保持到 ~1e-16）；环上的键不是桥，走
continuation。

### 1.3 在对称商空间上搜索

商掉冗余自由度是全套设计里唯一反复出现的抽象，已在六处独立救场：

1. `chemical_key` 的图部分用自同构不变哈希（否则转移三个等价甲基 H 中的哪一个会产生三个"产物"）
2. 构象去重用对称感知 RMSD（甲基转 120°：plain RMSD 0.914 Å vs symmetric 0.006 Å）
3. 角度方向枚举按正则色去重
4. 对称转子降级 + 基本域限制
5. reservoir 多样性度量排除转子扭转
6. 立体化学签名只取**构型**手性，绝不含可旋转单键的扭转符号（否则 gauche⁺/gauche⁻ 会变成两个
   化学节点）

对称转子**降级但保留**：它对盆地发现无用，却实测贡献了一个真实一阶鞍点（甲基旋转 TS，
+130.5 meV，−247.5 cm⁻¹）。所以 `U_basin ≈ 0` 而 `U_TS > 0`。

判据不是"幅度小于基本域宽度"，而是**到最近对称副本的距离**：2.09 rad = 119.75°，距 2π/3 只有
0.25°，形式上在域内、物理上是恒等操作。

### 1.4 mode-resolved 收敛判据

固定力容差只控制"力有多小"，不控制"离驻点多远"：`|δq| ≈ |g_q|/k_q`。乙醇上扭转比键伸缩软约
150 倍（0.21 eV/rad² vs 33–53 eV/Å²），所以按实测刚度分级：

```
k > 50 eV/rad²        stiff               全局 fmax 已管住
0.01 ≤ k ≤ 50         soft                要求 |g/k| < 1°，沿自身路径 Newton 步
|k| < 0.01            free                记录，不要求收敛
k < -0.01             negative curvature   交给极小检验，polish 不碰
```

实测效果：乙醇 gauche 角度的批次间散布从 ±4° 收到 ±0.15°，开销约每试验 2%。后来补了
`flat_biased` 一档 —— 曲率在地板下**而梯度不为零**的坐标不是平的（一个甲基转子因此冻住
9.97 meV）。

### 1.5 三层身份，不可互相折叠

| 层 | key | 商掉什么 |
| --- | --- | --- |
| 化学态 | `chemical_key` = 图哈希 + 碎片 + 配置立体化学 | 构象被商掉 |
| 反应事件 | `reaction_event_key` = 键差 `B⁻=E₀\E₁`、`B⁺=E₁\E₀` | 在 `Aut(G₀) ∪ Aut(G₁)` 下正则化，方向商掉 |
| 通道 | `reaction_channels`，允许自环 | 退化反应**就是**自环，不是构象变化 |

**准入裁决权归 `chemical_key`，不归分类标签**：等价 H 转移会得到 `reactive` 标签但相同的
key，此时绝不开新节点。第二层被折叠回第一层发生过一次（P1 的退化转移），后果见
`PRRS_STATUS.md` §8.19.5。

立体化学不推断唯一 Lewis 结构，而是推断**一个 admissible bond-order assignment 集合**
（Kekulé 对返回恰好 2 个解，环上键正确地不被 lock）。无法解析时抛 `IdentityUnavailable`，
**绝不偷猜**。

### 1.6 两层网络与双预算

```
G_chem = (V_chem, E_rxn)
  └── 化学节点 C_i（由 chemical_key 判定）
        └── 构象水库 {R_i1 … R_iK}
```

同一节点内新增构象**不扣 chemical 预算** —— 反应可及性依赖构象：`P(B|A,conformer)` 而不是
`P(B|A)`。水库评分 `S = w_E S_E + w_D S_diversity + w_R S_response`（能量玻尔兹曼因子、可旋
转扭转上的圆周距离、以及在线响应项，未搜索者给乐观先验）。

### 1.7 曲率：Hessian 是验证极限，不是主观测量

有限幅度响应是对 PES 的 **sinc² 低通滤波**（已定量确认），只有 `a → 0` 那个极限才有谱定理。
有限幅度的本征值**不得**称为"频率"、**不得**报 cm⁻¹。

MACE 的力本身是 autograd 梯度，所以解析二阶导本来就在图里 —— 解析 Hessian 已是默认
（`hessian_source = "auto"`），比 FD 快 2 倍且到机器精度。刚体模地板门（6 个平凡模的 Rayleigh
商，零额外力评估）是永久的无参考质量监控；转动曲率正比于残余梯度（`ẋᵀHẋ = ω²(g·x⊥)`），
扣掉后残差 ~1e-16。

### 1.8 证据分类记录，不合并

| 字段 | 含义 | 强度 |
| --- | --- | --- |
| `attempts` | 这次试验自己从源走到了目标 | 直接观测 |
| `continuations` | probe 被继续交付后到达 | 较弱 |
| `ts_support` | 鞍点两侧下降 —— **无向** | 结构证据 |

TS 连接是无向结构性证据，**不自动生成有向边**。`arrival_through_saddle` 用的是"一阶鞍点恰好
在两个盆地的边界上"：两侧同名 ⇒ 该点邻域属于那个态。min-mode following 必须给**目标**方向 ——
不给就爬全局最软模（3-氧代丁醛上返回 8 个 −65…−80 cm⁻¹ 的转子鞍点，而质子转移在 −3186）。

### 1.9 失败必须关闭（fail closed）

判不出来就拒绝并说明理由，绝不用默认值替代。**存在/生效必须由执行证实，不能由记录证实**：
`importlib.metadata` 对 conda 装的 `openmm-ml`/`openmm-torch` 假阴性，所以改用 `find_spec` 探
模块，结果诚实标为 `"installed; no version metadata"`；`create_mace_context` 不信
`removeCMMotion=False` 这个**参数**，而是建完 System 后清点 force 列表。

### 1.10 职责边界（用户 2026-09-02 界定）

**PRRS 只负责 TS 在哪里、产物是什么。** 势垒和虚频的绝对值是势能面的属性，不是搜索方法的属
性，交下游 QM。苯甲酰丙酮三模型对照支持这条界定：`chemical_key` / `event_key` / `broken` /
`formed` / saddle order **逐字段一致**，r(O···O)@TS 差 0.036 Å；而**势垒差 3.7 倍**
（+316.3 / +84.8 / +96 meV）、**虚频差 3.3 倍**（−3086 / −926 / −1436 cm⁻¹）。

"好不好用"的干净拆分：

```
A_θ[E_m]  vs  T[E_m]     ← 搜索方法好不好用（同一张 PES 上比）
T[E_m]    vs  T[E_QM]    ← 模型准不准（交下游 QM）
```

把 `A_θ[E_m]` 直接和 QM 势垒比，会把搜索误差、旋钮误差、模型误差混成一个数，无法归因。

---

## 2. 当前状态

实现 `5b69dfed41c4bdb2`（2026-09-05）· `233 passed, 1 skipped` · digest 绿（一条已审查迁移）
· `src/prrs` 17 个模块 / 约 8.4k 行，14 个测试文件。

### 2.1 已冻结的结论

| 阶段 | 体系 | 打的轴 | 结果 |
| --- | --- | --- | --- |
| **P0** | 乙醇 | 构象发现 + 架构，无反应 | 13/13，逐字段复现；频率由独立的 `ase.vibrations` 路径交叉验证（`validation.py`，**搜索侧从不 import**） |
| **P1** | 丙二醛 | reaction event，**退化**自环 | `runs/p1_round8` **13/13**，与 `p1_round7` 逐字段相同；`event_key a04410b76e16b117`；4 个 TS，−3189.6…−3186.7 cm⁻¹；势垒 **419.14 meV**；`reactions = 0` —— **这是对的**，退化转移两端是同一化学态 |
| **P2** | 3-氧代丁醛烯醇 | 网络扩张 A→B，`k_A ≠ k_B` | `runs/p2_round7` **16/16**；T=300 的 16 个 seed 全部完成、0 失败；Q1（有向边落在转移通道上）16/16，95% Wilson [0.81, 1.00]；势垒 **+361 meV**；ΔE(B−A) −11.567…−10.758 meV |
| **P3** | 苯甲酰丙酮 | 芳香图 + 螯合 + 竞争异构体 | 判据 2026-09-03 冻结并判定通过（`docs/P3_BENZOYLACETONE_ACCEPTANCE.md`）；势垒 +316 meV、虚频 −3086 cm⁻¹；TS 处 r(O···O) 从 2.5289 收缩到 2.3234 Å（重原子门控有实测值） |

两条贯穿的不变式：

```
不同极小  ≠  不同 chemical state
conformer discovery  ⇏  chemical budget consumption
```

**冻结时明确不成立的七条**（不要引用成结论）：n=16 的 16/16 不是"接近 1"（95% 下界 0.81）；
B1b 对 4 个 seed 记"不评"；Q1 从 12/16 升到 16/16 是旁证不是隔离归因；0.81 meV 是端点散布不是
能量分辨率；直接观测与续行是两类强度不同的证据不可合并；seed 71 的 C4/D1/D2 是真实失败且与
ceiling 无关；openmm 8.6.0.dev 下 backend parity 为何变红仍未查明（现在绿是因为回滚到 8.5.2）。

### 2.2 QM 复核（已通两例）

协议：**r2SCAN-3c 出几何与频率 → ωB97M-V/def2-TZVPD 在同一几何上补单点。** r2SCAN-3c 单独用
把 SN2 中心势垒错 11 kcal —— **几何好、能量必须补**。校准：SN2 上复现文献 0.4–0.9 kcal。

- **P1**（`docs/experiments/p1_qm/`）：双向 IRC 收敛、两端 0 虚频、正向端点对极小 RMSD
  0.0004 Å、`event_key` 复现冻结值。两张面的鞍点位置几乎重合（对称感知 RMSD **0.0151 Å**），
  所以相对 `T[E_QM] = 134.4 meV` 的那 285 meV 差是**纯能量的、落在模型上**，不是几何、不是
  方法。13/13 没有被推翻。
- **b04**（`docs/experiments/b04_qm/`）：从 PRRS 的 TS 出发精化到 2.3645/2.3645，与独立求得的
  QM 鞍点小数点后四位相同；双向 IRC 收敛、两端 0 虚频互为镜像；`event_key 847306d3f8c57bd4`
  与 PRRS 一致。

两例同一形状：**几何与事件身份属于方法且正确，能量与曲率属于模型。**

### 2.3 基准套件（`benchmarks/`）

一格 = 一个反应体系 + 一组**跑前冻结**的判据 + 一个机械评定它的脚本。判据先于数据冻结，
preflight 先于判据。

| id | 体系 | 打的轴 | 状态 |
| --- | --- | --- | --- |
| b00–b03 | 乙醇、丙二醛、3-氧代丁醛烯醇、苯甲酰丙酮 | 身份层与网络层的能力 | `fixture` —— 判据冻结过、判定通过过，验收职责移交 P4+ 的 QM 协议 |
| b04 | Cl⁻ + CH₃Cl | **电荷** · 分子间 · 背面进攻 | `preflight` —— 5/5 门通过，PRRS 找到的 TS 已 QM 复核通过；**仍缺正式判据与机械 `check.py`** |
| b05 | 叔卤代烷 + 亲核试剂 | 取代/消除竞争通道 | `blocked`（等 b04） |
| b06 | 亲电加成底物 | 极性加成 | `blocked` |
| b07 | 环己烯（逆 DA）/ 丁二烯+乙烯（DA） | **碎片数 1↔2** · 重原子坐标 | `blocked`（等 b03）；探针已绿 —— MACE-OFF24 在 C–C 断键区光滑，那个 −579 meV 台阶是扫描伪影 |

b00–b02 打的其实是**同一条化学轴**（分子内质子转移、碎片数恒为 1、势垒 419/361/316 meV 同一
带、虚频 −3186/−3010/−3086 cm⁻¹ 同一带）。让这种重复在加格之前就看得见，正是这个套件存在的
理由之一。

**从 P4 起默认模型是 `MACE-POLAR-1-M`**（用户 2026-09-04），依据是两组 QM 对照而不是偏好：
POLAR-1 在中性有机质子转移上偏 +12 meV、在离子 SN2 中心势垒上偏 +2.1 kcal，是唯一两边都成立
的候选。P0–P3 冻结在 MACE-OFF24 上，不重跑、不改标签，且**不自动**成为 P4+ 的参照。

### 2.4 已知限制与开着的缺口

- **离子化学被 MACE-OFF24 阻塞**：它对总电荷完全不响应（charge 0/−1/+1 逐位相同，模型里没有
  任何 charge/spin 模块），而域门会按构造拒绝每一个离子端点。charge-aware checkpoint 已在机器
  上，还剩三处代码缺口（`benchmarks/b04_sn2_chloride/README.md`）。
- **键拉伸的交付不是片段刚体的。** 拉 C–Cl 0.600 Å 时只有 C 与 Cl 动，每根 C–H 缩短
  0.0628 Å；CH₃ 若随 C 刚体平移则变化 0.000000 Å。扭转已用 `rotate_fragment` 避开，键没有等价
  物。**是位移实现，不是必要耦合** —— 这就是下一步。
- **bracket 不独立于 collateral 门。** 细化只比相邻样本，0.6 处一次 `collateral` 拒绝就隔断了
  唯一有效的 (0.3, 1.0) 对。次序是先修扰动交付，再评估扫描策略。
- **非键接近家族在密堆分子上交付不了**（P1：12/60 试验被 `collateral_bond` 挡掉，全部是
  `compress` 作用在非键对上，共价键被拖动 0.08–0.41 Å 而预算是 0.05 Å）。门是对的，交付满足
  不了它。
- **方向排序仍是启发式**，应改为 `U(P) = (U_basin, U_reaction, U_TS)`。
- **大体系没有 `n_-` 判定路径。** 完整 Hessian 是 6N 次力评估，projected Lanczos 只需数十个
  `Hv`、crossover 约 15 原子；但解析 Hessian 在 N ≲ 30 上便宜 3–4 倍**且严格**，所以第一件该做
  的是让 `hessian_spectrum` 改用 `calc.get_hessian()`，不是实现 Lanczos。
- **mode-resolved 判据只覆盖 torsion**，其他软模需要 projected-Hessian 路线。
- **`free` 与 `stiff` 两档只有解析测试覆盖** —— 乙醇的所有扭转都落在 soft 带
  （0.19–0.69 eV/rad²）。
- **`pulse` 家族未验收**（外功与撤力切换的连续性未测）。
- **`S_response` 只有乐观先验那一半在起作用。**
- **苯甲酰丙酮上 OFF24 的 `side +1` 到不了 A**，即使 `quench_steps = 5000`；记为**那张势能面的
  已知限制**（OFF24 在 A 处最软内坐标模 0.00262，POLAR-1 是 0.00470），并保留"两个模型差的不
  只是这一个数，仍非隔离归因"这句限定。
- **b04 仍缺正式判据与机械评定**（`CRITERIA.md` / `check.py` 是空壳）。

### 2.5 下一步（按用户定的顺序）

```
1  可复跑环境与版本记录、原版回归              ✅ env/BASELINE.md
2  三项正确性问题                              ✅ HANDOFF §12
3  单独修扰动交付，比较有效扰动比例与 bracket   ← 下一步
4  扫描策略 · b04 正式验收 · 扩展 B 组
5  效率用「每个经 QM 确认的独立通道成本」衡量
```

第 3 步的具体做法：对适用的桥键按分子图切分片段、整体移动，保持片段内部几何；环内键与跨碎片
接近分别处理。之后把扫描显式分为**返回 / 离开 / 被阻断 / 未知** —— 不得把拒绝当返回，也不得
无记录跳过。

---

## 3. 文件地图

```
src/prrs/                17 个模块、约 8.4k 行；搜索核心（前六个）不依赖积分器
  config.py              全部协议参数 + 校验（跑前拒绝不合法配置）
  internal.py            primitive internals、精确刚体旋转、continuation 位移、精确 ΔK 冲量
  chemistry.py           三层身份、自同构群、对称感知 RMSD、classify_transition
  perturbations.py       六家族提案、对称商、分层方向预算
  network.py             Registry：接纳、微观态水库、域门
  state.py               活性子集图编码、Kabsch RMSD、响应分类
  reliability.py         每次力评估的失败关闭门控
  runner.py              单次试验、mode-resolved 淬火、Hessian/极小检验、鞍点下降、路径记录
  search.py              预算调度、幅度细化、TS pool、manifest
  calculators.py         calculator 工厂（含 Committee —— MACE-OFF 单模型下不可用）
  openmm_backend.py      经实测验证的 openmm-ml/MACE 原语（交叉校验后端）
  validation.py          独立交叉校验（ase.vibrations）—— 搜索侧从不 import
  io.py / cli.py         持久化与命令行
tests/                   14 个文件、206 个测试函数
docs/                    状态、交接、验收文档、曲率路线、experiments/、figures/
benchmarks/              反应基准套件（装的是问题，不是数据）
scripts/                 p3_node_bundle/（分片重活）、扫描状态、provenance、实现哈希
runs/                    442 MB；哪些计入结论见 runs/README.md
models/MACE-OFF24_medium.model   随仓库走，sha e5ccf5837f685899
```

---

## 4. 怎么跑

```bash
make help          # 各 target 的作用
make check         # ruff + pytest + digest，按此顺序
make hash          # run manifest 会记的 implementation_sha256
make digest        # 黄金样本门：动 src/prrs 之前和之后都要跑
```

```bash
prrs demo   --output runs/demo --seed 42
prrs search input.extxyz --output runs/my_run --config config.json \
            --calculator mace --calculator-kwargs calc.json
prrs path   --run runs/my_run --attempt <id> --output seed.extxyz --images 9
```

`prrs path` 导出未精化的路径种子（`purpose = unrefined_path_seed`），供下游 NEB/dimer/string
精化 —— 淬火迭代被排除，因为它们不是动力学。

`make digest`（= `python docs/experiments/digest_regression.py`）对 332 文件 / 2136 帧的黄金样
本比对，并把差异**分类**：`graph_invariant` / `parity_only` / `key_only`。有意的改动用
`--allow parity_only` 之类放行 —— **但放行只确认差异的形状，不确认新值是对的**，那需要化学
判断。

### 两个环境，有意为之

| | |
| --- | --- |
| `openmm_dev`（conda） | 科学环境，它的版本是每个结果的一部分 |
| `.venv-tools` | 只放 ruff，这样格式化永远不会推动 openmm/torch/mace |

backend parity 曾因 openmm 的一个 patch 版本变红，所以**纯开发工具一律不许装进
`openmm_dev`** —— 别在里面 `pip install ruff`；缺失时 `make tools` 重建。

```
python 3.12.13   Linux-7.2.0-1-cachyos   mace-torch 0.3.16   torch 2.12.0   ase 3.29.0
numpy 2.4.3      openmm 8.5.2            openmm-ml / openmm-torch / nnpops：无 version metadata
本机 GPU  RTX 2080 Ti 11 GB   单点 52.2 ms（12 原子，OFF24）
节点 GPU  RTX 5080   16 GB    单点 23.8 ms
```

并发是实测出来的：2080 Ti 天花板 3.2×，而**更快的** 5080 只有 1.7× —— 这排除了 GPU 计算是争抢
源；绑核 + 单线程毫无变化，排除了 CPU 线程争抢。**不要开 MPS（实测 0.15×）。** 杠杆是多卡，
不是单卡堆并发。

---

## 5. 项目的规矩

这些是用户定的，不是建议。违反它们产出的结果会被作废，本项目已经作废过一轮。

- **判据跑前冻结，逐字不改。** 要改必须作为**显式修订**记录：旧文本原样保留、改的理由、以及
  **修订发生在看到哪些数据之后**。
- **不许"边跑边改"。** 绝不为了得到绿色结果而换 seed 或调参。
- **干预必须分别命名**，这样改进能归因到具体某一项。
- **未完成 ≠ 失败。** 分母是**发起数**；未完成记未完成（0 分因为没回答），崩溃记崩溃并附原
  因。`scripts/scan_status.sh` 机械地区分这三态。
- **旋钮先验证再设计实验。** 在确认被改的变量能推动被测机制之前设计方差实验，等于用算力买一
  个已经注定的答案。（这条是 T=0 扫描 8/8 那次教训的产物。）
- **门（gate）不能因为它挡住了想要的结果就调。** 要动必须先量清楚它挡掉的是什么。
- **测量优先于推断。** 本项目有大量"我推断 X，检验杀了它"的记录，那些记录是资产不是耻辱
  （`HANDOFF.md` §6）。
