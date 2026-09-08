# PRRS 新人导读

这套文档是给**第一次接触这个仓库**的人写的。目标：读完能回答三个问题——
这个项目在解决什么问题、它的数学对象是什么、怎么跑一次并正确地读结果。

仓库里既有的文档（`PRRS_STATUS.md`、`PRRS_CURVATURE_ROUTES.md`、各 `P*_ACCEPTANCE.md`）
是**工作记录**，按时间顺序写、假定读者已经在场。这套是**入口**，按逻辑顺序写。

## 阅读顺序

| 文件 | 内容 | 大概花多久 |
| --- | --- | --- |
| [01-问题与框架.md](01-问题与框架.md) | 要解决什么、为什么不用 NEB、PRRS 的五步循环 | 15 min |
| [02-内坐标与扰动.md](02-内坐标与扰动.md) | 扰动的数学：内坐标、梯度、守恒律、retraction、精确冲量、对称商 | 45 min |
| [03-曲率与响应算子.md](03-曲率与响应算子.md) | 曲率的数学：$Hv$、质量加权、投影、$\kappa_a$ 的滤波解释、Morse 指标、TS 认证 | 60 min |
| [04-化学身份与反应网络.md](04-化学身份与反应网络.md) | 图编码、颜色细化、Lewis 枚举、立体化学、两层网络、预算 | 40 min |
| [05-怎么跑与怎么读结果.md](05-怎么跑与怎么读结果.md) | CLI、配置、输出文件、门控、已有结果与**不成立的说法** | 30 min |

先读 01 和 05，再回头读 02–04，也是合理的路径。

## 30 秒版本

给一个**反应物结构** $A$ 和一个**机器学习势**（MACE 等），不给产物、不给反应坐标、
不给路径。系统被结构化地扰动，观察它的响应，自动发现它能做哪些反应：

$$A \xrightarrow{\ P_k\ } \text{短 MD} \longrightarrow \{B_1, B_2, \ldots\}$$

输出是一张**反应网络** $\mathcal G_{\rm rxn} = (V_{\rm basin}, E_{\rm reaction})$，
外加通过认证的一阶鞍点（TS 候选）。

传统方法为"每次力评估都很贵"而优化；MLP 让力评估变得极便宜且高度并行，所以
PRRS 反过来用**大量廉价的动力学实验**换信息。

## 代码地图

搜索核心（前六个模块）**不依赖积分器**，ASE 只在结构/分析/执行层出现。

```
src/prrs/
  config.py         冻结的协议参数 + 校验；全部进 manifest
  internal.py       primitive internal coordinate、精确刚体旋转、continuation 位移、精确 ΔK 冲量
  perturbations.py  六个扰动家族的提案、对称商、分层方向预算、二维交付验收
  chemistry.py      chemical_key（自同构不变图哈希 + 碎片 + 构型立体化学）、自同构群、对称感知 RMSD
  state.py          活性子集图编码、Kabsch RMSD、响应分类
  network.py        两层网络（chemical node / conformer microstate）、reservoir 评分、双预算
  reliability.py    每次力评估的 fail-closed 门控
  runner.py         单次试验编排、mode-resolved 淬火、Hessian/极小检验、min-mode 爬升、鞍点下降
  search.py         预算调度、幅度细化、TS pool、manifest
  openmm_backend.py 经实测的 openmm-ml/MACE 部署原语（交叉校验后端）
  validation.py     曲率层的独立交叉校验（ase.vibrations）；**搜索从不导入它**
  io.py / cli.py    持久化与命令行
```

## 一条纪律，先说在前面

这个仓库对**"测到了什么"和"能声称什么"分得很开**，几乎每个记录里都带一个
`meaning` 字段说明这条数字**不是**什么。新人最容易犯的错是把下面这些混起来：

- 找到一个一阶鞍点 $\ne$ 证明了 $A \to B$（还需要两侧下降各自落到 $A$ 和 $B$）；
- 临界幅度 $a_c$ $\ne$ 活化自由能 $\Delta G^\ddagger$；
- 搜索找到了这些反应 $\ne$ 这就是全部反应（提案基不张成整个 $3N$ 空间）；
- 曲率的 $\sigma$ 是**回归残差尺度**，不是标定过的不确定度，方向错了它看不见。

05 章有完整清单。
