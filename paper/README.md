# MLPATH 论文初稿

2026-11-13，v0.1。英文正文已经成文；这一版用截至 11-12 的开发批次建立文章叙事，供逐段修改。计算规则仍按原运行记录保存，11-12 的通用规则选项尚未采纳。

## 从这里读

- [英文正文](MLPATH_MainText.md)：摘要、引言、方法、四节结果、讨论、结论、参考文献。
- [SI 初稿](MLPATH_SI.md)：输入与参数、证据路径、评分、数值未决项、计费边界、QM 证据和复现命令。
- [横向验证设计](HORIZONTAL_VALIDATION_DESIGN_2026-11-13.md)：三开发＋三留出槽位、六体系局部/四体系端到端、共同 C2 与预算边界；设计稿，未开计算。
- [参考文献 BibTeX](references.bib)：五篇已核对的基础方法/模型来源。
- [机器可读证据](data/draft_evidence.json)：当前两张主表、原始来源和镜像修正。
- [图 1：局部定位](figures/Fig1_local.pdf)；[图 2：预算—召回](figures/Fig2_recall.pdf)。另有同名 PNG。

## 文章现在讲什么

从反应物姿态出发，利用约束扫描形成局部分支，再用反力括区与 E/F 优化定位 TS。b04、b05 的局部对照展示了求值成本优势；N1 漏掉协同 DA、Dimer 找到额外取代通道，展示了搜索覆盖的边界。

当前标题暂用 **MLPATH: Energy-and-force transition-state search through constrained reaction branches**。全文将“候选几何”“一阶性质”“两侧连接”“QM”分开报告。首篇聚焦 HF TS 搜索；PARS 选择收益、mature、b09 分岔动力学留在各自证据链里。

## 这一版已经处理的关键数字

| 项目 | 初稿口径 |
| --- | --- |
| b05 局部 HF 的 19 次参照匹配 | 其中 16 次一阶通过；正文严格命中数写 16 |
| N1 局部 Dimer | 采用已保存的镜像审计：3 次命中、覆盖 2 类 DA 参照 |
| b05 端到端 sB | 严格 2/3；第三个几何匹配但阶数未决 |
| b05 Dimer 的零命中 | 限于三个 E2 参照；同时写明发现额外取代鞍点 |
| N1 事件机制 | 后续下降与原 HC 有一项端点判法冲突，保留这一未决项 |
| 约 9.6 / 5.8 倍 | 两批局部搜索的总 forward 比值，包含失败；不写成整条流程墙钟加速 |
| 两个种子 | 逐个展示，未添加置信区间 |

## 按什么顺序修到下一版

1. **先改文章的主张和语言。** 这一版已有可读正文，可以开始讨论摘要是否准确、方法贡献是否清楚、结果是否足够支撑结论。作者、目标期刊与格式尚未指定。
2. **收口正式验证协议。** 一次确定通用数值测量、同键图脊规则、参考集及双方共同的认证成本边界。新的规则版本与现有开发批次分开。
3. **补最小验证。** 优先处理阶数重算不稳定与 Dimer fp32 数值适用性；再确定冻结后的重复、留出体系和关键消融。N1 负结果保留，不为 DA 增加参考位点特例。
4. **补化学与发布材料。** 对进入主张的新增通道补必要端点/QM 核查，完成相关工作检索、checkpoint 来源、环境与公开工件。正文 Methods 与结果随新证据同步更新。

上述是写作与证据收尾顺序，不是新计算的执行记录。本轮只写论文、提取既有汇总与作图。

## 图的使用

图 1、图 2 已由真实汇总生成，可以用于讨论与改稿。图 1 的命中是重复尝试数；图 2 的参照集遗漏参照外新通道。两个图都展示开发批次，图注已写明边界。后续可增加算法流程图与真实坐标的反应结构图；本版未用示意分子或虚构势能面代替计算数据。

## 复现

在项目根目录运行：

```sh
python paper/tools/build_draft_evidence.py
MPLCONFIGDIR=/tmp/mlpath-draft-mpl /home/ruigengji/miniforge3/envs/openmm_dev/bin/python paper/tools/build_draft_evidence.py --figures
```

脚本仅读取五份明确路径的汇总；每项费用求和与汇总值核对后生成表格，并将镜像修正作为显式叠加。它不调用模型，不更新历史结果。SI S1 列出全部生成文件。正文和 SI 的 Markdown 是本版主文件，未制作 Word/LaTeX/PDF 正文导出。

## 参考来源核对

基础文献与模型来源已在 2026-11-13 通过以下原始来源核对：

- [Dimer：作者机构记录](https://iris.hi.is/en/publications/a-dimer-method-for-finding-saddle-points-on-high-dimensional-pote/)。
- [CI-NEB：作者机构记录](https://iris.hi.is/en/publications/climbing-image-nudged-elastic-band-method-for-finding-saddle-poin/)。
- [FIRE：APS 论文页](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.97.170201)。
- [ASE：作者仓储记录](https://wrap.warwick.ac.uk/id/eprint/87141/)及[官方 Dimer 文档](https://docs.ase-lib.org/ase/dimer.html)。
- [MACE-POLAR-1：作者预印本 v1](https://arxiv.org/abs/2602.19411v1)。

局部 checkpoint 与公开发布版本的对应关系仍按已有项目证据记录，不用模型家族论文代替这项核验。
