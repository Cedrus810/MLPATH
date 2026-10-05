# geomeTRIC 局部几何与 TS 精修

2026-10-03：已提供可复用 Python 适配层与 `prrs optimize` 独立入口。
定位是 **PRRS / PARS 提供初猜 → geomeTRIC 局部优化 → 驻点与连接验证**。
PRRS 的 quench 和 PARS 的批量 `refine` 仍使用各自现有流程；本入口尚未自动串入它们。

## 放在哪一步

- **端点优化**：反应物、产物或构象的局部弛豫；随后按原来的身份和极小确认规则验收。
- **TS 精修**：PRRS 的 TS 候选或 PARS CI-NEB 的最高图像作为初猜，使用 `--transition`。
  这是局部鞍点搜索，不保证从任意几何发现指定反应的 TS。
- **验证**：优化收敛之后检查 Hessian 的负模数，再做双侧下降／IRC 与预期 A/B 核对。
  `result.json` 始终记录 `stationary_point_certified=false`、`connectivity_verified=false`，
  避免把优化终止当成反应连接证明。

## 安装与使用

geomeTRIC 是可选依赖：`pip install -e ".[geometric]"`。
OpenMM XML 路线还需要 OpenMM：`pip install -e ".[geometric,openmm]"`。
现有 `openmm_dev` 环境已经含 geomeTRIC 1.1.1，不需要重复安装。
MACE/OpenMM 路线复用项目已有的 openmm-ml/MACE 环境，两个 extras 不包含 openmm-ml。

已有 ASE calculator（包括 MACE、QM 或带可靠性守卫的 factory）：

```bash
prrs optimize endpoint.extxyz --output runs/endpoint_geometric \
  --calculator prrs.calculators:mace_factory --calculator-kwargs calc.json

prrs optimize ts_guess.extxyz --output runs/ts_geometric \
  --calculator prrs.calculators:mace_factory --calculator-kwargs calc.json \
  --transition --convergence-set GAU_TIGHT --maxiter 300
```

已有 MACE/OpenMM 部署入口：

```bash
prrs optimize ts_guess.extxyz --output runs/ts_openmm_geometric \
  --mace-model /absolute/path/model.model --device cuda --platform CUDA \
  --charge -1 --multiplicity 1 --transition
```

标准 OpenMM System／ForceField XML + PDB：

```bash
prrs optimize molecule.pdb --output runs/mm_geometric \
  --openmm-xml system.xml --platform CPU
```

XML 需要是本地文件；System XML 按原内容反序列化，ForceField XML 使用 PDB topology
创建 `NoCutoff`、无约束、无 CMMotionRemover 的有限体系。原子顺序必须与 System 一致。
普通固定拓扑 MM 力场的化学范围由该力场决定，不能据此验证 MLP 反应中的成键／断键。

支持 `--coordsys tric|dlc|cart`（默认 TRIC），`--convergence-set GAU|GAU_TIGHT|GAU_VERYTIGHT`。
输入必须只有一个几何；完整的 PARS 路径要先明确提取端点或 CI 图像，不能默认取最后一帧找 TS。
极小优化使用近似 Hessian；TS 优化首先用力差分计算笛卡尔 Hessian，约需 `6N` 次额外力评估。
大体系请根据这项成本设置预算；`--maxiter` 是优化迭代上限，不是力调用上限。

## 输出与状态

输出目录必须不存在；不会覆盖旧 run。输出包括：

- `input.extxyz`：输入快照。
- `optimization.extxyz`：优化器记录的迭代几何、能量和力，**不是动力学轨迹**，也不包括所有 Hessian 差分点。
- `geometry.extxyz`：最后一个优化器记录的几何；未收敛时也可能存在，必须连同状态读取。
- `geometric.log`、`scratch/`：优化与数值 Hessian 的日志／中间数据。
- `result.json`：`CONVERGED`、`NOT_CONVERGED` 或 `FAILED`；配置、版本、输入／实现指纹、
  backend 来源、优化迭代数、力评估数、最后几何的能量与最大原子力。

只有 `CONVERGED` 的 CLI 返回 0；其余返回 2。输出生成前的输入／后端构建错误打印到 stderr。
数值非法、calculator 拒绝等计算异常会记录为 `FAILED`；迭代预算结束是 `NOT_CONVERGED`。
保留 caller 的 charge/spin 数组与结构 info；ASE 路线的电荷配置遵循 calculator 自己的契约。
`--charge`／`--multiplicity` 向 ASE input info 或 MACE/OpenMM System 传递，必须与输入已有的
`charge`／`spin` 一致。MACE/OpenMM 省略时优先读输入，输入也未声明时采用 0／1。
XML 路线以 XML 定义的 Hamiltonian 为准，不接受这两个覆盖参数。

## Python 复用与边界

```python
from prrs.geometric_backend import GeometricConfig, openmm_evaluator, optimize_geometry

result = optimize_geometry(
    atoms, openmm_evaluator(atoms, context), "runs/local_ts",
    GeometricConfig(transition=True), metadata=backend_metadata,
)
```

ASE 路线用 `ase_evaluator(atoms, calculator)`。接口的 evaluator 返回 eV 和 eV/Å 的力；
适配层转换为 geomeTRIC 的 Bohr、Hartree 和 Hartree/Bohr 梯度，并处理 `gradient=-force`。
OpenMM 路线复用传入 Context，包括平台与动态 MLP 邻域，不强制 Reference，不做平台回退。
geomeTRIC 生成的距离拓扑仅用于优化坐标，不作为 PRRS 化学身份的权威。

当前仅支持至少两个真实原子的有限非周期体系；拒绝 ASE constraints、OpenMM constraints、
virtual sites 和周期 System。优化中只内置形状／有限值检查；不自动套用 PRRS 的 OOD 与
力阈值。需要这些门时，ASE factory 应返回配置好的 `GuardedCalculator`。
同一个 OpenMM Context 会被优化更新位置，不应并发共享；本适配层的日志也面向串行调用。

上游资料：[geomeTRIC OpenMM 与 CustomEngine 文档](https://geometric.readthedocs.io/en/latest/engines.html#openmm)。
官方原生 OpenMM engine 使用 Reference 平台；这里用 CustomEngine 桥接现有 Context，保留平台选择。

## 本次工程验证

2026-10-03 在本机 `openmm_dev` 中执行：最终适配层测试 **13 passed**，包含解析双阱的
极小／TS 搜索、原子单位梯度有限差分、OpenMM Context 与 System／ForceField XML 路线、
水的隐式刚性约束拒绝／关闭、预算失败、电荷冲突和完整路径误用拒绝。
根测试全量运行 **551 passed, 2 skipped**；全量收集之后的适配层改动另由最终 13 项测试覆盖。
ruff 格式与 lint 通过；冻结身份 digest **332 文件、2136 帧、4 对事件零变化**。

本地 `MACE-OFF24_medium.model` → OpenMM CPU Context → geomeTRIC 的水分子极小优化
实跑收敛，3 次优化迭代、4 次力评估。DLC／Cartesian 也各跑过一个 ASE/EMT 小分子检查。
这些是接口工程验证；尚未完成真实反应 TS、CUDA 优化或大体系性能验收。
