# 继续审查：三项大缺口与修复 · 2026-10-04

基线为 `741f708`。按用户“大缺口先修、小缺口先提出后修”的顺序，
本轮修复 PARS 的三项可执行反例；PRRS 源 run 与既有结果不覆写。
[上轮验收](MAJOR_FIXES_2026-10-04.md) 与
[两个延后身份问题](DEFERRED_FIXES_2026-10-04.md) 保留原证据边界。

## 1. [P1] 不同缺失指标形成循环支配，清空前沿

旧 `dominates` 只比较两个候选共同具备的指标，导致每对比较用不同集合：

```python
vectors = {
    "a": {"x": 0., "y": None, "z": 2.},
    "b": {"x": 1., "y": 0., "z": None},
    "c": {"x": None, "y": 1., "z": 1.},
}
```

a 经 x 支配 b，b 经 y 支配 c，c 经 z 支配 a。旧前沿为空，三者均进入
冷归档，margin 为 infinity，进一步触发严格 JSON 输出拒绝。
这是支持缺失值的通用 API 反例；尚未证明当前分子语料实际发生过这种循环。

现在只有可用指标集合相同的候选可比较；不同集合互不支配。
同一规则覆盖 dominance、directed margin 与 objective cover，避免前沿保留后
又因部分指标相同而被代表压缩掉。每个非空可比组均保留前沿，组内被支配者
有有限 margin。三个反例均保留；同组重复候选仍能正常压缩。
这是保守的策略变更，混合缺失模式下代表数可能增加，旧 selection 需重做。

## 2. [P1] 无关候选的 sigma 把缺测候选变成“确定被支配”

旧 tiering 以整个比较群是否存在任意 sigma 决定能否硬排除；随后把比较对
缺失的 sigma 当作零。`good=(0,0)`、`bad=(1,1)`、`noisy=(2,2)`，
`margin_gamma=10` 时，仅给 noisy 增加 `sigma_x=0.1`，就会把完全没有区间
信息的 bad 从 reserve 变成 `certainly_dominated`，禁止 gamma 恢复。
默认 V1 不生成 objective sigma；问题在调用者提供部分区间数据时触发。

现在每一对都要求双方每个可用指标有完整区间；一方缺测或某指标缺测，
不能形成区间确定支配。其他候选的数据不能为该比较对提供许可。
完整且分离的区间仍可硬排除；明确声明的零宽区间有效，不自动补零。
保留原有两项显式策略：nominal 模式可硬剪枝；双方均无区间时，正 tolerance
可按用户设定的点估计分离策略排除。后者不代表测得的统计确定性；部分区间
不会被 tolerance 绕过。

## 3. [P1] CI-NEB 成功标签接受高阶驻点和未解析的负曲率

解析外场 `E=(x²−1)²+(y²−1)²+z²/2`，端点 `(-1,-1,0)` 与
`(1,1,0)` 均为正曲率极小；均匀对角路径让 NEB/CI 完全驻定在原点。
旧实现返回 `ci_neb_converged`、势垒 2 eV；原点的 Hessian 谱是 `[-4,-4,1]`，
有两个负方向。经 `(1,-1,0)` 的显式路径最高仅 1 eV。
因此，优化器收敛与留在原 band 内不足以认定可用的一阶鞍点。

现在接受 CI 后独立做中心力差分 Hessian；按同一 SE(3) 设置投影，要求
恰有一个低于负噪声阈值的特征值。新增参数：

- `refinement_saddle_hessian_step=0.001` A。
- `refinement_saddle_curvature_tolerance=0.0001` eV/A²。

谱、负模数、移除刚体模数、步长、阈值、image index 写入 `saddle_check`。
零个或多个解析负模均 `FAILED`，`barrier_eV=None`，保留诊断链与能量；
不退回 plain NEB 来掩盖已经发现的驻点问题。
CI 未接受时仍沿用已有 plain NEB 的分辨率保护，仅报告离散 band 估计。
每个接受的 CI 增加 6N 次力计算；成本会进入现有 calculator 预算。
新参数和代码指纹属于 refinement conditions，旧结果需重新执行 `--redo`
或新建 selection；不把旧成功记录自动升级成新验收结果。

这只是给定 MLP 上、指定差分步长与噪声阈值下的局部阶数检查。
输出明确 `connectivity=not_checked`；没有执行负模双向下降/IRC、没有 QM 验证。
弱负曲率与平坦驻点可能拒绝，需独立分析，不能通过调阈值追逐验收率。

## 验收与重跑

- PARS `make check`：144 passed；格式、lint、独立模块边界通过。
- 根目录 `tests/test_two_saddles_reference.py`：4 passed；含硬力预算失败。
- 新回归覆盖循环、有限 JSON margin、不同掩码 cover、无关 sigma、部分区间、
  完整区间、一阶成功、二阶拒绝、零个解析负模拒绝。
- 三个隔离副本变异（恢复交集比较、把未知 sigma 当零、跳过阶数检查）均被
  对应反例测试击穿；没有改动工作区源码。结果见
  [`mutations.json`](../benchmarks/pars_two_saddles/results/2026-10-04-guards/mutations.json)。
- 原 CLI 的三原子独立双阱样例具有三个同步负模；回归现在检查失败原因与
  空势垒正常写入，不再把该样例误记为成功。
- 重新执行相同冻结协议的完整双鞍点基准：passed，12 个候选，20 个随机种子，
  每臂最多两个候选和 6000 次力计算；没有重复利用 oracle 的精修结果。
  PARS 与 channel-pruning 均 2/2、1926 次力计算；cheap-energy 为 1/2、2091 次；
  random 平均覆盖 0.85，14/20 全覆盖。总墙钟 36.04 s。
  相比上轮 1890 次，PARS 多 36 次，正好是两个三原子鞍点的额外差分成本。
  初始化 252、参考认证 2578、全池 oracle 18024 次力计算另列。
  这仍不是优于 channel-pruning 的覆盖证据，也不是分子/QM 证据。

新基准完整记录单独保存于
[`2026-10-04-guards`](../benchmarks/pars_two_saddles/results/2026-10-04-guards/report.json)，
包括 manifest、reference、oracle、arms 与 report；旧结果保持原样。

```sh
make -C pars check
PYTHONPATH=src:pars/src /home/ruigengji/miniforge3/envs/openmm_dev/bin/python \
  -m pytest -q tests/test_two_saddles_reference.py
/home/ruigengji/miniforge3/envs/openmm_dev/bin/python \
  benchmarks/pars_two_saddles/run.py --output runs/new_guard_validation
```

## 仍需完成的大门

1. 精修鞍点的负模双向下降/IRC，与映射后的端点身份验收。局部一阶曲率与
   整条 band 端点正确，仍不能代替该鞍点的实际连接证据。
2. 真实分子多机制、同预算覆盖/成本与独立 QM 势垒；解析双鞍点只验收工程行为。
3. 整网冻结判据下的搜索时长与停止规则验收；已有 b04 配对回放仅一条 scan，
   且时长敏感，不是整网/静态盆地边界认证。

本轮没有声称这些研究/连接验收已实现。未改动已有冻结协议或历史分子结果。
