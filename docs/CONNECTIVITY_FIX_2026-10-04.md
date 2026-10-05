# 精修鞍点连接验收 · 2026-10-04

基线 `14d8325` 的局部一阶 Hessian 检查不证明该鞍点连接指定两端。
本轮落实 [上轮未关闭的连接门](MAJOR_GUARDS_2026-10-04.md)，
保留 PRRS/PARS 模块边界、源 run、旧结果与冻结物理协议。
两个身份指纹小问题仍在 [延后清单](DEFERRED_FIXES_2026-10-04.md)。

## 实现与拒绝边界

CI 通过一阶曲率检查后，保存投影 Hessian 的负模并还原为归一化笛卡尔向量。
分别用正负号扰动鞍点；每个起点必须低于原鞍点能量。然后独立执行带回溯的
梯度下降，不沿原 band、不使用未来标签、也不把准牛顿驻定当作下降证明。
接受的每步必须严格降能，最大原子位移受限；不收敛或无法找到降能步即拒绝。

落点达到力阈值后，另做端点同规格的 Hessian 检查。两侧均无解析负模，且
能以唯一的一对一映射匹配 reactant/product，才返回 `CI_NEB_CONVERGED`。
比较使用固定原子映射；孤立分子按配置 Kabsch 对齐，外场不去除刚体运动。
同盆地、第三盆地、歧义映射、落点仍有负曲率、未完成的下降或力预算耗尽均
`FAILED`，不报告势垒，保留 band 与部分下降记录。不回退 plain NEB 掩盖失败。

显式不确定度阈值先检查整个 band，再应用于下降每次 trial；缺测/超阈值为
`UNRESOLVED_MLP`，非有限或负的不确定度为 `FAILED`。这些拒绝也不报告势垒。
无阈值不编造误差区间。plain NEB 保留原有分辨率保护与离散 band 估计，
明确 `connectivity_check.status=not_checked`，不冒充已验收的鞍点连接。

`connectivity_check` 保存方法、策略、正负扰动起点、接受步能量序列、回溯次数、
最后位置、力残差、落点谱、到两端 RMSD、匹配标签与不确定度覆盖。
有独立 source/quench 盆地 witness 时，scope 为 `source_basin_witnesses`；
旧种子没有 witness 时仅为 `refined_endpoints_only`，不能据此补认证旧源身份。

默认参数在完整预算基准前固定，没有根据覆盖率调整：

| 参数 | 默认值 |
| --- | --- |
| `refinement_descent_displacement` | 0.05 A，总笛卡尔模长 |
| `refinement_descent_step` | 0.05 A²/eV，初始梯度步长 |
| `refinement_descent_max_displacement` | 每原子最多 0.05 A/trial |
| `refinement_descent_fmax` | 0.001 eV/A，与端点阈值取更严格者 |
| `refinement_descent_max_steps` | 每侧最多 2000 个接受步 |
| `refinement_descent_max_backtracks` | 每步最多 20 次折半，另含初始 trial |

参数纳入 `ParsConfig`、refinement conditions 与结果策略回显；非法步长、非有限
阈值、非整数或负预算拒绝。零步预算可用于保留未完成的诊断，不能绕过验收。
所有下降 trial 和每落点 6N 次差分力计算复用原 calculator 及其硬预算。
预算中断不生成成功结果，保留已完成的分支诊断。

这是给定 MLP、固定参数与数值容差下的笛卡尔下降连接证据，不是质量加权 IRC，
也不是 QM 验证。弱曲率、平坦区域、刚体设置不适用或不够精确的 CI 点可能拒绝；
没有扩大扰动或调阈值以追逐成功率。多步反应经过中间盆地应拆分/重新分析，不能
把局部鞍点标为直接连接整个通道的证据。

## 实际反例与验收

解析外场 `E=x²(x²−1)²+(y²+z²)/2`：极小 x=-1,0,+1，
一阶鞍点 x=±1/√3，谱 `[-8/3,1,1]`，局部势垒 4/27 eV。
要求直接连接 -1 与 +1 的均匀种子，得到的局部鞍点实际连接其中一端和中间极小。
实际导出的 `14d8325` 源码返回 `ci_neb_converged`、0.148148135572 eV；
当前实现保留两侧落点并拒绝，势垒为 null。不是以模拟旧逻辑代替实际旧源码。
[反例重跑脚本](../benchmarks/pars_connectivity/run.py) 与
[旧/新完整记录](../benchmarks/pars_connectivity/results/2026-10-04/current.json) 已保存。

- PARS `make check`：160 passed；格式、lint、独立模块边界通过。
- 独立导出模块：159 passed、1 skipped；唯一跳过项依赖不随模块分发的历史 run。
- 根目录独立参考/硬预算测试：5 passed；新增下降阶段真实 calculator 预算耗尽回归。
- 回归覆盖真正连接、源 witness scope、负模符号交换、第三盆地、两侧同盆地、
  歧义匹配、低力负曲率落点、未完成下降、非降能步、不确定度缺测/超阈值与非法策略。
- 隔离变异跳过连接验收、接受歧义映射均被针对性回归击穿；工作区不改写。

双鞍点基准重新完整运行，冻结势能面、候选池、各臂选择计划、20 个随机种子、
每臂两个候选和 6000 次力计算上限均与上轮一致。manifest 在参考/oracle 前保存
所有有效配置、源码/runner 哈希与选择计划。12 个全池候选均具备两个单调下降的
正曲率落点，唯一匹配两个源 witness；谱、下降与落点信息随 oracle 单独保存。

| 策略 | 鞍点覆盖 | 力计算 |
| --- | --- | --- |
| PARS | 2/2 | 2616 |
| channel-pruning | 2/2 | 2616 |
| cheap-energy | 1/2 | 2805 |
| random，20 种子 | 平均 0.85，14/20 全覆盖 | 每臂最多 4827 |

PARS 相比上轮 1926 次增加 690 次，全部下降与落点检查计入原硬预算。
共享初始化 252、参考认证 2578、全池 oracle 22164 次另列；总墙钟 46.28 s。
这仍不是优于 channel-pruning 的覆盖证据，不是分子或 QM 验证。
新结果独立保存在
[`2026-10-04-connectivity`](../benchmarks/pars_two_saddles/results/2026-10-04-connectivity/report.json)，
旧结果保持原样。

```sh
make -C pars check
PYTHONPATH=src:pars/src /home/ruigengji/miniforge3/envs/openmm_dev/bin/python \
  -m pytest -q tests/test_two_saddles_reference.py
PYTHONPATH=pars/src /home/ruigengji/miniforge3/envs/openmm_dev/bin/python \
  benchmarks/pars_connectivity/run.py --output runs/new_connectivity_counterexample.json
/home/ruigengji/miniforge3/envs/openmm_dev/bin/python \
  benchmarks/pars_two_saddles/run.py --output runs/new_connectivity_validation
```

现有精修记录的 conditions 会因实现与新参数改变而拒绝普通 resume。
需重新 select 建立新验收目录，或在保留旧结果后明确 `pars refine --redo`。
本轮没有重写历史分子精修或改变原 source run。

仍待研究验收：真实分子多机制的同预算覆盖/成本、独立 QM 势垒，以及整网冻结
判据下的搜索时长/停止规则。当前工程连接门不能关闭这些研究门。
