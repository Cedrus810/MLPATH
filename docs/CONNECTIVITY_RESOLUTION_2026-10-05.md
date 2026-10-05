# 连接验收的步长敏感性 · 2026-10-05

继续审查基线 `7d299d5`，复现一个新的大缺口：**有限负模扰动与离散降能步骤
可以跳过窄中间盆地**。两侧最终匹配指定端点，并不足以排除此问题。
本轮补强 [上一轮连接验收](CONNECTIVITY_FIX_2026-10-04.md)；旧源 run/结果不改写。

## 实际反例

解析外场：

`E=(x²−1)²+(y²+z²)/2−0.02 exp(−((x−0.015)/0.002)²)`。

主鞍点在 x≈0，局部 Hessian 有一个负方向。独立力零点求解给出右侧先遇到的
极小 **x=0.01500600110285**，下一个鞍点 x=0.02007934646557；用能量差分确认
该中间极小正曲率，并检查主鞍点右侧到该极小之前的力方向。
力同时通过独立能量差分回归，未用“观察到的标签”构造选择指标。

旧 0.05 A 正扰动直接越过中间极小和下一个鞍点，后续每个接受步仍可降能，
最终到 x≈+1；负侧到 x≈−1。实际导出的 `7d299d5` 源码因此返回
`ci_neb_converged`、1 eV 势垒。不是模拟旧逻辑的替身。
新复核到达中间极小，拒绝该直接连接与势垒。
[反例脚本](../benchmarks/pars_connectivity/resolution.py) 和
[实际旧/新记录](../benchmarks/pars_connectivity/results/2026-10-05/current.json) 已归档。
这是解析外场的工程反例，尚未证明历史分子结果受影响。

## 修复契约

新增 `refinement_descent_resolution_factor=0.25`，必须在 (0,1) 内。
粗尺度双向下降均匹配后，再独立以更细配置执行双向下降：负模扰动、梯度步长、
每原子位移上限同时乘该因子；力/曲率阈值、端点判据和每侧接受步上限保持不变。
同一 calculator/硬力预算连续使用，不缓存精修结果，不重新获得预算额度。

更细检查未完成、到第三/同一盆地、配对歧义、与粗尺度符号标签不一致，均拒绝
CI 成功并清空势垒；缺测/超阈值保持 `UNRESOLVED_MLP`，不能继承粗尺度测量。
不确定度汇总包括已完成/部分完成的细尺度分支与原 band 的已知测量。
`resolution_check` 保存细尺度完整参数、分支、能量、落点谱、力残差、RMSD 和状态。
敏感或未完成标记 `sensitive_or_unresolved`，一致标记 `stable_on_tested_scales`。
新条件与实现指纹拒绝直接复用旧单尺度精修记录。

成功另写 **`continuous_connection_certified=false`**：两次有限尺度检验不能证明
连续梯度流、质量加权 IRC，不能排除更窄的未解析盆地，也不能认证 MLP/QM 精度。
这是固定策略下的数值敏感性门，不能把有限复核包装为数学连接证明。
因子在本轮预算比较前固定，没有为通过 6000 次预算而调参或扩大预算。

## 验收与实际成本

- PARS `make check`：170 passed，格式/lint/模块边界通过。
- 根独立参考/硬预算：6 passed；新增细尺度复核不能重置预算的真实 calculator 回归。
- 独立导出模块：169 passed、1 个依赖历史 run 的 smoke 跳过。
- 跳过细尺度复核的隔离变异被窄盆地反例击穿，工作区不改写。
- 回归包括窄盆地、力与能量一致性、普通连接尺度稳定、非法因子、细尺度缺测。

重新执行双鞍点全池 oracle 与 23 个预算臂。冻结势能面、12 候选、原选择计划、
20 随机种子、每臂两候选、**6000 次力上限**均与上轮一致；manifest 在参考/oracle
之前保存全部有效配置及源码/runner 哈希。所有未限预算 oracle 仍到两个已知鞍点，
各自具有粗/细两尺度、四条已验收单调下降分支与唯一端点匹配。

| 策略 | 预算内覆盖 | 力计算 |
| --- | --- | --- |
| PARS | 2/2 | 5464 |
| channel-pruning | 2/2 | 5464 |
| cheap-energy | 1/2 | 5781 |
| random，20 种子 | 平均 0.525，仅 1/20 全覆盖 | 多臂耗尽 6000 |

PARS 从上一轮 2616 次增至 5464 次，**额外 2848 次**全部计入原预算。
随机平均覆盖从 0.85 降至 0.525，因新增验收成本使许多臂无法完成第二个候选。
这是当前策略下的预算内结果，不能解释成搜索准确率改善或一般优越性。
PARS 与 channel-pruning 仍相同；固定解析池也不能代替分子发现召回。
共享初始化 252、参考认证 2578、全池 oracle 39246 次另列；总墙钟 71.57 s。

新记录单独保存于
[`2026-10-05-resolution`](../benchmarks/pars_two_saddles/results/2026-10-05-resolution/report.json)。
旧记录不升级为两尺度验收结果；历史单尺度连接仅保留当时有限分辨率的证据。

```sh
make -C pars check
PYTHONPATH=src:pars/src /home/ruigengji/miniforge3/envs/openmm_dev/bin/python \
  -m pytest -q tests/test_two_saddles_reference.py
PYTHONPATH=pars/src /home/ruigengji/miniforge3/envs/openmm_dev/bin/python \
  benchmarks/pars_connectivity/resolution.py --output runs/new_resolution_counterexample.json
/home/ruigengji/miniforge3/envs/openmm_dev/bin/python \
  benchmarks/pars_two_saddles/run.py --output runs/new_resolution_validation
```

后续大门仍包括真实分子多机制/同预算验证、独立 QM 势垒、整网时长与停止判据。
[已提出的小问题](DEFERRED_FIXES_2026-10-04.md) 本轮保持延后。
