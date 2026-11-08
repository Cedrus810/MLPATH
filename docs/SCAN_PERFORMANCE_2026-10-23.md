# 约束扫描性能挖掘（2026-10-23）

范围：用户明确本轮关心当前扫描速度；误发的“节省 VRAM”回复不作为需求。
检查对象是 `explore.probe.constrained_scan` 与 POLAR-1 求值链。该记录不是 M5
编排施工，也不改变发现层的科学参数、计分口径或冻结验收。

## 已核对的成本与现状

| 证据 | 已有观测 | 对扫描提速的含义 |
| --- | --- | --- |
| `benchmarks/explore_smoke/results/2026-10-13/run/ledger.jsonl` | b04 18 点、157 forward、段内 19.274 s，约 8.72 forward/点、0.123 s/forward | 已有点间几何延续；下一步应分别降低每次求值耗时与每点求值数 |
| `benchmarks/explore_smoke/results/2026-10-22/factory_reuse/factory_reuse.json` | 一次 b04 **认证** 39.49 → 31.27 s；199 forward 不变；factory 调用 13 次 | 是认证收益；不能给单条扫描直接套用 21%，扫描在循环外只调用一次 factory |
| `benchmarks/explore_smoke/results/2026-10-22/fp32_vs_fp64/compare.json` | 单次求值 b04 0.278 → 0.111 s、b09 0.266 → 0.191 s | 是固定几何求值收益；尚不能声称整条扫描快 2.5×/1.4× |
| `docs/PLAN_EXPLORE_V1_2026-10-13.md` 的 10-22 执行约束 | 已同意搜索 fp32、认证 fp64，长寿命进程与 calculator 复用 | 决策已落盘；旧 `b04_scan.py` 仍用默认 fp64/fresh factory，M5 编排尚未接入 |
| `benchmarks/explore_counts/results/2026-10-22/COUNTS.md` 与 DONE_LOG | b09 全枚举成本来自构象配对 × 姿态 × 意图的乘积，用户已同意预算内抽样 | 内核数倍提速不能单独解决组合规模；抽样属于 M5 的另一个工作项 |

上述 10-22 实测由 kasuga02 写入。当前读到了 NFS 上的既有工件与完成记录，
未重跑旧实验；新的冒烟状态和结果直接在 183 上读取。

## 源码发现与优先级

1. **无需求的 stress 计算：已做小样本测量，收益小，优先级下降。** 当前安装的
   `openmm_dev/lib/python3.12/site-packages/mace/calculators/mace.py`
   `calculate()` 按 `model_type` 无条件为 PolarMACE 设置 `compute_stress=True`，
   并不依照本次 ASE 请求是否包含 stress。扫描只读取能量和力。
   `PolarMACE.forward` 把该选项传入 `prepare_graph` 与 `get_outputs`，
   关闭它能去掉应变相关求导工作，但真实速度和力差异必须测量。
   本轮只做 benchmark 内的 forward 关键字包装，未改 site-packages 或生产 factory。
2. **量每点重置 BFGS 的额外 forward。** `probe.py` 在每个 delta 内重新
   `BFGS(work, ...)`；上个点的几何保留，曲率估计丢弃。当前 ASE 的 `initialize`
   使用对角 `H0=alpha*I`。可单列“携带 Hessian 预条件”的试验，保持目标网格、
   maxstep=0.05 Å、fmax、最大步数及 guard 不变。不能直接把跨目标的
   位移/投影力差当成同一目标下的 secant；换目标时至少重置该 secant 历史，
   并核对新约束切空间。收益未知，可能换支或更难收敛。
3. **再分解单次求值。** 分别计邻居图/ASE 转换、POLAR 自洽与反向求导、
   GPU→CPU 返回，以及优化器与 guard 的 CPU 时间，再决定线程数或后端优化。
   不根据低 GPU 利用率就断言“算力不足”或“Python 是瓶颈”。
4. **多意图同步批处理是吞吐方向。** 同一条扫描后一点依赖前一点的松弛结果，
   不能把这一串点当成独立输入一次并行求解；不同意图的当前优化步可以组合。
   需要独立的优化器/约束/门拒与逐结构 forward 计费，也要保存不同步完成和失败。
   上位 v1 暂不包含 batching，本轮只列候选，不暗中接入生产编排。
5. **compile/CuEq/OEq 后置。** 当前 MACE 提供入口，但安装、模型转换、首次编译、
   动态邻居图和 POLAR 长程支路会影响可用性与摊销。未测量，也未修改生产环境。

## 本轮有界冒烟

脚本：`benchmarks/explore_smoke/scan_performance.py`。一进程只用一种 dtype。
benchmark 内临时包装 model.forward 的 `compute_stress` 参数，退出上下文即恢复；
保留原 MACECalculator 的图构造、电子态映射、单位转换和结果处理。

- 固定几何：b04 电荷 −1、b09 电荷 0，多重度均为 1；固定 seed=0 的
  0.01 Å 坐标扰动；两路先预热，8 对求值交替 ON/OFF 与 OFF/ON 顺序；
  每次清 ASE 缓存，CUDA 同步，记录逐对能量/力差和逐次墙钟。
- 整条扫描：b04 C0–Cl5，18 点、每点 −0.1 Å，fmax=0.01 eV/Å、
  maxstep=0.05 Å、最多 300 步/点、max_rise=3 eV，与原冒烟相同。
  ABBA 顺序（stress on/off/off/on），每条预算 400 forward；
  全部点、未收敛、门拒、预算耗尽与已付成本都保留。
- 节点：183（kasuga02），2080 Ti，OMP/MKL 各 4 线程。
  启动前 GPU 0%、显存 1 MiB，无 GPU 计算进程。
- 单进程外层 `timeout 180s`；不可续跑，结果目录要求新建；
  每个固定几何体系、每条完成扫描后写 measurement.json，账本逐段写盘。
  若进程被硬终止，当前段的未闭合记录与超时状态单独保留，不能算成功。
- 输出：`benchmarks/explore_smoke/results/2026-10-23/scan_performance/`。
  JSON 记录实际包版本、源码内容指纹、模型路径、设备与线程。指纹只用于来源记录。

## 测量结果与下一步

在 183 上直接读取结果并生成 `scan_performance/SUMMARY.json`，确认两个冒烟的
supervisor/Python PID 均已退出，没有遗留 GPU 作业。

| 固定几何、fp32，每体系 8 对求值 | stress 开，中位秒/次 | stress 关，中位秒/次 | 比值 |
| --- | --- | --- | --- |
| b04（6 原子） | 0.07905 | 0.07521 | 1.051× |
| b09（29 原子） | 0.07670 | 0.07435 | 1.032× |

这只是一次小样本、固定几何的中位数比较。原始 stress-on 数据分别有 0.87/0.83 s
慢调用，全部保留；未删点，也未证明 3–5% 是可稳定复现的收益。
关闭 stress 的数值等价性也未闭合：b04 最大能量差 0、力差 3.19e-4 eV/Å；
b09 最大能量差 1.953e-3 eV、力差 1.609e-3 eV/Å。
本轮没有 b09 同后端、同几何的重复噪声对照，不能把全部差异归因于某个选项，
也不能据此宣布无行为变化。**暂不把关闭 stress 接入生产。**

整条 b04 stress-on 基线完成：**18/18 点收敛、158 forward、15.024 s**
（包含扫描调用内的账本开闭；probe 段自身 13.451 s）。每点 BFGS 步数
`6,7,7,8,8,8,8,8,8,8,6,5,8,8,10,11,10,6`，平均 7.78 步、
8.78 forward/点。与 10-13 fp64 的 19.274 s 不在同一次配对实验里，
不能把两者比值当成精度切换的净收益。

**启动与未完成保留：**

- 第一次在 180 s 窗口后退出；最后活跃时进程为 `DNl`、等待通道
  `rpc_wait_bit_killable`，元数据已写但 calculator 加载未完成；没有求值结果。
  未捕获到退出码，不能补写成已知 124，状态另存 `float32_startup_status.json`。
- 第二次同协议重试，calculator 构造/加载用 **105.961 s**，不含先前的顶层导入。
  两体系固定几何测量与第一条扫描完成；第二条 stress-off 扫描已启动但账本未闭合；
  余下两条未启动。外层超时退出码 **124**，完整 ABBA 比较没有完成。
- 这条未闭合扫描的实际 forward 数未知，不能记为 0；由该扫描预算给出上界 400。
  首版脚本未填段 reservation，因此不能依赖其账本重放恢复该数。
  未做第三次重试，也未延长三分钟上限。

**收口：** 当前最值得后续验证的是减少每点优化步数；随后量单次求值内部构成、
再决定独立意图 batching 或后端加速。stress 开关目前只有小幅 latency 差，
不作为主提速路线。M5 仍需把已经同意的搜索 fp32 与长寿命 worker 接到调用端，
才能摊销实测很大的冷启动。全枚举规模则按已经同意的预算内抽样解决。

**交付与验证边界：** 生产源码、环境、默认精度、科学阈值均未改。
新增 benchmark 与本文；ruff 格式/lint、语法检查通过。
实测脚本完整归档为 `scan_performance/measured_scan_performance.py`，
已核对与测量 JSON 的源码指纹一致。之后仅改善测量脚本的可观测性：
45 s 周期栈回溯、启动阶段日志、扫描段 reservation=400；这些改善经静态检查，
**没有重新上 GPU 验证**，因此不能用新版源码冒充这批结果的实际实现。
