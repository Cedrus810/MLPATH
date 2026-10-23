# 扫描优化：曲率复用、常驻模型与批量控制器（2026-10-23）

用户本轮明确指示修复每点丢失曲率，并允许模型常驻、多个小体系同时推进，
要求检查控制器耗时。该指示更新此前“跨候选 batching 留到 v2”的执行范围：
本轮实现并单独验证性能原语，不等同于已完成 M5 的发现调度、两账户或恢复。

已读发现层设计 §1–11、评估协议 §1–6；约束扫描仍遵守发现层 §4.1，
每结构计费遵守评估协议 §2.1。本轮不变更扫描幅度、fmax、maxstep、
每点步数上限、能量上限、门拒绝或认证算法。

## 实现

- `explore/src/explore/probe.py`：`ScanSettings.reuse_curvature` 默认开启；
  设为 `False` 可运行此前每点重置的基线。每点仍新建 BFGS，只传递上点
  已收敛时的曲率矩阵；不传递位移、力差或 secant 历史。
  新约束梯度单位向量为 n，P=I−nnᵀ，初值为
  `H0=P H_previous P + alpha nnᵀ`，alpha 使用 ASE 自己的初值。
  这保留新切空间的曲率，法向重新初始化；零步点无已学矩阵时自然退回原初始化。
- `explore/src/explore/parallel.py`：`ScanJob`、`scan_many` 与
  `MACEBatchEvaluator`。模型由一个线程独占，在同一对象内常驻，
  可以跨多次 `scan_many` 调用复用；各扫描由独立 CPU 控制器推进，
  每个有自己的原子、约束、优化器、calculator 缓存和账本。
  一轮收集每个活跃控制器的下一次请求，合成 disjoint graph 一次求值，
  完成者从待办队列补位，无定时凑批窗口。
- 电荷、spin、head 沿原 MACECalculator 的 KeySpecification 转入每张图，
  保留原单位转换与 stress 设置；只输出扫描所需的无偏能量和力。
  不支持 committee 或编译模型；C2/Hessian 仍走原 calculator，不伪装成批量认证。
- 测量：`batch_rounds`、请求等待墙钟、控制器线程 CPU 秒、evaluator 外关键路径墙钟，
  以及组图、模型计算、结果回传的同步计时。控制器 CPU 用 `time.thread_time`，
  仅主线程开 cProfile；不把阻塞等 GPU 的时间报成 CPU 忙。
- 费用在批次实际 dispatch 前按结构逐个扣除；预算耗尽项不送入模型。
  完成回调负责回收 worker，包括初始化或清理阶段的异常；协调器中断时唤醒所有等待者。

## 验证设计与状态

解析测试势是旋转/平移不变的三原子耦合键势，约束极小能量与反力有解析式。
验证了曲率复用在该势上减少真实调用，并保持解析极小、反力与投影力条件；
不同长度的并发扫描、预算耗尽、失败任务和尾批能独立闭合。
测试曾抓到客户端把约束一并送给通用 evaluator、导致返回投影力的问题；
现在发送模型的快照去掉约束，扫描端自行投影。

本机 `make -C explore check`：格式/lint 通过，**116 passed**。
日志：`benchmarks/explore_smoke/results/2026-10-23/scan_throughput_check.log`。

真实模型脚本 `benchmarks/explore_smoke/scan_throughput.py`：

1. 一个进程只加载一次 POLAR-1；启动导入、模型加载分别计时。
2. b04 阴离子、同几何阳离子诊断与 b09 中性结构，单独、合批、重排合批分别比较
   能量/力，检查不同原子数和电子态是否串扰。它们只是后端诊断，不是新反应结果。
3. b04 18 点旧优化器一条；曲率复用的同一组四条扫描先串行再四条合批。
   第 i 条仅对输入整体旋转 17i 度，串行与合批用相同坐标；不是独立科学重复。
4. 各扫描独立 500 forward 预算与段 reservation；保存全部点、成本和失败。
   比较逐点能量、坐标、反力和键图，不以快但换支冒充等价提速。
5. 主线程 profile 开启，耗时包含剖析开销；只代表这次性能冒烟。

节点 183，2080 Ti；OMP/MKL 各 4 线程。启动前 GPU 0%、显存 1 MiB。
整个进程 `timeout --kill-after=10s 300s`，不可续跑；逐组保存，超时保留未决。
初次输出：`benchmarks/explore_smoke/results/2026-10-23/scan_throughput/fp32/`；
修复后同协议重试：同目录下 `fp32_fixed/`。

## 实测结果

**首次诊断失败保留：** 导入 45.711 s，模型构造/加载 98.736 s，混合电荷/原子数
合批诊断完成；能量差均为 0，力差最大 1.116e-3 eV/Å（fp32，未独立量重复噪声）。
没有完成扫描组。线程栈显示 worker 已退出、主线程仍等 active 项，GPU 0%。
本机最小复现确认当前 Python 下多个 cProfile 同时 enable 抛
`ValueError: Another profiling tool is already active`；原实现把剖析初始化放在
worker try 之前，未获取其 Future 异常，故没有回收 active。

已按确切 PID 1034512 发送 SIGTERM，最终 supervisor 退出码 137，三个相关 PID
均已消失，GPU 已释放；未补写为正常退出或成功扫描。
初版相关源码归档在 `scan_throughput/first_source/`，保留原始日志与元数据。
修复为主线程单 profiler + 完成回调兜底；新回归覆盖启动异常和主协调器中断，
相关 26 项测试通过。修复后完整包检查与真实对比继续单列，不能用新源码追认旧结果。

### 修复后：完整冒烟通过（183，fp32）

`fp32_fixed/measurement.json` 与 `scan_throughput/SUMMARY.json` 已在 183 上读取。
退出码 **0**、总墙钟 174.165 s；其中导入 30.160 s、模型加载 69.078 s。
现场确认本轮 PID 已退出、GPU 0% / 1 MiB，没有残留作业。

| 组 | 扫描数 | 墙钟 / s | 结构 forward | 实际模型调用 | 控制器线程 CPU 合计 / s | evaluator 外墙钟 / s |
| --- | --- | --- | --- | --- | --- | --- |
| 每点重置 BFGS，串行 | 1 | 17.537 | 158 | 158 | 1.273 | 1.365 |
| 曲率复用，串行 | 4 | 37.161 | 84 / 83 / 81 / 80（共 328） | 328 | 2.814 | 3.080 |
| 曲率复用，四路合批 | 4 | 16.600 | 83 / 82 / 76 / 82（共 323） | 83 | 3.212 | 3.323 |

同一起点第 0 条扫描：forward **158→84（−46.8%）**，probe 账本段
**17.482→9.527 s（约 1.83×）**。四条相同输入的串行与合批墙钟比
**2.239×**；两组 forward 数略不同，按实际结构数归一化吞吐仍为 **2.205×**。
所以收益不只是减少了 5 次求值，也不是把整批冒充成一次 forward。
三组逐结构模型单位与账本总数分别为 158、328、323，全部闭合。
这是有主线程 cProfile 的一次性能冒烟，不是重复统计或新反应发现比较。

9 条扫描共 **162/162 点收敛**，比较中的每点键图完全一致。
曲率复用相对旧基线：最大能量差 1.953e-3 eV、原子位移 0.00347 Å、
反力差 0.00307 eV/Å。合批相对串行：最大能量差 3.906e-3 eV、
原子位移 0.00861 Å、反力差 0.00942 eV/Å。
混合电子态/原子数的固定几何诊断中能量差均为 0；单元素 batch 自身重复的力差
2.75e-4 eV/Å，混合/重排 batch 最大力差 9.88e-4 eV/Å。
保留这些差异，不宣布逐位等价；没有通过放宽 fmax 或门来取速度。

**控制器定位：** 原扫描约 7.8% 墙钟在 evaluator 之外；合批后这一比例约 20.0%。
合批组的 evaluator 内，组图/传输约 2.499 s、模型计算（同步完成）10.642 s、
结果回传约 0.0145 s。`outside_evaluator_wall_s` 含调度、优化器、guard、
账本 I/O 与等待，不等于纯 CPU 时间；线程 CPU 秒另列，不能与墙钟相加。
在这组小体系上，先降优化步数、再合批确有收益；控制器残余约 3.3 s 仍可优化，
但不能把整个 16.6 s 都归因于控制器慢。

当前交付为真实模型 b04 工程冒烟通过。尚未验证多通道搜索收益、
全局预算/恢复、长期负载或任意 batch 组成下的科学等价；批处理入口由调用方显式选择。
实测五个相关源码文件已归档 `scan_throughput/fixed_source/`，逐个核对与 manifest
来源一致。实测后仅补了 `reverse_check` 对 `reuse_curvature` 的透传，
确保显式旧基线在反扫时仍是旧基线；该设置透传不在本次 GPU 路径内。
最终 `make -C explore check`：格式/lint 通过、**118 passed in 15.35 s**，
日志 `scan_throughput_check_final.log`；新脚本语法与文档目标/空白检查通过。

## 使用入口

```python
from explore.parallel import MACEBatchEvaluator, ScanJob, scan_many

# calc 在进程启动时建立一次；jobs 为独立的 ScanJob，各自带 Ledger。
backend = MACEBatchEvaluator(calc)
outcomes, timings = scan_many(jobs, backend, max_active=4)
# 后续批次继续用 backend，不重建模型。
```

`reuse_curvature=False` 留作单独基线，合批规模由调用方显式指定。
当前接口是扫描执行原语，尚未接入 M5 未施工的全链控制器或精确全局预算。
