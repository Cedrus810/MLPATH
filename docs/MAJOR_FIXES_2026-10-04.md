# 大缺口施工与验收 · 2026-10-04

用户优先级：先处理大缺口；两个身份指纹小漏洞先提出、延后修复，见
[待修清单](DEFERRED_FIXES_2026-10-04.md)。依据为
[审查 §5–§6](REVIEW_2026-10-03.md)。保持旧 source run 与冻结化学判据不变。

## A. 搜索时窗：明确被测对象，并提供配对验收

普通搜索继续测量“给定 probe 与初始热种子，自由响应固定时长后淬火，落到哪里”。
它与从 transition frame 出发的首次进入 committor 不同；从源极小直接用首次进入
规则会先承诺回到源，不能不定义离开规则就替换。没有强行将首次进入搬入 trial。

新 manifest、network、每条 amplitude scan 写入 `window_semantics`：响应步数、
步长/时长、`horizon_stability=not_measured`、`static_basin_boundary=false`。
这次是产物解释契约改变，几何/能量/方向/停止规则保持原算法。

新增 `prrs audit-windows`：同一 source microstate、probe、种子、幅度集和淬火配置，
只改变 `response_steps`。所有选中 scan 的原时长端点必须先以绝对 1e-5 Å
复现，再执行任何额外时长；不复现立即写失败报告并拒绝继续。
原 run 的 manifest/network/attempts 文件哈希前后核对，不覆写源目录。
新 manifest 必须匹配 calculator identity；历史 run 依赖基线回放，绝不补写身份。
外部模型可用 `--calculator-metadata` 明确声明文件与内容哈希。
遗留 run 缺完整 calculator identity 时，端点几何复现只能核验这组有限观测，
不能补成历史 PES 的密码学身份保证；报告冻结的是本次明确声明的模型/环境。

验收对比全部 sample 的分类和到达的 microstate，不只比较 bracket 位置。
幅度点缺失、无样本、失败/瞬态/未识别端点均不认证；目标 B 变成 C 也判敏感。
仅完全一致的有限网格标记 `stable_on_tested_grid`，不能升级成静态盆地边界证明。
CLI 不认证或失败返回非零；失败仍保留带 provenance 的报告与已完成 trial。

```sh
PYTHONPATH=src:. /home/ruigengji/miniforge3/envs/openmm_dev/bin/python -m prrs audit-windows \
  --run runs/my_completed_run --output runs/my_new_audit \
  --calculator module:factory --calculator-kwargs /path/to/kwargs.json \
  --calculator-metadata /path/to/model_sources.json \
  --response-steps 100 400 --scans 0
```

scope 明确是所选 scan 的原幅度集；不能把局部 audit 当成整网冻结判据验收。
旧 b04 记录已展示 25/50/100 fs 时长敏感；本轮新入口的真实 CPU 回放与完整检查
结果见下方。要改变普通搜索停止规则、或者主张时长不变性，仍需独立冻结协议。

## B. 多鞍点：修正剪枝顺序并落地预算对照

实现缺陷：原 PARS 先在一个通道里全局 Pareto 剪枝，再对幸存者分管。
独特但廉价指标较差的路径可以在机理保留生效前变成冷归档。

现在默认 `preserve_mechanism_diversity=True` 时先对全部候选分管，
以同一通道归一化尺度在各管内计算 dominance/margin/tier，再做 active cover。
冷归档成员仍在该管里，可审计/作为后续备用，不进入当前代表集合。
`False` 保留原通道剪枝作为明确反事实。输出标注 `dominance_scope`。
`refinement_order` 先轮流访问每管一个代表，再访问更多代表；CLI 采用该顺序。
旧 selection 要重新 select，代表数可能增加；描述符本身仍是近似。

可运行的 [双鞍点基准](../benchmarks/pars_two_saddles/README.md) 与冻结协议已保存。
两个鞍点分别用力/Hessian/正负模下降确认连接同一对正曲率极小。
每策略在完全相同池、参数、候选上限与硬力计算预算下真实执行冷 NEB，
记录失败、端点验收成本、共享初始化、参考/全池 oracle 成本与总墙钟。
精修标签不进入任何选择规则；完整计划在精修前写入 manifest。
覆盖分母固定为两个已认证鞍点，不以选中或成功的结果重定义分母。

这关闭了“没有可运行多鞍点对照”和剪枝顺序的工程缺口。
真实分子多机制、发现召回、同预算普遍收益与 QM 势垒精度仍待验收；
解析势不能替代这些研究结论。

## C. 本轮证据与剩余门

CPU 环境 GPU 驱动不可用；未修改系统配置。

| 检查 | 结果 |
| --- | --- |
| 根 `make check` | 566 passed、2 skipped；格式/lint 通过；457.98 s |
| 新窗口/双鞍点相关测试 | 10 passed；含随后补入的严格力预算失败回归 |
| PARS `make check` | 136 passed；格式/lint/模块边界通过 |
| 历史 digest | 332 文件、2136 帧、4 事件对不变；254 次既有迁移比较 |
| 暂存快照导出 | 根新增 10 条通过；PARS 135 passed、1 个真实 connect run smoke 因不分发历史 runs 跳过 |

**真实 b04 回放完成**：POLAR-1、float64、CPU 两线程；源 `b04_frozen2` 的 scan 0，
7 幅度 × 25/50/100 fs = 21 次 trial，约 **2429.57 s**。
先完成 7 个原时长端点复现，再运行其余时长；原 manifest/network/attempts 哈希不变。
50 fs 的窗口 `[0.75, 0.80]` 在 100 fs 移为 `[0.30, 0.60]`，0.80/1.00 出现回越。
25 fs 的 0.80、100 fs 的 0.75 为未确认瞬态，不能认证为 confirmed departure。
总体 **not_certified**、label_changes=true，拒绝静态盆地边界解释。
25 fs 汇总与旧 geometry-only 探索脚本不同，原因是新门不把 unconfirmed 当 confirmed。
这是一条 scan 的配对测量，不能代替整网冻结判据验收。
完整小型 [真实审计报告](../benchmarks/b04_window_audit/results/2026-10-04/report.json)
及 [复现说明](../benchmarks/b04_window_audit/README.md) 已纳入版本管理。

解析双鞍点最终运行通过，约 47.53 s（含共享准备、参考认证、全池 oracle 和所有臂）。
两个参考势垒为约 0.600573512 / 0.675102168 eV；12 候选全池覆盖两鞍点。
所有策略上限为 2 候选、6000 次力计算；精修成本包括端点/参照 Hessian。

| 策略 | 鞍点覆盖 | 本次实际力计算 | 本次失败/未匹配 |
| --- | --- | --- | --- |
| PARS 管内剪枝 | 2/2 | 1890 | 0 |
| 原通道剪枝反事实 | 2/2 | 1890 | 0 |
| 廉价能量 excursion 排序 | 1/2 | 2055 | 0 |
| 随机 20 个预定种子 | 平均 85%；14/20 次覆盖两鞍点 | 每臂独立记录 | 逐臂记录 |

两候选 PARS 精修势垒约 0.600573512 / 0.675103274 eV；与独立参考匹配。
共享池生成 252 次力计算、参考认证 2578 次、全池 oracle 17808 次，分别列账，
没有把验证 oracle 的成本隐藏成在线方法收益。相同预算是共同上限，实际消耗不必相等。
本池原通道剪枝同样保住两鞍点，因此不能宣称本池证明剪枝修复提高了覆盖。
独特路径廉价指标完全较差的专门回归和隔离变异检查验证剪枝问题；
本池的对照只证明多鞍点验收链路可运行及这组有限工程结果。

小型完整记录纳入版本管理：
[manifest](../benchmarks/pars_two_saddles/results/2026-10-04/manifest.json)、
[reference](../benchmarks/pars_two_saddles/results/2026-10-04/reference.json)、
[oracle](../benchmarks/pars_two_saddles/results/2026-10-04/oracle.json)、
[arms](../benchmarks/pars_two_saddles/results/2026-10-04/arms.json)、
[report](../benchmarks/pars_two_saddles/results/2026-10-04/report.json)。

隔离临时副本中的三项变异均被抓住：取消管内剪枝、只比较窗口而忽略目标身份、
放宽基线端点复现到足以掩盖错误。源工作区没有写入变异。
原始日志 `/tmp/mlpath-major-root-check.log`、`/tmp/mlpath-major-pars-last-check.log`、
`/tmp/mlpath-major-new-tests-final.log`、`/tmp/mlpath-major-mutations.json` 为临时材料；
关键结果已保存于本记录与版本管理的基准 JSON。

仍需独立验收的真实分子门：

1. 在确定模型/参数上建立至少两个同端点、独立驻点/双侧连接验证的真实分子鞍点。
   b05 当前局部产物区驻点不能自动充当不同反应机理参照。
2. 冻结代表预算/力调用预算与端点、Hessian、射击、选择、精修的总成本，再比较多种策略。
   当前解析池不证明真实分子成本收益、发现召回或普遍机理保真。
3. 若普通搜索改成首次进入/首次离开，需要独立定义驻留与离开语义，并整网复测冻结判据。
4. MLP/QM 势垒误差仍在；上述几何和连接证据不认证 QM 能量或速率。
