# 已提出、暂不实施的小缺口 · 2026-10-04

用户顺序：先修搜索时窗与多鞍点验证的大缺口，小缺口先提出、后修。
本清单记录当前已复现的问题；本轮没有修改两个 provenance 实现。
没有证据表明历史 run 已因这两项污染。

1. **PRRS 绑定方法 factory 漏掉所属实例状态**。
   `calculator_identity(Factory(1).build, {})` 与 `Factory(2).build` 相同且
   `verifiable=True`，即便 `build` 用 `self.scale` 构造不同 PES。
   待修：同时记录 `__func__` 与 `__self__` 的可核验状态，处理循环引用；
   用真实 kill→resume 拒绝与同条件等价回归验收。
2. **PARS 模型 state_dict 之外的行为参数未覆盖**。
   同一 `torch.nn.Module` 类、相同 buffer/weight，但 `self.scale=1/2`，
   `forward(3)` 分别返回 3/6，calculator identity 仍相同。
   待修：记录模型类型/代码及影响 forward 的可核验非权重配置，无法核验则拒绝；
   加入权重相同、行为不同的续跑拒绝回归。不能宣称仅 state_dict 足以标识 PES。

外部资源/模型仍应显式声明来源和内容哈希。大缺口实验使用明确固定的解析势参数，
真实 b04 时长回放显式冻结 POLAR-1 checkpoint 内容；不使用上述特殊 factory/模型
形态作为条件相同的依据。
