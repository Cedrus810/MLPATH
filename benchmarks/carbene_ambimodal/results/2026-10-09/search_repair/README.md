# 搜索设计修复的验证记录（2026-10-09，Codex 产出，原放在 /tmp）

说明见 [docs/SEARCH_DESIGN_REPAIR_2026-10-09.md](../../../../../docs/SEARCH_DESIGN_REPAIR_2026-10-09.md)。
这些文件原本在 181 的 `/tmp/mlpath-*`，2026-10-09 15:39 原样拷进仓库（`cp -p`，保留时间戳）；`/tmp` 原件未删。

| 文件 | 内容 |
| --- | --- |
| `b09-audit-20261009.py` | 生成 `../b09_probe1_audit.json` 的核查脚本 |
| `coarse-first-{red,green,focused}.log`、`coarse-first-baseline-digest.log` | 方案 B 修复前失败 / 修复后通过 / 定向测试 / 历史 digest 基线 |
| `restricted-modes-red.log`、`constrained-delivery-{red,green}.log` | 受限软模、保键交付的失败 → 通过记录 |
| `mechanisms-core.diff` | 三处核心源码改动的 diff |
| `mechanisms-final-focused.log`、`mechanisms-full-check.log` | 定向 64 passed；`make check` 592 passed, 2 skipped |
| `mechanism-mutations.{py,json}`、`mutations/*.log` | 三组变异检验（log 按变异名改了文件名；json 里 `log_original` 是原 /tmp 路径）。脚本里的输出路径仍是 /tmp，未改 |
| `coarse-first-pars-check.json`、`pars_check_run/` | PARS 文件接口检查及其合成 run（`run_original` 是原 /tmp 路径） |
