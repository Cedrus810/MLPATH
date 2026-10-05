# b04 怎么跑

**2026-09-04 的这一次是探针，不是验收** —— `CRITERIA.md` 尚未冻结。
套件规矩允许判据冻结之前跑探针（`docs/P3_BENZOYLACETONE_PROBES.md` 就是这个形状，
378 行里没有一条是验收结果）；跑它是为了知道判据该怎么写。

## 输入

`input.extxyz` 由 `build.py` 确定性生成，与探针用的逐字节相同：

```
python benchmarks/b04_sn2_chloride/build.py benchmarks/b04_sn2_chloride/input.extxyz
```

Cl⁻ 在 C 的背面 3.6 Å（`pair_cutoff_A = 4.0` 是提案层能看到的最远距离；
更近就已经是复合物了，等于把产物递给搜索）。总电荷 −1 写在结构的 `info` 里，
`prrs.search._stamped_with_charge` 会与 config 核对，冲突则拒绝开跑。

## 命令

`config.json` 带 `//` 注释，`prrs` 的 `load_json` 不吃，先剥掉：

```bash
python - <<'PY'
import json
s = open("benchmarks/b04_sn2_chloride/config.json").read()
d = json.loads("\n".join(l for l in s.splitlines() if not l.lstrip().startswith("//")))
json.dump(d, open("/tmp/b04_config.json", "w"), indent=1)
PY

PYTHONPATH=$PWD:$PWD/src python -m prrs search \
    benchmarks/b04_sn2_chloride/input.extxyz \
    --output runs/b04_probe \
    --config /tmp/b04_config.json \
    --calculator benchmarks.b04_sn2_chloride.calculator:polar1
```

`PYTHONPATH` 要同时含仓库根（为了 `benchmarks.` 这个包路径）与 `src`
（`prrs` 未安装，见 `docs/SESSION_2026-09-04.md` 的技术债一节）。

**计算都在 183 上跑**（用户 2026-09-09），所以实际用的是这一行：

```bash
scripts/prrs_remote.sh search benchmarks/b04_sn2_chloride/input.extxyz \
    --output runs/b04_probe3_rigidbond \
    --config .scratch/b04_probe3_config.json \
    --calculator benchmarks.b04_sn2_chloride.calculator:polar1
```

它替你做的是 `cd`、`PYTHONPATH`、`PRRS_MACE_MODEL` 与 `setsid` ——
四样里有两样是先踩过才写进去的，见 `env/BASELINE.md` 的"计算放在哪"。

`.scratch/` 是草稿目录（不进版本），所以**那份 config 的存档是
`runs/<run>/config.json`** —— 它记的是解析后的完整配置加上环境与实现哈希，
probe2 与 probe3 的"逐字段相同"就是拿这两份比的。

## 为什么工厂在这个目录里

`prrs.calculators.mace_factory` 调 `MACECalculator` 时不传 `model_type`，
所以 POLAR-1（checkpoint 类是 `PolarMACE`）会被当成普通 MACE 载入，
**算错而不报错**。在有东西真正走过这条路之前不去拓宽 `mace_factory`，
工厂就放在 `calculator.py`，经 `--calculator module:function` 这个既有钩子进来。

## 三次探针

| run | 变量 | 结果 |
| --- | --- | --- |
| `runs/b04_probe` | 首次 | 64 trials · 9 次 collateral 拒绝 · brackets 0 · ts 0（`seed_above_ceiling`） |
| `runs/b04_probe2_openceiling` | `saddle_search_max_rise_eV` 1.5 → 10.0 | ts 0 → 1（−535.65 cm⁻¹）· brackets 仍 0 |
| `runs/b04_probe3_rigidbond` | 实现树（交付修复 + bracket 记录修复） | **0 拒绝 · brackets 1 `[0.75, 0.80]` · ts0000 复现（−535.82）** |

三次的 config 除第二列写明的那一项外逐字段相同。逐条数字在 `PROBE_2026-09-04.md`。

## 2026-09-04 那次探针

`runs/b04_probe`，MACE-POLAR-1-M，96 trials 预算，起于 21:23，
在 `korakuen@192.168.0.183` 后台运行（`setsid`，PPID=1）。
**结果尚未评定** —— 探针要回答的是三件事，都不是"通过没有"：

1. 搜索能不能在两碎片阴离子体系上跑完
2. 实际成本（每 trial 墙钟、显存）
3. `collateral_gate` 对跨碎片探针的行为 —— `fragments_after != fragments_before` 会拒绝，
   但**两碎片合成一个**的方向从没测过（`README.md` 的未测项）


---

## 冻结后的第一次运行 · `runs/b04_frozen1` · 2026-09-13

**这是 b04 第一次在冻结判据下跑的运行。** probe1–4 都在冻结之前，是探针不是验收。

跑在 kasuga02（`korakuen@192.168.0.183`），走 `env/BASELINE.md` 的跨机协议：先核两机
md5（两个坑各踩一次 —— `ssh 不 cd`、locale 排序差异 —— 核完两机逐字节一致），
再在那边跑 digest（332 文件 / 2136 帧 / 4 事件对，**零差异**）与 pytest（两机同为
`416 passed, 1 skipped, 9 xfailed`），最后才起跑。

```bash
python -c "剥掉 // 注释" > .scratch/b04_frozen_config.json
scripts/prrs_remote.sh search benchmarks/b04_sn2_chloride/input.extxyz \
    --output runs/b04_frozen1 \
    --config .scratch/b04_frozen_config.json \
    --calculator benchmarks.b04_sn2_chloride.calculator:polar1
```

`implementation_sha256 = c393da782a1350a1`（probe4 是 `fa15a829869cdc29`）。

### 机械评定

`check.py` 退出码 0（**可裁定**）：

```
取代通道 : observed_direct   D1-D4 全 PASS
           ts0000 与 ts0001 两个 order-1 鞍点，都支持 rc0000（自环）
           event_key = 847306d3f8c57bd4  —— 与 probe2/3/4 逐字节相同
分步侧   : unanswered   理由 candidate_not_selected, proposed_not_executed
           **不是** searched_no_hit：基被截断，逐微态 uncovered=14 / not_executed=5，
           `censored_by = direction_quota_fraction`，refusals=no
```

`event_key` 与三次探针逐字节相同，是 `reaction_event_key` 改走**生成群闭包**之后
在一次活运行里的确认 —— 与那次逐条重算 50 个通道的离线结果一致。

### 与 probe4 的对照：只能读什么，不能读什么

方向排序这次由**谱**给出，运行自己记着来源
（每条 scan 的 `proposal.budget.ranked_by`）：

```
ranked_by : spectral overlap with 4 soft internal mode(s) (analytic; 0 torsional removed)
```

族内头名的变化（族间交错是既有的 round-robin 配额，没动）：

| 族 | probe4（shuffle） | frozen1（谱） | overlap |
| --- | --- | --- | --- |
| stretch | `0-1:+` | `0-1:+` | 0.9651（shuffle 撞对了） |
| compress | `2-5:-` | **`0-5:-`** | 0.2558 → **0.6893** |
| kick | `2-5:-` | **`0-1:-`** | 0.2558 → **0.9651** |
| bend | `1-0-2:-` | `1-0-2:-` | 0.0856 |

**这张表是 A 在一次运行里唯一能证明的事**：同样 28 个候选、同样选中 28、同样被
`direction_quota_fraction` 截断，**改变的只是这有限预算先花在哪根坐标上**。

两处数字确实动了，**但不得归因给 A**：

| | probe4 | frozen1 |
| --- | --- | --- |
| order-1 鞍点 | 1 | **2** |
| `reactive` 分类 | 2 | **4** |
| `directions_executed` | 17 | 18 |

**为什么不能归因**：两次运行的 `implementation_sha256` 不同，而且这一轮同时动了七个
文件（生成群闭包、`(0,π)` 守卫、谱提案、碎片刚体转动）。**这不是受控对照**，
是"同一格子在新实现下重跑一次"。要隔离 A，需要一次只切 `spectral_rank_modes`
的 A/B —— 没做，也不在这次的声明范围内。

**没有变的**：基覆盖仍是 28 选 28、逐微态 `uncovered=14`，所以**负侧照旧 `unanswered`**。
A 不改覆盖率，只改顺序；把它读成"覆盖变好了"是错的。

---

## 冻结后的第二次运行 · `runs/b04_frozen2` · 2026-09-22

**为什么又跑一次**：`runs/b04_frozen1` 的评定跑在 `c393da782a1350a1`，之后 `src/` 移了四次
（`b4e92642` → `b1c4d440` → `1a6a2971` → `ec32bc3d`）。`CRITERIA.md` 只说实现树由每次运行的
`implementation_sha256` 记录，没有任何一条规矩说指纹移位后判定作废 —— 与其论证那个缺口，
不如跑一次（六个原子，44 分钟）。

```bash
scripts/prrs_remote.sh --wait search benchmarks/b04_sn2_chloride/input.extxyz \
    --output runs/b04_frozen2 \
    --config .scratch/b04_frozen_config.json \
    --calculator benchmarks.b04_sn2_chloride.calculator:polar1
python benchmarks/b04_sn2_chloride/check.py runs/b04_frozen2/network.json   # exit 0
```

`implementation_sha256 = ec32bc3df5d2bf14`。判定与 frozen1 相同：
取代通道 `observed_direct`（D1–D4 全 PASS，`event_key = 847306d3f8c57bd4`），
分步侧 `unanswered`（`candidate_not_selected` / `proposed_not_executed`）。

### 与 frozen1 逐字段对照

| | frozen1 (`c393da78`) | frozen2 (`ec32bc3d`) |
| --- | --- | --- |
| ts_candidates | 2 `['ts0000','ts0001']` | 同 |
| outcome_counts | `{'completed': 64}` | 同 |
| classification_counts | `{'elastic': 60, 'reactive': 4}` | 同 |
| admission_counts | `{'existing_microstate': 63, 'new_microstate': 1}` | 同 |
| termination | `configured_scope_exhausted` 64/96 | 同 |
| 覆盖（逐微态） | `uncovered=14 not_executed=5` | 同 |
| budgets | — | **`+ admissions_capped_by_depth: 0`** |

**四次指纹移位对这一格的编排层输出零影响。** 唯一的差是 09-20 那件随手件新加的记账字段
（`max_depth` 这第五种预算截断有了词），值为 0 说明这一格根本没撞到深度上限。

所以「判定要不要绑指纹」那个口径不必落笔：不是"digest 零差异所以推定不变"，是**实测不变**。
