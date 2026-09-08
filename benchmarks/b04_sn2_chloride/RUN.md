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

## 为什么工厂在这个目录里

`prrs.calculators.mace_factory` 调 `MACECalculator` 时不传 `model_type`，
所以 POLAR-1（checkpoint 类是 `PolarMACE`）会被当成普通 MACE 载入，
**算错而不报错**。在有东西真正走过这条路之前不去拓宽 `mace_factory`，
工厂就放在 `calculator.py`，经 `--calculator module:function` 这个既有钩子进来。

## 2026-09-04 那次探针

`runs/b04_probe`，MACE-POLAR-1-M，96 trials 预算，起于 21:23，
在 `korakuen@192.168.0.183` 后台运行（`setsid`，PPID=1）。
**结果尚未评定** —— 探针要回答的是三件事，都不是"通过没有"：

1. 搜索能不能在两碎片阴离子体系上跑完
2. 实际成本（每 trial 墙钟、显存）
3. `collateral_gate` 对跨碎片探针的行为 —— `fragments_after != fragments_before` 会拒绝，
   但**两碎片合成一个**的方向从没测过（`README.md` 的未测项）
