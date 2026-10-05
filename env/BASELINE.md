# 复跑基线

记录时间 2026-09-05 · 主机 kasuga01

`implementation_sha256` = **`99edddde34294b69…`**（2026-09-26；全量通过由用户报告，变异检验在 yayoi）
沿革：`69d53501…`（本文件初次记录）→ `5b69dfed41c4bdb2…` → `d7d53219baef17e2…` → `fa15a829869cdc29…`（b04 probe4）→ `c393da782a1350a1…`（2026-09-13）→ `b4e92642c66c9bab…`（2026-09-20 上午）→ `b1c4d44009f8785d…`（2026-09-20 下午）→ `1a6a297143deb645…`（2026-09-20 晚,resume）→ `ec32bc3df5d2bf14…`（2026-09-21 验收收尾）→ `63be7064f9e09f82…`（2026-09-22 committor 判归修复）→ `4014c4b267b25142…`（2026-09-23 committor 射击不再二次扰动）→ `ae47c4d9142203dc…`（2026-09-24 续跑回滚重做 + 三处小修）→ `b601b124e14538a2…`（2026-09-26 committor 首次进入）→ **`99edddde34294b69…`**（2026-09-26 C1 样本分类 + D1 下降几何）

**2026-09-26 第二次移位：C1 + D1。** scan 样本多 `outcome` / `outcome_reason` 两个键；鞍点下降的 `descent.extxyz` 现在两支都在（`descent.jsonl` 里 minus 支的 `frame` 编号因此从 plus 支之后接着数）。只动新 run 的输出，已有 run 与 digest 读的 `runs/p3_paths/` 不受影响。

**2026-09-26 这次移位：committor 首次进入（A5），digest 零差异。** 默认关（`committor_entry_steps = 0`），普通 trial 与固定规则射击的记录不变；committor 记录的 `p_B` 改名 `p_B_given_AB` 并多 `termination`（committor 本身默认关）。配置多两个字段 `committor_entry_steps` / `committor_max_steps`，新 run 的 `config.json` 因此多两个键。kasuga02：`make check` 绿，变异检验两条红，生产路径重打确认批 200 发逐发等价。

**2026-09-24 这次移位：续跑回滚重做 + 三处小修，digest 零差异。** 未中断运行的输出只多一个文件 `resume_point.json`（每个微态出队前写；`network.json` / `attempts.jsonl` 逐字节不变）。resume 多两条拒绝：`input_structure_sha256` 不符、无 `resume_point.json`。`completeness.cover` 拒非正交归一基底。committor 记录多三个键（默认关）。本机 kasuga01：476 passed / 1 skipped，digest 332 / 2136 / 4 零差异，变异检验五次全红。细节 `docs/REVIEW_2026-09-24.md` §2。

**2026-09-23 这次移位只动 committor 射击，digest 零差异：** `run_trial` 新增 `deliver=True` 形参，`shoot` 传 False（跳过位移 / kick / pulse，只留热动量），射击记录多一个 `probe_delivery` 键，普通 trial 记录逐字节不变。测试 +1 `test_a_shot_starts_from_the_frame_it_is_given`，**变异检验**：`shoot` 改回 `deliver=True` ⇒ 它红。kasuga02：466 passed / 7 skipped，digest 332 / 2136 / 4 零差异。

**2026-09-22 这次移位只动 committor 的判归，digest 零差异：**

| | |
| --- | --- |
| `Registry.identify()` 新增 | 只读的 `key_of()["key"] → find_node → _match_microstate`。committor 的射击判归此前读 `record["target_microstate"]`，而那个键只有扫描路径在 `admit` 之后才写 —— A/B 两个分支不可达，首次实跑四个 bracket 全报 `A 0 / B 0 / other 16` |
| `bracket["committor"]["landings"]` 新增 | 记下每发落在哪个微态；原来的 `other N` 分不出第三个盆、淬火没终止、判归有 bug |
| 测试 +1 | `test_a_shot_that_terminates_is_judged_a_basin_not_other`。**变异检验**：换回旧键 ⇒ 它红、其余 5 条绿 |

committor 默认关，所以语料里没有这条路径 —— digest 零差异是它应有的样子，不是它的证据。

**2026-09-21 这次移位是三件收尾，digest 零差异，行为改动只有一处（manifest 多一个字段）：**

| | |
| --- | --- |
| manifest 记 `host` | `search.py` 的 manifest 新增 `"host": platform.node()`。此前只有 `versions()["platform"]` 的内核串，跨主机区分是**碰巧**能用（HANDOFF §9）。它进 `experiment_sha256`（该值本来就含内核串，所以不是新性质）；**不影响续跑**，`_rebuild_for_resume` 比的是 `config` 与 `implementation_sha256` |
| `make fmt` 全量 | **15 个文件**（`src/` 5 + `tests/` 10）自 2026-09-13 起漂出格式而每道门都绿 —— 因为 `make check` 不含 `fmt-check` |
| `make check` 加 `fmt-check` | 上面那条的成因，补在 Makefile 里，`fmt-check` 排第一 |

照 `ROUNDS.md` 的规矩，全量格式化补了变异检验，靶点就用它点名的那两个 fail-closed 分支：
`internal_mode_recheck` 的 `"blocked": True → False`、`follow_unstable_mode` 的三条件判据 `→ if False`。
**两条安全性测试如期失败**（`test_a_recheck_that_could_not_run_fails_the_confirmation_closed`、
`test_a_contracting_soft_coordinate_no_longer_discards_the_descent_endpoint`），还原后 81 passed。
排版没有把断言变成空壳。

在 kasuga01 上跑（全栈，`test_backend_parity` 确实被收集）：

```
make check = fmt-check + lint + test + digest
44 files already formatted · All checks passed! · 470 passed, 1 skipped
digest：332 files, 2136 frames, 4 event pairs —— unchanged
```

**2026-09-20 下午的这次移位，digest 零差异。** 在 yayoi 上跑：

```
462 passed, 7 skipped（环境性：MACE checkpoint/CUDA 本机无）, 9 xfailed
digest：332 files, 2136 frames, 4 event pairs —— unchanged
新增 40 条测试（canonical 12 + committor 5 + completeness 7 + resume 3
  + internal 2 + soft_mode_tiers 1 + quotient 改写并摘 xfail + schema 上午已计）
```

这一轮（下午+晚,两移位）动了 `canonical.py(新) / chemistry / perturbations / network / runner /
search / internal / config / cli / completeness.py(新)`。晚上的 resume 移位：464 passed、
kill→续跑逐字节验收、digest 零差异。**关键的性质**：正则标号引擎落地、
`structural_key` 成为跨编号身份、轨道取代色类（先量：2136 帧零差异，改动免费）、
P1-3/P1-4/P1-7 关闭；`chemical_key` 的判据路径在语料上零行为改动（digest 为证）。

**2026-09-20 的这次移位，digest 零差异。** 在 yayoi 上跑：

```
332 files, 2136 frames, 4 event pairs —— digests unchanged
  + 254 comparisons on reviewed migration 2026-09-04-unknown-is-not-empty
pytest：420 passed, 7 skipped, 9 xfailed（7 个 skip 全部环境性：本机无 MACE
checkpoint / 无 CUDA；kasuga02 上同一套为 416 passed, 1 skipped 基线）
新增 10 条测试：schema 锚 9（tests/test_schema.py）+ P1-4 拒建 1
```

这一轮改了 `search / network / cli` 三个文件 + 新增 `schema.py`（schema 显式化、
`publish()` 末尾记录 `schema_warnings`、P1-4 溢出拒建），**没有一个持久化摘要移动**。

**2026-09-13 的那次移位，digest 零差异。** 在 kasuga02 上跑（`python -B`，先核过两机 md5）：

```
332 files, 2136 frames, 4 event pairs —— digests unchanged
  + 254 comparisons on reviewed migration 2026-09-04-unknown-is-not-empty
四条命名检查全绿（300 构象塌成一个 key、P1 两结构一个物质、P2 两节点两物质、P3 A≠B）
pytest 两机一致：416 passed, 1 skipped, 9 xfailed
```

这一轮改了 `chemistry / internal / config / network / perturbations / runner / search` 七个文件
（`reaction_event_key` 改走生成群闭包、`rotate_about_vertex` 的 (0,π) 守卫、谱提案、碎片刚体转动），
**没有一个持久化摘要移动** —— 与那次逐条重算 50 个 `event_key` 的结果一致。

## 环境（`versions()` 探测的九项，进 manifest）

| | |
| --- | --- |
| ase | 3.29.0 |
| numpy | 2.4.3 |
| scipy | 1.17.1 |
| torch | 2.12.0 |
| mace-torch | 0.3.16 |
| openmm | 8.5.2 |
| openmm-ml | installed; no version metadata |
| openmm-torch | installed; no version metadata |
| nnpops | installed; no version metadata |
| python | 3.12.13 |
| platform | Linux-7.2.0-1-cachyos-x86_64-with-glibc2.44 |

完整清单：`openmm_dev.yml`（conda）、`pip-freeze.txt`（pip）。

**不在 `_VERSION_PROBES` 里、因而不进 manifest 指纹的**：`geometric 1.1.1`、`opi 2.0.0`、
`rdkit`。三者都不被 `prrs` import，只在 preflight 与外部对照里用。

## 第二个环境：`.venv-tools`（只放 ruff）

**ruff 不在 `openmm_dev` 里，是有意的** —— 见 `Makefile` 头部：
backend parity 曾因 openmm 的一个 patch 版本变红，所以纯开发工具一律不许装进科学环境。

| | |
| --- | --- |
| 位置 | `.venv-tools/`（仓库内，不进结果） |
| ruff | 0.16.6 · python 3.12.13 |
| 缺失时 | `make tools` 重建（`venv` + `pip install ruff`） |

`make lint` / `fmt` / `fmt-check` 走 `.venv-tools/bin/ruff`，
`make test` / `digest` / `hash` 走 `openmm_dev`。**别在 `openmm_dev` 里 `pip install ruff`。**

## 复跑

```
make check          # ruff + pytest + digest
make hash           # 现算 implementation_sha256
```

ORCA 在 `korakuen@192.168.0.183`（`/home/korakuen/orca_6_1_1_..._openmpi418_avx2`）与
本机 `/home/ruigengji/ORCA611`。`/home/ruigengji` 两机 NFS 共享。

## 计算放在哪（用户 2026-09-09 定，当天两次调整）

**两台机器都可用**：kasuga01（本机）与 kasuga02（`korakuen@192.168.0.183`），各一块 RTX 2080 Ti。
`/home/ruigengji` 是两机同一个 NFS 挂载，所以 `runs/<name>` 两边指同一个目录，没有东西要拷。

| 地址 | 是谁 | 用法 |
| --- | --- | --- |
| `192.168.0.181` / `.182` | **本机 kasuga01** 的两块网卡（`eno1np0` / `eno2np1`） | 直接跑，不用 ssh |
| `192.168.0.183` | kasuga02 | `scripts/prrs_remote.sh` / `scripts/orca_remote.sh` |
| `192.168.0.180` | 通（0.25 ms）但**没有可用登录**（korakuen/ruigengji/kasuga/yayoi/gpu 全部 publickey 拒绝） | **暂不使用**（用户 2026-09-09） |

沿革：当天先定"全部丢 183"，随后改为"181 和 183 都可以"，180 探过之后搁置。
先前那句"本机只做编辑"已作废。

| | |
| --- | --- |
| MLP 搜索 | `scripts/prrs_remote.sh [--wait] <prrs 参数...>` |
| ORCA | `scripts/orca_remote.sh <input.inp> [nprocs]` |
| pytest / digest | `ssh korakuen@192.168.0.183 'cd /home/ruigengji/MLPATH; ...'`，见下 |

两处不随 NFS 过去，两处都先踩过：

1. **ssh 不 cd。** 不写 `cd /home/ruigengji/MLPATH` 就落在远端登录用户的家目录。
2. **`~/.cache/mace` 是每机各一份**，而 kasuga02 上登录的是 `korakuen`。
   所以要显式给 `PRRS_MACE_MODEL=/home/ruigengji/MLP/mace/MACE-OFF24_medium.model`
   —— 与本机 `~/.cache/mace/` 那份 **sha256 相同**（`e5ccf583…bc69`，已核对）。
   不给它，`test_backend_parity.py` 那 6 条会静默 skip：`231 passed, 7 skipped` 看着是绿的，
   而 backend parity 恰恰是这个项目上过一次颜色的那道网（openmm 小版本升级）。
   给了它，两机逐条一致：`237 passed, 1 skipped`。

liveness 不从本机断言：两机共享文件系统但**不共享进程表**，`scan_status.sh` 的注释
记着这一条是怎么学来的。

### 在一台机器上编辑、去另一台跑：先核 md5

2026-09-09 实测踩到：在 kasuga01 改完 `src/prrs/perturbations.py`，**立刻**在 kasuga02
上跑 pytest，得到 `KeyError: 'candidate_ids'` —— 那边还看着旧文件（NFS 客户端属性缓存，
也可能是 `__pycache__` 的 mtime 在两机之间比出了错的新旧）。同一条测试一分钟后再跑就绿。

**差一点把它读成代码缺陷。** 规则：跨机跑之前对一下

```bash
md5sum src/prrs/<改过的文件>
ssh korakuen@192.168.0.183 'cd /home/ruigengji/MLPATH; md5sum src/prrs/<同一个文件>'
```

两边不一致就等，或者 `python -B`（不写不读 pyc）再跑一次。
这条对**跑一个多小时的搜索**尤其要紧：起跑时读到旧代码，而 manifest 记的是
`make hash` 在本机现算的新哈希 —— 那就是一份自相矛盾却看不出来的记录。

