# 复跑基线

记录时间 2026-09-05 · 主机 kasuga01

`implementation_sha256` = **`5b69dfed41c4bdb2…`**（现算：`make hash`；`69d53501…` 是本文件初次记录时的值）

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

