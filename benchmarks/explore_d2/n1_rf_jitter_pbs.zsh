#!/bin/zsh
# PBS 模板。只改标了 ← 的几处，其余不动。来源：/home/ruigengji/E3D/learn_thermal.zsh 头部（跑通过）。
# 不加 set -u（mamba activate 遇到未定义变量会退出）；不夹带冒烟/测试步骤。
# 提交：qsub job.zsh
# 节点备注：yayoi07 物理 5×RTX A4000 16GB，PBS 只登记 ngpus=4、ncpus=32、mem≈126GB（2026-10-23 pbsnodes + nvidia-smi 核实）
# #PBS 行末尾不能带注释（qsub 报 directive error，2026-10-23 实测），要改的项说明写在上一行：
# ← -q：walltime > 72h 改成 week
#PBS -q default
# ← select：host 填空节点；ncpus/ngpus/mem 按节点改（纯 CPU 节点去掉 ngpus）
#PBS -l select=1:ncpus=32:ngpus=1:mem=64gb:host=yayoi38
# ← walltime：default 最多 72h；week 最多 144h
#PBS -l walltime=4:00:00
#PBS -j oe
# ← -N：作业名
#PBS -N n1_rfjitter
test $PBS_O_WORKDIR && cd $PBS_O_WORKDIR
. /home/apps/Modules/init/profile.sh
export MODULEPATH=/home/ruigengji/modulefiles:$MODULEPATH
# module load CUDA                                         # ← 只有 pyscf-env 需要（用户印象，未核实），其他环境别加

export MAMBA_EXE=/home/ruigengji/miniforge3/bin/mamba
export MAMBA_ROOT_PREFIX=/home/ruigengji/miniforge3
source /home/ruigengji/miniforge3/etc/profile.d/mamba.sh
mamba activate openmm_dev                                  # ← 环境名

cd /home/ruigengji/MLPATH                                  # ← 项目目录
export OMP_NUM_THREADS=16 MKL_NUM_THREADS=16                # ← yayoi38（EPYC 9354 32 物理核，1×RTX 5080，128 GB）：一个进程，16 线程（单进程用不满 32 核）
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python    # ← 和环境名一起改
echo "=== $(date '+%F %T')  host $(hostname) ==="
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

# ← 运行命令：HF-5 n1（docs/PLAN_HF5_CROSS_SYSTEM_2026-11-08.md §5）：N1 同支反力抖动 σ 与几何噪声，只取成键改变的脊（用户 11-09「可以」），4 条扫描 × 峰点 / 峰前一点。
# 同一脚本第一版 11-09 在 180 上 N1 全程跑通（46 min，取点规则不对，见 DONE_LOG）；本版只改了选扫描那一段。预计 ≤ 1 h（180 / 2080 Ti 实测 46 min）；walltime 4 h。不可续跑（中断重跑即从头）。
R=benchmarks/explore_d2/results/2026-11-09
$P -u benchmarks/explore_d2/rf_jitter_scan.py $R/rf_jitter_n1_v2n.json n1 20261023 2.1 8 4 4 > $R/rf_jitter_n1_v2n.log 2>&1
echo "=== done $(date '+%F %T') ==="
