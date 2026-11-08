#!/bin/zsh
# PBS 模板。只改标了 ← 的几处，其余不动。来源：/home/ruigengji/E3D/learn_thermal.zsh 头部（跑通过）。
# 不加 set -u（mamba activate 遇到未定义变量会退出）；不夹带冒烟/测试步骤。
# 提交：qsub job.zsh
# 节点备注：yayoi07 物理 5×RTX A4000 16GB，PBS 只登记 ngpus=4、ncpus=32、mem≈126GB（2026-10-23 pbsnodes + nvidia-smi 核实）
# #PBS 行末尾不能带注释（qsub 报 directive error，2026-10-23 实测），要改的项说明写在上一行：
# ← -q：walltime > 72h 改成 week
#PBS -q default
# ← select：host 填空节点；ncpus/ngpus/mem 按节点改（纯 CPU 节点去掉 ngpus）
#PBS -l select=1:ncpus=32:ngpus=1:mem=56gb:host=yayoi40
# ← walltime：default 最多 72h；week 最多 144h
#PBS -l walltime=24:00:00
#PBS -j oe
# ← -N：作业名
#PBS -N b05_v6val
test $PBS_O_WORKDIR && cd $PBS_O_WORKDIR
. /home/apps/Modules/init/profile.sh
export MODULEPATH=/home/ruigengji/modulefiles:$MODULEPATH
# module load CUDA                                         # ← 只有 pyscf-env 需要（用户印象，未核实），其他环境别加

export MAMBA_EXE=/home/ruigengji/miniforge3/bin/mamba
export MAMBA_ROOT_PREFIX=/home/ruigengji/miniforge3
source /home/ruigengji/miniforge3/etc/profile.d/mamba.sh
mamba activate openmm_dev                                  # ← 环境名

cd /home/ruigengji/MLPATH                                  # ← 项目目录
export OMP_NUM_THREADS=32 MKL_NUM_THREADS=32               # ← 等于 ncpus
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python    # ← 和环境名一起改
echo "=== $(date '+%F %T')  host $(hostname) ==="
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

$P -u benchmarks/explore_d2/midpoint_coarse_candidate_validate.py benchmarks/explore_d2/results/2026-11-07/midpoint_coarse_candidate_validate/full 0.025 22 20261102 0 2 benchmarks/explore_d2/results/2026-10-26/b05_hf_sA benchmarks/explore_d2/results/2026-10-26/b05_hf_sB benchmarks/explore_d2/results/2026-10-28/b05_hf_v11_sA benchmarks/explore_d2/results/2026-10-28/b05_hf_v11_sB benchmarks/explore_d2/results/2026-10-28/b05_hf_v11r_sA benchmarks/explore_d2/results/2026-10-28/b05_hf_v11r_sB benchmarks/explore_d2/results/2026-10-29/b05_hf_v4_sA benchmarks/explore_d2/results/2026-10-29/b05_hf_v4_sB   # ← 运行命令
echo "=== done $(date '+%F %T') ==="
