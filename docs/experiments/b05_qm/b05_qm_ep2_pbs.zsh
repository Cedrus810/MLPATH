#!/bin/zsh
# PBS 模板。只改标了 ← 的几处，其余不动。来源：/home/ruigengji/E3D/learn_thermal.zsh 头部（跑通过）。
# 不加 set -u（mamba activate 遇到未定义变量会退出）；不夹带冒烟/测试步骤。
# 提交：qsub job.zsh
# 节点备注：yayoi07 物理 5×RTX A4000 16GB，PBS 只登记 ngpus=4、ncpus=32、mem≈126GB（2026-10-23 pbsnodes + nvidia-smi 核实）
# #PBS 行末尾不能带注释（qsub 报 directive error，2026-10-23 实测），要改的项说明写在上一行：
# ← -q：walltime > 72h 改成 week
#PBS -q default
# ← select：host 填空节点；ncpus/ngpus/mem 按节点改（纯 CPU 节点去掉 ngpus）
#PBS -l select=1:ncpus=48:mpiprocs=48:mem=150gb:host=yayoi13
# ← walltime：default 最多 72h；week 最多 144h
#PBS -l walltime=24:00:00
#PBS -j oe
# ← -N：作业名
#PBS -N b05_qm_ep2
test $PBS_O_WORKDIR && cd $PBS_O_WORKDIR
. /home/apps/Modules/init/profile.sh
export MODULEPATH=/home/ruigengji/modulefiles:$MODULEPATH
module load orca/6.1.1                                     # ← yayoi13 纯 CPU（Xeon Gold 5220R，48 核，192 GB）：ORCA 6.1.1 模块

export MAMBA_EXE=/home/ruigengji/miniforge3/bin/mamba
export MAMBA_ROOT_PREFIX=/home/ruigengji/miniforge3
source /home/ruigengji/miniforge3/etc/profile.d/mamba.sh
mamba activate openmm_dev                                  # ← 环境名

cd /home/ruigengji/MLPATH/docs/experiments/b05_qm        # ← 项目目录（QM 链目录）
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1                  # ← ORCA 用 MPI 并行；三个端点并行各 16 进程 = 48 核
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python    # ← 和环境名一起改
echo "=== $(date '+%F %T')  host $(hostname) ==="
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

command -v mpirun >/dev/null || module load openmpi5.0.10  # ← 运行命令：模块未带 mpirun 时补 Open MPI 5（183 上 ORCA 6.1.1 + Open MPI 5 实测可用）
echo "orca: $(command -v orca)  mpirun: $(command -v mpirun)"
./run_endpoint2.sh new_10-12_react 16 > endpoint2_new_10-12_react.log 2>&1 &
./run_endpoint2.sh new_10-12_prod 16 > endpoint2_new_10-12_prod.log 2>&1 &
./run_endpoint2.sh ctrl_6-9_prod 16 > endpoint2_ctrl_6-9_prod.log 2>&1 &
wait
echo "=== done $(date '+%F %T') ==="
