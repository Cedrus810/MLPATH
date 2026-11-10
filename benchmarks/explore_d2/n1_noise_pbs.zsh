#!/bin/zsh
# PBS 模板。只改标了 ← 的几处，其余不动。来源：/home/ruigengji/E3D/learn_thermal.zsh 头部（跑通过）。
# 不加 set -u（mamba activate 遇到未定义变量会退出）；不夹带冒烟/测试步骤。
# 提交：qsub job.zsh
# 节点备注：yayoi07 物理 5×RTX A4000 16GB，PBS 只登记 ngpus=4、ncpus=32、mem≈126GB（2026-10-23 pbsnodes + nvidia-smi 核实）
# #PBS 行末尾不能带注释（qsub 报 directive error，2026-10-23 实测），要改的项说明写在上一行：
# ← -q：walltime > 72h 改成 week
#PBS -q default
# ← select：host 填空节点；ncpus/ngpus/mem 按节点改（纯 CPU 节点去掉 ngpus）
#PBS -l select=1:ncpus=32:ngpus=1:mem=56gb:host=yayoi41
# ← walltime：default 最多 72h；week 最多 144h
#PBS -l walltime=2:00:00
#PBS -j oe
# ← -N：作业名
#PBS -N n1_noise
test $PBS_O_WORKDIR && cd $PBS_O_WORKDIR
. /home/apps/Modules/init/profile.sh
export MODULEPATH=/home/ruigengji/modulefiles:$MODULEPATH
# module load CUDA                                         # ← 只有 pyscf-env 需要（用户印象，未核实），其他环境别加

export MAMBA_EXE=/home/ruigengji/miniforge3/bin/mamba
export MAMBA_ROOT_PREFIX=/home/ruigengji/miniforge3
source /home/ruigengji/miniforge3/etc/profile.d/mamba.sh
mamba activate openmm_dev                                  # ← 环境名

cd /home/ruigengji/MLPATH                                  # ← 项目目录
export OMP_NUM_THREADS=16 MKL_NUM_THREADS=16                # ← yayoi41（16 物理核）：一个进程
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python    # ← 和环境名一起改
echo "=== $(date '+%F %T')  host $(hostname) ==="
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

# ← 运行命令：HF-5 n1：N1 扫描判读噪声 noise_eV（run_d2.py --measure-noise，定义同 run_b05_d2.py / M5-D5）。未冒烟：180 上那次卡在 NFS I/O（D 状态 10 min，未开始计算），用户 11-09 要求直接交节点。
R=benchmarks/explore_d2/results/2026-11-09
$P -u benchmarks/explore_d2/run_d2.py --system n1 --measure-noise $R/n1_noise_pbs --seed 20261023 --shell 2.1 --n-samples 8 --n-far 4 > $R/n1_noise_pbs.log 2>&1
echo "=== done $(date '+%F %T') ==="
