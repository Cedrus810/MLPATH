#!/bin/zsh
# PBS 模板。只改标了 ← 的几处，其余不动。来源：/home/ruigengji/E3D/learn_thermal.zsh 头部（跑通过）。
# 不加 set -u（mamba activate 遇到未定义变量会退出）；不夹带冒烟/测试步骤。
# 提交：qsub job.zsh
# 节点备注：yayoi07 物理 5×RTX A4000 16GB，PBS 只登记 ngpus=4、ncpus=32、mem≈126GB（2026-10-23 pbsnodes + nvidia-smi 核实）
# #PBS 行末尾不能带注释（qsub 报 directive error，2026-10-23 实测），要改的项说明写在上一行：
# ← -q：walltime > 72h 改成 week
#PBS -q default
# ← select：host 填空节点；ncpus/ngpus/mem 按节点改（纯 CPU 节点去掉 ngpus）
#PBS -l select=1:ncpus=32:ngpus=1:mem=48gb:host=yayoi31
# ← walltime：default 最多 72h；week 最多 144h
#PBS -l walltime=2:00:00
#PBS -j oe
# ← -N：作业名
#PBS -N b04_signres
test $PBS_O_WORKDIR && cd $PBS_O_WORKDIR
. /home/apps/Modules/init/profile.sh
export MODULEPATH=/home/ruigengji/modulefiles:$MODULEPATH
# module load CUDA                                         # ← 只有 pyscf-env 需要（用户印象，未核实），其他环境别加

export MAMBA_EXE=/home/ruigengji/miniforge3/bin/mamba
export MAMBA_ROOT_PREFIX=/home/ruigengji/miniforge3
source /home/ruigengji/miniforge3/etc/profile.d/mamba.sh
mamba activate openmm_dev                                  # ← 环境名

cd /home/ruigengji/MLPATH                                  # ← 项目目录
export OMP_NUM_THREADS=16 MKL_NUM_THREADS=16                # ← yayoi31（Ryzen 9950X 16 物理核，1×RTX 4080 SUPER，64 GB）：一个进程，线程 = 物理核数
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python    # ← 和环境名一起改
echo "=== $(date '+%F %T')  host $(hostname) ==="
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

# ← 运行命令：HF-5 a2（docs/PLAN_HF5_CROSS_SYSTEM_2026-11-08.md §4）：b04 符号分辨界线测量，S 组（a1 的收敛中点按 |RF| 三层各 8 个 = 24 点）× 8 次重复 = 192 次松弛；a1 没有失败中点，不跑 F 组。
# 冒烟 2026-11-08 180 通过（smoke/，1 点 × 2 次，19–23 forward、1.4–3.3 s / 次）。预计 5–11 min + 模型加载；walltime 2 h。可续跑：已有结果文件跳过。
R=benchmarks/explore_d2/results/2026-11-08
$P -u benchmarks/explore_d2/midpoint_sign_resolution.py $R/b04_sign_resolution/full S 8 8 20261108 0 $R/b04_hf5_v12_sA > $R/b04_sign_resolution/full.log 2>&1
echo "=== done $(date '+%F %T') ==="
