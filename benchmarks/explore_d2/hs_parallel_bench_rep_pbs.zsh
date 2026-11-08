#!/bin/zsh
# PBS 模板。只改标了 ← 的几处，其余不动。来源：/home/ruigengji/E3D/learn_thermal.zsh 头部（跑通过）。
# 不加 set -u（mamba activate 遇到未定义变量会退出）；不夹带冒烟/测试步骤。
# 提交：qsub job.zsh
# 节点备注：yayoi07 物理 5×RTX A4000 16GB，PBS 只登记 ngpus=4、ncpus=32、mem≈126GB（2026-10-23 pbsnodes + nvidia-smi 核实）
# #PBS 行末尾不能带注释（qsub 报 directive error，2026-10-23 实测），要改的项说明写在上一行：
# ← -q：walltime > 72h 改成 week
#PBS -q default
# ← select：host 填空节点；ncpus/ngpus/mem 按节点改（纯 CPU 节点去掉 ngpus）
#PBS -l select=1:ncpus=32:ngpus=2:mem=120gb:host=yayoi03
# ← walltime：default 最多 72h；week 最多 144h
#PBS -l walltime=24:00:00
#PBS -j oe
# ← -N：作业名
#PBS -N hs_par_rep
test $PBS_O_WORKDIR && cd $PBS_O_WORKDIR
. /home/apps/Modules/init/profile.sh
export MODULEPATH=/home/ruigengji/modulefiles:$MODULEPATH
# module load CUDA                                         # ← 只有 pyscf-env 需要（用户印象，未核实），其他环境别加

export MAMBA_EXE=/home/ruigengji/miniforge3/bin/mamba
export MAMBA_ROOT_PREFIX=/home/ruigengji/miniforge3
source /home/ruigengji/miniforge3/etc/profile.d/mamba.sh
mamba activate openmm_dev                                  # ← 环境名

cd /home/ruigengji/MLPATH                                  # ← 项目目录
export OMP_NUM_THREADS=16 MKL_NUM_THREADS=16                # ← yayoi03（EPYC 7452 32 物理核，2×RTX 3090）：两次重复各占一张卡，各 16 线程，合计 32
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python    # ← 和环境名一起改
echo "=== $(date '+%F %T')  host $(hostname) ==="
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

# ← 运行命令：并行 HS 测量重复两次（用户 11-08「BC 我全给你3090 跑」），与 11-08 yayoi43（4080S）第一次同项同臂；分清 fp32 波动与并发度对结局的影响。可续跑（每臂一个文件）。
CUDA_VISIBLE_DEVICES=0 $P -u benchmarks/explore_d2/hs_parallel_bench.py benchmarks/explore_d2/results/2026-11-08/hs_parallel_bench/full_r1 benchmarks/explore_d2/results/2026-11-04/midpoint_fire_validate/full/items.json benchmarks/explore_d2/results 0 serial,p4,p8,p16 > benchmarks/explore_d2/results/2026-11-08/hs_parallel_bench/full_r1.log 2>&1 &
CUDA_VISIBLE_DEVICES=1 $P -u benchmarks/explore_d2/hs_parallel_bench.py benchmarks/explore_d2/results/2026-11-08/hs_parallel_bench/full_r2 benchmarks/explore_d2/results/2026-11-04/midpoint_fire_validate/full/items.json benchmarks/explore_d2/results 0 serial,p4,p8,p16 > benchmarks/explore_d2/results/2026-11-08/hs_parallel_bench/full_r2.log 2>&1 &
wait
echo "=== done $(date '+%F %T') ==="
