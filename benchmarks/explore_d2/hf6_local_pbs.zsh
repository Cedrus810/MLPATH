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
#PBS -l walltime=48:00:00
#PBS -j oe
# ← -N：作业名
#PBS -N hf6_local
test $PBS_O_WORKDIR && cd $PBS_O_WORKDIR
. /home/apps/Modules/init/profile.sh
export MODULEPATH=/home/ruigengji/modulefiles:$MODULEPATH
# module load CUDA                                         # ← 只有 pyscf-env 需要（用户印象，未核实），其他环境别加

export MAMBA_EXE=/home/ruigengji/miniforge3/bin/mamba
export MAMBA_ROOT_PREFIX=/home/ruigengji/miniforge3
source /home/ruigengji/miniforge3/etc/profile.d/mamba.sh
mamba activate openmm_dev                                  # ← 环境名

cd /home/ruigengji/MLPATH                                  # ← 项目目录
export OMP_NUM_THREADS=5 MKL_NUM_THREADS=5                  # ← yayoi31（Ryzen 9950X 16 物理核，1×RTX 4080 SUPER，64 GB）：3 个进程共用一张卡，各 5 线程 = 15；换节点按物理核数 / 3 改
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python    # ← 和环境名一起改
echo "=== $(date '+%F %T')  host $(hostname) ==="
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

# ← 运行命令：HF-6 第一阶段（docs/PLAN_HF6_DIMER_2026-11-10.md §3）：三个体系的全部 HS 冻结括区，两臂 hf / dimer，单项上限 = 各运行 hs_cap（3000）。
# 冒烟：2026-11-10 180 b04 第 1 项两臂跑通（hf 242 forward 出 SN2 候选；dimer 1353 forward 收敛但漂离起点 12.3 Å、低 0.17 eV）。可续跑（每项每臂一个 json，已有跳过；同一节点重提即接着跑）。耗时见交 qsub 时给的计划表；walltime 48 h。
R=benchmarks/explore_d2/results
O=$R/2026-11-10/hf6_local
$P -u benchmarks/explore_d2/hf6_local.py $O/b04 $R/2026-11-08/b04_hf5_v12_sA 0 hf,dimer > ${O}_b04.log 2>&1 &
$P -u benchmarks/explore_d2/hf6_local.py $O/b05 $R/2026-11-08/b05_hf_v14_sA 0 hf,dimer > ${O}_b05.log 2>&1 &
$P -u benchmarks/explore_d2/hf6_local.py $O/n1 $R/2026-11-09/n1_hf5_v12_sA_n0020 0 hf,dimer > ${O}_n1.log 2>&1 &
wait
echo "=== done $(date '+%F %T') ==="
