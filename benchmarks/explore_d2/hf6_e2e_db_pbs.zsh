#!/bin/zsh
# PBS 模板。只改标了 ← 的几处，其余不动。来源：/home/ruigengji/E3D/learn_thermal.zsh 头部（跑通过）。
# 不加 set -u（mamba activate 遇到未定义变量会退出）；不夹带冒烟/测试步骤。
# 提交：qsub job.zsh
# 节点备注：yayoi07 物理 5×RTX A4000 16GB，PBS 只登记 ngpus=4、ncpus=32、mem≈126GB（2026-10-23 pbsnodes + nvidia-smi 核实）
# #PBS 行末尾不能带注释（qsub 报 directive error，2026-10-23 实测），要改的项说明写在上一行：
# ← -q：walltime > 72h 改成 week
#PBS -q default
# ← select：host 填空节点；ncpus/ngpus/mem 按节点改（纯 CPU 节点去掉 ngpus）
#PBS -l select=1:ncpus=32:ngpus=1:mem=64gb:host=yayoi39
# ← walltime：default 最多 72h；week 最多 144h
#PBS -l walltime=48:00:00
#PBS -j oe
# ← -N：作业名
#PBS -N hf6_e2e_db
test $PBS_O_WORKDIR && cd $PBS_O_WORKDIR
. /home/apps/Modules/init/profile.sh
export MODULEPATH=/home/ruigengji/modulefiles:$MODULEPATH
# module load CUDA                                         # ← 只有 pyscf-env 需要（用户印象，未核实），其他环境别加

export MAMBA_EXE=/home/ruigengji/miniforge3/bin/mamba
export MAMBA_ROOT_PREFIX=/home/ruigengji/miniforge3
source /home/ruigengji/miniforge3/etc/profile.d/mamba.sh
mamba activate openmm_dev                                  # ← 环境名

cd /home/ruigengji/MLPATH                                  # ← 项目目录
export OMP_NUM_THREADS=10 MKL_NUM_THREADS=10                # ← yayoi39（EPYC 9354 32 物理核，1×RTX 5080，128 GB）：3 个进程共用一张卡，各 10 线程 = 30 ≤ 32
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python    # ← 和环境名一起改
echo "=== $(date '+%F %T')  host $(hostname) ==="
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

# ← 运行命令：HF-6 第二阶段 DB 臂（docs/PLAN_HF6_DIMER_E2E_2026-11-11.md §2.2）：ASE Dimer 从同一 L0 姿态原样出发，默认随机位移与随机初始模式，单次 3000 forward，累计到 B* = 对应 D2 运行的 ledger（b04 26129 / b05 408888 / N1 376296）。种子 20261111。
# 冒烟：2026-11-11 180 b04 2 次尝试跑通（均用满 3000、约 0.06 s / forward）。预计 b05 / N1 各约 4–8 h（第一阶段 3 进程共卡约 0.07 s / forward × 约 40 万）；walltime 48 h。可续跑（每次尝试一个 json，已有的计入预算并跳过；同一节点重提即接着跑）。
R=benchmarks/explore_d2/results
O=$R/2026-11-11/hf6_e2e_db
$P -u benchmarks/explore_d2/hf6_e2e_db.py $O/b04 $R/2026-11-09/b04_hf5_v14_sA 26129 20261111 > ${O}_b04.log 2>&1 &
$P -u benchmarks/explore_d2/hf6_e2e_db.py $O/b05 $R/2026-11-08/b05_hf_v14_sA 408888 20261111 > ${O}_b05.log 2>&1 &
$P -u benchmarks/explore_d2/hf6_e2e_db.py $O/n1 $R/2026-11-09/n1_hf5_v12_sA_n0020 376296 20261111 > ${O}_n1.log 2>&1 &
wait
echo "=== done $(date '+%F %T') ==="
