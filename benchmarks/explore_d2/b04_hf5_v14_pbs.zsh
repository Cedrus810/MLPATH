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
#PBS -l walltime=8:00:00
#PBS -j oe
# ← -N：作业名
#PBS -N b04_hf5_v14
test $PBS_O_WORKDIR && cd $PBS_O_WORKDIR
. /home/apps/Modules/init/profile.sh
export MODULEPATH=/home/ruigengji/modulefiles:$MODULEPATH
# module load CUDA                                         # ← 只有 pyscf-env 需要（用户印象，未核实），其他环境别加

export MAMBA_EXE=/home/ruigengji/miniforge3/bin/mamba
export MAMBA_ROOT_PREFIX=/home/ruigengji/miniforge3
source /home/ruigengji/miniforge3/etc/profile.d/mamba.sh
mamba activate openmm_dev                                  # ← 环境名

cd /home/ruigengji/MLPATH                                  # ← 项目目录
export OMP_NUM_THREADS=3 MKL_NUM_THREADS=3                  # ← yayoi31（Ryzen 9950X 16 物理核，1×RTX 4080 SUPER，64 GB）：5 个进程共用一张卡，各 3 线程，合计 15 ≤ 16 物理核；换节点时按物理核数 / 5 改
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python    # ← 和环境名一起改
echo "=== $(date '+%F %T')  host $(hostname) ==="
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

# ← 运行命令：HF-5 a3 + a4（docs/PLAN_HF5_CROSS_SYSTEM_2026-11-08.md §4）：b04 D2-HF-v1.4（B_b04 = 0.0237，用户 11-09）。
# a3：种子 sA / sB；a4：sA 同命令重复（fp32 重复噪声）、sA 加 --rotate 1、sA 加 --magnitude-step 0.05（扫描点数加倍，search 额度相应加倍）。
# L0、cap、额度同 a1（v1.2 sA：search 22.8k / 47k、refine 4.9k / 200k，4080 SUPER 上约 10 min）。--max-active 8 同 b05 v1.4。
# 预计：每个运行约 10–20 min（a1 实测；5 个共用一张卡会更慢，未测），walltime 8 h。可续跑：同一节点重提即接着跑。
R=benchmarks/explore_d2/results/2026-11-09
V14=(--hf --start-cluster --same-graph-ridge --hs-bracket --sample-parallel)
L0=(--shell 2.0 --n-samples 8 --n-far 4 --noise-eV 0.001953125 --per-pose-quota 1000)
CAP=(--certify-budget 200000 --scan-cap 800 --reverse-cap 200 --release-cap 400 --certify-cap 1500
     --audit-fraction 0.1 --max-active 8 --wave-size 16)
D2=benchmarks/explore_d2/run_d2.py
$P -u $D2 --system b04 --run $R/b04_hf5_v14_sA $V14 --seed 20261023 $L0 --search-budget 47000 $CAP > $R/b04_hf5_v14_sA.log 2>&1 &
$P -u $D2 --system b04 --run $R/b04_hf5_v14_sB $V14 --seed 20261024 $L0 --search-budget 47000 $CAP > $R/b04_hf5_v14_sB.log 2>&1 &
$P -u $D2 --system b04 --run $R/b04_hf5_v14_sA_rep $V14 --seed 20261023 $L0 --search-budget 47000 $CAP > $R/b04_hf5_v14_sA_rep.log 2>&1 &
$P -u $D2 --system b04 --run $R/b04_hf5_v14_sA_rot1 $V14 --rotate 1 --seed 20261023 $L0 --search-budget 47000 $CAP > $R/b04_hf5_v14_sA_rot1.log 2>&1 &
$P -u $D2 --system b04 --run $R/b04_hf5_v14_sA_step005 $V14 --magnitude-step 0.05 --seed 20261023 $L0 --search-budget 94000 $CAP > $R/b04_hf5_v14_sA_step005.log 2>&1 &
wait
echo "=== done $(date '+%F %T') ==="
