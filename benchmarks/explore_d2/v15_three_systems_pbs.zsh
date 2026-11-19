#!/bin/zsh
# PBS 模板。只改标了 ← 的几处，其余不动。来源：/home/ruigengji/E3D/learn_thermal.zsh 头部（跑通过）。
# 不加 set -u（mamba activate 遇到未定义变量会退出）；不夹带冒烟/测试步骤。
# 提交：qsub job.zsh
# 节点备注：yayoi07 物理 5×RTX A4000 16GB，PBS 只登记 ngpus=4、ncpus=32、mem≈126GB（2026-10-23 pbsnodes + nvidia-smi 核实）
# #PBS 行末尾不能带注释（qsub 报 directive error，2026-10-23 实测），要改的项说明写在上一行：
# ← -q：walltime > 72h 改成 week
#PBS -q default
# ← select：host 填空节点；ncpus/ngpus/mem 按节点改（纯 CPU 节点去掉 ngpus）
#PBS -l select=1:ncpus=32:ngpus=1:mem=60gb:host=yayoi29
# ← walltime：default 最多 72h；week 最多 144h
#PBS -l walltime=48:00:00
#PBS -j oe
# ← -N：作业名
#PBS -N d2_v15_x3
test $PBS_O_WORKDIR && cd $PBS_O_WORKDIR
. /home/apps/Modules/init/profile.sh
export MODULEPATH=/home/ruigengji/modulefiles:$MODULEPATH
# module load CUDA                                         # ← 只有 pyscf-env 需要（用户印象，未核实），其他环境别加

export MAMBA_EXE=/home/ruigengji/miniforge3/bin/mamba
export MAMBA_ROOT_PREFIX=/home/ruigengji/miniforge3
source /home/ruigengji/miniforge3/etc/profile.d/mamba.sh
mamba activate openmm_dev                                  # ← 环境名

cd /home/ruigengji/MLPATH                                  # ← 项目目录
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4                  # ← yayoi29（Ryzen 5950X 16 物理核，1×RTX 3090 24 GB，64 GB）：同时 4 个进程 × 4 线程 = 16 = 物理核数；同卡 4 个约 6 GB / 进程（11-09 16 GB 卡 3 个约 5.3 GB 正常、5 个约 3.2 GB OOM）
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python    # ← 和环境名一起改
echo "=== $(date '+%F %T')  host $(hostname) ==="
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

# ← 运行命令：D2-HF-v1.5（docs/PLAN_HF_D2_V15_2026-11-13.md；用户 11-13 定 G3 / M2 / N1）在 b04 / b05 / N1 上统一重跑，
# 各两个种子（sA 20261023、sB 20261024）。σ、B、noise_eV、几何噪声全部按 v1.5 规则（--v15），不再逐体系给。
# L0 与额度同上一轮各体系运行（b04 v1.4、b05 v1.4、N1 v1.2），只换规则。单卡：b05 / N1 各两种子同时 4 个进程，b04 排在同种子 b05 之后。可续跑：同一节点重提即接着跑（NFS 属性缓存 6 h）。
R=benchmarks/explore_d2/results/2026-11-13
COMMON=(--v15 --n-samples 8 --n-far 4 --per-pose-quota 1000 --scan-cap 800 --reverse-cap 200
        --release-cap 400 --certify-cap 1500 --audit-fraction 0.1 --max-active 8 --wave-size 16)
B04=(--system b04 --shell 2.0 --search-budget 47000 --certify-budget 200000)
B05=(--system b05 --shell 2.4 --search-budget 340000 --certify-budget 240000)
N1=(--system n1 --shell 2.1 --search-budget 340000 --certify-budget 240000)
for s in 20261023 20261024; do
  t=$([[ $s == 20261023 ]] && echo sA || echo sB)
  ( $P -u benchmarks/explore_d2/run_d2.py --run $R/b05_v15_$t $B05 $COMMON --seed $s > $R/b05_v15_$t.log 2>&1
    $P -u benchmarks/explore_d2/run_d2.py --run $R/b04_v15_$t $B04 $COMMON --seed $s > $R/b04_v15_$t.log 2>&1 ) &
  $P -u benchmarks/explore_d2/run_d2.py --run $R/n1_v15_$t $N1 $COMMON --seed $s > $R/n1_v15_$t.log 2>&1 &
done
wait
echo "=== done $(date '+%F %T') ==="
