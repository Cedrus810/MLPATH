#!/bin/zsh
# PBS 模板。只改标了 ← 的几处，其余不动。来源：/home/ruigengji/E3D/learn_thermal.zsh 头部（跑通过）。
# 不加 set -u（mamba activate 遇到未定义变量会退出）；不夹带冒烟/测试步骤。
# 提交：qsub job.zsh
# 节点备注：yayoi07 物理 5×RTX A4000 16GB，PBS 只登记 ngpus=4、ncpus=32、mem≈126GB（2026-10-23 pbsnodes + nvidia-smi 核实）
# #PBS 行末尾不能带注释（qsub 报 directive error，2026-10-23 实测），要改的项说明写在上一行：
# ← -q：walltime > 72h 改成 week
#PBS -q default
# ← select：host 填用户给的节点；ncpus/ngpus 填用户给的数，给多少占多少，不往下压；mem 不能填满（全占直接报错），非特殊情况留 10 GB，即用户给的内存减 10 GB（纯 CPU 节点去掉 ngpus）
#PBS -l select=1:ncpus=32:ngpus=2:mem=118gb:host=yayoi10
# ← walltime：default 最多 72h；week 最多 144h
#PBS -l walltime=72:00:00
#PBS -j oe
# ← -N：作业名
#PBS -N w1_regress_sB
test $PBS_O_WORKDIR && cd $PBS_O_WORKDIR
. /home/apps/Modules/init/profile.sh
export MODULEPATH=/home/ruigengji/modulefiles:$MODULEPATH
# module load CUDA                                         # ← 只有 pyscf-env 需要（用户印象，未核实），其他环境别加

export MAMBA_EXE=/home/ruigengji/miniforge3/bin/mamba
export MAMBA_ROOT_PREFIX=/home/ruigengji/miniforge3
source /home/ruigengji/miniforge3/etc/profile.d/mamba.sh
mamba activate openmm_dev                                  # ← 环境名

cd /home/ruigengji/MLPATH                                  # ← 项目目录
export OMP_NUM_THREADS=16 MKL_NUM_THREADS=16               # ← yayoi10 EPYC 7452 32 物理核；两路并行各 16，合计不超过物理核数
export AFIR_MAX_ATOMS_PER_CALL=768                         # ← RTX 3090 24 GB，每路独占一张卡（PBS 29729 同值）
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python    # ← 和环境名一起改
echo "=== $(date '+%F %T')  host $(hostname) ==="
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

# ← 运行命令：AFIR-W1 开发集回归的 sB 种子（yayoi31 上的 w1_regress 跑 sA + b09 两种子）。
# 种子、壳层、预算同 PBS 29729 sB。两张卡各跑一路，每路内部串行；可续跑：同一命令重提。
D=benchmarks/afir_mc/results/2026-11-27/w1_regress_sB
mkdir -p $D
R=benchmarks/afir_mc/run_mc_afir_batch.py
W1=(--mode sc --sc-poses 64 --wave 32 --sc-split 1 --skip-start-min 1 --confirm-dedup 1 \
    --double-xx first --double-order class --reflow 32)
run() { local gpu=$1 name=$2; shift 2; echo "=== $(date '+%F %T') start $name (gpu $gpu)"
  CUDA_VISIBLE_DEVICES=$gpu $P -u $R "$@" --out $D/$name > $D/$name.log 2>&1 || echo "=== $name FAILED"
  echo "=== $(date '+%F %T') end $name"; }
( run 0 n1_sB  --system n1  --shell 2.1 --seed 20261024 --budget 1761340 $W1
  run 0 b04_sB --system b04 --shell 2.0 --seed 20261024 --budget 165930  $W1 ) &
( run 1 b05_sB --system b05 --shell 2.4 --seed 20261024 --budget 2075845 $W1 ) &
wait
$P -u benchmarks/afir_mc/afir_validate.py $D/validate.json $D/b04_sB:b04:sB $D/b05_sB:b05:sB $D/n1_sB:n1:sB > $D/validate.log 2>&1 || echo "=== validate FAILED"
echo "=== done $(date '+%F %T') ==="
