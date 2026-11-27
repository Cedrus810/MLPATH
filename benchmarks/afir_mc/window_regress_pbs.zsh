#!/bin/zsh
# PBS 模板。只改标了 ← 的几处，其余不动。来源：/home/ruigengji/E3D/learn_thermal.zsh 头部（跑通过）。
# 不加 set -u（mamba activate 遇到未定义变量会退出）；不夹带冒烟/测试步骤。
# 提交：qsub job.zsh
# 节点备注：yayoi07 物理 5×RTX A4000 16GB，PBS 只登记 ngpus=4、ncpus=32、mem≈126GB（2026-10-23 pbsnodes + nvidia-smi 核实）
# #PBS 行末尾不能带注释（qsub 报 directive error，2026-10-23 实测），要改的项说明写在上一行：
# ← -q：walltime > 72h 改成 week
#PBS -q default
# ← select：host 填用户给的节点；ncpus/ngpus 填用户给的数，给多少占多少，不往下压；mem 不能填满（全占直接报错），非特殊情况留 10 GB，即用户给的内存减 10 GB（纯 CPU 节点去掉 ngpus）
#PBS -l select=1:ncpus=32:ngpus=1:mem=118gb:host=yayoi15
# ← walltime：default 最多 72h；week 最多 144h
#PBS -l walltime=72:00:00
#PBS -j oe
# ← -N：作业名
#PBS -N window_regress
test $PBS_O_WORKDIR && cd $PBS_O_WORKDIR
. /home/apps/Modules/init/profile.sh
export MODULEPATH=/home/ruigengji/modulefiles:$MODULEPATH
# module load CUDA                                         # ← 只有 pyscf-env 需要（用户印象，未核实），其他环境别加

export MAMBA_EXE=/home/ruigengji/miniforge3/bin/mamba
export MAMBA_ROOT_PREFIX=/home/ruigengji/miniforge3
source /home/ruigengji/miniforge3/etc/profile.d/mamba.sh
mamba activate openmm_dev                                  # ← 环境名

cd /home/ruigengji/MLPATH                                  # ← 项目目录
export OMP_NUM_THREADS=16 MKL_NUM_THREADS=16               # ← yayoi15 7950X3D 16 物理核；串行单进程
export AFIR_MAX_ATOMS_PER_CALL=768                         # ← 4090 24 GB 单进程（PBS 29729 的 3090 同值）
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python    # ← 和环境名一起改
echo "=== $(date '+%F %T')  host $(hostname) ==="
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

# ← 运行命令：--locator window 回归（DONE_LOG 2026-11-25；用户「可以」）。与 PBS 29729（b04/b05/N1）和 b09 第一阶段
# 同体系、同种子、同壳层、同预算、同 32 个朝向；b04/b05/N1 只多开 --locator window（其余 11-15 规则），
# b09 沿用第一阶段的 A/C/D 并把定位器换成 window（peak-dedup 在 window 模式不生效）。yayoi15 单卡、串行，最后检验。
# 可续跑：同一命令重提。
D=benchmarks/afir_mc/results/2026-11-25/window_regress
mkdir -p $D
R=benchmarks/afir_mc/run_mc_afir_batch.py
SC=(--mode sc --sc-poses 32 --wave 32 --locator window)
B9=(--shell 2.5 --budget 1700000 --sc-split 1 --skip-start-min 1 --confirm-dedup 1)
run() { local name=$1; shift; echo "=== $(date '+%F %T') start $name"
  $P -u $R "$@" --out $D/$name > $D/$name.log 2>&1 || echo "=== $name FAILED"
  echo "=== $(date '+%F %T') end $name"; }
# 串行（用户 2026-11-25「串行」）：先 b04 sA 快速确认，再 b09，之后 N1、b05、b04 sB
run b04_sA --system b04 --shell 2.0 --seed 20261023 --budget 161275 $SC
run b09_sB --system b09 --seed 20261024 $SC $B9
run n1_sA  --system n1  --shell 2.1 --seed 20261023 --budget 1676490 $SC
run b05_sA --system b05 --shell 2.4 --seed 20261023 --budget 2084515 $SC
run b04_sB --system b04 --shell 2.0 --seed 20261024 --budget 165930 $SC
$P -u benchmarks/afir_mc/afir_validate.py $D/validate.json $D/b04_sA:b04:sA $D/b04_sB:b04:sB $D/b05_sA:b05:sA $D/n1_sA:n1:sA > $D/validate.log 2>&1 || echo "=== validate FAILED"
echo "=== done $(date '+%F %T') ==="
