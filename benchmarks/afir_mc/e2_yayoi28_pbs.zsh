#!/bin/zsh
# PBS 模板。只改标了 ← 的几处，其余不动。来源：/home/ruigengji/E3D/learn_thermal.zsh 头部（跑通过）。
# 不加 set -u（mamba activate 遇到未定义变量会退出）；不夹带冒烟/测试步骤。
# 提交：qsub job.zsh
# E2（docs/TS_LOCATOR_DESIGN_2026-11-19.md §5）：b04 / b05 两种子，PBS 29729 的 84 条括区路径，短 String + 固定方向爬坡。
# 第 2 版（作业 32178 失败后）：四块依次跑（并行时四个进程挤一张 3090 OOM）；zsh 里 $k:$sys 被当成修饰符，变量全加花括号。
#PBS -q default
#PBS -l select=1:ncpus=32:ngpus=1:mem=54gb:host=yayoi28
#PBS -l walltime=72:00:00
#PBS -j oe
#PBS -N e2_b04b05
test $PBS_O_WORKDIR && cd $PBS_O_WORKDIR
. /home/apps/Modules/init/profile.sh
export MODULEPATH=/home/ruigengji/modulefiles:$MODULEPATH

export MAMBA_EXE=/home/ruigengji/miniforge3/bin/mamba
export MAMBA_ROOT_PREFIX=/home/ruigengji/miniforge3
source /home/ruigengji/miniforge3/etc/profile.d/mamba.sh
mamba activate openmm_dev

cd /home/ruigengji/MLPATH
export OMP_NUM_THREADS=32 MKL_NUM_THREADS=32
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python
echo "=== $(date '+%F %T')  host $(hostname) ==="
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

R=benchmarks/afir_mc/results/2026-11-15/pbs_pair_double
D=benchmarks/afir_mc/results/2026-11-19/e2_pbs
mkdir -p $D
for run in b04_sA b04_sB b05_sA b05_sB; do
  sys=${run%_*}
  $P -u benchmarks/afir_mc/e2_string_climb.py --run ${R}/${run} --system ${sys} --out ${D}/${run} > ${D}/${run}.log 2>&1 || echo "${run} FAILED"
done
for k in 0 2; do
  for sys in b04 b05; do
    $P -u benchmarks/afir_mc/afir_validate.py ${D}/validate_${sys}_k${k}.json "${D}/${sys}_sA/k${k}:${sys}:sA" "${D}/${sys}_sB/k${k}:${sys}:sB" > ${D}/validate_${sys}_k${k}.log 2>&1 || echo "validate ${sys} k${k} FAILED"
  done
done
echo "=== done $(date '+%F %T') ==="
