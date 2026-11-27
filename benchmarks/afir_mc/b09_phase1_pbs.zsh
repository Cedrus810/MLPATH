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
#PBS -N b09_afir1
test $PBS_O_WORKDIR && cd $PBS_O_WORKDIR
. /home/apps/Modules/init/profile.sh
export MODULEPATH=/home/ruigengji/modulefiles:$MODULEPATH
# module load CUDA                                         # ← 只有 pyscf-env 需要（用户印象，未核实），其他环境别加

export MAMBA_EXE=/home/ruigengji/miniforge3/bin/mamba
export MAMBA_ROOT_PREFIX=/home/ruigengji/miniforge3
source /home/ruigengji/miniforge3/etc/profile.d/mamba.sh
mamba activate openmm_dev                                  # ← 环境名

cd /home/ruigengji/MLPATH                                  # ← 项目目录
export OMP_NUM_THREADS=16 MKL_NUM_THREADS=16               # ← 两路并行各占 ncpus 的一半（ncpus 改了这里跟着改成 ncpus/2）
export AFIR_MAX_ATOMS_PER_CALL=768                         # ← 24 GB 卡（PBS 29729 同值）；卡更小按显存比例压
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python    # ← 和环境名一起改
echo "=== $(date '+%F %T')  host $(hostname) ==="
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

# ← 运行命令：b09 第一次测量（docs/B09_DISTANCE_2026-11-20.md §5 C）。基础规则同 PBS 29729（11-15 定稿：
# 单对推力 + 双对层、TS 坐标解耦、gamma 106.9、输入 + 32 个随机朝向），另加下面的试行规则；都对所有体系同一套，不为 b09 单独改参数。
# 壳层 2.5 = 无偏输入复合物实测接触比（build.json，同 b05 的取法）。预算每种子 1 700 000 forward（N1 同量级；b09 无 B*）。
# 两路并行各占一张 GPU（未设 CUDA_VISIBLE_DEVICES 时 0,1）。可续跑：同一命令重提，已完成波次跳过。
D=benchmarks/afir_mc/results/2026-11-22/b09_phase1
mkdir -p $D
R=benchmarks/afir_mc/run_mc_afir_batch.py
# 2026-11-22 加试行规则（用户 11-20「可以动手试试」、11-22「可以 做」，方案 PLAN_AFIR_PAIR_SPLIT §7.2 A/C/D、
# PLAN_AFIR_PEAK_DEDUP d* = 0.25 Å）：片段内原子对只在输入复合物推、跳过未变键图的起点极小化、确认前按事件判重、括区前跳过离已知鞍点 ≤ 0.25 Å 的峰。
# 2026-11-26：默认定位器改为 window 后补 --locator hf，保证本脚本续跑 / 复现仍用原 HF 定位器（不改任何结果）。
SC=(--mode sc --sc-poses 32 --wave 32 --shell 2.5 --budget 1700000 --locator hf \
    --sc-split 1 --skip-start-min 1 --confirm-dedup 1 --peak-dedup 0.25)
run() { local name=$1; shift; echo "=== $(date '+%F %T') start $name"
  $P -u $R "$@" --out $D/$name > $D/$name.log 2>&1 || echo "=== $name FAILED"
  echo "=== $(date '+%F %T') end $name"; }
GPUS=(${(s:,:)${CUDA_VISIBLE_DEVICES:-0,1}})
echo "GPUs for the two streams: ${GPUS[1]} / ${GPUS[2]}"
( export CUDA_VISIBLE_DEVICES=${GPUS[1]}; run b09_sA --system b09 --seed 20261023 $SC ) &
( export CUDA_VISIBLE_DEVICES=${GPUS[2]}; run b09_sB --system b09 --seed 20261024 $SC ) &
wait
echo "=== done $(date '+%F %T') ==="
