#!/bin/zsh
# PBS 模板。只改标了 ← 的几处，其余不动。来源：/home/ruigengji/E3D/learn_thermal.zsh 头部（跑通过）。
# 不加 set -u（mamba activate 遇到未定义变量会退出）；不夹带冒烟/测试步骤。
# 提交：qsub job.zsh
# 节点备注：yayoi07 物理 5×RTX A4000 16GB，PBS 只登记 ngpus=4、ncpus=32、mem≈126GB（2026-10-23 pbsnodes + nvidia-smi 核实）
# #PBS 行末尾不能带注释（qsub 报 directive error，2026-10-23 实测），要改的项说明写在上一行：
# ← -q：walltime > 72h 改成 week
#PBS -q default
# ← select：host 填用户给的节点；ncpus/ngpus 填用户给的数，给多少占多少，不往下压；mem 不能填满（全占直接报错），非特殊情况留 10 GB，即用户给的内存减 10 GB（纯 CPU 节点去掉 ngpus）
#PBS -l select=1:ncpus=32:ngpus=1:mem=54gb:host=yayoi32
# ← walltime：default 最多 72h；week 最多 144h
#PBS -l walltime=24:00:00
#PBS -j oe
# ← -N：作业名
#PBS -N b05_kscan
test $PBS_O_WORKDIR && cd $PBS_O_WORKDIR
. /home/apps/Modules/init/profile.sh
export MODULEPATH=/home/ruigengji/modulefiles:$MODULEPATH
# module load CUDA                                         # ← 只有 pyscf-env 需要（用户印象，未核实），其他环境别加

export MAMBA_EXE=/home/ruigengji/miniforge3/bin/mamba
export MAMBA_ROOT_PREFIX=/home/ruigengji/miniforge3
source /home/ruigengji/miniforge3/etc/profile.d/mamba.sh
mamba activate openmm_dev                                  # ← 环境名

cd /home/ruigengji/MLPATH                                  # ← 项目目录
export OMP_NUM_THREADS=16 MKL_NUM_THREADS=16               # ← yayoi32 是 Ryzen 9950X，PBS 登记 32 = 16 物理核 × 2 线程；按全局规则线程数不超过物理核数
export AFIR_MAX_ATOMS_PER_CALL=512                         # ← RTX 4080 SUPER 16 GB：按 24 GB → 768 的比例取 512；183 的 11 GB 卡未设此值、默认 1024 时 OOM（DONE_LOG 11-20）
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python    # ← 和环境名一起改
echo "=== $(date '+%F %T')  host $(hostname) ==="
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

# ← 运行命令：b05 String 步数诊断重跑（原 183 作业 PID 24862 于 2026-11-20 15:25 CUDA OOM 退出；用户 2026-11-23
# 「给你一个pbs 节点 慢慢跑」）。命令与 k183_e2_b05_kscan.sh 相同：b05 sB 中 HF 内核匹配、E2 k=2 不匹配、起点相同的
# 19 条路径，String 快照 k = 0,2,10,50,200,1000，各自爬坡，再逐 k 检验。重放帧从 11-20 目录复制（不重放）。
R=benchmarks/afir_mc/results/2026-11-15/pbs_pair_double
D=benchmarks/afir_mc/results/2026-11-23/e2_b05_kscan
mkdir -p $D/b05_sB/frames
cp -n benchmarks/afir_mc/results/2026-11-20/e2_b05_kscan/b05_sB/frames/*.npz $D/b05_sB/frames/
L="n1|2-12+13-14,n1|3-14,n2|2-12+13-14,n2|3-14,n4|2-7+11-14,n5|1-7+3-14,n5|3-10+5-14,n5|3-14,n5|3-14+6-11,n5|3-14+7-8,n6|2-8+7-14,n6|5-14,n7|5-14,n8|0-11+13-14,n8|1-11+13-14,n8|1-2+11-14,n8|11-14,n8|2-12+11-14,n8|3-10+11-14"
$P -u benchmarks/afir_mc/e2_string_climb.py --run $R/b05_sB --system b05 --out $D/b05_sB \
  --labels "$L" --checkpoints 0,2,10,50,200,1000 > $D/b05_sB.log 2>&1 || echo "=== e2 FAILED"
for k in 0 2 10 50 200 1000; do
  $P -u benchmarks/afir_mc/afir_validate.py $D/validate_k${k}.json "$D/b05_sB/k${k}:b05:sB" > $D/validate_k${k}.log 2>&1 || echo "=== validate k${k} FAILED"
done
echo "=== done $(date '+%F %T') ==="
