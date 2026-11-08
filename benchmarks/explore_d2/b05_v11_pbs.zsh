#!/bin/zsh
# PBS 模板。只改标了 ← 的几处，其余不动。来源：/home/ruigengji/E3D/learn_thermal.zsh 头部（跑通过）。
# 不加 set -u（mamba activate 遇到未定义变量会退出）；不夹带冒烟/测试步骤。
# 提交：qsub job.zsh
# 节点备注：yayoi07 物理 5×RTX A4000 16GB，PBS 只登记 ngpus=4、ncpus=32、mem≈126GB（2026-10-23 pbsnodes + nvidia-smi 核实）
# #PBS 行末尾不能带注释（qsub 报 directive error，2026-10-23 实测），要改的项说明写在上一行：
# ← -q：walltime > 72h 改成 week
#PBS -q default
# ← select：host 填空节点；ncpus/ngpus/mem 按节点改（纯 CPU 节点去掉 ngpus）
#PBS -l select=1:ncpus=16:ngpus=1:mem=60gb:host=yayoi10
# ← walltime：default 最多 72h；week 最多 144h
#PBS -l walltime=48:00:00
#PBS -j oe
# ← -N：作业名
#PBS -N b05_v11
test $PBS_O_WORKDIR && cd $PBS_O_WORKDIR
. /home/apps/Modules/init/profile.sh
export MODULEPATH=/home/ruigengji/modulefiles:$MODULEPATH
# module load CUDA                                         # ← 只有 pyscf-env 需要（用户印象，未核实），其他环境别加

export MAMBA_EXE=/home/ruigengji/miniforge3/bin/mamba
export MAMBA_ROOT_PREFIX=/home/ruigengji/miniforge3
source /home/ruigengji/miniforge3/etc/profile.d/mamba.sh
mamba activate openmm_dev                                  # ← 环境名

cd /home/ruigengji/MLPATH                                  # ← 项目目录
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8                  # ← 两个块并行，各 8 线程（共 16 = ncpus）
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python    # ← 和环境名一起改
echo "=== $(date '+%F %T')  host $(hostname) ==="
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

# ← 运行命令：b05 D2-v1.1（F3 判读修复 + F6 扫描切向传入精修），全程一个 fp32 模型；与 b05_fp32_s{A,B} 同种子同参数，做 A/B。可续跑：在同一节点重提即接着跑。
#   （续跑必须在写出 run 的同一主机上，NFS 属性缓存 6 h）。noise_eV 来自 b05_noise_r2（2026-10-23，M5-D5 规则）。
R=benchmarks/explore_d2/results/2026-10-24
ARGS=(--shell 2.4 --n-samples 8 --n-far 4 --noise-eV 3.90625e-3 --per-pose-quota 1000
      --search-budget 340000 --certify-budget 240000 --scan-cap 800 --reverse-cap 200
      --release-cap 400 --certify-cap 1500 --audit-fraction 0.1 --max-active 4 --wave-size 16)
$P -u benchmarks/explore_d2/run_b05_d2.py --run $R/b05_v11_sA --seed 20261023 $ARGS > $R/b05_v11_sA.log 2>&1 &
$P -u benchmarks/explore_d2/run_b05_d2.py --run $R/b05_v11_sB --seed 20261024 $ARGS > $R/b05_v11_sB.log 2>&1 &
wait
echo "=== done $(date '+%F %T') ==="
