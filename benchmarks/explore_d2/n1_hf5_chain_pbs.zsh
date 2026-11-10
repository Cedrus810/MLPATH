#!/bin/zsh
# PBS 模板。只改标了 ← 的几处，其余不动。来源：/home/ruigengji/E3D/learn_thermal.zsh 头部（跑通过）。
# 不加 set -u（mamba activate 遇到未定义变量会退出）；不夹带冒烟/测试步骤。
# 提交：qsub job.zsh
# 节点备注：yayoi07 物理 5×RTX A4000 16GB，PBS 只登记 ngpus=4、ncpus=32、mem≈126GB（2026-10-23 pbsnodes + nvidia-smi 核实）
# #PBS 行末尾不能带注释（qsub 报 directive error，2026-10-23 实测），要改的项说明写在上一行：
# ← -q：walltime > 72h 改成 week
#PBS -q default
# ← select：host 填空节点；ncpus/ngpus/mem 按节点改（纯 CPU 节点去掉 ngpus）
#PBS -l select=1:ncpus=32:ngpus=1:mem=56gb:host=yayoi41
# ← walltime：default 最多 72h；week 最多 144h
#PBS -l walltime=48:00:00
#PBS -j oe
# ← -N：作业名
#PBS -N n1_hf5_chain
test $PBS_O_WORKDIR && cd $PBS_O_WORKDIR
. /home/apps/Modules/init/profile.sh
export MODULEPATH=/home/ruigengji/modulefiles:$MODULEPATH
# module load CUDA                                         # ← 只有 pyscf-env 需要（用户印象，未核实），其他环境别加

export MAMBA_EXE=/home/ruigengji/miniforge3/bin/mamba
export MAMBA_ROOT_PREFIX=/home/ruigengji/miniforge3
source /home/ruigengji/miniforge3/etc/profile.d/mamba.sh
mamba activate openmm_dev                                  # ← 环境名

cd /home/ruigengji/MLPATH                                  # ← 项目目录
export OMP_NUM_THREADS=16 MKL_NUM_THREADS=16                # ← 节点物理核数（yayoi41 / 31：16；EPYC 节点按实际改）
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python    # ← 和环境名一起改
echo "=== $(date '+%F %T')  host $(hostname) ==="
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

# ← 运行命令：HF-5 N1 一次交完（docs/PLAN_HF5_CROSS_SYSTEM_2026-11-08.md §5 n1–n2）：
#   1. noise_eV（run_d2.py --measure-noise，定义同 b05 M5-D5）；已有 noise.json 就跳过
#   2. v1.2 基线 sA，noise_eV 直接取上一步的 noise_eV_rule；L0 / cap / 额度同 b05 v1.4 整跑（N1 19 原子，b05 16 原子）
#   3. 用 2 的中点测 B_N1（midpoint_sign_resolution.py，S 组 24 点 × 8 次，再跑 F 组）
# 2 与 3 可续跑（同一节点重提即接着跑）。耗时不知道：b05 v1.3 串行整跑 5.3–9 h 可作量级参考；walltime 48 h。
set -e
R=benchmarks/explore_d2/results/2026-11-09
N=$R/n1_noise_pbs
[ -f $N/noise.json ] || $P -u benchmarks/explore_d2/run_d2.py --system n1 --measure-noise $N --seed 20261023 --shell 2.1 --n-samples 8 --n-far 4 > $R/n1_noise_pbs.log 2>&1
NOISE=$($P -c "import json; print(json.load(open('$N/noise.json'))['noise_eV_rule'])")
echo "noise_eV $NOISE"
$P -u benchmarks/explore_d2/run_d2.py --system n1 --run $R/n1_hf5_v12_sA --hf --start-cluster --same-graph-ridge \
  --seed 20261023 --shell 2.1 --n-samples 8 --n-far 4 --noise-eV $NOISE --per-pose-quota 1000 \
  --search-budget 340000 --certify-budget 240000 --scan-cap 800 --reverse-cap 200 --release-cap 400 \
  --certify-cap 1500 --audit-fraction 0.1 --max-active 8 --wave-size 16 > $R/n1_hf5_v12_sA.log 2>&1
echo "=== v1.2 done $(date '+%F %T') ==="
for G in S F; do
  $P -u benchmarks/explore_d2/midpoint_sign_resolution.py $R/n1_sign_resolution/full $G 8 8 20261109 0 $R/n1_hf5_v12_sA > $R/n1_sign_resolution_$G.log 2>&1
done
echo "=== done $(date '+%F %T') ==="
