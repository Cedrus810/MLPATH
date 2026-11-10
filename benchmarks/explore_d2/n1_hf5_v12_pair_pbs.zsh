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
#PBS -N n1_v12_pair
test $PBS_O_WORKDIR && cd $PBS_O_WORKDIR
. /home/apps/Modules/init/profile.sh
export MODULEPATH=/home/ruigengji/modulefiles:$MODULEPATH
# module load CUDA                                         # ← 只有 pyscf-env 需要（用户印象，未核实），其他环境别加

export MAMBA_EXE=/home/ruigengji/miniforge3/bin/mamba
export MAMBA_ROOT_PREFIX=/home/ruigengji/miniforge3
source /home/ruigengji/miniforge3/etc/profile.d/mamba.sh
mamba activate openmm_dev                                  # ← 环境名

cd /home/ruigengji/MLPATH                                  # ← 项目目录
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8                  # ← yayoi41（Ryzen 9950X3D 16 物理核，1×RTX 5080，64 GB）：两个运行共用一张卡，各 8 线程 = 16
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python    # ← 和环境名一起改
echo "=== $(date '+%F %T')  host $(hostname) ==="
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

# ← 运行命令：HF-5 n2（docs/PLAN_HF5_CROSS_SYSTEM_2026-11-08.md §5）：N1 v1.2 基线 sA，noise_eV 两种取法并跑（用户未定，DONE_LOG 11-09）：
#   n0527：0.052734375（run_d2 --measure-noise 原定义，含构象盆能量差）；n0020：0.001953125（2 个 fp32 量子，b05 同为 2 个量子）。
#   各自跑完后用其中点测符号分辨界线（S、F 组）。冒烟：180 上 noise + v1.2 小规模 2026-11-09 18:04 通过。
#   L0 / cap / 额度同 b05 v1.4（冒烟实测约 238 search forward / 意图）。耗时不知道（b05 v1.3 串行 5.3–9 h 作量级）；walltime 48 h。可续跑：同一节点重提即接着跑。
R=benchmarks/explore_d2/results/2026-11-09
arm() {
  $P -u benchmarks/explore_d2/run_d2.py --system n1 --run $R/n1_hf5_v12_sA_$1 --hf --start-cluster --same-graph-ridge \
    --seed 20261023 --shell 2.1 --n-samples 8 --n-far 4 --noise-eV $2 --per-pose-quota 1000 \
    --search-budget 340000 --certify-budget 240000 --scan-cap 800 --reverse-cap 200 --release-cap 400 \
    --certify-cap 1500 --audit-fraction 0.1 --max-active 8 --wave-size 16 > $R/n1_hf5_v12_sA_$1.log 2>&1
  echo "=== $1 v1.2 done $(date '+%F %T') ==="
  for G in S F; do
    $P -u benchmarks/explore_d2/midpoint_sign_resolution.py $R/n1_sign_resolution_$1 $G 8 8 20261109 0 $R/n1_hf5_v12_sA_$1 > $R/n1_sign_resolution_$1_$G.log 2>&1
  done
  echo "=== $1 done $(date '+%F %T') ==="
}
arm n0527 0.052734375 &
arm n0020 0.001953125 &
wait
echo "=== done $(date '+%F %T') ==="
