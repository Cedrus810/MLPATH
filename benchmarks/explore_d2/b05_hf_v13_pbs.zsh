#!/bin/zsh
# PBS 模板。只改标了 ← 的几处，其余不动。来源：/home/ruigengji/E3D/learn_thermal.zsh 头部（跑通过）。
# 不加 set -u（mamba activate 遇到未定义变量会退出）；不夹带冒烟/测试步骤。
# 提交：qsub job.zsh
# 节点备注：yayoi07 物理 5×RTX A4000 16GB，PBS 只登记 ngpus=4、ncpus=32、mem≈126GB（2026-10-23 pbsnodes + nvidia-smi 核实）
# #PBS 行末尾不能带注释（qsub 报 directive error，2026-10-23 实测），要改的项说明写在上一行：
# ← -q：walltime > 72h 改成 week
#PBS -q default
# ← select：host 填空节点；ncpus/ngpus/mem 按节点改（纯 CPU 节点去掉 ngpus）
#PBS -l select=1:ncpus=32:ngpus=1:mem=48gb:host=yayoi32
# ← walltime：default 最多 72h；week 最多 144h
#PBS -l walltime=36:00:00
#PBS -j oe
# ← -N：作业名
#PBS -N b05_hf_v13
test $PBS_O_WORKDIR && cd $PBS_O_WORKDIR
. /home/apps/Modules/init/profile.sh
export MODULEPATH=/home/ruigengji/modulefiles:$MODULEPATH
# module load CUDA                                         # ← 只有 pyscf-env 需要（用户印象，未核实），其他环境别加

export MAMBA_EXE=/home/ruigengji/miniforge3/bin/mamba
export MAMBA_ROOT_PREFIX=/home/ruigengji/miniforge3
source /home/ruigengji/miniforge3/etc/profile.d/mamba.sh
mamba activate openmm_dev                                  # ← 环境名

cd /home/ruigengji/MLPATH                                  # ← 项目目录
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8                  # ← yayoi32（Ryzen 9950X 16 物理核，1×RTX 4080 SUPER，64 GB）：两块共用一张卡，各 8 线程，合计 16 = 物理核数
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python    # ← 和环境名一起改
echo "=== $(date '+%F %T')  host $(hostname) ==="
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

# ← 运行命令：b05 D2-HF-v1.3（docs/PLAN_HF_D2_V13_2026-11-07.md）= v1.1 起点分簇 + 同键图的脊最后派发（v1.2）+ HS 括区规则 HF-EF-bracket-v7（两阶段中点、精阶段 FIRE、粗阶段出候选、端点出候选；粗阈值 2e-3、B 0.025）。与 10-29 yayoi32 b05_hf_v4_s{A,B} 同种子同参数，对比用。yayoi32 那次两块共用一张 4080S 实耗 5.4 h，walltime 放宽到 36 h。可续跑：在同一节点重提即接着跑（NFS 属性缓存 6 h，续跑必须在写出 run 的同一主机上）。noise_eV 来自 b05_noise_r2（2026-10-23，M5-D5 规则）。
R=benchmarks/explore_d2/results/2026-11-07
ARGS=(--shell 2.4 --n-samples 8 --n-far 4 --noise-eV 3.90625e-3 --per-pose-quota 1000
      --search-budget 340000 --certify-budget 240000 --scan-cap 800 --reverse-cap 200
      --release-cap 400 --certify-cap 1500 --audit-fraction 0.1 --max-active 4 --wave-size 16)
$P -u benchmarks/explore_d2/run_b05_d2.py --run $R/b05_hf_v13_sA --hf --start-cluster --same-graph-ridge --hs-bracket --seed 20261023 $ARGS > $R/b05_hf_v13_sA.log 2>&1 &
$P -u benchmarks/explore_d2/run_b05_d2.py --run $R/b05_hf_v13_sB --hf --start-cluster --same-graph-ridge --hs-bracket --seed 20261024 $ARGS > $R/b05_hf_v13_sB.log 2>&1 &
wait
echo "=== done $(date '+%F %T') ==="
