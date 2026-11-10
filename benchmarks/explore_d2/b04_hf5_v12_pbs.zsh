#!/bin/zsh
# PBS 模板。只改标了 ← 的几处，其余不动。来源：/home/ruigengji/E3D/learn_thermal.zsh 头部（跑通过）。
# 不加 set -u（mamba activate 遇到未定义变量会退出）；不夹带冒烟/测试步骤。
# 提交：qsub job.zsh
# 节点备注：yayoi07 物理 5×RTX A4000 16GB，PBS 只登记 ngpus=4、ncpus=32、mem≈126GB（2026-10-23 pbsnodes + nvidia-smi 核实）
# #PBS 行末尾不能带注释（qsub 报 directive error，2026-10-23 实测），要改的项说明写在上一行：
# ← -q：walltime > 72h 改成 week
#PBS -q default
# ← select：host 填空节点；ncpus/ngpus/mem 按节点改（纯 CPU 节点去掉 ngpus）
#PBS -l select=1:ncpus=20:ngpus=1:mem=64gb:host=yayoi42
# ← walltime：default 最多 72h；week 最多 144h
#PBS -l walltime=6:00:00
#PBS -j oe
# ← -N：作业名
#PBS -N b04_hf5_v12
test $PBS_O_WORKDIR && cd $PBS_O_WORKDIR
. /home/apps/Modules/init/profile.sh
export MODULEPATH=/home/ruigengji/modulefiles:$MODULEPATH
# module load CUDA                                         # ← 只有 pyscf-env 需要（用户印象，未核实），其他环境别加

export MAMBA_EXE=/home/ruigengji/miniforge3/bin/mamba
export MAMBA_ROOT_PREFIX=/home/ruigengji/miniforge3
source /home/ruigengji/miniforge3/etc/profile.d/mamba.sh
mamba activate openmm_dev                                  # ← 环境名

cd /home/ruigengji/MLPATH                                  # ← 项目目录
export OMP_NUM_THREADS=12 MKL_NUM_THREADS=12                # ← yayoi42（i7-12700K 8P+4E = 12 物理核 / 20 线程，1×RTX 4080 SUPER，128 GB）：一个进程，线程 = 物理核数
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python    # ← 和环境名一起改
echo "=== $(date '+%F %T')  host $(hostname) ==="
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

# ← 运行命令：HF-5 a1（docs/PLAN_HF5_CROSS_SYSTEM_2026-11-08.md §4）：b04 D2-HF-v1.2 基线，种子 sA；不加括区规则（B_b04 未测）。冒烟 2026-11-08 180 通过（b04_hf5_smoke_180，与 10-28 v1.1 冒烟结局逐项相同）。
# L0 与 cap 同 b05 v1.4 整跑；壳层 2.0（b04 实测接触比 2.02）；certify-max-rise 用默认 3.0（M5-R10，冒烟为了对照用的是 1.0）。
# 额度：search 47000 = 184 意图上限 × 冒烟 127 forward / 意图 × 2；refine 200000（冒烟升级比例 6/42 推约 26 项 × 6000，取余量）。
# 预计：约 27k forward；冒烟在 180（2080 Ti）实测 0.047 s / forward（含开销），约 20–30 min + 模型加载；walltime 6 h。可续跑：同一节点重提即接着跑。
R=benchmarks/explore_d2/results/2026-11-08
$P -u benchmarks/explore_d2/run_d2.py --system b04 --run $R/b04_hf5_v12_sA --hf --start-cluster --same-graph-ridge \
  --seed 20261023 --shell 2.0 --n-samples 8 --n-far 4 --noise-eV 0.001953125 --per-pose-quota 1000 \
  --search-budget 47000 --certify-budget 200000 --scan-cap 800 --reverse-cap 200 --release-cap 400 \
  --certify-cap 1500 --audit-fraction 0.1 --max-active 8 --wave-size 16 > $R/b04_hf5_v12_sA.log 2>&1
echo "=== done $(date '+%F %T') ==="
