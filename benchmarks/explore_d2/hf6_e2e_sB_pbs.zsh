#!/bin/zsh
# PBS 模板。只改标了 ← 的几处，其余不动。来源：/home/ruigengji/E3D/learn_thermal.zsh 头部（跑通过）。
# 不加 set -u（mamba activate 遇到未定义变量会退出）；不夹带冒烟/测试步骤。
# 提交：qsub job.zsh
# 节点备注：yayoi07 物理 5×RTX A4000 16GB，PBS 只登记 ngpus=4、ncpus=32、mem≈126GB（2026-10-23 pbsnodes + nvidia-smi 核实）
# #PBS 行末尾不能带注释（qsub 报 directive error，2026-10-23 实测），要改的项说明写在上一行：
# ← -q：walltime > 72h 改成 week
#PBS -q default
# ← select：host 填空节点；ncpus/ngpus/mem 按节点改（纯 CPU 节点去掉 ngpus）
#PBS -l select=1:ncpus=32:ngpus=1:mem=64gb:host=yayoi35
# ← walltime：default 最多 72h；week 最多 144h
#PBS -l walltime=48:00:00
#PBS -j oe
# ← -N：作业名
#PBS -N hf6_e2e_sB
test $PBS_O_WORKDIR && cd $PBS_O_WORKDIR
. /home/apps/Modules/init/profile.sh
export MODULEPATH=/home/ruigengji/modulefiles:$MODULEPATH
# module load CUDA                                         # ← 只有 pyscf-env 需要（用户印象，未核实），其他环境别加

export MAMBA_EXE=/home/ruigengji/miniforge3/bin/mamba
export MAMBA_ROOT_PREFIX=/home/ruigengji/miniforge3
source /home/ruigengji/miniforge3/etc/profile.d/mamba.sh
mamba activate openmm_dev                                  # ← 环境名

cd /home/ruigengji/MLPATH                                  # ← 项目目录
export OMP_NUM_THREADS=10 MKL_NUM_THREADS=10                # ← yayoi35（EPYC 9354 32 物理核，1×RTX 5080，128 GB）：3 个进程共用一张卡，各 10 线程 = 30 ≤ 32
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python    # ← 和环境名一起改
echo "=== $(date '+%F %T')  host $(hostname) ==="
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

# ← 运行命令：HF-6 第二阶段第二个种子（20261024；用户 11-11「要不先把第二份 DB D2 做了」）。D2：b04 / b05 用已有 v1.4 sB（ledger 28071 / 405075）；N1 新跑 v1.2 sB（noise 0.001953125，其余同 sA，冒烟过的路径），跑完读其 ledger 作 N1 的 B*，接着跑 N1 DB。DB 设定同 sA（hf6_e2e_db.py，2026-11-11 冒烟通过）。
# 预计：b04 / b05 DB 4–8 h；N1 D2 约 3 h（sA 在 yayoi41 共卡实测）+ N1 DB 4–8 h；walltime 48 h。可续跑（D2 与 DB 都是同一节点重提接着跑）。
set -e
R=benchmarks/explore_d2/results
O=$R/2026-11-11/hf6_e2e_db_sB
N1B=$R/2026-11-11/n1_hf5_v12_sB_n0020
$P -u benchmarks/explore_d2/hf6_e2e_db.py $O/b04 $R/2026-11-09/b04_hf5_v14_sB 28071 20261024 > ${O}_b04.log 2>&1 &
$P -u benchmarks/explore_d2/hf6_e2e_db.py $O/b05 $R/2026-11-08/b05_hf_v14_sB 405075 20261024 > ${O}_b05.log 2>&1 &
{
  $P -u benchmarks/explore_d2/run_d2.py --system n1 --run $N1B --hf --start-cluster --same-graph-ridge \
    --seed 20261024 --shell 2.1 --n-samples 8 --n-far 4 --noise-eV 0.001953125 --per-pose-quota 1000 \
    --search-budget 340000 --certify-budget 240000 --scan-cap 800 --reverse-cap 200 --release-cap 400 \
    --certify-cap 1500 --audit-fraction 0.1 --max-active 8 --wave-size 16 > ${N1B}.log 2>&1
  B=$($P -c "import json; print(json.load(open('$N1B/report.json'))['ledger_total'])")
  echo "=== n1 D2 sB done $(date '+%F %T'), ledger $B ==="
  $P -u benchmarks/explore_d2/hf6_e2e_db.py $O/n1 $N1B $B 20261024 > ${O}_n1.log 2>&1
} &
wait
echo "=== done $(date '+%F %T') ==="
