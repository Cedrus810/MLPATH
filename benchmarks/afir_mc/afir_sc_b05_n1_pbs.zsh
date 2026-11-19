#!/bin/zsh
# PBS 模板。只改标了 ← 的几处，其余不动。来源：/home/ruigengji/E3D/learn_thermal.zsh 头部（跑通过）。
# 不加 set -u（mamba activate 遇到未定义变量会退出）；不夹带冒烟/测试步骤。
# 提交：qsub job.zsh
# 节点备注：yayoi07 物理 5×RTX A4000 16GB，PBS 只登记 ngpus=4、ncpus=32、mem≈126GB（2026-10-23 pbsnodes + nvidia-smi 核实）
# #PBS 行末尾不能带注释（qsub 报 directive error，2026-10-23 实测），要改的项说明写在上一行：
# ← -q：walltime > 72h 改成 week
#PBS -q default
# ← select：host 填空节点；ncpus/ngpus/mem 按节点改（纯 CPU 节点去掉 ngpus）
#PBS -l select=1:ncpus=20:ngpus=1:mem=100gb:host=yayoi42
# ← walltime：default 最多 72h；week 最多 144h
#PBS -l walltime=24:00:00
#PBS -j oe
# ← -N：作业名
#PBS -N afir_sc
test $PBS_O_WORKDIR && cd $PBS_O_WORKDIR
. /home/apps/Modules/init/profile.sh
export MODULEPATH=/home/ruigengji/modulefiles:$MODULEPATH
# module load CUDA                                         # ← 只有 pyscf-env 需要（用户印象，未核实），其他环境别加

export MAMBA_EXE=/home/ruigengji/miniforge3/bin/mamba
export MAMBA_ROOT_PREFIX=/home/ruigengji/miniforge3
source /home/ruigengji/miniforge3/etc/profile.d/mamba.sh
mamba activate openmm_dev                                  # ← 环境名

cd /home/ruigengji/MLPATH                                  # ← 项目目录
export OMP_NUM_THREADS=12 MKL_NUM_THREADS=12               # ← yayoi42（i7-12700K 12 物理核 / 20 线程，RTX 4080 SUPER 16 GB）：一个进程，12 线程 = 物理核数
export AFIR_MAX_ATOMS_PER_CALL=512                         # ← 16 GB 卡：fulda 5090（32 GB）实测约 2 000 原子 OOM，按显存减半
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python    # ← 和环境名一起改
echo "=== $(date '+%F %T')  host $(hostname) ==="
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

# ← 运行命令：Non-Hessian AFIR，SC 模式（docs/NONHESSIAN_AFIR_DESIGN_2026-11-14.md、DONE_LOG 11-14 / 11-15），
# b05 / N1 各两个种子依次跑，早弃关闭（默认），反应物 MIN = 输入 + 32 个随机朝向，预算 5 x B*（D2-HF-v1.5 同体系同种子
# ledger）；最后统一检验（阶数 + 参照匹配 + 不同 TS / 事件 + 与 D2 v1.5 的 R(b)）。可续跑：同一节点重提即接着跑。
D=benchmarks/afir_mc/results/2026-11-15/pbs
mkdir -p $D
R=benchmarks/afir_mc/run_mc_afir_batch.py
SC=(--mode sc --sc-poses 32 --wave 32)
run() { local name=$1; shift; echo "=== $(date '+%F %T') start $name"
  $P -u $R "$@" --out $D/$name > $D/$name.log 2>&1 || echo "=== $name FAILED"
  echo "=== $(date '+%F %T') end $name"; }
run b05_sc_sA --system b05 --shell 2.4 --seed 20261023 --budget 2084515 $SC
run b05_sc_sB --system b05 --shell 2.4 --seed 20261024 --budget 2075845 $SC
run n1_sc_sA  --system n1  --shell 2.1 --seed 20261023 --budget 1676490 $SC
run n1_sc_sB  --system n1  --shell 2.1 --seed 20261024 --budget 1761340 $SC
$P -u benchmarks/afir_mc/afir_validate.py $D/validate.json $D/b05_sc_sA:b05:sA $D/b05_sc_sB:b05:sB \
  $D/n1_sc_sA:n1:sA $D/n1_sc_sB:n1:sB > $D/validate.log 2>&1 || echo "=== validate FAILED"
echo "=== done $(date '+%F %T') ==="
