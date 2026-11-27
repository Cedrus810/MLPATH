#!/bin/zsh
# PBS 模板。只改标了 ← 的几处，其余不动。来源：/home/ruigengji/E3D/learn_thermal.zsh 头部（跑通过）。
# 不加 set -u（mamba activate 遇到未定义变量会退出）；不夹带冒烟/测试步骤。
# 提交：qsub job.zsh
# 节点备注：yayoi07 物理 5×RTX A4000 16GB，PBS 只登记 ngpus=4、ncpus=32、mem≈126GB（2026-10-23 pbsnodes + nvidia-smi 核实）
# #PBS 行末尾不能带注释（qsub 报 directive error，2026-10-23 实测），要改的项说明写在上一行：
# ← -q：walltime > 72h 改成 week
#PBS -q default
# ← select：host 填空节点；ncpus/ngpus/mem 按节点改（纯 CPU 节点去掉 ngpus）
#PBS -l select=1:ncpus=32:ngpus=2:mem=120gb:host=yayoi08
# ← walltime：default 最多 72h；week 最多 144h
#PBS -l walltime=48:00:00
#PBS -j oe
# ← -N：作业名
#PBS -N afir_pair2
test $PBS_O_WORKDIR && cd $PBS_O_WORKDIR
. /home/apps/Modules/init/profile.sh
export MODULEPATH=/home/ruigengji/modulefiles:$MODULEPATH
# module load CUDA                                         # ← 只有 pyscf-env 需要（用户印象，未核实），其他环境别加

export MAMBA_EXE=/home/ruigengji/miniforge3/bin/mamba
export MAMBA_ROOT_PREFIX=/home/ruigengji/miniforge3
source /home/ruigengji/miniforge3/etc/profile.d/mamba.sh
mamba activate openmm_dev                                  # ← 环境名

cd /home/ruigengji/MLPATH                                  # ← 项目目录
export OMP_NUM_THREADS=16 MKL_NUM_THREADS=16               # ← yayoi08（EPYC 7452 32 核，2 x RTX 3090 24 GB）：两路各占一张卡、各 16 线程 = 物理核数
export AFIR_MAX_ATOMS_PER_CALL=768                         # ← 每张 24 GB 卡一个进程（fulda 32 GB 实测约 2 000 原子 OOM，按显存比例再留余量）
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python    # ← 和环境名一起改
echo "=== $(date '+%F %T')  host $(hostname) ==="
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

# ← 运行命令：Non-Hessian AFIR SC，2026-11-15 定稿的默认（设计稿 §2.2 / §2.3；DONE_LOG 11-15）：
# 单对推力 alpha r_kl + 双对层 alpha (r_kl + r_mn)/sqrt 2（单独记账），TS 坐标与推力解耦（峰帧冻结权重的原子对集合），
# 不按两个极大停止，括区从峰向两侧逐轮展开，gamma 106.9；反应物 MIN = 输入 + 32 个随机朝向；预算 5 x B*
# （D2-HF-v1.5 同体系同种子 ledger；用户 11-15：保留 5xB* 上限，N1 双对按固定预算部分覆盖）。
# 两路并行，各占一张 GPU（取 PBS 给的 CUDA_VISIBLE_DEVICES，未设时 0,1）：A 路 N1 sA -> b05 sA -> b04 sA；
# B 路 N1 sB -> b05 sB -> b04 sB。最后统一检验。
# 可续跑：同一节点重提即接着跑（已完成的波次跳过）。
D=benchmarks/afir_mc/results/2026-11-15/pbs_pair_double
mkdir -p $D
R=benchmarks/afir_mc/run_mc_afir_batch.py
# 2026-11-26：默认定位器改为 window 后补 --locator hf，保证本脚本续跑 / 复现仍用原 HF 定位器（不改任何结果）。
SC=(--mode sc --sc-poses 32 --wave 32 --locator hf)
run() { local name=$1; shift; echo "=== $(date '+%F %T') start $name"
  $P -u $R "$@" --out $D/$name > $D/$name.log 2>&1 || echo "=== $name FAILED"
  echo "=== $(date '+%F %T') end $name"; }
GPUS=(${(s:,:)${CUDA_VISIBLE_DEVICES:-0,1}})
echo "GPUs for the two streams: $GPUS[1] / $GPUS[2]"
(
  export CUDA_VISIBLE_DEVICES=$GPUS[1]
  run n1_sA  --system n1  --shell 2.1 --seed 20261023 --budget 1676490 $SC
  run b05_sA --system b05 --shell 2.4 --seed 20261023 --budget 2084515 $SC
  run b04_sA --system b04 --shell 2.0 --seed 20261023 --budget 161275  $SC
) &
(
  export CUDA_VISIBLE_DEVICES=$GPUS[2]
  run n1_sB  --system n1  --shell 2.1 --seed 20261024 --budget 1761340 $SC
  run b05_sB --system b05 --shell 2.4 --seed 20261024 --budget 2075845 $SC
  run b04_sB --system b04 --shell 2.0 --seed 20261024 --budget 165930  $SC
) &
wait
$P -u benchmarks/afir_mc/afir_validate.py $D/validate.json \
  $D/b04_sA:b04:sA $D/b04_sB:b04:sB $D/b05_sA:b05:sA $D/b05_sB:b05:sB \
  $D/n1_sA:n1:sA $D/n1_sB:n1:sB > $D/validate.log 2>&1 || echo "=== validate FAILED"
echo "=== done $(date '+%F %T') ==="
