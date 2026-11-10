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
#PBS -l walltime=4:00:00
#PBS -j oe
# ← -N：作业名
#PBS -N n1_rfjit_seeds
test $PBS_O_WORKDIR && cd $PBS_O_WORKDIR
. /home/apps/Modules/init/profile.sh
export MODULEPATH=/home/ruigengji/modulefiles:$MODULEPATH
# module load CUDA                                         # ← 只有 pyscf-env 需要（用户印象，未核实），其他环境别加

export MAMBA_EXE=/home/ruigengji/miniforge3/bin/mamba
export MAMBA_ROOT_PREFIX=/home/ruigengji/miniforge3
source /home/ruigengji/miniforge3/etc/profile.d/mamba.sh
mamba activate openmm_dev                                  # ← 环境名

cd /home/ruigengji/MLPATH                                  # ← 项目目录
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4                  # ← yayoi41（Ryzen 9950X3D 16 物理核，1×RTX 5080 16 GB，64 GB）：每批 4 个进程共用一张卡，各 4 线程 = 16 物理核
P=/home/ruigengji/miniforge3/envs/openmm_dev/bin/python    # ← 和环境名一起改
echo "=== $(date '+%F %T')  host $(hostname) ==="
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

# ← 运行命令：HF-5 n1 续（DONE_LOG 11-09）：N1 第一个姿态上成键改变的脊为 0。用已跑通的同一脚本（rf_jitter_scan.py，yayoi38 作业 25432 跑通）换 7 个种子（各自 L0 的第一个姿态朝向不同）找成键改变的脊并测 σ / 几何噪声。代码未改。
# 两批各 4 / 3 个进程（11-09 五进程共用 16 GB 卡出过一次 OOM，这次少放）。单个运行 yayoi38（5080）实测 12 min；预计 30–60 min；walltime 4 h。不可续跑（已出 json 的种子可不重跑，脚本本身不跳过）。
R=benchmarks/explore_d2/results/2026-11-09
run() { $P -u benchmarks/explore_d2/rf_jitter_scan.py $R/rf_jitter_n1_s$1.json n1 $1 2.1 8 4 4 > $R/rf_jitter_n1_s$1.log 2>&1 }
for s in 20261024 20261025 20261026 20261027; do run $s & done
wait
for s in 20261028 20261029 20261030; do run $s & done
wait
echo "=== done $(date '+%F %T') ==="
