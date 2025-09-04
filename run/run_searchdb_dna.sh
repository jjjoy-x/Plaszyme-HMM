#!/bin/bash
#SBATCH --job-name=HMM  # 作业名
#SBATCH --partition=gpu3090         # gpu3090队列
#SBATCH --qos=1gpu                  # gpu3090Qos
#SBATCH --nodes=1                   # 节点数量
#SBATCH --ntasks-per-node=1         # 每节点进程数
#SBATCH --cpus-per-task=4           # 1:4 的 GPU:CPU 配比 
#SBATCH --gres=gpu:1                # 2 块 GPU
#SBATCH --output=logs/%j.out        # 标准输出
#SBATCH --error=logs/%j.err         # 错误输出
#SBATCH --mail-user=yueqing.xing22@student.xjtlu.edu.cn
#SBATCH --mail-type=ALL

module load hmmer/3.3.2-gcc-8.5.0-7vk4bd2
module load mafft/7.505-gcc-8.5.0-omqojgg
module load prodigal/2.6.3-gcc-8.5.0-pziudqi

PROJECT_DIR="/gpfs/work/bio/yueqingxing22/igem/HMM_functional_annotator/scripts"
SEQ_DB="/gpfs/work/bio/yueqingxing22/igem/HMM_functional_annotator/data/DNA/random_dna_db.fasta"
HMM_DIR="/gpfs/work/bio/yueqingxing22/igem/HMM_functional_annotator/hmm_profiles/plastic_families/hmm"                           
OUT_DIR="/gpfs/work/bio/yueqingxing22/igem/HMM_functional_annotator/DNA_results/searchdb_results/plastic_families"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

cd "$PROJECT_DIR"
python functional_annotator.py searchdb \
  --seq_db "$SEQ_DB" \
  --seq_db_is_dna \
  --prodigal_mode meta \
  --genetic_code 11 \
  --hmm_dir "$HMM_DIR" \
  --out_dir "$OUT_DIR" \
  --threads "$THREADS" \
  --iE 1e-6 --bits 30 --qcov 0.35
