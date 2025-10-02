module load hmmer/3.3.2-gcc-8.5.0-7vk4bd2
module load mafft/7.505-gcc-8.5.0-omqojgg
module load prodigal/2.6.3-gcc-8.5.0-pziudqi

PROJECT_DIR="your_project_dir_path"
UNKNOWN_FASTA="your_unknown_fasta_path"
HMM_DIR="your_hmm_dir_path"
OUT_DIR="your_out_dir_path"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

cd "$PROJECT_DIR"
python functional_annotator.py annotate \
    --unknown_fasta "$UNKNOWN_FASTA" \
    --is_dna \
    --prodigal_mode meta \
    --genetic_code 11 \
    --hmm_dir "$HMM_DIR" \
    --out_dir "$OUT_DIR" \
    --threads "$THREADS" \
    --iE 1e-5 --bits 25 --qcov 0.30
