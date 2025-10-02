module load hmmer/3.3.2-gcc-8.5.0-7vk4bd2
module load mafft/7.505-gcc-8.5.0-omqojgg

PROJECT_DIR="your_project_dir_path"
FAMILIES_DIR="your_families_dir_path"
HMM_DIR="your_hmm_dir_path"
ALN_DIR="your_aln_dir_path"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

cd "$PROJECT_DIR"
python functional_annotator.py build\
  --families_dir "$FAMILIES_DIR" \
  --aln_dir "$ALN_DIR" \
  --hmm_dir "$HMM_DIR" \
  --threads "$THREADS"
