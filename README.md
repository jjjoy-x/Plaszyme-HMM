# Plaszyme HMM

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./LICENSE)
[![WebApp](https://img.shields.io/badge/WebApp-Online-brightgreen?)](http://plaszyme.org/plaszyme)
[![iGEM](https://img.shields.io/badge/iGEM-XJTLU--AI--China-blue?logo=)](https://teams.igem.org/5580)
[![iGEM GitLab](https://img.shields.io/badge/GitLab-XJTLU--AI--China-orange?logo=gitlab)](https://gitlab.igem.org/2025/software-tools/xjtlu-ai-china)

## Overview

A lightweight, reproducible pipeline for **protein functional annotation** using **Prodigal**,**Multiple Sequence Alignment (MAFFT)** and **profile HMMs (HMMER)**.

---

## Features

- **Build per-family HMMs** from protein FASTA files (≥2 sequences per family), using `hmmbuild`.
- **Annotate unknown proteins** against a combined HMM database using `hmmscan`.
- **Search downloaded databases** (e.g., UniProt subsets) **with HMMs** using `hmmsearch`.
- Ready-to-run **Slurm scripts** for HPC clusters.

---

## Directory Layout

Recommended project layout (create at any path):

```
functional_annotator/
├─ scripts/
│  └─ functional_annotator.py         # main pipeline (MAFFT + HMMER)
├─ run/
│  ├─ logs/                            # Slurm logs
│  ├─ run_annotation_dna.sh            # annotate unknown sequences with hmmscan (use DNA sequences)
│  ├─ run_searchdb_dna.sh              # search a downloaded DB with hmmsearch (use DNA sequences)
│  ├─ run_build_hmm.sh                 # build per-family HMMs and combined DB
│  ├─ run_annotation.sh                # annotate unknown sequences with hmmscan
│  └─ run_searchdb.sh                  # search a downloaded DB with hmmsearch
├─ data/                               # raw inputs (UNaligned)
│  ├─ families/                        # one FASTA per family (≥2 sequences)
│  │  ├─ plastic_families
│  │  ├─ pfam_families
│  │  └─ clan_families
│  ├─ unknown_sequences                # sequences to annotate
   |   └─ proteins_converted.fasta
│  └─ database/                              # downloaded protein DBs for searchdb
│     └─ uniprot_swissprot_2025_08.fasta
├─ hmm_profiles/                       # outputs from build (auto)
├─ annotation_results/                 # outputs from annotate (auto)
└─ search_results/                     # outputs from searchdb (auto)
```


## Dependencies

Install with conda (recommended):

```bash
conda create -n hmm python=3.10 -y
conda activate hmm

# core tools
conda install -c bioconda hmmer mafft -y

# python libs for report
conda install -c conda-forge pandas openpyxl -y
```

Verify:
```bash
mafft --version
hmmbuild -h
hmmpress -h
hmmscan -h
hmmsearch -h
```

> **Input type:** All inputs are **protein FASTA** unless you have DNA-specific HMMs (not covered here).

---

## Usage (CLI)

> Change the path to your own path

### 1) Build per-family HMMs and combined DB

```bash
cd ../HMM_functional_annotator/run
sbatch run_build_hmm.sh

```

---

### 2) Annotate unknown proteins with the HMM DB

```bash
cd ../HMM_functional_annotator/run
sbatch run_annotation.sh

```

---

### 3) Search a downloaded protein database with your HMMs 

```bash
cd ../HMM_functional_annotator/run
sbatch run_annotation.sh

```

### 4) Annotate unknown proteins with the HMM DB (use DNA sequences)

```bash
cd ../HMM_functional_annotator/run
sbatch run_annotation_dna.sh

```

---

### 5) Search a downloaded protein database with your HMMs (use DNA sequences)

```bash
cd ../HMM_functional_annotator/run
sbatch run_annotation_dna.sh

```

## Parameters & Tips

| Parameter | Meaning | Typical |
|---|---|---|
| `--threads` | CPU threads for MAFFT/HMMER | 8–32 |
| `--iE` | domain independent E-value cutoff | `1e-5` (strict) |
| `--bits` | domain bit-score cutoff | `20–30` |
| `--qcov` | query coverage threshold (0–1) | `0.30–0.50` |

- Ensure **protein** FASTA inputs for protein HMMs.
- Each family FASTA should have **≥2 sequences**; otherwise `hmmbuild` is not meaningful.

---

## License

This project is licensed under the MIT License.  
See the [LICENSE](./LICENSE) file for details.

---

## Acknowledgments

- Thanks to the developers of **HMMER**, **MAFFT**, and **Prodigal**, which form the core of this pipeline.  
- We acknowledge the use of publicly available databases such as **UniProt** and **Pfam**, which provide essential protein family and sequence data.  
- This project was developed as part of the **iGEM XJTLU-AI-China** team efforts, with contributions from both wet-lab and dry-lab members.  
- We are grateful to the open-source community for providing the computational biology tools and resources that made this work possible.
- Inspiration for this work was partly drawn from the **mparty model**, which motivated the design of the functional annotation workflow.   
- Special thanks to collaborators and advisors who provided feedback and guidance throughout the development of this pipeline.  
