# 16S rRNA Amplicon Sequencing Pipeline

QIIME2-based pipeline for 16S rRNA amplicon analysis, including denoising, taxonomic classification, diversity analysis, and functional prediction.

## Repository Structure

```
16S/
├── main.sh              # Main QIIME2 pipeline (Steps 1–8)
├── submit_16S.pbs       # PBS job submission wrapper (HPC)
├── ggpicrust.R          # Functional pathway abundance analysis (PICRUSt2 output)
├── mob_analysis.R       # Methane-oxidising bacteria stacked bar plot
├── table_combine.R      # Merge feature table with taxonomic ranks
├── Data/
│   ├── raw_data/        # Input FASTQ files + manifest.tsv
│   ├── processed_data/  # Intermediate QIIME2 artifacts (.qza/.qzv)
│   ├── reference_dbs/   # Classifier and reference sequences (not tracked in git)
│   └── metadata/        # Sample metadata (metadata.tsv, decontam-metadata.tsv)
├── Results/
│   ├── denoise_mode/    # DADA2 outputs (table, rep-seqs, taxonomy, diversity)
│   └── cluster_mode/    # OTU clustering outputs
└── Logs/
    ├── checkpoints/     # Step completion checkpoints
    └── qiime2_pipeline.log
```

## Pipeline Steps

| Step | Function |
|------|----------|
| 1 | Environment setup — create/verify `qiime2` conda environment |
| 2 | Import paired-end FASTQ files via manifest |
| 3 | Visualize demultiplexed reads (quality check) |
| 4 | Primer removal with Cutadapt (optional) |
| 5 | **Denoising** (DADA2) or **OTU clustering** (vsearch) |
| 6 | Decontamination (optional, requires control samples) |
| 7 | Taxonomic classification (SILVA classifier) |
| 8 | Phylogenetic tree + alpha/beta diversity |

## Quick Start

### Local
```bash
# Run all steps
bash main.sh

# Resume from step 5
bash main.sh 5

# OTU clustering mode instead of DADA2
bash main.sh -m cluster

# Custom DADA2 truncation lengths and threads
TRUNC_LEN_F=230 TRUNC_LEN_R=200 N_THREADS=16 bash main.sh
```

### HPC (PBS)
```bash
qsub submit_16S.pbs

# Pass arguments (e.g. start from step 5)
qsub -- submit_16S.pbs 5
```

PBS resources: 16 cores, 64 GB RAM, 24 h walltime.

## Required Input

**`Data/raw_data/manifest.tsv`** — tab-separated file with columns:
```
sample-idforward-absolute-filepathreverse-absolute-filepath
```

**`Data/metadata/metadata.tsv`** — QIIME2 metadata file (required for diversity step). If running decontamination, add a `control_status` column (`sample` / `control`).

## Reference Databases

Place files in `Data/reference_dbs/`. The classifier is auto-detected in this order:

1. `silva-138-99-nb-classifier.qza`
2. `classifier.qza`
3. `silva-v4-classifier.qza`
4. `custom-classifier.qza`

Download pre-trained SILVA 138 classifier:
```bash
wget -P Data/reference_dbs/ \
  https://data.qiime2.org/2024.10/common/silva-138-99-nb-classifier.qza
```

## Key Options

| Variable | Default | Description |
|----------|---------|-------------|
| `MODE` | `denoise` | `denoise` (DADA2) or `cluster` (OTU) |
| `ENV_NAME` | `qiime2` | Conda environment name |
| `TRUNC_LEN_F` | `250` | DADA2 forward truncation length |
| `TRUNC_LEN_R` | `250` | DADA2 reverse truncation length |
| `N_THREADS` | `0` (auto) | Number of threads for DADA2 |
| `PRIMER_CHOICE` | `3` (skip) | `1`=515F/806R, `2`=custom, `3`=skip |
| `DECONTAMINATION_CHOICE` | `2` (skip) | `1`=run decontam, `2`=skip |

## Checkpoint System

Completed steps are saved as checkpoints in `Logs/checkpoints/`. The pipeline skips already-completed steps automatically.

```bash
bash main.sh -s              # Show completion status
bash main.sh -r 7            # Re-run step 7
bash main.sh -c              # Reset all checkpoints (start fresh)
bash main.sh -d              # Delete step 4/5 intermediate files
```

## Downstream R Scripts

- **`table_combine.R`** — Merges the exported feature table with taxonomy ranks into a single TSV (`final-table-with-ranks.tsv`)
- **`mob_analysis.R`** — Generates stacked bar plots of methane-oxidising bacteria (MOB) from taxonomy results
- **`ggpicrust.R`** — Extracts and visualises targeted functional pathways from PICRUSt2 output

## Dependencies

- [QIIME2](https://qiime2.org/) ≥ 2024.10 (amplicon distribution)
- Conda / Miniconda3
- R ≥ 4.0 with packages: `tidyverse`, `ggplot2`, `dplyr`, `readr`

Install QIIME2 environment:
```bash
conda env create -n qiime2 \
  --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2025.10/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml
```
