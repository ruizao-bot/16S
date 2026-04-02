# 16S rRNA Amplicon Sequencing Pipeline

This repository contains a checkpointed QIIME2 Bash workflow for paired-end 16S amplicon data. The pipeline in `main.sh` covers environment setup, FASTQ import, demux QC, optional primer trimming, DADA2 denoising or OTU clustering, optional decontamination, taxonomic classification, and diversity analysis.

## Repository Structure

```text
16S/
├── main.sh                   # Main pipeline (steps 1–8)
├── README.md                 # Usage notes
├── Data/
│   ├── raw_data/             # FASTQ files + manifest.tsv
│   ├── processed_data/       # Imported / trimmed QIIME2 artifacts
│   ├── metadata/             # metadata.tsv and optional decontam metadata
│   └── reference_dbs/        # Classifiers or reference reads/taxonomy
├── Results/
│   ├── denoise_mode/         # DADA2 outputs
│   └── cluster_mode/         # vsearch clustering outputs
└── Logs/
    ├── checkpoints/          # Per-step checkpoint files
    └── qiime2_pipeline.log   # Pipeline log
```

## Pipeline Steps

| Step | Action | Notes |
|------|--------|-------|
| 1 | Environment setup | Verifies Conda and creates/uses `ENV_NAME` (default `qiime2`) |
| 2 | Import data | Imports paired-end FASTQ from `Data/raw_data/manifest.tsv` |
| 3 | Visualize demux | Produces `demux-paired-end.qzv` for quality review |
| 4 | Remove primers | Optional; `PRIMER_CHOICE=1` uses 515F/806R, `2` uses custom primers, default skips |
| 5 | Denoise or cluster | `MODE=denoise` runs DADA2; `MODE=cluster` runs the vsearch OTU workflow |
| 6 | Decontamination | Enabled by default; requires metadata with a `control_status` column |
| 7 | Taxonomic classification | Uses an existing classifier or builds one with RESCRIPt |
| 8 | Tree + diversity | Builds a phylogeny and runs Shannon / unweighted UniFrac analyses |

## Requirements

- Conda / Miniconda available on `PATH`
- QIIME2 amplicon environment compatible with the pipeline
- Paired-end manifest at `Data/raw_data/manifest.tsv`
- Sample metadata at `Data/metadata/metadata.tsv` for grouped diversity visualizations
- Classifier at `16S/Data/reference_dbs`, default is `silva-138-99-nb-classifier.qza`

Step 1 uses the QIIME2 amplicon distribution below if the environment does not already exist:

```bash
mamba env create -f qiime_env.yml
```

## Required Input Files

### `Data/raw_data/manifest.tsv`
Expected tab-separated columns:

```text
sample-id	forward-absolute-filepath	reverse-absolute-filepath
```

### `Data/metadata/metadata.tsv`
Used in step 8 for group significance and Emperor plots.

### Optional decontam metadata
If `DECONTAMINATION_CHOICE=1`, provide either:

- `Data/metadata/decontam-metadata.tsv`, or
- `Data/metadata/metadata.tsv`

with a `control_status` column containing values such as `sample` and `control`.

## Local Usage

```bash
bash main.sh [OPTIONS] [START_STEP]
```

### Common examples

```bash
# Run the full pipeline
bash main.sh

# Resume from step 5
bash main.sh 5

# Run clustering mode instead of DADA2
bash main.sh -m cluster

# Custom DADA2 truncation lengths and threads
TRUNC_LEN_F=230 TRUNC_LEN_R=200 N_THREADS=16 bash main.sh

# Use the built-in 515F/806R primers for trimming
PRIMER_CHOICE=1 bash main.sh 4

# Rebuild a classifier in the current environment and rerun taxonomy
bash main.sh -r 7
CLASSIFIER_CHOICE=2 bash main.sh 7
```

### Command-line options

| Option | Meaning |
|--------|---------|
| `-h`, `--help` | Show help |
| `-c`, `--clean` | Remove all checkpoints and the pipeline log |
| `-s`, `--status` | Show step completion status |
| `-d`, `--delete-intermediate` | Delete selected step 4/5 intermediate outputs |
| `-r`, `--remove <N|name>` | Remove a specific checkpoint |
| `-m`, `--mode <denoise|cluster>` | Select pipeline mode |
| `-e`, `--env <name>` | Set the Conda environment name |
| `START_STEP` | Resume from a step number `1`–`8` |

## SLURM Usage

An example job wrapper is provided in `run_qiime.slurm`:

```bash
sbatch run_qiime.slurm
```

Current defaults in that wrapper:

- partition: `large_336`
- CPUs: `8`
- memory: `64G`
- walltime: `24:00:00`
- starts the pipeline from step `2`

## Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `ENV_NAME` | `qiime2` | Conda environment used by the pipeline |
| `MODE` | `denoise` | `denoise` for DADA2 or `cluster` for vsearch OTU mode |
| `TRUNC_LEN_F` | `250` | Forward truncation length for DADA2 |
| `TRUNC_LEN_R` | `250` | Reverse truncation length for DADA2 |
| `N_THREADS` | `0` | Threads for DADA2 (`0` lets QIIME decide) |
| `PRIMER_CHOICE` | `3` | `1` = built-in 515F/806R, `2` = custom primers, other = skip |
| `PRIMER_F` / `PRIMER_R` | unset | Custom primer sequences when `PRIMER_CHOICE=2` |
| `CLUSTER_CHOICE` | `1` | Cluster mode: `1` de novo, `2` closed reference, `3` skip extra clustering |
| `DECONTAMINATION_CHOICE` | `1` | `1` runs decontamination by default; set another value to skip |
| `DECONTAM_METADATA` | unset | Optional explicit metadata path for decontam |
| `DECONTAM_THRESHOLD` | `0.5` | Threshold for contaminant calling |
| `REMOVE_CONTROLS` | `y` | Remove negative controls after decontam |
| `CLASSIFIER_CHOICE` | `skip` | `1` = use `USER_CLASSIFIER`, `2` = build with RESCRIPt, otherwise skip if none found |
| `USER_CLASSIFIER` | unset | Path to a classifier when `CLASSIFIER_CHOICE=1` |
| `EXTRACT_PRIMERS` | `n` | When building a classifier, optionally extract reads matching `PRIMER_F` / `PRIMER_R` |

## Reference Databases and Classifiers

Step 7 auto-detects the first existing classifier in this order:

1. `Data/reference_dbs/silva-138-99-nb-classifier.qza`
2. `Data/reference_dbs/classifier.qza`
3. `Data/reference_dbs/silva-v4-classifier.qza`
4. `Data/reference_dbs/custom-classifier.qza`

If none is found, the script can:

- use a user-supplied classifier with `CLASSIFIER_CHOICE=1 USER_CLASSIFIER=/path/to/file.qza`, or
- build a new classifier with `CLASSIFIER_CHOICE=2` using RESCRIPt.

> `silva-138-99-seqs-515-806.qza` is a reference-reads artifact, not a trained classifier. Step 7 needs a compatible `*-classifier.qza` or a rebuild via `CLASSIFIER_CHOICE=2`.

> Pretrained classifiers must match the current QIIME2 / `scikit-learn` version. If you hit a version-mismatch error, rebuild the classifier in the active environment.

## Checkpoints and Re-running Steps

Each step writes a checkpoint to `Logs/checkpoints/`. On later runs, completed steps are skipped automatically.

```bash
bash main.sh -s      # show status
bash main.sh -r 7    # remove step 7 checkpoint
bash main.sh -c      # clear all checkpoints and log
bash main.sh -d      # remove selected intermediate outputs
```

If you change parameters, metadata, or classifiers, remove the relevant checkpoint before rerunning that step.

## Outputs

Results are written under `Results/<mode>_mode/`, including:

- `table.qza` / `rep-seqs.qza`
- `taxonomy.qza`, `taxonomy.qzv`, `taxa-bar-plots.qzv`
- `rooted-tree.qza` and exported Newick tree
- `diversity/` outputs such as Shannon vectors, UniFrac distances, PCoA, and Emperor plots

QIIME2 visualizations (`.qzv`) can be viewed at <https://view.qiime2.org>.
