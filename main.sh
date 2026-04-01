#!/usr/bin/env bash

# QIIME2 Pipeline with Checkpoint System
# This script processes demultiplexed paired-end FASTQ files using QIIME2 to generate an OTU table and representative sequences
# If unfamiliar with QIIME, see https://amplicon-docs.qiime2.org/en/latest/explanations/getting-started.html

set -euo pipefail  # Exit on error, undefined variables, and pipe failures

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
CHECKPOINT_DIR="${PROJECT_DIR}/Logs/checkpoints"
LOG_FILE="${PROJECT_DIR}/Logs/qiime2_pipeline.log"
ENV_NAME="qiime"
# Non-interactive / batch mode for HPC (Slurm) runs. When true, prompts are skipped
# and environment variables or sensible defaults are used instead.
NON_INTERACTIVE="false"
# Mode: 'denoise' (DADA2) or 'cluster' (OTU clustering from demux)
# Users can override with -m|--mode when running the script
MODE="denoise"

# Create checkpoint directory
mkdir -p "${CHECKPOINT_DIR}"

# Logging function
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "${LOG_FILE}"
}

# Error handling function
error_exit() {
    log "ERROR: $1"
    exit 1
}

# Checkpoint functions
create_checkpoint() {
    local step_name="$1"
    echo "$(date '+%Y-%m-%d %H:%M:%S')" > "${CHECKPOINT_DIR}/${step_name}.checkpoint"
    log "Checkpoint created for step: ${step_name}"
}

check_checkpoint() {
    local step_name="$1"
    if [[ -f "${CHECKPOINT_DIR}/${step_name}.checkpoint" ]]; then
        log "Step ${step_name} already completed. Skipping..."
        return 0
    else
        return 1
    fi
}

# Function to pause and wait for user input
pause_script() {
    local step_name="$1"
    local message="${2:-Press Enter to continue to the next step, or Ctrl+C to exit}"
    
    log "Completed step: ${step_name}"
    echo ""
    echo "=================================================="
    echo "Step completed: ${step_name}"
    echo "Checkpoint saved at: ${CHECKPOINT_DIR}/${step_name}.checkpoint"
    echo "=================================================="
    echo "${message}"
    # In non-interactive mode (batch), do not wait for user input
    if [[ "${NON_INTERACTIVE}" == "true" ]]; then
        log "Non-interactive mode: continuing without pause"
    else
        read -r
    fi
}

# Step 1: Environment Setup
step1_environment_setup() {
    if check_checkpoint "step1_environment_setup"; then
        return 0
    fi
    
    log "Starting Step 1: Environment Setup"
    
    # Check if conda is available
    if ! command -v conda &> /dev/null; then
        error_exit "Conda is not installed or not in PATH"
    fi
    
    # Initialize conda for this shell session
    source "$(conda info --base)/etc/profile.d/conda.sh"
    
    # Check if environment already exists
    if conda env list | grep -q "^${ENV_NAME} "; then
        log "QIIME2 environment '${ENV_NAME}' already exists"
        echo ""
        echo "=================================================="
        echo "Existing environment detected: ${ENV_NAME}"
        echo "=================================================="
        echo "1) Use existing environment"
        echo "2) Remove and recreate environment"
        echo "3) Exit and use a different environment name"
        echo ""
        read -p "Enter choice (1-3): " ENV_CHOICE
        
        case ${ENV_CHOICE} in
            1)
                log "Using existing QIIME2 environment: ${ENV_NAME}"
                ;;
            2)
                log "Removing existing environment and recreating..."
                conda env remove --name "${ENV_NAME}" -y || error_exit "Failed to remove existing environment"
                
                log "Creating new QIIME2 environment..."
                if ! conda env create \
                    --name "${ENV_NAME}" \
                    --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2025.10/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml; then
                    
                    log "Initial environment creation failed. Trying with flexible channel priority..."
                    conda config --set channel_priority flexible
                    
                    conda env create \
                        --name "${ENV_NAME}" \
                        --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2025.10/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml || \
                        error_exit "Failed to create QIIME2 environment"
                fi
                ;;
            3)
                log "Exiting. Please set ENV_NAME variable in the script to use a different environment name."
                exit 0
                ;;
            *)
                log "Invalid choice. Using existing environment."
                ;;
        esac
    else
        # Environment doesn't exist, create it
        log "Creating QIIME2 environment..."
        if ! conda env create \
            --name "${ENV_NAME}" \
            --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2025.10/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml; then
            
            log "Initial environment creation failed. Trying with flexible channel priority..."
            conda config --set channel_priority flexible
            
            conda env create \
                --name "${ENV_NAME}" \
                --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2025.10/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml || \
                error_exit "Failed to create QIIME2 environment"
        fi
    fi
    
    log "Activating QIIME2 environment..."
    conda activate "${ENV_NAME}" || error_exit "Failed to activate QIIME2 environment"
    
    # Verify installation
    qiime --version || error_exit "QIIME2 installation verification failed"
    
    create_checkpoint "step1_environment_setup"
    pause_script "Environment Setup" "Environment is ready. Make sure your manifest.tsv file is properly formatted before continuing."
}
# Step 2: Import Data
step2_import_data() {
    if check_checkpoint "step2_import_data"; then
        return 0
    fi
    
    log "Starting Step 2: Import demultiplexed paired-end FASTQ files"
    
    # Ensure environment is activated
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "${ENV_NAME}" || error_exit "Failed to activate QIIME2 environment"
    
    # Check if manifest file exists - look in Data/raw_data/16S first, then current directory
    MANIFEST_PATH=""
    if [[ -f "Data/raw_data/16S/manifest.tsv" ]]; then
        MANIFEST_PATH="Data/raw_data/16S/manifest.tsv"
        log "Using manifest file: Data/raw_data/16S/manifest.tsv"
    elif [[ -f "manifest.tsv" ]]; then
        MANIFEST_PATH="manifest.tsv"
        log "Using manifest file: manifest.tsv (in working directory)"
    else
        error_exit "manifest.tsv file not found. Please create this file in Data/raw_data/16S/ or current directory according to QIIME2 specifications."
    fi
    
    log "Importing paired-end data from ${MANIFEST_PATH}..."
    qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path "${MANIFEST_PATH}" \
        --output-path Data/processed_data/demux-paired-end.qza \
        --input-format PairedEndFastqManifestPhred33V2 || \
        error_exit "Failed to import data"
    
    # Verify output file was created
    if [[ ! -f "Data/processed_data/demux-paired-end.qza" ]]; then
        error_exit "Output file Data/processed_data/demux-paired-end.qza was not created"
    fi
    
    create_checkpoint "step2_import_data"
    pause_script "Data Import" "Data imported successfully. Review the Data/processed_data/demux-paired-end.qza file before proceeding."
}

# Step 3: Visualize Demux Data
step3_visualize_demux() {
    if check_checkpoint "step3_visualize_demux"; then
        return 0
    fi
    
    log "Starting Step 3: Visualize demux data for quality check"
    
    # Ensure environment is activated
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "${ENV_NAME}" || error_exit "Failed to activate QIIME2 environment"
    
    # Check if input file exists - look in Data/processed_data first, then current directory
    DEMUX_PATH=""
    if [[ -f "Data/processed_data/demux-paired-end.qza" ]]; then
        DEMUX_PATH="Data/processed_data/demux-paired-end.qza"
    elif [[ -f "demux-paired-end.qza" ]]; then
        DEMUX_PATH="demux-paired-end.qza"
    else
        error_exit "Input file demux-paired-end.qza not found in Data/processed_data/ or current directory. Please run step 2 first."
    fi
    
    log "Creating demux visualization from ${DEMUX_PATH}..."
    qiime demux summarize \
        --i-data "${DEMUX_PATH}" \
        --o-visualization Data/processed_data/demux-paired-end.qzv || \
        error_exit "Failed to create demux visualization"
    
    # Verify output file was created
    if [[ ! -f "Data/processed_data/demux-paired-end.qzv" ]]; then
        error_exit "Output file Data/processed_data/demux-paired-end.qzv was not created"
    fi
    
    create_checkpoint "step3_visualize_demux"
    pause_script "Demux Visualization" "Visualization created at Data/processed_data/demux-paired-end.qzv. Please drag this file to https://view.qiime2.org to check quality plots and determine truncation lengths for the next step."
}

# Step 4: Remove Primers/Adapters (Optional)
step4_remove_primers() {
    if check_checkpoint "step4_remove_primers"; then
        return 0
    fi
    
    log "Starting Step 4: Remove primers/adapters (optional)"
    
    # Ensure environment is activated
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "${ENV_NAME}" || error_exit "Failed to activate QIIME2 environment"
    
    echo ""
    echo "=================================================="
    echo "Primer/Adapter Removal"
    echo "=================================================="
    echo ""
    echo "Do you need to remove primers/adapters from your sequences?"
    echo ""
    echo "Common 16S primer sets:"
    echo "1) 515F/806R (V4 region) - Earth Microbiome Project"
    echo "   Forward: GTGCCAGCMGCCGCGGTAA"
    echo "   Reverse: GGACTACHVGGGTWTCTAAT"
    echo ""
    echo "2) Custom primers (you will enter them)"
    echo "3) Skip primer removal (already removed or not needed)"
    echo ""
    if [[ "${NON_INTERACTIVE}" == "true" && -n "${PRIMER_CHOICE:-}" ]]; then
        log "Non-interactive: using PRIMER_CHOICE=${PRIMER_CHOICE}"
    elif [[ "${NON_INTERACTIVE}" == "true" ]]; then
        # default: skip primer removal in batch mode unless PRIMER_CHOICE is provided
        PRIMER_CHOICE=3
        log "Non-interactive: default PRIMER_CHOICE=3 (skip primer removal)"
    else
        read -p "Enter choice (1-3): " PRIMER_CHOICE
    fi
    
    case ${PRIMER_CHOICE} in
        1)
            PRIMER_F="GTGCCAGCMGCCGCGGTAA"
            PRIMER_R="GGACTACHVGGGTWTCTAAT"
            log "Using 515F/806R primers (V4 region)"
            ;;
        2)
            echo ""
            read -p "Enter forward primer sequence (5' to 3'): " PRIMER_F
            read -p "Enter reverse primer sequence (5' to 3'): " PRIMER_R
            log "Using custom primers: F=${PRIMER_F}, R=${PRIMER_R}"
            ;;
        3)
            log "Skipping primer removal"
            create_checkpoint "step3b_remove_primers"
            return 0
            ;;
        *)
            log "Invalid choice. Skipping primer removal"
            create_checkpoint "step3b_remove_primers"
            return 0
            ;;
    esac

    # Check if input file exists
    DEMUX_PATH=""
    if [[ -f "Data/processed_data/demux-paired-end.qza" ]]; then
        DEMUX_PATH="Data/processed_data/demux-paired-end.qza"
    elif [[ -f "demux-paired-end.qza" ]]; then
        DEMUX_PATH="demux-paired-end.qza"
    else
        error_exit "Input file demux-paired-end.qza not found. Please run step 2 first."
    fi

    log "Removing primers using cutadapt..."
    log "Forward primer: ${PRIMER_F}"
    log "Reverse primer: ${PRIMER_R}"

    qiime cutadapt trim-paired \
        --i-demultiplexed-sequences "${DEMUX_PATH}" \
        --p-front-f "${PRIMER_F}" \
        --p-front-r "${PRIMER_R}" \
        --p-match-read-wildcards \
        --p-match-adapter-wildcards \
        --o-trimmed-sequences Data/processed_data/demux-trimmed.qza \
        --verbose || error_exit "Primer removal failed"

    # Verify output file was created
    if [[ ! -f "Data/processed_data/demux-trimmed.qza" ]]; then
        error_exit "Output file Data/processed_data/demux-trimmed.qza was not created"
    fi

    # Create visualization of trimmed data
    log "Creating visualization of trimmed sequences..."
    qiime demux summarize \
        --i-data Data/processed_data/demux-trimmed.qza \
        --o-visualization Data/processed_data/demux-trimmed.qzv || \
        log "Warning: Failed to create trimmed data visualization"

    create_checkpoint "step4_remove_primers"

    echo ""
    echo "=================================================="
    echo "Primer Removal Complete!"
    echo "=================================================="
    echo ""
    echo "Generated files:"
    echo "- Data/processed_data/demux-trimmed.qza (primer-trimmed sequences)"
    echo "- Data/processed_data/demux-trimmed.qzv (visualization)"
    echo ""
    echo "IMPORTANT: Use demux-trimmed.qza for downstream analysis (Step 5)"
    echo ""

    pause_script "Primer Removal" "Primers removed. Review demux-trimmed.qzv to verify primer removal before proceeding."
}

# Step 5: DADA2 Denoising
step5_dada2_denoising() {
    if check_checkpoint "step5_dada2_denoising"; then
        return 0
    fi
    
    log "Starting Step 5: DADA2 denoising and feature table generation"
    
    # Ensure environment is activated
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "${ENV_NAME}" || error_exit "Failed to activate QIIME2 environment"
    
    # Check if input file exists - prefer trimmed version if available
    DEMUX_PATH=""
    if [[ -f "Data/processed_data/demux-trimmed.qza" ]]; then
        DEMUX_PATH="Data/processed_data/demux-trimmed.qza"
        log "Using primer-trimmed sequences: ${DEMUX_PATH}"
    elif [[ -f "Data/processed_data/demux-paired-end.qza" ]]; then
        DEMUX_PATH="Data/processed_data/demux-paired-end.qza"
        log "Using original sequences (primers not removed): ${DEMUX_PATH}"
    elif [[ -f "demux-paired-end.qza" ]]; then
        DEMUX_PATH="demux-paired-end.qza"
    else
        error_exit "Input file not found. Please run step 2 (and optionally step 3b) first."
    fi
    
    # Get truncation lengths from user
    echo ""
    echo "Based on the quality plots from Data/processed_data/demux-paired-end.qzv:"
    if [[ "${NON_INTERACTIVE}" == "true" && -n "${TRUNC_LEN_F:-}" ]]; then
        log "Non-interactive: using TRUNC_LEN_F=${TRUNC_LEN_F}"
    elif [[ "${NON_INTERACTIVE}" == "true" ]]; then
        TRUNC_LEN_F=250
        log "Non-interactive: default TRUNC_LEN_F=${TRUNC_LEN_F}"
    else
        read -p "Enter forward truncation length (default: 250): " TRUNC_LEN_F
    fi

    if [[ "${NON_INTERACTIVE}" == "true" && -n "${TRUNC_LEN_R:-}" ]]; then
        log "Non-interactive: using TRUNC_LEN_R=${TRUNC_LEN_R}"
    elif [[ "${NON_INTERACTIVE}" == "true" ]]; then
        TRUNC_LEN_R=250
        log "Non-interactive: default TRUNC_LEN_R=${TRUNC_LEN_R}"
    else
        read -p "Enter reverse truncation length (default: 250): " TRUNC_LEN_R
    fi

    TRUNC_LEN_F=${TRUNC_LEN_F:-250}
    TRUNC_LEN_R=${TRUNC_LEN_R:-250}
    
    log "Using truncation lengths: Forward=${TRUNC_LEN_F}, Reverse=${TRUNC_LEN_R}"
    
    # Get number of threads to use
    echo ""
    echo "DADA2 can use multiple CPU cores to speed up processing."
    if [[ "${NON_INTERACTIVE}" == "true" && -n "${N_THREADS:-}" ]]; then
        log "Non-interactive: using N_THREADS=${N_THREADS}"
    elif [[ "${NON_INTERACTIVE}" == "true" ]]; then
        N_THREADS=0
        log "Non-interactive: default N_THREADS=0 (all available)"
    else
        read -p "Enter number of threads (0 = all available, default: 0): " N_THREADS
        N_THREADS=${N_THREADS:-0}
    fi
    
    log "Using ${N_THREADS} threads (0 = all available cores)"
    
    # Create Results directory for denoise mode
    mkdir -p Results/denoise_mode
    
    log "Running DADA2 denoising (this may take a while)..."
    qiime dada2 denoise-paired \
        --i-demultiplexed-seqs "${DEMUX_PATH}" \
        --p-trunc-len-f "${TRUNC_LEN_F}" \
        --p-trunc-len-r "${TRUNC_LEN_R}" \
        --p-trim-left-f 0 \
        --p-trim-left-r 0 \
        --p-n-threads "${N_THREADS}" \
        --o-table Results/denoise_mode/table.qza \
        --o-representative-sequences Results/denoise_mode/rep-seqs.qza \
        --o-denoising-stats Results/denoise_mode/denoising-stats.qza \
        --o-base-transition-stats Results/denoise_mode/base-transition-stats.qza || \
        error_exit "DADA2 denoising failed"
    
    # Verify output files were created
    for file in Results/denoise_mode/table.qza Results/denoise_mode/rep-seqs.qza Results/denoise_mode/denoising-stats.qza Results/denoise_mode/base-transition-stats.qza; do
        if [[ ! -f "${file}" ]]; then
            error_exit "Output file ${file} was not created"
        fi
    done
    
    # Create visualizations of DADA2 results
    log "Creating visualizations of DADA2 results..."
    
    log "Creating denoising stats visualization..."
    qiime metadata tabulate \
        --m-input-file Results/denoise_mode/denoising-stats.qza \
        --o-visualization Results/denoise_mode/denoising-stats.qzv || \
        log "Warning: Failed to create denoising stats visualization"
    
    log "Creating feature table summary..."
    qiime feature-table summarize \
        --i-table Results/denoise_mode/table.qza \
        --o-visualization Results/denoise_mode/table.qzv || \
        log "Warning: Failed to create feature table visualization"
    
    log "Creating representative sequences summary..."
    qiime feature-table tabulate-seqs \
        --i-data Results/denoise_mode/rep-seqs.qza \
        --o-visualization Results/denoise_mode/rep-seqs.qzv || \
        log "Warning: Failed to create rep-seqs visualization"
    
    create_checkpoint "step5_dada2_denoising"
    
    echo ""
    echo "=================================================="
    echo "DADA2 Denoising Complete!"
    echo "=================================================="
    echo ""
    echo "Generated files:"
    echo "1. Results/denoise_mode/table.qza - Feature table (ASV counts)"
    echo "2. Results/denoise_mode/rep-seqs.qza - Representative sequences"
    echo "3. Results/denoise_mode/denoising-stats.qza - Denoising statistics"
    echo ""
    echo "Visualizations (open at https://view.qiime2.org):"
    echo "1. Results/denoise_mode/denoising-stats.qzv - See how many reads passed filters"
    echo "2. Results/denoise_mode/table.qzv - Feature table summary (# features, sampling depth)"
    echo "3. Results/denoise_mode/rep-seqs.qzv - Browse representative sequences"
    echo ""
    
    pause_script "DADA2 Denoising" "DADA2 processing completed. Review the .qzv files at https://view.qiime2.org before proceeding."
}

# Step 6: Decontamination

step6_decontamination() {
    if check_checkpoint "step6_decontamination"; then
        return 0
    fi

    log "Starting Step 6: Decontamination"

    # Ensure environment is activated
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "${ENV_NAME}" || error_exit "Failed to activate QIIME2 environment"

    # Determine which mode we're in and set paths
    if [[ "${MODE}" == "denoise" ]]; then
        TABLE_PATH="Results/denoise_mode/table.qza"
        REP_SEQS_PATH="Results/denoise_mode/rep-seqs.qza"
        OUTPUT_DIR="Results/denoise_mode"
    elif [[ "${MODE}" == "cluster" ]]; then
        TABLE_PATH="Results/cluster_mode/table.qza"
        REP_SEQS_PATH="Results/cluster_mode/rep-seqs.qza"
        OUTPUT_DIR="Results/cluster_mode"
    else
        error_exit "Unknown mode: ${MODE}"
    fi

    # Check if required files exist
    if [[ ! -f "${TABLE_PATH}" ]]; then
        error_exit "Feature table not found. Please run step 5 first."
    fi

    echo ""
    echo "=================================================="
    echo "Decontamination Step"
    echo "=================================================="
    echo ""
    echo "Do you want to perform decontamination to remove negative controls?"
    echo "1) Yes"
    echo "2) No (skip decontamination)"
    echo ""
    if [[ "${NON_INTERACTIVE}" == "true" && -n "${DECONTAMINATION_CHOICE:-}" ]]; then
        log "Non-interactive: using DECONTAMINATION_CHOICE=${DECONTAMINATION_CHOICE}"
    elif [[ "${NON_INTERACTIVE}" == "true" ]]; then
        DECONTAMINATION_CHOICE=2
        log "Non-interactive: default DECONTAMINATION_CHOICE=2 (skip decontamination)"
    else
        read -p "Enter choice (1-2): " DECONTAMINATION_CHOICE
    fi

    case ${DECONTAMINATION_CHOICE} in
        1)
            # Non-interactive: the pipeline no longer prompts for example IDs or generates metadata automatically.
            # Instead, the user must provide a prepared metadata file with a 'sample-id' column and a
            # 'control_status' column (values: 'control' or 'not_control').
            
            # Priority for metadata file:
            # 1) Environment variable DECONTAM_METADATA (full path)
            # 2) Data/metadata/decontam-metadata.tsv
            # 3) Data/metadata/metadata.tsv
            
            if [[ -n "${DECONTAM_METADATA:-}" && -f "${DECONTAM_METADATA}" ]]; then
                META_FILE="${DECONTAM_METADATA}"
                log "Using metadata from DECONTAM_METADATA=${META_FILE}"
            elif [[ -f "Data/metadata/decontam-metadata.tsv" ]]; then
                META_FILE="Data/metadata/decontam-metadata.tsv"
                log "Using metadata from Data/metadata/decontam-metadata.tsv"
            elif [[ -f "Data/metadata/metadata.tsv" ]]; then
                META_FILE="Data/metadata/metadata.tsv"
                log "Using metadata from Data/metadata/metadata.tsv"
            else
                log "No decontamination metadata file found. Please provide one via DECONTAM_METADATA or place Data/metadata/decontam-metadata.tsv"
                log "Skipping decontamination step."
                create_checkpoint "step6_decontamination"
            fi

            if [[ -n "${META_FILE:-}" ]]; then
                log "Running decontam identification with prevalence method..."
                qiime quality-control decontam-identify \
                    --i-table "${TABLE_PATH}" \
                    --m-metadata-file "${META_FILE}" \
                    --p-method prevalence \
                    --p-prev-control-column control_status \
                    --p-prev-control-indicator control \
                    --o-decontam-scores "${OUTPUT_DIR}/decontam-scores.qza" || \
                    error_exit "Decontam identification failed"

                log "Creating visualization of decontam scores..."
                qiime quality-control decontam-score-viz \
                    --i-decontam-scores "${OUTPUT_DIR}/decontam-scores.qza" \
                    --i-table "${TABLE_PATH}" \
                    --o-visualization "${OUTPUT_DIR}/decontam-scores.qzv" || \
                    error_exit "Decontam visualization failed"

                log "Exporting decontam scores to TSV format..."
                mkdir -p "${OUTPUT_DIR}/decontam_scores_export"
                qiime tools export \
                    --input-path "${OUTPUT_DIR}/decontam-scores.qza" \
                    --output-path "${OUTPUT_DIR}/decontam_scores_export" || \
                    error_exit "Failed to export decontam scores"

                # Find the actual scores file (QIIME2 may name it differently)
                SCORES_FILE=$(find "${OUTPUT_DIR}/decontam_scores_export" -name "*.tsv" | head -1)
                if [[ -z "${SCORES_FILE}" ]]; then
                    error_exit "No TSV file found after export in ${OUTPUT_DIR}/decontam_scores_export"
                fi
                log "Found scores file: ${SCORES_FILE}"

                log "Extracting contaminant feature IDs based on score threshold..."
                # Get score threshold from user or use default (0.5 for prevalence method)
                if [[ "${NON_INTERACTIVE}" == "true" && -n "${DECONTAM_THRESHOLD:-}" ]]; then
                    SCORE_THRESHOLD="${DECONTAM_THRESHOLD}"
                    log "Non-interactive: using DECONTAM_THRESHOLD=${SCORE_THRESHOLD}"
                elif [[ "${NON_INTERACTIVE}" == "true" ]]; then
                    SCORE_THRESHOLD=0.5
                    log "Non-interactive: using default SCORE_THRESHOLD=0.5"
                else
                    echo ""
                    echo "The decontam score indicates probability of being a contaminant (0-1)."
                    echo "Higher scores = more likely to be a contaminant."
                    echo "Recommended threshold for prevalence method: 0.5"
                    read -p "Enter score threshold for identifying contaminants (default: 0.5): " SCORE_THRESHOLD
                    SCORE_THRESHOLD=${SCORE_THRESHOLD:-0.5}
                fi

                log "Using score threshold: ${SCORE_THRESHOLD}"

                # Debug: show the file structure
                log "First 10 lines of ${SCORES_FILE}:"
                head -10 "${SCORES_FILE}" | tee -a "${LOG_FILE}"

                # Extract contaminant feature IDs using Python
                python3 - "${SCORES_FILE}" "${SCORE_THRESHOLD}" "${OUTPUT_DIR}" << 'PYTHON_SCRIPT'
import csv
import sys

scores_file = sys.argv[1]
threshold = float(sys.argv[2])
output_dir = sys.argv[3]

contaminant_ids_file = f"{output_dir}/decontam_scores_export/contaminant-ids.txt"

contaminant_ids = []
try:
    with open(scores_file, 'r') as f:
        # Read the file and print column headers for debugging
        reader = csv.DictReader(f, delimiter='\t')
        headers = reader.fieldnames
        print(f"DEBUG: Found columns: {headers}", file=sys.stderr)
        
        row_count = 0
        for row in reader:
            row_count += 1
            try:
                # Print first row for debugging
                if row_count == 1:
                    print(f"DEBUG: First row: {row}", file=sys.stderr)
                
                # Get the score from column 'p' (decontam p-value)
                # Also try 'score' and other variations as fallback
                score = None
                score_col_found = None
                
                for col in ['p', 'score', 'Score', 'p-value', 'p_value', 'decontam score']:
                    if col in row:
                        try:
                            score = float(row[col])
                            score_col_found = col
                            break
                        except ValueError:
                            pass
                
                if row_count == 1 and score is not None:
                    print(f"DEBUG: Score column: '{score_col_found}', value: {score}", file=sys.stderr)
                
                if score is not None and score > threshold:
                    # Get feature ID from '#OTU ID' or similar
                    feature_id = None
                    for col in ['#OTU ID', '#OTU id', 'OTU ID', 'feature-id', 'feature_id', 'id']:
                        if col in row:
                            feature_id = row[col]
                            break
                    
                    if feature_id:
                        contaminant_ids.append(feature_id)
                        if len(contaminant_ids) <= 3:
                            print(f"DEBUG: Found contaminant: {feature_id} with score {score}", file=sys.stderr)
            except Exception as e:
                if row_count == 1:
                    print(f"DEBUG: Error processing row 1: {e}", file=sys.stderr)
                continue
    
    print(f"DEBUG: Total rows processed: {row_count}", file=sys.stderr)
    print(f"DEBUG: Contaminants found: {len(contaminant_ids)}", file=sys.stderr)
    
    # Write contaminant IDs to file in metadata format
    with open(contaminant_ids_file, 'w') as f:
        f.write('feature-id\n')
        for fid in contaminant_ids:
            f.write(f'{fid}\n')
    
    print(f"Identified {len(contaminant_ids)} contaminant features with score > {threshold}")
except Exception as e:
    print(f"Error: {e}", file=sys.stderr)
    import traceback
    traceback.print_exc()
    sys.exit(1)
PYTHON_SCRIPT

                if [[ ! -f "${OUTPUT_DIR}/decontam_scores_export/contaminant-ids.txt" ]]; then
                    error_exit "Failed to extract contaminant feature IDs"
                fi

                # Count contaminants
                CONTAM_COUNT=$(tail -n +2 "${OUTPUT_DIR}/decontam_scores_export/contaminant-ids.txt" 2>/dev/null | wc -l)
                log "Identified ${CONTAM_COUNT} contaminant features"

                # Only proceed with filtering if contaminants were found
                if [[ ${CONTAM_COUNT} -eq 0 ]]; then
                    log "No contaminants identified with score threshold ${SCORE_THRESHOLD}"
                    log "Using original table for downstream analysis"
                    FINAL_TABLE="${TABLE_PATH}"
                    FINAL_REP_SEQS="${REP_SEQS_PATH}"
                    FINAL_DESC="Original feature table (no contaminants identified)"
                else
                    log "Filtering contaminant features from the table..."
                    qiime feature-table filter-features \
                        --i-table "${TABLE_PATH}" \
                        --m-metadata-file "${OUTPUT_DIR}/decontam_scores_export/contaminant-ids.txt" \
                        --p-exclude-ids \
                        --o-filtered-table "${OUTPUT_DIR}/table-no-contam.qza" || \
                        error_exit "Failed to filter contaminant features from table"

                    log "Filtering representative sequences to match filtered table..."
                    qiime feature-table filter-seqs \
                        --i-data "${REP_SEQS_PATH}" \
                        --i-table "${OUTPUT_DIR}/table-no-contam.qza" \
                        --o-filtered-data "${OUTPUT_DIR}/rep-seqs-no-contam.qza" || \
                        error_exit "Failed to filter representative sequences"
                    
                    FINAL_TABLE="${OUTPUT_DIR}/table-no-contam.qza"
                    FINAL_REP_SEQS="${OUTPUT_DIR}/rep-seqs-no-contam.qza"
                    FINAL_DESC="Feature table with ${CONTAM_COUNT} contaminant features removed"
                fi

                # Ask about removing negative control samples
                echo ""
                echo "Do you want to remove negative control samples from the analysis?"
                echo "This is recommended as they were only used for contamination identification."
                if [[ "${NON_INTERACTIVE}" == "true" ]]; then
                    REMOVE_CONTROLS=${REMOVE_CONTROLS:-y}
                    log "Non-interactive: REMOVE_CONTROLS=${REMOVE_CONTROLS}"
                else
                    read -p "Remove negative controls? (y/n): " REMOVE_CONTROLS
                    REMOVE_CONTROLS=${REMOVE_CONTROLS:-y}
                fi

                if [[ "${REMOVE_CONTROLS}" =~ ^[Yy]$ ]]; then
                    log "Removing negative control samples from the table..."
                    qiime feature-table filter-samples \
                        --i-table "${FINAL_TABLE}" \
                        --m-metadata-file "${META_FILE}" \
                        --p-where "control_status!='control'" \
                        --o-filtered-table "${OUTPUT_DIR}/table-clean.qza" || \
                        error_exit "Failed to remove negative control samples"

                    # Update table path for final output
                    FINAL_TABLE="${OUTPUT_DIR}/table-clean.qza"
                    FINAL_DESC="${FINAL_DESC} (negative controls removed)"
                else
                    log "Keeping negative control samples in the final table..."
                fi

                # Verify output files were created
                for file in "${OUTPUT_DIR}/decontam-scores.qza" "${OUTPUT_DIR}/decontam-scores.qzv" "${FINAL_TABLE}" "${FINAL_REP_SEQS}"; do
                    if [[ ! -f "${file}" ]]; then
                        error_exit "Output file ${file} was not created"
                    fi
                done

                create_checkpoint "step6_decontamination"

                echo ""
                echo "=================================================="
                echo "Decontamination Complete!"
                echo "=================================================="
                echo ""
                echo "Generated files:"
                echo "1. ${OUTPUT_DIR}/decontam-scores.qza - Decontam scores for each feature"
                echo "2. ${OUTPUT_DIR}/decontam-scores.qzv - Visualization of decontam scores"
                echo "3. ${OUTPUT_DIR}/decontam_scores_export/scores.tsv - Decontam scores in TSV format"
                echo "4. ${OUTPUT_DIR}/decontam_scores_export/contaminant-ids.txt - Identified contaminants"
                echo "5. ${FINAL_TABLE} - ${FINAL_DESC}"
                echo "6. ${FINAL_REP_SEQS} - Representative sequences matching filtered table"
                echo ""
                echo "Decontamination statistics:"
                echo "- Score threshold used: ${SCORE_THRESHOLD}"
                echo "- Contaminant features removed: ${CONTAM_COUNT}"
                echo ""

                pause_script "Decontamination" "Decontamination completed. Cleaned data ready for taxonomic classification."
            fi
            ;;
        2)
            log "Skipping decontamination step."
            create_checkpoint "step6_decontamination"
            ;;
        *)
            log "Invalid choice. Skipping decontamination step."
            create_checkpoint "step6_decontamination"
            ;;
    esac
}

# Step 5 Alternative: Build OTU table and rep-seqs directly from demux (cluster mode)
step5_cluster_from_demux() {
    if check_checkpoint "step5_cluster_from_demux"; then
        return 0
    fi

    log "Starting Step 5 (cluster mode): Generate OTU table and representative sequences from demux"

    # Ensure environment is activated
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "${ENV_NAME}" || error_exit "Failed to activate QIIME2 environment"

    # Check if input file exists - prefer trimmed version if available
    DEMUX_PATH=""
    if [[ -f "Data/processed_data/demux-trimmed.qza" ]]; then
        DEMUX_PATH="Data/processed_data/demux-trimmed.qza"
        log "Using primer-trimmed sequences: ${DEMUX_PATH}"
    elif [[ -f "Data/processed_data/demux-paired-end.qza" ]]; then
        DEMUX_PATH="Data/processed_data/demux-paired-end.qza"
        log "Using original sequences (primers not removed): ${DEMUX_PATH}"
    elif [[ -f "demux-paired-end.qza" ]]; then
        DEMUX_PATH="demux-paired-end.qza"
    else
        error_exit "Input file not found. Please run step 2 (and optionally step 3b) first."
    fi

    # Create Results directory for cluster mode
    mkdir -p Results/cluster_mode

    log "Merging paired reads using QIIME2 vsearch..."
    qiime vsearch merge-pairs \
        --i-demultiplexed-seqs "${DEMUX_PATH}" \
        --o-merged-sequences Results/cluster_mode/joined.qza \
        --o-unmerged-sequences Results/cluster_mode/unmerged.qza || error_exit "merge-pairs failed"

    log "Quality-filtering merged reads..."
    qiime quality-filter q-score \
        --i-demux Results/cluster_mode/joined.qza \
        --o-filtered-sequences Results/cluster_mode/filtered-seqs.qza \
        --o-filter-stats Results/cluster_mode/filtered-stats.qza || error_exit "quality-filter failed"

    log "Dereplicating sequences..."
    qiime vsearch dereplicate-sequences \
        --i-sequences Results/cluster_mode/filtered-seqs.qza \
        --o-dereplicated-table Results/cluster_mode/derep-table.qza \
        --o-dereplicated-sequences Results/cluster_mode/derep-seqs.qza || error_exit "dereplicate-sequences failed"
    
    log "Removing chimeras from dereplicated sequences..."
    qiime vsearch uchime-denovo \
        --i-table Results/cluster_mode/derep-table.qza \
        --i-sequences Results/cluster_mode/derep-seqs.qza \
        --o-chimeras Results/cluster_mode/chimeras.qza \
        --o-nonchimeras Results/cluster_mode/rep-seqs.qza \
        --o-stats Results/cluster_mode/chimera-stats.qza || error_exit "chimera removal failed"
    
    log "Filtering feature table to remove chimeric sequences..."
    qiime feature-table filter-features \
        --i-table Results/cluster_mode/derep-table.qza \
        --m-metadata-file Results/cluster_mode/rep-seqs.qza \
        --o-filtered-table Results/cluster_mode/table.qza || error_exit "table filtering failed"

    # At this point we have table.qza and rep-seqs.qza; offer clustering options
    echo ""
    echo "Choose clustering method to produce OTUs (from dereplicated sequences):"
    echo "1) de novo clustering (vsearch cluster-features-de-novo, 97%)"
    echo "2) closed-reference clustering (requires reference database)"
    echo "3) skip clustering (keep dereplicated sequences as features)"
    read -p "Enter choice (1-3): " CLUSTER_CHOICE

    case ${CLUSTER_CHOICE} in
        1)
            log "Running de novo clustering (97%)..."
            qiime vsearch cluster-features-de-novo \
                --i-table Results/cluster_mode/table.qza \
                --i-sequences Results/cluster_mode/rep-seqs.qza \
                --p-perc-identity 0.97 \
                --o-clustered-table Results/cluster_mode/table-dn-97.qza \
                --o-clustered-sequences Results/cluster_mode/rep-seqs-dn-97.qza || error_exit "de-novo clustering failed"
            # Promote outputs to standard names for downstream steps
            mv Results/cluster_mode/table-dn-97.qza Results/cluster_mode/table.qza
            mv Results/cluster_mode/rep-seqs-dn-97.qza Results/cluster_mode/rep-seqs.qza
            ;;
        2)
            if [[ ! -f "Data/reference_dbs/silva_97_otus.qza" ]] && [[ ! -f "silva_97_otus.qza" ]]; then
                log "Reference database not found in Data/reference_dbs/ or current directory. Cannot run closed-reference clustering."
            else
                REF_DB_PATH=""
                if [[ -f "Data/reference_dbs/silva_97_otus.qza" ]]; then
                    REF_DB_PATH="Data/reference_dbs/silva_97_otus.qza"
                else
                    REF_DB_PATH="silva_97_otus.qza"
                fi
                log "Running closed-reference clustering with ${REF_DB_PATH}..."
                qiime vsearch cluster-features-closed-reference \
                    --i-table Results/cluster_mode/table.qza \
                    --i-sequences Results/cluster_mode/rep-seqs.qza \
                    --i-reference-sequences "${REF_DB_PATH}" \
                    --p-perc-identity 0.97 \
                    --o-clustered-table Results/cluster_mode/table-cr-97.qza \
                    --o-clustered-sequences Results/cluster_mode/rep-seqs-cr-97.qza \
                    --o-unmatched-sequences Results/cluster_mode/unmatched.qza || error_exit "Closed-reference clustering failed"
                mv Results/cluster_mode/table-cr-97.qza Results/cluster_mode/table.qza
                mv Results/cluster_mode/rep-seqs-cr-97.qza Results/cluster_mode/rep-seqs.qza
            fi
            ;;
        3)
            log "Skipping additional clustering; using dereplicated features"
            ;;
        *)
            log "Invalid choice. Using dereplicated features"
            ;;
    esac

    # Verify outputs exist
    for file in Results/cluster_mode/table.qza Results/cluster_mode/rep-seqs.qza; do
        if [[ ! -f "${file}" ]]; then
            error_exit "Expected output ${file} was not created"
        fi
    done

    create_checkpoint "step5_cluster_from_demux"
    pause_script "Cluster-from-demux" "Clustering-from-demux completed. Outputs in Results/cluster_mode/: table.qza, rep-seqs.qza"
}

# Step 7: Taxonomic Classification
step7_taxonomic_classification() {
    if check_checkpoint "step7_taxonomic_classification"; then
        return 0
    fi
    
    log "Starting Step 7: Taxonomic Classification"
    
    # Ensure environment is activated
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "${ENV_NAME}" || error_exit "Failed to activate QIIME2 environment"
    
    # Determine which mode we're in and set paths
    if [[ "${MODE}" == "denoise" ]]; then
        TABLE_PATH="Results/denoise_mode/table.qza"
        REP_SEQS_PATH="Results/denoise_mode/rep-seqs.qza"
        OUTPUT_DIR="Results/denoise_mode"
        
        # Check if decontaminated files exist and prefer them
        if [[ -f "Results/denoise_mode/table-clean.qza" ]]; then
            TABLE_PATH="Results/denoise_mode/table-clean.qza"
            REP_SEQS_PATH="Results/denoise_mode/rep-seqs-no-contam.qza"
            log "Using decontaminated table and sequences"
        elif [[ -f "Results/denoise_mode/table-no-contam.qza" ]]; then
            TABLE_PATH="Results/denoise_mode/table-no-contam.qza"
            REP_SEQS_PATH="Results/denoise_mode/rep-seqs-no-contam.qza"
            log "Using contaminant-filtered table and sequences"
        fi
    elif [[ "${MODE}" == "cluster" ]]; then
        TABLE_PATH="Results/cluster_mode/table.qza"
        REP_SEQS_PATH="Results/cluster_mode/rep-seqs.qza"
        OUTPUT_DIR="Results/cluster_mode"
        
        # Check if decontaminated files exist and prefer them
        if [[ -f "Results/cluster_mode/table-clean.qza" ]]; then
            TABLE_PATH="Results/cluster_mode/table-clean.qza"
            REP_SEQS_PATH="Results/cluster_mode/rep-seqs-no-contam.qza"
            log "Using decontaminated table and sequences"
        elif [[ -f "Results/cluster_mode/table-no-contam.qza" ]]; then
            TABLE_PATH="Results/cluster_mode/table-no-contam.qza"
            REP_SEQS_PATH="Results/cluster_mode/rep-seqs-no-contam.qza"
            log "Using contaminant-filtered table and sequences"
        fi
    else
        error_exit "Unknown mode: ${MODE}"
    fi
    
    # Check if required files exist
    if [[ ! -f "${TABLE_PATH}" ]] || [[ ! -f "${REP_SEQS_PATH}" ]]; then
        error_exit "Required files not found. Please run step 5 first."
    fi
    
    # Check for classifier (detect but always ask user whether to use it)
    CLASSIFIER_PATH=""
    if [[ -f "Data/reference_dbs/silva-138-99-nb-classifier.qza" ]]; then
        CLASSIFIER_PATH="Data/reference_dbs/silva-138-99-nb-classifier.qza"
    elif [[ -f "silva-138-99-nb-classifier.qza" ]]; then
        CLASSIFIER_PATH="silva-138-99-nb-classifier.qza"
    elif [[ -f "Data/reference_dbs/classifier.qza" ]]; then
        CLASSIFIER_PATH="Data/reference_dbs/classifier.qza"
    elif [[ -f "classifier.qza" ]]; then
        CLASSIFIER_PATH="classifier.qza"
    elif [[ -f "Data/reference_dbs/silva-v4-classifier.qza" ]]; then
        CLASSIFIER_PATH="Data/reference_dbs/silva-v4-classifier.qza"
    elif [[ -f "silva-v4-classifier.qza" ]]; then
        CLASSIFIER_PATH="silva-v4-classifier.qza"
    fi

    # If a classifier was found, ask the user whether to use it or choose other options
    if [[ -n "${CLASSIFIER_PATH}" ]]; then
        echo ""
        echo "Detected existing classifier: ${CLASSIFIER_PATH}"
        if [[ "${NON_INTERACTIVE}" == "true" ]]; then
            log "Non-interactive: automatically using detected classifier: ${CLASSIFIER_PATH}"
        else
            read -p "Use this classifier? (y = use it / n = show options to build/download/choose another) [Y/n]: " USE_FOUND
            USE_FOUND=${USE_FOUND:-Y}
            if [[ "${USE_FOUND}" =~ ^[Yy]$ ]]; then
                log "Using detected classifier: ${CLASSIFIER_PATH}"
            else
                # Clear so the interactive menu below is shown
                CLASSIFIER_PATH=""
            fi
        fi
    fi
    
    if [[ -z "${CLASSIFIER_PATH}" ]]; then
        echo ""
        echo "=================================================="
        echo "Taxonomic Classifier Not Found"
        echo "=================================================="
        echo ""
        echo "Choose an option:"
        echo "1) Use an existing classifier file (you provide the path)"
        echo "2) Build a custom classifier (requires reference sequences)"
        echo "3) Download a pre-trained classifier"
        echo "4) Skip taxonomic classification"
        echo ""
        read -p "Enter choice (1-4): " CLASSIFIER_CHOICE
        
        # In non-interactive mode, allow CLASSIFIER_CHOICE env var or BUILD_CLASSIFIER to control behavior
        if [[ "${NON_INTERACTIVE}" == "true" && -n "${CLASSIFIER_CHOICE:-}" ]]; then
            log "Non-interactive: using CLASSIFIER_CHOICE=${CLASSIFIER_CHOICE}"
        elif [[ "${NON_INTERACTIVE}" == "true" && -n "${BUILD_CLASSIFIER:-}" ]]; then
            # map BUILD_CLASSIFIER values to numeric choices: silva->2, download->3, skip->4, custom->2 (handled later)
            case "${BUILD_CLASSIFIER}" in
                silva|Silva|SILVA)
                    CLASSIFIER_CHOICE=2
                    ;;
                download|pretrained)
                    CLASSIFIER_CHOICE=3
                    ;;
                skip|none)
                    CLASSIFIER_CHOICE=4
                    ;;
                custom)
                    CLASSIFIER_CHOICE=2
                    ;;
                *)
                    CLASSIFIER_CHOICE=4
                    ;;
            esac
            log "Non-interactive: mapped BUILD_CLASSIFIER=${BUILD_CLASSIFIER} to CLASSIFIER_CHOICE=${CLASSIFIER_CHOICE}"
        fi

        case ${CLASSIFIER_CHOICE} in
            1)
                read -p "Enter the full path to your classifier file: " USER_CLASSIFIER
                if [[ -f "${USER_CLASSIFIER}" ]]; then
                    CLASSIFIER_PATH="${USER_CLASSIFIER}"
                    log "Using classifier: ${CLASSIFIER_PATH}"
                else
                    error_exit "Classifier file not found: ${USER_CLASSIFIER}"
                fi
                ;;
            2)
                echo ""
                echo "Building a custom classifier using RESCRIPt and SILVA."
                echo "This will download, filter, extract, and train automatically."
                echo ""
                # Check for RESCRIPt installation and offer interactive install into the active env
                if ! qiime rescript --help &> /dev/null; then
                    echo ""
                    echo "RESCRIPt plugin is not installed in the active environment: ${ENV_NAME}"
                    if [[ "${NON_INTERACTIVE}" == "true" ]]; then
                        # In batch mode use INSTALL_RESCRIPT env var (yes/no). Default: no.
                        INSTALL_RESCRIPT=${INSTALL_RESCRIPT:-n}
                        log "Non-interactive: INSTALL_RESCRIPT=${INSTALL_RESCRIPT}"
                    else
                        read -p "Would you like to attempt to install RESCRIPt into '${ENV_NAME}' now? (y/n): " INSTALL_RESCRIPT
                    fi
                    if [[ "${INSTALL_RESCRIPT}" =~ ^[Yy]$ ]]; then
                        log "Attempting to install RESCRIPt into environment: ${ENV_NAME}"
                        # Try installing into the named environment. Use -n so it works whether or not env is active.
                        if ! conda install -n "${ENV_NAME}" -c conda-forge -c bioconda qiime2-rescript -y; then
                            log "Automatic installation failed. Please install manually with: conda activate ${ENV_NAME} && conda install -c conda-forge -c bioconda qiime2-rescript"
                            error_exit "Automatic RESCRIPt installation failed"
                        fi
                        # Reactivate environment to ensure new plugins are available
                        source "$(conda info --base)/etc/profile.d/conda.sh"
                        conda activate "${ENV_NAME}" || error_exit "Failed to reactivate environment after installing RESCRIPt"
                        if ! qiime rescript --help &> /dev/null; then
                            error_exit "RESCRIPt still not available after installation. Please check conda output or install manually in '${ENV_NAME}'."
                        fi
                        log "RESCRIPt successfully installed and available in '${ENV_NAME}'"
                    else
                        error_exit "RESCRIPt plugin is required to build a classifier with RESCRIPt. Install it in '${ENV_NAME}' and re-run the script."
                    fi
                fi
                mkdir -p Data/reference_dbs
                log "Downloading SILVA reference data (this may take several minutes)..."
                # Recent RESCRIPt versions use --o-silva-sequences and --o-silva-taxonomy
                qiime rescript get-silva-data \
                    --p-version '138.1' \
                    --p-target 'SSURef_NR99' \
                    --o-silva-sequences Data/reference_dbs/silva-seqs.qza \
                    --o-silva-taxonomy Data/reference_dbs/silva-tax.qza || error_exit "Failed to download SILVA data"
                log "Culling low-quality sequences..."
                qiime rescript cull-seqs \
                    --i-sequences Data/reference_dbs/silva-seqs.qza \
                    --o-clean-sequences Data/reference_dbs/silva-seqs-clean.qza || error_exit "Failed to cull sequences"
                log "Filtering sequences by length..."
                qiime rescript filter-seqs-length-by-taxon \
                    --i-sequences Data/reference_dbs/silva-seqs-clean.qza \
                    --i-taxonomy Data/reference_dbs/silva-tax.qza \
                    --p-min-lens 1000 \
                    --o-filtered-seqs Data/reference_dbs/silva-seqs-filtered.qza || error_exit "Failed to filter sequences by length"
                echo ""
                echo "Do you want to extract a specific region using primers?"
                echo "This is recommended if you used specific primers (e.g., 515F/806R)"
                if [[ "${NON_INTERACTIVE}" == "true" ]]; then
                    EXTRACT_PRIMERS=${EXTRACT_PRIMERS:-n}
                    log "Non-interactive: EXTRACT_PRIMERS=${EXTRACT_PRIMERS}"
                else
                    read -p "Extract primers? (y/n): " EXTRACT_PRIMERS
                fi
                if [[ "${EXTRACT_PRIMERS}" =~ ^[Yy]$ ]]; then
                    if [[ "${NON_INTERACTIVE}" == "true" ]]; then
                        # In batch mode require PRIMER_F/PRIMER_R env vars to be set
                        if [[ -z "${PRIMER_F:-}" || -z "${PRIMER_R:-}" ]]; then
                            error_exit "Non-interactive: PRIMER_F and PRIMER_R must be set when EXTRACT_PRIMERS=yes"
                        fi
                        log "Non-interactive: using primers F=${PRIMER_F} R=${PRIMER_R}"
                    else
                        read -p "Enter forward primer sequence: " PRIMER_F
                        read -p "Enter reverse primer sequence: " PRIMER_R
                    fi
                    log "Extracting amplicon region from reference sequences..."
                    qiime feature-classifier extract-reads \
                        --i-sequences Data/reference_dbs/silva-seqs-filtered.qza \
                        --p-f-primer "${PRIMER_F}" \
                        --p-r-primer "${PRIMER_R}" \
                        --o-reads Data/reference_dbs/silva-extracted.qza || error_exit "Failed to extract reads"
                    REF_SEQS="Data/reference_dbs/silva-extracted.qza"
                else
                    REF_SEQS="Data/reference_dbs/silva-seqs-filtered.qza"
                fi
                log "Training classifier (this may take 30+ minutes)..."
                qiime feature-classifier fit-classifier-naive-bayes \
                    --i-reference-reads "${REF_SEQS}" \
                    --i-reference-taxonomy Data/reference_dbs/silva-tax.qza \
                    --o-classifier Data/reference_dbs/custom-classifier.qza || error_exit "Classifier training failed"
                CLASSIFIER_PATH="Data/reference_dbs/custom-classifier.qza"
                log "Custom classifier created: ${CLASSIFIER_PATH}"
                ;;
            3)
                echo ""
                echo "To perform taxonomic classification, you need a pre-trained classifier."
                echo ""
                echo "You can:"
                echo "1. Download a pre-trained classifier from:"
                echo "   https://docs.qiime2.org/2024.10/data-resources/"
                echo ""
                echo "2. For 16S data, common classifiers:"
                echo "   - Silva 138 99% OTUs"
                echo "   - Greengenes2 2022.10"
                echo "   - GTDB"
                echo ""
                echo "3. Place the classifier in Data/reference_dbs/classifier.qza"
                echo ""
                echo "Example download command:"
                echo "  wget -O Data/reference_dbs/classifier.qza \\"
                echo "    https://data.qiime2.org/classifiers/sklearn-1.4.2/silva/silva-138-99-nb-classifier.qza"
                echo ""
                log "Skipping taxonomic classification. Please download a classifier and run step 7 again."
                create_checkpoint "step7_taxonomic_classification"
                return 0
                ;;
            4)
                log "Skipping taxonomic classification."
                create_checkpoint "step7_taxonomic_classification"
                return 0
                ;;
            *)
                log "Invalid choice. Skipping taxonomic classification."
                create_checkpoint "step7_taxonomic_classification"
                return 0
                ;;
        esac
    fi
    
    log "Using classifier: ${CLASSIFIER_PATH}"
    
    log "Running taxonomic classification (this may take a while)..."
    qiime feature-classifier classify-sklearn \
        --i-classifier "${CLASSIFIER_PATH}" \
        --i-reads "${REP_SEQS_PATH}" \
        --o-classification "${OUTPUT_DIR}/taxonomy.qza" || \
        error_exit "Taxonomic classification failed"
    
    log "Creating taxonomy visualization..."
    qiime metadata tabulate \
        --m-input-file "${OUTPUT_DIR}/taxonomy.qza" \
        --o-visualization "${OUTPUT_DIR}/taxonomy.qzv" || \
        error_exit "Failed to create taxonomy visualization"
    
    log "Creating taxonomic bar plots..."
    qiime taxa barplot \
        --i-table "${TABLE_PATH}" \
        --i-taxonomy "${OUTPUT_DIR}/taxonomy.qza" \
        --o-visualization "${OUTPUT_DIR}/taxa-bar-plots.qzv" || \
        error_exit "Failed to create taxonomic bar plots"
    
    # Optional: Export taxonomy table to TSV
    log "Exporting taxonomy table to TSV..."
    qiime tools export \
        --input-path "${OUTPUT_DIR}/taxonomy.qza" \
        --output-path "${OUTPUT_DIR}/exported-taxonomy" || \
        log "Warning: Failed to export taxonomy table"
    
    create_checkpoint "step7_taxonomic_classification"
    
    echo ""
    echo "=================================================="
    echo "Taxonomic Classification Complete!"
    echo "=================================================="
    echo ""
    echo "Generated files:"
    echo "1. ${OUTPUT_DIR}/taxonomy.qza"
    echo "   - Taxonomic assignments for each feature"
    echo ""
    echo "2. ${OUTPUT_DIR}/taxonomy.qzv"
    echo "   - Taxonomy table visualization"
    echo ""
    echo "3. ${OUTPUT_DIR}/taxa-bar-plots.qzv"
        # export feature table to biom and tsv
        echo "Exporting feature table to biom and tsv..."
        mkdir -p "${OUTPUT_DIR}/exported-table"
        qiime tools export \
            --input-path "${OUTPUT_DIR}/table.qza" \
            --output-path "${OUTPUT_DIR}/exported-table"
        if command -v biom &> /dev/null; then
            biom convert \
                -i "${OUTPUT_DIR}/exported-table/feature-table.biom" \
                -o "${OUTPUT_DIR}/exported-table/feature-table.tsv" \
                --to-tsv
        else
            echo "Warning: biom tool not found, cannot convert biom to tsv."
        fi
    echo "   - Interactive taxonomic bar plots (like your example image!)"
    echo "   - View at https://view.qiime2.org"
    echo ""
    if [[ -f "${OUTPUT_DIR}/exported-taxonomy/taxonomy.tsv" ]]; then
        echo "4. ${OUTPUT_DIR}/exported-taxonomy/taxonomy.tsv"
        echo "   - Taxonomy table in TSV format for further analysis"
        echo ""
    fi
    echo "The taxa-bar-plots.qzv file will show relative abundances"
    echo "of different taxa across your samples, just like your example!"
    echo ""
    
    create_checkpoint "step7_taxonomic_classification"
    pause_script "Taxonomic Classification" "Taxonomic classification complete. Open taxa-bar-plots.qzv to see your taxonomic composition!"
}

# Step 8: Phylogenetic Tree and Diversity Analysis
step8_diversity_analysis() {
    if check_checkpoint "step8_diversity_analysis"; then
        return 0
    fi
    
    log "Starting Step 8: Phylogenetic Tree Generation and Diversity Analysis"
    
    # Ensure environment is activated
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "${ENV_NAME}" || error_exit "Failed to activate QIIME2 environment"
    
    # Determine which mode we're in and set paths
    if [[ "${MODE}" == "denoise" ]]; then
        TABLE_PATH="Results/denoise_mode/table.qza"
        REP_SEQS_PATH="Results/denoise_mode/rep-seqs.qza"
        OUTPUT_DIR="Results/denoise_mode"
        
        # Check if decontaminated files exist and prefer them
        if [[ -f "Results/denoise_mode/table-clean.qza" ]]; then
            TABLE_PATH="Results/denoise_mode/table-clean.qza"
            REP_SEQS_PATH="Results/denoise_mode/rep-seqs-no-contam.qza"
            log "Using decontaminated table and sequences"
        elif [[ -f "Results/denoise_mode/table-no-contam.qza" ]]; then
            TABLE_PATH="Results/denoise_mode/table-no-contam.qza"
            REP_SEQS_PATH="Results/denoise_mode/rep-seqs-no-contam.qza"
            log "Using contaminant-filtered table and sequences"
        fi
    elif [[ "${MODE}" == "cluster" ]]; then
        TABLE_PATH="Results/cluster_mode/table.qza"
        REP_SEQS_PATH="Results/cluster_mode/rep-seqs.qza"
        OUTPUT_DIR="Results/cluster_mode"
        
        # Check if decontaminated files exist and prefer them
        if [[ -f "Results/cluster_mode/table-clean.qza" ]]; then
            TABLE_PATH="Results/cluster_mode/table-clean.qza"
            REP_SEQS_PATH="Results/cluster_mode/rep-seqs-no-contam.qza"
            log "Using decontaminated table and sequences"
        elif [[ -f "Results/cluster_mode/table-no-contam.qza" ]]; then
            TABLE_PATH="Results/cluster_mode/table-no-contam.qza"
            REP_SEQS_PATH="Results/cluster_mode/rep-seqs-no-contam.qza"
            log "Using contaminant-filtered table and sequences"
        fi
    else
        error_exit "Unknown mode: ${MODE}"
    fi
    
    # Check if required files exist
    if [[ ! -f "${TABLE_PATH}" ]] || [[ ! -f "${REP_SEQS_PATH}" ]]; then
        error_exit "Required files not found. Please run step 5 first."
    fi
    
    # Create diversity output directory
    mkdir -p "${OUTPUT_DIR}/diversity"
    
    echo ""
    echo "=================================================="
    echo "Phylogenetic Tree Generation"
    echo "=================================================="
    echo ""
    
    log "Generating phylogenetic tree (required for UniFrac metrics)..."
    log "This involves: alignment → masking → tree building → rooting"
    
    qiime phylogeny align-to-tree-mafft-fasttree \
        --i-sequences "${REP_SEQS_PATH}" \
        --o-alignment "${OUTPUT_DIR}/aligned-rep-seqs.qza" \
        --o-masked-alignment "${OUTPUT_DIR}/masked-aligned-rep-seqs.qza" \
        --o-tree "${OUTPUT_DIR}/unrooted-tree.qza" \
        --o-rooted-tree "${OUTPUT_DIR}/rooted-tree.qza" || \
        error_exit "Phylogenetic tree generation failed"
    
    log "Phylogenetic tree generated successfully"
    
    echo ""
    echo "=================================================="
    echo "Diversity Analysis"
    echo "=================================================="
    echo ""
    
    # Check for metadata file
    METADATA_PATH=""
    if [[ -f "Data/metadata/metadata.tsv" ]]; then
        METADATA_PATH="Data/metadata/metadata.tsv"
        log "Using metadata file: Data/metadata/metadata.tsv"
    elif [[ -f "metadata.tsv" ]]; then
        METADATA_PATH="metadata.tsv"
        log "Using metadata file: metadata.tsv (in working directory)"
    else
        log "Warning: No metadata.tsv file found. Diversity metrics will be calculated but group comparisons will not be available."
        log "To enable group comparisons, create a metadata file at Data/metadata/metadata.tsv"
    fi
    
    # Determine sampling depth
    echo ""
    echo "For diversity analysis, you need to choose a sampling depth (rarefaction depth)."
    echo "This should be based on the feature table summary (table.qzv)."
    echo ""
    echo "View ${OUTPUT_DIR}/table.qzv at https://view.qiime2.org to see:"
    echo "- Sampling depth per sample"
    echo "- Recommended rarefaction depth"
    echo ""
    
    if [[ "${NON_INTERACTIVE}" == "true" && -n "${SAMPLING_DEPTH:-}" ]]; then
        log "Non-interactive: using SAMPLING_DEPTH=${SAMPLING_DEPTH}"
    elif [[ "${NON_INTERACTIVE}" == "true" ]]; then
        SAMPLING_DEPTH=1000
        log "Non-interactive: default SAMPLING_DEPTH=${SAMPLING_DEPTH}"
    else
        read -p "Enter sampling depth for rarefaction (e.g., 1000): " SAMPLING_DEPTH
    fi
    
    SAMPLING_DEPTH=${SAMPLING_DEPTH:-1000}
    log "Using sampling depth: ${SAMPLING_DEPTH}"
    
    # Only calculate and export alpha diversity (shannon) and unweighted unifrac
    log "Calculating alpha diversity (shannon)..."
    qiime diversity alpha \
        --i-table "${TABLE_PATH}" \
        --p-metric shannon \
        --o-alpha-diversity "${OUTPUT_DIR}/diversity/shannon_vector.qza" || \
        log "Warning: Failed to calculate Shannon index"

    log "Exporting shannon vector..."
    qiime tools export \
        --input-path "${OUTPUT_DIR}/diversity/shannon_vector.qza" \
        --output-path "${OUTPUT_DIR}/diversity/shannon_exported" || \
        log "Warning: Failed to export Shannon vector"

    log "Visualizing shannon diversity..."
    qiime diversity alpha-group-significance \
        --i-alpha-diversity "${OUTPUT_DIR}/diversity/shannon_vector.qza" \
        --m-metadata-file "${METADATA_PATH}" \
        --o-visualization "${OUTPUT_DIR}/diversity/shannon_group_significance.qzv" || \
        log "Warning: Failed to visualize Shannon diversity"

    log "Calculating unweighted UniFrac distance..."
    qiime diversity beta-phylogenetic \
        --i-table "${TABLE_PATH}" \
        --i-phylogeny "${OUTPUT_DIR}/rooted-tree.qza" \
        --p-metric unweighted_unifrac \
        --o-distance-matrix "${OUTPUT_DIR}/diversity/unweighted_unifrac_distance_matrix.qza" || \
        log "Warning: Failed to calculate unweighted UniFrac distance"

    log "Exporting unweighted UniFrac distance matrix..."
    qiime tools export \
        --input-path "${OUTPUT_DIR}/diversity/unweighted_unifrac_distance_matrix.qza" \
        --output-path "${OUTPUT_DIR}/diversity/unweighted_unifrac_exported" || \
        log "Warning: Failed to export unweighted UniFrac distance matrix"

    log "PCoA and visualization for unweighted UniFrac..."
    qiime diversity pcoa \
        --i-distance-matrix "${OUTPUT_DIR}/diversity/unweighted_unifrac_distance_matrix.qza" \
        --o-pcoa "${OUTPUT_DIR}/diversity/unweighted_unifrac_pcoa.qza" || \
        log "Warning: Failed to run PCoA for unweighted UniFrac"

    qiime emperor plot \
        --i-pcoa "${OUTPUT_DIR}/diversity/unweighted_unifrac_pcoa.qza" \
        --m-metadata-file "${METADATA_PATH}" \
        --o-visualization "${OUTPUT_DIR}/diversity/unweighted_unifrac_emperor.qzv" || \
        log "Warning: Failed to generate Emperor plot for unweighted UniFrac"
    
    log "Exporting unweighted UniFrac PCoA results..."
    qiime tools export \
        --input-path "${OUTPUT_DIR}/diversity/unweighted_unifrac_pcoa.qza" \
        --output-path "${OUTPUT_DIR}/diversity/exported-unweighted_unifrac-pcoa" || \
        log "Warning: Failed to export unweighted UniFrac PCoA"
    
    # Create visualizations for additional metrics
    if [[ -n "${METADATA_PATH}" ]]; then
        if [[ -f "${OUTPUT_DIR}/diversity/chao1-vector.qza" ]]; then
            qiime diversity alpha-group-significance \
                --i-alpha-diversity "${OUTPUT_DIR}/diversity/chao1-vector.qza" \
                --m-metadata-file "${METADATA_PATH}" \
                --o-visualization "${OUTPUT_DIR}/diversity/chao1-significance.qzv" || \
                log "Warning: Failed to create Chao1 significance visualization"
        fi
        
        if [[ -f "${OUTPUT_DIR}/diversity/simpson-vector.qza" ]]; then
            qiime diversity alpha-group-significance \
                --i-alpha-diversity "${OUTPUT_DIR}/diversity/simpson-vector.qza" \
                --m-metadata-file "${METADATA_PATH}" \
                --o-visualization "${OUTPUT_DIR}/diversity/simpson-significance.qzv" || \
                log "Warning: Failed to create Simpson significance visualization"
        fi
    fi
    
    # Export alpha diversity metrics to TSV
    log "Exporting alpha diversity metrics to TSV format..."
    
    for metric in shannon observed_features evenness faith_pd chao1 simpson; do
        METRIC_FILE=""
        if [[ -f "${OUTPUT_DIR}/diversity/core-metrics/${metric}_vector.qza" ]]; then
            METRIC_FILE="${OUTPUT_DIR}/diversity/core-metrics/${metric}_vector.qza"
        elif [[ -f "${OUTPUT_DIR}/diversity/${metric}-vector.qza" ]]; then
            METRIC_FILE="${OUTPUT_DIR}/diversity/${metric}-vector.qza"
        fi
        
        if [[ -n "${METRIC_FILE}" ]]; then
            mkdir -p "${OUTPUT_DIR}/diversity/exported-${metric}"
            qiime tools export \
                --input-path "${METRIC_FILE}" \
                --output-path "${OUTPUT_DIR}/diversity/exported-${metric}" || \
                log "Warning: Failed to export ${metric}"
        fi
    done
    
    # Export beta diversity distance matrices
    log "Exporting beta diversity distance matrices to TSV format..."
    
    for metric in bray_curtis jaccard weighted_unifrac unweighted_unifrac; do
        if [[ -f "${OUTPUT_DIR}/diversity/core-metrics/${metric}_distance_matrix.qza" ]]; then
            mkdir -p "${OUTPUT_DIR}/diversity/exported-${metric}"
            qiime tools export \
                --input-path "${OUTPUT_DIR}/diversity/core-metrics/${metric}_distance_matrix.qza" \
                --output-path "${OUTPUT_DIR}/diversity/exported-${metric}" || \
                log "Warning: Failed to export ${metric} distance matrix"
        fi
    done
    
    # Export PCoA results to TSV format
    log "Exporting PCoA results to TSV format..."
    
    for metric in bray_curtis jaccard weighted_unifrac unweighted_unifrac; do
        if [[ -f "${OUTPUT_DIR}/diversity/core-metrics/${metric}_pcoa_results.qza" ]]; then
            mkdir -p "${OUTPUT_DIR}/diversity/exported-${metric}-pcoa"
            qiime tools export \
                --input-path "${OUTPUT_DIR}/diversity/core-metrics/${metric}_pcoa_results.qza" \
                --output-path "${OUTPUT_DIR}/diversity/exported-${metric}-pcoa" || \
                log "Warning: Failed to export ${metric} PCoA results"
        fi
    done
    
    # Export rarefied table
    log "Exporting rarefied feature table..."
    if [[ -f "${OUTPUT_DIR}/diversity/core-metrics/rarefied_table.qza" ]]; then
        mkdir -p "${OUTPUT_DIR}/diversity/exported-rarefied-table"
        qiime tools export \
            --input-path "${OUTPUT_DIR}/diversity/core-metrics/rarefied_table.qza" \
            --output-path "${OUTPUT_DIR}/diversity/exported-rarefied-table" || \
            log "Warning: Failed to export rarefied table"
        
        # Convert rarefied BIOM to TSV if biom tool is available
        if command -v biom &> /dev/null && [[ -f "${OUTPUT_DIR}/diversity/exported-rarefied-table/feature-table.biom" ]]; then
            biom convert \
                -i "${OUTPUT_DIR}/diversity/exported-rarefied-table/feature-table.biom" \
                -o "${OUTPUT_DIR}/diversity/exported-rarefied-table/feature-table.tsv" \
                --to-tsv || \
                log "Warning: Failed to convert rarefied table to TSV"
        fi
    fi
    
    # Export phylogenetic tree in Newick format
    log "Exporting phylogenetic trees..."
    
    if [[ -f "${OUTPUT_DIR}/rooted-tree.qza" ]]; then
        mkdir -p "${OUTPUT_DIR}/exported-tree"
        qiime tools export \
            --input-path "${OUTPUT_DIR}/rooted-tree.qza" \
            --output-path "${OUTPUT_DIR}/exported-tree" || \
            log "Warning: Failed to export rooted tree"
        
        # Rename to more descriptive filename
        if [[ -f "${OUTPUT_DIR}/exported-tree/tree.nwk" ]]; then
            mv "${OUTPUT_DIR}/exported-tree/tree.nwk" "${OUTPUT_DIR}/exported-tree/rooted-tree.nwk"
        fi
    fi
    
    if [[ -f "${OUTPUT_DIR}/unrooted-tree.qza" ]]; then
        mkdir -p "${OUTPUT_DIR}/exported-tree"
        qiime tools export \
            --input-path "${OUTPUT_DIR}/unrooted-tree.qza" \
            --output-path "${OUTPUT_DIR}/exported-tree/unrooted" || \
            log "Warning: Failed to export unrooted tree"
        
        # Move and rename
        if [[ -f "${OUTPUT_DIR}/exported-tree/unrooted/tree.nwk" ]]; then
            mv "${OUTPUT_DIR}/exported-tree/unrooted/tree.nwk" "${OUTPUT_DIR}/exported-tree/unrooted-tree.nwk"
            rmdir "${OUTPUT_DIR}/exported-tree/unrooted" 2>/dev/null || true
        fi
    fi
    
    create_checkpoint "step8_diversity_analysis"
    
    echo ""
    echo "=================================================="
    echo "Diversity Analysis Complete!"
    echo "=================================================="
    echo ""
    echo "Generated files:"
    echo ""
    echo "Phylogenetic Tree:"
    echo "1. ${OUTPUT_DIR}/rooted-tree.qza - Rooted phylogenetic tree"
    echo "2. ${OUTPUT_DIR}/unrooted-tree.qza - Unrooted phylogenetic tree"
    echo "3. ${OUTPUT_DIR}/exported-tree/rooted-tree.nwk - Tree in Newick format"
    echo ""
    echo "Alpha Diversity Metrics (QZA + TSV):"
    echo "1. Shannon Index - ${OUTPUT_DIR}/diversity/core-metrics/shannon_vector.qza"
    echo "   Exported: ${OUTPUT_DIR}/diversity/exported-shannon/alpha-diversity.tsv"
    echo "2. Observed Features (Richness) - ${OUTPUT_DIR}/diversity/core-metrics/observed_features_vector.qza"
    echo "   Exported: ${OUTPUT_DIR}/diversity/exported-observed_features/alpha-diversity.tsv"
    echo "3. Chao1 - ${OUTPUT_DIR}/diversity/chao1-vector.qza"
    echo "   Exported: ${OUTPUT_DIR}/diversity/exported-chao1/alpha-diversity.tsv"
    echo "4. Simpson Index - ${OUTPUT_DIR}/diversity/simpson-vector.qza"
    echo "   Exported: ${OUTPUT_DIR}/diversity/exported-simpson/alpha-diversity.tsv"
    echo "5. Faith's Phylogenetic Diversity - ${OUTPUT_DIR}/diversity/core-metrics/faith_pd_vector.qza"
    echo "   Exported: ${OUTPUT_DIR}/diversity/exported-faith_pd/alpha-diversity.tsv"
    echo "6. Pielou's Evenness - ${OUTPUT_DIR}/diversity/core-metrics/evenness_vector.qza"
    echo "   Exported: ${OUTPUT_DIR}/diversity/exported-evenness/alpha-diversity.tsv"
    echo ""
    echo "Beta Diversity Distance Matrices (QZA + TSV):"
    echo "1. Bray-Curtis - ${OUTPUT_DIR}/diversity/core-metrics/bray_curtis_distance_matrix.qza"
    echo "   Exported: ${OUTPUT_DIR}/diversity/exported-bray_curtis/distance-matrix.tsv"
    echo "2. Jaccard - ${OUTPUT_DIR}/diversity/core-metrics/jaccard_distance_matrix.qza"
    echo "   Exported: ${OUTPUT_DIR}/diversity/exported-jaccard/distance-matrix.tsv"
    echo "3. Weighted UniFrac - ${OUTPUT_DIR}/diversity/core-metrics/weighted_unifrac_distance_matrix.qza"
    echo "   Exported: ${OUTPUT_DIR}/diversity/exported-weighted_unifrac/distance-matrix.tsv"
    echo "4. Unweighted UniFrac - ${OUTPUT_DIR}/diversity/core-metrics/unweighted_unifrac_distance_matrix.qza"
    echo "   Exported: ${OUTPUT_DIR}/diversity/exported-unweighted_unifrac/distance-matrix.tsv"
    echo ""
    echo "PCoA Results (coordinates for plotting):"
    echo "- ${OUTPUT_DIR}/diversity/exported-bray_curtis-pcoa/ordination.txt"
    echo "- ${OUTPUT_DIR}/diversity/exported-jaccard-pcoa/ordination.txt"
    echo "- ${OUTPUT_DIR}/diversity/exported-weighted_unifrac-pcoa/ordination.txt"
    echo "- ${OUTPUT_DIR}/diversity/exported-unweighted_unifrac-pcoa/ordination.txt"
    echo ""
    echo "PCoA Visualizations (interactive plots):"
    echo "- ${OUTPUT_DIR}/diversity/core-metrics/bray_curtis_emperor.qzv"
    echo "- ${OUTPUT_DIR}/diversity/core-metrics/jaccard_emperor.qzv"
    echo "- ${OUTPUT_DIR}/diversity/core-metrics/weighted_unifrac_emperor.qzv"
    echo "- ${OUTPUT_DIR}/diversity/core-metrics/unweighted_unifrac_emperor.qzv"
    echo ""
    echo "Rarefied Feature Table:"
    echo "- ${OUTPUT_DIR}/diversity/exported-rarefied-table/feature-table.tsv"
    echo ""
    if [[ -n "${METADATA_PATH}" ]]; then
        echo "Alpha Diversity Significance Tests:"
        echo "- ${OUTPUT_DIR}/diversity/core-metrics/shannon-group-significance.qzv"
        echo "- ${OUTPUT_DIR}/diversity/core-metrics/faith-pd-group-significance.qzv"
        echo "- ${OUTPUT_DIR}/diversity/core-metrics/evenness-group-significance.qzv"
        echo "- ${OUTPUT_DIR}/diversity/chao1-significance.qzv"
        echo "- ${OUTPUT_DIR}/diversity/simpson-significance.qzv"
        echo ""
    fi
    echo "All .qzv files can be viewed at https://view.qiime2.org"
    echo "All TSV files are ready for analysis in R, Python, or Excel"
    echo ""
    
    # Export representative sequences for downstream analysis (e.g., PICRUSt2, phylogenetic placement)
    log "Exporting representative sequences for downstream analysis..."
    echo ""
    echo "Exporting representative sequences for downstream analysis..."
    
    if [[ -f "Results/denoise_mode/rep-seqs-no-contam.qza" ]]; then
        qiime tools export \
            --input-path Results/denoise_mode/rep-seqs-no-contam.qza \
            --output-path Results/denoise_mode/exported-rep-seqs || \
            log "Warning: Failed to export representative sequences"
        
        if [[ -f "Results/denoise_mode/exported-rep-seqs/dna-sequences.fasta" ]]; then
            echo "Representative sequences exported to: Results/denoise_mode/exported-rep-seqs/dna-sequences.fasta"
            log "Representative sequences exported successfully"
        fi
    else
        log "Warning: rep-seqs-no-contam.qza not found. Skipping export."
    fi
    echo ""
    
    pause_script "Diversity Analysis" "Diversity analysis complete. All results exported to TSV format!"
}

# Main execution function
main() {
    log "Starting QIIME2 Pipeline with Checkpoint System"
    log "Script directory: ${SCRIPT_DIR}"
    log "Project directory: ${PROJECT_DIR}"
    log "Checkpoint directory: ${CHECKPOINT_DIR}"
    log "Log file: ${LOG_FILE}"
    
    echo ""
    echo "=============================================="
    echo "QIIME2 Pipeline with Checkpoint System"
    echo "=============================================="
    echo ""
    
    # Ask user for environment name if not already set via command line
    if [[ -z "${ENV_NAME_SET:-}" ]]; then
        echo "Current conda environment name: ${ENV_NAME}"
        read -p "Enter your QIIME2 environment name (or press Enter to use '${ENV_NAME}'): " USER_ENV_NAME
        if [[ -n "${USER_ENV_NAME}" ]]; then
            ENV_NAME="${USER_ENV_NAME}"
            log "Using user-specified environment name: ${ENV_NAME}"
        else
            log "Using default environment name: ${ENV_NAME}"
        fi
        echo ""
    fi
    
    echo "This script will run the following steps:"
    echo "1. Environment Setup"
    echo "2. Import Data"
    echo "3. Visualize Demux Data"
    echo "4. Remove Primers/Adapters (optional)"
    echo "5. Denoising or Clustering (mode: denoise or cluster)"
    echo "6. Decontamination (optional)"
    echo "7. Taxonomic Classification"
    echo "8. Diversity Analysis (phylogenetic tree + alpha/beta diversity)"
    echo ""
    echo "You can:"
    echo "- Run all steps: Press Enter"
    echo "- Run from specific step: Type step number (1-8)"
    echo "- Exit: Ctrl+C"
    echo ""
    if [[ -z "${START_STEP:-}" ]]; then
        if [[ "${NON_INTERACTIVE}" == "true" ]]; then
            # In batch mode default to step 1 unless START_STEP exported
            START_STEP=${START_STEP:-1}
            log "Non-interactive: starting from step ${START_STEP}"
        else
            read -p "Enter step number to start from (or press Enter for all steps): " START_STEP
            START_STEP=${START_STEP:-1}
        fi
    else
        log "Starting from step ${START_STEP} (provided via CLI)"
    fi

    # Validate start step
    if ! [[ "${START_STEP}" =~ ^[1-8]$ ]]; then
        error_exit "Invalid step number. Please enter 1-8."
    fi
    
    # Run steps based on start step
    if [[ ${START_STEP} -le 1 ]]; then
        step1_environment_setup
    fi
    
    if [[ ${START_STEP} -le 2 ]]; then
        step2_import_data
    fi
    
    if [[ ${START_STEP} -le 3 ]]; then
        step3_visualize_demux
    fi
    
    if [[ ${START_STEP} -le 4 ]]; then
        step4_remove_primers
    fi
    
    if [[ ${START_STEP} -le 5 ]]; then
        if [[ "${MODE}" == "denoise" ]]; then
            step5_dada2_denoising
        elif [[ "${MODE}" == "cluster" ]]; then
            step5_cluster_from_demux
        else
            error_exit "Unknown mode: ${MODE}. Supported: denoise, cluster"
        fi
    fi
    
    if [[ ${START_STEP} -le 6 ]]; then
        step6_decontamination
    fi
    
    if [[ ${START_STEP} -le 7 ]]; then
        step7_taxonomic_classification
    fi
    
    if [[ ${START_STEP} -le 8 ]]; then
        step8_diversity_analysis
    fi

    log "Pipeline completed successfully!"
    echo ""
    echo "=============================================="
    echo "Pipeline completed successfully!"
    echo "=============================================="
    echo ""
    echo "Generated files:"
    echo "- Data/processed_data/demux-paired-end.qza (imported data)"
    echo "- Data/processed_data/demux-paired-end.qzv (quality visualization)"
    if [[ "${MODE}" == "denoise" ]]; then
        # Check for decontaminated files
        if [[ -f "Results/denoise_mode/table-clean.qza" ]]; then
            echo "- Results/denoise_mode/table-clean.qza (decontaminated feature table, controls removed)"
            echo "- Results/denoise_mode/rep-seqs-no-contam.qza (filtered representative sequences)"
        elif [[ -f "Results/denoise_mode/table-no-contam.qza" ]]; then
            echo "- Results/denoise_mode/table-no-contam.qza (contaminant-filtered feature table)"
            echo "- Results/denoise_mode/rep-seqs-no-contam.qza (filtered representative sequences)"
        else
            echo "- Results/denoise_mode/table.qza (feature table)"
            echo "- Results/denoise_mode/rep-seqs.qza (representative sequences)"
        fi
        echo "- Results/denoise_mode/denoising-stats.qza (DADA2 statistics)"
    elif [[ "${MODE}" == "cluster" ]]; then
        # Check for decontaminated files
        if [[ -f "Results/cluster_mode/table-clean.qza" ]]; then
            echo "- Results/cluster_mode/table-clean.qza (decontaminated feature table, controls removed)"
            echo "- Results/cluster_mode/rep-seqs-no-contam.qza (filtered representative sequences)"
        elif [[ -f "Results/cluster_mode/table-no-contam.qza" ]]; then
            echo "- Results/cluster_mode/table-no-contam.qza (contaminant-filtered feature table)"
            echo "- Results/cluster_mode/rep-seqs-no-contam.qza (filtered representative sequences)"
        else
            echo "- Results/cluster_mode/table.qza (feature table)"
            echo "- Results/cluster_mode/rep-seqs.qza (representative sequences)"
        fi
        echo "- Results/cluster_mode/joined.qza, filtered-seqs.qza (intermediate files)"
    fi
    if [[ -f "Results/${MODE}_mode/table-cr-97.qza" || -f "Results/${MODE}_mode/table-or-97.qza" ]]; then
        echo "- Results/${MODE}_mode/table-*-97.qza (additional clustered tables, if VSEARCH was used)"
    fi
    echo ""
    echo "Checkpoints saved in: ${CHECKPOINT_DIR}"
    echo "Log file: ${LOG_FILE}"
    echo ""
    echo "Next steps:"
    echo "1. Review the quality plots and visualizations at https://view.qiime2.org"
    echo "2. Analyze alpha and beta diversity results"
    echo "3. Perform additional statistical tests if needed"
    echo "4. Export data for publication or further analysis"
    echo ""
}

# Cleanup function for graceful exit
cleanup() {
    log "Script interrupted. Checkpoints preserved for resuming."
    exit 130
}

# Set up signal handlers
trap cleanup SIGINT SIGTERM

# Usage function
usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  -h, --help         Show this help message"
    echo "  -c, --clean        Remove all checkpoints and start fresh"
    echo "  -s, --status       Show status of completed steps"
    echo "  -r, --remove       Remove a specific checkpoint"
    echo "                     Example: -r step7_taxonomic_classification"
    echo "  -d, --delete-intermediate"
    echo "                     Delete intermediate outputs from step 4 and 5"
    echo "  -m, --mode         Pipeline mode: 'denoise' (DADA2) or 'cluster' (generate OTUs from demux). Default: denoise"
    echo "                     Example: -m cluster"
    echo "  -e, --env          Conda environment name to use. Default: qiime"
    echo "                     Example: -e qiime2-2024.5"
    echo ""
    echo "Examples:"
    echo "  $0                 Run the pipeline (resume from last checkpoint)"
    echo "  $0 --clean         Start the pipeline from scratch"
    echo "  $0 --status        Check which steps have been completed"
    echo "  $0 -r 6            Remove checkpoint for step 6"
    echo "  $0 -r 7            Remove checkpoint for step 7"
    echo "  $0 -d              Delete intermediate outputs and re-run step 4/5"
    echo "  $0 -e myqiime      Use conda environment named 'myqiime'"
    echo "  $0 -m cluster -e qiime2-amplicon  Use cluster mode with custom environment"
    echo ""
    echo "Available checkpoints to remove:"
    echo "  - step1_environment_setup"
    echo "  - step2_import_data"
    echo "  - step3_visualize_demux"
    echo "  - step4_remove_primers"
    echo "  - step5_dada2_denoising"
    echo "  - step5_cluster_from_demux"
    echo "  - step6_decontamination"
    echo "  - step7_taxonomic_classification"
    echo "  - step8_diversity_analysis"
}

# Status function
show_status() {
    echo "Pipeline Status:"
    echo "================"
    
    local steps=("step1_environment_setup" "step2_import_data" "step3_visualize_demux" "step4_remove_primers" "step5_dada2_denoising" "step5_cluster_from_demux" "step6_decontamination" "step7_taxonomic_classification" "step8_diversity_analysis")
    local step_names=("Step 1: Environment Setup" "Step 2: Import Data" "Step 3: Visualize Demux" "Step 4: Remove Primers" "Step 5: DADA2 Denoising" "Step 5: Cluster from Demux" "Step 6: Decontamination" "Step 7: Taxonomic Classification" "Step 8: Diversity Analysis")
    
    for i in "${!steps[@]}"; do
        local step="${steps[$i]}"
        local name="${step_names[$i]}"
        if [[ -f "${CHECKPOINT_DIR}/${step}.checkpoint" ]]; then
            local timestamp=$(cat "${CHECKPOINT_DIR}/${step}.checkpoint")
            echo "✓ ${name} (completed: ${timestamp})"
        else
            echo "✗ ${name} (not completed)"
        fi
    done
    echo ""
}

# Clean function
clean_checkpoints() {
    echo "Removing all checkpoints..."
    rm -rf "${CHECKPOINT_DIR}"
    rm -f "${LOG_FILE}"
    echo "Checkpoints and log file removed. Pipeline will start from the beginning."
}

# Remove specific checkpoint function
remove_checkpoint() {
    local step_name="$1"
    local checkpoint_file="${CHECKPOINT_DIR}/${step_name}.checkpoint"
    
    if [[ -f "${checkpoint_file}" ]]; then
        rm "${checkpoint_file}"
        echo "Checkpoint removed: ${step_name}"
        log "Checkpoint manually removed: ${step_name}"
    else
        echo "Checkpoint not found: ${step_name}"
        echo "Available checkpoints:"
        if [[ -d "${CHECKPOINT_DIR}" ]]; then
            for checkpoint in "${CHECKPOINT_DIR}"/*.checkpoint; do
                if [[ -f "${checkpoint}" ]]; then
                    basename "${checkpoint}" .checkpoint
                fi
            done
        else
            echo "  (none)"
        fi
    fi
}

# Remove intermediate outputs function
remove_intermediate_outputs() {
    echo ""
    echo "=================================================="
    echo "Remove Intermediate Outputs"
    echo "=================================================="
    echo ""
    echo "This will delete intermediate files from:"
    echo ""
    echo "Step 4 (Remove Primers):"
    echo "  - demux-trimmed.qza"
    echo "  - demux-trimmed.qzv"
    echo ""
    echo "Step 5 (DADA2 Denoising):"
    echo "  - table.qza"
    echo "  - rep-seqs.qza"
    echo "  - denoising-stats.qza"
    echo "  - base-transition-stats.qza"
    echo "  - rep-seqs.qzv"
    echo ""
    echo "And reset checkpoints for these steps."
    echo ""
    read -p "Are you sure you want to delete these files? (y/n): " CONFIRM_DELETE
    
    if [[ "${CONFIRM_DELETE}" =~ ^[Yy]$ ]]; then
        echo "Deleting intermediate outputs..."
        
        # Remove Step 4 files
        if [[ -f "demux-trimmed.qza" ]]; then
            rm "demux-trimmed.qza"
            echo "✓ Deleted: demux-trimmed.qza"
        fi
        if [[ -f "demux-trimmed.qzv" ]]; then
            rm "demux-trimmed.qzv"
            echo "✓ Deleted: demux-trimmed.qzv"
        fi
        
        # Remove Step 5 files
        if [[ -f "table.qza" ]]; then
            rm "table.qza"
            echo "✓ Deleted: table.qza"
        fi
        if [[ -f "rep-seqs.qza" ]]; then
            rm "rep-seqs.qza"
            echo "✓ Deleted: rep-seqs.qza"
        fi
        if [[ -f "denoising-stats.qza" ]]; then
            rm "denoising-stats.qza"
            echo "✓ Deleted: denoising-stats.qza"
        fi
        if [[ -f "base-transition-stats.qza" ]]; then
            rm "base-transition-stats.qza"
            echo "✓ Deleted: base-transition-stats.qza"
        fi
        if [[ -f "rep-seqs.qzv" ]]; then
            rm "rep-seqs.qzv"
            echo "✓ Deleted: rep-seqs.qzv"
        fi
        
        # Remove checkpoints
        echo "Removing checkpoints..."
        remove_checkpoint "step4_remove_primers"
        remove_checkpoint "step5_dada2_denoising"
        
        echo ""
        echo "Intermediate outputs removed successfully."
        echo "You can now re-run these steps with different parameters."
        log "Intermediate outputs from step4 and step5 removed by user"
    else
        echo "Cancelled. No files were deleted."
    fi
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help)
            usage
            exit 0
            ;;
        -c|--clean)
            clean_checkpoints
            exit 0
            ;;
        -s|--status)
            show_status
            exit 0
            ;;
        -d|--delete-intermediate)
            remove_intermediate_outputs
            exit 0
            ;;
        -r|--remove)
            if [[ -n "${2:-}" ]]; then
                arg="$2"
                case "${arg}" in
                    [1-7])
                        case "${arg}" in
                            1)
                                remove_checkpoint "step1_environment_setup"
                                ;;
                            2)
                                remove_checkpoint "step2_import_data"
                                ;;
                            3)
                                remove_checkpoint "step3_visualize_demux"
                                ;;
                            4)
                                remove_checkpoint "step4_remove_primers"
                                ;;
                            5)
                                remove_checkpoint "step5_dada2_denoising"
                                remove_checkpoint "step5_cluster_from_demux"
                                ;;
                            6)
                                remove_checkpoint "step6_decontamination"
                                ;;
                            7)
                                remove_checkpoint "step8_taxonomic_classification"
                                ;;
                        esac
                        exit 0
                        ;;
                    *)
                        # Accept explicit checkpoint name as before
                        remove_checkpoint "${arg}"
                        exit 0
                        ;;
                esac
            else
                echo "Error: --remove requires an argument (checkpoint number 1-7 or checkpoint name)"
                echo "Example: $0 --remove 7"
                exit 1
            fi
            ;;
        -b|--batch|--non-interactive)
            # Run in non-interactive batch mode (suitable for HPC). Prompts will be skipped
            # and environment variables will be used where available.
            NON_INTERACTIVE="true"
            shift
            ;;
        -m|--mode)
            if [[ -n "${2:-}" ]]; then
                MODE="$2"
                shift 2
            else
                echo "Error: --mode requires an argument (denoise|cluster)"
                exit 1
            fi
            ;;
        -e|--env)
            if [[ -n "${2:-}" ]]; then
                ENV_NAME="$2"
                ENV_NAME_SET="true"
                shift 2
            else
                echo "Error: --env requires an argument (environment name)"
                exit 1
            fi
            ;;
        # positional: start step (1-9)
        [1-7])
            START_STEP="$1"
            shift
            ;;
        *)
            echo "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

main

