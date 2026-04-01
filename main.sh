#!/usr/bin/env bash
set -euo pipefail

# ── Paths ─────────────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
DATA_DIR="${PROJECT_DIR}/Data"
PROCESSED_DIR="${DATA_DIR}/processed_data"
RESULTS_DIR="${PROJECT_DIR}/Results"
CHECKPOINT_DIR="${PROJECT_DIR}/Logs/checkpoints"
LOG_FILE="${PROJECT_DIR}/Logs/qiime2_pipeline.log"

# ── Configuration ─────────────────────────────────────────────────────────────
ENV_NAME="${ENV_NAME:-qiime2}"
MODE="${MODE:-denoise}"

mkdir -p "${CHECKPOINT_DIR}" "${PROCESSED_DIR}" "${RESULTS_DIR}" \
         "${DATA_DIR}/metadata" "${DATA_DIR}/reference_dbs"

# ── Core helpers ──────────────────────────────────────────────────────────────
log()         { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "${LOG_FILE}"; }
error_exit()  { log "ERROR: $1"; exit 1; }

create_checkpoint() {
    date '+%Y-%m-%d %H:%M:%S' > "${CHECKPOINT_DIR}/${1}.checkpoint"
    log "Checkpoint: ${1}"
}

check_checkpoint() {
    [[ -f "${CHECKPOINT_DIR}/${1}.checkpoint" ]] && { log "Already done: ${1}. Skipping."; return 0; }
    return 1
}

# Activate conda env (called once per step that needs qiime)
activate_env() {
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "${ENV_NAME}" || error_exit "Failed to activate ${ENV_NAME}"
}

# Wrap qiime tools export with error handling
qiime_export() {   # usage: qiime_export <input.qza> <output_dir>
    qiime tools export --input-path "$1" --output-path "$2" \
        || log "Warning: export failed for $1"
}

# Return the output dir for the current MODE
mode_outdir() {
    case "${MODE}" in
        denoise) echo "${RESULTS_DIR}/denoise_mode" ;;
        cluster) echo "${RESULTS_DIR}/cluster_mode" ;;
        *)       error_exit "Unknown MODE: ${MODE}" ;;
    esac
}

# Resolve the best available table + rep-seqs for downstream steps
# Sets globals: OUT_TABLE, OUT_REP_SEQS
resolve_outputs() {
    local base
    base="$(mode_outdir)"
    if [[ -f "${base}/table-clean.qza" ]]; then
        OUT_TABLE="${base}/table-clean.qza"
        OUT_REP_SEQS="${base}/rep-seqs-no-contam.qza"
        log "Using decontaminated + control-filtered table"
    elif [[ -f "${base}/table-no-contam.qza" ]]; then
        OUT_TABLE="${base}/table-no-contam.qza"
        OUT_REP_SEQS="${base}/rep-seqs-no-contam.qza"
        log "Using contaminant-filtered table"
    else
        OUT_TABLE="${base}/table.qza"
        OUT_REP_SEQS="${base}/rep-seqs.qza"
    fi
    [[ -f "${OUT_TABLE}" && -f "${OUT_REP_SEQS}" ]] \
        || error_exit "table/rep-seqs not found in ${base}. Run step 5 first."
}

# Pick the best demux input (trimmed preferred)
resolve_demux() {
    if   [[ -f "${PROCESSED_DIR}/demux-trimmed.qza"    ]]; then echo "${PROCESSED_DIR}/demux-trimmed.qza"
    elif [[ -f "${PROCESSED_DIR}/demux-paired-end.qza" ]]; then echo "${PROCESSED_DIR}/demux-paired-end.qza"
    else error_exit "Demux file not found. Run step 2 first."
    fi
}

# ── Step 1: Environment Setup ─────────────────────────────────────────────────
step1_environment_setup() {
    check_checkpoint "step1_environment_setup" && return 0
    log "Step 1: Environment Setup"

    command -v conda &>/dev/null || error_exit "Conda not found in PATH"
    source "$(conda info --base)/etc/profile.d/conda.sh"

    if conda env list | grep -q "^${ENV_NAME}[[:space:]]"; then
        log "Using existing environment: ${ENV_NAME}"
    else
        log "Creating environment '${ENV_NAME}' from QIIME2 distribution..."
        conda env create --name "${ENV_NAME}" \
            --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2025.10/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml \
            || error_exit "Failed to create environment '${ENV_NAME}'"
    fi

    conda activate "${ENV_NAME}" || error_exit "Failed to activate ${ENV_NAME}"
    qiime --version              || error_exit "QIIME2 not functional in ${ENV_NAME}"

    create_checkpoint "step1_environment_setup"
}

# ── Step 2: Import Data ───────────────────────────────────────────────────────
step2_import_data() {
    check_checkpoint "step2_import_data" && return 0
    log "Step 2: Import paired-end FASTQ files"
    activate_env

    local manifest="${DATA_DIR}/raw_data/manifest.tsv"
    [[ -f "${manifest}" ]] || error_exit "manifest.tsv not found: ${manifest}"

    qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path "${manifest}" \
        --output-path "${PROCESSED_DIR}/demux-paired-end.qza" \
        --input-format PairedEndFastqManifestPhred33V2 \
        || error_exit "Data import failed"

    create_checkpoint "step2_import_data"
}

# ── Step 3: Visualize Demux ───────────────────────────────────────────────────
step3_visualize_demux() {
    check_checkpoint "step3_visualize_demux" && return 0
    log "Step 3: Visualize demux"
    activate_env

    local demux="${PROCESSED_DIR}/demux-paired-end.qza"
    [[ -f "${demux}" ]] || error_exit "${demux} not found. Run step 2 first."

    qiime demux summarize \
        --i-data "${demux}" \
        --o-visualization "${PROCESSED_DIR}/demux-paired-end.qzv" \
        || error_exit "demux summarize failed"

    log "View ${PROCESSED_DIR}/demux-paired-end.qzv at https://view.qiime2.org"
    create_checkpoint "step3_visualize_demux"
}

# ── Step 4: Remove Primers ────────────────────────────────────────────────────
step4_remove_primers() {
    check_checkpoint "step4_remove_primers" && return 0
    log "Step 4: Primer removal (PRIMER_CHOICE=${PRIMER_CHOICE:-3})"

    local choice="${PRIMER_CHOICE:-3}"
    local primer_f primer_r

    case "${choice}" in
        1) primer_f="GTGCCAGCMGCCGCGGTAA"; primer_r="GGACTACHVGGGTWTCTAAT"
           log "515F/806R primers (V4)" ;;
        2) primer_f="${PRIMER_F:?PRIMER_F required}"; primer_r="${PRIMER_R:?PRIMER_R required}" ;;
        *) log "Skipping primer removal"; create_checkpoint "step4_remove_primers"; return 0 ;;
    esac

    activate_env

    qiime cutadapt trim-paired \
        --i-demultiplexed-sequences "${PROCESSED_DIR}/demux-paired-end.qza" \
        --p-front-f "${primer_f}" --p-front-r "${primer_r}" \
        --p-match-read-wildcards --p-match-adapter-wildcards \
        --o-trimmed-sequences "${PROCESSED_DIR}/demux-trimmed.qza" \
        --verbose || error_exit "Primer removal failed"

    qiime demux summarize \
        --i-data "${PROCESSED_DIR}/demux-trimmed.qza" \
        --o-visualization "${PROCESSED_DIR}/demux-trimmed.qzv" \
        || log "Warning: trimmed demux visualization failed"

    create_checkpoint "step4_remove_primers"
}

# ── Step 5a: DADA2 Denoising ──────────────────────────────────────────────────
step5_dada2_denoising() {
    check_checkpoint "step5_dada2_denoising" && return 0
    log "Step 5 (denoise): DADA2"
    activate_env

    local demux; demux="$(resolve_demux)"
    local out="${RESULTS_DIR}/denoise_mode"
    local trunc_f="${TRUNC_LEN_F:-250}" trunc_r="${TRUNC_LEN_R:-250}" threads="${N_THREADS:-0}"
    log "DADA2: trunc_f=${trunc_f} trunc_r=${trunc_r} threads=${threads}"
    mkdir -p "${out}"

    qiime dada2 denoise-paired \
        --i-demultiplexed-seqs "${demux}" \
        --p-trunc-len-f "${trunc_f}" --p-trunc-len-r "${trunc_r}" \
        --p-trim-left-f 0 --p-trim-left-r 0 \
        --p-n-threads "${threads}" \
        --o-table "${out}/table.qza" \
        --o-representative-sequences "${out}/rep-seqs.qza" \
        --o-denoising-stats "${out}/denoising-stats.qza" \
        --o-base-transition-stats "${out}/base-transition-stats.qza" \
        || error_exit "DADA2 failed"

    qiime metadata tabulate \
        --m-input-file "${out}/denoising-stats.qza" \
        --o-visualization "${out}/denoising-stats.qzv" \
        || log "Warning: denoising-stats.qzv failed"
    qiime feature-table summarize \
        --i-table "${out}/table.qza" \
        --o-visualization "${out}/table.qzv" \
        || log "Warning: table.qzv failed"
    qiime feature-table tabulate-seqs \
        --i-data "${out}/rep-seqs.qza" \
        --o-visualization "${out}/rep-seqs.qzv" \
        || log "Warning: rep-seqs.qzv failed"

    create_checkpoint "step5_dada2_denoising"
}

# ── Step 5b: OTU Clustering ───────────────────────────────────────────────────
step5_cluster_from_demux() {
    check_checkpoint "step5_cluster_from_demux" && return 0
    log "Step 5 (cluster): OTU clustering"
    activate_env

    local demux; demux="$(resolve_demux)"
    local out="${RESULTS_DIR}/cluster_mode"
    mkdir -p "${out}"

    qiime vsearch merge-pairs \
        --i-demultiplexed-seqs "${demux}" \
        --o-merged-sequences "${out}/joined.qza" \
        --o-unmerged-sequences "${out}/unmerged.qza" \
        || error_exit "merge-pairs failed"

    qiime quality-filter q-score \
        --i-demux "${out}/joined.qza" \
        --o-filtered-sequences "${out}/filtered-seqs.qza" \
        --o-filter-stats "${out}/filtered-stats.qza" \
        || error_exit "quality-filter failed"

    qiime vsearch dereplicate-sequences \
        --i-sequences "${out}/filtered-seqs.qza" \
        --o-dereplicated-table "${out}/derep-table.qza" \
        --o-dereplicated-sequences "${out}/derep-seqs.qza" \
        || error_exit "dereplicate-sequences failed"

    qiime vsearch uchime-denovo \
        --i-table "${out}/derep-table.qza" \
        --i-sequences "${out}/derep-seqs.qza" \
        --o-chimeras "${out}/chimeras.qza" \
        --o-nonchimeras "${out}/rep-seqs.qza" \
        --o-stats "${out}/chimera-stats.qza" \
        || error_exit "chimera removal failed"

    qiime feature-table filter-features \
        --i-table "${out}/derep-table.qza" \
        --m-metadata-file "${out}/rep-seqs.qza" \
        --o-filtered-table "${out}/table.qza" \
        || error_exit "table filtering failed"

    case "${CLUSTER_CHOICE:-1}" in
        1)
            log "De novo clustering at 97%..."
            qiime vsearch cluster-features-de-novo \
                --i-table "${out}/table.qza" --i-sequences "${out}/rep-seqs.qza" \
                --p-perc-identity 0.97 \
                --o-clustered-table "${out}/table-dn-97.qza" \
                --o-clustered-sequences "${out}/rep-seqs-dn-97.qza" \
                || error_exit "de-novo clustering failed"
            mv "${out}/table-dn-97.qza" "${out}/table.qza"
            mv "${out}/rep-seqs-dn-97.qza" "${out}/rep-seqs.qza"
            ;;
        2)
            local ref="${DATA_DIR}/reference_dbs/silva_97_otus.qza"
            [[ -f "${ref}" ]] || error_exit "Reference DB not found: ${ref}"
            qiime vsearch cluster-features-closed-reference \
                --i-table "${out}/table.qza" --i-sequences "${out}/rep-seqs.qza" \
                --i-reference-sequences "${ref}" --p-perc-identity 0.97 \
                --o-clustered-table "${out}/table-cr-97.qza" \
                --o-clustered-sequences "${out}/rep-seqs-cr-97.qza" \
                --o-unmatched-sequences "${out}/unmatched.qza" \
                || error_exit "closed-reference clustering failed"
            mv "${out}/table-cr-97.qza" "${out}/table.qza"
            mv "${out}/rep-seqs-cr-97.qza" "${out}/rep-seqs.qza"
            ;;
        3) log "Skipping additional clustering; using dereplicated features." ;;
    esac

    create_checkpoint "step5_cluster_from_demux"
}

# ── Step 6: Decontamination ───────────────────────────────────────────────────
step6_decontamination() {
    check_checkpoint "step6_decontamination" && return 0
    log "Step 6: Decontamination"

    if [[ "${DECONTAMINATION_CHOICE:-2}" != "1" ]]; then
        log "Skipping decontamination"
        create_checkpoint "step6_decontamination"; return 0
    fi

    activate_env
    local out; out="$(mode_outdir)"
    local table="${out}/table.qza"
    [[ -f "${table}" ]] || error_exit "${table} not found. Run step 5 first."

    # Locate metadata
    local meta=""
    if   [[ -n "${DECONTAM_METADATA:-}" && -f "${DECONTAM_METADATA}" ]]; then meta="${DECONTAM_METADATA}"
    elif [[ -f "${DATA_DIR}/metadata/decontam-metadata.tsv" ]];           then meta="${DATA_DIR}/metadata/decontam-metadata.tsv"
    elif [[ -f "${DATA_DIR}/metadata/metadata.tsv" ]];                    then meta="${DATA_DIR}/metadata/metadata.tsv"
    else log "No decontam metadata found. Skipping."; create_checkpoint "step6_decontamination"; return 0
    fi
    log "Decontam metadata: ${meta}"

    qiime quality-control decontam-identify \
        --i-table "${table}" --m-metadata-file "${meta}" \
        --p-method prevalence \
        --p-prev-control-column control_status --p-prev-control-indicator control \
        --o-decontam-scores "${out}/decontam-scores.qza" \
        || error_exit "decontam-identify failed"

    qiime quality-control decontam-score-viz \
        --i-decontam-scores "${out}/decontam-scores.qza" --i-table "${table}" \
        --o-visualization "${out}/decontam-scores.qzv" \
        || error_exit "decontam-score-viz failed"

    local export_dir="${out}/decontam_scores_export"
    qiime_export "${out}/decontam-scores.qza" "${export_dir}"

    local scores_file; scores_file=$(find "${export_dir}" -name "*.tsv" | head -1)
    [[ -n "${scores_file}" ]] || error_exit "No TSV found in ${export_dir}"

    local threshold="${DECONTAM_THRESHOLD:-0.5}"
    python3 - "${scores_file}" "${threshold}" "${export_dir}" << 'PYTHON'
import csv, sys
scores, threshold, outdir = sys.argv[1], float(sys.argv[2]), sys.argv[3]
ids = []
with open(scores) as f:
    for row in csv.DictReader(f, delimiter='\t'):
        for col in ['p', 'score', 'Score', 'p-value', 'p_value']:
            if col in row:
                try:
                    if float(row[col]) > threshold:
                        for id_col in ['#OTU ID', '#OTU id', 'feature-id', 'feature_id', 'id']:
                            if id_col in row: ids.append(row[id_col]); break
                    break
                except ValueError: pass
with open(f"{outdir}/contaminant-ids.txt", 'w') as f:
    f.write('feature-id\n')
    f.writelines(f'{i}\n' for i in ids)
print(f"Contaminants identified: {len(ids)} (threshold={threshold})")
PYTHON

    local contam_ids="${export_dir}/contaminant-ids.txt"
    [[ -f "${contam_ids}" ]] || error_exit "contaminant-ids.txt not created"
    local n; n=$(tail -n +2 "${contam_ids}" | wc -l)
    log "Contaminants: ${n}"

    local final_table="${table}" final_rep_seqs="${out}/rep-seqs.qza"

    if [[ "${n}" -gt 0 ]]; then
        qiime feature-table filter-features \
            --i-table "${table}" --m-metadata-file "${contam_ids}" --p-exclude-ids \
            --o-filtered-table "${out}/table-no-contam.qza" \
            || error_exit "Feature filtering failed"
        qiime feature-table filter-seqs \
            --i-data "${out}/rep-seqs.qza" --i-table "${out}/table-no-contam.qza" \
            --o-filtered-data "${out}/rep-seqs-no-contam.qza" \
            || error_exit "Seq filtering failed"
        final_table="${out}/table-no-contam.qza"
        final_rep_seqs="${out}/rep-seqs-no-contam.qza"
    fi

    if [[ "${REMOVE_CONTROLS:-y}" =~ ^[Yy]$ ]]; then
        qiime feature-table filter-samples \
            --i-table "${final_table}" --m-metadata-file "${meta}" \
            --p-where "control_status!='control'" \
            --o-filtered-table "${out}/table-clean.qza" \
            || error_exit "Sample filtering failed"
        log "Negative controls removed -> ${out}/table-clean.qza"
    fi

    create_checkpoint "step6_decontamination"
}

# ── Step 7: Taxonomic Classification ─────────────────────────────────────────
step7_taxonomic_classification() {
    check_checkpoint "step7_taxonomic_classification" && return 0
    log "Step 7: Taxonomic Classification"
    activate_env

    resolve_outputs   # sets OUT_TABLE, OUT_REP_SEQS
    local out; out="$(mode_outdir)"

    # Auto-detect classifier
    local classifier=""
    for c in \
        "${DATA_DIR}/reference_dbs/silva-138-99-nb-classifier.qza" \
        "${DATA_DIR}/reference_dbs/classifier.qza" \
        "${DATA_DIR}/reference_dbs/silva-v4-classifier.qza" \
        "${DATA_DIR}/reference_dbs/custom-classifier.qza"; do
        [[ -f "${c}" ]] && { classifier="${c}"; break; }
    done

    if [[ -z "${classifier}" ]]; then
        case "${CLASSIFIER_CHOICE:-skip}" in
            1)
                classifier="${USER_CLASSIFIER:?USER_CLASSIFIER required}"
                [[ -f "${classifier}" ]] || error_exit "Classifier not found: ${classifier}"
                ;;
            2)
                log "Building classifier via RESCRIPt..."
                local db="${DATA_DIR}/reference_dbs"
                qiime rescript get-silva-data \
                    --p-version '138.1' --p-target 'SSURef_NR99' \
                    --o-silva-sequences "${db}/silva-seqs.qza" \
                    --o-silva-taxonomy "${db}/silva-tax.qza" \
                    || error_exit "SILVA download failed"
                qiime rescript cull-seqs \
                    --i-sequences "${db}/silva-seqs.qza" \
                    --o-clean-sequences "${db}/silva-seqs-clean.qza" \
                    || error_exit "cull-seqs failed"
                qiime rescript filter-seqs-length-by-taxon \
                    --i-sequences "${db}/silva-seqs-clean.qza" \
                    --i-taxonomy "${db}/silva-tax.qza" --p-min-lens 1000 \
                    --o-filtered-seqs "${db}/silva-seqs-filtered.qza" \
                    || error_exit "filter-seqs failed"
                local ref="${db}/silva-seqs-filtered.qza"
                if [[ "${EXTRACT_PRIMERS:-n}" =~ ^[Yy]$ ]]; then
                    qiime feature-classifier extract-reads \
                        --i-sequences "${ref}" \
                        --p-f-primer "${PRIMER_F:?}" --p-r-primer "${PRIMER_R:?}" \
                        --o-reads "${db}/silva-extracted.qza" \
                        || error_exit "extract-reads failed"
                    ref="${db}/silva-extracted.qza"
                fi
                qiime feature-classifier fit-classifier-naive-bayes \
                    --i-reference-reads "${ref}" \
                    --i-reference-taxonomy "${db}/silva-tax.qza" \
                    --o-classifier "${db}/custom-classifier.qza" \
                    || error_exit "Classifier training failed"
                classifier="${db}/custom-classifier.qza"
                ;;
            *)
                log "No classifier. Skipping step 7. Place a classifier in ${DATA_DIR}/reference_dbs/ and re-run."
                create_checkpoint "step7_taxonomic_classification"; return 0 ;;
        esac
    fi
    log "Classifier: ${classifier}"

    qiime feature-classifier classify-sklearn \
        --i-classifier "${classifier}" --i-reads "${OUT_REP_SEQS}" \
        --o-classification "${out}/taxonomy.qza" \
        || error_exit "Classification failed"

    qiime metadata tabulate \
        --m-input-file "${out}/taxonomy.qza" \
        --o-visualization "${out}/taxonomy.qzv" \
        || log "Warning: taxonomy.qzv failed"

    qiime taxa barplot \
        --i-table "${OUT_TABLE}" --i-taxonomy "${out}/taxonomy.qza" \
        --o-visualization "${out}/taxa-bar-plots.qzv" \
        || error_exit "taxa barplot failed"

    qiime_export "${out}/taxonomy.qza" "${out}/exported-taxonomy"

    # Export feature table -> biom + tsv
    mkdir -p "${out}/exported-table"
    qiime_export "${OUT_TABLE}" "${out}/exported-table"
    if command -v biom &>/dev/null && [[ -f "${out}/exported-table/feature-table.biom" ]]; then
        biom convert \
            -i "${out}/exported-table/feature-table.biom" \
            -o "${out}/exported-table/feature-table.tsv" --to-tsv \
            || log "Warning: biom conversion failed"
    fi

    create_checkpoint "step7_taxonomic_classification"
    log "Step 7 complete. View ${out}/taxa-bar-plots.qzv at https://view.qiime2.org"
}

# ── Step 8: Phylogenetic Tree & Diversity ─────────────────────────────────────
step8_diversity_analysis() {
    check_checkpoint "step8_diversity_analysis" && return 0
    log "Step 8: Phylogenetic Tree & Diversity"
    activate_env

    resolve_outputs   # sets OUT_TABLE, OUT_REP_SEQS
    local out; out="$(mode_outdir)"
    mkdir -p "${out}/diversity"

    qiime phylogeny align-to-tree-mafft-fasttree \
        --i-sequences "${OUT_REP_SEQS}" \
        --o-alignment "${out}/aligned-rep-seqs.qza" \
        --o-masked-alignment "${out}/masked-aligned-rep-seqs.qza" \
        --o-tree "${out}/unrooted-tree.qza" \
        --o-rooted-tree "${out}/rooted-tree.qza" \
        || error_exit "Tree generation failed"

    local meta=""
    [[ -f "${DATA_DIR}/metadata/metadata.tsv" ]] && meta="${DATA_DIR}/metadata/metadata.tsv"
    [[ -z "${meta}" ]] && log "Warning: No metadata.tsv; group comparisons unavailable."

    # Alpha: Shannon
    qiime diversity alpha \
        --i-table "${OUT_TABLE}" --p-metric shannon \
        --o-alpha-diversity "${out}/diversity/shannon_vector.qza" \
        || log "Warning: Shannon failed"
    qiime_export "${out}/diversity/shannon_vector.qza" "${out}/diversity/exported-shannon"

    if [[ -n "${meta}" ]]; then
        qiime diversity alpha-group-significance \
            --i-alpha-diversity "${out}/diversity/shannon_vector.qza" \
            --m-metadata-file "${meta}" \
            --o-visualization "${out}/diversity/shannon_group_significance.qzv" \
            || log "Warning: Shannon group-significance failed"
    fi

    # Beta: unweighted UniFrac
    local unifrac="${out}/diversity/unweighted_unifrac_distance_matrix.qza"
    qiime diversity beta-phylogenetic \
        --i-table "${OUT_TABLE}" --i-phylogeny "${out}/rooted-tree.qza" \
        --p-metric unweighted_unifrac \
        --o-distance-matrix "${unifrac}" \
        || log "Warning: UniFrac failed"
    qiime_export "${unifrac}" "${out}/diversity/exported-unweighted_unifrac"

    local pcoa="${out}/diversity/unweighted_unifrac_pcoa.qza"
    qiime diversity pcoa --i-distance-matrix "${unifrac}" --o-pcoa "${pcoa}" \
        || log "Warning: PCoA failed"
    qiime_export "${pcoa}" "${out}/diversity/exported-unweighted_unifrac_pcoa"

    if [[ -n "${meta}" ]]; then
        qiime emperor plot \
            --i-pcoa "${pcoa}" --m-metadata-file "${meta}" \
            --o-visualization "${out}/diversity/unweighted_unifrac_emperor.qzv" \
            || log "Warning: Emperor plot failed"
    fi

    # Export tree (Newick)
    mkdir -p "${out}/exported-tree"
    qiime_export "${out}/rooted-tree.qza" "${out}/exported-tree"
    [[ -f "${out}/exported-tree/tree.nwk" ]] \
        && mv "${out}/exported-tree/tree.nwk" "${out}/exported-tree/rooted-tree.nwk"

    create_checkpoint "step8_diversity_analysis"
    log "Step 8 complete. Outputs in ${out}/diversity/"
}

# ── Utility commands ──────────────────────────────────────────────────────────
show_status() {
    echo "Pipeline Status  (PROJECT_DIR=${PROJECT_DIR})"
    echo "================================================"
    local steps=("step1_environment_setup" "step2_import_data" "step3_visualize_demux"
                 "step4_remove_primers" "step5_dada2_denoising" "step5_cluster_from_demux"
                 "step6_decontamination" "step7_taxonomic_classification" "step8_diversity_analysis")
    for s in "${steps[@]}"; do
        local f="${CHECKPOINT_DIR}/${s}.checkpoint"
        printf "  %s %s\n" "$([[ -f $f ]] && echo checkmark || echo x)" \
               "$s$([[ -f $f ]] && echo "  ($(cat "$f"))" || echo '')"
    done
}

clean_checkpoints() {
    rm -rf "${CHECKPOINT_DIR}" "${LOG_FILE}"
    echo "Checkpoints and log removed."
}

remove_checkpoint() {
    local f="${CHECKPOINT_DIR}/${1}.checkpoint"
    if [[ -f "${f}" ]]; then rm "${f}"; echo "Removed: ${1}"
    else echo "Not found: ${1}"; fi
}

remove_intermediate_outputs() {
    log "Removing step 4/5 intermediate outputs..."
    local files=(
        "${PROCESSED_DIR}/demux-trimmed.qza"
        "${PROCESSED_DIR}/demux-trimmed.qzv"
        "${RESULTS_DIR}/denoise_mode/table.qza"
        "${RESULTS_DIR}/denoise_mode/rep-seqs.qza"
        "${RESULTS_DIR}/denoise_mode/denoising-stats.qza"
        "${RESULTS_DIR}/denoise_mode/base-transition-stats.qza"
        "${RESULTS_DIR}/denoise_mode/rep-seqs.qzv"
    )
    for f in "${files[@]}"; do [[ -f "${f}" ]] && { rm "${f}"; log "Deleted: ${f}"; }; done
    remove_checkpoint "step4_remove_primers"
    remove_checkpoint "step5_dada2_denoising"
}

usage() {
    cat <<EOF
Usage: $0 [OPTIONS] [START_STEP]

Options:
  -h, --help                Show this help
  -c, --clean               Remove all checkpoints (start fresh)
  -s, --status              Show step completion status
  -d, --delete-intermediate Delete step 4/5 intermediate outputs
  -r, --remove <N|name>     Remove specific checkpoint (1-8 or full name)
  -m, --mode <denoise|cluster>  Pipeline mode (default: denoise)
  -e, --env  <name>         Conda environment name (default: qiime2)

Positional:
  START_STEP  Step number 1-8 to resume from

Examples:
  bash code/main.sh                   # run all steps
  bash code/main.sh 5                 # resume from step 5
  bash code/main.sh -m cluster        # OTU clustering mode
  bash code/main.sh -r 7              # re-run step 7
  TRUNC_LEN_F=230 TRUNC_LEN_R=200 N_THREADS=16 bash code/main.sh
EOF
}

# ── Signal handler ────────────────────────────────────────────────────────────
trap 'log "Interrupted. Checkpoints preserved."; exit 130' SIGINT SIGTERM

# ── Main ──────────────────────────────────────────────────────────────────────
main() {
    log "====== QIIME2 16S Pipeline  MODE=${MODE}  ENV=${ENV_NAME} ======"
    local start="${START_STEP:-1}"
    [[ "${start}" =~ ^[1-8]$ ]] || error_exit "Invalid START_STEP: ${start}"

    [[ ${start} -le 1 ]] && step1_environment_setup
    [[ ${start} -le 2 ]] && step2_import_data
    [[ ${start} -le 3 ]] && step3_visualize_demux
    [[ ${start} -le 4 ]] && step4_remove_primers

    if [[ ${start} -le 5 ]]; then
        [[ "${MODE}" == "denoise" ]] && step5_dada2_denoising \
            || step5_cluster_from_demux
    fi

    [[ ${start} -le 6 ]] && step6_decontamination
    [[ ${start} -le 7 ]] && step7_taxonomic_classification
    [[ ${start} -le 8 ]] && step8_diversity_analysis

    log "====== Pipeline complete. Results: ${RESULTS_DIR}/${MODE}_mode/ ======"
}

# ── Argument parsing ──────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help)            usage; exit 0 ;;
        -c|--clean)           clean_checkpoints; exit 0 ;;
        -s|--status)          show_status; exit 0 ;;
        -d|--delete-intermediate) remove_intermediate_outputs; exit 0 ;;
        -r|--remove)
            [[ -n "${2:-}" ]] || { echo "Error: --remove needs an argument"; exit 1; }
            case "$2" in
                1) remove_checkpoint "step1_environment_setup" ;;
                2) remove_checkpoint "step2_import_data" ;;
                3) remove_checkpoint "step3_visualize_demux" ;;
                4) remove_checkpoint "step4_remove_primers" ;;
                5) remove_checkpoint "step5_dada2_denoising"
                   remove_checkpoint "step5_cluster_from_demux" ;;
                6) remove_checkpoint "step6_decontamination" ;;
                7) remove_checkpoint "step7_taxonomic_classification" ;;
                8) remove_checkpoint "step8_diversity_analysis" ;;
                *) remove_checkpoint "$2" ;;
            esac
            shift 2; exit 0 ;;
        -m|--mode) MODE="$2"; shift 2 ;;
        -e|--env)  ENV_NAME="$2"; shift 2 ;;
        [1-8]) START_STEP="$1"; shift ;;
        *) echo "Unknown option: $1"; usage; exit 1 ;;
    esac
done

main
