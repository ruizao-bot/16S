# PICRUSt2 Targeted Abundance Extraction
library(readr)
library(tidyverse)
library(dplyr)

# ---------------------------------------------------------
# 1. Setup & Data Loading
# ---------------------------------------------------------
# Set base directory
base_dir <- "/Users/jiayi/Desktop/metagenomic_pipeline/QIIME"
setwd(base_dir)

# Load metadata and pathway abundance
abundance_file <- file.path(base_dir, "picrust2_output_final/pathways_out/path_abun_unstrat.tsv")
metadata_file <- file.path(base_dir, "Data/metadata/metadata.tsv")

# Read data
metadata <- read_delim(metadata_file, delim = "\t", trim_ws = TRUE)
abundance <- read_delim(abundance_file, delim = "\t", trim_ws = TRUE) %>% 
  column_to_rownames(var = colnames(.)[1])

# Filter samples to match metadata
clean_metadata <- metadata %>% filter(control_status == "sample")
common_samples <- intersect(colnames(abundance), clean_metadata$`sample-id`)
abundance <- abundance %>% dplyr::select(dplyr::all_of(common_samples))

cat("Loaded data: ", ncol(abundance), "samples.\n")

# ---------------------------------------------------------
# 2. Define Targets
# ---------------------------------------------------------
methane_pathways <- list(
  methanogenesis = c("METHANOGENESIS-PWY", "PWY-7084", "PWY-5173", "PWY-7255", "METHGLYUT-PWY"),
  methanotroph = c("METHANE-OXIDATION-PWY", "PWY-7357", "FORMALDEHYDE-OXIDATION-PWY", "PWY-7013")
)

methane_enzymes <- list(
  methanogenesis = list(
    EC = c("2.8.4.1", "1.12.98.1", "2.1.1.86", "1.2.7.4", "2.3.1.101"),
    KO = c("K00399", "K00401", "K00402", "K00577", "K00578", "K00579", "K00580", "K00581", "K00582", "K00583", "K00584", "K00200", "K00201", "K00202", "K00672", "K01499")
  ),
  methanotroph = list(
    EC = c("1.14.13.25", "1.14.18.3", "1.1.1.1", "1.2.1.46", "1.1.2.7"),
    KO = c("K10944", "K10945", "K10946", "K16157", "K16158", "K16159", "K14028", "K14029", "K00121", "K00122", "K08691")
  )
)

# ---------------------------------------------------------
# 3. Annotation Database (for descriptions)
# ---------------------------------------------------------
ensure_file <- function(url, file_path) {
  if (!file.exists(file_path)) {
    cat("Downloading:", file_path, "...\n")
    tryCatch(download.file(url, file_path, method = "auto"), error = function(e) cat("Download failed:", url, "\n"))
  }
}

ensure_file("https://raw.githubusercontent.com/picrust/picrust2/master/picrust2/default_files/description_mapfiles/ko_info.tsv.gz", "picrust2_output_final/ko_info.tsv.gz")
ensure_file("https://raw.githubusercontent.com/picrust/picrust2/master/picrust2/default_files/description_mapfiles/ec_level4_info.tsv.gz", "picrust2_output_final/ec_level4_info.tsv.gz")

load_desc <- function(path) if(file.exists(path)) read_tsv(path, col_names = c("ID", "Description"), col_types = "cc") else data.frame(ID=character(), Description=character())
enzyme_descriptions <- bind_rows(load_desc("picrust2_output_final/ko_info.tsv.gz"), load_desc("picrust2_output_final/ec_level4_info.tsv.gz")) %>%
  mutate(ID = gsub("^EC:", "", ID), ID = gsub("^ko:", "", ID))

cat("Loaded", nrow(enzyme_descriptions), "enzyme descriptions (KO + EC)\n")

# Manual additions
manual_desc <- data.frame(
  ID = c("2.1.1.21"),
  Description = c("N5-methyltetrahydromethanopterin:coenzyme M methyltransferase (transferred to 2.1.1.86)"),
  stringsAsFactors = FALSE
)
enzyme_descriptions <- bind_rows(enzyme_descriptions, manual_desc)

# ---------------------------------------------------------
# 4. Processing Functions
# ---------------------------------------------------------
save_abundance <- function(data, features, type_name, id_col, output_prefix) {
  # Filter data
  found_data <- data %>% filter(!!sym(id_col) %in% features)
  
  if (nrow(found_data) > 0) {
    # Add descriptions if available
    found_data_desc <- found_data %>%
      left_join(enzyme_descriptions, by = setNames("ID", id_col)) %>%
      dplyr::select(dplyr::all_of(id_col), Description, everything())
    
    # Report matching statistics
    n_matched <- sum(!is.na(found_data_desc$Description))
    cat("  ", n_matched, "/", nrow(found_data_desc), id_col, "entries matched with descriptions\n")
    
    # Save
    outfile <- paste0("picrust2_output_final/", output_prefix, "_abundance.csv")
    write.csv(found_data_desc, outfile, row.names = FALSE)
    cat("Saved:", outfile, "(", nrow(found_data), "features )\n")
  } else {
    cat("No features found for", output_prefix, "\n")
  }
}

process_enzyme_file <- function(file_path, id_col, target_list) {
  if (file.exists(file_path)) {
    cat("\nProcessing", file_path, "...\n")
    data <- read_delim(file_path, delim = "\t", trim_ws = TRUE) %>% rename(!!id_col := 1)
    if (id_col == "KO") data$KO <- gsub("ko:", "", data$KO)
    
    for (group in names(target_list)) {
      save_abundance(data, target_list[[group]][[id_col]], group, id_col, paste0(group, "_", id_col))
    }
  } else {
    cat("File not found:", file_path, "\n")
  }
}

# ---------------------------------------------------------
# 5. Execute
# ---------------------------------------------------------

# 1. Pathways
cat("\n--- Extracting Pathways ---\n")
pathway_df <- abundance %>% rownames_to_column("Pathway")
for (group in names(methane_pathways)) {
  found_pathways <- pathway_df %>% filter(Pathway %in% methane_pathways[[group]])
  if (nrow(found_pathways) > 0) {
    outfile <- paste0("picrust2_output_final/", group, "_pathways_abundance.csv")
    write.csv(found_pathways, outfile, row.names = FALSE)
    cat("Saved:", outfile, "\n")
  }
}

# 2. Enzymes (KO & EC)
cat("\n--- Extracting Enzymes ---\n")
process_enzyme_file("picrust2_output_final/KO_metagenome_out/pred_metagenome_unstrat.tsv", "KO", methane_enzymes)
process_enzyme_file("picrust2_output_final/EC_metagenome_out/pred_metagenome_unstrat.tsv", "EC", methane_enzymes)

cat("\nDone.\n")

