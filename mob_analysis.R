library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

setwd("/Users/jiayi/Desktop/metagenomic_pipeline/QIIME")

input_file <- "Results/denoise_mode/final-table-with-ranks.tsv"
output_plot <- "Results/denoise_mode/stacked_barplot_MOB.pdf"

# Load and preprocess data
header <- scan(input_file, what = "character", nlines = 1, sep = "\t", quiet = TRUE)
data <- read.table(input_file, sep = "\t", skip = 1, header = FALSE, stringsAsFactors = FALSE)

if (ncol(data) == length(header) + 1) {
  colnames(data) <- c("FeatureID", header)
} else {
  colnames(data) <- header
}

metadata_cols <- c("FeatureID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
sample_cols <- setdiff(colnames(data), metadata_cols)
sample_cols <- setdiff(sample_cols, c("A1", "A2"))

# Save full sample list (including samples with zero counts) for downstream completeness
ALL_SAMPLE_IDS <- sample_cols
cat("Detected sample columns (", length(ALL_SAMPLE_IDS), "):\n")
print(ALL_SAMPLE_IDS)

long_data <- data %>%
  pivot_longer(cols = all_of(sample_cols), names_to = "SampleID", values_to = "Count") %>%
  filter(Count > 0) %>%
  mutate(
    Domain = trimws(str_replace_all(Domain, "[a-z]__", "")),
    Phylum = trimws(str_replace_all(Phylum, "[a-z]__", "")),
    Class = trimws(str_replace_all(Class, "[a-z]__", "")),
    Order = trimws(str_replace_all(Order, "[a-z]__", "")),
    Family = trimws(str_replace_all(Family, "[a-z]__", "")),
    Genus = trimws(str_replace_all(Genus, "[a-z]__", ""))
  )

# Remove Cyanobacteria and Chloroplast
cat("Filtering out Cyanobacteria and Chloroplast...\n")
cat("  Features before filtering:", nrow(long_data), "\n")

long_data <- long_data %>%
  filter(
    !grepl("Cyanobacteria", Phylum, ignore.case = TRUE),
    !grepl("Cyanobacteriia", Class, ignore.case = TRUE),
    !grepl("Chloroplast", Order, ignore.case = TRUE),
    !grepl("Chloroplast", Family, ignore.case = TRUE),
    !grepl("Chloroplast", Genus, ignore.case = TRUE)
  )

cat("  Features after filtering:", nrow(long_data), "\n\n")

# Filter datasets by domain
# Dataset 1: Bacteria only
long_data_bacteria <- long_data %>%
  filter(grepl("Bacteria", Domain, ignore.case = TRUE))

# Dataset 2: Bacteria and Archaea
long_data_bacteria_archaea <- long_data %>%
  filter(grepl("Bacteria|Archaea", Domain, ignore.case = TRUE))

cat("Dataset summary:\n")
cat("  Total features:", nrow(long_data), "\n")
cat("  Bacteria only:", nrow(long_data_bacteria), "\n")
cat("  Bacteria + Archaea:", nrow(long_data_bacteria_archaea), "\n\n")

# Create output directories for each dataset
dir.create("Results/denoise_mode/bacteria_only", showWarnings = FALSE, recursive = TRUE)
dir.create("Results/denoise_mode/bacteria_archaea", showWarnings = FALSE, recursive = TRUE)

# Function to run MOB analysis on a dataset
run_mob_analysis <- function(data_subset, output_prefix, dataset_name) {
  cat("\n=== Processing", dataset_name, "===\n")
  
  # MOB identification and filtering
  target_genera <- c(
    "Methylacidiphilum", "Crenothrix", "Methylibium", "Methylobacillus", 
    "Methylobacterium", "Methylocella", "Methylomonas", "Methylonatrum", 
    "Methylophaga", "Methylopila", "Methylosinus", "Methylovirgula"
  )
  
  target_pattern <- paste(target_genera, collapse = "|")
  
  filtered_data <- data_subset %>%
    mutate(
      MatchedGenus = str_extract(Genus, regex(target_pattern, ignore_case = TRUE)),
      IsUnclassifiedMethylophilaceae = grepl("Methylophilaceae", Family, ignore.case = TRUE) & 
                                       (Genus == "" | Genus == "uncultured" | Genus == "Unassigned"),
      IsUnclassifiedMethylacidiphilaceae = grepl("Methylacidiphilaceae", Family, ignore.case = TRUE) & 
                                       (Genus == "" | Genus == "uncultured" | Genus == "Unassigned")
    ) %>%
    filter(!is.na(MatchedGenus) | IsUnclassifiedMethylophilaceae | IsUnclassifiedMethylacidiphilaceae) %>%
    mutate(
      Taxon = case_when(
        !is.na(MatchedGenus) ~ str_to_title(MatchedGenus),
        IsUnclassifiedMethylophilaceae ~ "Unclassified Methylophilaceae",
        IsUnclassifiedMethylacidiphilaceae ~ "Unclassified Methylacidiphilaceae"
      ),
      Taxon = ifelse(Taxon == "Methylacidiphilum", "Candidatus Methylacidiphilum", Taxon)
    )
  
  if (nrow(filtered_data) == 0) {
    cat("WARNING: No MOB species found in", dataset_name, "\n")
    return(NULL)
  }
  
  # MOB relative abundance calculation and visualization
  plot_data <- filtered_data %>%
    group_by(SampleID, Taxon) %>%
    summarise(Count = sum(Count), .groups = "drop") %>%
    group_by(SampleID) %>%
    mutate(RelativeAbundance = Count / sum(Count) * 100) %>%
    ungroup()
  
  plot_data <- plot_data %>%
    mutate(Group = case_when(
      grepl("B$", SampleID) ~ "Bark",
      grepl("S$", SampleID) ~ "Soil",
      TRUE ~ "Wood"
    ))
  
  mob_colors <- c(
    "Candidatus Methylacidiphilum" = "#CD5C5C",
    "Crenothrix" = "black",
    "Methylibium" = "#FF8C00",
    "Methylobacillus" = "#98FB98",
    "Methylobacterium" = "#7FFF00",
    "Methylocella" = "#8B4513",
    "Methylomonas" = "#E6E6FA",
    "Methylonatrum" = "#483D8B",
    "Methylophaga" = "#800080",
    "Methylopila" = "#C0C0C0",
    "Methylosinus" = "#DA70D6",
    "Unclassified Methylophilaceae" = "#AFEEEE"
  )
  
  plot_data$Taxon <- factor(plot_data$Taxon, levels = sort(unique(plot_data$Taxon)))
  
  p <- ggplot(plot_data, aes(x = SampleID, y = RelativeAbundance, fill = Taxon)) +
    geom_bar(stat = "identity", width = 0.8) +
    facet_grid(~ Group, scales = "free_x", space = "free_x") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = mob_colors) +
    labs(x = "Tree", y = "Relative abundance (%)", fill = "MOB - 16S", 
         title = paste("MOB Distribution -", dataset_name)) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5), 
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "white", colour = "black"),
      strip.text = element_text(size = 12),
      legend.position = "right",
      legend.text = element_text(size = 10, face = "italic")
    )
  
  print(p)
  ggsave(paste0(output_prefix, "_stacked_barplot_MOB.pdf"), p, width = 10, height = 6)
  
  # MOB total relative abundance per sample
  # Ensure all samples (including those with zero total counts) are represented
  samples_df <- tibble(SampleID = ALL_SAMPLE_IDS)

  mob_totals <- filtered_data %>%
    group_by(SampleID) %>%
    summarise(MOB_Total = sum(Count), .groups = "drop") %>%
    right_join(samples_df, by = "SampleID") %>%
    mutate(MOB_Total = replace_na(MOB_Total, 0))

  # Use original wide 'data' to compute total abundance per sample (includes zeros)
  # safe colSums conversion
  total_abundance_vec <- rep(0, length(ALL_SAMPLE_IDS))
  names(total_abundance_vec) <- ALL_SAMPLE_IDS
  if (exists("data") && length(ALL_SAMPLE_IDS) > 0) {
    # coerce to numeric if necessary
    mat <- data[ALL_SAMPLE_IDS]
    mat2 <- as.data.frame(lapply(mat, function(x) as.numeric(as.character(x))))
    total_abundance_vec <- colSums(mat2, na.rm = TRUE)
  }

  all_totals <- tibble(SampleID = ALL_SAMPLE_IDS, Total_Abundance = as.numeric(total_abundance_vec[ALL_SAMPLE_IDS]))

  mob_rel_abundance <- mob_totals %>%
    left_join(all_totals, by = "SampleID") %>%
    mutate(
      MOB_RelAbundance = ifelse(Total_Abundance > 0, MOB_Total / Total_Abundance * 100, 0),
      Species = str_extract(SampleID, "^[A-Z]"),
      Surrounding = case_when(
        grepl("B$", SampleID) ~ "Bark",
        grepl("S$", SampleID) ~ "Soil",
        TRUE ~ "Wood"
      )
    )
  
  mob_per_sample <- mob_rel_abundance %>%
    select(SampleID, Species, Surrounding, MOB_RelAbundance) %>%
    rename(`MOB_RelAbundance (%)` = MOB_RelAbundance)
  
  print("=== MOB Relative Abundance Per Sample ===")
  print(mob_per_sample)
  # write.csv(mob_per_sample, paste0(output_prefix, "_MOB_relative_abundance_per_sample.csv"), row.names = FALSE)

  # --- Additional: Count unique Genus+Species per sample and compare MOB% across surroundings ---
  # Build a combined Genus+Species identifier to avoid double-counting same species
  taxon_counts <- data_subset %>%
    mutate(
      Genus_clean = trimws(Genus),
      Species_clean = trimws(Species),
      GenusSpecies = case_when(
        !is.na(Genus_clean) & Genus_clean != "" & !tolower(Genus_clean) %in% c("uncultured", "unassigned") &
          !is.na(Species_clean) & Species_clean != "" & !tolower(Species_clean) %in% c("uncultured", "unassigned") ~ paste(Genus_clean, Species_clean),
        !is.na(Genus_clean) & Genus_clean != "" & !tolower(Genus_clean) %in% c("uncultured", "unassigned") ~ Genus_clean,
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(GenusSpecies)) %>%
    group_by(SampleID) %>%
    summarise(Genus_Species_Count = n_distinct(GenusSpecies), .groups = "drop")

  # Join combined taxon counts to MOB per-sample table
  mob_per_sample_counts <- mob_per_sample %>%
    left_join(taxon_counts, by = "SampleID") %>%
    mutate(Genus_Species_Count = replace_na(Genus_Species_Count, 0))

  write.csv(mob_per_sample_counts, paste0(output_prefix, "_MOB_per_sample_with_genus_species_counts.csv"), row.names = FALSE)

  # --- Statistical test: Kruskal-Wallis for Surrounding ---
  surrounding_levels <- unique(mob_per_sample_counts$Surrounding)
  if (length(surrounding_levels) >= 2) {
    kruskal_surrounding <- tryCatch({
      kruskal.test(`MOB_RelAbundance (%)` ~ Surrounding, data = mob_per_sample_counts)
    }, error = function(e) e)
    
    print("=== Kruskal-Wallis Test: MOB abundance by Surrounding ===")
    print(kruskal_surrounding)
    
    cat("Kruskal-Wallis test for Surrounding completed.\n")
    
    # Pairwise Wilcoxon tests for Surrounding
    if (length(surrounding_levels) >= 2) {
      pairwise_surrounding <- tryCatch({
        pairwise.wilcox.test(mob_per_sample_counts$`MOB_RelAbundance (%)`, 
                            mob_per_sample_counts$Surrounding, 
                            p.adjust.method = "BH", 
                            exact = FALSE)
      }, error = function(e) e)
      
      if (!inherits(pairwise_surrounding, "error")) {
        print("\n=== Pairwise Wilcoxon Test: MOB abundance by Surrounding (BH adjusted) ===")
        print(pairwise_surrounding)
        cat("\n")
      }
    }
  } else {
    cat("Not enough surrounding levels for Kruskal-Wallis test.\n")
  }

  # --- Statistical test: Kruskal-Wallis for Species ---
  species_levels <- unique(mob_per_sample_counts$Species)
  if (length(species_levels) >= 2) {
    kruskal_species <- tryCatch({
      kruskal.test(`MOB_RelAbundance (%)` ~ Species, data = mob_per_sample_counts)
    }, error = function(e) e)
    
    print("=== Kruskal-Wallis Test: MOB abundance by Species ===")
    print(kruskal_species)
    
    cat("Kruskal-Wallis test for Species completed.\n")
    
    # Pairwise Wilcoxon tests for Species
    if (length(species_levels) >= 2) {
      pairwise_species <- tryCatch({
        pairwise.wilcox.test(mob_per_sample_counts$`MOB_RelAbundance (%)`, 
                            mob_per_sample_counts$Species, 
                            p.adjust.method = "BH", 
                            exact = FALSE)
      }, error = function(e) e)
      
      if (!inherits(pairwise_species, "error")) {
        print("\n=== Pairwise Wilcoxon Test: MOB abundance by Species (BH adjusted) ===")
        print(pairwise_species)
        cat("\n")
      }
    }
  } else {
    cat("Not enough species levels for Kruskal-Wallis test.\n")
  }

  # Pairwise comparisons for Species
  # species_levels <- unique(mob_per_sample_counts$Species)
  # pairwise_species <- NULL
  # if (length(species_levels) >= 2) {
  #   pairwise_species <- tryCatch({
  #     pw2 <- pairwise.wilcox.test(mob_per_sample_counts$`MOB_RelAbundance (%)`, mob_per_sample_counts$Species, p.adjust.method = "BH", exact = FALSE)
  #     pw2
  #   }, error = function(e) e)
  # 
  #   if (!inherits(pairwise_species, "error")) {
  #     write.csv(as.data.frame.matrix(pairwise_species$p.value), paste0(output_prefix, "_pairwise_species_pvalues.csv"), row.names = TRUE)
  #   } else {
  #     writeLines(as.character(pairwise_species), con = paste0(output_prefix, "_pairwise_species_pvalues.txt"))
  #   }
  # } else {
  #   writeLines("Not enough species levels for pairwise comparisons", con = paste0(output_prefix, "_pairwise_species_pvalues.txt"))
  # }
  # 
  # # If exactly two species, also perform standard Wilcoxon test and save statistic
  # if (length(species_levels) == 2) {
  #   try({
  #     wres <- wilcox.test(`MOB_RelAbundance (%)` ~ Species, data = mob_per_sample_counts, exact = FALSE)
  #     sink(paste0(output_prefix, "_species_wilcox_test.txt"))
  #     print(wres)
  #     sink()
  #   }, silent = TRUE)
  # }

  # Boxplot: MOB percentage by surrounding
  p_mob_box <- ggplot(mob_per_sample_counts, aes(x = Surrounding, y = `MOB_RelAbundance (%)`)) +
    geom_boxplot(fill = "#87CEEB") +
    geom_jitter(width = 0.2, height = 0, size = 1, alpha = 0.8) +
    labs(x = "Surrounding", y = "MOB Relative Abundance (%)", title = paste("MOB % by Surrounding -", dataset_name)) +
    theme_bw()

  ggsave(paste0(output_prefix, "_MOB_by_surrounding_boxplot.pdf"), p_mob_box, width = 6, height = 4)

  # Plot Genus+Species count distribution by surrounding
  genus_by_sample <- mob_per_sample_counts %>%
    select(SampleID, Surrounding, Genus_Species_Count)

  p_genus_box <- ggplot(genus_by_sample, aes(x = Surrounding, y = Genus_Species_Count)) +
    geom_boxplot(fill = "#FFD700") +
    geom_jitter(width = 0.2, height = 0, size = 1, alpha = 0.8) +
    labs(x = "Surrounding", y = "Unique Genus+Species Count", title = paste("Genus+Species counts by Surrounding -", dataset_name)) +
    theme_bw()

  ggsave(paste0(output_prefix, "_genus_species_counts_by_surrounding_boxplot.pdf"), p_genus_box, width = 6, height = 4)

  cat("Wrote MOB per-sample table with genus+species counts and summary/plots for Surrounding comparison.\n")
  
  # Dominant vs rare MOB group classification
  mob_genus_totals <- plot_data %>%
    group_by(Taxon) %>%
    summarise(Total_RelAbundance = sum(RelativeAbundance)) %>%
    arrange(desc(Total_RelAbundance)) %>%
    ungroup()
  
  key_groups <- mob_genus_totals %>%
    slice(1:4) %>%
    pull(Taxon)
  
  mob_classification <- mob_genus_totals %>%
    mutate(
      Group_Type = case_when(
        Taxon %in% key_groups ~ "Dominant/Key Group",
        TRUE ~ "Rare Group"
      ),
      Rank = row_number()
    )
  
  print("=== MOB Group Classification ===")
  print(mob_classification)
  write.csv(mob_classification, paste0(output_prefix, "_MOB_group_classification.csv"), row.names = FALSE)
  
  # Key species identification by site
  mob_by_site <- plot_data %>%
    group_by(Group, Taxon) %>%
    summarise(Total_RelAbundance = sum(RelativeAbundance), .groups = "drop") %>%
    arrange(Group, desc(Total_RelAbundance))
  
  key_species_by_site <- mob_by_site %>%
    group_by(Group) %>%
    slice(1:3) %>%
    mutate(Rank = row_number()) %>%
    ungroup()
  
  print("\n=== Key MOB Species by Site (Top 3 per Surrounding) ===")
  print(key_species_by_site)
  write.csv(key_species_by_site, paste0(output_prefix, "_MOB_key_species_by_site.csv"), row.names = FALSE)
  
  cat("\n", dataset_name, "analysis complete\n")
}

# Run analysis for both datasets
run_mob_analysis(long_data_bacteria, "Results/denoise_mode/bacteria_only/bacteria_only", "Bacteria Only")
run_mob_analysis(long_data_bacteria_archaea, "Results/denoise_mode/bacteria_archaea/bacteria_archaea", "Bacteria + Archaea")

