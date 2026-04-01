library(dplyr)
# Set working directory if needed (optional, can be removed if run from project root)
setwd("/Users/jiayi/Desktop/metagenomic_pipeline/QIIME")
# Read feature table (ASV counts per sample)
# biom->TSV files usually have:
#  line 1: "# Constructed from biom file"
#  line 2: "#OTU ID\tSample1\tSample2..."
ft_path <- 'Results/denoise_mode/exported-table/feature-table.tsv'
feature_table_lines <- readLines(ft_path)

# Locate the header line that starts with "#OTU ID"
header_idx <- which(grepl('^#OTU ID', feature_table_lines))[1]
if (is.na(header_idx)) {
    stop('Could not find header line beginning with "#OTU ID" in feature-table.tsv')
}

# Extract column names from that header (remove leading #)
header_raw <- sub('^#', '', feature_table_lines[header_idx])
col_names <- strsplit(header_raw, '\t')[[1]]
# Fallback: if not tab-delimited for some reason, split on whitespace
if (length(col_names) == 1L) {
    col_names <- strsplit(header_raw, '\\s+')[[1]]
}

# Read the data lines after the header line
feature_table <- read.table(ft_path,
                                                        sep='\t', header=FALSE, row.names=1,
                                                        skip=header_idx, check.names=FALSE,
                                                        comment.char='')

# Assign correct column names (skip first column name which is for row names)
if (length(col_names) - 1L != ncol(feature_table)) {
    stop(sprintf('Column name count mismatch: header has %d data columns, table has %d',
                             length(col_names) - 1L, ncol(feature_table)))
}
colnames(feature_table) <- col_names[-1]

# Read taxonomy (Feature ID -> Taxon mapping)
taxonomy <- read.table('Results/denoise_mode/exported-taxonomy/taxonomy.tsv', 
                      sep='\t', header=TRUE, row.names=1)

# Extract only the Taxon column from taxonomy
taxonomy_only <- data.frame(Taxon = taxonomy$Taxon, row.names = rownames(taxonomy))

# Keep only features that have taxonomy assigned
# Find common feature IDs (preserve feature table order)
common_features <- intersect(rownames(feature_table), rownames(taxonomy_only))

# Subset matched rows
feature_table_matched <- feature_table[common_features, , drop=FALSE]
taxonomy_matched <- taxonomy_only[common_features, , drop=FALSE]

# Merge: add Taxon column as the first column, then all sample counts
merged <- cbind(taxonomy_matched, feature_table_matched)

# Save with proper column names
write.table(merged, 'Results/denoise_mode/taxonomy-abundance-table.tsv', 
            sep='\t', quote=FALSE, row.names=TRUE)


parse_taxonomy <- function(tax_string) {
    ranks <- c(Domain=NA, Phylum=NA, Class=NA, Order=NA, Family=NA, Genus=NA, Species=NA)
    if (is.na(tax_string) || tax_string == '') {
        ranks[is.na(ranks)] <- 'Unassigned'
        return(ranks)
    }
    parts <- unlist(strsplit(tax_string, ';'))
    for (p in parts) {
        p <- trimws(p)
        if (p == '' ) next
        if (!grepl('^[dpcofgs]__', p)) next  # skip unexpected tokens
        prefix <- substr(p, 1, 1)
        value <- sub('^[dpcofgs]__', '', p)
        # Normalize empty or placeholder values
        if (value == '' || value == 'uncultured' || value == 'metagenome') value <- 'Unassigned'
        rankName <- switch(prefix,
                                             d='Domain', p='Phylum', c='Class', o='Order',
                                             f='Family', g='Genus', s='Species')
        ranks[rankName] <- ifelse(value == '', 'Unassigned', value)
    }
    ranks[is.na(ranks)] <- 'Unassigned'
    return(ranks)
}

rank_matrix <- t(vapply(taxonomy_matched$Taxon, parse_taxonomy, character(7)))
rank_df <- as.data.frame(rank_matrix, stringsAsFactors = FALSE)
rownames(rank_df) <- rownames(taxonomy_matched)


merged_ranks <- cbind(rank_df, feature_table_matched)

write.table(merged_ranks, 'Results/denoise_mode/final-table-with-ranks.tsv',
                        sep='\t', quote=FALSE, row.names=TRUE)
