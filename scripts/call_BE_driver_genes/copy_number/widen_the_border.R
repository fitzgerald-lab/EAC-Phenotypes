#=========================================================#
# Script: Widen GISTIC Peaks and Extract Nearby Genes
# Author: Lianlian Wu
# Date: 2024-05-01
# Description:
#   This script expands GISTIC amplification peak regions
#   and retrieves all genes falling into the widened regions
#   using Ensembl GRCh37 (hg19) annotation.
# A peak would often sit next to, but not overlap, a well-characterized oncogene or tumor suppressor. 
# To account for this tendency, we widened the amplification peak sizes upstream and downstream
# by twice the size of each peak to ensure that we captured all possible drivers.
#=========================================================#

# Load required libraries
library(dplyr)
library(stringr)
library(tidyr)
library(GenomicRanges)
library(biomaRt)

#---------------------------------------------------------#
# Step 1: Define file paths
#---------------------------------------------------------#
source_path <- "./Genomic_analysis/results/0_call_BE_driver_genes/copy_number/BE_samples_205/conf_75_q_25_205"

amp_peak_path <- file.path(source_path, 'amp_genes.conf_75.txt')
del_peak_path <- file.path(source_path, 'del_genes.conf_75.txt')  # (defined but not used yet)

#---------------------------------------------------------#
# Step 2: Read and process amplification peak regions
#---------------------------------------------------------#
cat("Reading amplification peaks...\n")

amp_peak_region <- read.table(amp_peak_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
  filter(cytoband == "wide peak boundaries")

# Parse genomic coordinates
amp_peak_region_parsed <- as.data.frame(t(amp_peak_region[, -1])) %>%
  separate(V1, into = c("chromosome", "start", "end"), sep = "[:-]", convert = TRUE) %>%
  filter(!is.na(start)) %>%
  mutate(
    length = end - start,
    start_new = ifelse(start - 2 * length <= 1, 1, start - 2 * length),
    end_new = end + 2 * length,
    length_widen = end_new - start_new
  )

# Save widened peak regions
write.csv(amp_peak_region_parsed, file.path(source_path, 'amp_peak_region_new.csv'), row.names = FALSE)

cat("Amplification peak regions processed and widened.\n")

#---------------------------------------------------------#
# Step 3: Connect to Ensembl GRCh37 and prepare query
#---------------------------------------------------------#
cat("Connecting to Ensembl GRCh37...\n")

ensembl_grch37 <- useEnsembl(biomart = "ensembl",
                             dataset = "hsapiens_gene_ensembl",
                             host = "https://grch37.ensembl.org")

# Check available datasets (optional)
# datasets <- listDatasets(ensembl_grch37)

# Define a query function
query_genes_in_range <- function(chromosome, start, end) {
  chromosome <- gsub("chr", "", chromosome)  # Remove "chr" prefix
  getBM(
    attributes = c('hgnc_symbol', 'chromosome_name', 'start_position', 'end_position'),
    filters = c('chromosome_name', 'start', 'end'),
    values = list(chromosome, start, end),
    mart = ensembl_grch37
  )
}

#---------------------------------------------------------#
# Step 4: Query genes within widened peak regions
#---------------------------------------------------------#
cat("Querying genes within widened peak regions...\n")

genes_list <- lapply(seq_len(nrow(amp_peak_region_parsed)), function(i) {
  cat(sprintf("Processing region %d / %d: %s:%d-%d\n",
              i, nrow(amp_peak_region_parsed),
              amp_peak_region_parsed$chromosome[i],
              amp_peak_region_parsed$start_new[i],
              amp_peak_region_parsed$end_new[i]))
  
  query_genes_in_range(
    chromosome = amp_peak_region_parsed$chromosome[i],
    start = amp_peak_region_parsed$start_new[i],
    end = amp_peak_region_parsed$end_new[i]
  )
})

# Combine results into a single data frame
genes_df <- do.call(rbind, genes_list)

# Extract non-empty gene symbols
genes_extracted <- genes_df$hgnc_symbol
filtered_genes_list <- genes_extracted[genes_extracted != "" & !is.na(genes_extracted)]

# Save the final gene list
writeLines(filtered_genes_list, file.path(source_path, "amp_gene_list_widen.txt"))

cat("Gene extraction completed and saved.\n")