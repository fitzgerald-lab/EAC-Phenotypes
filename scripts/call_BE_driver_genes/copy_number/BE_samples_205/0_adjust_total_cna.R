#================================================================#
# Script: Adjust CNA Values called by ASCAT to Seg.CN (log2() -1 of copy number) required by GISTIC2
# Author: Lianlian Wu
# Date: 2024-04-30
# Description:
#   This script adjusts raw copy number values to log2
#   copy number values for downstream analysis.
#================================================================#

library(dplyr)

# Define the source directory 
base_dir <- "./Genomic_analysis/data/0_call_BE_driver_genes/copy_number/BE_samples_205/"

# List all CSV files recursively
ploidy_files <- list.files(path = base_dir, pattern = "samplestatistics.txt$", full.names = TRUE, recursive = TRUE)
cna_files <- list.files(path = base_dir, pattern = "copynumber.csv$", full.names = TRUE, recursive = TRUE)

print(ploidy_files[1])
print(cna_files[1])

for (i in seq_along(ploidy_files)) {
    ploidy_data_path <- ploidy_files[i]
    cna_data_path <- cna_files[i]
    print(paste0("=============", i, "=============="))
    # Read and process ploidy data
    ploidy_data <- read.table(ploidy_data_path, header = FALSE, sep = " ")
    names(ploidy_data) <- c("Parameter", "Value")
    ploidy_value <- as.numeric(ploidy_data$Value[ploidy_data$Parameter == "Ploidy"])
    print(ploidy_value)
  
    # Read and process CNA data
    cna_data <- read.csv(cna_data_path)
    #cna_data$ploidy_adjusted_Total_CN <- cna_data$Total_CN / ploidy_value
    #cna_data$ploidy_adjusted_log2_minus_1_CN <- log2(cna_data$ploidy_adjusted_Total_CN) - 1
    cna_data$cna_log2_minus_1 <- log2(cna_data[[7]]) - 1
    print(cna_data[[7]][1])
    # Extract sample ID from the filename and adjust for saving
    sample_id <- gsub(".copynumber.caveman.csv$", "", basename(cna_data_path))
    cna_data$sample_id <- sample_id
    
    # Set custom column names
    column_names <- c("index", "Chromosome", "Start", "End", "Total_CN_normal", "Minor_CN_normal", 
                      "Total_CN", "Minor_CN", "cna_log2_minus_1", "sample_id")
    names(cna_data) <- column_names

    # Define new filename for saving the adjusted CNA data
    adjusted_filename <- sub("copynumber.csv$", "copynumber.caveman_loged.csv", cna_data_path)

    # Write the adjusted CNA data back to the same directory
    write.csv(cna_data, adjusted_filename, row.names = FALSE)
}
