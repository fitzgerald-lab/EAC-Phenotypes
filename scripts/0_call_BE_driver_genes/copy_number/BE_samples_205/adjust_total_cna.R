#================================================================#
# Script: Adjust CNA Values to Ploidy-adjusted CNA
# Author: Lianlian Wu
# Date: 2024-04-30
# Description:
#   This script adjusts raw copy number values to ploidy-adjusted
#   copy number values for downstream analysis.
#================================================================#

# Load required libraries
library(dplyr)

# Define the source directory 
base_dir <- "/path/to/output/folder/of/CAVAMEN/"

# List all CSV files recursively
ploidy_files <- list.files(path = base_dir, pattern = "samplestatistics.csv$", full.names = TRUE, recursive = TRUE)
cna_files <- list.files(path = base_dir, pattern = "copynumber.caveman.csv$", full.names = TRUE, recursive = TRUE)
cat("Number of CNA files found:", length(cna_files), "\n")

# Generate new columns with ploidy_adjusted_Total_CN and ploidy_adjusted_log2_minus_1_CN, and save files for further analysis
# Loop over all files and process
for (i in seq_along(ploidy_files)) {
    ploidy_data_path <- ploidy_files[i]
    cna_data_path <- cna_files[i]
    
    print(paste("==============Processing iteration:", i))
    print(cna_data_path)
    print(ploidy_data_path)

    # Read and process ploidy data
    ploidy_data <- read.csv(ploidy_data_path, header = FALSE, sep = " ")
    names(ploidy_data) <- c("Parameter", "Value")
    ploidy_value <- as.numeric(ploidy_data$Value[ploidy_data$Parameter == "Ploidy"])
    print(ploidy_value)

    # Read and process CNA data
    cna_data <- read.csv(cna_data_path)
    cna_data$ploidy <- ploidy_value
    
    cna_data$ploidy_adjusted_Total_CN <- cna_data$Total_CN / ploidy_value
    cna_data$ploidy_adjusted_log2_minus_1_CN <- log2(cna_data$ploidy_adjusted_Total_CN) - 1
    #print(cna_data$ploidy_adjusted_log2_minus_1_CN)
    # Extract sample ID from the filename and adjust for saving
    pattern_ <- "__pv.*__rg\\.grch37.*"
    file_name <- basename(cna_data_path)
    sample_id <- sub(pattern_, "", file_name)
    cna_data$sample_id <- sample_id
    print(file_name)
    # Define new filename for saving the adjusted CNA data
    adjusted_filename <- sub("copynumber.caveman.csv$", "copynumber.caveman_adjusted.csv", cna_data_path)
    print(adjusted_filename)
    # Write the adjusted CNA data back to the same directory
    write.csv(cna_data, adjusted_filename, row.names = FALSE)
}
cat("\n=========== Adjustment Completed ===========\n")