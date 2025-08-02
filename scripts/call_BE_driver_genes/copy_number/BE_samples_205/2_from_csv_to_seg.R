#======================================================================#
#======================================================================#
######                                                            ######
######   Script convert .csv to .seg format required by GISTIC2   ######
######                                                            ######
#======================================================================#
#======================================================================#

# Author: Lianlian Wu
# Date: 2024-04-30

library(dplyr)
library(stringr)
library(readr)

data_dir <- "./Genomic_analysis/data/0_call_BE_driver_genes/copy_number/BE_samples_205/"
output_dir <- '/scratchb/stlab-icgc/users/wu04/project/bo_gene_list/data/aks_trio_ross_noora/cna/'

csv_files <- list.files(data_dir, pattern = "copynumber_loged_with_marker\\.csv$", full.names = TRUE, recursive = TRUE)
#print(csv_files)

process_and_save <- function(file_path) {
  # Read CSV file
  data_df <- read.csv(file_path)
  
  # Rename and select specific columns
  data_df <- data_df %>%
    rename(
      Sample = sample_id,
      `Start Position` = Start,
      `End Position` = End,
      `Num Markers` = Unique_Markers,
      #`Seg.CN` = ploidy_adjusted_log2_minus_1_CN
      `Seg.CN` = cna_log2_minus_1
    ) %>%
    select(Sample, Chromosome, `Start Position`, `End Position`, `Num Markers`, `Seg.CN`) %>%
    mutate(
      `Start Position` = as.integer(`Start Position`),
      `End Position` = as.integer(`End Position`),
      `Num Markers` = as.integer(`Num Markers`),
      `Seg.CN` = as.numeric(`Seg.CN`),
      Chromosome = as.character(Chromosome),
      Sample = as.character(Sample))
  
  # Print structure and head of DataFrame (optional, for verification)
  #print(str(data_df))
  #print(head(data_df))
  
  # Construct output filename based on the input file's name
  output_file_name <- sub(".csv$", ".seg", basename(file_path))
  
  # Save the DataFrame to a .seg file in the same directory
  write.table(data_df, file = file.path(dirname(file_path), output_file_name), sep = "\t", row.names = FALSE, quote = FALSE)
}

lapply(csv_files, process_and_save)

print('===============Process for each .csv to .seg completed================')

seg_files <- list.files(data_dir, pattern = "\\.copynumber_loged_with_marker\\.seg$", full.names = TRUE, recursive = TRUE)

#all_seg_data <- lapply(seg_files, read.table, sep = "\t", header = TRUE) %>% 
 # bind_rows()

all_seg_data <- lapply(seg_files, function(file) {
  read.table(file, sep = "\t", header = TRUE, colClasses = c(Chromosome = "character"))
}) %>% 
  bind_rows()

# Exclude six GM sample from the cohort analysis 
excluded_samples <- c("list of excluded samples here")  # Replace with actual sample names
filtered_seg_data <- all_seg_data %>%
  filter(!grepl(paste(excluded_samples, collapse="|"), Sample))

output_file_path <- paste0(data_dir, 'cna_for_gistic2_205.seg')
write.table(filtered_seg_data, file = output_file_path, sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

