#=====================================================================#
#=====================================================================#
######                                                           
######      Script process CNA data to calculate fraction of aberrant genome, LOH, and WGD
######                                                           
#=====================================================================#
#=====================================================================#

# Author: Lianlian Wu
# Date: 2025-05-22

library(dplyr)
library(ggplot2)

# Define the source directory 
base_dir <- "/mnt/scratchc/stlab-icgc/users/wu04/project/bo_gene_list/data/aks_trio_ross_noora/cna/sample_211"
result_dir <- "/mnt/scratchc/stlab-icgc/users/wu04/project/bo_gene_list/results/samples_211/cna"

# List all CSV files recursively
# cna_files <- list.files(path = base_dir, pattern = "copynumber.caveman_adjusted.csv$", full.names = TRUE, recursive = TRUE)

# Combine all CNA data and save the combined one
# combined_cna_data <- cna_files %>%
#     lapply(function(file) {
#         data <- read.csv(file)
#         data$Chromosome <- as.character(data$Chromosome)
#         return(data)
#     }) %>%
#     bind_rows()

#write.csv(combined_cna_data, paste0(base_dir, "/combined_cna.csv"), row.names = FALSE)

# Load data 
combined_cna_data <- read.csv(paste0(base_dir, "/combined_cna.csv"))

data_df <- combined_cna_data %>%
    mutate(Length = End - Start,
    Major_CN = Total_CN - Minor_CN) %>% 
    select(sample_id, Chromosome, Start, End, Length, Total_CN, Major_CN, Minor_CN, ploidy) %>%
    filter(Chromosome != "X", Chromosome != "Y")
head(data_df)

# Calculate the fraction of aberrant genome==========================================================
fraction_aberrant <- data_df  %>%  # Exclude X and Y chromosomes
  group_by(sample_id) %>%
  summarize(Total_Length = sum(Length), .groups = 'drop') %>%  # Summing total length of all regions
  full_join(
    data_df %>%
      filter(round(Total_CN) != round(ploidy)) %>%  # Aberrant regions based on CN != ploidy or LOH
      group_by(sample_id) %>%
      summarize(aberrant_Length = sum(Length), .groups = 'drop'),  # Summing aberrant region lengths
    by = "sample_id"
  ) %>%
  mutate(fraction_aberrant_genome = ifelse(is.na(aberrant_Length), 0, aberrant_Length / Total_Length))  # Calculate fraction of aberrant length

write.csv(fraction_aberrant, paste0(result_dir, "/fraction_aberrant.csv"), row.names = FALSE)

# Calculate the fraction of LOH genome==========================================================
fraction_LOH <- data_df %>%
  group_by(sample_id) %>%
  summarize(Total_Length = sum(Length), .groups = 'drop') %>%  # Summing total length of all regions
  inner_join(
    data_df %>%
      filter(Major_CN>0 & Minor_CN==0) %>% 
      group_by(sample_id) %>%
      summarize(aberrant_LOH_Length = sum(Length), .groups = 'drop'),  # Summing aberrant region lengths
    by = "sample_id"
  ) %>%
  mutate(fraction_LOH = aberrant_LOH_Length / Total_Length)

write.csv(fraction_LOH, paste0(result_dir, "/fraction_LOH.csv"), row.names = FALSE)

# WGD analysis (Definition: More than half of the genome with Major_CN >= 2)========================================
WGD_analysis <- data_df %>%
  group_by(sample_id) %>%
  summarize(Total_Length = sum(Length), .groups = 'drop') %>%  # Summing total length of all regions
  inner_join(
    data_df %>%
      filter(Major_CN >= 2) %>%  # Aberrant regions with WGD
      group_by(sample_id) %>%
      summarize(WGD_Length = sum(Length), .groups = 'drop'),  # Summing aberrant region lengths
    by = "sample_id"
  ) %>%
  mutate(fraction_WGD = WGD_Length / Total_Length,
      WGD = ifelse(fraction_WGD>=0.5, TRUE, FALSE)) 

write.csv(WGD_analysis, paste0(result_dir, "/WGD_analysis.csv"), row.names = FALSE)
