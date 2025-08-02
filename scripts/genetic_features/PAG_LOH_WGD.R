#=====================================================================#
#=====================================================================#
######                                                           
######      Script to process CNA data and calculate genomic features:
######      - Fraction of aberrant genome (PAG)
######      - Fraction of genome with loss of heterozygosity (LOH)
######      - Whole genome duplication (WGD) status
######                                                           
#=====================================================================#
#=====================================================================#

# Author: Lianlian Wu
# Date: 2024-09-22

library(dplyr)
library(ggplot2)

# Define directories for input data and results
base_dir <- "path/to/input/data"  # Replace with the actual path to input data
result_dir <- "path/to/output/results"  # Replace with the actual path to save results

# Load data 
combined_cna_data <- read.csv(file.path(base_dir, "combined_cna.csv"))
data_df <- combined_cna_data %>%
    mutate(Start = Start.Position,
    End = End.Position,
    Total_CN = CN,
    Length = End - Start,
    Minor_CN = Total_CN - Major_CN) %>% 
    select(sample_id, Chromosome, Start, End, Length, Total_CN, Major_CN, Minor_CN, ploidy) %>%
    filter(Chromosome != "X", Chromosome != "Y") # Exclude X and Y chromosomes
head(data_df)


# Calculate the fraction of aberrant genome==========================================================
fraction_aberrant <- data_df  %>%  
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

#=======================================#
######                                                           
######      Summarize the data          #
######                                                           
#=======================================#
## Calculate the median of each feature
median_PAG <- median(fraction_aberrant$fraction_aberrant_genome, na.rm = TRUE)
median_LOH <- median(fraction_LOH$fraction_LOH, na.rm = TRUE)
WGD_proportion = mean(WGD_analysis$WGD == TRUE)
print(paste("Median PAG:", median_PAG))
print(paste("Median LOH:", median_LOH))
print(paste("WGD proportion:", WGD_proportion))

# mean_PAG <- mean(fraction_aberrant$fraction_aberrant_genome, na.rm = TRUE)
# mean_LOH <- mean(fraction_LOH$fraction_LOH, na.rm = TRUE)
# mean_WGD <- mean(WGD_analysis$fraction_WGD, na.rm = TRUE)
