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
sample_211_data_path <- "/mnt/scratchc/stlab-icgc/users/wu04/project/bo_gene_list/data/aks_trio_ross_noora/cna/sample_211"
result_dir <- "/mnt/scratchc/stlab-icgc/users/wu04/project/bo_gene_list/results/samples_211/cna"

# Load data 
combined_cna_data <- read.csv(paste0(base_dir, "/combined_cna.csv"))

data_df <- combined_cna_data %>%
    mutate(Length = End - Start,
    Major_CN = Total_CN - Minor_CN) %>% 
    select(sample_id, Chromosome, Start, End, Length, Total_CN, Major_CN, Minor_CN, ploidy) %>%
    filter(Chromosome != "X", Chromosome != "Y")
head(data_df)

# Prepare data separating p and q chromosome arms =========================================================
wget("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/gap.txt.gz",
     destfile = paste0(base_dir, "/arm_level_cna.csv")

# Calculate arm-level CNA =================================================================================
# Step 1: Define centromere positions for hg19 (example subset, expand as needed)
centromere_hg19 <- data.frame(
  Chromosome = as.character(1:22),
  Centromere = c(
    121535434, 92326171, 90504854, 49660117, 46405641,
    58830166, 58054331, 43838887, 47367679, 39254935,
    51644205, 34856694, 16000000, 16000000, 17000000,
    36000000, 22000000, 17000000, 26500000, 27500000,
    13000000, 13700000
  )
)

# Step 2: Assign p or q arm
data_df_with_arm <- data_df %>%
  left_join(centromere_hg19, by = "Chromosome") %>%
  mutate(
    arm = case_when(
      End <= Centromere ~ "p",
      Start >= Centromere ~ "q",
      Start < Centromere & End > Centromere ~ "pq",  # split interval
      TRUE ~ NA_character_
    )
  ) %>%
  filter(arm %in% c("p", "q"))  # Exclude overlapping regions

# Step 3: Weighted average CN values per sample, chromosome, and arm
weighted_CN <- data_df_with_arm %>%
  group_by(sample_id, Chromosome, arm) %>%
  summarise(
    total_length = sum(Length, na.rm = TRUE),
    weighted_Total_CN = sum(Total_CN * Length, na.rm = TRUE) / total_length,
    weighted_Major_CN = sum(Major_CN * Length, na.rm = TRUE) / total_length,
    weighted_Minor_CN = sum(Minor_CN * Length, na.rm = TRUE) / total_length,
    ploidy = first(ploidy),
    .groups = "drop"
  ) %>%
  mutate(
    event = case_when(
      round(weighted_Total_CN) > round(ploidy) ~ "gain",
      round(weighted_Total_CN) < round(ploidy) ~ "loss",
      TRUE ~ "neutral"
    ),
    LOH = ifelse(weighted_Minor_CN < 0.1, TRUE, FALSE)
  )
