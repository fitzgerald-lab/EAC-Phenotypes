#================================================================#
#================================================================#
######    Script filter out insignificant types of SNV and duplicates
######    Input: combined_snv_consequences.csv, filter_list.txt
######    Output: consequence_snv_211_samples_filtered_unique_2.csv
#================================================================#
#================================================================#

# Author: Lianlian Wu
# Date: 2024-05-01

library(stringr)
library(openxlsx)
library(purrr)
library(dplyr)
library(tidyr)

file_path <- "/scratchb/stlab-icgc/users/wu04/project/bo_gene_list/results/samples_211/combined_snv_consequences.csv"
included_variation_file <- "/scratchb/stlab-icgc/users/wu04/project/bo_gene_list/data/filter_list.txt"
output_dir <- "/scratchb/stlab-icgc/users/wu04/project/bo_gene_list/results/samples_211/"
sequence_id_file_path <- "/scratchb/stlab-icgc/users/wu04/project/bo_gene_list/data/aks_trio_ross_noora/sequence_id_211.txt"

sequence_ids <- readLines(sequence_id_file_path)

data_df <- read.csv(file_path)

# Filter out the interested samples
data_df$sequence_id <- str_extract(data_df$sample_id, "^(.+?)_vs") 
data_df$sequence_id <- gsub("_vs", "", data_df$sequence_id)

filtered_data_df <- data_df[data_df$sequence_id %in% sequence_ids, ]

print('Before filtering interested samples: The unique number of sequence_id are:')
print(length(unique(data_df$sequence_id)))
print('Before filtering interested samples: The unique number of sequence_id_vs_sequence_id are:')
print(length(unique(data_df$sample_id)))

print('After filtering interested samples: The unique number of sequence_id are:')
print(length(unique(filtered_data_df$sequence_id)))
print('After filtering interested samples: The unique number of sequence_id_vs_sequence_id are:')
print(length(unique(filtered_data_df$sample_id)))

# Select included consequences
included_consequences <- readLines(included_variation_file)
included_consequences <- trimws(included_consequences)
included_pattern <- paste(included_consequences, collapse = "|")
print(included_pattern)

filtered_data_df <- filtered_data_df %>%
  filter(str_detect(Consequence, included_pattern))

print('The variation after filtering consequences types')
print(unique(filtered_data_df$Consequence))

consequences_after_filtering <- unique(filtered_data_df$Consequence)

# Remove the undesired variations in a new column
keep_matching_variations <- function(consequence, pattern) {
  # Split consequence into individual variations
  variations <- str_split(consequence, pattern = ",", simplify = FALSE)[[1]]
  #print(variations)
  # Filter variations based on pattern
  matched_variations <- variations[str_detect(variations, pattern)]
  
  # Return filtered variations as a single string
  return(paste(matched_variations, collapse = ","))
}

filtered_data_df <- filtered_data_df %>%
  rowwise() %>%
  mutate(Consequence_undesired_removed = keep_matching_variations(Consequence, included_pattern)) %>%
  ungroup()
#write.xlsx(filtered_data_df, file = paste0(output_dir, "consequence_snv_146_samples_filtered.xlsx"))

# Drop the rows with the combination of Location, symbol, sequence_id, Consequence_undesired_removed are the same
filtered_data_df_unique <- filtered_data_df %>% distinct(Location, symbol, sequence_id, Consequence_undesired_removed, .keep_all = TRUE)
filtered_data_df_unique$variation <- 'snv'
write.csv(filtered_data_df_unique, file = paste0(output_dir, "consequence_snv_211_samples_filtered_unique_2.csv"))
