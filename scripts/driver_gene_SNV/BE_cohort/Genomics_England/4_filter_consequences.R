#================================================================#
#================================================================#
######    Script filter out insignificant types of mutations and duplicates
######    Input: combined_consequences.csv, filter_list.txt
######    Output: consequence_gel_samples_filtered_unique.csv
#================================================================#
#================================================================#

library(openxlsx)
library(stringr)
library(purrr)
library(dplyr)

file_path <- "/scratchb/stlab-icgc/users/wu04/project/bo_gene_list/results/gel/combined_consequences.csv"
included_variation_file <- "/scratchb/stlab-icgc/users/wu04/project/bo_gene_list/data/filter_list.txt"
output_dir <- "/scratchb/stlab-icgc/users/wu04/project/bo_gene_list/results/gel/"

data_df <- read.csv(file_path)

# Select included consequences
included_consequences <- readLines(included_variation_file)
included_consequences <- trimws(included_consequences)
included_pattern <- paste(included_consequences, collapse = "|")
print(included_pattern)
print(names(data_df))

filtered_data_df <- data_df %>%
  filter(str_detect(.data[["Consequence"]], included_pattern))
print(names(filtered_data_df))

print('The variation after filtering consequences types')
print(unique(filtered_data_df$Consequence))

consequences_after_filtering <- unique(filtered_data_df$Consequence)

# Remove the undesired variations in a new column
keep_matching_variations <- function(consequence, pattern) {
  # Split consequence into individual variations
  variations <- str_split(consequence, pattern = ",", simplify = FALSE)[[1]]
  print(variations)
  # Filter variations based on pattern
  matched_variations <- variations[str_detect(variations, pattern)]
  
  # Return filtered variations as a single string
  return(paste(matched_variations, collapse = ","))
}

filtered_data_df <- data_df %>%
  rowwise() %>%
  mutate(Consequence_undesired_removed = keep_matching_variations(Consequence, included_pattern)) %>%
  ungroup()
write.csv(filtered_data_df, file = paste0(output_dir, "consequence_gel_samples_filtered.csv"))

# Drop the rows with the combination of Location, symbol, sequence_id, Consequence_undesired_removed are the same
filtered_data_df_unique <- filtered_data_df %>% 
  rename(sequence_id=sample_id) %>%
  filter(Consequence_undesired_removed != "") %>%
  distinct(Location, symbol, sequence_id, Consequence_undesired_removed, .keep_all = TRUE)

filtered_data_df_unique <- filtered_data_df_unique %>%
  mutate(variation = if_else(
    str_detect(Consequence_undesired_removed, "inframe_insertion|inframe_deletion|frameshift"),
    'indel',
    'snp'
  ))

write.csv(filtered_data_df_unique, file = paste0(output_dir, "consequence_gel_samples_filtered_unique.csv"))
