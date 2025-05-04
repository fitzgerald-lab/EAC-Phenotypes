#===========================================================================================
#===========================================================================================      
######   Script extract and combine SNV consequences of the cohort   
######   Input: .vep files          
######   Output: combined_snv_consequences.csv, all_snv_count.csv                                                                      
#===========================================================================================
#===========================================================================================

# Author: Lianlian Wu
# Date: 2024-05-01

library(dplyr)
library(tidyr)
library(stringr)
library(openxlsx)
library(ggplot2)

vep_files_dir <- "/scratchb/stlab-icgc/users/wu04/project/bo_gene_list/data/aks_trio_ross_noora/snp/"
output_dir <- "/scratchb/stlab-icgc/users/wu04/project/bo_gene_list/results/samples_211/"

vep_files <- list.files(vep_files_dir, pattern = "\\.vep$", full.names = TRUE)


count_variance_each_gene <- function(vep_file_path) {
  # Load VEP and skip the head
  lines <- readLines(vep_file_path)
  data_start_line <- which(grepl("^#Uploaded_variation", lines))-1
  vep_data <- read.delim(vep_file_path, skip = data_start_line, header = TRUE, stringsAsFactors = FALSE)
  
  # Extract symbol from the Extra column and generate a new df with four columns: sample_id, gene, Consequence, symbol
  pattern_ <- "^(.+?)__rg.grch37"
  file_name <- basename(vep_file_path)
  #print(regexpr(pattern_, file_name))
  sample_id <- regmatches(file_name, regexpr(pattern_, file_name))
  
  vep_data$symbol <- str_extract(vep_data$Extra, "(?<=SYMBOL=)[^;]+")
  
  cleaned_data <- vep_data %>% select(Location, Gene, Consequence, symbol) %>% 
      mutate(sample_id = sample_id) %>% 
      distinct(sample_id, Location, .keep_all = TRUE)
  # Calculate snv count
  snv_count <- nrow(cleaned_data)
  print(sample_id)
  print(snv_count)
  snv_count_df <- data.frame(sample_id = sample_id, snv_count = snv_count)
  
  return(list(cleaned_data = cleaned_data, snv_count_df = snv_count_df))
  #df_count <- cleaned_data %>%  group_by(sample_id, symbol) %>%  
   # mutate(count = n()) %>% ungroup() %>% select(symbol, count, sample_id) %>% 
    #distinct(sample_id, symbol, .keep_all = TRUE)
  #return(list(cleaned_data = cleaned_data, df_count = df_count))
}

# Apply the function to each VEP file and combine the results
all_consequences <- vep_files %>% 
  lapply(function(x) count_variance_each_gene(x)$cleaned_data) %>% 
  bind_rows()

all_snv_count <- vep_files %>% 
  lapply(function(x) count_variance_each_gene(x)$snv_count_df) %>% 
  bind_rows()

# Save the combined result to a xlsx file (CSV file also works, but Microsoft Excel tries to automatically detect and format data looks like date to date when opening CSV files.)
write.csv(all_consequences, file = paste0(output_dir, "combined_snv_consequences.csv"))
write.csv(all_snv_count, file = paste0(output_dir, "all_snv_count.csv"))