#===========================================================================================
#===========================================================================================      
######   Script extract and combine consequences of the cohort   
######   Input: ../snv_hg19_filter_pass/*.vcf files          
######   Output: combined_consequences.csv                                                                     
#===========================================================================================
#===========================================================================================

# Author: Lianlian Wu
# Date: 2024-05-01

library(dplyr)
library(tidyr)
library(stringr)
library(openxlsx)

vcf_files_dir <- "path/to/vcf_files_directory"
output_dir <- "path/to/output_directory"

vcf_files <- list.files(vcf_files_dir, pattern = "\\.vcf$", full.names = TRUE)

# Create a named vector for renaming columns
column_names <- c(V1 = "chromosome",
                  V2 = "location",
                  V3 = "snp_id",
                  V4 = "ori_allele",
                  V5 = "mut_allele",
                  V6 = "V6",  # If V6 is to keep its name
                  V7 = "filter",
                  V8 = "mut_details",
                  V9 = "indexs",
                  V10 = "values_1",
                  V11 = "values_2")

process_vcf_data <- function(vcf_data, column_names) {
  # Apply column names
  names(vcf_data) <- column_names[names(vcf_data)]
  
  # Perform mutations
  vcf_data <- vcf_data %>%
    mutate(
      qc = if_else(str_detect(mut_details, "CSQT=1\\|"), 
                   sub("CSQT=1\\|.*$", "", mut_details), 
                   mut_details),
      mut_details_all = if_else(str_detect(mut_details, "CSQT=1\\|"), 
                                sub(".*?CSQT=1\\|", "", mut_details), 
                                NA)
    ) %>%
    mutate(mut_details_all = sub(";CSQR=1.*", "", mut_details_all)) %>%
    ungroup() %>%
    select(-mut_details)
  
  # Splitting the mut_details_all column dynamically
  vcf_data <- vcf_data %>%
    mutate(split_count = str_count(mut_details_all, "\\|") + 1)
  
  max_splits <- max(vcf_data$split_count, na.rm = TRUE)
  
  # Prepare column names based on the dynamically calculated maximum split count
  base_names <- c("symbol", "Gene", "Consequence")
  num_repeats <- ceiling(max_splits / length(base_names))
  full_names <- unlist(sapply(1:num_repeats, function(i) {
    paste(base_names, i, sep = "_")
  }))
  
  # Separate the mut_details_all column
  vcf_data <- vcf_data %>%
    separate(mut_details_all, into = full_names, sep = "\\|", remove = FALSE, extra = "merge", fill = "right") %>%
    mutate(across(contains("consequences"), ~gsub(",1", "", .x, fixed = TRUE)))
  
  # Reshape the extended columns
  vcf_data <- vcf_data %>%
    pivot_longer(
      cols = matches("^(symbol|Gene|Consequence)_\\d+$"),
      names_to = c(".value", "set"),
      names_sep = "_\\d"
    )
  
  return(vcf_data)
}

# Process each file and combine the consequences====================================================================
all_consequences <- data.frame()
for (vcf_file_path in vcf_files) {
  # Read VCF data, skipping comments
  vcf_data <- read.delim(vcf_file_path, header = FALSE, stringsAsFactors = FALSE, comment.char = "#")
  
  # Apply the predefined processing function
  processed_data <- process_vcf_data(vcf_data, column_names)  # Ensure this function and column_names are correctly defined
  
  # Reformat the data frame to variance count matrix------------------------------------------------
  consequence_data <- processed_data %>%
    select(chromosome, location, symbol, Gene, Consequence) %>%
    mutate(
      chromosome = str_replace_all(chromosome, "chr", ""),
      chromosome = paste0(chromosome, ":"),
      Location = paste0(chromosome, location)
    ) %>%
    select(-chromosome, -location) %>%
    select(Location, Gene, Consequence, symbol) %>%
    mutate(sample_id = sub(".somatic_pass", "", tools::file_path_sans_ext(basename(vcf_file_path)))) %>%
    mutate(Consequence = str_replace_all(Consequence, ",1", "")) %>%
    filter(!is.na(Gene))
  
  # Accumulate results
  all_consequences <- bind_rows(all_consequences, consequence_data)
}

# Save the combined data to a file
write.csv(all_consequences, paste0(output_dir, "/combined_consequences.csv"), row.names = FALSE)

# End of script
#===========================================================================================================================
