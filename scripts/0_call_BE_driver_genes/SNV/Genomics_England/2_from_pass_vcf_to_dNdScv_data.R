#===========================================================================================
#===========================================================================================      
######   Script generate dNdS input   
######   Input: .snp.pass.vep files          
######   Output: data_extracted_gel.csv                                                                    
#===========================================================================================
#===========================================================================================

# Author: Lianlian Wu
# Date: 2024-05-01

library(stringr)
library(dplyr)

# Specify the source_dir containing the VCF files
source_dir <- "/scratchb/stlab-icgc/users/wu04/project/bo_gene_list/data/gel/snv_hg19_filter_pass/"
output_dir <- "/scratchb/stlab-icgc/users/wu04/project/bo_gene_list/results/gel/for_dNdScv/"

# List all VCF files
vcf_files <- list.files(path = source_dir, pattern = "\\.somatic_lifted_hg19_pass\\.vcf$", full.names = TRUE)

process_vcf <- function(file_path) {
  vcf_data <- read.delim(file_path, header = FALSE, stringsAsFactors = FALSE, comment.char = "#")
  
  # Specify the new column names
  column_names <- c("chr", "pos", "snp_id", "ref", "mut", "V6", "filter", "mut_details", "indexs", "values_1", "values_2")
  
  # Ensure the column names match up with the existing names if the file has fewer columns
  names(vcf_data) <- column_names[seq_along(vcf_data)]
  
  # Extract the sample ID from the file name
  sample_id <- sub(".somatic_lifted_hg19_pass.vcf", "", basename(file_path))
  
  # Clean up the chromosome number and add sampleID
  vcf_data <- vcf_data %>%
    mutate(sampleID = sample_id, chr = str_remove(chr, "^chr")) %>%
    select(sampleID, chr, pos, ref, mut)
  
  return(vcf_data)
}

# Process each file and combine results
combined_vcf_data <- lapply(vcf_files, process_vcf) %>% bind_rows()

# Check the first few rows of the combined DataFrame
head(combined_vcf_data)

# View the structure of the DataFrame
str(combined_vcf_data)

# Save to CSV 
write.csv(combined_vcf_data, paste0(output_dir, "data_extracted_gel.csv"), row.names = FALSE)
