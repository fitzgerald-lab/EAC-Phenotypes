#===========================================================================================
#===========================================================================================      
######   Script generate indel consequences for dNdS input   
######   Input: .indel.pass.vep files          
######   Output: combined_indel_for_dNdS.csv                                                                    
#===========================================================================================
#===========================================================================================

# Author: Lianlian Wu
# Date: 2024-05-01

library(tidyr)
library(dplyr)
library(stringr)

# Specify the source_dir containing the VCF files
source_dir <- "/scratchb/stlab-icgc/users/wu04/project/bo_gene_list/data/aks_trio_ross_noora/indel/"
output_dir <- "/scratchb/stlab-icgc/users/wu04/project/bo_gene_list/results/samples_211/"

vep_files <- list.files(source_dir, pattern = "\\.indel.pass.vep$", full.names = TRUE)

process_vep_file <- function(file_path) {
  lines <- readLines(file_path)
  data_start_line <- which(grepl("^#Uploaded_variation", lines))-1
  vep_data <- read.delim(file_path, skip = data_start_line, header = TRUE, stringsAsFactors = FALSE)
  
  # Extract sample ID from the filename
  pattern_ <- "__pv.*__rg.grch37_g1k__al.bwa_mem__.indel.pass.vep"
  sample_id <- sub(pattern_, "", basename(file_path))
  
print(colnames(vep_data))
  # Modify data frame
  vep_data <- vep_data %>%
    mutate(sampleID = sample_id) %>%
    separate(X.Uploaded_variation, into = c("chr", "pos", "ref_mut"), sep = "_", extra = "merge", fill = "right") %>%
    separate(ref_mut, into = c("ref", "mut"), sep = "/", extra = "merge", fill = "right") %>%
    select(sampleID, chr, pos, ref, mut)
  
  return(vep_data)
}

all_vep_data <- lapply(vep_files, process_vep_file) %>%
  bind_rows()

write.csv(all_vep_data, paste0(output_dir, "combined_indel_for_dNdS.csv"), row.names = FALSE)

# End of script