#================================================================#
#================================================================#
######                                                      ######
######   Script exclude 6 GM and combine 205 IM samples     ######
######                                                      ######
#================================================================#
#================================================================#

# Author: Lianlian Wu
# Date: 2024-04-30

library(dplyr)

# Path to the .seg files
file1_path <- "/scratchb/stlab-icgc/users/wu04/project/bo_gene_list/data/aks_trio_ross_noora/cna/211_cna_for_gistic_2.seg"

seg_data1 <- read.table(file1_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#seg_data2 <- read.table(file2_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

#combined_seg_data <- bind_rows(seg_data1, seg_data2)

# Exclude six GM sample from the cohort analysis 
excluded_samples <- c("LP6008269-DNA_H01", "SLX-21762.idtUDP0094", "SLX-21762.idtUDP0096", "SLX-22426.idtUDP0109", "SLX-22426.idtUDP0119", "SLX-22427.idtUDP0130")

filtered_seg_data <- seg_data1 %>%
  filter(!grepl(paste(excluded_samples, collapse="|"), Sample))

# Two samples with fragmentized CNA from ASCAT
#excluded_samples_2 <- c("SLX-18929.UDP0028", "LP6008341-DNA_A03")

#filtered_seg_data <- filtered_seg_data %>%
#  filter(!grepl(paste(excluded_samples_2, collapse="|"), Sample))

# Save the filtered .seg file
output_file_path <- "/scratchb/stlab-icgc/users/wu04/project/bo_gene_list/data/aks_trio_ross_noora/cna/"
#write.table(combined_seg_data, file = paste0(output_file_path, '290_samples_cna.seg'), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
write.table(filtered_seg_data, file = paste0(output_file_path, '205_samples_cna_excludeGM_2.seg'), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

print(paste(".seg files saved at:", output_file_path))

