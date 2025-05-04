library(dplyr)

source_dir <- '/scratchb/stlab-icgc/users/wu04/project/bo_gene_list/results/'
file1 <- read.csv(paste0(source_dir, "gel/for_dNdScv/data_extracted_gel.csv"))
file2 <- read.csv(paste0(source_dir, "samples_211/combined_snv_for_dNdS.csv"))
file3 <- read.csv(paste0(source_dir, "samples_211/combined_indel_for_dNdS.csv"))

# Convert the 'pos' column in all files to character 
file1$pos <- as.character(file1$pos)
file2$pos <- as.character(file2$pos)
file3$pos <- as.character(file3$pos)

combined_data <- bind_rows(file1, file2, file3)

write.csv(combined_data, paste0(source_dir, "290_samples_combined_for_dNdS.csv"), row.names = FALSE)
