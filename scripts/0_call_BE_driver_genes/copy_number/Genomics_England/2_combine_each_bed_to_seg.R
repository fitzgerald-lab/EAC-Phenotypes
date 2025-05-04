#==============================================================================#
#==============================================================================#
######                                                                    ######
######   Script combining CNA mapped with hg19 for downstream analysis    ######
######       -- 79 BE samples sequenced by Genomics England               ######
#==============================================================================#
#==============================================================================#

# Author: Lianlian Wu
# Date: 2024-09-02

setwd("/scratchb/stlab-icgc/users/wu04/project/bo_gene_list/scripts/gel/cna/")

####################################################
#### Source required functions & load libraries ####
####################################################

library(stringr)
library(dplyr)
library(tidyr)

# Specify the source_dir containing the VCF files
source_dir <- "/scratchb/stlab-icgc/users/wu04/project/bo_gene_list/data/gel/cna_canvas_pass_hg19/"

bed_files <- list.files(path = source_dir, pattern = "somatic.SV_Canvas_only_pass_cleaned_hg19\\.bed$", full.names = TRUE, recursive = TRUE)

process_bed <- function(file_path) {
  bed_data <- read.delim(file_path, header = FALSE, stringsAsFactors = FALSE)
  
  # Specify the new column names
  column_names <- c('Chromosome', 'Start.Position', 'End.Position', 'Num Markers', 'Seg.CN')
  names(bed_data) <- column_names[seq_along(bed_data)]
  
  sample_id <- sub(".somatic.SV_Canvas_only_pass_cleaned_hg19.bed", "", basename(file_path))
  
  # Clean vcf data
  bed_data <- bed_data %>%
    mutate(Sample = sample_id,
           Chromosome = as.character(gsub("chr", "", Chromosome, fixed = TRUE)),
           `Seg.CN` = as.numeric(`Seg.CN`),
           `Start.Position` = as.numeric(`Start.Position`),
           `End.Position` = as.numeric(`End.Position`),
           ) %>%
    select(Sample, Chromosome, 'Start.Position', 'End.Position', 'Num Markers', 'Seg.CN')
  
  return(bed_data)
}

combined_bed_data <- lapply(bed_files, process_bed) %>% bind_rows()

output_path <- paste0(source_dir, "combined_gel_cna.seg")
write.table(combined_bed_data, file = output_path, sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

print(paste("Combined .seg file saved at:", output_path))
