#==========================================================================#
#==========================================================================#
######                                                                ######
######  Script convert the output of CANVAS to the input of liftover  ######
######                        hg38 -> hg19                            ######
#==========================================================================#
#==========================================================================#

# Author: Lianlian Wu
# Date: 2024-09-02

setwd("~/Genomic_analysis")

####################################################
####   Source required libraries & load data   ####
####################################################

library(stringr)
library(dplyr)
library(tidyr)

# Specify the source_dir containing the VCF files
source_dir <- "./data/0_call_BE_driver_genes/copy_number/Genomics_England"

# List all VCF files
vcf_files <- list.files(path = source_dir, pattern = "somatic.SV_Canvas_only_pass\\.vcf$", full.names = TRUE, recursive = TRUE)

#print(vcf_files)

####################################################
process_vcf <- function(file_path) {
  lines <- readLines(file_path)
  
  # Find the line with the overall ploidy
  ploidy_line_index <- grep("^##OverallPloidy=", lines)
  ploidy_line <- lines[ploidy_line_index]
  ploidy_value <- sub("^##OverallPloidy=", "", ploidy_line)
  print(ploidy_value)
  
  vcf_data <- read.delim(file_path, header = FALSE, stringsAsFactors = FALSE, comment.char = "#")
  # vcf_data <- read.delim(vcf_files[1], header = FALSE, stringsAsFactors = FALSE, comment.char = "#")
  
  # Specify the new column names
  column_names <- c("Chromosome", "pos", "id", "ref", "alt", "qual", "filter", "info", "format", "proband_normal", "proband_tumor")
  names(vcf_data) <- column_names[seq_along(vcf_data)]
  
  output_file_path <- str_replace(file_path, "somatic.SV_Canvas_only_pass\\.vcf$", "somatic.SV_Canvas_only_pass_cleaned.bed")
  print(output_file_path)

  sample_id <- str_extract(basename(file_path), "^[^.]+")
  print(sample_id)

  # Clean vcf data
  vcf_data_ <- vcf_data %>%
    separate(id, into = c("Canvas", "type", "Chromosome", "range"), sep = ":") %>%
    separate(range, into = c("Start.Position", "End.Position"), sep = "-") %>%
    rowwise() %>%
    mutate(
      CN = {
        format_list <- strsplit(format, ":", fixed = TRUE)[[1]]
        tumor_list <- strsplit(proband_tumor, ":", fixed = TRUE)[[1]]
        cn_index <- match("CN", format_list)
        as.numeric(tumor_list[cn_index])
        },
      Major_CN = {
        format_list <- strsplit(format, ":", fixed = TRUE)[[1]]
        tumor_list <- strsplit(proband_tumor, ":", fixed = TRUE)[[1]]
        Major_index <- match("MCC", format_list)
        # tumor_list[Major_index]
        ifelse(tumor_list[Major_index]==".", 0, as.numeric(tumor_list[Major_index]))
      },
      sample_id = sample_id
    ) %>%
    ungroup()
  
  # Caculate CNA
  cleaned_data <- vcf_data_ %>%
    mutate(
      `Start.Position` = as.numeric(`Start.Position`),
      `End.Position` = as.numeric(`End.Position`),
      CN = as.numeric(CN),
      Major_CN = as.numeric(Major_CN),
      ploidy = ploidy_value,
      ploidy = as.numeric(ploidy),
      `Num Markers` = ifelse(is.na(`Start.Position`) | is.na(`End.Position`), NA,
        as.integer((`End.Position` - `Start.Position`) / 1e+05)),
      Seg.CN = log2(CN) - 1) %>%
      #Seg.CN = log2(CN/ploidy) - 1) %>%
    select(sample_id, Chromosome, 'Start.Position', 'End.Position', 'Num Markers', 'Seg.CN', ploidy, CN, Major_CN) %>%
    mutate(`Seg.CN` = as.numeric(`Seg.CN`))
  
  write.table(cleaned_data, 
              file = output_file_path, 
              quote = FALSE, 
              sep = "\t", 
              row.names = FALSE, 
              col.names = TRUE)
  
  return(cleaned_data)
}

# Process each file and combine results
combined_vcf_data <- lapply(vcf_files, process_vcf)
