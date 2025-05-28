#================================================================#
#================================================================#
######                                                      ######
######   Script count markers in segments across samples    ######
######                                                      ######
#================================================================#
#================================================================#

# Author: Lianlian Wu
# Date: 2024-04-30

# Load required libraries
library(dplyr)
library(readr)
library(stringr)

source_dir <- "/path/to/output/folder/of/CAVAMEN/"

# List all caveman adjusted CSV files
adjusted_files <- list.files(source_dir, pattern = "\\.copynumber\\.caveman_loged\\.csv$", full.names = TRUE, recursive = TRUE)

# Process each file
for (seg_path in adjusted_files) {
  # Construct the path for the corresponding marker file
  marker_path <- str_replace(seg_path, "caveman_loged\\.csv$", "txt")
  
  # Check if marker file exists
  if (file.exists(marker_path)) {
    # Load data
    segments <- read_csv(seg_path)
    markers <- read_delim(marker_path, delim = "\t", escape_double = FALSE, trim_ws = TRUE)
    
    # Define a function to count markers within a segment
    count_markers <- function(chromosome, start, end) {
      sum(markers$Chromosome == chromosome & markers$Position >= start & markers$Position <= end)
    }
    
    # Add the count of unique markers to the segments data frame
    segments <- segments %>%
      rowwise() %>%
      mutate(Unique_Markers = count_markers(Chromosome, Start, End))
    
    # Define the output path
    output_file_path <- str_replace(seg_path, "caveman_loged\\.csv$", "caveman_loged_with_marker.csv")
    
    # Write the updated data to a new CSV file
    write_csv(segments, output_file_path)
    print(paste("Processed and saved:", output_file_path))
  } else {
    print(paste("Marker file not found for:", seg_path))
  }
}

