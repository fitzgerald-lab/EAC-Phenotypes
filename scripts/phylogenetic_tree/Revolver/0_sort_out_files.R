#=====================================================================#
#=====================================================================#
######                                                           
######  Script sort out files and storage files by patient names  
######                                                           
#=====================================================================#
#=====================================================================#

# Author: Lianlian Wu
# Date: 2025-05-05

# Set your working directory
folder_path <- "/Users/wu04/Library/CloudStorage/OneDrive-UniversityofCambridge/projects/BO_gene_list/results/trees/Trees_structure/20250110_BO_green"

# List all files in the directory (excluding folders)
all_files <- list.files(folder_path, full.names = TRUE)
file_names <- list.files(folder_path)

# Filter only files (not directories)
all_files <- all_files[file.info(all_files)$isdir == FALSE]
file_names <- file_names[file.info(all_files)$isdir == FALSE]

# Extract the first 6 characters of each file name
prefixes <- substr(basename(all_files), 1, 6)

# Group and move files
for (prefix in unique(prefixes)) {
  # Create subfolder
  subfolder_path <- file.path(folder_path, prefix)
  if (!dir.exists(subfolder_path)) {
    dir.create(subfolder_path)
  }
  
  # Find files with the matching prefix
  matching_files <- all_files[prefixes == prefix]
  
  # Move each file into the subfolder
  for (file_path in matching_files) {
    file.rename(from = file_path, to = file.path(subfolder_path, basename(file_path)))
  }
}
