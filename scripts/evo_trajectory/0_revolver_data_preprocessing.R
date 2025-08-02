#=====================================================================#
#=====================================================================#
######                                                           
######      Script to Process PyClone Output for Revolver Input   
######       
#=====================================================================#
#=====================================================================#

# Author: Lianlian Wu
# Date: 2025-05-07

# Description:
# This script processes PyClone output files to generate input data for 
# the Revolver package. It includes functions for identifying clonal 
# clusters, formatting data, and handling manual adjustments.

# Set working directory
setwd("/path/to/project/scripts/phylogenetic_tree")

####################################################
#### Load Required Libraries and Functions 
####################################################
library(dplyr)
library(tidyr)
library(stringr)
source("functions_phy_tress.R") # Custom functions for phylogenetic tree analysis

####################################################
#### Define Paths and Load Data
####################################################
# Define the directory containing PyClone output files
data_dir <- "/path/to/project/results/phylogenetic_tree/Trees_structure"

# List all relevant files in the directory
file_list <- list.files(data_dir, pattern = "*_final_py_loci.csv$", full.names = TRUE, recursive = TRUE)

# Define the directory for saving results
result_dir <- "/path/to/project/results/phylogenetic_tree/Revolver"

# Filter the file list to keep only the most recent file for each patient
clean_file_list <- get_latest_files(file_list)

####################################################
#### Process the data into required format for Revolver
####################################################
# The function identifying the dominant cluster for each sample and assigning the clonal flag
assign_clonality_flag <- function(data_df) {
    # Step 1: Identify CCF columns
    ccf_cols <- grep("CCF", names(data_df), value = TRUE)
    print(paste0("CCF columns: ", paste(ccf_cols, collapse = ", ")))
    # Step 2: Compute mean CCF per cluster per sample
    avg_ccf_per_cluster <- data_df %>%
        group_by(cluster) %>%
        summarise(across(all_of(ccf_cols), mean, na.rm = TRUE), .groups = "drop")
    
    # Step 3: Find the cluster with highest avg CCF for each sample
    highest_avg_cluster_per_sample <- pivot_longer(
        avg_ccf_per_cluster,
        cols = all_of(ccf_cols),
        names_to = "sample",
        values_to = "mean_ccf"
    ) %>%
        group_by(sample) %>%
        slice_max(order_by = mean_ccf, n = 1, with_ties = FALSE) %>%
        ungroup()
    
    # Print results
    print(highest_avg_cluster_per_sample)
    
    # Step 4: Check if all top clusters are the same
    if (length(unique(highest_avg_cluster_per_sample$cluster)) == 1) {
        # Consistent dominant cluster
        dominant_cluster <- unique(highest_avg_cluster_per_sample$cluster)
        message("✔ Dominant cluster is consistent: cluster ", dominant_cluster)
        
        data_df <- data_df %>%
            mutate(is.clonal = ifelse(cluster == dominant_cluster, "True", "False"))
        
    } else {
        # Inconsistent top clusters, print results
        message("❌ Dominant cluster is not consistent across samples. Here's the summary:")
        print(highest_avg_cluster_per_sample)
        # Optionally, set is.clonal to NA or "False" for all
        data_df <- data_df %>%
            mutate(is.clonal = NA)
    }
    return(data_df)
}

process_single_dataframe_for_revolver <- function(data_df_path) {
    print(paste0("Processing file: ", data_df_path))
    # Read the data
    data_df <- read.csv(data_df_path)
    data_df <- assign_clonality_flag(data_df)
    # Mutate and select required columns and change column formats
    data_df <- data_df %>% 
        mutate(
            # is.clonal = ifelse(cluster == 1, "True", "False"),
                     cluster = as.character(cluster),
                     patientID = str_extract(mutation_id, "^[A-Z]{2}_[0-9]{3}")
                     ) %>%
        select(Misc = mutation_id, patientID, variantID = gene_name, cluster, is.driver = Driver, 
                     is.clonal, 
                     contains("CCF"))
    
    # Merge the CCF values of multiple samples into one column into the example format: R1:0.86;R2:1;R3:1
    data_df_2 <- data_df %>%
        mutate(across(contains("CCF"), as.numeric) / 100) %>%
        rowwise() %>%
        mutate(CCF = paste(across(contains("CCF"), ~ paste0(str_remove(cur_column(), "\\.CCF"), ":", .x)),
                                            collapse = ";"
        )) %>%
        ungroup() %>%
        select(-contains(".CCF")) 
    
    return(data_df_2)
}

# Call the function with clean_file_list
processed_dfs <- lapply(clean_file_list, process_single_dataframe_for_revolver)
# Combine all processed dataframes into one
combined_df <- do.call(rbind, processed_dfs)

# Print the patients without clonal clusters flaged
print(paste0("Patients without clonal clusters flagged: ", 
             paste(unique(combined_df$patientID[is.na(combined_df$is.clonal)]), collapse = ", ")))

# Manually check and asign the clonal falg=======================================================
# Define manual clonal flags for specific patients
assign_manual_clonal_flag <- function(df, manual_flag) {
    df %>%
        mutate(is.clonal = ifelse(
            patientID %in% names(manual_flag) & 
                cluster == manual_flag[patientID],
            "True",
            "False"
        ),
        is.clonal = as.logical(is.clonal))
}

# Apply manual clonal flag assignments
manual_clonal_flag <- list(
    "Patient1" = "Cluster1",
    "Patient2" = "Cluster2"
    # Add more patient-cluster mappings as needed
)

combined_df <- assign_manual_clonal_flag(combined_df, manual_clonal_flag)

# Display the first few rows of the processed dataframe
head(combined_df)

# Print the number of unique patients in the dataset
print(paste0("Number of unique patients: ", length(unique(combined_df$patientID))))

# Save the processed dataframe to a CSV file
output_file <- file.path(result_dir, "combined_revolver_input.csv")
write.csv(combined_df, file = output_file, row.names = FALSE)
message("✔ Processed data saved to: ", output_file)
#=====================================================================#
# End of script