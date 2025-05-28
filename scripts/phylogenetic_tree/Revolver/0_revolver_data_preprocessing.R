#=====================================================================#
#=====================================================================#
######                                                           
######      Script process pyclone output to Revolver input   
######       
#=====================================================================#
#=====================================================================#

# Author: Lianlian Wu
# Date: 2025-05-07

setwd("/mnt/scratchc/stlab-icgc/users/wu04/project/bo_gene_list/scripts/phylogenetic_tree")

####################################################
#### Source required functions & load libraries 
####################################################
library(dplyr)
library(tidyr)
library(stringr)
source("functions_phy_tress.R")

####################################################
#### Set path and load data
####################################################
# data_dir <- "/mnt/scratchc/stlab-icgc/users/wu04/project/bo_gene_list/results/phylogenetic_tree/Trees_structure/20250110_BO_green"
data_dir <- "/mnt/scratchc/stlab-icgc/users/wu04/project/bo_gene_list/results/phylogenetic_tree/Trees_structure/20250110_Non_BO_blue"
file_list <- list.files(data_dir, pattern = "*_final_py_loci.csv$", full.names = TRUE, recursive = TRUE)
result_dir <- "/mnt/scratchc/stlab-icgc/users/wu04/project/bo_gene_list/results/phylogenetic_tree/Revolver"

clean_file_list <- get_latest_files(file_list) # keep the most recent one for each patient

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
manual_clonal_flag_BO <- c(
  "AH_135" = 1,
  "AH_412" = 1, # The dominant cluster can also be 2, the CNAs results could be wrong
  "ED_136" = 1,
  "ED_143" = 1,
  "NT_302" = 1,
  "PL_176" = 1, # No clonal CCF in PL_176_D_T, the CNAs results could be wrong
  "QE_077" = 1,
  "QE_130" = 1,
  "QE_173" = 1
)

manual_clonal_flag_Non_BO <- c(
  "ED_139" = 1,
  "NT_139" = 1,
  "NT_314" = 1,
  "NT_350" = 1,
)

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

# assign_manual_clonal_flag
# combined_df <- assign_manual_clonal_flag(combined_df, manual_clonal_flag_BO)
combined_df <- assign_manual_clonal_flag(combined_df, manual_clonal_flag_Non_BO)

head(combined_df)
print(length(unique(combined_df$patientID)))
# write.csv(combined_df, file = file.path(result_dir, "combined_revolver_input.csv"), row.names = FALSE)
# write.csv(combined_df, file = file.path(result_dir, "combined_revolver_input_BO.csv"), row.names = FALSE)

