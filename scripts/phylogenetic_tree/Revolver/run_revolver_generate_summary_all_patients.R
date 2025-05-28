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
# devtools::install_github("caravagnalab/revolver")
require(revolver)
library(mtree)
library(ctree)
library(easypar)

library(dplyr)
library(tidyr)
library(stringr)
library(tidyverse)

library(ggdendro)
library(dendextend)
library(ggdendro)

library(ggplot2)

# Get full paths of all .R files in the directory
# r_files <- list.files(path = "./revolver-master/R/", pattern = "\\.R$", full.names = TRUE)

# Source each file
# sapply(r_files, source)

####################################################
#### Set path and load data
####################################################
result_dir <- "/mnt/scratchc/stlab-icgc/users/wu04/project/bo_gene_list/results/phylogenetic_tree/Revolver"

input_BO_df <- read.csv(file.path(result_dir, "combined_revolver_input_BO.csv")) %>% mutate(cluster = as.character(cluster)) 
# input_BO_df <- input_BO_df %>% filter(patientID != "PL_176")

input_Non_BO_df <- read.csv(file.path(result_dir, "combined_revolver_input_Non_BO.csv")) %>% mutate(cluster = as.character(cluster)) 

# ####################################################
# #### Remove duplicate variantIDs for the same patient
# ####################################################
remove_duplicate_variants <- function(input_df) {
    duplicate_rows <- input_df %>%
        group_by(variantID, patientID) %>%
        filter(n() > 1) %>%
        ungroup()

    deduplicated_df <- input_df %>%
        mutate(cluster = as.numeric(cluster)) %>%  # Ensure cluster is numeric for removing duplicates
        group_by(patientID, variantID) %>%
        slice_min(order_by = cluster, n = 1, with_ties = FALSE) %>%
        ungroup() %>%
        mutate(cluster = as.character(cluster))  # Convert cluster back to character

    print(paste("Original number of rows:", nrow(input_df)))
    print(paste("Number of duplicate rows:", nrow(duplicate_rows)))
    print(paste("Number of deduplicated rows:", nrow(deduplicated_df)))

    return(deduplicated_df)
}

input_Non_BO_df_deduplicated <- remove_duplicate_variants(input_Non_BO_df) 

input_BO_df_deduplicated <- remove_duplicate_variants(input_BO_df) 


combined_deduplicated_df <- rbind(input_BO_df_deduplicated %>% mutate(patientID = paste0("BE_", patientID)), input_Non_BO_df_deduplicated %>% mutate(patientID = paste0("Non-BE_", patientID)))

print(length(unique(combined_deduplicated_df$patientID)))

# ####################################################
# #### Create REVOLVER cohort and visualize
# ####################################################
# print("Starting to create REVOLVER cohort------------------------")
# # my_BO_cohort = revolver_cohort(input_BO_df_deduplicated, CCF_parser = CCF_parser, ONLY.DRIVER = FALSE, MIN.CLUSTER.SIZE = 0, annotation = "BE+ve patients" )
# # my_Non_BO_cohort = revolver_cohort(input_Non_BO_df_deduplicated, CCF_parser = CCF_parser, ONLY.DRIVER = FALSE, MIN.CLUSTER.SIZE = 0, annotation = "BE-ve patients" )
my_FULL_cohort = revolver_cohort(combined_deduplicated_df, CCF_parser = CCF_parser, ONLY.DRIVER = FALSE, MIN.CLUSTER.SIZE = 0, annotation = "FULL dataset" )

summary_data_all_patients_minCluster_0 <- Stats_cohort(my_FULL_cohort) # Check cohort summary statistics
write.csv(summary_data_all_patients_minCluster_0, file = file.path(result_dir, "summary_data_all_patients_minCluster_0.csv"), row.names = FALSE)
