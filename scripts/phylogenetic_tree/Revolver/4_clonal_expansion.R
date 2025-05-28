#=====================================================================#
#=====================================================================#
######                                                           
######      Script process pyclone output to calculate recent subclonal expansion
######                                                           
#=====================================================================#
#=====================================================================#

# Author: Lianlian Wu
# Date: 2025-05-22

setwd("/mnt/scratchc/stlab-icgc/users/wu04/project/bo_gene_list/scripts/phylogenetic_tree")

library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
source("functions_phy_tress.R")

# Load the data and get the list of patient directories
data_path <- "/mnt/scratchc/stlab-icgc/users/wu04/project/bo_gene_list/results/phylogenetic_tree/Trees_structure"
tree_files <- list.files(data_path, pattern = "tree_structure_models\\.csv$", full.names = TRUE, recursive = TRUE)

# Extract unique folder paths
tree_folders <- unique(dirname(tree_files))

# Exclude "PL_176" & "AH_412"
tree_folders <- tree_folders[!grepl("PL_176|AH_412", tree_folders)]


#======================================================================#
#####                                                               ####
#####  Calculate recent subclonal expansion score for each patient  ####
#####                                                               ####
#======================================================================#
#' Compute Recent Subclonal Expansion Score for a Patient
#'
#' This function reads a clone tree and associated CCF file for a given patient,
#' extracts CCF values for leaf nodes (terminal clusters), and calculates
#' the recent subclonal expansion score: the largest average CCF among the leaf clusters
#' across all tumour regions.
#'
#' @param patient_path Path to the directory for a given patient.
#' @return A numeric value representing the recent subclonal expansion score.
#' @examples
#' score <- compute_subclonal_expansion_score("/path/to/patient/folder")
compute_subclonal_expansion_score <- function(patient_path) {
  # Load tree structure (get the first .csv file that matches)
  tree_data_path <- list.files(
    patient_path,
    pattern = "tree_structure_models\\.csv$",
    full.names = TRUE,
    recursive = TRUE
  )[1]

  if (is.na(tree_data_path)) stop("Tree structure file not found.")

  # Load latest CCF file
  CCF_file_path <- get_latest_files(
    list.files(
      patient_path,
      pattern = "final_py_loci.*\\.csv$",
      full.names = TRUE,
      recursive = TRUE
    )
  )

  if (is.na(CCF_file_path)) stop("CCF file not found.")

  # Read files
  tree_df <- read.csv(tree_data_path, stringsAsFactors = FALSE)
  CCF_df <- read.csv(CCF_file_path, stringsAsFactors = FALSE)

  # Identify leaf nodes from tree (terminal clones)
  leaf_nodes <- tree_df %>%
    filter(is.term == TRUE) %>%
    pull(lab)

  # Subset CCF values for leaf nodes
  CCF_leaf_df <- CCF_df %>%
    filter(cluster %in% leaf_nodes) %>%
    select(cluster, contains("CCF"))

  # Compute average CCF per cluster per sample
  avg_ccf_df <- CCF_leaf_df %>%
    group_by(cluster) %>%
    summarise(across(ends_with("CCF"), mean, na.rm = TRUE), .groups = "drop")

  # Calculate the recent subclonal expansion score:
  # the maximum mean CCF across all regions among all leaf clusters
  Recent_subclonal_expansion_score <- avg_ccf_df %>%
    pivot_longer(
      cols = -cluster,
      names_to = "sample",
      values_to = "mean_ccf"
    ) %>%
    slice_max(order_by = mean_ccf, n = 1) %>%
    pull(mean_ccf)

  return(Recent_subclonal_expansion_score)
}

# Apply to multiple patients
subclonal_scores_df <- lapply(tree_folders, function(patient_path) {
  patient_id <- basename(patient_path)

  # Use tryCatch to handle errors gracefully per patient
  score <- tryCatch({
    compute_subclonal_expansion_score(patient_path)
  }, error = function(e) {
    message("❌ Error for ", patient_id, ": ", e$message)
    return(NA)  # Return NA if there was an error
  })

  data.frame(
    patient_id = patient_id,
    subclonal_expansion_score = score,
    stringsAsFactors = FALSE
  )
}) %>% do.call(rbind, .)

# Preview result
print(subclonal_scores_df)
write.csv(subclonal_scores_df, file.path(data_path, "subclonal_expansion_scores.csv"), 
          row.names = FALSE)


#======================================================================#
#####                                                               ####
#####                Plot the score between groups                  ####
#####                                                               ####
#======================================================================#
phenotype_data <- read.csv("/mnt/scratchc/stlab-icgc/users/wu04/project/bo_gene_list/data/WES/20250126_Stages_Leanne.csv")
phenotype_data$OCCAMS_ID <- gsub("/", "_", phenotype_data$OCCAMS_ID)

subclonal_scores_df <- read.csv(file.path(data_path, "subclonal_expansion_scores.csv"))
RSE_df <- subclonal_scores_df %>% left_join(phenotype_data, by = c("patient_id" = "OCCAMS_ID")) %>%
  select(patient_id, subclonal_expansion_score, group) %>%
  mutate(group = ifelse(group == "pos", "BE+ve", "BE-ve"))
RSE_df$group <- factor(RSE_df$group, levels = c("BE+ve", "BE-ve"))

# Perform Wilcoxon test for RSE values between the two groups
wilcox_test <- wilcox.test(subclonal_expansion_score ~ group, data = RSE_df)
p_valule <- round(wilcox_test$p.value, 2)

# Create box plot with scatters
ggplot(RSE_df, aes(x = group, y = subclonal_expansion_score, fill = group)) +
    geom_boxplot() +
    geom_jitter(width = 0.15, height = 0, alpha = 0.5) +
    labs(x = "", y = "Recent subclonal expansion score") +
    ggtitle("") +
    theme_minimal() +
    theme(
      panel.background = element_blank(),
      plot.background = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.text = element_text(color = "black", size = 16),
      axis.title = element_text(color = "black", size = 18),
      legend.position = "none"
    ) +
    scale_y_continuous(breaks = seq(0, 100, 20), limits = c(0, 106)) +
    scale_fill_manual(values = c("BE+ve" = "lightgreen", "BE-ve" = "lightblue")) +
    annotate("text", x = 1.5, y = 105, label = paste("p-value:", p_valule), hjust = 0.5, vjust = 0.5, size = 5)

ggsave(file.path(data_path, "plot", paste0("subclonal_expansion_score.jpg")), width = 4, height = 6, dpi = 300)



