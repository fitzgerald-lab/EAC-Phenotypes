#=====================================================================#
#=====================================================================#
######                                                           
######  Script to analyze the metrics including
######  the number of samples and clusters per patient
######  
#=====================================================================#
#=====================================================================#

# Author: Lianlian Wu
# Date: 2025-05-05

setwd("/path/to/project/scripts/evo_trajectory") # Adjust the path as needed

# Load required libraries
library(dplyr)
library(readr)
library(stringr)

# Set file_path and load data
result_dir <- "./results/phylogenetic_tree/Trees_structure"
file_dir_BO <- file.path(result_dir, "20250110_BO_green")
file_dir_NonBO <- file.path(result_dir, "20250110_Non_BO_blue")

# Get the list of files of BE+ve, and filter(patientID != "PL_176" & patientID != "AH_412"), 
# because these two patients don't have clonal CCF -- the CNA looks not right, probably from insufficient ASCAT QC
file_list_BO <- list.files(file_dir_BO, pattern = "*cluster_ccf.csv$", full.names = TRUE, recursive = TRUE)
file_list_BO <- file_list_BO[!grepl("PL_176|AH_412", file_list_BO)]

# Get the list of files of BE-ve
file_list_NonBO <- list.files(file_dir_NonBO, pattern = "*cluster_ccf.csv$", full.names = TRUE, recursive = TRUE)

#=====================================================================# 
######  There are multiple .final_py_loci.csv files in the list
######  The codes will keep the most recent one for each patient   
#=====================================================================#
# Extract patient ID and date from each path
file_list <- file_list_NonBO
# file_list <- file_list_BO
file_df <- data.frame(
  file_path = file_list,
  stringsAsFactors = FALSE
) %>%
  mutate(
    file_name = basename(file_path),
    patient_id = str_extract(file_name, "^[A-Z]{2}_[0-9]{3}"),
    file_date = str_extract(file_name, "\\d{4}-\\d{2}-\\d{2}"),
    file_date = as.Date(file_date)
  )

# Keep the latest file per patient
latest_files <- file_df %>%
  group_by(patient_id) %>%
  slice_max(order_by = file_date, n = 1, with_ties = FALSE) %>%
  ungroup()

# Result: only the most recent .csv per patient
cleaned_file_list <- latest_files$file_path

#=====================================================================# 
######  Calculate the number of samples and clusters per patient
######  Save the results for BO and Non-BO patients separetely
#=====================================================================#
# Function to calculate number of samples and clusters for a single patient
calculate_samples_clusters <- function(file_path) {
    test_df <- read.csv(file_path)
    n_samples <- length(unique(test_df$sample))
    n_clusters <- length(unique(test_df$cluster))
    patient_id <- substr(basename(file_path), 1, 6)
    
    return(list(n_samples = n_samples, n_clusters = n_clusters, patient_id = patient_id))
}

result_list <- lapply(cleaned_file_list, calculate_samples_clusters)
result_df <- do.call(rbind, lapply(result_list, function(x) {
    data.frame(
        patient_id = x$patient_id,
        n_samples = x$n_samples,
        n_clusters = x$n_clusters
    )
}))

# write.csv(result_df, file.path(result_dir, "number_of_samples_clusters_BO.csv"), row.names = FALSE)
write.csv(result_df, file.path(result_dir, "number_of_samples_clusters_Non_BO.csv"), row.names = FALSE)

#============================================================================================================#
#=====================================================PLOT===================================================# 
#============================================================================================================#
# Load required libraries
library(ggplot2)

BO_n_result_df <- read.csv(file.path(result_dir, "number_of_samples_clusters_BO.csv")) %>% mutate(group = "BE+ve")
NonBO_n_result_df <- read.csv(file.path(result_dir, "number_of_samples_clusters_Non_BO.csv")) %>% mutate(group = "BE-ve")
sum(BO_n_result_df$n_samples)
sum(NonBO_n_result_df$n_samples)

n_result_df <- rbind(BO_n_result_df, NonBO_n_result_df)
head(n_result_df)
write.csv(n_result_df, file.path(result_dir, "number_of_samples_clusters.csv"), row.names = FALSE)
print(length(unique(n_result_df$patient_id)))
print(sum(n_result_df$n_samples))
print(sum(n_result_df$n_clusters))
############################ N_SAMPLES############################
# Perform Wilcoxon test for ITH values between the two groups
wilcox_test <- wilcox.test(n_samples ~ group, data = n_result_df)
p_valule <- round(wilcox_test$p.value, 2)

n_result_df$group <- factor(n_result_df$group, levels = c("BE+ve", "BE-ve"))

# Create box plot with scatters
ggplot(n_result_df, aes(x = group, y = n_samples, fill = group)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, height = 0, alpha = 0.5) +
    labs(x = "", y = "Number of samples") +
    ggtitle("") +
    theme_minimal() +
    theme(
      panel.background = element_blank(),
      axis.line = element_line(color = "black"),
      axis.text = element_text(color = "black", size = 16),
      axis.title = element_text(color = "black", size = 18),
      legend.position = "none"
    ) +
    scale_y_continuous(breaks = seq(2, 7, 1), limits = c(2, 7)) +
    scale_fill_manual(values = c("BE+ve" = "lightgreen", "BE-ve" = "lightblue")) +
    annotate("text", x = 1.5, y = 7, label = paste("p-value:", p_valule), hjust = 0.5, vjust = 0.5, size = 5)

ggsave(file.path(result_dir, "plot", "number_samples_boxplot.jpg"), width = 8, height = 6, dpi = 300)

############################ N_CLUSTERS ############################
# Perform Wilcoxon test for ITH values between the two groups
wilcox_test <- wilcox.test(n_clusters ~ group, data = n_result_df)
p_valule <- round(wilcox_test$p.value, 2)

# Create box plot with scatters
ggplot(n_result_df, aes(x = group, y = n_clusters, fill = group)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, height = 0, alpha = 0.5) +
    labs(x = "", y = "Number of clusters") +
    ggtitle("") +
    theme_minimal() +
    theme(
      panel.background = element_blank(),
      axis.line = element_line(color = "black"),
      axis.text = element_text(color = "black", size = 16),
      axis.title = element_text(color = "black", size = 18),
      legend.position = "none"
    ) +
    scale_y_continuous(breaks = seq(2, 10, 1), limits = c(2, 11)) +
    scale_fill_manual(values = c("BE+ve" = "lightgreen", "BE-ve" = "lightblue")) +
    annotate("text", x = 1.5, y = 11, label = paste("p-value:", p_valule), hjust = 0.5, vjust = 0.5, size = 5)

ggsave(file.path(result_dir, "plot", "number_clusters_boxplot.jpg"), width = 8, height = 6, dpi = 300)

