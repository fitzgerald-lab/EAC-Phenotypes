#=====================================================================#
#=====================================================================#
######                                                           
######  Script calculate mutational ITH per sample (# subclonal SNVs / # total SNVs)  
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
library(ggplot2)
source("functions_phy_tress.R")

# Set file_path and load data
data_dir <- "./data/phylogenetic_tree/Trees_structure"
file_dir_BO <- file.path(data_dir, "20250110_BO_green")
file_dir_NonBO <- file.path(data_dir, "20250110_Non_BO_blue")

# Get the list of files of BE+ve
file_list_BO <- list.files(file_dir_BO, pattern = "*final_py_loci.csv$", full.names = TRUE, recursive = TRUE)
file_list_NonBO <- list.files(file_dir_NonBO, pattern = "*final_py_loci.csv$", full.names = TRUE, recursive = TRUE)

#  keep the most recent one for each patient  
file_list <- file_list_NonBO
cleaned_file_list <- get_latest_files(file_list)

#=====================================================================# 
######  Calculate mutational ITH per patient (# subclonal SNVs / # total SNVs)
######  Save the results for BO and Non-BO patients separetely
#=====================================================================#

# Function to calculate ITH for a single patient
calculate_ith_from_file <- function(file_path) {
  df <- read.csv(file_path)
  
  subclonal_count <- nrow(df[df$cluster != 1, ])
  total_count <- nrow(df)
  ith <- ifelse(total_count > 0, subclonal_count / total_count, NA)
  
  patient_id <- substr(basename(file_path), 1, 6)
  
  return(data.frame(
    patient_id = patient_id,
    Subclonal_count = subclonal_count,
    Total_count = total_count,
    Mutational_ITH = ith,
    stringsAsFactors = FALSE
  ))
}

# Apply to all files
ith_results_Non_BO_df <- do.call(rbind, lapply(cleaned_file_list, calculate_ith_from_file))
ith_results_BO_df <- do.call(rbind, lapply(file_list_BO, calculate_ith_from_file))

# Save output
write.csv(ith_results_Non_BO_df, file.path(result_dir, "mutational_ITH_per_patient_Non_BO.csv"), row.names = FALSE)
write.csv(ith_results_BO_df, file.path(result_dir, "mutational_ITH_per_patient_BO.csv"), row.names = FALSE)

#============================================================================================================#
#=====================================================PLOT===================================================# 
#============================================================================================================#

BO_ith_result_df <- read.csv(file.path(result_dir, "mutational_ITH_per_patient_BO.csv")) %>% mutate(group = "BE+ve")
NonBO_ith_result_df <- read.csv(file.path(result_dir, "mutational_ITH_per_patient_Non_BO.csv")) %>% mutate(group = "BE-ve")
ith_result_df <- rbind(BO_ith_result_df, NonBO_ith_result_df)
# ith_result_df <- rbind(ith_results_Non_BO_df %>% mutate(group = "BE-ve"), ith_results_BO_df %>% mutate(group = "BE+ve"))
ith_result_df$group <- factor(ith_result_df$group, levels = c("BE+ve", "BE-ve"))

# Perform Wilcoxon test for ITH values between the two groups
wilcox_test <- wilcox.test(Mutational_ITH ~ group, data = ith_result_df)
p_valule <- round(wilcox_test$p.value, 2)

# Create box plot with scatters
ggplot(ith_result_df, aes(x = group, y = Mutational_ITH, fill = group)) +
    geom_boxplot() +
    geom_jitter(width = 0.15, height = 0, alpha = 0.5) +
    labs(x = "", y = "Mutational ITH") +
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
    scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1.1)) +
    scale_fill_manual(values = c("BE+ve" = "lightgreen", "BE-ve" = "lightblue")) +
    annotate("text", x = 1.5, y = 1.1, label = paste("p-value:", p_valule), hjust = 0.5, vjust = 0.5, size = 5)

ggsave(file.path(result_dir, "plot", paste0("mutational_ITH_boxplot.jpg")), width = 4, height = 6, dpi = 300)

#=====================================================================#
### End of script
