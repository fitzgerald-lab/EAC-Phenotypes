#=====================================================================#
# Script to analyze Revolver cohort summary data                     
# Includes metrics such as truncal SNVs, number of driver SNVs,      
# and number of clones with drivers.                                 
#=====================================================================#

# Author: Lianlian Wu
# Date: 2025-05-05

setwd("/path/to/project/scripts/evo_trajectory") # Adjust the path as needed

# Load required libraries
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)

# Set file_path and load data
result_dir <- "./results/phylogenetic_tree/Revolver"
data_df <- read.csv(file.path(result_dir, "summary_data_all_patients_minCluster_0.csv"))

#=====================================================================# 
######  Generate the group info for each patient
######  
#=====================================================================#
# print(colnames(data_df))
data_df$group <- sub("_.*", "", data_df$patientID)
data_df$group <- ifelse(data_df$group == "BE", "BE+ve", "BE-ve")
data_df$group <- factor(data_df$group, levels = c("BE+ve", "BE-ve"))
#=====================================================PLOT===================================================# 
#============================================================================================================#
# Perform Wilcoxon test 
# Perform Wilcoxon test for all numeric columns by group and generate a summary data frame
numeric_cols <- sapply(data_df, is.numeric)
numeric_colnames <- names(data_df)[numeric_cols]

wilcox_results <- lapply(numeric_colnames, function(col) {
  test <- wilcox.test(data_df[[col]] ~ data_df$group)
  data.frame(
    variable = col,
    p_value = test$p.value,
    statistic = test$statistic
  )
})

wilcox_df <- do.call(rbind, wilcox_results)
wilcox_df$p_value <- round(wilcox_df$p_value, 4)
print(wilcox_df)
write.csv(wilcox_df, file.path(result_dir, "wilcox_test_summary_data.csv"), row.names = FALSE)

# Create box plot with scatters
plot_box_with_stats <- function(data_df, column, breaks = breaks, limits = limits, group_col = "group", stat_df = wilcox_df) {
  # Convert column to symbol for tidy evaluation
  col_sym <- rlang::ensym(column)
  group_sym <- rlang::ensym(group_col)
  
  # Get maximum value for setting y-axis limit
  max_y <- ceiling(max(data_df[[as.character(col_sym)]], na.rm = TRUE))
  
  # Extract p-value from the stat_df
  p_val <- stat_df$p_value[stat_df$variable == as.character(col_sym)]
  p_label <- paste("p-value:", format(p_val, digits = 3, scientific = TRUE))

  # Plot
  ggplot(data_df, aes(x = !!group_sym, y = !!col_sym, fill = !!group_sym)) +
    geom_boxplot() +
    geom_jitter(width = 0.15, height = 0, alpha = 0.5) +
    labs(x = "", y = as.character(col_sym)) +
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
    scale_y_continuous(breaks = breaks, limits = limits) +
    scale_fill_manual(values = c("BE+ve" = "lightgreen", "BE-ve" = "lightblue")) +
    annotate("text", x = 1.5, y = 1.05*max_y, label = p_label, hjust = 0.5, size = 5)
}

# Plot for numMutations
plot_box_with_stats(data_df, column = numMutations, breaks = seq(0, 1200, 200), limits = c(0, 1200)) 
ggsave(file.path(result_dir, "plots", "boxplot_numMutations.jpg"), width = 8, height = 6, dpi = 300)

# Plot for numDriverMutations
plot_box_with_stats(data_df, column = numDriverMutations, breaks = seq(0, 20, 5), limits = c(0, 20))
ggsave(file.path(result_dir, "plots", "boxplot_numDriverMutations.jpg"), width = 4, height = 6, dpi = 300)

# Plot for numClonesWithDriver
max(data_df$numClonesWithDriver, na.rm = TRUE)
plot_box_with_stats(data_df, column = numClonesWithDriver, breaks = seq(0, 8, 2), limits = c(0, 8))
ggsave(file.path(result_dir, "plots", "boxplot_numClonesWithDriver.jpg"), width = 8, height = 6, dpi = 300)

# Plot for numTruncalMutations
max(data_df$numTruncalMutations, na.rm = TRUE)
plot_box_with_stats(data_df, column = numTruncalMutations, breaks = seq(0, 300, 50), limits = c(0, 300))
ggsave(file.path(result_dir, "plots", "boxplot_numTruncalMutations.jpg"), width = 8, height = 6, dpi = 300)

# End of script