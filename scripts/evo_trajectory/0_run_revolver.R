#=====================================================================#
#=====================================================================#
######                                                           
######      Script process pyclone output to Revolver input   
######                                                           
#=====================================================================#
#=====================================================================#

# Author: Lianlian Wu
# Date: 2025-05-07

setwd("./scripts/phylogenetic_tree")

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

# Get full paths of all .R files in the directory. Notable!! This can only be loaded after processing the cohort, and used in plotting functions, otherwise the argument type will be wrong.
# r_files <- list.files(path = "./revolver-master/R/", pattern = "\\.R$", full.names = TRUE)
# sapply(r_files, source)

####################################################
#### Set path and load data
####################################################
result_dir <- "./results/phylogenetic_tree/Revolver"

# BO cohort---------------
# input_BO_df <- read.csv(file.path(result_dir, "combined_revolver_input_BO.csv")) %>% mutate(cluster = as.character(cluster)) 
# print(length(unique(input_BO_df$patientID)))

# # Non-BO cohort---------------
# input_Non_BO_df <- read.csv(file.path(result_dir, "combined_revolver_input_Non_BO.csv")) %>% mutate(cluster = as.character(cluster)) 

########################################################
#### Remove duplicate variantIDs for the same patient
########################################################
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

# combined_revolver_input <- rbind(input_BO_df_deduplicated %>% mutate(group = "BE+ve"), input_Non_BO_df_deduplicated %>% mutate(group = "BE-ve"))
# write.csv(combined_revolver_input, file.path(result_dir, "combined_revolver_input.csv"), row.names = FALSE)

combined_deduplicated_df <- rbind(input_BO_df_deduplicated %>% mutate(patientID = paste0("BE_", patientID)), input_Non_BO_df_deduplicated %>% mutate(patientID = paste0("Non-BE_", patientID)))

# # print(length(unique(combined_deduplicated_df$patientID)))

# # # ####################################################
# # # #### Create REVOLVER cohort and visualize
# # # ####################################################
# # print("Starting to create REVOLVER cohort------------------------")
my_BO_cohort = revolver_cohort(input_BO_df_deduplicated, CCF_parser = CCF_parser, ONLY.DRIVER = FALSE, MIN.CLUSTER.SIZE = 0, annotation = "BE+ve patients" )
my_Non_BO_cohort = revolver_cohort(input_Non_BO_df_deduplicated, CCF_parser = CCF_parser, ONLY.DRIVER = FALSE, MIN.CLUSTER.SIZE = 0, annotation = "BE-ve patients" )
my_FULL_cohort = revolver_cohort(combined_deduplicated_df, CCF_parser = CCF_parser, ONLY.DRIVER = FALSE, MIN.CLUSTER.SIZE = 10, annotation = "FULL dataset" )


# print(my_cohort)
cohort_summary_plot <- plot(my_FULL_cohort)
# cohort_summary_BO_plot <- plot(my_BO_cohort)
# cohort_summary_Non_BO_plot <- plot(my_Non_BO_cohort)

ggsave(cohort_summary_plot, file = file.path(result_dir, "plots/cohort_full_summary_plot.png"), width = 10, height = 8, dpi = 300, units = "in")
ggsave(cohort_summary_BO_plot, file = file.path(result_dir, "cohort_summary_BO_plot.png"), width = 10, height = 8, dpi = 300, units = "in")
ggsave(cohort_summary_Non_BO_plot, file = file.path(result_dir, "cohort_summary_Non_BO_plot.png"), width = 10, height = 8, dpi = 300, units = "in")

###################################################
### Clonal/Subclonal drivers
###################################################
plot_drivers_clonality_with_order <- function(x, variant_order = NULL) {
  st <- Stats_drivers(x) %>% arrange(desc(numClonal), desc(numSubclonal))
  
  if (is.null(variant_order)) {
    variant_order <- st$variantID
  } else {
    st <- st %>% filter(variantID %in% variant_order)
  }

  st$numSubclonal <- -st$numSubclonal
  st <- st %>%
    select(variantID, numClonal, numSubclonal) %>%
    rename(Clonal = numClonal, Subclonal = numSubclonal) %>%
    reshape2::melt(id = "variantID")

  st$variantID <- factor(st$variantID, levels = variant_order)

  N <- length(x$patients)

  ggplot(st, aes(x = variantID, y = value, fill = variable)) + 
    geom_bar(stat = "identity", position = "identity") + 
    coord_flip(clip = "off") +
    scale_fill_manual(values = c(Clonal = "steelblue", Subclonal = "darkorange3")) +
    labs(
      title = paste("Driver burden"),
      y = paste0("Occurrences (n = ", N, " patients)"),
      x = "Driver",
      subtitle = paste(x$annotation)
    ) +
    guides(fill = guide_legend("")) +
    theme_minimal()
}

variant_order <- Stats_drivers(FULL_clone_jackknife) %>%
  arrange(numClonal, numSubclonal) %>%
  pull(variantID)

clonal_subclonal_full_cohort_plot <- plot_drivers_clonality_with_order(FULL_clone_jackknife, variant_order)
clonal_subclonal_BO_plot <- plot_drivers_clonality_with_order(my_BO_cohort, variant_order)
clonal_subclonal_Non_BO_plot <- plot_drivers_clonality_with_order(my_Non_BO_cohort, variant_order)
# plot_drivers_clonality(my_Non_BO_cohort)


# st <- Stats_drivers(Non_BO_clone_jackknife) %>% arrange(desc(numClonal), desc(numSubclonal))

ggsave(clonal_subclonal_full_cohort_plot, file = file.path(result_dir, "clonal_subclonal_full_cohort_plot.jpg"), width = 4, height = 8, dpi = 300, units = "in")
ggsave(clonal_subclonal_BO_plot, file = file.path(result_dir, "clonal_subclonal_BO_plot.jpg"), width = 4, height = 8, dpi = 300, units = "in")
ggsave(clonal_subclonal_Non_BO_plot, file = file.path(result_dir, "clonal_subclonal_Non_BO_plot.jpg"), width = 4, height = 8, dpi = 300, units = "in")


####################################################
#### Fit Trees, Fit via TL, REVOLVER clustering and Jackknife statistics
####################################################
print("Starting to process REVOLVER cohort------------------------")
# Non_BO cohort
print("Fit Trees----------------------------------------------------------")
Non_BO_clone_trees <- compute_clone_trees(my_Non_BO_cohort)
saveRDS(Non_BO_clone_trees, file.path(result_dir, "Non_BO_clone_trees.rds"))

print("Remove drivers with only one occurrence in the cohort--------------")
Stats_drivers_non_BO <- Stats_drivers(my_Non_BO_cohort) %>% filter(N_tot == 1)
NEW_Non_BO_cohort = remove_drivers(Non_BO_clone_trees, Stats_drivers_non_BO$variantID)
saveRDS(NEW_Non_BO_cohort, file.path(result_dir, "Non_BO_rm_one_occur_diver_cohort.rds"))

print("Fit via TL----------------------------------------------------------")
NEW_Non_BO_cohort <- readRDS(file.path(result_dir, "Non_BO_rm_one_occur_diver_cohort.rds"))
Non_BO_clone_fit <- revolver_fit(NEW_Non_BO_cohort)
saveRDS(Non_BO_clone_fit, file.path(result_dir, "Non_BO_clone_fit.rds"))

print("Clustering----------------------------------------------------------")
Non_BO_clone_cluster <- revolver_cluster(Non_BO_clone_fit)
saveRDS(Non_BO_clone_cluster, file.path(result_dir, "Non_BO_clone_cluster.rds"))

print("Jackknife statistics------------------------------------------------")
Non_BO_clone_cluster <- readRDS(file.path(result_dir, "Non_BO_clone_cluster.rds"))
Non_BO_clone_jackknife <- revolver_jackknife(Non_BO_clone_cluster)
saveRDS(Non_BO_clone_jackknife, "Non_BO_clone_jackknife.rds")
Load and continue if interrupted
Non_BO_clone_jackknife <- readRDS(file.path(result_dir, "Non_BO_all.rds"))

BO cohort
Stats_trees(my_cohort)  # Check tree status per patient
BO_clone_trees <- compute_clone_trees(my_BO_cohort)
# Remove drivers with only one occurrence in the cohort
Stats_drivers_BO <- Stats_drivers(my_BO_cohort) %>% filter(N_tot == 1)
NEW_BO_cohort = remove_drivers(BO_clone_trees, Stats_drivers_BO$variantID)
BO_clone_fit <- revolver_fit(NEW_BO_cohort)
BO_clone_cluster <- revolver_cluster(BO_clone_fit)
BO_clone_jackknife <- revolver_jackknife(BO_clone_cluster)

saveRDS(BO_clone_jackknife, "BO_cohort_late_stages.rds")
# Load and continue if interrupted
# BO_clone_jackknife <- readRDS("BO_cohort_late_stages.rds")
options(width = 200)
str(my_FULL_cohort)

# Full cohort
print("Fit Trees----------------------------------------------------------")

#' Safely compute clone trees for all patients in a REVOLVER cohort
#'
#' @param cohort A REVOLVER cohort object.
#' @param overwrite Logical; if TRUE, overwrite existing trees (default = FALSE).
#' @param ... Additional arguments passed to ctree::ctrees().
#'
#' @return The same cohort object with $phylogenies populated.
#' @export

safe_compute_trees <- function(cohort, overwrite = TRUE) {
  # Get all patient IDs from the cohort
  all_patients <- Stats_cohort(cohort)$patientID

  # Initialize results
  phylo_list <- list()
  failed_patients <- c()

  # Loop through patients
  for (patient in all_patients) {
    message(">>> Processing patient: ", patient)

    # Skip if tree already exists and overwrite is FALSE
    if (!overwrite && has_patient_trees(cohort, patient)) {
      message(" - Skipping ", patient, " (already has tree, overwrite = FALSE)")
      phylo_list[[patient]] <- cohort$phylogenies[[patient]]
      next
    }

    # Safely compute tree for the patient
    tryCatch({
      result <- ctree::ctrees(
        CCF_clusters(cohort, patient),
        Drivers(cohort, patient),
        Samples(cohort, patient),
        patient = patient
      )
      phylo_list[[patient]] <- result
      message(" - Successfully computed tree.")
    }, error = function(e) {
      message(" !! Failed for ", patient, ": ", e$message)
      failed_patients <<- c(failed_patients, patient)
    })
  }

  # Assign successful results to the cohort
  cohort$phylogenies <- phylo_list

  # Summary
  message("\n✅ Tree construction completed.")
  message("✔ Success: ", length(phylo_list), " patients.")

  # Print all patient IDs with successful trees
  message("✔ Success: ", length(cohort$phylogenies), " patients.")
  message("Patient IDs: ", paste(names(cohort$phylogenies), collapse = ", "))

  return(cohort)
}

# # Example call
FULL_clone_trees <- safe_compute_trees(my_FULL_cohort, overwrite = TRUE)

FULL_clone_trees <- compute_clone_trees(my_FULL_cohort)
saveRDS(FULL_clone_trees, file.path(result_dir, "FULL_clone_trees_10.rds"))

print("Remove drivers with only one occurrence in the cohort--------------")
# FULL_clone_trees <- readRDS(file.path(result_dir, "FULL_clone_trees_10.rds"))
Stats_drivers_FULL <- Stats_drivers(FULL_clone_trees) %>% filter(N_tot == 1)
NEW_FULL_cohort = remove_drivers(FULL_clone_trees, Stats_drivers_FULL$variantID)
saveRDS(NEW_FULL_cohort, file.path(result_dir, "FULL_rm_one_occur_diver_cohort_10.rds"))

print("Fit via TL----------------------------------------------------------")
NEW_FULL_cohort <- readRDS(file.path(result_dir, "FULL_rm_one_occur_diver_cohort_10.rds"))
FULL_clone_fit <- revolver_fit(NEW_FULL_cohort)
saveRDS(FULL_clone_fit, file.path(result_dir, "FULL_clone_fit_10.rds"))

print("Clustering----------------------------------------------------------")
FULL_clone_fit <- readRDS(file.path(result_dir, "FULL_clone_fit_10.rds"))
FULL_clone_cluster <- revolver_cluster(FULL_clone_fit)
saveRDS(FULL_clone_cluster, file.path(result_dir, "FULL_clone_cluster_10.rds"))

print("Jackknife statistics------------------------------------------------")
FULL_clone_cluster <- readRDS(file.path(result_dir, "FULL_clone_cluster_10.rds"))
FULL_clone_jackknife <- revolver_jackknife(FULL_clone_cluster)

print("Save the REVOLVER full cohort----------------------------------------")
saveRDS(FULL_clone_jackknife, file.path(result_dir, "FULL_clone_jackknife_10.rds"))
print("The REVOLVER full cohort has been built and saved as Full_cohort.rds")
FULL_clone_jackknife <- readRDS(file.path(result_dir, "FULL_clone_jackknife_10.rds"))

####################################################
#### Plot the graph-alike summary statistics for both phenotypes
####################################################
# Plot the  graph-alike summary statistics for the cohort drivers
# Non_BO_driver_graph <- plot_drivers_graph(Non_BO_clone_jackknife)
# ggsave(Non_BO_driver_graph, file = file.path(result_dir, "driver_graph_Non_BO.jpg"), width = 14, height = 12, dpi = 200, units = "in")

# BO_driver_graph <- plot_drivers_graph(BO_clone_jackknife)
# ggsave(BO_driver_graph, file = file.path(result_dir, "driver_graph_BO.jpg"), width = 14, height = 12, dpi = 200, units = "in")

####################################################
#### Plot the heatmaps of REVOLVER"s clusters -- using FULL cohort
####################################################
# trajectories_drivers_heatmap <- plot_clusters(Non_BO_clone_jackknife) 
# ggsave(trajectories_drivers_heatmap, file = file.path(result_dir, "plots/trajectories_drivers_heatmap_Non_BO.jpg"), width = 8, height = 6, dpi = 200, units = "in")

# trajectories_drivers_heatmap_BO <- plot_clusters(BO_clone_jackknife) 
# ggsave(trajectories_drivers_heatmap_BO, file = file.path(result_dir, "plots/trajectories_drivers_heatmap_BO.jpg"), width = 8, height = 6, dpi = 200, units = "in")

# trajectories_drivers_heatmap_FULL <- plot_clusters(FULL_clone_jackknife, cutoff_drivers = 10, cutoff_trajectories =8) 
trajectories_drivers_heatmap_FULL <- plot_clusters(FULL_clone_jackknife, cutoff_drivers = 0, cutoff_trajectories = 2) 
ggsave(trajectories_drivers_heatmap_FULL, file = file.path(result_dir, "plots/trajectories_drivers_heatmap_full_cohort_filtered_cutoff_drivers_0_cutoff_trajectories_2.jpg"), width = 8, height = 20, dpi = 300, units = "in")

####################################################
#### Plot dendrogram clusters plot
####################################################
# plot_dendrogram(FULL_clone_jackknife)
class(FULL_clone_jackknife$cluster$fits$hc)
# str(FULL_clone_jackknife$cluster$fits$hc)

  
# Dendrogram - it gives the ordering of the patients which are displayed on the x-axcis
hc = FULL_clone_jackknife$cluster$fits$hc

# Patient ordering - from the dedrogram, thi defines the levels of the factors used
factors_patient_level = hc$order.lab

# Get colors for the clusters
clusters_colors = get_cluster_colors(FULL_clone_jackknife, distinct_palette_few)

# Assign the colors following factors_patient_level
patients_factors_colors = sapply(factors_patient_level,
                                function(y)
                                  clusters_colors[Cluster(FULL_clone_jackknife, y) %>% pull(cluster)])
names(patients_factors_colors) = factors_patient_level

bars_separation = Cluster(FULL_clone_jackknife, factors_patient_level)
bars_separation$cluster = factor(bars_separation$cluster, levels = unique(bars_separation$cluster))
bars_separation = bars_separation %>% pull(cluster) %>% table %>% cumsum + 0.5

# Number of clusters
nclusters = Cluster(FULL_clone_jackknife) %>% pull(cluster) %>% unique %>% length

x = FULL_clone_jackknife
# Dendrogram plot
hc_ <- as.hclust(FULL_clone_jackknife$cluster$fits$hc)
cluster_plot <- ggdendro::ggdendrogram(hc_,
                      rotate = FALSE, 
                      size = 2) +
my_ggplot_theme() +
theme(axis.text.x = element_text(
  angle = 90,
  size = 8,
  color = patients_factors_colors
)) +
geom_vline(
  xintercept = bars_separation,
  size = .3,
  color = 'darkred',
  linetype = 'dashed'
) +
labs(
  x = 'Patient',
  y = "REVOLVER evolutionary distance",
  title = 'REVOLVER cluster dendrogram',
  subtitle = x$annotation,
  caption = paste0('k = ', nclusters, ' clusters, ',
                    'n = ', x$n$patients, ' patients.')
) +
scale_color_manual(values = clusters_colors) 

ggsave(cluster_plot, file = file.path(result_dir, "plots/cluster_plot_full_cohort.jpg"), width = 8, height = 3, dpi = 200, units = "in")

####################################################
#### Change the branch colors
####################################################

# Use the hclust object from REVOLVER
hc <- FULL_clone_jackknife$cluster$fits$hc

# Get dendrogram data
dendro_data <- ggdendro::dendro_data(hc_, type = "rectangle")

# Map patient colors to branch segments
# Assume that `patients_factors_colors` is a named vector with patient names as names
# First, extract the order of the patients (i.e., the leaf order)
leaf_order <- hc_$order
leaf_labels <- hc_$labels[leaf_order]

# Create a named vector of colors in order of appearance
ordered_colors <- patients_factors_colors[leaf_labels]

# Assign colors to segments based on the xend of the segment (tip positions)
segment_data <- dendro_data$segments
label_data <- dendro_data$labels

# Match xend to patient positions and assign color
segment_data$branch_color <- NA
for (i in seq_len(nrow(segment_data))) {
  # If the segment ends at a leaf (i.e., connects to a label), color it
  xend <- segment_data$xend[i]
  match_idx <- which(label_data$x == xend)
  if (length(match_idx) == 1) {
    segment_data$branch_color[i] <- ordered_colors[match_idx]
  } else {
    segment_data$branch_color[i] <- NA
  }
}


# Step 1: Compute max x per color
color_bounds <- segment_data %>%
  filter(!is.na(branch_color)) %>%
  group_by(branch_color) %>%
  summarise(max_x = max(x, na.rm = TRUE), .groups = "drop") %>%
  arrange(max_x)

# Step 2: Build intervals (prepend 0 to the left)
color_ranges <- tibble(
  left = c(0, color_bounds$max_x[-nrow(color_bounds)]),
  right = color_bounds$max_x,
  assign_color = color_bounds$branch_color
)
# Step 3: Initialize new column
segment_data$branch_color_filled <- segment_data$branch_color

# Step 4: Assign colors to NA rows based on x falling in each interval
for (i in seq_len(nrow(color_ranges))) {
  range_left <- color_ranges$left[i]
  range_right <- color_ranges$right[i]
  color_to_assign <- color_ranges$assign_color[i]

  segment_data$branch_color_filled[
    is.na(segment_data$branch_color_filled) &
    segment_data$x > range_left &
    segment_data$x <= range_right
  ] <- color_to_assign
}

	# 1.	Define a function to get the color range for any given x value.
	# 2.	For each row in segment_data, check:
	# •	color_x  = get_color_for(x)
	# •	color_xend = get_color_for(xend)
	# •	If not equal, set branch_color_filled = "black"
# Helper function to assign color based on x position
get_color_for_x <- function(x_value) {
  for (i in seq_len(nrow(color_ranges))) {
    if (x_value > color_ranges$left[i] && x_value <= color_ranges$right[i]) {
      return(color_ranges$assign_color[i])
    }
  }
  return(NA)  # fallback
}

# Loop through each row and override color if x and xend fall in different intervals
for (i in seq_len(nrow(segment_data))) {
  if (is.na(segment_data$branch_color_filled[i])) next  # already NA, skip

  color_x <- get_color_for_x(segment_data$x[i])
  color_xend <- get_color_for_x(segment_data$xend[i])

  if (!is.na(color_x) && !is.na(color_xend) && color_x != color_xend) {
    segment_data$branch_color_filled[i] <- "black"
  }
}

segment_data_ <- segment_data %>% select(-branch_color) %>% rename(branch_color = branch_color_filled)

cluster_plot <- ggplot() +
  geom_segment(
    data = segment_data_,
    aes(x = x, y = y, xend = xend, yend = yend, color = branch_color),
    size = 0.8
  ) +
  scale_color_identity() +  # Use exact colors from vector
  geom_text(
    data = label_data,
    aes(x = x, y = y - 2, label = label, color = ordered_colors),
    angle = 90,
    hjust = 1,
    size = 2.5
  ) +
  my_ggplot_theme() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  ) +
  geom_vline(
    xintercept = bars_separation,
    size = 0.3,
    color = 'darkred',
    linetype = 'dashed'
  ) +
  labs(
    x = 'Patient',
    y = "REVOLVER evolutionary distance",
    title = 'REVOLVER cluster dendrogram',
    subtitle = x$annotation,
    caption = paste0('k = ', nclusters, ' clusters, ',
                     'n = ', x$n$patients, ' patients.')
  )

cluster_plot

ggsave(cluster_plot, file = file.path(result_dir, "plots/cluster_plot_full_cohort_colored.jpg"), width = 8, height = 3, dpi = 200, units = "in")

####################################################
#### Plot annotation bars
####################################################
# Get the orders of patients to plot labels
str(label_data)
patient_order <- label_data$label
cleaned_patient_order <- gsub("^(BE_|Non-BE_)", "", patient_order)

metadata_df <- read.csv(file.path(result_dir, "/20250126_Stages_Leanne.csv"))
metadata_df$OCCAMS_ID <- gsub("/", "_", metadata_df$OCCAMS_ID)

# Step 2: Filter metadata for patients in the list
matched_df <- metadata_df[metadata_df$OCCAMS_ID %in% cleaned_patient_order, ]

# Step 3: Set factor levels based on patient_order
matched_df$OCCAMS_ID <- factor(matched_df$OCCAMS_ID, levels = cleaned_patient_order)

# Step 4: Arrange the data frame to match the order in patient_order
matched_df <- matched_df[order(matched_df$OCCAMS_ID), ]

# Ensure OCCAMS_ID is a factor with levels in the original row order
matched_df$OCCAMS_ID <- factor(matched_df$OCCAMS_ID, levels = matched_df$OCCAMS_ID)

# Plot groups
label_bar <- ggplot(matched_df, aes(x = OCCAMS_ID, fill = group)) +
  geom_bar(stat = "count") +
  scale_fill_manual(values = c("pos" = "lightgreen", "neg" = "#329ca8")) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),              # remove grid lines
    axis.title = element_blank(),              # remove axis titles
    axis.text.y = element_text(size = 6),      # show y-axis labels (patient IDs)
    axis.text.x = element_text(size = 6, angle = 90, hjust = 0.5),
    axis.ticks = element_blank(),              # remove ticks
    legend.position = "right",
    plot.margin = margin(5, 5, 5, 5),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    plot.title = element_text(size = 10, hjust = 0.5)
  ) +
  labs(
    fill = "Group",
    title = "Patients grouped by status"
  )

ggsave(label_bar,file = file.path(result_dir, "plots/label_bar_phenotype.jpg"), width = 8, height = 1, dpi = 300, units = "in")


# Plot stages
label_bar <- ggplot(matched_df, aes(x = OCCAMS_ID, fill = ClinStageShort)) +
  geom_bar(stat = "count") +
  scale_fill_manual(values = c("I" = "#C1DBB3", "II" = "#FAEDCA", "III" = "#F2C078", "IV" = "#FE5D26")) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),              # remove grid lines
    axis.title = element_blank(),              # remove axis titles
    axis.text.y = element_text(size = 6),      # show y-axis labels (patient IDs)
    axis.text.x = element_text(size = 6, angle = 90, hjust = 0.5),
    axis.ticks = element_blank(),              # remove ticks
    legend.position = "right",
    plot.margin = margin(5, 5, 5, 5),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    plot.title = element_text(size = 10, hjust = 0.5)
  ) +
  labs(
    fill = "Group",
    title = "Patients grouped by status"
  )

ggsave(label_bar,file = file.path(result_dir, "plots/label_bar_stages.jpg"), width = 8, height = 1, dpi = 200, units = "in")

#### End of the script
#=====================================================================#