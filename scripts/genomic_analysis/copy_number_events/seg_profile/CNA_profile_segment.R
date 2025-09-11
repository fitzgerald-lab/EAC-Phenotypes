#=====================================================================#
#=====================================================================#
######                                                           
######  Script process EAC and BE CNA data to profile the segmented CNAs
######                                                           
#=====================================================================#
#=====================================================================#

# Author: Lianlian Wu
# Date: 2025-06-03

library(dplyr)
library(ggplot2)
library(scales)
library(stringr)
library(patchwork) 

# Define the source directory 
EAC_cna_path <- "/path/to/case_caveman_adjusted_combined.csv"
hg19_data_dir <- "/path/to/data"
EAC_phenotype_path <- "/path/to/case_phenotype.csv"
result_dir <- "/path/to/results/cna/seg_profile/"

#=================================#
######                                                           
######    Load data
######                                                           
#=================================#
EAC_cna_data <- read.csv(EAC_cna_path, header=TRUE)

head(EAC_cna_data)
print(nrow(EAC_cna_data))
print(unique(EAC_cna_data$sample_id))

sample_order_pos <- read.csv(file.path(sample_order_dir, "sample_id_order_fig_2_EAC_pos.csv"), stringsAsFactors = FALSE)$sample_id
sample_order_neg <- read.csv(file.path(sample_order_dir, "sample_id_order_fig_2_EAC_neg.csv"), stringsAsFactors = FALSE)$sample_id
sample_order_un <- read.csv(file.path(sample_order_dir, "sample_id_order_fig_2_EAC_un.csv"), stringsAsFactors = FALSE)$sample_id

EAC_cna_data_pos <- EAC_cna_data %>% filter(sample_id %in% sample_order_pos) 
print(length(sample_order_pos))
print(length(unique(EAC_cna_data_pos$sample_id)))
# unmatched_samples <- setdiff(sample_order_pos, unique(EAC_cna_data$sample_id))

EAC_cna_data_neg <- EAC_cna_data %>% filter(sample_id %in% sample_order_neg)
print(length(sample_order_neg))
print(length(unique(EAC_cna_data_neg$sample_id)))

EAC_cna_data_un <- EAC_cna_data %>% filter(sample_id %in% sample_order_un)
print(length(sample_order_un))
print(length(unique(EAC_cna_data_un$sample_id)))

#=====================================================================#
######                                                           
######  Prepare data for CNA profile segment plot
######                                                           
#=====================================================================#
## Chromosome info
gap_hg19_df <- read.table(file.path(hg19_data_dir, "gap.txt"),
  header = FALSE, sep = "\t", stringsAsFactors = FALSE)

gap_hg19_df_filtered <- gap_hg19_df %>%
  filter(V8 == "telomere") %>%
  select(V2, V3, V4, V8) %>%
  rename(Chromosome = V2, Centromere_Start = V3, Centromere_End = V4, type = V8) %>%
  mutate(
    Chromosome = gsub("chr", "", Chromosome),
    Centromere_Start = as.numeric(Centromere_Start),
    Centromere_End = as.numeric(Centromere_End)
  ) %>%
  filter(Chromosome %in% as.character(1:22)) %>% 
  filter(Centromere_Start != 0) %>%
  arrange(Chromosome)

## Ensure chromosome order is numeric (1–22), sort before cumsum
hg19_chrom_sizes <- gap_hg19_df_filtered %>%
  mutate(Chromosome = as.character(Chromosome)) %>%
  filter(Chromosome %in% as.character(1:22)) %>%  # only autosomes
  arrange(as.numeric(Chromosome)) %>%  # proper numeric order
  mutate(chr_length = Centromere_End) %>%
  mutate(cum_start = c(0, cumsum(chr_length)[-n()]))

## Define a function to process CNA data by adding chromosome position info
process_cna_data <- function(data_df, hg19_chrom_sizes) {
  # Join and calculate new positions
  processed_data_df <- data_df %>%
    left_join(hg19_chrom_sizes, by = "Chromosome") %>%
    mutate(
      Start_cum = Start + cum_start,
      End_cum = End + cum_start
    )
  
  # Convert to factor to keep chromosome order
  processed_data_df$Chromosome <- factor(processed_data_df$Chromosome, levels = as.character(1:22))
  
  # Prepare chromosome label positions
  chrom_labels <- processed_data_df %>%
    mutate(Chromosome = as.integer(as.character(Chromosome))) %>%  # ensure numeric sorting
    group_by(Chromosome) %>%
    summarise(chr_mid = (min(Start_cum) + max(End_cum)) / 2) %>%
    arrange(Chromosome)  # ensures correct 1-22 order
  
  chrom_labels$Chromosome <- factor(chrom_labels$Chromosome, levels = as.character(1:22))
  
  return(list(processed_data_df = processed_data_df, chrom_labels = chrom_labels))
}

# Process the CNA data for positive, negative, and unknown samples
EAC_cna_data_pos_chrom_info_added <- process_cna_data(EAC_cna_data_pos, hg19_chrom_sizes)[["processed_data_df"]]
chrom_labels_pos <- process_cna_data(EAC_cna_data_pos, hg19_chrom_sizes)[["chrom_labels"]]

EAC_cna_data_neg_chrom_info_added <- process_cna_data(EAC_cna_data_neg, hg19_chrom_sizes)[["processed_data_df"]]
chrom_labels_neg <- process_cna_data(EAC_cna_data_neg, hg19_chrom_sizes)[["chrom_labels"]]

EAC_cna_data_un_chrom_info_added <- process_cna_data(EAC_cna_data_un, hg19_chrom_sizes)[["processed_data_df"]]
chrom_labels_un <- process_cna_data(EAC_cna_data_un, hg19_chrom_sizes)[["chrom_labels"]]

#=====================================================================#
######                                                           
######    PLOT the CNA profile segment
######                                                           
#=====================================================================#
## If plot by sample_id order in Fig 2A============================================================================
EAC_cna_data_pos_chrom_info_added <- EAC_cna_data_pos_chrom_info_added %>%
  mutate(sample_id = factor(sample_id, levels = sample_order_pos))

max_cn <- max(EAC_cna_data_pos_chrom_info_added$ploidy_adjusted_CN, na.rm = TRUE)

plot_cna_profile_segment <- function(EAC_cna_data_chrom_info_added, chrom_labels) {
  ggplot(EAC_cna_data_chrom_info_added) +
    geom_rect(aes(xmin = as.numeric(factor(sample_id)) - 0.5,
                  xmax = as.numeric(factor(sample_id)) + 0.5,
                  ymin = Start_cum,
                  ymax = End_cum,
                  fill = ploidy_adjusted_CN)) +
    scale_fill_gradientn(colors = c("navy", "white", "firebrick", "darkred"),
                         values = rescale(c(0, 1, 2, 4)),
                         limits = c(0, 4),
                         name = "Ploidy-adjusted copy number") +
    scale_x_continuous(breaks = unique(as.numeric(factor(EAC_cna_data_chrom_info_added$sample_id))),
                       labels = unique(EAC_cna_data_chrom_info_added$sample_id),
                       expand = c(0, 0)) +
    scale_y_continuous(
      breaks = chrom_labels$chr_mid,
      labels = chrom_labels$Chromosome,
      expand = c(0, 0),
      trans = "reverse"
    ) +
    geom_hline(data = hg19_chrom_sizes, aes(yintercept = cum_start),
               linetype = "dashed", color = "grey60", linewidth = 0.3) +
    labs(y = "Chromosomes", x = "Sample ID") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 3),
      axis.text.y = element_text(size = 8),
      panel.grid = element_blank(),
      legend.position = "right"
    )
}

cna_profile_segment_pos <- plot_cna_profile_segment(EAC_cna_data_pos_chrom_info_added, chrom_labels_pos)
cna_profile_segment_neg <- plot_cna_profile_segment(EAC_cna_data_neg_chrom_info_added, chrom_labels_neg)
cna_profile_segment_un <- plot_cna_profile_segment(EAC_cna_data_un_chrom_info_added, chrom_labels_un)

## Save the plot
ggsave(file.path(result_dir, "cna_profile_segment_same_order_Fig2_pos.png"), cna_profile_segment_pos, width = 8, height = 8, dpi = 300)
ggsave(file.path(result_dir, "cna_profile_segment_same_order_Fig2_neg.png"), cna_profile_segment_neg, width = 8, height = 8, dpi = 300)
ggsave(file.path(result_dir, "cna_profile_segment_same_order_Fig2_un.png"), cna_profile_segment_un, width = 8, height = 8, dpi = 300)

## If plot by phenotypes========================================================================================================
EAC_cohort_phenotype <- readxl::read_excel(EAC_phenotype_path) 
print(colnames(EAC_cohort_phenotype))
EAC_cohort_phenotype <- EAC_cohort_phenotype %>% 
  dplyr::select(sample_id = DNA_ID, Phenotype, tStage) 

EAC_cohort_phenotype$sample_id <- gsub("SLX-18929\\.", "SLX-18929_", EAC_cohort_phenotype$sample_id)
EAC_cohort_phenotype$sample_id <- gsub("SLX-18928\\.", "SLX-18928_", EAC_cohort_phenotype$sample_id)

print(unique(EAC_cohort_phenotype$tStage))
EAC_cohort_phenotype_ordered <- EAC_cohort_phenotype %>%
  mutate(tStage = factor(tStage, levels = c("T1", "T2", "T3", "T4", "NA"))) %>%
  arrange(tStage)

## Define a function to order the CNA in satges info
order_cna_data_by_phenotype <- function(cna_data, phenotype_data) {
  ordered_data <- cna_data %>%
    left_join(phenotype_data, by = "sample_id") %>%
    mutate(sample_id = factor(sample_id, levels = phenotype_data$sample_id))
  print(length(unique(ordered_data$sample_id)))
  print(length(unique(cna_data$sample_id)))
  return(ordered_data)
}

EAC_cna_data_pos_ordered <- order_cna_data_by_phenotype(EAC_cna_data_pos_chrom_info_added, EAC_cohort_phenotype_ordered)
EAC_cna_data_neg_ordered <- order_cna_data_by_phenotype(EAC_cna_data_neg_chrom_info_added, EAC_cohort_phenotype_ordered)
EAC_cna_data_un_ordered <- order_cna_data_by_phenotype(EAC_cna_data_un_chrom_info_added, EAC_cohort_phenotype_ordered)

## Order the samples in most variable to less and within phenotype groups===========
# Step 1: Assign group labels (combine HGD and IMC)

# Step 2: Calculate absolute diff and weighted average per sample
# Function to compute instability scores and order sample IDs by tStage and score
get_ordered_sample_ids <- function(cna_data, tstage_levels = c("T1", "T2", "T3", "T4", "NA")) {
  instability_scores <- cna_data %>%
    mutate(Length = End - Start) %>%
    group_by(sample_id, tStage) %>%
    summarize(
      score = sum(abs(ploidy_adjusted_CN - 1) * Length, na.rm = TRUE) / sum(Length, na.rm = TRUE),
      .groups = "drop"
    )
  
  ordered_sample_ids <- instability_scores %>%
    arrange(factor(tStage, levels = tstage_levels), score) %>%
    pull(sample_id)
  
  return(ordered_sample_ids)
}

# Example usage:
ordered_sample_ids_pos <- get_ordered_sample_ids(EAC_cna_data_pos_ordered)
ordered_sample_ids_neg <- get_ordered_sample_ids(EAC_cna_data_neg_ordered)
ordered_sample_ids_un <- get_ordered_sample_ids(EAC_cna_data_un_ordered)

write.csv(ordered_sample_ids_pos, file = file.path(result_dir, "ordered_sample_ids_pos.csv"), row.names = FALSE)
write.csv(ordered_sample_ids_neg, file = file.path(result_dir, "ordered_sample_ids_neg.csv"), row.names = FALSE)
write.csv(ordered_sample_ids_un, file = file.path(result_dir, "ordered_sample_ids_un.csv"), row.names = FALSE)

# Step 4. Set factor levels
EAC_cna_data_pos_ordered <- EAC_cna_data_pos_ordered %>%
  mutate(sample_id = factor(sample_id, levels = ordered_sample_ids_pos))

EAC_cna_data_neg_ordered <- EAC_cna_data_neg_ordered %>%
  mutate(sample_id = factor(sample_id, levels = ordered_sample_ids_neg))

EAC_cna_data_un_ordered <- EAC_cna_data_un_ordered %>%
  mutate(sample_id = factor(sample_id, levels = ordered_sample_ids_un))

phenotype_colors <- c("T1" = "green", "T2" = "orange", "T3" = "blue","T4" = "red2", "NA" = "grey")

plot_cna_profile_segment_with_phenotype <- function(EAC_cna_data, chrom_labels, phenotype_colors = c("T1" = "green", "T2" = "orange", "T3" = "blue", "T4" = "red2", "NA" = "grey")) {
  # phenotype_bar
  phenotype_bar <- ggplot(
    distinct(EAC_cna_data, sample_id, tStage)
  ) +
    geom_tile(aes(
      x = as.numeric(sample_id),
      y = 1,
      fill = tStage
    )) +
    scale_fill_manual(values = phenotype_colors, name = "Phenotype") +
    scale_x_continuous(
      expand = c(0, 0),
      breaks = seq_along(levels(EAC_cna_data$sample_id)),
      labels = levels(EAC_cna_data$sample_id)
    ) +
    theme_void() +
    theme(
      legend.position = "top",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(0, 5, 0, 5)
    )
  
  # cna_profile_segment_ordered
  cna_profile_segment_ordered <- ggplot(EAC_cna_data) +
    geom_rect(aes(
      xmin = as.numeric(sample_id) - 0.5,
      xmax = as.numeric(sample_id) + 0.5,
      ymin = Start_cum,
      ymax = End_cum,
      fill = ploidy_adjusted_CN
    )) +
    scale_fill_gradientn(
      colors = c("navy", "white", "firebrick", "darkred"),
      values = scales::rescale(c(0, 1, 2, 4)),
      limits = c(0, 4),
      name = "Ploidy-adjusted copy number"
    ) +
    scale_x_continuous(
      breaks = seq_along(levels(EAC_cna_data$sample_id)),
      labels = levels(EAC_cna_data$sample_id),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = chrom_labels$chr_mid,
      labels = chrom_labels$Chromosome,
      expand = c(0, 0),
      trans = "reverse"
    ) +
    geom_hline(
      data = hg19_chrom_sizes,
      aes(yintercept = cum_start),
      linetype = "dashed", color = "grey60", linewidth = 0.3
    ) +
    labs(y = "Chromosomes", x = "Sample ID") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 3),
      axis.text.y = element_text(size = 8),
      panel.grid = element_blank(),
      legend.position = "right",
      plot.margin = margin(0, 5, 0, 5)
    )
  
  # Combine plots
  final_plot <- phenotype_bar / cna_profile_segment_ordered + patchwork::plot_layout(heights = c(0.02, 1))
  return(final_plot)
}

final_plot_pos <- plot_cna_profile_segment_with_phenotype (EAC_cna_data_pos_ordered, chrom_labels_pos)
final_plot_neg <- plot_cna_profile_segment_with_phenotype (EAC_cna_data_neg_ordered, chrom_labels_neg)
final_plot_un <- plot_cna_profile_segment_with_phenotype (EAC_cna_data_un_ordered, chrom_labels_un)

ggsave(
  file.path(result_dir, "cna_profile_segment_ordered_by_phenotype_pos.png"),
  final_plot_pos,
  width = 8,
  height = 8,
  dpi = 300)

ggsave(
  file.path(result_dir, "cna_profile_segment_ordered_by_phenotype_neg.png"),
  final_plot_neg,
  width = 8,
  height = 8,
  dpi = 300)

ggsave(
  file.path(result_dir, "cna_profile_segment_ordered_by_phenotype_un.png"),
  final_plot_un,
  width = 8,
  height = 8,
  dpi = 300)
