#=====================================================================#
#=====================================================================#
######                                                           
######    Script calculate arm level CNA data by aberrant proportions -- EAC
######                                                           
#=====================================================================#
#=====================================================================#

# Author: Lianlian Wu
# Date: 2025-05-22

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(patchwork)
library(purrr)
library(broom)

# Define the source directory 
data_dir <- "/path/to/data"
sample_id_order_dir <- "/path/to/ordered_sample_ids"
EAC_cna_path <- "/path/to/EAC_710_caveman_adjusted_combined.csv"
EAC_phenotype_path <- "/path/to/phenotype_data.xlsx"
result_dir <- "/path/to/EAC_results"

# Load data 
EAC_cna_data <- read.csv(EAC_cna_path)
head(EAC_cna_data)
print(nrow(EAC_cna_data))
EAC_cna_data_major_less_than_minor <- EAC_cna_data %>%
  filter(Major_CN < Minor_CN)
print(nrow(EAC_cna_data_major_less_than_minor))
#=====================================================================#
######                                                           
######  Prepare data for arm-level CNA analysis 
######                                                           
#=====================================================================#

# Calculate arm-level CNA ============================================
# Step 1: Define centromere positions for hg19
# download.file(
#   url = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/gap.txt.gz",
#   destfile = file.path(data_dir, "gap.txt.gz"),
#   method = "auto"  # Or "wget" if you're on Linux and have wget installed
# )
gap_hg19_df <- read.table(file.path(data_dir, "gap.txt"),
  header = FALSE, sep = "\t", stringsAsFactors = FALSE)

gap_hg19_df_filtered <- gap_hg19_df %>%
  filter(V8 == "centromere") %>%
  select(V2, V3, V4, V8) %>%
  rename(Chromosome = V2, Centromere_Start = V3, Centromere_End = V4, type = V8) %>%
  mutate(
    Chromosome = gsub("chr", "", Chromosome),
    Centromere_Start = as.numeric(Centromere_Start),
    Centromere_End = as.numeric(Centromere_End)
  ) %>%
  filter(Chromosome %in% as.character(1:22)) %>% 
  # mutate(Chromosome = as.numeric(Chromosome)) %>%
  arrange(Chromosome)

str(gap_hg19_df_filtered)

# Step 2: Assign p or q arm
# Join centromere info
cna_with_centromere <- EAC_cna_data %>%
  left_join(gap_hg19_df_filtered, by = "Chromosome")

# Split into three parts as needed
split_cna_by_arm <- function(cna_with_centromere) {
  bind_rows(
    # p-arm: Start to Centromere_Start
    cna_with_centromere %>%
      filter(Start < Centromere_Start) %>%
      mutate(
        arm = "p",
        End = pmin(End, Centromere_Start),
        Length = End - Start
      ) %>%
      filter(Length > 0),
    # centromere region
    cna_with_centromere %>%
      filter(Start < Centromere_End & End > Centromere_Start) %>%
      mutate(
        Start = pmax(Start, Centromere_Start),
        End = pmin(End, Centromere_End),
        arm = "centromere",
        Length = End - Start
      ) %>%
      filter(Length > 0),
    # q-arm: Centromere_End to End
    cna_with_centromere %>%
      filter(End > Centromere_End) %>%
      mutate(
        arm = "q",
        Start = pmax(Start, Centromere_End),
        Length = End - Start
      ) %>%
      filter(Length > 0)
  ) %>%
    select(sample_id, Chromosome, Start, End, Length, Total_CN, Major_CN, Minor_CN, ploidy, arm, Centromere_Start, Centromere_End)
}

data_df_with_arm <- split_cna_by_arm(cna_with_centromere)


# Step 3: Weighted average CN values per sample, chromosome, and arm
data_df_with_arm_p_q <- data_df_with_arm %>% filter(arm %in% c("p", "q"))
print(length(unique(data_df_with_arm_p_q$sample_id)))

CN_assigned <- data_df_with_arm_p_q %>%
  mutate(
    event = case_when(
      round(Total_CN) > round(ploidy) ~ "gain",
      round(Total_CN) < round(ploidy) ~ "loss",
      round(Total_CN) == round(ploidy) & Minor_CN == 0 ~ "neutral LOH",
      TRUE ~ "neutral"
    ))

CN_assigned_arm <- CN_assigned %>%
  group_by(sample_id, Chromosome, arm) %>%
  summarise(
    total_length = max(End, na.rm = TRUE) - min(Start, na.rm = TRUE),
    gain_length = sum(Length[event == "gain"], na.rm = TRUE),
    loss_length = sum(Length[event == "loss"], na.rm = TRUE),
    loh_length  = sum(Length[event == "neutral LOH"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    arm_event = case_when(
      gain_length / total_length > 0.5 ~ "gain",
      loss_length / total_length > 0.5 ~ "loss",
      loh_length  / total_length > 0.5 ~ "neutral LOH",
      TRUE ~ "neutral"
    )
  )

print(colnames(CN_assigned_arm))

arm_cout_df <- CN_assigned_arm %>%
  group_by(sample_id) %>%
  summarise(n_rows = n(), .groups = "drop")

#=====================================================================#
######                                                           
######  PLOT arm-level CNA heatmap and bar plot
######                                                           
#=====================================================================#
# Prepare data for heatmap and bar plot, and seperate the samples into different phenotypes==============================
## load phenotype data
EAC_cohort_phenotype <- readxl::read_excel(EAC_phenotype_path) 
print(colnames(EAC_cohort_phenotype))
EAC_cohort_phenotype <- EAC_cohort_phenotype %>% 
  dplyr::select(sample_id = DNA_ID, Phenotype, tStage) 

CN_assigned_arm_with_phenotype <- CN_assigned_arm %>%
  left_join(EAC_cohort_phenotype, by = "sample_id")

options(width = 200)
head(CN_assigned_arm)
head(CN_assigned_arm_with_phenotype)
print(length(unique(CN_assigned_arm_with_phenotype$sample_id)))
print(unique(CN_assigned_arm_with_phenotype$Phenotype))

CN_assigned_arm_pos <- CN_assigned_arm_with_phenotype %>% filter(Phenotype == "BE/IM EAC")
CN_assigned_arm_neg <- CN_assigned_arm_with_phenotype %>% filter(Phenotype == "Non-BE/IM EAC")
CN_assigned_arm_un <- CN_assigned_arm_with_phenotype %>% filter(Phenotype == "Unknown BE/IM EAC")

n_samples_pos <- length(unique(CN_assigned_arm_pos$sample_id))
n_samples_neg <- length(unique(CN_assigned_arm_neg$sample_id))
n_samples_un <- length(unique(CN_assigned_arm_un$sample_id))

# Define chr_arm levels and events level in desired order
chr_arms <- unlist(lapply(1:22, function(x) c(paste0(x, "p"), paste0(x, "q"))))

# Save all CNA events
CN_EAC_all_events <- CN_assigned_arm_with_phenotype %>%
    mutate(chr_arm = paste0(Chromosome, arm)) %>%
    mutate(
      chr_arm = factor(chr_arm, levels = chr_arms)) %>%
    select(-c(Chromosome, arm))
write.csv(CN_EAC_all_events, file.path(result_dir, "EAC_arm_level_CNA_all_events.csv"), row.names = FALSE)


prepare_arm_level_cna_data <- function(CN_assigned_arm, chr_arms, n_samples) {
  events_df <- CN_assigned_arm %>%
    mutate(chr_arm = paste0(Chromosome, arm)) %>%
    select(sample_id, chr_arm, arm_event) %>%
    mutate(
      chr_arm = factor(chr_arm, levels = chr_arms))
  
  summary_df <- events_df %>%
    filter(arm_event %in% c("gain", "loss", "neutral LOH")) %>%
    group_by(chr_arm, arm_event) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(chr_arm) %>%
    mutate(percentage = n / n_samples * 100) %>%
    mutate(
      chr_arm = factor(chr_arm, levels = rev(chr_arms)),
      arm_event = factor(arm_event, levels = c("neutral LOH", "loss", "gain")))
  
  list(events_df = events_df, summary_df = summary_df)
}

# Apply the function to each group
arm_level_cna_results_all <- prepare_arm_level_cna_data(CN_assigned_arm_with_phenotype, chr_arms, 710)

## BE positive
arm_level_cna_results_pos <- prepare_arm_level_cna_data(CN_assigned_arm_pos, chr_arms, n_samples_pos)
events_df_pos <- arm_level_cna_results_pos$events_df
summary_df_pos <- arm_level_cna_results_pos$summary_df

## BE negative
arm_level_cna_results_neg <- prepare_arm_level_cna_data(CN_assigned_arm_neg, chr_arms, n_samples_neg)
events_df_neg <- arm_level_cna_results_neg$events_df
summary_df_neg <- arm_level_cna_results_neg$summary_df

## BE unknown
arm_level_cna_results_un <- prepare_arm_level_cna_data(CN_assigned_arm_un, chr_arms, n_samples_un)
events_df_un <- arm_level_cna_results_un$events_df
summary_df_un <- arm_level_cna_results_un$summary_df

# Set the order of sample_id the same as Fig 2A==============================
sample_order_pos <- read.csv(file.path(sample_id_order_dir, "sample_id_order_fig_2_EAC_pos.csv"), stringsAsFactors = FALSE)$sample_id
sample_order_neg <- read.csv(file.path(sample_id_order_dir, "sample_id_order_fig_2_EAC_neg.csv"), stringsAsFactors = FALSE)$sample_id
sample_order_un <- read.csv(file.path(sample_id_order_dir, "sample_id_order_fig_2_EAC_un.csv"), stringsAsFactors = FALSE)$sample_id


# Create the heatmap 
# Function to create arm-level CNA heatmap and bar plot-------------------------------------
plot_arm_level_cna <- function(events_df, summary_df, sample_order, chr_arms) {
  # Prepare events_df with correct sample order and chr_arm factor levels
  events_df_ordered <- events_df %>%
    mutate(
      sample_id = factor(sample_id, levels = sample_order),
      chr_arm = factor(chr_arm, levels = rev(chr_arms))
    )
  
  # Prepare summary_df with correct chr_arm factor levels
  summary_df <- summary_df %>%
    mutate(
      chr_arm = factor(chr_arm, levels = rev(chr_arms)),
      arm_event = factor(arm_event, levels = c("neutral LOH", "loss", "gain"))
    )
  
  # Heatmap plot
  heatmap_plot <- ggplot(events_df_ordered, aes(x = sample_id, y = chr_arm, fill = arm_event)) +
    geom_tile(color = "white") +
    scale_fill_manual(
      values = c(gain = "#EA2E49", loss = "#1d7cdb", 'neutral LOH' = "orange", neutral = "grey90"),
      name = "Event"
    ) +
    labs(x = "", y = "") +
    theme_minimal() +
    theme(
      axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=0.5, size = 4.5),
      axis.text.x = element_text(angle = 90, vjust = 0, hjust=0.5, size = 2),
      legend.position = "none"
    )
  
  # Bar plot
  bar_plot <- ggplot(summary_df, aes(x = chr_arm, y = percentage, fill = arm_event)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = c(gain = "#EA2E49", loss = "#1d7cdb", 'neutral LOH' = "orange")) +
    theme_minimal() +
    coord_flip() +
    scale_y_continuous(limits = c(0, 80), breaks = seq(0, 80, by = 20)) +
    theme(
      axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5, size = 3),
      axis.ticks.x = element_line(color = "black", linewidth = 0.1),
      axis.line.x = element_line(color = "black", linewidth = 0.1),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      panel.grid = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "", y = "")
  
  # Combine plots
  final_plot <- (heatmap_plot | bar_plot) + plot_layout(widths = c(1, 0.25))
  
  return(final_plot)
}

# Example usage for BE/IM EAC positive group:
final_plot_pos <- plot_arm_level_cna(events_df_pos, summary_df_pos, sample_order_pos, chr_arms)
final_plot_neg <- plot_arm_level_cna(events_df_neg, summary_df_neg, sample_order_neg, chr_arms)
final_plot_un <- plot_arm_level_cna(events_df_un, summary_df_un, sample_order_un, chr_arms)

ggsave(file.path(result_dir, "arm_level_CNA_bar_plot_ordered_Fig2_pos.png"),
  final_plot_pos, width = 4.2, height = 3.5, dpi = 300)

ggsave(file.path(result_dir, "arm_level_CNA_bar_plot_ordered_Fig2_neg.png"),
  final_plot_neg, width = 4.2, height = 3.5, dpi = 300)

ggsave(file.path(result_dir, "arm_level_CNA_bar_plot_ordered_Fig2_un.png"),
  final_plot_un, width = 4.2, height = 3.5, dpi = 300)
#=====================================================================#
######                                                           
######  PLOT arm-level CNA by phenotypes
######                                                           
#=====================================================================#
## Read the sample id orders from seg_profile data
ordered_sample_ids_seg_pos <- read.csv(file.path(sample_id_order_dir, "ordered_sample_ids_pos.csv"), stringsAsFactors = FALSE)$x
ordered_sample_ids_seg_neg <- read.csv(file.path(sample_id_order_dir, "ordered_sample_ids_neg.csv"), stringsAsFactors = FALSE)$x
ordered_sample_ids_seg_un <- read.csv(file.path(sample_id_order_dir, "ordered_sample_ids_un.csv"), stringsAsFactors = FALSE)$x

# final_plot_pos_ordered_by_stages <- plot_arm_level_cna(events_df_pos, summary_df_pos, ordered_sample_ids_seg_pos, chr_arms)


plot_arm_level_cna_with_phenotype <- function(events_df, ordered_sample_ids, phenotype_df) {
  ## Define the colors for TStages
  phenotype_colors <- c("T1" = "green", "T2" = "orange", "T3" = "blue","T4" = "red2", "NA" = "grey")
  # Merge phenotype info and set factor levels
  events_df_ordered <- events_df %>%
    left_join(phenotype_df, by = "sample_id") %>%
    mutate(
      sample_id = factor(sample_id, levels = ordered_sample_ids),
      chr_arm = factor(chr_arm, levels = rev(chr_arms))
    )
  print(colnames(events_df_ordered))
  
  # Phenotype bar plot
  phenotype_bar <- ggplot(
    distinct(events_df_ordered, sample_id, tStage)
  ) +
    geom_tile(aes(
      x = as.numeric(sample_id),
      y = 1,
      fill = tStage
    )) +
    scale_fill_manual(values = phenotype_colors, name = "Phenotype") +
    scale_x_continuous(
      expand = c(0, 0),
      breaks = seq_along(levels(events_df_ordered$sample_id)),
      labels = levels(events_df_ordered$sample_id)
    ) +
    theme_void() +
    theme(
      legend.position = "top",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(0, 5, 0, 5)
    )
  
  # CNA profile heatmap
  cna_profile_arm_ordered <- ggplot(events_df_ordered, aes(x = sample_id, y = chr_arm, fill = arm_event)) +
    geom_tile(color = "white") +
    scale_fill_manual(
      values = c(gain = "#EA2E49", loss = "#1d7cdb", 'neutral LOH' = "orange", neutral = "grey90"),
      name = "Event"
    ) +
    labs(x = "", y = "") +
    theme_minimal() +
    theme(
      axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=0.5, size = 6),
      axis.text.x = element_text(angle = 90, vjust = 0, hjust=0.5, size = 2),
      legend.position = "none"
    )
  
  # Combine
  final_plot <- phenotype_bar / cna_profile_arm_ordered + patchwork::plot_layout(heights = c(0.02, 1))
  return(final_plot)
}

final_arm_plot_pos <- plot_arm_level_cna_with_phenotype(events_df_pos, ordered_sample_ids_seg_pos, EAC_cohort_phenotype)
final_arm_plot_neg <- plot_arm_level_cna_with_phenotype(events_df_neg, ordered_sample_ids_seg_neg, EAC_cohort_phenotype)
final_arm_plot_un <- plot_arm_level_cna_with_phenotype(events_df_un, ordered_sample_ids_seg_un, EAC_cohort_phenotype)

ggsave(file.path(result_dir, "arm_level_CNA_bar_plot_ordered_by_phenotype_pos.png"),
  final_plot, width = 4, height = 6, dpi = 300)

ggsave(file.path(result_dir, "arm_level_CNA_bar_plot_ordered_by_phenotype_neg.png"),
  final_arm_plot_neg, width = 4, height = 6, dpi = 300)

ggsave(file.path(result_dir, "arm_level_CNA_bar_plot_ordered_by_phenotype_un.png"),
  final_arm_plot_un, width = 4, height = 6, dpi = 300)


#=====================================================================#
######                                                           
######  calculate the diffrence of arm-level CNA between BE/IM EAC and Non-BE/IM EAC
######                                                           
#=====================================================================#
# Function to calculate arm-level CNA proportions for a given events_df and sample count
calculate_arm_level_cna_proportions <- function(events_df, n_samples, phenotype_label) {
  proportions <- events_df %>%
    group_by(chr_arm, arm_event) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(
      percentage = n / n_samples * 100,
      Phenotype = phenotype_label
    )
  
  proportions_aberrant <- events_df %>%
    mutate(aberrant = ifelse(arm_event %in% c("gain", "loss", "neutral LOH"), 1, 0)) %>%
    group_by(chr_arm, aberrant) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(
      percentage = n / n_samples * 100,
      Phenotype = phenotype_label
    )
  
  list(proportions = proportions, proportions_aberrant = proportions_aberrant)
}

# Apply the function to BE/IM EAC and Non-BE/IM EAC
pos_results <- calculate_arm_level_cna_proportions(events_df_pos, n_samples_pos, "BE/IM EAC")
proportions_pos <- pos_results$proportions
proportions_aberrant_pos <- pos_results$proportions_aberrant

neg_results <- calculate_arm_level_cna_proportions(events_df_neg, n_samples_neg, "Non-BE/IM EAC")
proportions_neg <- neg_results$proportions
proportions_aberrant_neg <- neg_results$proportions_aberrant

# Perform chi-squared test for each chr_arm
# Function to perform chi-squared test for each chr_arm between two phenotypes
perform_arm_level_cna_chisq_test <- function(proportions_aberrant_pos, proportions_aberrant_neg) {
  combined_df <- bind_rows(proportions_aberrant_pos, proportions_aberrant_neg)
  wide_df <- combined_df %>%
    pivot_wider(
      id_cols = chr_arm,
      names_from = c(Phenotype, aberrant),
      values_from = n,
      values_fill = 0
    )
  
  test_results <- wide_df %>%
    mutate(
      chisq_result = pmap(
        list(`BE/IM EAC_0`, `BE/IM EAC_1`, `Non-BE/IM EAC_0`, `Non-BE/IM EAC_1`),
        ~ {
          mat <- matrix(c(..1, ..2, ..3, ..4), nrow = 2, byrow = TRUE)
          chisq.test(mat)
        }
      ),
      p_value = map_dbl(chisq_result, "p.value")
    ) %>%
    select(chr_arm, p_value)
  
  # Apply multiple testing correction
  test_results <- test_results %>%
    mutate(
      p_adj_fdr = p.adjust(p_value, method = "fdr"),
      p_adj_bonferroni = p.adjust(p_value, method = "bonferroni")
    )
  
  return(test_results)
}

test_results <- perform_arm_level_cna_chisq_test(proportions_aberrant_pos, proportions_aberrant_neg)
write.csv(test_results, file.path(result_dir, "arm_level_CNA_chisq_test_results.csv"), row.names = FALSE)

# Only count gain and loss------------------------------
calculate_arm_level_cna_proportions_gain_loss_only <- function(events_df, n_samples, phenotype_label) {
  proportions_aberrant <- events_df %>%
    mutate(aberrant = ifelse(arm_event %in% c("gain", "loss"), 1, 0)) %>%
    group_by(chr_arm, aberrant) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(
      percentage = n / n_samples * 100,
      Phenotype = phenotype_label
    )
  
  return(proportions_aberrant)
}

# Apply the function to BE/IM EAC and Non-BE/IM EAC for gain and loss only
pos_results_gain_gain <- calculate_arm_level_cna_proportions_gain_loss_only(events_df_pos, n_samples_pos, "BE/IM EAC")
neg_results_gain_loss <- calculate_arm_level_cna_proportions_gain_loss_only(events_df_neg, n_samples_neg, "Non-BE/IM EAC")

test_results_gain_loss <- perform_arm_level_cna_chisq_test(pos_results_gain_gain, neg_results_gain_loss)
write.csv(test_results, file.path(result_dir, "arm_level_CNA_chisq_test_results.csv"), row.names = FALSE)

# For comparing mutiple types of events ------------------------------
combined_df <- bind_rows(proportions_pos, proportions_neg)
wide_df <- combined_df %>%
  select(chr_arm, Phenotype, arm_event, n) %>%
  pivot_wider(
    names_from = Phenotype,
    values_from = n,
    values_fill = 0
  )

# Run chi-square test for each arm
chi_results <- wide_df %>%
  group_by(chr_arm) %>%
  summarise(
    test_result = list({
      tbl <- matrix(c(
        `BE/IM EAC`[arm_event == "gain"],
        `BE/IM EAC`[arm_event == "loss"],
        `BE/IM EAC`[arm_event == "neutral"],
        `BE/IM EAC`[arm_event == "neutral LOH"],
        `Non-BE/IM EAC`[arm_event == "gain"],
        `Non-BE/IM EAC`[arm_event == "loss"],
        `Non-BE/IM EAC`[arm_event == "neutral"],
        `Non-BE/IM EAC`[arm_event == "neutral LOH"]
      ), nrow = 2, byrow = TRUE)

      # Decide test type
      test_obj <- chisq.test(tbl)
      if (any(test_obj$expected < 5)) {
        fisher.test(tbl)
      } else {
        test_obj
      }
    }),
    method_used = list({
      tbl <- matrix(c(
        `BE/IM EAC`[arm_event == "gain"],
        `BE/IM EAC`[arm_event == "loss"],
        `BE/IM EAC`[arm_event == "neutral"],
        `BE/IM EAC`[arm_event == "neutral LOH"],
        `Non-BE/IM EAC`[arm_event == "gain"],
        `Non-BE/IM EAC`[arm_event == "loss"],
        `Non-BE/IM EAC`[arm_event == "neutral"],
        `Non-BE/IM EAC`[arm_event == "neutral LOH"]
      ), nrow = 2, byrow = TRUE)

      test_obj <- chisq.test(tbl)
      if (any(test_obj$expected < 5)) {
        "Fisher"
      } else {
        "Chi-squared"
      }
    }),
    .groups = "drop"
  ) %>%
  mutate(
    p_value = map_dbl(test_result, ~ .x$p.value),
    method = unlist(method_used),
    p_adj_fdr = p.adjust(p_value, method = "fdr"),
    p_adj_bonferroni = p.adjust(p_value, method = "bonferroni")
  )

# Clean output
chi_results_cleaned <- chi_results %>%
  select(chr_arm, p_value, p_adj_fdr, p_adj_bonferroni)

# View results
print(chi_results_cleaned)
write.csv(chi_results_cleaned, file.path(result_dir, "arm_level_CNA_chisq_test_results_multiple_events.csv"), row.names = FALSE)

#================================================================================================================================#


#=====================================================================#
######                                                           
######  Unsuprised clustering of arm-level CNA
######                                                           
#=====================================================================#
library(ComplexHeatmap)
library(circlize)

CN_assigned_arm_pos_clean <- CN_assigned_arm_pos %>% 
  mutate(chr_arm = paste0(Chromosome, arm)) %>%
  select(sample_id, chr_arm, arm_event, Phenotype, tStage)

CN_assigned_arm_neg_clean <- CN_assigned_arm_neg %>% 
  mutate(chr_arm = paste0(Chromosome, arm)) %>%
  select(sample_id, chr_arm, arm_event, Phenotype, tStage)

combine_pos_neg_events <- rbind(CN_assigned_arm_pos_clean, CN_assigned_arm_neg_clean)
write.csv(combine_pos_neg_events, file.path(result_dir, "arm_level_CNA_pos_neg_events.csv"), row.names = FALSE)

# End of the script
#================================================================================================================================#

