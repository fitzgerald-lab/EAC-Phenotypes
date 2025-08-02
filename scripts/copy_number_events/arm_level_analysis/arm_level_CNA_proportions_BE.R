#=====================================================================#
#=====================================================================#
######                                                           
######    Script calculate arm level CNA data by aberrant proportions -- BE
######                                                           
#=====================================================================#
#=====================================================================#

# Author: Lianlian Wu
# Date: 2025-05-22

library(dplyr)
library(ggplot2)
library(stringr)
library(patchwork) 

# Define the source directory 
data_dir <- "/path/to/data"
combined_cna_path <- file.path(data_dir, "BE_samples/alelle_cna_data.csv")
sample_order_df <- read.csv("/path/to/sample_id_order_fig_2.csv", stringsAsFactors = FALSE)
result_dir <- "/path/to/output_directory"

# Load data 
combined_cna_data <- read.csv(combined_cna_path)
head(combined_cna_data)
print(nrow(combined_cna_data))
combined_cna_data_major_less_than_minor <- combined_cna_data %>%
  filter(Major_CN < Minor_CN)

print(length(unique(combined_cna_data_major_less_than_minor$sample_id)))
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
cna_with_centromere <- combined_cna_data %>%
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
# Prepare data for heatmap and bar plot==============================
events_df <- CN_assigned_arm %>%
  mutate(chr_arm = paste0(Chromosome, arm)) %>%
  select(sample_id, chr_arm, arm_event)

summary_df <- events_df %>%
  filter(arm_event %in% c("gain", "loss", "neutral LOH")) %>%
  group_by(chr_arm, arm_event) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(chr_arm) %>%
  mutate(percentage = n / 284 * 100)

# Define chr_arm levels and events level in desired order
chr_arms <- unlist(lapply(1:22, function(x) c(paste0(x, "p"), paste0(x, "q"))))

events_df <- events_df %>%
  mutate(
    chr_arm = factor(chr_arm, levels = chr_arms))

summary_df <- summary_df %>%
  mutate(
    chr_arm = factor(chr_arm, levels = rev(chr_arms)),
    arm_event = factor(arm_event, levels = c("neutral LOH", "loss", "gain")))

write.csv(events_df, file.path(result_dir, "BE/BE_arm_level_CNA_events.csv"), row.names = FALSE)
write.csv(summary_df, file.path(result_dir, "BE/BE_arm_level_CNA_summary.csv"), row.names = FALSE)
# Set the order of sample_id the same as Fig 2A==============================
sample_order <- sample_order_df$sample_id
events_df_ordered <- events_df %>%
  mutate(matched_sample = sapply(sample_id, function(x) {
    match <- sample_order[str_detect(x, sample_order)]
    if (length(match) > 0) match[1] else NA
  }))
print(unique(events_df_ordered$matched_sample))

events_df_ordered <- events_df_ordered %>%
  mutate(matched_sample = factor(matched_sample, levels = sample_order)) %>%
  mutate(chr_arm = factor(chr_arm, levels = rev(chr_arms)))


# Create the heatmap 
heatmap_plot <- ggplot(events_df_ordered, aes(x = matched_sample, y = chr_arm, fill = arm_event)) +
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

ggsave(
  heatmap_plot,
  filename = file.path(result_dir, paste0("arm_level_CNA_heatmap_newly_liftover_allele_corrected", Sys.Date(), ".pdf")),
  width = 10, height = 8, dpi = 300
)

# Create the bar plot
bar_plot <- ggplot(summary_df, aes(x = chr_arm, y = percentage, fill = arm_event)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c(gain = "#EA2E49", loss = "#1d7cdb", 'neutral LOH' = "orange")) +
  theme_minimal() +
  coord_flip() +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5, size = 3),
    axis.ticks.x = element_line(color = "black", linewidth = 0.1),
    axis.line.x = element_line(color = "black", linewidth = 0.1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid = element_blank(),           # Remove all grid lines
    panel.grid.major = element_blank(),     # Remove major grid lines
    panel.grid.minor = element_blank()      # Remove minor grid lines
  ) +
  labs(x = "", y = "")

# ggsave(bar_plot,
#   filename = file.path(result_dir, paste0("arm_level_CNA_bar_plot_with_LOH_newly_liftover_allele_corrected", Sys.Date(), ".pdf")),
#   width = 12, height = 6)

final_plot <- (heatmap_plot | bar_plot) + plot_layout(widths = c(1, 0.25))

ggsave(
  file.path(result_dir, "arm_level_CNA_bar_plot_ordered_Fig2.png"),
  final_plot,
  width = 3.6,
  height = 3.8,
  dpi = 300
)

#=====================================================================#
######                                                           
######  PLOT arm-level CNA by phenotypes
######                                                           
#=====================================================================#
BO_cohort_phenotype <- read.csv("/path/to/case_phenotype.csv") 
BO_cohort_phenotype <- BO_cohort_phenotype %>% 
  dplyr::select(sequence_id, Phenotype = specimen_pathology) %>%
  filter(!is.na(sequence_id))
BO_cohort_phenotype_ordered <- BO_cohort_phenotype %>%
  mutate(Phenotype = factor(Phenotype, levels = c("NDBE", "LGD", "HGD", "IMC"))) %>%
  arrange(Phenotype)

events_df_ordered_phenotype <- events_df_ordered %>%
  left_join(BO_cohort_phenotype_ordered, by = c("matched_sample" = "sequence_id")) %>%
  mutate(matched_sample = factor(matched_sample, levels = BO_cohort_phenotype_ordered$sequence_id))

## Order the samples in most variable to less and within phenotype groups===========
# Step 1: Assign group labels (combine HGD and IMC)
events_df_ordered_phenotype <- events_df_ordered_phenotype %>%
  mutate(
    phenotype_group = case_when(
      Phenotype == "NDBE" ~ "NDBE",
      Phenotype == "LGD" ~ "LGD",
      Phenotype %in% c("HGD", "IMC") ~ "HGD/IMC",
      TRUE ~ NA_character_
    )
  )

# Step 2/3. Define final sample order
ordered_sample_ids <- read.csv(file.path("/path/to/ordered_sample_ids.csv"))$x

# Step 4. Set factor levels
events_df_ordered_phenotype <- events_df_ordered_phenotype %>%
  mutate(matched_sample = factor(matched_sample, levels = ordered_sample_ids))


phenotype_colors <- c("NDBE" = "green", "LGD" = "orange", "HGD/IMC" = "red")

phenotype_bar <- ggplot(
  distinct(events_df_ordered_phenotype, matched_sample, phenotype_group)
) +
  geom_tile(aes(
    x = as.numeric(matched_sample),
    y = 1,
    fill = phenotype_group
  )) +
  scale_fill_manual(values = phenotype_colors, name = "Phenotype") +
  scale_x_continuous(
    expand = c(0, 0),
    breaks = seq_along(levels(cna_data_cum_ordered$matched_sample)),
    labels = levels(cna_data_cum_ordered$matched_sample)
  ) +
  theme_void() +
  theme(
    legend.position = "top",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(0, 5, 0, 5)
  )

cna_profile_arm_ordered <- ggplot(events_df_ordered_phenotype, aes(x = matched_sample, y = chr_arm, fill = arm_event)) +
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

final_plot <- phenotype_bar / cna_profile_arm_ordered + plot_layout(heights = c(0.02, 1))
ggsave(
  file.path(result_dir, "arm_level_CNA_bar_plot_ordered_by_phenotype.png"),
  final_plot,
  width = 4,
  height = 6,
  dpi = 300
)

# End of script
#=====================================================================#