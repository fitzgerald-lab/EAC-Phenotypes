#===========================================================================================
#===========================================================================================      
######   Script combine SNV and indel consequences for all BE samples and calculate their proportions   
######   Input = ../gel/consequence_gel_samples_filtered_unique.csv, ../samples_211/"consequence_indel_211_samples_filtered_unique.csv, ../samples_211/consequence_snv_211_samples_filtered_unique.csv
######   Output = consequence_gel_211_samples_filtered_unique.csv, proportions_gel_and_211.csv, proportions_top_gel_and_211.jpeg                                                  
#===========================================================================================
#===========================================================================================

# Author: Lianlian Wu
# Date: 2024-05-04

library(ggplot2)
library(dplyr)
library(tidyr)
#library(openxlsx)

file_path_211 <- "/scratchb/stlab-icgc/users/wu04/project/bo_gene_list/results/"
file_path_gel <- "/scratchb/stlab-icgc/users/wu04/project/bo_gene_list/results/gel/consequence_gel_samples_filtered_unique.csv"
output_dir <- "/scratchb/stlab-icgc/users/wu04/project/bo_gene_list/results/gel/"

indel_data <- read.csv(paste0(file_path_211,"samples_211/", "consequence_indel_211_samples_filtered_unique.csv"))
snv_df <- read.csv(paste0(file_path_211,"samples_211/", "consequence_snv_211_samples_filtered_unique.csv"))
gel_df <- read.csv(file_path_gel)

all_variations <- bind_rows(indel_data, snv_df, gel_df)
write.csv(all_variations, paste0(output_dir, "consequence_gel_211_samples_filtered_unique.csv"))

print("The number of samples:")
sample_number <- length(unique(all_variations$sequence_id))
print(sample_number)
print("The number of genes:")
print(length(unique(all_variations$symbol)))

gene_variation_proportions <- all_variations %>% group_by(sequence_id, symbol) %>% 
  summarise(variation_combo = paste(sort(unique(variation)), collapse = "_")) %>%
  group_by(symbol, variation_combo) %>% summarise(count = n()) %>% ungroup() %>%
  pivot_wider(names_from = variation_combo, values_from = count, values_fill = list(count = 0)) %>%
  mutate(total_samples = length(unique(all_variations$sequence_id)),
         prop_snv_only = `snv` / total_samples,
         prop_indel_only = `indel` / total_samples,
         prop_both = `indel_snv`/ total_samples) %>%
  select(symbol, prop_snv_only, prop_indel_only, prop_both) %>%
  mutate(sum_props = prop_snv_only + prop_indel_only + prop_both)

gene_variation_proportions <- gene_variation_proportions %>%
  arrange(desc(sum_props))

proportions_df_long <- gene_variation_proportions %>%
  pivot_longer(cols = starts_with("prop"), names_to = "prop_type", values_to = "value") 

proportions_df_long <- proportions_df_long %>%
  filter(symbol != "TTN")

write.csv(proportions_df_long, paste0(output_dir, "proportions_gel_and_211.csv"))

proportions_top <- proportions_df_long %>% arrange(desc(sum_props)) %>% 
  slice_head(n = 150)
print("top_genes:")
print(proportions_top$symbol)

top_genes_plot <- ggplot(proportions_top, aes(x = reorder(symbol, sum_props), y = value, fill = prop_type)) +
  geom_col() +
  coord_flip() + # Flip coordinates for horizontal layout
  labs(title = "Proportions of Mutation Types per Gene",
       x = "Gene Symbol",
       y = "Proportion") +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal()

output_file_path <- paste0(output_dir, "proportions_top_gel_and_211.jpeg")
ggsave(output_file_path, plot = top_genes_plot, width = 8, height = 16, dpi = 400)

