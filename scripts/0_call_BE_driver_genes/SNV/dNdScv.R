library(dplyr)
library(tidyr)
library(stringr)
library(dndscv)

data_path <- "C:/Users/lw748/OneDrive - University of Cambridge/projects/BO_gene_list/data/290_samples_combined_for_dNdS.csv"
output_dir <- "C:/Users/lw748/OneDrive - University of Cambridge/projects/BO_gene_list/results/"
dNdS_data <- read.csv(data_path)
#sample_id <- unique(dNdS_data$sampleID)
#write.csv(sample_id, paste0(output_dir, 'samples_with_snv_indel.csv'))

# Exclude GM
dNdS_data_exclude_GM <- dNdS_data %>%
  filter(!str_detect(sampleID, "LP6008269-DNA_H01|SLX-21762.idtUDP0094|SLX-21762.idtUDP0096|SLX-22426.idtUDP0109|SLX-22426.idtUDP0119|SLX-22427.idtUDP0130"))
print(length(unique(dNdS_data_exclude_GM$sampleID)))

#dndsout = dndscv(dNdS_data)
dndsout_exclude_GM = dndscv(dNdS_data_exclude_GM)

#sel_cv = dndsout$sel_cv
sel_cv_exclude_GM = dndsout_exclude_GM$sel_cv

#print(head(sel_cv), digits = 3)
#print(dndsout$globaldnds)
#signif_genes = sel_cv[sel_cv$qglobal_cv<0.1, c("gene_name","qglobal_cv")]
#print(signif_genes)

signif_genes_exclude_GM = sel_cv_exclude_GM[sel_cv_exclude_GM$qglobal_cv<0.1, c("gene_name","qglobal_cv")]
print(signif_genes_exclude_GM)

save(dndsout_exclude_GM, file = paste0(output_dir, "dndsout_284_samples_excluding_GM.RData"))
write.csv(signif_genes_exclude_GM, paste0(output_dir, "dndsout_significant_genes_from_284_samples_excluding_GM.csv"))
