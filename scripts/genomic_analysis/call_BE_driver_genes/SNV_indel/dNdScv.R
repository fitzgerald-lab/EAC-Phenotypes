library(dplyr)
library(tidyr)
library(stringr)
library(dndscv)

data_path <- "./data/genomic/call_BE_driver_genes/input_for_dnds_example.csv"
output_dir <- "./results/"
dNdS_data <- read.csv(data_path)
options(width = 200)
head(dNdS_data, n=200)
#sample_id <- unique(dNdS_data$sampleID)
#write.csv(sample_id, paste0(output_dir, 'samples_with_snv_indel.csv'))

#dndsout = dndscv(dNdS_data)
dNdS_data = dndscv(dNdS_data)

#sel_cv = dndsout$sel_cv
dNdS_data = dNdS_data$sel_cv

#print(head(sel_cv), digits = 3)
#print(dndsout$globaldnds)
#signif_genes = sel_cv[sel_cv$qglobal_cv<0.1, c("gene_name","qglobal_cv")]
#print(signif_genes)

dNdS_data = dNdS_data[dNdS_data$qglobal_cv<0.1, c("gene_name","qglobal_cv")]
print(dNdS_data)

save(dNdS_data, file = paste0(output_dir, "dndsout.RData"))
write.csv(dNdS_data, paste0(output_dir, "dndsout_significant_genes_from_BE_samples.csv"))
