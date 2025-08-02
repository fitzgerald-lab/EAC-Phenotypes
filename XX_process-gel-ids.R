# NB: the GEL file names use the normal_abnormal convention rather than the abnormal_vs_normal used in our convention.

# read in the TMB data that i processed on the cluster
gel_bo_oac_snv_indel_purity <-
  vroom::vroom('data/original/wgs/gel_oac_tmb_sz_20230401.csv', col_names = F) %>%
  rename(
    DNA_ID = X1,
    snv_count = X2,
    indel_count = X3,
    estimated_purity = X4,
  ) %>%
  mutate(
    normal_DNA_ID = str_extract(DNA_ID, '.+_L'),
    normal_DNA_ID = str_replace(normal_DNA_ID, '_L', ''),
    abnormal_DNA_ID = str_extract(DNA_ID, '[:digit:]_.+'),
    abnormal_DNA_ID = str_replace(abnormal_DNA_ID, '[:digit:]_', ''),
  ) %>% 
  relocate(c(abnormal_DNA_ID, normal_DNA_ID), .after = DNA_ID)

gel_oac_cohort <-
  vroom::vroom('data/original/wgs/gel_submission_metadata_ginny_BO_OAC.csv.csv', col_names = TRUE) %>% 
  rename(abnormal_DNA_ID = abnormal_ID) %>% 
  filter(grepl('OCCAMS', ID))
  
gel_oac_total_snv_indel_TMB <- inner_join(gel_oac_cohort,gel_bo_oac_snv_indel_purity, "abnormal_DNA_ID") %>%
  select(ID, snv_count, indel_count, estimated_purity) %>%
  group_by(ID) %>%
  mutate(snv_count = mean(snv_count, na.rm = FALSE)) %>% 
  mutate(indel_count = mean(indel_count, na.rm = FALSE)) %>% 
  mutate(total_count = snv_count + indel_count) %>% 
  mutate(estimated_purity = mean(estimated_purity, na.rm = FALSE)) %>% 
  distinct(ID, .keep_all = TRUE)
  

