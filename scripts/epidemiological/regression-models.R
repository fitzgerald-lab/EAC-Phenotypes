library(tidyverse)
library(finalfit)


dependent <- "phenotype"
independent_base <- c('dem_age_group', 'dem_gender')
independent_full <-
  c(
    "dem_age_group",
    "dem_gender",
    'dem_bmi_group',
    'smk_ever_smoked',
    'med_taken_nsaids',
    'sym_heartburn_inferred',
    'pth_pTNM'
  )
independent_no_nsaid <-
  c(
    "dem_age_group",
    "dem_dem_gender",
    'dem_bmi_group',
    'smk_ever_smoked',
    # 'med_taken_nsaids',
    'sym_heartburn_inferred',
    'pth_pTNM'
  )

independent_permute <-
  c("dem_bmi_group",
    "smk_ever_smoked",
    "med_med_taken_nsaids",
    'sym_heartburn_inferred',
    'pth_pTNM')
independent_permute_sensitivity <-
  c("dem_bmi_group", 
    "smk_ever_smoked", 
    "med_taken_nsaids", 
    'pth_pTNM')

independent_no_TNMStage <-
  c(
    "dem_age_group",
    "dem_gender",
    'dem_bmi_group',
    'smk_ever_smoked',
    'med_taken_nsaids',
    'sym_heartburn_inferred'
    )
independent_permute_no_TNMStage <-
  c("dem_bmi_group",
    "smk_ever_smoked",
    "med_taken_nsaids",
    'sym_heartburn_inferred'
)

BE_IM_EAC_vs_Non_BE_IM_EAC_logit_model_data %>%   
summary_factorlist(dependent,
                     c(independent_base, independent),
                     column = TRUE,
                     fit_id = TRUE) %>%
  # Univariable
  ff_merge(
    glmuni(
      BE_IM_EAC_vs_Non_BE_IM_EAC_logit_model_data,
      dependent,
      c(independent)
    ) %>%
      fit2df(estimate_suffix = " (Univariable)")
  ) %>%
  # Base
  ff_merge(
    glmmulti(
      BE_IM_EAC_vs_Non_BE_IM_EAC_logit_model_data,
      dependent,
      independent_base
    ) %>%
      fit2df(estimate_suffix = " (Base model [age+sex])")
  ) %>%
  # BMI Model
  ff_merge(
    glmmulti(
      BE_IM_EAC_vs_Non_BE_IM_EAC_logit_model_data,
      dependent,
      c(independent_base, independent_permute[1])
    ) %>%
      fit2df(estimate_suffix = " (BMI model)")
  ) %>%
    # Smoking Model
  ff_merge(
    glmmulti(
      BE_IM_EAC_vs_Non_BE_IM_EAC_logit_model_data,
      dependent,
      c(independent_base, independent_permute[2])
    ) %>%
      fit2df(estimate_suffix = " (Smoking model)")
  ) %>%
    # NSAID Model
  ff_merge(
    glmmulti(
      BE_IM_EAC_vs_Non_BE_IM_EAC_logit_model_data,
      dependent,
      c(independent_base, independent_permute[3])
    ) %>%
      fit2df(estimate_suffix = " (NSAID model)")
  ) %>%
    # Reflux Model
  ff_merge(
    glmmulti(
      BE_IM_EAC_vs_Non_BE_IM_EAC_logit_model_data,
      dependent,
      c(independent_base, independent_permute[4])
    ) %>%
      fit2df(estimate_suffix = " (Reflux model)")
  ) %>%
      # TNM Model
  ff_merge(
    glmmulti(
      BE_IM_EAC_vs_Non_BE_IM_EAC_logit_model_data,
      dependent,
      c(independent_base, independent_permute[5])
    ) %>%
      fit2df(estimate_suffix = " (TNM model)")
  ) %>%
    ff_merge(
    glmmulti(
      BE_IM_EAC_vs_Non_BE_IM_EAC_logit_model_data,
      dependent,
      c(independent_base, independent_permute)
    ) %>%
      fit2df(estimate_suffix = " (Fully adjusted model)")
  ) %>%
    ff_merge(
    glmmulti(
      BE_IM_EAC_vs_Non_BE_IM_EAC_logit_model_data,
      dependent, 
      c(independent_base, independent_no_nsaid)
    ) %>%
      fit2df(estimate_suffix = " (Fully adjusted model sensitivity: remove aspirin)")
  ) |> 
    ff_merge(
    glmmulti(
      BE_IM_EAC_vs_Non_BE_IM_EAC_logit_model_data,
      dependent, 
      c(independent_base, independent_permute_sensitivity)
    ) %>%
      fit2df(estimate_suffix = " (Fully adjusted model sensitivity: remove reflux)")
  ) %>%
  select(-c(fit_id, index)) %>%    
  kableExtra::kable(row.names = F) %>%
  kableExtra::kable_styling(full_width = F)
