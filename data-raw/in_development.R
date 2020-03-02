library(tidyverse)



# Isobologram

df <- calc_combo(combo_1_2, drug1, drug2)

df2 <- ncr_calc_combo(ncr_combo, drug1, drug2)

df %>%
  mutate(drug = paste0('drug_', drug)) %>%
  select(-D_combo, -label) %>%
  pivot_wider(id_cols = c(id, fa), names_from = drug, values_from = c(dose_single, dose_combo))
