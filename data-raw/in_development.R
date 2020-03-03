library(tidyverse)

# > calc_combo(combo_1_2, drug1, drug2)
#
#      id    fa drug   label      dose_combo dose_single     m      Dm
#
# 1     1 0.45  drug_1 Paclitaxel      0.001     0.00185  1.25 0.00217
# 2     1 0.45  drug_2 cis-pt          0.1       0.279    1.46 0.320
# 3     2 0.701 drug_1 Paclitaxel      0.002     0.00429  1.25 0.00217
# 4     2 0.701 drug_2 cis-pt          0.2       0.574    1.46 0.320
# 5     3 0.91  drug_1 Paclitaxel      0.005     0.0138   1.25 0.00217
# 6     3 0.91  drug_2 cis-pt          0.5       1.56     1.46 0.320
# 7     4 0.968 drug_1 Paclitaxel      0.01      0.0333   1.25 0.00217
# 8     4 0.968 drug_2 cis-pt          1         3.31     1.46 0.320
#
#
# > ncr_calc_combo(ncr_combo, drug1, drug2)
#
#      id    fa drug   label      dose_combo dose_single     m      Dm
#
# 1     1 0.45  drug_1 Paclitaxel      0.001     0.00185  1.25 0.00217
# 2     1 0.45  drug_2 cis-pt          0.1       0.279    1.46 0.320
# 3     2 0.701 drug_1 Paclitaxel      0.002     0.00429  1.25 0.00217
# 4     2 0.701 drug_2 cis-pt          0.2       0.574    1.46 0.320
# 5     3 0.91  drug_1 Paclitaxel      0.005     0.0138   1.25 0.00217
# 6     3 0.91  drug_2 cis-pt          0.5       1.56     1.46 0.320
# 7     4 0.968 drug_1 Paclitaxel      0.01      0.0333   1.25 0.00217
# 8     4 0.968 drug_2 cis-pt          1         3.31     1.46 0.320

# Isobologram

df <- calc_combo(combo_1_2, drug1, drug2) %>%
  select(-label, -m, -Dm)

# Generate the points on the axis to be connected by lines:
# - brute force separate the drug_1 and drug_2 single doses to separate rows, with coordinates (drug_1, 0) and (0, drug_2)
# - keeping fa, as each line will represent the equipotency line for the isobologram
df2 <- df %>%
  select(id, fa, drug, dose_single) %>%
  pivot_wider(names_from = drug, values_from = dose_single)

df_axis <-bind_rows(df2 %>% select(id, fa, drug_1), df2 %>% select(id, fa, drug_2)) %>%
  replace_na(list(drug_1 = 0, drug_2 = 0)) %>%
  mutate(fa = as.character(fa))



df_points <- df %>%
  select(id, fa, drug, dose_combo) %>%
  pivot_wider(names_from = drug, values_from = dose_combo) %>%
  mutate(fa = as.character(fa))


df_axis %>% ggplot(aes(drug_1, drug_2, color = fa)) +
  geom_line() +
  geom_point(data = df_points)


# normalized isobologram
df_points %>%
  left_join(df2 %>% transmute(id, d1 = drug_1, d2 = drug_2), by = 'id') %>%
  mutate(drug_1 = drug_1 / d1, drug_2 = drug_2 / d2) %>%
  ggplot(aes(drug_1, drug_2)) +
    geom_point(aes(color = fa)) +
    geom_line(data = data.frame(drug_1 = c(1, 0), drug_2 = c(0, 1)), aes(drug_1, drug_2)) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0,1))


# pivot_wider(id_cols = c(id, fa), names_from = drug, values_from = c(dose_single, dose_combo))
