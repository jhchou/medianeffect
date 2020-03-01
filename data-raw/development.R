library(tidyverse)
library(medianeffect)

# Chou, T. C., R. J. Motzer, Y. Tong, and G. J. Bosl. 1994.
# “Computerized Quantitation of Synergism and Antagonism of Taxol, Topotecan, and Cisplatin against Human Teratocarcinoma Cell Growth: A Rational Approach to Clinical Protocol Design.”
# Journal of the National Cancer Institute 86 (20): 1517–24. https://doi.org/10.1093/jnci/86.20.1517.


# Paclitaxel
drug1 <- drug_effects(
  D = c(0.002, 0.004, 0.005, 0.01, 0.02),
  fa = c(0.429, 0.708, 0.761, 0.882, 0.932),
  name = "Paclitaxel",
  label = "Paclitaxel")
drug1


# Cis-Pt
drug2 <- drug_effects(
  D = c(0.05, 0.1, 0.2, 0.5, 1, 2),
  fa = c(0.055, 0.233, 0.301, 0.559, 0.821, 0.953),
  name = "Cisplatin",
  label = "cis-pt")
drug2


# Topotecan
drug3 <- drug_effects(
  D = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5),
  fa = c(0.069, 0.213, 0.373, 0.785, 0.94, 0.991),
  name = 'Topotecan',
  label = 'Topotecan'
)
drug3


# Paclitaxel:Cis-Pt Combo, ratio = 0.01
combo_1_2 <- drug_effects(
  D = c(0.001, 0.002, 0.005, 0.01) + c(0.1, 0.2, 0.5, 1),
  fa = c(0.45, 0.701, 0.910, 0.968),
  name = 'Paclitaxel - Cisplatin',
  label = 'combo_1_2',
  ratio = c(0.001, 0.1))
combo_1_2


# Cis-Pt:Topotecan Combo, ratio = 0.1
combo_2_3 <- drug_effects(
  D = c(0.05, 0.1, 0.2, 0.5, 1) + c(0.005, 0.01, 0.02, 0.05, 0.1),
  fa = c(0.304, 0.413, 0.675, 0.924, 0.977),
  name = 'Cisplatin - Topotecan',
  label = 'combo_2_3',
  ratio = c(0.1, 0.01))
combo_2_3


# Paclitaxel:Topotecan Combo, ratio = 0.1
combo_1_3 <- drug_effects(
  D = c(0.001, 0.002, 0.005, 0.01) + c(0.01, 0.02, 0.05, 0.1),
  fa = c(0.274, 0.579, 0.901, 0.965),
  name = 'Paclitaxel - Topotecan',
  label = 'combo_1_3',
  ratio = c(0.01, 0.1))
combo_1_3


# Paclitaxel:Cis-Pt:Topotecan Combo, ratio = 0.001 : 0.1 : 0.01
combo_1_2_3 <- drug_effects(
  D = c(0.001, 0.002, 0.003, 0.005) + c(0.1, 0.2, 0.3, 0.5) + c(0.01, 0.02, 0.03, 0.05),
  fa = c(0.456, 0.806, 0.947, 0.995),
  name = 'Paclitaxel - Cisplatin - Topotecan',
  label = 'combo_1_2_3',
  ratio = c(0.001, 0.1, 0.01))
combo_1_2_3

median_effect_plot(drug1, drug2, drug3, combo_1_2, combo_1_3, combo_2_3, combo_1_2_3)

median_effect_plot(drug1, drug2, combo_1_2)

dose_effect_plot(drug1, drug2, drug3, combo_1_2, combo_1_3, combo_2_3,combo_1_2_3)
dose_effect_plot(drug1, drug2, drug3, combo_1_2, combo_1_3, combo_2_3,combo_1_2_3) +
  coord_cartesian(xlim = c(0, 2.1))

fa_ci_plot(combo_1_2_3, drug1, drug2, drug3)

fa_dri_plot(combo_1_2_3, drug1, drug2, drug3)

drug1

drug2

drug3

combo_1_2
calc_CI(combo_1_2, drug1, drug2)

combo_2_3
calc_CI(combo_2_3, drug2, drug3)

combo_1_3
calc_CI(combo_1_3, drug1, drug3)

combo_1_2_3
calc_CI(combo_1_2, drug1, drug2, fa = c(0.5, 0.75, 0.9, 0.95))

calc_CI(combo_1_2_3, drug1, drug2, drug3)
calc_CI(combo_1_2_3, drug1, drug2, drug3, fa = c(0.5, 0.75, 0.9, 0.95))

calc_DRI(combo_1_2_3, drug1, drug2, drug3, fa = c(0.5, 0.75, 0.9, 0.95))



f <- function(...) {
  all(map_lgl(list(...), inherits, "drug_effects"))
}
