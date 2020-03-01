#' Drug combination calculations
#'
#' Given drug_combo and arbitrary number of single-drug effect objects, calculate
#' all parameters needed for later calculation of combination or dose reduction index,
#' at either a vector of specified fa values or using actual doses / fa from drug_combo.
#' Returns unique id, total dose in combination, fa, single drug doses needed for fa,
#' and dose in the combination for each drug.
#'
#' @param drug_combo Drug effect fixed-ratio combination object
#' @param ...  Drug effect objects
#' @param fa Vector of fraction affected (fa) at which calculations will be made (optional)
#' @export
#'
calc_combo <- function(drug_combo, ..., fa = double()) {

  if (!inherits(drug_combo, "combo_drug_effects")) { stop("Requires a fixed-ratio combination drug effects object", call. = FALSE) }
  # if (!all(inherits(c(...), "drug_effects"))) { stop("Requires drug effects objects", call. = FALSE) }

  drug_effect_list <- function(...) {
    map_dfr(list(...), function(x) { data.frame(m = x$m, Dm = x$Dm) }) %>%
      mutate(drug = row_number())
  }

  ratio <- drug_combo$ratio

  # Generate dataframe of drug combination total doses and fa
  # - columns: D_combo, fa, id
  # - including id in case fa is non-unique for the actual dose / Fa
  if (length(fa) == 0) { # use the actual doses and effects
    D_combo <- drug_combo$D # actual doses
    fa <- drug_combo$fa     # the actual fa observed in the combo
  } else { # calculate the predicted doses for the list of given fa
    D_combo <- calc_D(drug_combo, fa)
  }
  df_combo <- data.frame(D_combo = D_combo, fa = fa) %>% mutate(id = row_number())

  # Generate dataframe of each drug of combination and m / Dm
  df_drugs <- drug_effect_list(...)

  df <- crossing(df_combo, df_drugs) %>%
    left_join( # generate column of proportion of each drug in the combination
      data.frame(proportion = ratio / sum(ratio)) %>% mutate(drug = row_number()),
      by = 'drug') %>%
    mutate(
      dose_single = (Dm*(fa / (1 - fa))^(1/m)), # calculate predicted single drug alone dose
      dose_combo = D_combo * proportion
    ) %>%
    select(
      id, D_combo, fa, drug, dose_single, dose_combo
    )

  return(df)
}



calc_CI <- function(drug_combo, ..., fa = double()) {
   calc_combo(drug_combo = drug_combo, ..., fa = fa) %>%
      group_by(D_combo, fa, id) %>%
      summarize(CI = sum(dose_combo / dose_single)) %>%
      ungroup() %>%
      select(-id)
}
# To generate fa-CI plot for both curve and actual points
# calc_CI(drug_combo, drug1, drug2, fa = seq(0.05, 0.95, 0.01)) %>%
#   ggplot(aes(fa, CI)) +
#   geom_line() +
#   geom_point(data = calc_CI(drug_combo, drug1, drug2), aes(fa, CI))



calc_DRI <- function(drug_combo, ..., fa = double()) {
   calc_combo(drug_combo = drug_combo, ..., fa = fa) %>%
    group_by(fa, id, drug) %>%
    summarise(dri = dose_single / dose_combo) %>%
    ungroup() %>%
    pivot_wider(names_from = drug, values_from = dri, names_prefix = 'dri_drug_') %>%
    select(-id)
}

# To generate fa-DRI plot for both curve and actual points
# calc_DRI(drug_combo, drug1, drug2, fa = seq(0.05, 0.95, 0.01)) %>%
#   pivot_longer(names_to = 'drug', values_to = 'DRI', cols = contains('dri_drug_'), names_prefix = 'dri_') %>%
#   ggplot(aes(fa, DRI, color = drug)) + geom_line() +
#   geom_point(data = calc_DRI(drug_combo, drug1, drug2) %>%
#                pivot_longer(names_to = 'drug', values_to = 'DRI', cols = contains('dri_drug_'), names_prefix = 'dri_'),
#              aes(fa, DRI, color = drug)
#   )


