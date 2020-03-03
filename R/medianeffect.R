# R package for Median Effect principle calculations and plots
#
# Requires:
# - knitr (for kable)
# - dplyr, ggplot2
#
# Class: drug_effects
#
# Single drug series
# - vector of doses (D)
# - vector of fraction affected (fa)
# - drug name (name) - full name of tested drug
# - label (label)    - label to use on plots
# - info (info)      - additional information
# - ratio: vector of ratio of each drug in constant ratio combinations
# - Calculated: Dm, m, R^2, R, log10(D), log10(fa / (1-fa))
#
# Combination drug series, NON-constant ratio
# - separate object type?
# - vector of drug 1 doses (D1)
# - vector of drug 2 doses (D2)
# - vector of fraction affected (fa)
# - drug name (name) - full name of tested drugs
# - label (label)    - label to use on plots
# - info (info)      - additional information


#' Drug effects object constructor
#'
#' Internal constructor function to generate a drug_effects object, for single drugs and constant ratio combinations.
#'
#' @param D Vector of doses (single drug dose or sum of doses for fixed-ratio combinations)
#' @param fa Vector of fraction affected
#' @param name Long name of drug (optional)
#' @param label Short label of drug (for plot labels, optional)
#' @param info Longer description (optional)
#' @param ratio If present, specifies vector of ratio of each drug in constant ratio combinations
#'
new_drug_effects <- function(D = double(), fa = double(), name = "", label = "", info = "", ratio = double()) {
  stopifnot(is.double(D), is.double(fa), is.double(ratio))
  class <- c("drug_effects")
  if (length(ratio) != 0) {
    # ratio should either be length 0 for single-drug, or else length = # of drugs
    if (length(ratio) == 1) {
      # allow shorthand of entering just ratio to refer to ratio:1 for 2-drug combination
      ratio <- c(ratio, 1.0)
    }
    class <- append("combo_drug_effects", class)
  }
  values <- list(D = D, fa = fa, name = name, label = label, info = info, ratio = ratio)
  structure(values, class = class)
}



#' Drug effects object validator + calculator
#'
#' Internal validator + calculator function.
#'
#' Checks that doses are non-negative, fraction affected is between (0,1),
#' equal length vector of at least two doses and effects, ratio positive (if present).
#'
#' Also calculates log(D) and log(fa / (1-fa)), fits a linear model, and calculates b, m, Dm (median effect dose), and
#' R^2 statistic (as well as R).
#'
#' Returns the drug_effects object augmented with the additional calculated information.
#'
#' @param x drug_effects object
#'
validate_drug_effects <- function(x) {
  D  <- x$D
  fa <- x$fa
  ratio <- x$ratio

  if (!all(!is.na(D) & D > 0.0)) { stop("All doses `D` must be non-missing and greater than zero", call. = FALSE) }
  if (!all(!is.na(fa) & fa > 0.0 & fa < 1.0)) { stop("All fraction affected `fa` must be non-missing, greater than zero, and less than one", call. = FALSE) }
  if (length(D) != length(fa)) { stop("The number of doses must equal the number of fraction affected", call. = FALSE) }
  if (length(D) < 2) { stop("There must be at least two dose / fraction affected points", call. = FALSE) }

  if (length(ratio) != 0) {
    if (!all(ratio > 0)) { stop("The drug ratio components must all be positive", call. = FALSE) }
  }

  # Do the statistical calculations; add more elements to x
  x$log_D     <- log10(x$D)
  x$log_fa_fu <- log10(x$fa / (1-x$fa))

  fit  <- stats::lm(log_fa_fu ~ log_D, data = data.frame(list(log_D = x$log_D, log_fa_fu = x$log_fa_fu)))

  x$b  <- unname(fit$coefficients[1])
  x$m  <- unname(fit$coefficients[2])
  x$Dm <- unname(10^(-x$b/x$m))
  x$R2 <- unname(summary(fit)$r.squared)
  x$R  <- unname(sqrt(x$R2)) * sign(x$m) # m should always be positive, but just in case

  x # return the validated and augmented object
}




#' Drug object generator
#'
#' Function to create a drug_effects object, either for a single drug or fixed-ratio combination of 2 drugs
#'
#' @param D Vector of doses (single drug dose or sum of doses for fixed-ratio combinations)
#' @param fa Vector of fraction affected
#' @param name Long name of drug (optional)
#' @param label Short label of drug (for plot labels, optional)
#' @param info Longer description (optional)
#' @param ratio If present, specifies vector of ratio of each drug in constant ratio combinations
#' @export
#' @examples
#' D  <- c(0.10, 0.15, 0.20, 0.25, 0.35)
#' fa <- c(0.24, 0.44, 0.63, 0.81, 0.90)
#' drug_effects(D = D, fa = fa, name = 'Example', label = 'drug1')
#'
drug_effects <- function(D = double(), fa = double(), name = "", label = "", info = "", ratio = double()) {
  D <- as.double(D)
  fa <- as.double(fa)
  name <- as.character(name)
  label <- as.character(label)
  info <- as.character(info)
  ratio <- as.double(ratio)
  validate_drug_effects(new_drug_effects(D, fa, name, label, info, ratio))
}



#' Calculate predicted fraction affected (fa)
#'
#' Given a drug_effects object and vector of doses, return predicted fraction affected
#'
#' @param drug Drug effect object
#' @param D Vector of doses
#' @export
#'
calc_fa <- function(drug, D) {
  if (!inherits(drug, "drug_effects")) { stop("Requires a `drug_effects` object", call. = FALSE) }
  m <- drug$m
  Dm <- drug$Dm
  return (1 / (1 + ((Dm / D)^m)))
}



#' Calculate predicted doses
#'
#' Given a drug_effects object and vector of fraction affected (fa), return predicted doses required
#'
#' @param drug Drug effect object
#' @param fa Vector of fraction affected (fa)
#' @export
#'
#
# [ ] To consider: what about drug combination objects? break down by ratio and return dataframe?
#
calc_d <- function(drug, fa) {
  if (!inherits(drug, "drug_effects")) { stop("Requires a `drug_effects` object", call. = FALSE) }
  m <- drug$m
  Dm <- drug$Dm
  return (Dm*(fa / (1 - fa))^(1/m))
}



#' Function to decompose a drug_combo drug effects object into its component drug doses
#'
#' @param drug_combo A fixed-ratio combination drug-effects object
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @export
drug_combo_decompose <- function(drug_combo) {
  if (!inherits(drug_combo, "combo_drug_effects")) { stop("Requires a fixed-ratio combination drug effects object", call. = FALSE) }

  D <- drug_combo$D
  fa <- drug_combo$fa
  ratio <- drug_combo$ratio
  df <- data.frame(list(D = D, fa = fa))

  f <- function(D, fa, ratio) {
    # Function for pmap to iterate over dataframe of rows of D and fa
    # - returns a WIDE format dataframe row with D, fa, and a column for each drug portion
    # - D and fa are a single number (because receiving one row at a time)
    # - ratio is a vector
    data.frame(D = D, fa = fa, dose_portion = D * ratio / sum(ratio)) %>%
      dplyr::mutate(drug = dplyr::row_number()) %>%
      tidyr::pivot_wider(names_from = .data$drug, values_from = .data$dose_portion, names_prefix = "drug_")
  }
  purrr::pmap_dfr(df, f, ratio)
}


#' Drug combination calculations
#'
#' Given drug combination object (either fixed or non-fixed ratio) and arbitrary number
#' of single-drug effect objects, calculate
#' all parameters needed for later calculation of combination or dose reduction index,
#' at either a vector of specified fa values or using actual doses / fa from drug_combo.
#' Returns unique id, total dose in combination, fa, single drug doses needed for fa,
#' and dose in the combination for each drug.
#'
#' @param drug_combo Drug effect fixed-ratio combination object
#' @param ...  Drug effect objects
#' @param fa Vector of fraction affected (fa) at which calculations will be made (optional)
#' @importFrom rlang .data
#' @export
#'
calc_combo <- function(drug_combo, ..., fa = double()) {

  # `inherits` seems to do an boolean OR for the vector of classes
  if (!inherits(drug_combo, c("combo_drug_effects", "ncr_drug_effects"))) { stop("Requires a combination drug effects object", call. = FALSE) }
  if (!all(purrr::map_lgl(list(...), inherits, "drug_effects"))) { stop("Some objects not drug effects objects", call. = FALSE) }

  if (inherits(drug_combo, "combo_drug_effects")) {
    if (length(list(...)) != length(drug_combo$ratio)) { stop("Different number of drugs in ratio and single drug effects objects", call. = FALSE) }
  } else {
    if (length(list(...)) != ncol(drug_combo$dose_df)) { stop("Different number of drugs in combination and single drug effects objects", call. = FALSE) }
    if (length(fa) != 0) { stop("Can not specify arbitrary fa levels for non-constant ratio drug combinations", call. = FALSE) }
  }

  # Generate dataframe of each single drug's m / Dm / label
  df_drugs <- list(...) %>%
    purrr::map_dfr(function(x) { data.frame(m = x$m, Dm = x$Dm, label = x$label, stringsAsFactors = FALSE) }) %>%
    dplyr::mutate(drug = paste0('drug_', dplyr::row_number())) %>%
    dplyr::mutate(label = ifelse(.data$label == '', .data$drug, .data$label)) # fill in missing labels


  # need id column in case fa is non-unique
  if (inherits(drug_combo, "combo_drug_effects")) { # constant-ratio specific code
    ratio <- drug_combo$ratio
    if (length(fa) == 0) { # if no fa values specified...
      dose_total <- drug_combo$D # ... use the actual doses
      fa <- drug_combo$fa        # ... use the actual fa
    } else {
      dose_total <- calc_d(drug_combo, fa) # ... otherwise, predict doses for the given fa
    }
    df <- data.frame(dose_total = dose_total, fa = fa) %>% dplyr::mutate(id = dplyr::row_number()) %>%
      tidyr::crossing(df_drugs) %>%
      dplyr::left_join( # generate column of proportion of each drug in the combination
        data.frame(proportion = ratio / sum(ratio)) %>% dplyr::mutate(drug = paste0('drug_', dplyr::row_number())),
        by = 'drug') %>%
      dplyr::mutate(
        dose_single = (.data$Dm*(.data$fa / (1 - .data$fa))^(1/.data$m)), # calculate predicted single drug alone dose
        dose_combo = .data$dose_total * .data$proportion
      )

  } else { # non-constant ratio specific code

    df <- cbind(id = seq(1, nrow(drug_combo$dose_df)), fa = drug_combo$fa, drug_combo$dose_df) %>%
      tidyr::pivot_longer(tidyselect::contains('drug_'), names_to = 'drug', values_to = 'dose_combo') %>%
      dplyr::left_join(df_drugs, by = 'drug') %>%
      dplyr::mutate(dose_single = (.data$Dm*(.data$fa / (1 - .data$fa))^(1/.data$m))) # calculate predicted single drug alone dose
  }

  df <- df %>% dplyr::select(.data$id, .data$fa, .data$drug, .data$label, .data$dose_combo, .data$dose_single, .data$m, .data$Dm)
  return(df)
}


#' Calculate combination index (CI)
#'
#' Given drug_combo and arbitrary number of single-drug effect objects, calculate
#' combination index at vector of fa values. If fa is empty, then use the actual fa
#' and actual doses in drug_combo.
#'
#' @param drug_combo Drug effect fixed-ratio combination object
#' @param ...  Drug effect objects
#' @param fa Vector of fraction affected (fa) at which calculations will be made (optional)
#' @importFrom rlang .data
#' @export
#'
calc_ci <- function(drug_combo, ..., fa = double()) {
  if (!inherits(drug_combo, "combo_drug_effects")) { stop("Requires a fixed-ratio combination drug effects object", call. = FALSE) }
  if (!all(purrr::map_lgl(list(...), inherits, "drug_effects"))) { stop("Some objects not drug effects objects", call. = FALSE) }
  if (length(list(...)) != length(drug_combo$ratio)) { stop("Different number of drugs in ratio and single drug effects objects", call. = FALSE) }

  calc_combo(drug_combo = drug_combo, ..., fa = fa) %>%
    dplyr::group_by(.data$fa, .data$id) %>%
    dplyr::summarize(CI = sum(.data$dose_combo / .data$dose_single)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.data$id)
}





#' Calculate dose reduction index (DRI)
#'
#' Given drug_combo and arbitrary number of single-drug effect objects, calculate
#' the dose reduction index for each drug at a vector of fa values. If fa is empty,
#' then use the actual fa and actual doses in drug_combo.
#'
#' @param drug_combo Drug effect fixed-ratio combination object
#' @param ...  Drug effect objects
#' @param fa Vector of fraction affected (fa) at which calculations will be made (optional)
#' @importFrom rlang .data
#' @export
#'
calc_dri <- function(drug_combo, ..., fa = double()) {
  if (!inherits(drug_combo, "combo_drug_effects")) { stop("Requires a fixed-ratio combination drug effects object", call. = FALSE) }
  if (!all(purrr::map_lgl(list(...), inherits, "drug_effects"))) { stop("Some objects not drug effects objects", call. = FALSE) }
  if (length(list(...)) != length(drug_combo$ratio)) { stop("Different number of drugs in ratio and single drug effects objects", call. = FALSE) }

  calc_combo(drug_combo = drug_combo, ..., fa = fa) %>%
    dplyr::group_by(.data$fa, .data$id, .data$drug) %>%
    dplyr::summarise(dri = .data$dose_single / .data$dose_combo) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_wider(names_from = .data$drug, values_from = .data$dri, names_prefix = 'dri_') %>%
    dplyr::select(-.data$id)
}



#' Print method for objects of drug_effects class
#'
#' @param x A drug_effects object
#' @param ... (not used)
#' @param stats Boolean for whether to display calculated m, Dm, R^2, and R
#' @importFrom dplyr %>%
#' @export
#'
print.drug_effects <- function(x, ..., stats = TRUE) {

  D <- x$D
  fa <- x$fa
  name <- x$name
  label <- x$label
  if (label != '') { label <- paste0(' (', label, ')') }
  info <- x$info
  if (info != '') { info <- paste0('\n\n', info) }
  ratio <- x$ratio
  Dm <- x$Dm

  cat(name, label, info, sep = '')

  df <- data.frame(list(D = D, fa = fa))

  # Additional processing for constant ratio combination drug effects
  if (inherits(x, "combo_drug_effects")) {
    cat('\nDrug ratio:', paste(signif(ratio / min(ratio), 4), collapse = " : "))
    f <- function(D, fa, ratio) {
      # Function for pmap to iterate over dataframe of rows of D and fa
      # - returns a WIDE format dataframe row with D, fa, and a column for each drug portion
      # - D and fa are a single number (because receiving one row at a time)
      # - ratio is a vector
      data.frame(D = D, fa = fa, dose_portion = signif(D * ratio / sum(ratio), 4)) %>%
        dplyr::mutate(drug = dplyr::row_number()) %>%
        tidyr::pivot_wider(names_from = .data$drug, values_from = .data$dose_portion, names_prefix = "drug_")
    }
    df <- purrr::pmap_dfr(df, f, ratio)
    Dm <- paste0(signif(Dm, 4), ' = ', paste0(signif(Dm * ratio / sum(ratio), 4), collapse = ' + '))
  }

  print(knitr::kable(df))
  if (stats) { cat('\nm: ', x$m, '\nDm: ', Dm, '\nR2: ', x$R2, '\nR: ', x$R, sep = '') }
}




#' Plot method for objects of drug_effects class
#'
#' @param x A drug_effects object
#' @param y (not used)
#' @param ... (not used)
#' @param color Line color
#' @importFrom rlang .data
#' @export
#'
plot.drug_effects <- function(x, y, ..., color = 'blue') {
  title <- "Median-Effect Plot"
  if (x$name != "") { title <- paste0(title, ": ", x$name) }
  if (x$label != "") { title <- paste0(title, " (", x$label, ")") }

  df <- data.frame(list(D = x$D, fa = x$fa, log_D = x$log_D, log_fa_fu = x$log_fa_fu))

  g <- ggplot2::ggplot(data = df, ggplot2::aes(.data$log_D, .data$log_fa_fu)) +
    ggplot2::geom_point() +
    ggplot2::stat_smooth(method = "lm", color = color) +
    ggplot2::xlab("log (D)") +
    ggplot2::ylab("log (fa / fu)") +
    ggplot2::ggtitle(title)

  g
}





#' Generate Median-Effect Plot
#'
#' Generates Median-Effect plots (as ggplot objects) for an arbitrary number of drug_effects objects
#'
#' Plots linearized log(fa / (1-fa)) versus log(D)
#'
#' @param ... drug_effects objects to plot
#' @importFrom rlang .data
#' @export
median_effect_plot <- function(...) {
  if (!all(purrr::map_lgl(list(...), inherits, "drug_effects"))) { stop("Some objects not drug effects objects", call. = FALSE) }

  df <- data.frame()
  i <- 0
  for (d in list(...)) {
    i <- i + 1
    if (d$label == '') { d$label <- paste('Drug', i) } # assign labels if empty
    df <- rbind(df, data.frame(list(log_D = d$log_D, log_fa_fu = d$log_fa_fu, label = d$label), stringsAsFactors = FALSE))
  }
  df$label <- factor(df$label, levels = unique(df$label)) # prevent re-ordering of labels; `unique` appears to maintain order
  ggplot2::ggplot(data = df, ggplot2::aes(.data$log_D, .data$log_fa_fu, shape = .data$label, color = .data$label)) +
    ggplot2::geom_point() +
    ggplot2::stat_smooth(method = "lm", se = FALSE, fullrange = TRUE) + # fullrange allows extrapolation of line beyond the data points
    ggplot2::xlab("log (D)") +
    ggplot2::ylab("log (fa / fu)") +
    ggplot2::labs(color = 'Drug', shape = 'Drug')
}



#' Generate Dose-Effect Plot
#'
#' Generates Dose-Effect plots (as ggplot objects) for an arbitrary number of drug_effects objects
#'
#' Plots curves of fa versus D.
#'
#' @param ... drug_effects objects to plot
#' @param from Start of range of fraction affected (fa)
#' @param to End of range of fa
#' @param by Step size of fa
#' @importFrom rlang .data
#' @importFrom dplyr %>%
#' @export
dose_effect_plot <- function(..., from = 0.01, to = 0.99, by = 0.01) {
  if (!all(purrr::map_lgl(list(...), inherits, "drug_effects"))) { stop("Some objects not drug effects objects", call. = FALSE) }

  df_lines <- data.frame()
  df_points <- data.frame()

  i <- 0
  for (d in list(...)) {
    i <- i + 1
    if (d$label == '') { d$label <- paste('Drug', i) } # assign labels if empty
    df_lines  <- rbind(df_lines,  data.frame(list(label = d$label, m = d$m, Dm = d$Dm))) # dataframe with label, m, Dm
    df_points <- rbind(df_points, data.frame(list(label = d$label, D = d$D, fa = d$fa))) # dataframe with label, D, fa
  }
  labels <- unique(df_lines$label) # prevent re-ordering of labels; `unique` appears to maintain order
  df_lines$label  <- factor(df_lines$label,  levels = labels)
  df_points$label <- factor(df_points$label, levels = labels)

  # Generate curves from m, Dm in df_lines
  # - generate a cartesian product with all fa levels, for each drug
  df_lines <- df_lines %>%
    tidyr::crossing(data.frame(fa = seq(from = from, to = to, by = by))) %>%
    dplyr::mutate(D = .data$Dm*(.data$fa / (1 - .data$fa))^(1/.data$m))

  g <- ggplot2::ggplot(df_lines, ggplot2::aes(.data$D, .data$fa, color = .data$label)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(data = df_points, ggplot2::aes(.data$D, .data$fa, color = .data$label, shape = .data$label)) +
    ggplot2::xlab("Dose") +
    ggplot2::ylab("Fraction affected (fa)") +
    ggplot2::labs(color = 'Drug', shape = 'Drug')
  g
}



#' Generate fa-CI Plot
#'
#' Generates fa-CI plots (as ggplot objects) for a fixed-ratio combination
#'
#' @param drug_combo Drug effect fixed-ratio combination object
#' @param ...  Drug effect objects
#' @param from Start of range of fraction affected (fa)
#' @param to End of range of fa
#' @param by Step size of fa
#' @importFrom rlang .data
#' @export
fa_ci_plot <- function(drug_combo, ..., from = 0.01, to = 0.99, by = 0.01) {
  if (!inherits(drug_combo, "combo_drug_effects")) { stop("Requires a fixed-ratio combination drug effects object", call. = FALSE) }
  if (!all(purrr::map_lgl(list(...), inherits, "drug_effects"))) { stop("Some objects not drug effects objects", call. = FALSE) }
  if (length(list(...)) != length(drug_combo$ratio)) { stop("Different number of drugs in ratio and single drug effects objects", call. = FALSE) }

  g <- calc_ci(drug_combo, ..., fa = seq(from, to, by)) %>%
    ggplot2::ggplot(ggplot2::aes(.data$fa, .data$CI)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(data = calc_ci(drug_combo, ...), ggplot2::aes(.data$fa, .data$CI)) +
    ggplot2::geom_hline(yintercept = 1.0, linetype = 'dotted') +
    ggplot2::xlab("fa") +
    ggplot2::ylab("CI") +
    ggplot2::coord_cartesian(xlim = c(0,1))
  g
}


#' Generate fa-DRI Plot
#'
#' Generates fa-DRI plots (as ggplot objects) for a fixed-ratio combination
#'
#' @param drug_combo Drug effect fixed-ratio combination object
#' @param ...  Drug effect objects
#' @param from Start of range of fraction affected (fa)
#' @param to End of range of fa
#' @param by Step size of fa
#' @importFrom rlang .data
#' @export
fa_dri_plot <- function(drug_combo, ..., from = 0.01, to = 0.99, by = 0.01) {
  if (!inherits(drug_combo, "combo_drug_effects")) { stop("Requires a fixed-ratio combination drug effects object", call. = FALSE) }
  if (!all(purrr::map_lgl(list(...), inherits, "drug_effects"))) { stop("Some objects not drug effects objects", call. = FALSE) }
  if (length(list(...)) != length(drug_combo$ratio)) { stop("Different number of drugs in ratio and single drug effects objects", call. = FALSE) }

  # This next chunk of code is an ugly mess and should be refactored
  # Create labels (if missing) and maintain order
  # Store in a datframe to join the results of calc_dri from, which don't have labels
  df_labels <- data.frame()
  i <- 0
  for (d in list(...)) {
    i <- i + 1
    if (d$label == '') { d$label <- paste0('drug_', i) } # assign labels if empty
    df_labels  <- rbind(df_labels,  data.frame(list(drug = paste0('drug_', i), label = d$label), stringsAsFactors = FALSE))
  }
  labels <- unique(df_labels$label) # prevent re-ordering of labels; `unique` appears to maintain order
  df_labels$label  <- factor(df_labels$label,  levels = labels)


  df_points <- calc_dri(drug_combo, ...) %>%
    tidyr::pivot_longer(names_to = 'drug', values_to = 'DRI', cols = tidyselect::contains('dri_drug_'), names_prefix = 'dri_') %>%
    dplyr::left_join(df_labels, by = 'drug')

  df_lines <- calc_dri(drug_combo, ..., fa = seq(from, to, by)) %>%
    tidyr::pivot_longer(names_to = 'drug', values_to = 'DRI', cols = tidyselect::contains('dri_drug_'), names_prefix = 'dri_') %>%
    dplyr::left_join(df_labels, by = 'drug')

  g <- df_lines %>%
    ggplot2::ggplot(ggplot2::aes(.data$fa, .data$DRI, color = .data$label)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(data = df_points, ggplot2::aes(.data$fa, .data$DRI, color = .data$label, shape = .data$label)) +
    ggplot2::xlab("fa") +
    ggplot2::ylab("DRI") +
    ggplot2::labs(color = 'Drug', shape = 'Drug')
  g
}



#' Non-constant ratio drug effects object generator
#'
#' Function to generate a non-constant ratio drug_effects object
#'
#' @param ... Vectors of doses, one for each drug, with length equal to number of data points
#' @param fa Vector of fraction affected
#' @param name Name to describe set (optional)
#' @param label Short label (optional)
#' @param info Longer description (optional)
#' @export
#'
ncr_drug_effects <- function(..., fa = double(), name = "", label = "", info = "") {

  ### Validation ###
  # - not checking for edge case of all doses equal to 0...
  if (length(list(...)) == 0) { stop("Require at least one dose and effect point", call. = FALSE) }
  if (!all(purrr::map_lgl(list(...), is.double))) { stop("Doses must be numeric", call. = FALSE) }
  if (!all(purrr::map_lgl(list(...), function(x) {all(x >= 0)} ))) { stop("All doses must be greater than or equal to zero", call. = FALSE) }

  # Check that all the dose vectors are of equal length
  lengths <- purrr::map_int(list(...), length)
  n <- min(lengths)
  if (n != max(lengths)) { stop("Number of doses for each drug not all equal", call. = FALSE) }

  # Check fa
  if (!is.double(fa)) { stop("Fraction affected must be numeric", call. = FALSE) }
  if (length(fa) != n) { stop("Number of fraction affected must equal number of doses", call. = FALSE) }
  if (!all(fa > 0.0 & fa < 1.0)) { stop("All fraction affected `fa` must be greater than zero, and less than one", call. = FALSE) }

  ### Generation ###

  # Convert vectors of doses to dataframe
  # - one row per point
  # - one column per drug, named drug_#
  dose_df <- data.frame(matrix(NA, nrow = n, ncol = length(list(...))))
  i <- 0
  for (d in list(...)) {
    i <- i + 1
    dose_df[i] <- d
  }
  names(dose_df) <- paste0('drug_', seq(1, length(list(...))))

  values <- list(dose_df = dose_df, fa = fa, name = name, label = label, info = info)
  structure(values, class = "ncr_drug_effects")
}



#' Print method for non-constant ratio drug_effects class
#'
#' @param x A non-constant ratio drug_effects object
#' @param ... (not used)
#' @export
#'
print.ncr_drug_effects <- function(x, ...) {
  name <- x$name
  label <- x$label
  if (label != '') { label <- paste0(' (', label, ')') }
  info <- x$info
  if (info != '') { info <- paste0('\n\n', info) }

  cat(name, label, info, sep = '')
  print(knitr::kable(cbind(fa = x$fa, x$dose_df)))
}




#' Calculate combination index (CI) for non-constant ratio
#'
#' Given non-constant ratio drug effects object and arbitrary number of single-drug
#' effect objects, calculate combination index observed fa values.
#'
#' @param ncr_combo Non-constant ratio combination drug effects object
#' @param ...  Drug effect objects
#' @importFrom rlang .data
#' @export
#'
ncr_calc_ci <- function(ncr_combo, ...) {

  if (!inherits(ncr_combo, "ncr_drug_effects")) { stop("Requires a non-constant ratio combination drug effects object", call. = FALSE) }
  if (!all(purrr::map_lgl(list(...), inherits, "drug_effects"))) { stop("Some objects not drug effects objects", call. = FALSE) }
  if (length(list(...)) != ncol(ncr_combo$dose_df)) { stop("Different number of drugs in combination and single drug effects objects", call. = FALSE) }

  calc_combo(ncr_combo, ...) %>%
    dplyr::group_by(.data$fa, .data$id) %>%
    dplyr::summarize(CI = sum(.data$dose_combo / .data$dose_single)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.data$id)
}
