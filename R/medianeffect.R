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
# - ratio: drug 1 / drug 2 (ratio), only defined for constant ratio drug combinations
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
#' @param ratio If present, specify constant ratio of drug1/drug2 doses
#'
new_drug_effects <- function(D = double(), fa = double(), name = "", label = "", info = "", ratio = double()) {
  stopifnot(is.double(D), is.double(fa), is.double(ratio))
  values <- list(D = D, fa = fa, name = name, label = label, info = info, ratio = ratio)
  class <- c("drug_effects")
  if (length(ratio) != 0) { class <- append("combo_drug_effects", class) }
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
    if (length(ratio) != 1) { stop("Only a single constant ratio can be used", call. = FALSE) }
    if (ratio < 0) { stop("The ratio of drug1 / drug2 must be non-negative", call. = FALSE) }
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
#' @param ratio If present, specify constant ratio of drug1/drug2 doses
#' @export
#' @examples
#' D  <- c(0.10, 0.15, 0.20, 0.25, 0.35)
#' fa <- c(0.24, 0.44, 0.63, 0.81, 0.90)
#' drug_effects(D = D, fa = fa, name = 'Example', label = 'drug1')
#'
drug_effects <- function(D = double(), fa = double(), name = "", label = "", info = "", ratio = double()) {
  D <- as.double(D)
  fa = as.double(fa)
  name = as.character(name)
  label = as.character(label)
  info = as.character(info)
  ratio = as.double(ratio)
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
calc_D <- function(drug, fa) {
  if (!inherits(drug, "drug_effects")) { stop("Requires a `drug_effects` object", call. = FALSE) }
  m <- drug$m
  Dm <- drug$Dm
  return (Dm*(fa / (1 - fa))^(1/m))
}



#' Calculate combination index (CI)
#'
#' Given drug 1, drug 2, and drug_combo objects, calculate combination index CI at vector of fa values.
#' If fa is empty, then use the actual fa and actual doses in drug_combo.
#'
#' @param drug1 Drug effect object
#' @param drug2 Drug effect object
#' @param drug_combo Drug effect fixed-ratio combination object
#' @param fa Vector of fraction affected (fa) at which combination indices (CI) will be calculated (optional)
#' @export
#'
calc_CI <- function(drug1, drug2, drug_combo, fa = double()) {
  # [ ] should do class checking here
  ratio <- drug_combo$ratio

  if (length(fa) == 0) {
    fa <- drug_combo$fa # the actual fa observed in the combo
    D_combo <- drug_combo$D # actual doses
  } else {
    D_combo <- calc_D(drug_combo, fa) # calculate predicted doses
  }
  D1 <- D_combo * (ratio / (ratio + 1))
  D2 <- D_combo / (ratio + 1)

  Dx1 <- calc_D(drug1, fa)
  Dx2 <- calc_D(drug2, fa)

  CI <- (D1 / Dx1) + (D2 / Dx2)
  return(CI)
}


#' Print method for objects of drug_effects class
#'
#' @param x A drug_effects object
#' @param ... (not used)
#' @param stats Boolean for whether to display calculated m, Dm, R^2, and R
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

  # If want to also include columns of log(D) and log(fa / fu)
  # df <- data.frame(list(D = D, fa = fa, log_D = x$log_D, log_fa_fu = x$log_fa_fu))

  df <- data.frame(list(D = D, fa = fa))
  if (length(ratio) != 0) {
    cat("\n\nRatio (drug 1/drug 2): ", ratio, sep = '')
    D1 <- (ratio / (ratio+1)) * D
    D2 <- D / (ratio + 1)
    df$D1 <- D1
    df$D2 <- D2
    # For combination, also show the contributing doses from drug 1 + 2, and not just the sum
    Dm <- paste0(signif(Dm, 4), ' = ', signif((ratio / (ratio+1)) * Dm, 4), ' + ', signif(Dm / (ratio + 1), 4))
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

  g <- ggplot2::ggplot(
    data = df,
    ggplot2::aes(.data$log_D, .data$log_fa_fu)
  ) +
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
  df <- data.frame()
  i <- 0
  for (d in list(...)) {
    if (inherits(d, "drug_effects")) {
      i <- i + 1
      if (d$label == '') { d$label <- paste('Drug', i) } # assign labels if empty
      df <- rbind(df, data.frame(list(log_D = d$log_D, log_fa_fu = d$log_fa_fu, label = d$label), stringsAsFactors = FALSE))
    }
  }
  df$label <- factor(df$label, levels = unique(df$label)) # prevent re-ordering of labels; `unique` appears to maintain order
  ggplot2::ggplot(data = df, ggplot2::aes(.data$log_D, .data$log_fa_fu, shape = .data$label, color = .data$label)) +
    ggplot2::geom_point() +
    ggplot2::stat_smooth(method = "lm", se = FALSE, fullrange = TRUE) # fullrange allows extrapolation of line beyond the data points
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
#' @export
dose_effect_plot <- function(..., from = 0.01, to = 0.99, by = 0.01) {
  # takes an arbitrary number of drug_effect objects
  df_lines <- data.frame()
  df_points <- data.frame()

  i <- 0
  for (d in list(...)) {
    if (inherits(d, "drug_effects")) {
      i <- i + 1
      if (d$label == '') { d$label <- paste('Drug', i) } # assign labels if empty
      df_lines  <- rbind(df_lines,  data.frame(list(label = d$label, m = d$m, Dm = d$Dm))) # dataframe with label, m, Dm
      df_points <- rbind(df_points, data.frame(list(label = d$label, D = d$D, fa = d$fa))) # dataframe with label, D, fa
    }
  }
  labels <- unique(df_lines$label) # prevent re-ordering of labels; `unique` appears to maintain order
  df_lines$label  <- factor(df_lines$label,  levels = labels)
  df_points$label <- factor(df_points$label, levels = labels)

  # Generate curves from m, Dm in df_lines
  # - generate a cartesian product with all fa levels, for each drug, by doing a join on a dummy variable (probably a better way)
  df_lines$dummy <- 'dummy'
  df2 <- data.frame(list(fa = seq(from = from, to = to, by = by), dummy = 'dummy'), stringsAsFactors = FALSE)
  df_lines <- dplyr::left_join(df_lines, df2, by = 'dummy')
  df_lines <- dplyr::mutate(df_lines, D = .data$Dm*(.data$fa / (1 - .data$fa))^(1/.data$m))

  g <- ggplot2::ggplot(df_lines, ggplot2::aes(.data$D, .data$fa, color = .data$label)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(data = df_points, ggplot2::aes(.data$D, .data$fa, color = .data$label, shape = .data$label))
  g
}



#' Generate fa-CI Plot
#'
#' Generates fa-CI plots (as ggplot objects) for a two-drug fixed-ratio combination
#'
#' @param drug1 Drug effect object
#' @param drug2 Drug effect object
#' @param drug_combo Drug effect fixed-ratio combination object
#' @param from Start of range of fraction affected (fa)
#' @param to End of range of fa
#' @param by Step size of fa
#' @export
fa_ci_plot <- function(drug1, drug2, drug_combo, from = 0.01, to = 0.99, by = 0.01) {

  fa <- seq(from = from, to = to, by = by)
  CI = calc_CI(drug1, drug2, drug_combo, fa)
  df <- data.frame(list(fa = fa, CI = CI))

  df_points <- data.frame(list(fa = drug_combo$fa, CI = calc_CI(drug1, drug2, drug_combo)))

  g <- ggplot2::ggplot(df, ggplot2::aes(fa, CI)) +
    ggplot2::geom_line() +
    ggplot2::geom_hline(yintercept = 1.0, linetype = 'dotted') +
    ggplot2::geom_point(data = df_points, ggplot2::aes(fa, CI)) +
    ggplot2::coord_cartesian(xlim = c(0,1))
  g
}
