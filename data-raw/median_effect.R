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


# Drug effects object constructor: for single drugs and constant ratio combinations
new_drug_effects <- function(D = double(), fa = double(), name = "", label = "", info = "", ratio = double()) {
  stopifnot(is.double(D), is.double(fa), is.double(ratio))
  values <- list(D = D, fa = fa, name = name, label = label, info = info, ratio = ratio)
  class <- c("drug_effects")
  if (length(ratio) != 0) { class <- append("combo_drug_effects", class) }
  structure(values, class = class)
}

# Drug effects object validator + statistics calculator
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


# Single drug object helper
drug_effects <- function(D = double(), fa = double(), name = "", label = "", info = "", ratio = double()) {
  D <- as.double(D)
  fa = as.double(fa)
  name = as.character(name)
  label = as.character(label)
  info = as.character(info)
  ratio = as.double(ratio)
  validate_drug_effects(new_drug_effects(D, fa, name, label, info, ratio))
}


# Given a drug_effects object and vector of doses, return predicted fraction affected
calc_fa <- function(drug, D) {
  if (!inherits(drug, "drug_effects")) { stop("Requires a `drug_effects` object", call. = FALSE) }
  m <- drug$m
  Dm <- drug$Dm
  return (1 / (1 + ((Dm / D)^m)))
}


# Given a drug_effects object and vector of fa, return predicted doses
calc_D <- function(drug, fa) {
  if (!inherits(drug, "drug_effects")) { stop("Requires a `drug_effects` object", call. = FALSE) }
  m <- drug$m
  Dm <- drug$Dm
  return (Dm*(fa / (1 - fa))^(1/m))
}


# Given drug 1, drug 2, and drug_combo objects, calculate combination index CI at vector of fa values
# - if fa is empty, then use the actual fa and actual doses in drug_combo
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



# Add print for drug_effects display
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


plot.drug_effects <- function(x, y, ..., color = 'blue') {
  title <- "Median-Effect Plot"
  if (x$name != "") { title <- paste0(title, ": ", x$name) }
  if (x$label != "") { title <- paste0(title, " (", x$label, ")") }

  df <- data.frame(list(D = x$D, fa = x$fa, log_D = x$log_D, log_fa_fu = x$log_fa_fu))

  g <- ggplot2::ggplot(
    data = df,
    ggplot2::aes(log_D, log_fa_fu)
  ) +
    ggplot2::geom_point() +
    ggplot2::stat_smooth(method = "lm", color = color) +
    ggplot2::xlab("log (D)") +
    ggplot2::ylab("log (fa / fu)") +
    ggplot2::ggtitle(title)

  g
}

median_effect_plot <- function(...) {
  # takes an arbitrary number of drug_effect objects
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
  ggplot2::ggplot(data = df, ggplot2::aes(log_D, log_fa_fu, shape = label, color = label)) +
    ggplot2::geom_point() +
    ggplot2::stat_smooth(method = "lm", se = FALSE, fullrange=TRUE) # fullrange allows extrapolation of line beyond the data points
}

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
  df_lines <- dplyr::mutate(df_lines, D = Dm*(fa / (1 - fa))^(1/m))

  g <- ggplot2::ggplot(df_lines, ggplot2::aes(D, fa, color = label)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(data = df_points, ggplot2::aes(D, fa, color = label, shape = label))
  g
}

#   m  <- drug$m
#   Dm <- drug$Dm
#
#   fa <- seq(from = from, to = to, by = by)
#   df2 <- data.frame(list(D = calc_d(fa, m, Dm), fa = fa))
#   g <- ggplot2::ggplot(data = df2, ggplot2::aes(D, fa)) +
#     ggplot2::geom_line() +
#     ggplot2::geom_point(data = df, ggplot2::aes(D, fa)) +
#     ggplot2::coord_cartesian(ylim = c(0,1))
#
#   g
# }


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

# df <- dplyr::tribble(
#   ~D, ~fa,
#   0.1, 0.24,
#   0.15, 0.44,
#   0.2, 0.63,
#   0.25, 0.81,
#   0.35, 0.9
# )
#
#
# df <- dplyr::tribble(
#   ~D, ~fa,
#   0.002, 0.429,
#   0.004, 0.708,
#   0.005, 0.761,
#   0.01, 0.882,
#   0.02, 0.932
# )
#
# df <- dplyr::tribble(
#   ~D, ~fa,
#   0.05, 0.055,
#   0.1, 0.233,
#   0.2, 0.301,
#   0.5, 0.559,
#   1, 0.821,
#   2, 0.953
# )


# Standard test example
D  <- c(0.10, 0.15, 0.20, 0.25, 0.35)
fa <- c(0.24, 0.44, 0.63, 0.81, 0.90)
d <- drug_effects(D = D, fa = fa)
d # print(d)
plot(d)

# str(d)
# class(d)
# inherits(d, "drug_effects")
# unclass(d)


# Taxol
drug1 <- drug_effects(D = c(0.002, 0.004, 0.005, 0.01, 0.02), fa = c(0.429, 0.708, 0.761, 0.882, 0.932), name = "Taxol", label = "taxol")

# Cis-Pt
drug2 <- drug_effects(D = c(0.05, 0.1, 0.2, 0.5, 1, 2), fa = c(0.055, 0.233, 0.301, 0.559, 0.821, 0.953), name = "Cisplatin", label = "cis-pt")

# Taxol:Cis-Pt Combo, ratio = 0.1
# MISTAKE: drug_combo <- drug_effects(D = c(0.11, 0.22, 0.55, 1.1), fa = c(0.45, 0.671, 0.921, 0.958), ratio = 0.1)
drug_combo <- drug_effects(
  D = c(0.001, 0.002, 0.005, 0.01) + c(0.1, 0.2, 0.5, 1),
  fa = c(0.45, 0.671, 0.921, 0.958),
  name = 'Taxol - Cisplatin',
  label = 'combo',
  ratio = 0.01)

calc_CI(drug1, drug2, drug_combo, fa = c(0.5, 0.75, 0.9, 0.95))

fa_ci_plot(drug1, drug2, drug_combo, from = 0.05, to = 0.95)

median_effect_plot(drug1, drug2, drug_combo)

dose_effect_plot(drug1, drug2, drug_combo, to = 0.98)
