#
# R package for Median Effect principle calculations and plots
#
# Requires
# - knitr (for kable)
# - ggplot2
#
#
# Class: drug_effects
#
# Single drug series
# - vector of doses (D)
# - vector of fraction affected (fa)
# - drug name  (name)
# - drug units (units)
# - Calculated: Dm, m, R^2, R, log10(D), log10(fa / (1-fa))
#
#
# Combination drug series, constant ratio
# - ADD: ratio: drug 1 / drug 2 (ratio)
# - ADD: drug 2 name  (name2)
# - ADD: drug 2 units (units2)
#
#
# Combination drug series, non-constant ratio
# - vector of drug 1 doses (D1)
# - vector of drug 2 doses (D2)
# - vector of fraction affected (fa)
# - drug 1 name  (name)
# - drug 1 units (units)
# - drug 2 name  (name2)
# - drug 2 units (units2)
#
#
# Methods:
#
# - print() / export()
# - calculate fa from D
# - calculate D from fa
# - sigmoidal: plot fa versus D
# - linear:    plot log10(fa / (1-fa)) versus log10(D)
#
# - fa-CI plot from drug 1, drug 2, constant / non-constant ratio combo
# - isobologram?
# - print() / export() combo analysis


# Drug effects object constructor: for single drugs and constant ratio combinations
new_drug_effects <- function(D = double(), fa = double(), name = "", units = "", ratio = double(), name2 = "", units2 = "") {
  stopifnot(is.double(D), is.double(fa), is.double(ratio))
  
  values <- list(D = D, fa = fa, name = name, units = units, ratio = ratio, name2 = name2, units2 = units2)
  class <- c("drug_effects")
  
  if (length(ratio) != 0) {
    # values <- c(values, list(ratio = ratio, name2 = name2, units2 = units2))
    class <- append("combo_drug_effects", class)
  }
  
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
drug_effects <- function(D = double(), fa = double(), name = "", units = "", ratio = double(), name2 = "", units2 = "") {
  D <- as.double(D)
  fa = as.double(fa)
  name = as.character(name)
  units = as.character(units)
  ratio = as.double(ratio)
  name2 = as.character(name2)
  units2 = as.character(units2)
  validate_drug_effects(new_drug_effects(D, fa, name, units, ratio, name2, units2))
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


# Given drug 1, drug 2, and drug combination drug_effects objects, calculate combination index CI at vector of fa values
calc_CI <- function(drug1, drug2, drug_combo, fa) {
  ratio <- drug_combo$ratio
  
  Dx1 <- calc_D(drug1, fa)
  Dx2 <- calc_D(drug2, fa)
  
  D_combo <- calc_D(drug_combo, fa)
  D1 <- (ratio / (ratio + 1)) * D_combo
  D2 <- D_combo / (ratio + 1)
  
  CI <- (D1 / Dx1) + (D2 / Dx2)
  return(CI)
}


# Add print for drug_effects display
print.drug_effects <- function(x, ..., stats = TRUE) {
  
  D <- x$D
  fa <- x$fa
  name <- x$name
  units <- x$units
  
  ratio <- x$ratio
  name2 <- x$name2
  units2 <- x$units2
  
  Dm <- x$Dm
  
  # log_D <- log10(D)
  # log_fa_fu <- log10(fa / (1-fa))
  
  if (x$name  != "") { cat("Drug: ", name, "\n", sep = '') }
  if (x$units != "") { cat("Units: ", units, "\n", sep = '')}
  
  if (x$name2  != "") { cat("Drug 2: ", name2, "\n", sep = '') }
  if (x$units2 != "") { cat("Units 2: ", units2, "\n", sep = '')}

  # df <- data.frame(list(D = D, fa = fa, log_D = log_D, log_fa_fu = log_fa_fu))
  df <- data.frame(list(D = D, fa = fa))
  if (length(ratio) != 0) {
    cat("\nRatio of drug 1 / drug 2: ", ratio, sep = '')
    D1 <- (ratio / (ratio+1)) * D
    D2 <- D / (ratio + 1)
    df$D1 <- D1
    df$D2 <- D2
    
    Dm <- paste0(signif(Dm, 4), ' = ', signif((ratio / (ratio+1)) * Dm, 4), ' + ', signif(Dm / (ratio + 1), 4))
  }
  
  print(knitr::kable(df))
  
  if (stats) { cat('\nm: ', x$m, '\nDm: ', Dm, '\nR2: ', x$R2, '\nR: ', x$R, sep = '') }
}


plot.drug_effects <- function(x, y, ..., color = 'blue') {
  title <- "Median-Effect Plot"
  if (x$name != "") { title <- paste0(title, ": ", x$name) }
  if (x$units != "") { title <- paste0(title, " (", x$units, ")") }

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


dose_effect_plot <- function(drug, from = 0.01, to = 0.99, by = 0.01) {
  if (!inherits(drug, "drug_effects")) {
    stop(
      "This requires a drug dose/fa object",
      call. = FALSE
    )
  }
  
  calc_d <- function(fa, m, Dm) { Dm*(fa / (1 - fa))^(1/m) }
  
  m  <- drug$m
  Dm <- drug$Dm
  
  df <- data.frame(list(D = drug$D, fa = drug$fa))
  
  fa <- seq(from = from, to = to, by = by)
  df2 <- data.frame(list(D = calc_d(fa, m, Dm), fa = fa))
  g <- ggplot2::ggplot(data = df2, ggplot2::aes(D, fa)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(data = df, ggplot2::aes(D, fa)) +
    ggplot2::coord_cartesian(ylim = c(0,1))
  
  g
}


fa_ci_plot <- function(drug1, drug2, drug_combo, from = 0.01, to = 0.99, by = 0.01) {
  fa <- seq(from = from, to = to, by = by)
  CI = calc_CI(drug1, drug2, drug_combo, fa)
  
  df <- data.frame(list(fa = fa, CI = CI))
  g <- ggplot2::ggplot(df, ggplot2::aes(fa, CI)) +
    ggplot2::geom_line() +
    ggplot2::geom_hline(yintercept = 1.0, linetype = 'dotted') +
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


D  <- c(0.10, 0.15, 0.20, 0.25, 0.35)
fa <- c(0.24, 0.44, 0.63, 0.81, 0.90)

d <- drug_effects(D = D, fa = fa, name = 'Drug name', units = 'mg')
d

d2 <- drug_effects(D = D, fa = fa, name = 'DrugOne', units = 'mg', ratio = 3, name2 = "DrugTwo")

plot(d)

dose_effect_plot(d)

# str(d)
# class(d)
# inherits(d, "drug_effects")
# unclass(d)


# Taxol
drug1 <- drug_effects(D = c(0.002, 0.004, 0.005, 0.01, 0.02), fa = c(0.429, 0.708, 0.761, 0.882, 0.932), name = "Taxol", units = "uM")

# Cis-Pt
drug2 <- drug_effects(D = c(0.05, 0.1, 0.2, 0.5, 1, 2), fa = c(0.055, 0.233, 0.301, 0.559, 0.821, 0.953), name = "Cis-Pt", units = "uM")

# Taxol:Cis-Pt Combo, ratio = 0.1
# drug_combo <- drug_effects(D = c(0.11, 0.22, 0.55, 1.1), fa = c(0.45, 0.671, 0.921, 0.958), ratio = 0.1)

drug_combo <- drug_effects(D = c(0.001, 0.002, 0.005, 0.01) + c(0.1, 0.2, 0.5, 1), fa = c(0.45, 0.671, 0.921, 0.958), ratio = 0.01)

calc_CI(drug1, drug2, drug_combo, fa = c(0.5, 0.75, 0.9, 0.95))

fa_ci_plot(drug1, drug2, drug_combo, from = 0.1, to = 0.9)
