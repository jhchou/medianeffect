#
# R package for Median Effect principle calculations and plots
#
# Requires
# - knitr (for kable)
# - ggplot2
#
#
# Classes:
#
# - single drug series
# - vector of doses (D)
# - vector of fraction affected (fa)
# - drug name  (name)
# - drug units (units)
# - Dm, m, R^2, R, log10(D), log10(fa / (1-fa))
#
#
# - combination drug series, constant ratio
# - subclass from single drug series?
# - vector of sum of drug 1 and drug 2 doses (D)
# - vector of fraction affected (fa)
# - ratio: drug 1 / drug 2 (ratio)
# - drug 1 name  (name1)
# - drug 1 units (units1)
# - drug 2 name  (name2)
# - drug 2 units (units2)
#
# - combination drug series, non-constant ratio
# - subclass from combination drug series? adds vector of drug 2 doses and removes ratio
# - vector of drug 1 doses (D1)
# - vector of drug 2 doses (D2)
# - vector of fraction affected (fa)
# - drug 1 name  (name1)
# - drug 1 units (units1)
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
#


# Single drug object constructor
# - a bit sloppy, as does not include the calculated values
new_single_drug <- function(D = double(), fa = double(), name = character(0), units = character(0)) {
  stopifnot(is.double(D), is.double(fa))
  structure(list(D = D, fa = fa, name = name, units = units), class = 'single_drug')
}

# Single drug object validator + statistics calculator
validate_single_drug <- function(x) {
  D <- x$D
  fa <- x$fa
  
  if (!all(!is.na(D) & D > 0.0)) {
    stop(
      "All doses `D` must be non-missing and greater than zero",
      call. = FALSE
    )
  }
  
  if (!all(!is.na(fa) & fa > 0.0 & fa < 1.0)) {
    stop(
      "All fraction affected `fa` must be non-missing, greater than zero, and less than one",
      call. = FALSE
    )
  }
  
  if (length(D) != length(fa)) {
    stop(
      "The number of doses must equal the number of fraction affected",
      call. = FALSE
    )
  }
  
  if (length(D) == 0) {
    stop(
      "There must be at least one dose / fraction affected point",
      call. = FALSE
    )
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
single_drug <- function(D = double(), fa = double(), name = "", units = "") {
  D <- as.double(D)
  fa = as.double(fa)
  name = as.character(name)
  units = as.character(units)
  validate_single_drug(new_single_drug(D, fa, name, units))
}



# Add print for single_drug display
print.single_drug <- function(x, ..., stats = TRUE) {
  
  D <- x$D
  fa <- x$fa
  name <- x$name
  units <- x$units
  
  log_D <- log10(D)
  log_fa_fu <- log10(fa / (1-fa))
  
  df <- data.frame(list(D = D, fa = fa, log_D = log_D, log_fa_fu = log_fa_fu))
  
  if (x$name  != "") { cat("Drug: ", name, "\n", sep = '') }
  if (x$units != "") { cat("Units: ", units, "\n", sep = '')}
  
  print(knitr::kable(df))
  
  if (stats) { cat('\nDm: ', x$Dm, '\nm: ', x$m, '\nR2: ', x$R2, '\nR: ', x$R, sep = '') }
}


plot.single_drug <- function(x, y, ..., color = 'blue') {
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
  if (!inherits(drug, "single_drug")) {
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

d <- single_drug(D = D, fa = fa, name = 'Drug name', units = 'mg')
d

plot(d)

# str(d)
# class(d)
# inherits(d, "single_drug")
# unclass(d)




calc_fa <- function(D, m, Dm) {
  1 / (1 + ((Dm / D)^m))
}

calc_d <- function(fa, m, Dm) {
  Dm*(fa / (1 - fa))^(1/m)
}



df %>% ggplot(aes(D, fa)) + geom_point() + theme_bw()

df %>% ggplot(aes(log_D, log_fa_fu)) + geom_point() + stat_smooth(method = "lm", col = "blue") + theme_bw()

fa <- seq(from = 0.01, to = 0.99, by = 0.01)
df2 <- as.data.frame(list(D = calc_d(fa, m, Dm), fa = fa))
df2 %>% ggplot(aes(D, fa)) + geom_line() + theme_bw() + coord_cartesian(ylim = c(0,1))

# df2 %>% mutate(log_D = log10(D), log_fa_fu = log10(fa / (1-fa))) %>% ggplot(aes(log_D, log_fa_fu)) + geom_line()






