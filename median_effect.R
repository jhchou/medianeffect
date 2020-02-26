#
# R package for Median Effect principle calculations and plots
#
# Classes:
#
# - single drug series
# - vector of doses (D)
# - vector of fraction affected (fa)
# - drug name  (name)
# - drug units (units)
# - CONSIDER: Dm, m, R^2, R, log10(D), log10(fa / (1-fa))
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


library(tidyverse)
library(broom)

calc_fa <- function(D, m, Dm) {
  1 / (1 + ((Dm / D)^m))
}

calc_d <- function(fa, m, Dm) {
  Dm*(fa / (1 - fa))^(1/m)
}

df <- tribble(
  ~D, ~fa,
  0.1, 0.24,
  0.15, 0.44,
  0.2, 0.63,
  0.25, 0.81,
  0.35, 0.9
)


df <- tribble(
  ~D, ~fa,
  0.002, 0.429,
  0.004, 0.708,
  0.005, 0.761,
  0.01, 0.882,
  0.02, 0.932
)

df <- tribble(
  ~D, ~fa,
  0.05, 0.055,
  0.1, 0.233,
  0.2, 0.301,
  0.5, 0.559,
  1, 0.821,
  2, 0.953
)




# Single drug object constructor
new_single_drug <- function(D = double(), fa = double(), name = character(), units = character()) {
  stopifnot(is.double(D), is.double(fa))
  structure(list(D = D, fa = fa, name = name, units = units), class = 'single_drug')
}

# Single drug object validator
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
  
  x # return the input
}


# Single drug object helper
single_drug <- function(D = double(), fa = double(), name = character(), units = character()) {
  D <- as.double(D)
  fa = as.double(fa)
  name = as.character(name)
  units = as.character(units)
  validate_single_drug(new_single_drug(D, fa, name, units))
}

df <- tribble(
  ~D, ~fa,
  0.1, 0.24,
  0.15, 0.44,
  0.2, 0.63,
  0.25, 0.81,
  0.35, 0.9
)


# Add print for single_drug display
print.single_drug <- function(x, ..., verbose = FALSE, stats = TRUE) {
  df <- data.frame(list(D = x$D, fa = x$fa))
  if (verbose) {
    df$logD <- log10(x$D)
    df$log_fa_fu <- log10(df$fa / (1-df$fa))
  }
  print(knitr::kable(df))
  cat("\nboog\n\n\n")
}


d <- single_drug(D =df$D, fa = df$fa)
class(d)
inherits(d, "single_drug")



df %>% ggplot(aes(D, fa)) + geom_point() + theme_bw()

df <- df %>%
  mutate(
    log_D = log10(D),
    log_fa_fu = log10(fa / (1-fa))
  )

df %>% ggplot(aes(log_D, log_fa_fu)) + geom_point() + stat_smooth(method = "lm", col = "blue") + theme_bw()

fit <- lm(log_fa_fu ~ log_D, df)
summary(fit)

summary(fit)$r.squared
sqrt(summary(fit)$r.squared)

fit$coefficients

g <- glance(fit)
g

g$r.squared %>% sqrt

t <- tidy(fit)
t

a <- augment(fit)

m <- t$estimate[2]
b <- t$estimate[1]

Dm <- 10^(-b/m)

m
Dm

fa <- seq(from = 0.01, to = 0.99, by = 0.01)
df2 <- as.data.frame(list(D = calc_d(fa, m, Dm), fa = fa))
df2 %>% ggplot(aes(D, fa)) + geom_line() + theme_bw() + coord_cartesian(ylim = c(0,1))

# df2 %>% mutate(log_D = log10(D), log_fa_fu = log10(fa / (1-fa))) %>% ggplot(aes(log_D, log_fa_fu)) + geom_line()
