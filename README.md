# medianeffect

## Dose-Effect Analysis with the Median Effect Principle

The aim of the `medianeffect` package is to facilitate single and multiple drug dose-effect analysis
via the Median Effect Principle, including quantitation of synergism and antagonism in drug combinations.

This work is neither approved nor endorsed by Dr. Ting-Chao Chou, the developer of this analytic approach.

## Installation

```r
# Install from GitHub:
devtools::install_github("jhchou/medianeffect", build_vignettes = TRUE)
```

## Example usage

A package vignette provides a detailed example of most functionality and includessome technical background.

```r
# View package vignette
vignette('medianeffect')
```

## Functions

### Object generation

* `drug_effects()` Drug effects object generator, either single drug or constant ratio combination
* `ncr_drug_effects()` Non-constant ratio drug effects object generator

### Calculators

* `calc_ci()` Calculate combination index (CI)
* `calc_dri()` Calculate dose reduction index (DRI)
* `ncr_calc_ci()` Calculate combination index (CI) for non-constant ratio

### Plots

* `median_effect_plot()` Generate Median-Effect Plot
* `dose_effect_plot()` Generate Dose-Effect Plot
* `fa_ci_plot()` Generate fa-CI Plot
* `fa_dri_plot()` Generate fa-DRI Plot
* `isobologram_plot()` Generate isobologram plot

### Methods

* `print.drug_effects()` Print method for objects of drug_effects class
* `print.ncr_drug_effects()` Print method for non-constant ratio drug_effects class

### Internal helper functions

* `calc_d()` Calculate predicted doses
* `calc_fa()` Calculate predicted fraction affected (fa)
* `calc_combo()` Drug combination calculations
* `drug_combo_decompose()` Function to decompose a drug_combo drug effects object into its component drug doses
