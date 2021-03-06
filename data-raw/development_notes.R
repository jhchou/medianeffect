devtools::document()
devtools::build_vignettes()

# Installation with vignettes:
devtools::install(build_vignettes = TRUE)
devtools::install_github("jhchou/medianeffect", build_vignettes = TRUE)

library(tidyverse)
library(medianeffect)
vignette('medianeffect')
browseVignettes('medianeffect')
help(package = medianeffect)


# Package vignettes:
# - http://r-pkgs.had.co.nz/vignettes.html#vignette-workflow-2
# - usethis::use_vignette("medianeffect")
# - not sure if need to do this: in RStudio configure build tools, check the use roxygen to generate vignettes box -- unchecked still works?
# - https://stackoverflow.com/questions/33614660/knitr-rmd-vignettes-do-not-appear-with-vignette
#
# Note devtools does not build vignettes by default when you devtools::install() (same thing for some install_* functions like install_github()) a package from a directory. You have to specify the argument build_vignettes = TRUE when you install the package. Currently there is no way to build vignettes using devtools if you just use the RStudio button Build & Reload. You have to Build Source Package, and run R CMD INSTALL on the tarball. Or run devtools::install(build_vignettes = TRUE) in the R console.


# To do:
# [ ] write readme
# [ ] make vignette available online

# Consider:
# [ ] isobologram
# [ ] load from Excel
# [ ] split code to multiple files
# [ ] export raw data calculations as CSV

