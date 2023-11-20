install.packages("covr")
# Install BIOS625 from GitHub
devtools::install_github("Jiayuny54/BIOS625")

library(covr)
library(BIOS625)

# If run with no arguments implicitly calls `package_coverage()`
covr::report()

covr::package_coverage()
