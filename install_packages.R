# Install all R packages required by the final Silwood project.
# Run once from the project root:
#   Rscript install_packages.R

required_packages <- c(
  "tidyverse",
  "here",
  "lme4",
  "Matrix",
  "MASS",
  "patchwork",
  "lubridate",
  "scales",
  "ragg",
  "digest"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]

if (length(missing_packages) == 0L) {
  message("All required packages are already installed.")
} else {
  install.packages(
    missing_packages,
    dependencies = TRUE,
    repos = c(CRAN = "https://cloud.r-project.org")
  )
  message("Installed: ", paste(missing_packages, collapse = ", "))
}
