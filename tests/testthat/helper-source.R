# testthat runs from tests/testthat, so the project code is two levels up
root <- normalizePath(file.path(testthat::test_path(), "..", ".."))

suppressPackageStartupMessages({
  library(dplyr); library(admiral); library(admiralonco); library(survival)
})
for (f in c("simulate", "write_sdtm", "derive_adsl", "derive_adrs", "derive_adtte", "methods"))
  source(file.path(root, "trial", "R", paste0(f, ".R")))
