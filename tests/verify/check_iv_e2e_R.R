#!/usr/bin/env Rscript
# End-to-end test_sharp_null with IV reg_formula: get R pvals for
# 3 stress cases to compare with Stata.

suppressPackageStartupMessages({
  library(haven); library(dplyr); library(fixest); library(osqp); library(Matrix)
})

setwd("~/Library/CloudStorage/OneDrive-Personal/Brown/26 RA with Jon/stata-testmechs")

df <- read_dta("data/mother_data_extended.dta")
needed <- c("treat", "grandmother", "relationship_husb", "motherfinancial", "uc",
            "age_baseline", "edu_mo_baseline", "wealth_baseline", "treat3")
df <- df[complete.cases(df[, needed]), ]
cat("Rows:", nrow(df), "\n\n")

# Load R TestMechs
r_dir <- "r_reference/TestMechs-master/R"
r_files <- list.files(r_dir, pattern = "\\.R$", full.names = TRUE)
for (f in r_files) {
  src <- readLines(f)
  src <- gsub("verbose = TRUE", "verbose = FALSE", src, fixed = TRUE)
  tmp <- tempfile(fileext = ".R")
  writeLines(src, tmp)
  source(tmp, local = FALSE)
}

reg_formula_iv <- "~ age_baseline + edu_mo_baseline + wealth_baseline | 0 | treat ~ treat3"

for (m_var in list(c("grandmother"),
                   c("relationship_husb"),
                   c("relationship_husb", "grandmother"))) {
  cat("========== IV treat3 for m =", paste(m_var, collapse=","), "==========\n")
  result <- tryCatch(
    test_sharp_null(
      df = df, d = "treat", m = m_var, y = "motherfinancial",
      method = "CS", num_Ybins = 5, cluster = "uc",
      reg_formula = reg_formula_iv),
    error = function(e) { cat("R errored:", conditionMessage(e), "\n"); NULL })

  if (!is.null(result)) {
    cat(sprintf("  T_CC  = %.7f\n", result$T_CC))
    cat(sprintf("  cv_CC = %.7f\n", result$cv_CC))
    cat(sprintf("  pval  = %.7f\n", result$pval))
    cat(sprintf("  df    = %d\n\n", result$df))
  }
}
