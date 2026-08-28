#!/usr/bin/env Rscript
# Dump IV IF for all 3 instruments (treat2, treat3, treat4).
# Cell fixed: y_bin=1, m_val=0, treated transform.

suppressPackageStartupMessages({
  library(haven); library(dplyr); library(fixest); library(sandwich)
})

setwd("~/Library/CloudStorage/OneDrive-Personal/Brown/26 RA with Jon/stata-testmechs")

df <- read_dta("data/mother_data_extended.dta")
needed <- c("treat", "grandmother", "motherfinancial", "uc",
            "age_baseline", "edu_mo_baseline", "wealth_baseline",
            "treat2", "treat3", "treat4")
df <- df[complete.cases(df[, needed]), ]
cat("Rows:", nrow(df), "\n\n")

df$y_bin <- ntile(df$motherfinancial, 5)

Y_VAL <- 1; M_VAL <- 0
df$lhs_tmp <- df$treat * as.numeric(df$y_bin == Y_VAL & df$grandmother == M_VAL)

for (instr in c("treat2", "treat3", "treat4")) {
  cat("========== IV with", instr, "==========\n")

  fml <- as.formula(sprintf(
    "lhs_tmp ~ age_baseline + edu_mo_baseline + wealth_baseline | 0 | treat ~ %s",
    instr))

  reg <- tryCatch(fixest::feols(fml = fml, data = df),
                  error = function(e) { cat("R errored:", conditionMessage(e), "\n"); NULL })
  if (is.null(reg)) next

  cat("Coef names:", paste(rownames(reg$coeftable), collapse=", "), "\n")
  cat("IV estimate for fit_treat:",
      reg$coeftable[grep("fit_treat|^treat$", rownames(reg$coeftable))[1], "Estimate"], "\n")

  S <- sandwich::estfun(reg)
  B <- sandwich::bread(reg)
  IF <- (S %*% t(B))
  treat_pos <- if ("treat" %in% colnames(IF)) which(colnames(IF) == "treat")
               else grep("fit_treat", colnames(IF))[1]

  if_treat <- as.numeric(IF[, treat_pos])
  cat("IF first 5:", head(if_treat, 5), "\n\n")

  out_df <- data.frame(row = seq_along(if_treat),
                       y_bin = df$y_bin,
                       if_treat = if_treat)
  write.csv(out_df, sprintf("/tmp/r_iv_if_%s.csv", instr), row.names = FALSE)
}

cat("Dumped 3 R IF files to /tmp/r_iv_if_*.csv\n")
