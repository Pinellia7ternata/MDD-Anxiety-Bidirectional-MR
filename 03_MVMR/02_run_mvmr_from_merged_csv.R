#!/usr/bin/env Rscript
# ============================================================
# MVMR from a pre-merged harmonised CSV matrix
# Use this when OpenGWAS download fails or when all covariate GWAS files are local.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(MVMR)
})

args <- commandArgs(trailingOnly = TRUE)
INFILE <- ifelse(length(args) >= 1, args[1], "D:/2026年/文章/因果推断/06_MVMR/MVMR_merged_harmonised.csv")
OUTDIR <- ifelse(length(args) >= 2, args[2], dirname(INFILE))
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# Direction can be MDD_to_Anxiety or Anxiety_to_MDD. If not supplied, both are attempted.
DIRECTION <- ifelse(length(args) >= 3, args[3], "both")

safe_num <- function(x) as.numeric(as.character(x))

run_direction <- function(dt, direction) {
  if (direction == "MDD_to_Anxiety") {
    exposure_names <- c("MDD", "BMI", "Education", "Smoking")
    betaY <- "beta_Anxiety"; seY <- "se_Anxiety"
  } else if (direction == "Anxiety_to_MDD") {
    exposure_names <- c("Anxiety", "BMI", "Education", "Smoking")
    betaY <- "beta_MDD"; seY <- "se_MDD"
  } else {
    stop("Unknown direction: ", direction)
  }
  needed <- c("SNP", betaY, seY, paste0("beta_", exposure_names), paste0("se_", exposure_names))
  missing <- setdiff(needed, names(dt))
  if (length(missing) > 0) stop("Missing columns for ", direction, ": ", paste(missing, collapse = ", "))
  d <- copy(dt)
  for (cc in needed[needed != "SNP"]) d[[cc]] <- safe_num(d[[cc]])
  d <- d[complete.cases(d[, ..needed])]
  if (nrow(d) < length(exposure_names) + 5) stop("Too few complete SNPs for ", direction, ": ", nrow(d))

  BX <- as.matrix(d[, paste0("beta_", exposure_names), with = FALSE])
  seBX <- as.matrix(d[, paste0("se_", exposure_names), with = FALSE])
  BY <- d[[betaY]]
  seBY <- d[[seY]]

  r_input <- format_mvmr(BXGs = BX, BYG = BY, seBXGs = seBX, seBYG = seBY, RSID = d$SNP)
  ivw <- ivw_mvmr(r_input, model = "random")
  strength <- tryCatch(strength_mvmr(r_input, gencov = 0), error = function(e) e)
  pleio <- tryCatch(pleiotropy_mvmr(r_input, gencov = 0), error = function(e) e)

  # Manual weighted regression backup - useful if MVMR package output is hard to parse.
  W <- diag(1 / (seBY^2))
  beta_hat <- solve(t(BX) %*% W %*% BX, t(BX) %*% W %*% BY)
  residual <- BY - BX %*% beta_hat
  Q <- as.numeric(t(residual) %*% W %*% residual)
  df <- nrow(BX) - ncol(BX)
  phi <- max(1, Q / df)
  vcov_beta <- phi * solve(t(BX) %*% W %*% BX)
  se_hat <- sqrt(diag(vcov_beta))
  z <- as.numeric(beta_hat) / se_hat
  p <- 2 * pnorm(abs(z), lower.tail = FALSE)
  out <- data.table(
    direction = direction,
    exposure = exposure_names,
    beta = as.numeric(beta_hat),
    se = se_hat,
    OR = exp(as.numeric(beta_hat)),
    CI_lower = exp(as.numeric(beta_hat) - 1.96 * se_hat),
    CI_upper = exp(as.numeric(beta_hat) + 1.96 * se_hat),
    p_value = p,
    n_snps = nrow(d),
    Q = Q,
    Q_df = df,
    Q_p = pchisq(Q, df = df, lower.tail = FALSE)
  )

  fwrite(out, file.path(OUTDIR, paste0(direction, "_MVMR_manual_IVW_random.csv")))
  sink(file.path(OUTDIR, paste0(direction, "_MVMR_package_output.txt")))
  cat("Direction:", direction, "\n")
  cat("Complete SNPs:", nrow(d), "\n\n")
  cat("--- MVMR package IVW ---\n"); print(ivw)
  cat("\n--- Conditional F / strength ---\n"); print(strength)
  cat("\n--- Pleiotropy / heterogeneity ---\n"); print(pleio)
  cat("\n--- Manual weighted regression backup ---\n"); print(out)
  sink()
  return(out)
}

dt <- fread(INFILE)
results <- list()
if (DIRECTION %in% c("both", "MDD_to_Anxiety")) results[["MDD_to_Anxiety"]] <- run_direction(dt, "MDD_to_Anxiety")
if (DIRECTION %in% c("both", "Anxiety_to_MDD")) results[["Anxiety_to_MDD"]] <- run_direction(dt, "Anxiety_to_MDD")
allres <- rbindlist(results, fill = TRUE)
fwrite(allres, file.path(OUTDIR, "MVMR_all_directions_summary.csv"))
message("Done. Summary: ", file.path(OUTDIR, "MVMR_all_directions_summary.csv"))
