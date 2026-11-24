# Stage 4: full multivariable IVW coefficient matrices for both instrument designs,
# plus the sensitivity run that supplies the cross-trait covariance (gencov) for the
# primary-exposure pair. Written because the main text now reports only the direct
# estimate and the supplementary table has to carry every coefficient.
#
# Output: mvmr_independent/mvmr_full_models.csv
#   model | exposure | nsnp | beta | se | or | lo | hi | p | conditional_F

.libPaths(c("C:/Rlibs", "C:/Users/\u9648\u8d85/Rlibs", .libPaths()))
suppressPackageStartupMessages(library(MVMR))

OUT <- "08_revision_analyses/mvmr/mvmr_full_models.csv"
d <- read.csv("08_revision_analyses/mvmr/mvmr_matrix_harmonised.csv",
              check.names = FALSE)

COV  <- c("beta_BMI", "beta_Education", "beta_SmokingInit")
CSE  <- c("se_BMI", "se_Education", "se_SmokingInit")

fit <- function(rows, expc, ob, os, label, primary_name) {
  NAMES <- c(primary_name, "BMI", "Education", "Smoking initiation")
  expse <- paste0("se_", sub("beta_", "", expc))
  BX <- as.matrix(rows[, c(expc, COV)]); seBX <- as.matrix(rows[, c(expse, CSE)])
  r  <- format_mvmr(BXGs = BX, BYG = rows[[ob]], seBXGs = seBX, seBYG = rows[[os]], RSID = rows$snp)
  iv <- suppressWarnings(ivw_mvmr(r))
  st <- suppressWarnings(strength_mvmr(r))
  cf <- as.numeric(as.vector(unlist(st)))
  data.frame(model = label, exposure = NAMES, nsnp = nrow(rows),
             beta = as.numeric(iv[, "Estimate"]), se = as.numeric(iv[, "Std. Error"]),
             or_ = exp(as.numeric(iv[, "Estimate"])),
             lo = exp(as.numeric(iv[, "Estimate"]) - 1.96 * as.numeric(iv[, "Std. Error"])),
             hi = exp(as.numeric(iv[, "Estimate"]) + 1.96 * as.numeric(iv[, "Std. Error"])),
             p = as.numeric(iv[, "Pr(>|t|)"]),
             conditional_F = cf[1:length(NAMES)],
             stringsAsFactors = FALSE)
}

res <- list()
for (dir in c("forward", "reverse")) {
  if (dir == "forward") {
    expc <- "beta_MDD";     ob <- "beta_Anxiety"; os <- "se_Anxiety"; flag <- "is_MDD_iv"
    lab_exp <- "MDD";       lab_dir <- "MDD -> anxiety"
  } else {
    expc <- "beta_Anxiety"; ob <- "beta_MDD";     os <- "se_MDD";     flag <- "is_Anx_iv"
    lab_exp <- "Anxiety";   lab_dir <- "Anxiety -> MDD"
  }
  primary  <- d[d[[flag]] == 1, ]
  exp_only <- primary[primary$indep_BMI == 0 & primary$indep_Education == 0 & primary$indep_SmokingInit == 0, ]
  cat("\n=====", lab_dir, "=====\n")
  cat("  primary-harmonised rows:", nrow(primary), " of which complete with no independent covariate hit:",
      nrow(exp_only), " union set:", nrow(d), "\n")
  res[[length(res) + 1]] <- fit(exp_only, expc, ob, os, paste0(lab_dir, " | covariates instrumented by the exposure"), lab_exp)
  res[[length(res) + 1]] <- fit(d,      expc, ob, os, paste0(lab_dir, " | covariate instruments selected independently"), lab_exp)
}
out <- do.call(rbind, res)
write.csv(out, OUT, row.names = FALSE)
print(out[out$exposure %in% c("MDD", "Anxiety"), c("model", "nsnp", "beta", "se", "or_", "p", "conditional_F")], digits = 4)
cat("\ncovariate coefficients (independent design only):\n")
print(out[grepl("independently", out$model) & !out$exposure %in% c("MDD", "Anxiety"),
          c("model", "exposure", "beta", "se", "or_", "p", "conditional_F")], digits = 4)
cat("\n->", OUT, "\n")
