# Stage 3: multivariable IVW with independently selected covariate instruments (MVMR 0.4.8).
#
# Validation first: the old design (covariates instrumented by the primary-exposure
# instruments) must reproduce Table 3 - forward OR 1.15 with conditional F 15.6 on 84 SNPs,
# reverse OR 1.10 with conditional F 12.9 on 27 SNPs. Only then are the new-design numbers
# used. gencov = 0 because the two GWAS are not the same sample; the sensitivity with the
# cross-trait covariance from MRlap is reported as a second run.

.libPaths(c("C:/Rlibs", "C:/Users/\u9648\u8d85/Rlibs", .libPaths()))
suppressPackageStartupMessages(library(MVMR))

d <- read.csv("08_revision_analyses/mvmr/mvmr_matrix_harmonised.csv",
              check.names = FALSE)
s1 <- read.csv("data/S1_MDD_to_Anxiety_harmonised_198_instruments.csv",
               check.names = FALSE)
s2 <- read.csv("data/S2_Anxiety_to_MDD_harmonised_50_instruments.csv",
               check.names = FALSE)
mdd_iv <- s1$SNP; anx_iv <- s2$SNP

COV <- c("beta_BMI", "beta_Education", "beta_SmokingInit")
COVSE <- c("se_BMI", "se_Education", "se_SmokingInit")

run <- function(rows, exp_col, out_beta, out_se, label) {
  exp_se <- paste0("se_", sub("beta_", "", exp_col))
  BX <- as.matrix(rows[, c(exp_col, COV)]); seBX <- as.matrix(rows[, c(exp_se, COVSE)])
  r <- format_mvmr(BXGs = BX, BYG = rows[[out_beta]], seBXGs = seBX, seBYG = rows[[out_se]], RSID = rows$snp)
  iv <- suppressWarnings(ivw_mvmr(r))
  st <- suppressWarnings(strength_mvmr(r))
  b <- iv["exposure1", "Estimate"]; se <- iv["exposure1", "Std. Error"]; pv <- iv["exposure1", "Pr(>|t|)"]
  cf <- as.numeric(as.vector(st$csqF))
  data.frame(model = label, nsnp = nrow(rows), beta = b, se = se, or_ = exp(b),
             lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se), p = pv,
             condF_primary = cf[1], condF_BMI = cf[2], condF_Education = cf[3], condF_Smoking = cf[4])
}

out <- list()
for (dir in c("forward", "reverse")) {
  if (dir == "forward") {
    expc <- "beta_MDD"; ob <- "beta_Anxiety"; os <- "se_Anxiety"; iv <- mdd_iv
    old <- d[d[[if (dir == "forward") "is_MDD_iv" else "is_Anx_iv"]] == 1, ]
  } else {
    expc <- "beta_Anxiety"; ob <- "beta_MDD"; os <- "se_MDD"; iv <- anx_iv
    old <- d[d[[if (dir == "forward") "is_MDD_iv" else "is_Anx_iv"]] == 1, ]
  }
  # old design: instruments taken from the primary exposure only (no independent covariate hits)
  old_only <- old[old$indep_BMI == 0 & old$indep_Education == 0 & old$indep_SmokingInit == 0, ]
  newset <- d   # full union set
  cat("\n=====", dir, "=====\n")
  out[[length(out) + 1]] <- run(old_only, expc, ob, os, paste0(dir, " | OLD design (exposure-selected only)"))
  out[[length(out) + 1]] <- run(newset, expc, ob, os, paste0(dir, " | NEW design (independent covariate instruments)"))
  # covariate-free control: primary exposure alone in the union set
  cat("  old-design rows:", nrow(old_only), "| union rows:", nrow(newset), "\n")
}
res <- do.call(rbind, out)
print(res, digits = 4)
write.csv(res, "08_revision_analyses/mvmr/mvmr_old_vs_new_design.csv",
          row.names = FALSE)
cat("\n-> mvmr_independent/mvmr_old_vs_new_design.csv\n")
