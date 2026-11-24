# Stage 5: multivariable IVW on instrument sets selected on the EXPOSURES only.
#
# Why: the union set used in stage 4 contains the 27 variants that were selected as genome-wide significant
# for anxiety disorders and the 84 selected for MDD, and each of those appears as an instrument in the model
# where the trait that selected it is the OUTCOME. Their outcome-side betas carry winner's curse from the
# outcome GWAS, which pulls the multivariable slope up. The stage-4 forward direct estimate of 2.82 against a
# univariable 2.34 is mostly that, not conditioning: the univariable estimate on the same 264 variants is
# already 2.79. Here the outcome-selected variants are removed from each direction's instrument set, and the
# univariable estimate is computed on exactly the same set so the two are comparable.
#
# Output: mvmr_independent/mvmr_exposure_only.csv
#   set | nsnp | beta | se | or | lo | hi | p | scale | condF_1 | condF_BMI | condF_Edu | condF_Smoke

.libPaths(c("C:/Rlibs", "C:/Users/\u9648\u8d85/Rlibs", .libPaths()))
suppressPackageStartupMessages(library(MVMR))

d <- read.csv("08_revision_analyses/mvmr/mvmr_matrix_harmonised.csv",
              check.names = FALSE)
COV <- c("beta_BMI", "beta_Education", "beta_SmokingInit")
CSE <- c("se_BMI", "se_Education", "se_SmokingInit")

fit <- function(rows, xb, yb, ys, label) {
  xs <- paste0("se_", sub("beta_", "", xb))
  BX <- as.matrix(rows[, c(xb, COV)]); seBX <- as.matrix(rows[, c(xs, CSE)])
  r  <- format_mvmr(BXGs = BX, BYG = rows[[yb]], seBXGs = seBX, seBYG = rows[[ys]], RSID = rows$snp)
  iv <- suppressWarnings(ivw_mvmr(r))
  st <- suppressWarnings(strength_mvmr(r))
  cf <- as.numeric(as.vector(unlist(st)))
  b <- as.numeric(iv["exposure1", "Estimate"]); se <- as.numeric(iv["exposure1", "Std. Error"])
  # ivw_mvmr prints the lm summary but returns only the coefficient matrix, so the residual dispersion that
  # scales its standard errors is recomputed here from the fit
  Xd <- as.matrix(rows[, c(xb, COV)])
  wY <- 1 / rows[[ys]] ^ 2
  resid <- rows[[yb]] - Xd %*% as.numeric(iv[, "Estimate"])
  sig <- sqrt(sum(wY * resid ^ 2) / (nrow(rows) - ncol(Xd)))
  data.frame(set = label, nsnp = nrow(rows), beta = b, se = se, or_ = exp(b),
             lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se),
             p = as.numeric(iv["exposure1", "Pr(>|t|)"]),
             scale = sig,
             condF_primary = cf[1], condF_BMI = cf[2], condF_Education = cf[3], condF_Smoking = cf[4],
             stringsAsFactors = FALSE)
}

out <- list()
# forward: exposures are MDD, BMI, education, smoking; anxiety disorders is the outcome
fwd_exp  <- d[d$is_Anx_iv == 0, ]                                     # drop anxiety-selected instruments
fwd_des  <- d[d$is_MDD_iv == 1 & d$indep_BMI == 0 & d$indep_Education == 0 & d$indep_SmokingInit == 0, ]
# reverse: exposures are anxiety disorders, BMI, education, smoking; MDD is the outcome
rev_exp  <- d[d$is_MDD_iv == 0, ]
rev_des  <- d[d$is_Anx_iv == 1 & d$indep_BMI == 0 & d$indep_Education == 0 & d$indep_SmokingInit == 0, ]

out[[1]] <- fit(fwd_exp, "beta_MDD",     "beta_Anxiety", "se_Anxiety", "forward | exposure-only union")
out[[2]] <- fit(fwd_des, "beta_MDD",     "beta_Anxiety", "se_Anxiety", "forward | covariates from exposure instruments")
out[[3]] <- fit(rev_exp, "beta_Anxiety", "beta_MDD",     "se_MDD",     "reverse | exposure-only union")
out[[4]] <- fit(rev_des, "beta_Anxiety", "beta_MDD",     "se_MDD",     "reverse | covariates from exposure instruments")

res <- do.call(rbind, out)
write.csv(res, "08_revision_analyses/mvmr/mvmr_exposure_only.csv",
          row.names = FALSE)
print(res, digits = 4)
cat("\n-> mvmr_independent/mvmr_exposure_only.csv\n")
