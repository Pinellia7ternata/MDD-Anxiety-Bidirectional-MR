# Stage 8b: conditional F-statistics for the LD-overlap sensitivity.
#
# Position clumping at 1 Mb leaves covariate instruments that sit close to a primary instrument: 14 of 154
# within 1 Mb and 4 within 100 kb. If BMI's conditional F of 16.6 is partly borrowed strength from instruments
# selected for the disorder, removing those pairs should move it. This refits the two exposure-only models on
# the reduced sets and reports the coefficients and the conditional F for all four exposures.
#
# Output: mvmr_independent/mvmr_ld_sensitivity.csv

.libPaths(c("C:/Rlibs", "C:/Users/\u9648\u8d85/Rlibs", .libPaths()))
suppressPackageStartupMessages(library(MVMR))

d <- read.csv("08_revision_analyses/mvmr/mvmr_matrix_harmonised.csv",
              check.names = FALSE)
COV <- c("beta_BMI", "beta_Education", "beta_SmokingInit")
CSE <- c("se_BMI", "se_Education", "se_SmokingInit")

prim <- d[d$is_MDD_iv == 1 | d$is_Anx_iv == 1, c("snp", "chrom", "pos")]
covr <- d[d$indep_BMI == 1 | d$indep_Education == 1 | d$indep_SmokingInit == 1, c("snp", "chrom", "pos")]
mind <- sapply(seq_len(nrow(covr)), function(i) {
  p <- prim$pos[prim$chrom == covr$chrom[i]]
  if (length(p) == 0) return(Inf)
  min(abs(p - covr$pos[i]))
})
near1  <- covr$snp[mind < 1e6]
near01 <- covr$snp[mind < 1e5]
cat("covariate instruments within 1 Mb of a primary instrument:", length(near1),
    "\n  within 100 kb:", length(near01), "\n")

fit <- function(rows, xb, yb, ys, label) {
  xs <- paste0("se_", sub("beta_", "", xb))
  BX <- as.matrix(rows[, c(xb, COV)]); seBX <- as.matrix(rows[, c(xs, CSE)])
  r  <- format_mvmr(BXGs = BX, BYG = rows[[yb]], seBXGs = seBX, seBYG = rows[[ys]], RSID = rows$snp)
  iv <- suppressWarnings(ivw_mvmr(r))
  st <- suppressWarnings(strength_mvmr(r))
  cf <- as.numeric(as.vector(unlist(st)))
  data.frame(model = label, nsnp = nrow(rows),
             beta_pkg = as.numeric(iv["exposure1", "Estimate"]),
             se_pkg = as.numeric(iv["exposure1", "Std. Error"]),
             or_pkg = exp(as.numeric(iv["exposure1", "Estimate"])),
             condF_primary = cf[1], condF_BMI = cf[2], condF_Education = cf[3], condF_Smoking = cf[4],
             stringsAsFactors = FALSE)
}

out <- list()
for (dir in c("forward", "reverse")) {
  if (dir == "forward") {
    xb <- "beta_MDD"; yb <- "beta_Anxiety"; ys <- "se_Anxiety"; keep <- d$is_Anx_iv == 0
  } else {
    xb <- "beta_Anxiety"; yb <- "beta_MDD"; ys <- "se_MDD"; keep <- d$is_MDD_iv == 0
  }
  for (lab in c("as reported", "minus pairs within 100 kb", "minus all within 1 Mb")) {
    drop <- if (lab == "as reported") character(0) else if (grepl("100 kb", lab)) near01 else near1
    rows <- d[keep & !(d$snp %in% drop), ]
    out[[length(out) + 1]] <- fit(rows, xb, yb, ys, paste0(dir, " | ", lab))
  }
}
res <- do.call(rbind, out)
write.csv(res, "08_revision_analyses/mvmr/mvmr_ld_sensitivity.csv",
          row.names = FALSE)
print(res, digits = 4)
cat("\n-> mvmr_independent/mvmr_ld_sensitivity.csv\n")
