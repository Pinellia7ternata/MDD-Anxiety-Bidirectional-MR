#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
})

BASE <- "D:/2026年/文章/因果推断"
OUT  <- file.path(BASE, "04_univariable_MR")
setwd(BASE)

cat("==================================================\n")
cat("MR Full Results Export (Direct MR, no re-harmonise)\n")
cat("==================================================\n\n")

# ============================================================
# DIRECTION 1: MDD -> Anxiety
# ============================================================
cat("=== MDD -> Anxiety (318 SNPs) ===\n")

h1 <- fread("03_harmonised/MDD_to_Anxiety_harmonised.csv")
d1 <- h1[h1$remove == FALSE]
cat(sprintf("  SNPs: %d\n", nrow(d1)))

# Build dat_mr format directly from already-harmonised data
# exposures/outcomes are already aligned so we skip harmonise_data
dat1 <- data.frame(
  SNP                    = d1$SNP,
  exposure               = "MDD",
  id.exposure            = "MDD",
  beta.exposure          = d1$beta.exposure,
  se.exposure            = d1$se.exposure,
  pval.exposure         = d1$pval.exposure,
  effect_allele.exposure = d1$effect_allele.exposure,
  other_allele.exposure  = d1$other_allele.exposure,
  eaf.exposure          = d1$eaf.exposure,
  samplesize.exposure    = d1$samplesize.exposure,
  outcome               = "Anxiety",
  id.outcome            = "Anxiety",
  beta.outcome          = d1$beta.outcome,
  se.outcome            = d1$se.outcome,
  pval.outcome          = d1$pval.outcome,
  effect_allele.outcome = d1$effect_allele.outcome,
  other_allele.outcome  = d1$other_allele.outcome,
  eaf.outcome          = d1$eaf.outcome,
  samplesize.outcome    = d1$samplesize.outcome,
  mr_keep               = TRUE,
  data.type.exposure    = "easiro",
  data.type.outcome     = "easiro"
)
class(dat1) <- c("data.frame", "dat.trich")

# All MR methods
r1 <- mr(dat1)
r1$OR     <- exp(r1$b)
r1$OR_LCI <- exp(r1$b - 1.96 * r1$se)
r1$OR_UCI <- exp(r1$b + 1.96 * r1$se)
fwrite(r1, file.path(OUT, "MDD_to_Anxiety_MR_methods_full.csv"))
cat("\n--- MR Results ---\n")
print(r1[, c("method", "b", "se", "pval", "OR", "OR_LCI", "OR_UCI")])

# Pleiotropy test (Egger intercept)
p1 <- mr_pleiotropy_test(dat1)
fwrite(p1, file.path(OUT, "MDD_to_Anxiety_pleiotropy.csv"))
cat("\n--- Pleiotropy Test (Egger Intercept) ---\n")
print(p1)

# Heterogeneity
h1 <- mr_heterogeneity(dat1)
fwrite(h1, file.path(OUT, "MDD_to_Anxiety_heterogeneity.csv"))
cat("\n--- Heterogeneity ---\n")
print(h1)

# Egger regression details
cat("\n--- Egger Detailed ---\n")
eg1 <- mr_egger_regression(
  bx = dat1$beta.exposure, bxse = dat1$se.exposure,
  by = dat1$beta.outcome, byse = dat1$se.outcome
)
cat("  Intercept: ", eg1[[1]]$intercept, " (SE=", eg1[[1]]$intercept_se, ")\n")
cat("  Slope:     ", eg1[[1]]$slope,    " (SE=", eg1[[1]]$slope_se,    ", pval=", eg1[[1]]$pval, ")\n")
cat("  Directional pleiotropy (|intercept|/SE > 2): ",
    ifelse(abs(eg1[[1]]$intercept / eg1[[1]]$intercept_se) > 2, "YES", "NO"), "\n")

# ============================================================
# DIRECTION 2: Anxiety -> MDD
# ============================================================
cat("\n=== Anxiety -> MDD (87 SNPs) ===\n")

h2 <- fread("03_harmonised/Anxiety_to_MDD_harmonised.csv")
d2 <- h2[h2$remove == FALSE]
cat(sprintf("  SNPs: %d\n", nrow(d2)))

dat2 <- data.frame(
  SNP                    = d2$SNP,
  exposure               = "Anxiety",
  id.exposure            = "Anxiety",
  beta.exposure          = d2$beta.exposure,
  se.exposure            = d2$se.exposure,
  pval.exposure         = d2$pval.exposure,
  effect_allele.exposure = d2$effect_allele.exposure,
  other_allele.exposure  = d2$other_allele.exposure,
  eaf.exposure          = d2$eaf.exposure,
  samplesize.exposure    = d2$samplesize.exposure,
  outcome               = "MDD",
  id.outcome            = "MDD",
  beta.outcome          = d2$beta.outcome,
  se.outcome            = d2$se.outcome,
  pval.outcome          = d2$pval.outcome,
  effect_allele.outcome = d2$effect_allele.outcome,
  other_allele.outcome  = d2$other_allele.outcome,
  eaf.outcome          = d2$eaf.outcome,
  samplesize.outcome    = d2$samplesize.outcome,
  mr_keep               = TRUE,
  data.type.exposure    = "easiro",
  data.type.outcome     = "easiro"
)
class(dat2) <- c("data.frame", "dat.trich")

r2 <- mr(dat2)
r2$OR     <- exp(r2$b)
r2$OR_LCI <- exp(r2$b - 1.96 * r2$se)
r2$OR_UCI <- exp(r2$b + 1.96 * r2$se)
fwrite(r2, file.path(OUT, "Anxiety_to_MDD_MR_methods_full.csv"))
cat("\n--- MR Results ---\n")
print(r2[, c("method", "b", "se", "pval", "OR", "OR_LCI", "OR_UCI")])

p2 <- mr_pleiotropy_test(dat2)
fwrite(p2, file.path(OUT, "Anxiety_to_MDD_pleiotropy.csv"))
cat("\n--- Pleiotropy Test (Egger Intercept) ---\n")
print(p2)

h2 <- mr_heterogeneity(dat2)
fwrite(h2, file.path(OUT, "Anxiety_to_MDD_heterogeneity.csv"))
cat("\n--- Heterogeneity ---\n")
print(h2)

cat("\n--- Egger Detailed ---\n")
eg2 <- mr_egger_regression(
  bx = dat2$beta.exposure, bxse = dat2$se.exposure,
  by = dat2$beta.outcome, byse = dat2$se.outcome
)
cat("  Intercept: ", eg2[[1]]$intercept, " (SE=", eg2[[1]]$intercept_se, ")\n")
cat("  Slope:     ", eg2[[1]]$slope,    " (SE=", eg2[[1]]$slope_se,    ", pval=", eg2[[1]]$pval, ")\n")
cat("  Directional pleiotropy (|intercept|/SE > 2): ",
    ifelse(abs(eg2[[1]]$intercept / eg2[[1]]$intercept_se) > 2, "YES", "NO"), "\n")

cat("\n==================================================\n")
cat("All files exported successfully.\n")
