#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(TwoSampleMR)
})

cat("==================================================\n")
cat("STRICT MR: Full Results Export\n")
cat("==================================================\n\n")

BASE <- "D:/2026年/文章/因果推断"
OUT_UMR <- file.path(BASE, "04_univariable_MR")
setwd(BASE)

# ============================================================
# DIRECTION 1: MDD -> Anxiety
# ============================================================
cat("=== MDD -> Anxiety ===\n")

harm1 <- fread("03_harmonised/MDD_to_Anxiety_harmonised.csv")
dat1  <- harm1[harm1$remove == FALSE, ]
cat(sprintf("  Harmonised SNPs: %d\n", nrow(dat1)))

# TwoSampleMR format
exp1 <- dat1[, .(SNP, beta.exposure, se.exposure, pval.exposure,
                  effect_allele.exposure, other_allele.exposure,
                  eaf.exposure, samplesize.exposure)]
out1 <- dat1[, .(SNP, beta.outcome, se.outcome, pval.outcome, samplesize.outcome)]

dat1_mr <- harmonise_data(exp1, out1)

# MR results (all methods)
mr_res1 <- mr(dat1_mr)
mr_res1$OR <- exp(mr_res1$b)
mr_res1$OR_LCI <- exp(mr_res1$b - 1.96 * mr_res1$se)
mr_res1$OR_UCI <- exp(mr_res1$b + 1.96 * mr_res1$se)
fwrite(mr_res1, file.path(OUT_UMR, "MDD_to_Anxiety_MR_methods_full.csv"))
print(mr_res1)

# MR-Egger regression (detailed)
cat("\n--- MR-Egger Detailed (MDD->Anxiety) ---\n")
egger1 <- mr_egger_regression(
  bx   = dat1_mr$beta.exposure,
  bxse = dat1_mr$se.exposure,
  by   = dat1_mr$beta.outcome,
  byse = dat1_mr$se.outcome
)
print(egger1)
saveRDS(egger1, file.path(OUT_UMR, "MDD_to_Anxiety_egger_regression.rds"))

# Pleiotropy test (Egger intercept)
pleio1 <- mr_pleiotropy_test(dat1_mr)
pleio1$direction <- "MDD->Anxiety"
fwrite(pleio1, file.path(OUT_UMR, "MDD_to_Anxiety_pleiotropy.csv"))
cat("\n--- Pleiotropy Test (MDD->Anxiety) ---\n")
print(pleio1)

# Heterogeneity
het1 <- mr_heterogeneity(dat1_mr)
fwrite(het1, file.path(OUT_UMR, "MDD_to_Anxiety_heterogeneity.csv"))
cat("\n--- Heterogeneity (MDD->Anxiety) ---\n")
print(het1)

# Leave-one-out
loo1 <- mr_leaveoneout(dat1_mr)
fwrite(loo1, file.path(OUT_UMR, "MDD_to_Anxiety_leaveoneout_full.csv"))

# Steiger
steiger1 <- mr_steiger(
  exp_b = dat1_mr$beta.exposure, out_b = dat1_mr$beta.outcome,
  exp_se = dat1_mr$se.exposure, out_se = dat1_mr$se.outcome,
  exp_n = dat1_mr$samplesize.exposure[1],
  out_n = dat1_mr$samplesize.outcome[1]
)
fwrite(as.data.table(steiger1), file.path(OUT_UMR, "MDD_to_Anxiety_steiger.csv"))
cat("\n--- Steiger (MDD->Anxiety) ---\n")
print(steiger1)

cat("\n")

# ============================================================
# DIRECTION 2: Anxiety -> MDD
# ============================================================
cat("=== Anxiety -> MDD ===\n")

harm2 <- fread("03_harmonised/Anxiety_to_MDD_harmonised.csv")
dat2  <- harm2[harm2$remove == FALSE, ]
cat(sprintf("  Harmonised SNPs: %d\n", nrow(dat2)))

exp2 <- dat2[, .(SNP, beta.exposure, se.exposure, pval.exposure,
                  effect_allele.exposure, other_allele.exposure,
                  eaf.exposure, samplesize.exposure)]
out2 <- dat2[, .(SNP, beta.outcome, se.outcome, pval.outcome, samplesize.outcome)]

dat2_mr <- harmonise_data(exp2, out2)

mr_res2 <- mr(dat2_mr)
mr_res2$OR <- exp(mr_res2$b)
mr_res2$OR_LCI <- exp(mr_res2$b - 1.96 * mr_res2$se)
mr_res2$OR_UCI <- exp(mr_res2$b + 1.96 * mr_res2$se)
fwrite(mr_res2, file.path(OUT_UMR, "Anxiety_to_MDD_MR_methods_full.csv"))
print(mr_res2)

cat("\n--- MR-Egger Detailed (Anxiety->MDD) ---\n")
egger2 <- mr_egger_regression(
  bx   = dat2_mr$beta.exposure,
  bxse = dat2_mr$se.exposure,
  by   = dat2_mr$beta.outcome,
  byse = dat2_mr$se.outcome
)
print(egger2)
saveRDS(egger2, file.path(OUT_UMR, "Anxiety_to_MDD_egger_regression.rds"))

pleio2 <- mr_pleiotropy_test(dat2_mr)
pleio2$direction <- "Anxiety->MDD"
fwrite(pleio2, file.path(OUT_UMR, "Anxiety_to_MDD_pleiotropy.csv"))
cat("\n--- Pleiotropy Test (Anxiety->MDD) ---\n")
print(pleio2)

het2 <- mr_heterogeneity(dat2_mr)
fwrite(het2, file.path(OUT_UMR, "Anxiety_to_MDD_heterogeneity.csv"))

loo2 <- mr_leaveoneout(dat2_mr)
fwrite(loo2, file.path(OUT_UMR, "Anxiety_to_MDD_leaveoneout_full.csv"))

steiger2 <- mr_steiger(
  exp_b = dat2_mr$beta.exposure, out_b = dat2_mr$beta.outcome,
  exp_se = dat2_mr$se.exposure, out_se = dat2_mr$se.outcome,
  exp_n = dat2_mr$samplesize.exposure[1],
  out_n = dat2_mr$samplesize.outcome[1]
)
fwrite(as.data.table(steiger2), file.path(OUT_UMR, "Anxiety_to_MDD_steiger.csv"))
cat("\n--- Steiger (Anxiety->MDD) ---\n")
print(steiger2)

cat("\n==================================================\n")
cat("Done.\n")
cat("==================================================\n")
