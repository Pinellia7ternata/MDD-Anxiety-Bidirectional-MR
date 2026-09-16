#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
})

BASE <- "D:/2026年/文章/因果推断"
OUT  <- file.path(BASE, "04_univariable_MR")
setwd(BASE)

h2 <- fread("03_harmonised/Anxiety_to_MDD_harmonised.csv")
d2 <- h2[h2$remove == FALSE]
cat(sprintf("Anxiety->MDD: %d SNPs\n", nrow(d2)))

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

cat("\nDone.\n")
