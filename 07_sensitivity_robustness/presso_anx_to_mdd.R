#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(data.table); library(MRPRESSO) })
set.seed(20260522)

BASE <- "D:/2026年/文章/因果推断"
setwd(BASE)

h <- fread("03_harmonised/Anxiety_to_MDD_harmonised.csv")
# Filter kept SNPs
d <- h[h$remove == FALSE]

cat(sprintf("Running MR-PRESSO: %s SNPs\n", nrow(d)))

pdat <- data.frame(
  BetaOutcome  = d$beta.outcome,
  BetaExposure = d$beta.exposure,
  SdOutcome    = d$se.outcome,
  SdExposure   = d$se.exposure,
  SNP          = d$SNP
)

res <- mr_presso(
  BetaOutcome  = "BetaOutcome",
  BetaExposure = "BetaExposure",
  SdOutcome    = "SdOutcome",
  SdExposure   = "SdExposure",
  OUTLIERtest    = TRUE,
  DISTORTIONtest = TRUE,
  data             = pdat,
  NbDistribution  = 10000,
  SignifThreshold = 0.05
)

sink("05_MR_PRESSO/Anxiety_to_MDD_MRPRESSO_full.txt")
print(res)
sink()

gp <- res$`MR-PRESSO Results`$`P-value`[1]
cat(sprintf("Global P = %s\n", gp))
cat("Done.\n")
