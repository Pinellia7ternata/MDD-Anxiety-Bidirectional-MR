#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(data.table); library(MRPRESSO) })
set.seed(20260522)
BASE <- "D:/2026年/文章/因果推断"
setwd(BASE)

OUT_HARM   <- file.path(BASE, "03_harmonised")
OUT_PRESSO <- file.path(BASE, "05_MR_PRESSO")
dir.create(OUT_PRESSO, showWarnings=FALSE)

for (f in c("MDD_to_Anxiety", "Anxiety_to_MDD")) {
  cat(sprintf("[%s] Running PRESSO...\n", f))
  h <- fread(file.path(OUT_HARM, sprintf("%s_harmonised.csv", f)))
  d <- h[!h$remove, ]
  pdat <- data.frame(BetaOutcome=d$beta.outcome, BetaExposure=d$beta.exposure,
                     SdOutcome=d$se.outcome, SdExposure=d$se.exposure, SNP=d$SNP)
  res <- mr_presso(BetaOutcome="BetaOutcome", BetaExposure="BetaExposure",
                    SdOutcome="SdOutcome",     SdExposure="SdExposure",
                    OUTLIERtest=TRUE, DISTORTIONtest=FALSE,
                    data=pdat, NbDistribution=100, SignifThreshold=0.05)
  sink(file.path(OUT_PRESSO, sprintf("%s_MRPRESSO_quick.txt", f)))
  print(res)
  sink()
  gp <- res$`MR-PRESSO Results`$`P-value`[1]
  cat(sprintf("  Global P=%.3e, Outliers=%d\n", gp, nrow(res$`MR-PRESSO Outlier Test`)))
}
cat("[Done]\n")
