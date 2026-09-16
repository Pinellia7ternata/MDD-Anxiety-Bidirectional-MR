#!/usr/bin/env Rscript
# -*- encoding: utf-8 -*-
# ============================================================
# STEP 3: Harmonization
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
})

BASE <- "D:/2026年/文章/因果推断"
setwd(BASE)

OUT_INSTR <- file.path(BASE, "02_instruments")
OUT_HARM <- file.path(BASE, "03_harmonised")

cat("[Step 3] Harmonization...\n")

# Load clumped instruments
mdd_clump <- fread(file.path(OUT_INSTR, "MDD_clumped_r2_0.001_10Mb.tsv"))
anx_clump <- fread(file.path(OUT_INSTR, "Anxiety_clumped_r2_0.001_10Mb.tsv"))

# Load full formatted data
exp_mdd <- fread(file.path(BASE, "01_formatted/MDD_PGC_2025_EUR_no23andMe_formatted.tsv.gz"))
exp_anx <- fread(file.path(BASE, "01_formatted/Anxiety_PGC_2026_EUR_formatted.tsv.gz"))

# --- Direction 1: MDD -> Anxiety ---
cat("\n  Direction 1: MDD -> Anxiety...\n")
exp_mdd_dat <- read_exposure_data(
  filename = file.path(OUT_INSTR, "MDD_clumped_r2_0.001_10Mb.tsv"),
  sep = "\t",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "ea.exposure",
  other_allele_col = "oa.exposure",
  eaf_col = "eaf.exposure",
  pval_col = "pval.exposure",
  samplesize_col = "samplesize.exposure"
)

out_anx_dat <- read_outcome_data(
  snps = exp_mdd_dat$SNP,
  filename = file.path(BASE, "01_formatted/Anxiety_PGC_2026_EUR_formatted.tsv.gz"),
  sep = "\t",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  samplesize_col = "samplesize"
)

harm1 <- harmonise_data(exp_mdd_dat, out_anx_dat, action = 2)
cat(sprintf("  Harmonized: %s SNPs\n", nrow(harm1)))
cat(sprintf("  Removed (strand flip/ambig): %s\n", sum(harm1$remove)))
fwrite(as.data.frame(harm1), file.path(OUT_HARM, "MDD_to_Anxiety_harmonised.csv"))

# --- Direction 2: Anxiety -> MDD ---
cat("\n  Direction 2: Anxiety -> MDD...\n")
exp_anx_dat <- read_exposure_data(
  filename = file.path(OUT_INSTR, "Anxiety_clumped_r2_0.001_10Mb.tsv"),
  sep = "\t",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "ea.exposure",
  other_allele_col = "oa.exposure",
  eaf_col = "eaf.exposure",
  pval_col = "pval.exposure",
  samplesize_col = "samplesize.exposure"
)

out_mdd_dat <- read_outcome_data(
  snps = exp_anx_dat$SNP,
  filename = file.path(BASE, "01_formatted/MDD_PGC_2025_EUR_no23andMe_formatted.tsv.gz"),
  sep = "\t",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  samplesize_col = "samplesize"
)

harm2 <- harmonise_data(exp_anx_dat, out_mdd_dat, action = 2)
cat(sprintf("  Harmonized: %s SNPs\n", nrow(harm2)))
cat(sprintf("  Removed: %s\n", sum(harm2$remove)))
fwrite(as.data.frame(harm2), file.path(OUT_HARM, "Anxiety_to_MDD_harmonised.csv"))

cat("\n[Step 3] Harmonization done.\n")
