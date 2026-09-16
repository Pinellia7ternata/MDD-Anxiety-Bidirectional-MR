#!/usr/bin/env Rscript
# -*- encoding: utf-8 -*-
# ============================================================
# MDD ↔ Anxiety 严格 MR 重跑
# 严格LD Clumping: P<5e-8, r²<0.001, window=10Mb
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(TwoSampleMR)
  library(MRPRESSO)
  library(MendelianRandomization)
})

cat("==================================================\n")
cat("STRICT MR ANALYSIS: MDD <-> Anxiety\n")
cat("Clumping: P<5e-8, r2=0.001, window=10Mb\n")
cat("==================================================\n\n")

BASE <- "D:/2026年/文章/因果推断"
setwd(BASE)

OUT_INSTR <- file.path(BASE, "02_instruments")
OUT_HARM <- file.path(BASE, "03_harmonised")
OUT_MR <- file.path(BASE, "04_univariable_MR")
OUT_PRESSO <- file.path(BASE, "05_MR_PRESSO")
OUT_MVMR <- file.path(BASE, "06_MVMR")

dir.create(OUT_INSTR, showWarnings=FALSE)
dir.create(OUT_HARM, showWarnings=FALSE)
dir.create(OUT_MR, showWarnings=FALSE)
dir.create(OUT_PRESSO, showWarnings=FALSE)
dir.create(OUT_MVMR, showWarnings=FALSE)

# ============================================================
# STEP 1: 加载格式化数据
# ============================================================
cat("[Step 1] Loading formatted GWAS data...\n")

exp_mdd <- fread(file.path(BASE, "01_formatted/MDD_PGC_2025_EUR_no23andMe_formatted.tsv.gz"))
exp_anx <- fread(file.path(BASE, "01_formatted/Anxiety_PGC_2026_EUR_formatted.tsv.gz"))

cat(sprintf("  MDD: %s SNPs\n", format(nrow(exp_mdd), big.mark=",")))
cat(sprintf("  Anxiety: %s SNPs\n", format(nrow(exp_anx), big.mark=",")))

# ============================================================
# STEP 2: LD Clumping (严格: r²<0.001, window=10Mb)
# ============================================================
cat("\n[Step 2] LD Clumping (r2=0.001, window=10Mb, EUR)...\n")
cat("  NOTE: TwoSampleMR::clump_data uses 1000G EUR reference\n")
cat("  NOTE: Very strict r2=0.001 may yield very few SNPs\n\n")

# --- MDD clump (for MDD -> Anxiety) ---
cat("  Clumping MDD instruments...\n")
mdd_sig <- exp_mdd[pval < 5e-8]
cat(sprintf("  MDD P<5e-8: %s SNPs\n", nrow(mdd_sig)))

mdd_for_clump <- data.frame(
  SNP = mdd_sig$SNP,
  pval.exposure = mdd_sig$pval,
  chr_name = mdd_sig$chr,
  chrom_start = mdd_sig$pos,
  id.exposure = "MDD_PGC2025"
)

tryCatch({
  mdd_clumped <- clump_data(mdd_for_clump, clump_kb = 10000, clump_r2 = 0.001, pop = "EUR")
  cat(sprintf("  MDD after clumping: %s SNPs\n", nrow(mdd_clumped)))
  fwrite(mdd_clumped, file.path(OUT_INSTR, "MDD_clumped_r2_0.001_10Mb.tsv"), sep = "\t")
}, error = function(e) {
  cat("  ERROR in clump_data:", conditionMessage(e), "\n")
  cat("  Trying alternative: save for PLINK clumping...\n")
  write.table(mdd_for_clump[, c("SNP","pval.exposure","chr_name","chrom_start")], 
              file.path(OUT_INSTR, "MDD_for_clump.txt"), sep="\t", row.names=FALSE, quote=FALSE)
})

# --- Anxiety clump (for Anxiety -> MDD) ---
cat("\n  Clumping Anxiety instruments...\n")
anx_sig <- exp_anx[pval < 5e-8]
cat(sprintf("  Anxiety P<5e-8: %s SNPs\n", nrow(anx_sig)))

anx_for_clump <- data.frame(
  SNP = anx_sig$SNP,
  pval.exposure = anx_sig$pval,
  chr_name = anx_sig$chr,
  chrom_start = anx_sig$pos,
  id.exposure = "Anxiety_PGC2026"
)

tryCatch({
  anx_clumped <- clump_data(anx_for_clump, clump_kb = 10000, clump_r2 = 0.001, pop = "EUR")
  cat(sprintf("  Anxiety after clumping: %s SNPs\n", nrow(anx_clumped)))
  fwrite(anx_clumped, file.path(OUT_INSTR, "Anxiety_clumped_r2_0.001_10Mb.tsv"), sep = "\t")
}, error = function(e) {
  cat("  ERROR in clump_data:", conditionMessage(e), "\n")
})

cat("\n[Step 2] Clumping done.\n")
cat("==================================================\n")
