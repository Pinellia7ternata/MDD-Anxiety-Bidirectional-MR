#!/usr/bin/env Rscript
# run_coloc_v2.R
# Colocalization Analysis using IEU OpenGWAS
# Author: Chen Chao | 2026-05-21

library(coloc)
library(TwoSampleMR)
library(dplyr)

cat("=== COLOCALIZATION ANALYSIS (IEU OpenGWAS) ===\n")

# ==================================================
# Step 1: Extract MDD instruments from IEU
# ==================================================
cat("\n[Step 1] Extracting MDD instruments from IEU...\n")

# MDD exposure (ieu-a-7: PGC MDD 2018)
mdd_exp <- extract_instruments(
  outcomes = 'ieu-a-7',
  p1 = 5e-8,
  clump = TRUE,
  r2 = 0.001,
  kb = 10000
)

cat(sprintf("MDD instruments extracted: %d SNPs\n", nrow(mdd_exp)))

# ==================================================
# Step 2: Extract Anxiety outcome data
# ==================================================
cat("\n[Step 2] Extracting Anxiety outcome data...\n")

# Anxiety outcome (ieu-a-117: PGC Anxiety 2018)
anxiety_out <- extract_outcome_data(
  snps = mdd_exp$SNP,
  outcomes = 'ieu-a-117'
)

cat(sprintf("Anxiety outcome data: %d SNPs\n", nrow(anxiety_out)))

# ==================================================
# Step 3: Harmonize data
# ==================================================
cat("\n[Step 3] Harmonizing data...\n")

harmonized <- harmonise_data(
  exposure_dat = mdd_exp,
  outcome_dat = anxiety_out
)

cat(sprintf("Harmonized pairs: %d SNPs\n", nrow(harmonized)))

# ==================================================
# Step 4: Format for coloc
# ==================================================
cat("\n[Step 4] Formatting for coloc...\n")

mdd_coloc <- list(
  pvalues = harmonized$pval.exposure,
  beta = harmonized$beta.exposure,
  N = harmonized$samplesize.exposure,
  snp = harmonized$SNP,
  gene = NA
)

anxiety_coloc <- list(
  pvalues = harmonized$pval.outcome,
  beta = harmonized$beta.outcome,
  N = harmonized$samplesize.outcome,
  snp = harmonized$SNP,
  gene = NA
)

# ==================================================
# Step 5: Run coloc.abf
# ==================================================
cat("\n[Step 5] Running coloc.abf...\n")

coloc_res <- coloc.abf(dataset1 = mdd_coloc, dataset2 = anxiety_coloc)

print(coloc_res$summary)

# ==================================================
# Step 6: Save results
# ==================================================
cat("\n[Step 6] Saving results...\n")

sink("D:/2026年/文章/因果推断/MR_results/coloc_IEU_summary.txt")
print(coloc_res$summary)
sink()

# Save detailed results
write.csv(
  coloc_res$results,
  "D:/2026年/文章/因果推断/MR_results/coloc_IEU_results.csv",
  row.names = FALSE
)

# Save harmonized data
write.csv(
  harmonized,
  "D:/2026年/文章/因果推断/MR_results/coloc_IEU_harmonized.csv",
  row.names = FALSE
)

cat("\n=== COLOCALIZATION COMPLETE ===\n")
cat("Results saved to: MR_results/coloc_IEU_*\n")
cat(sprintf("\nSummary:\n"))
cat(sprintf("  PP.H0 (no association): %.3f\n", coloc_res$summary[1]))
cat(sprintf("  PP.H1 (MDD only):      %.3f\n", coloc_res$summary[2]))
cat(sprintf("  PP.H2 (Anxiety only):  %.3f\n", coloc_res$summary[3]))
cat(sprintf("  PP.H3 (two distinct):  %.3f\n", coloc_res$summary[4]))
cat(sprintf("  PP.H4 (shared causal): %.3f\n", coloc_res$summary[5]))
