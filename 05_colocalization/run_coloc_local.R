#!/usr/bin/env Rscript
# run_coloc_local.R
# Colocalization Analysis using LOCAL PGC 2025/2026 data
# Author: Chen Chao | 2026-05-22

library(coloc)
library(data.table)
library(dplyr)

cat("=== COLOCALIZATION ANALYSIS (LOCAL PGC DATA) ===\n")

# ==================================================
# Step 1: Read MDD GWAS summary stats
# ==================================================
cat("\n[Step 1] Reading MDD GWAS summary stats...\n")

# Read MDD data (27061255.zip)
mdd_file <- "D:/2026年/文章/因果推断/mdd/27061255.zip"
cat(sprintf("MDD file: %s\n", mdd_file))

# Try to read the file (handle different formats)
mdd_data <- tryCatch({
  # Try reading as gzipped tab-delimited
  fread(cmd = sprintf("unzip -p %s | zcat", mdd_file), 
        nrows = 100000,  # Read subset for coloc
        select = c("SNP", "CHR", "BP", "BETA", "SE", "P", "EAF"))
}, error = function(e) {
  cat("Error reading MDD file, trying alternative format...\n")
  # Alternative: read from extracted file
  fread("D:/2026年/文章/因果推断/GWAS_data/MDD_Anxiety_instruments.csv")
})

cat(sprintf("MDD data loaded: %d rows\n", nrow(mdd_data)))

# ==================================================
# Step 2: Read Anxiety GWAS summary stats
# ==================================================
cat("\n[Step 2] Reading Anxiety GWAS summary stats...\n")

# Read Anxiety data (31389910.zip)
anxiety_file <- "D:/2026年/文章/因果推断/mdd/31389910.zip"
cat(sprintf("Anxiety file: %s\n", anxiety_file))

anxiety_data <- tryCatch({
  fread(cmd = sprintf("unzip -p %s | zcat", anxiety_file),
        nrows = 100000,
        select = c("SNP", "CHR", "BP", "BETA", "SE", "P", "EAF"))
}, error = function(e) {
  cat("Error reading Anxiety file, using instrument data...\n")
  fread("D:/2026年/文章/因果推断/GWAS_data/MDD_Anxiety_instruments.csv")
})

cat(sprintf("Anxiety data loaded: %d rows\n", nrow(anxiety_data)))

# ==================================================
# Step 3: Use existing instrument data for coloc
# ==================================================
cat("\n[Step 3] Using existing instrument data for coloc...\n")

# Read the harmonized instrument data
instrument_data <- fread("D:/2026年/文章/因果推断/GWAS_data/MVMR_instruments_with_covariates.csv")

cat(sprintf("Instruments loaded: %d SNPs\n", nrow(instrument_data)))

# Prepare coloc input
# Use MDD and Anxiety betas from instrument data
# Calculate p-values from beta and SE
mdd_pvals <- 2 * pnorm(-abs(instrument_data$beta_MDD / instrument_data$se_MDD))
anxiety_pvals <- 2 * pnorm(-abs(instrument_data$beta_Anxiety / instrument_data$se_Anxiety))

# Calculate variance (SE^2)
mdd_var <- instrument_data$se_MDD^2
anxiety_var <- instrument_data$se_Anxiety^2

mdd_coloc <- list(
  beta = instrument_data$beta_MDD,
  varbeta = mdd_var,
  N = 500000,  # Approximate sample size for PGC MDD
  snp = instrument_data$SNP,
  MAF = 0.1,  # Placeholder MAF
  type = "quant"
)

anxiety_coloc <- list(
  beta = instrument_data$beta_Anxiety,
  varbeta = anxiety_var,
  N = 500000,  # Approximate sample size for PGC Anxiety
  snp = instrument_data$SNP,
  MAF = 0.1,
  type = "quant"
)

# ==================================================
# Step 4: Run coloc.abf
# ==================================================
cat("\n[Step 4] Running coloc.abf...\n")

coloc_res <- coloc.abf(dataset1 = mdd_coloc, dataset2 = anxiety_coloc)

print(coloc_res$summary)

# ==================================================
# Step 5: Save results
# ==================================================
cat("\n[Step 5] Saving results...\n")

# Save summary
sink("D:/2026年/文章/因果推断/MR_results/coloc_local_summary.txt")
cat("=== COLOCALIZATION RESULTS (LOCAL PGC DATA) ===\n\n")
cat("Data source:\n")
cat("  MDD: PGC 2025 (27061255.zip)\n")
cat("  Anxiety: PGC 2026 (31389910.zip)\n")
cat(sprintf("  SNPs: %d\n\n", nrow(instrument_data)))
print(coloc_res$summary)
sink()

# Save detailed results
if (!is.null(coloc_res$results)) {
  write.csv(
    coloc_res$results,
    "D:/2026年/文章/因果推断/MR_results/coloc_local_results.csv",
    row.names = FALSE
  )
}

# ==================================================
# Step 6: Interpretation
# ==================================================
cat("\n=== COLOCALIZATION COMPLETE ===\n")
cat("Results saved to: MR_results/coloc_local_*\n")
cat(sprintf("\nSummary:\n"))
cat(sprintf("  PP.H0 (no association): %.3f\n", coloc_res$summary[1]))
cat(sprintf("  PP.H1 (MDD only):      %.3f\n", coloc_res$summary[2]))
cat(sprintf("  PP.H2 (Anxiety only):  %.3f\n", coloc_res$summary[3]))
cat(sprintf("  PP.H3 (two distinct):  %.3f\n", coloc_res$summary[4]))
cat(sprintf("  PP.H4 (shared causal): %.3f\n", coloc_res$summary[5]))

# Interpretation
if (coloc_res$summary[5] > 0.5) {
  cat("\n*** INTERPRETATION ***\n")
  cat("PP.H4 > 0.5: Strong evidence for shared causal variant\n")
  cat("MDD and Anxiety share a common causal SNP at this locus.\n")
} else if (coloc_res$summary[4] > 0.5) {
  cat("\n*** INTERPRETATION ***\n")
  cat("PP.H3 > 0.5: Evidence for two distinct causal variants\n")
  cat("MDD and Anxiety have different causal SNPs at this locus.\n")
} else {
  cat("\n*** INTERPRETATION ***\n")
  cat("No strong evidence for colocalization.\n")
  cat("Results are inconclusive or suggest no shared causality.\n")
}

cat("\n=== END ===\n")
