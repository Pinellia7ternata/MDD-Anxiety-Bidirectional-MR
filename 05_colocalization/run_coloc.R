#!/usr/bin/env Rscript
# run_coloc.R
# Colocalization Analysis for MDD and Anxiety
# Author: Chen Chao | 2026-05-21

library(coloc)
library(TwoSampleMR)
library(dplyr)
library(data.table)

cat("=== COLOCALIZATION ANALYSIS ===\n")

# ==================================================
# Step 1: Download GWAS data from IEU OpenGWAS
# ==================================================
cat("\n[Step 1] Downloading GWAS data from IEU OpenGWAS...\n")

# MDD exposure ID (PGC MDD 2023)
mdd_id <- "ieu-a-7"

# Anxiety outcome ID (PGC Anxiety 2023)  
anxiety_id <- "ieu-a-117"

# Download data
cat("Downloading MDD GWAS summary stats...\n")
mdd_dat <- extract_outcome_data(
  snps = NULL,
  outcomes = mdd_id,
  plink_ld = NULL
)

cat("Downloading Anxiety GWAS summary stats...\n")
anxiety_dat <- extract_outcome_data(
  snps = NULL,
  outcomes = anxiety_id,
  plink_ld = NULL
)

# ==================================================
# Step 2: Format data for coloc
# ==================================================
cat("\n[Step 2] Formatting data for coloc...\n")

# Format for coloc
mdd_coloc <- list(
  pvalues = mdd_dat$pval,
  beta = mdd_dat$beta.outcome,
  N = mdd_dat$n.outcome,
  snp = mdd_dat$SNP,
  gene = NA
)

anxiety_coloc <- list(
  pvalues = anxiety_dat$pval,
  beta = anxiety_dat$beta.outcome,
  N = anxiety_dat$n.outcome,
  snp = anxiety_dat$SNP,
  gene = NA
)

# ==================================================
# Step 3: Run coloc.abf (Approximate Bayes Factor)
# ==================================================
cat("\n[Step 3] Running coloc.abf...\n")

coloc_res <- coloc.abf(dataset1 = mdd_coloc, dataset2 = anxiety_coloc)

print(coloc_res$summary)

# ==================================================
# Step 4: Run coloc.susie (if susieR available)
# ==================================================
cat("\n[Step 4] Running coloc.susie (fine-mapping)...\n")

if (requireNamespace("susieR", quietly = TRUE)) {
  coloc_susie <- coloc.susie(
    dataset1 = mdd_coloc,
    dataset2 = anxiety_coloc,
    p1 = 1e-4,
    p2 = 1e-4
  )
  print(coloc_susie$summary)
} else {
  cat("susieR not installed. Skipping coloc.susie.\n")
  cat("Install with: install.packages('susieR')\n")
}

# ==================================================
# Step 5: Save results
# ==================================================
cat("\n[Step 5] Saving results...\n")

# Save summary
sink("D:/2026年/文章/因果推断/MR_results/coloc_summary.txt")
print(coloc_res$summary)
sink()

# Save full results
write.csv(
  coloc_res$results,
  "D:/2026年/文章/因果推断/MR_results/coloc_results.csv",
  row.names = FALSE
)

cat("\n=== COLOCALIZATION COMPLETE ===\n")
cat("Results saved to: MR_results/coloc_summary.txt\n")
