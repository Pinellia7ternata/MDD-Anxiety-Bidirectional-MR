#!/usr/bin/env Rscript
# UKB External Validation Analysis (V3)
# Using CHR:POS coordinate matching

suppressPackageStartupMessages(library(data.table))

BASE <- "D:/2026年/文章/因果推断"
OUT  <- "D:/2026年/文章/因果推断/10_External_Validation"
set.seed(2026)

cat("=== UKB External Validation ===\n\n")

# ============================================================
# Step 1: Load PGC IVs with coordinates
# ============================================================
cat("Loading PGC IVs with coordinates...\n")
pggc_mdd <- fread(file.path(BASE, "01_formatted/MDD_PGC_2025_EUR_no23andMe_formatted.tsv.gz"),
                  select = c("SNP", "chr", "pos", "effect_allele", "other_allele", "beta", "se", "eaf"))
pggc_anx <- fread(file.path(BASE, "01_formatted/Anxiety_PGC_2026_EUR_formatted.tsv.gz"),
                 select = c("SNP", "chr", "pos", "effect_allele", "other_allele", "beta", "se", "eaf"))

# Load instrument SNPs
ivs_mdd <- fread(file.path(BASE, "02_instruments_realLD/MDD_realLD_clumped.tsv"))$SNP
ivs_anx <- fread(file.path(BASE, "02_instruments_realLD/Anxiety_realLD_clumped.tsv"))$SNP

# Filter to IVs
pggc_mdd_ivs <- pggc_mdd[SNP %in% ivs_mdd]
pggc_anx_ivs <- pggc_anx[SNP %in% ivs_anx]

# Create coordinate key
pggc_mdd_ivs[, coord := paste(chr, pos, sep = ":")]
pggc_anx_ivs[, coord := paste(chr, pos, sep = ":")]

cat(sprintf("PGC MDD IVs with coords: %d\n", nrow(pggc_mdd_ivs)))
cat(sprintf("PGC Anxiety IVs with coords: %d\n", nrow(pggc_anx_ivs)))

# ============================================================
# Step 2: Load UKB data and create coordinate key
# ============================================================
cat("\nLoading UKB 20002_1287 (~582 MB)...\n")
ukb_1287 <- fread(file.path(OUT, "20002_1287.gwas.imputed_v3.both_sexes.tsv.bgz"),
                  select = c("variant", "beta", "se", "pval", "minor_AF", "n_complete_samples"))
# Parse variant: CHR:POS:REF:ALT -> CHR:POS
ukb_1287[, coord := gsub(":([^:]+)$", "", variant)]  # Remove REF:ALT suffix
ukb_1287[, minor_allele := NULL]  # Free memory

cat(sprintf("UKB 1287 total: %s variants\n", format(nrow(ukb_1287), big.mark = ",")))

cat("\nLoading UKB 20002_1286 (~582 MB)...\n")
ukb_1286 <- fread(file.path(OUT, "20002_1286.gwas.imputed_v3.both_sexes.tsv.bgz"),
                  select = c("variant", "beta", "se", "pval", "minor_AF", "n_complete_samples"))
ukb_1286[, coord := gsub(":([^:]+)$", "", variant)]

cat(sprintf("UKB 1286 total: %s variants\n", format(nrow(ukb_1286), big.mark = ",")))

# ============================================================
# Step 3: Coordinate-based matching
# ============================================================
setkey(pggc_mdd_ivs, coord)
setkey(pggc_anx_ivs, coord)
setkey(ukb_1287, coord)
setkey(ukb_1286, coord)

# Matches: PGC MDD IVs in UKB 1287
matches_mdd_1287 <- pggc_mdd_ivs[ukb_1287, nomatch = 0]
cat(sprintf("\nPGC MDD IVs found in UKB 1287: %d/%d (%.1f%%)\n", 
           nrow(matches_mdd_1287), length(ivs_mdd), 
           100 * nrow(matches_mdd_1287) / length(ivs_mdd)))

# Matches: PGC Anxiety IVs in UKB 1286
matches_anx_1286 <- pggc_anx_ivs[ukb_1286, nomatch = 0]
cat(sprintf("PGC Anxiety IVs found in UKB 1286: %d/%d (%.1f%%)\n", 
           nrow(matches_anx_1286), length(ivs_anx),
           100 * nrow(matches_anx_1286) / length(ivs_anx)))

# ============================================================
# Step 4: Harmonize alleles and calculate MR
# ============================================================
run_mr_direction <- function(exp_dat, out_dat, exp_name, out_name) {
  cat(sprintf("\n=== %s -> %s ===\n", exp_name, out_name))
  
  # Need to check allele alignment - simplified version
  # For now, use the beta directly assuming aligned direction
  
  n <- nrow(exp_dat)
  if (n < 10) {
    cat("Too few SNPs!\n")
    return(NULL)
  }
  
  # Simple IVW
  w <- 1 / exp_dat$se ^ 2
  beta_ivw <- sum(w * exp_dat$beta * exp_dat$i.beta) / sum(w * exp_dat$beta ^ 2)
  se_ivw <- sqrt(1 / sum(w * exp_dat$beta ^ 2))
  p_ivw <- 2 * pnorm(-abs(beta_ivw / se_ivw))
  
  # Weighted median
  ratios <- exp_dat$i.beta / exp_dat$beta
  b_wm <- median(ratios)
  se_wm <- mad(ratios) / sqrt(n)
  p_wm <- 2 * pnorm(-abs(b_wm / se_wm))
  
  res <- data.frame(
    Direction = paste(exp_name, "->", out_name),
    Method = c("IVW", "Weighted Median"),
    N_SNPs = c(n, n),
    OR = c(exp(beta_ivw), exp(b_wm)),
    LCI = c(exp(beta_ivw - 1.96 * se_ivw), NA),
    UCI = c(exp(beta_ivw + 1.96 * se_ivw), NA),
    P = c(p_ivw, p_wm)
  )
  
  print(res, digits = 4)
  return(res)
}

# Note: i.beta is the inner join variable from UKB
# Since we can't reliably align alleles without full harmonization,
# we'll just report the overlap statistics for now

cat("\n\n=== COORDINATE MATCHING COMPLETE ===\n")
cat("For proper MR analysis, need allele harmonization between:\n")
cat("  - PGC effect allele (effect_allele) vs UKB effect (minor_allele)\n")
cat("This requires detailed allele flipping logic.\n")

# Save match statistics for later
match_stats <- data.frame(
  Dataset = c("UKB_1287", "UKB_1286"),
  Phenotype_ID = c("20002_1287", "20002_1286"),
  PGC_Exposure = c("MDD", "Anxiety"),
  N_Available = c(length(ivs_mdd), length(ivs_anx)),
  N_Matched = c(nrow(matches_mdd_1287), nrow(matches_anx_1286)),
  Match_Rate = c(100 * nrow(matches_mdd_1287) / length(ivs_mdd),
                100 * nrow(matches_anx_1286) / length(ivs_anx)),
  stringsAsFactors = FALSE
)

print(match_stats)
write.csv(match_stats, file.path(OUT, "ukb_match_statistics.csv"), row.names = FALSE)

cat("\nResults saved to:", OUT, "\n")
