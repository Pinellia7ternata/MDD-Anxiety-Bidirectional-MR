#!/usr/bin/env Rscript
# UKB External Validation - Standardized Beta MR
# Key fix: UKB beta is per-unit of the trait (standardized, ~1 SD),
# while PGC beta is log-OR. Use ratio method.

suppressPackageStartupMessages(library(data.table))

BASE <- "D:/2026年/文章/因果推断"
OUT  <- file.path(BASE, "10_External_Validation")

cat("=== UKB External Validation (Standardized) ===\n\n")

# Load IVs
ivs_mdd <- fread(file.path(BASE, "02_instruments_realLD/MDD_realLD_clumped.tsv"))$SNP
ivs_anx <- fread(file.path(BASE, "02_instruments_realLD/Anxiety_realLD_clumped.tsv"))$SNP

# Load PGC data (exposure)
pggc_mdd <- fread(file.path(BASE, "01_formatted/MDD_PGC_2025_EUR_no23andMe_formatted.tsv.gz"),
                  select = c("SNP", "effect_allele", "other_allele", "beta", "se", "eaf"))
pggc_anx <- fread(file.path(BASE, "01_formatted/Anxiety_PGC_2026_EUR_formatted.tsv.gz"),
                  select = c("SNP", "effect_allele", "other_allele", "beta", "se", "eaf"))

# Load UKB with rsID
ukb_dep <- fread(file.path(OUT, "UKB_1287_depression_with_rsid.tsv.gz"),
                  select = c("rsid", "beta", "se", "pval", "minor_AF"))
setnames(ukb_dep, c("SNP", "beta_ukb", "se_ukb", "pval_ukb", "eaf_ukb"))

cat(sprintf("PGC MDD IVs: %d | Anxiety IVs: %d\n", length(ivs_mdd), length(ivs_anx)))

# ============================================================
# MR using RATIO METHOD (b_y/b_x, se = b_y/(b_x*sqrt(n)))
# This handles different beta scales naturally
# ============================================================
run_ratio_mr <- function(exp_snps, out_dat, exp_name, out_name, direction) {
  cat(sprintf("\n=== %s: %s -> %s ===\n", direction, exp_name, out_name))
  
  # Get exposure SNPs data
  if (exp_name == "PGC_MDD") {
    exp_dat <- pggc_mdd[SNP %in% exp_snps]
  } else if (exp_name == "PGC_Anxiety") {
    exp_dat <- pggc_anx[SNP %in% exp_snps]
  } else {
    exp_dat <- out_dat  # placeholder, handled below
  }
  
  # Get outcome data for these SNPs
  out_snp_dat <- out_dat[SNP %in% exp_dat$SNP]
  
  # Merge
  merged <- merge(exp_dat, out_snp_dat, by = "SNP", all = FALSE)
  n <- nrow(merged)
  cat(sprintf("Overlapping SNPs: %d\n", n))
  
  if (n < 10) return(NULL)
  
  # Method 1: IVW ratio (b_y/b_x)
  # b_ratio = b_out / b_exp; SE via delta method
  merged[, ratio := beta_ukb / beta]
  merged[, ratio_se := sqrt(se_ukb^2 / beta^2 + beta_ukb^2 * se^2 / beta^4)]
  
  b_ivw <- sum(merged$ratio / merged$ratio_se^2) / sum(1 / merged$ratio_se^2)
  se_ivw <- sqrt(1 / sum(1 / merged$ratio_se^2))
  p_ivw <- 2 * pnorm(-abs(b_ivw / se_ivw))
  
  # Method 2: Weighted median ratio
  b_wm <- median(merged$ratio)
  se_wm <- mad(merged$ratio) / sqrt(n)
  p_wm <- 2 * pnorm(-abs(b_wm / se_wm))
  
  # Method 3: Egger on ratios
  x <- merged$beta
  y <- merged$beta_ukb
  w <- 1 / merged$ratio_se^2
  X <- cbind(1, x)
  W <- diag(w)
  eg <- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% y
  b_eg <- as.numeric(eg[2])
  int_eg <- as.numeric(eg[1])
  vcov <- solve(t(X) %*% W %*% X)
  se_eg <- sqrt(as.numeric(vcov[2,2]))
  se_int <- sqrt(as.numeric(vcov[1,1]))
  p_eg <- 2 * pnorm(-abs(b_eg / se_eg))
  p_int <- 2 * pnorm(-abs(int_eg / se_int))
  
  # Report results (ratio = b_ukb / b_pgc = standardized-to-logOR ratio)
  # OR = exp(b_ratio * sd_x) where sd_x ~ 1 for MDD binary
  # Since b_pgc is log-OR per SD of liability, interpret as:
  # b_ratio = SD_change_in_outcome_per_SD_change_in_exposure
  res <- data.frame(
    Direction = direction,
    Method = c("IVW", "Weighted Median", "Egger"),
    N_SNPs = n,
    Beta_ratio = c(b_ivw, b_wm, b_eg),
    SE = c(se_ivw, se_wm, se_eg),
    P = c(p_ivw, p_wm, p_eg),
    Int_P = c(NA, NA, p_int),
    stringsAsFactors = FALSE
  )
  
  print(res, digits = 4)
  return(res)
}

# ============================================================
# Run analyses
# ============================================================
all_results <- list()

# 1. PGC MDD IVs -> UKB Depression
all_results[[1]] <- run_ratio_mr(ivs_mdd, ukb_dep, "PGC_MDD", "UKB_Depression", "MDD->UKB_DEP")

# 2. PGC Anxiety IVs -> UKB Depression
all_results[[2]] <- run_ratio_mr(ivs_anx, ukb_dep, "PGC_Anxiety", "UKB_Depression", "Anxiety->UKB_DEP")

# ============================================================
# Compare with primary analysis
# ============================================================
cat("\n\n========== EXTERNAL VALIDATION SUMMARY ==========\n")
primary <- data.frame(
  Direction = c("MDD->Anxiety", "Anxiety->MDD"),
  Primary_Beta = c(log(2.343), log(1.682)),
  Primary_P = c(0, 2.99e-113),
  stringsAsFactors = FALSE
)

comparison <- data.frame(
  Analysis = c("Primary: PGC MDD -> PGC Anxiety", 
               "Validation: PGC MDD -> UKB Depression",
               "Primary: PGC Anxiety -> PGC MDD",
               "Validation: PGC Anxiety -> UKB Depression"),
  Direction = c("MDD->Anxiety", "MDD->UKB_DEP", 
                "Anxiety->MDD", "Anxiety->UKB_DEP"),
  Beta = c(primary$Primary_Beta[1], all_results[[1]]$Beta_ratio[1],
           primary$Primary_Beta[2], all_results[[2]]$Beta_ratio[1]),
  SE = c(NA, all_results[[1]]$SE[1], NA, all_results[[2]]$SE[1]),
  P = c(primary$Primary_P[1], all_results[[1]]$P[1],
        primary$Primary_P[2], all_results[[2]]$P[1]),
  N_SNPs = c(198, all_results[[1]]$N_SNPs[1],
             57, all_results[[2]]$N_SNPs[1]),
  stringsAsFactors = FALSE
)

print(comparison, digits = 4)

# Save
all_res <- do.call(rbind, all_results)
write.csv(all_res, file.path(OUT, "ukb_external_validation_ratio.csv"), row.names = FALSE)
write.csv(comparison, file.path(OUT, "ukb_validation_comparison.csv"), row.names = FALSE)

cat("\nResults saved to:", OUT, "\n")
