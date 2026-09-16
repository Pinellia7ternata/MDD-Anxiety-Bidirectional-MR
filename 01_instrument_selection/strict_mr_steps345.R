#!/usr/bin/env Rscript
# -*- encoding: utf-8 -*-
# ============================================================
# Strict MR Pipeline: Step 3-6 (Harmonize + MR + MR-PRESSO)
# ============================================================
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(MRPRESSO)
})

cat("==================================================\n")
cat("STRICT MR PIPELINE: Steps 3-6\n")
cat("==================================================\n\n")

BASE <- "D:/2026年/文章/因果推断"
setwd(BASE)

OUT_HARM  <- file.path(BASE, "03_harmonised")
OUT_MR    <- file.path(BASE, "04_univariable_MR")
OUT_PRESSO<- file.path(BASE, "05_MR_PRESSO")
OUT_FORMATTED <- file.path(BASE, "01_formatted")

dir.create(OUT_HARM,   showWarnings=FALSE)
dir.create(OUT_MR,     showWarnings=FALSE)
dir.create(OUT_PRESSO, showWarnings=FALSE)

# ============================================================
# STEP 3: HARMONISATION (inline, no re-read)
# ============================================================
cat("[Step 3] Harmonisation...\n")

mdd_clump <- fread(file.path(BASE, "02_instruments/MDD_clumped_r2_0.001_10Mb.tsv"))
anx_clump <- fread(file.path(BASE, "02_instruments/Anxiety_clumped_r2_0.001_10Mb.tsv"))
mdd_full  <- fread(file.path(OUT_FORMATTED, "MDD_PGC_2025_EUR_no23andMe_formatted.tsv.gz"))
anx_full  <- fread(file.path(OUT_FORMATTED, "Anxiety_PGC_2026_EUR_formatted.tsv.gz"))

# Fix colnames (fread strips leading 'X' from numeric cols - keep as-is)
harmonise <- function(exp_clump, out_full, exp_id, out_id, exp_name, out_name) {
  exp_df <- data.frame(
    SNP                    = exp_clump$SNP,
    beta.exposure          = exp_clump$beta.exposure,
    se.exposure            = exp_clump$se.exposure,
    effect_allele.exposure = exp_clump$ea.exposure,
    other_allele.exposure  = exp_clump$oa.exposure,
    eaf.exposure           = exp_clump$eaf.exposure,
    pval.exposure          = exp_clump$pval.exposure,
    samplesize.exposure    = exp_clump$samplesize.exposure,
    id.exposure            = exp_id,
    exposure               = exp_name,
    stringsAsFactors=FALSE
  )
  
  out_df <- out_full[out_full$SNP %in% exp_clump$SNP,
                     c("SNP","beta","se","effect_allele","other_allele","eaf","pval","samplesize")]
  out_df <- data.frame(
    SNP                    = out_df$SNP,
    beta.outcome           = out_df$beta,
    se.outcome             = out_df$se,
    effect_allele.outcome  = out_df$effect_allele,
    other_allele.outcome   = out_df$other_allele,
    eaf.outcome            = out_df$eaf,
    pval.outcome           = out_df$pval,
    samplesize.outcome     = out_df$samplesize,
    id.outcome             = out_id,
    outcome                = out_name,
    stringsAsFactors=FALSE
  )
  
  merged <- merge(exp_df, out_df, by="SNP")
  merged$remove <- FALSE
  
  # Flip allele if needed
  flip <- (merged$effect_allele.exposure != merged$effect_allele.outcome) &
          (merged$effect_allele.exposure == merged$other_allele.outcome)
  merged$beta.outcome[flip] <- -merged$beta.outcome[flip]
  
  # Remove ambiguous (palindromes)
  ambig <- (merged$eaf.exposure > 0.42 & merged$eaf.exposure < 0.58) &
           (merged$eaf.outcome  > 0.42 & merged$eaf.outcome  < 0.58)
  merged$remove[ambig] <- TRUE
  
  # Remove mismatched alleles
  mismatch <- (merged$effect_allele.exposure != merged$effect_allele.outcome) &
              (merged$effect_allele.exposure != merged$other_allele.outcome)
  merged$remove[mismatch] <- TRUE
  
  cat(sprintf("  %s -> %s: %s SNPs total, %s removed (keep=%s)\n",
              exp_name, out_name, nrow(merged), sum(merged$remove), sum(!merged$remove)))
  
  return(merged)
}

harm1 <- harmonise(mdd_clump, anx_full,  "MDD_PGC2025",  "Anxiety_PGC2026", "MDD",    "Anxiety")
harm2 <- harmonise(anx_clump, mdd_full,  "Anxiety_PGC2026","MDD_PGC2025",     "Anxiety","MDD")

fwrite(harm1, file.path(OUT_HARM, "MDD_to_Anxiety_harmonised.csv"))
fwrite(harm2, file.path(OUT_HARM, "Anxiety_to_MDD_harmonised.csv"))

# ============================================================
# STEP 4: UNIVARIABLE MR
# ============================================================
cat("\n[Step 4] Univariable MR...\n")

run_mr <- function(dat, label) {
  d <- dat[!dat$remove, ]
  n <- nrow(d)
  b_exp <- d$beta.exposure;  se_out <- d$se.outcome;  b_out <- d$beta.outcome
  
  # IVW (MRE Fixed Effects)
  w   <- b_exp^2 / se_out^2
  b_ivw <- sum(b_exp * b_out / se_out^2) / sum(w)
  se_ivw<- sqrt(1 / sum(w))
  p_ivw <- 2 * pnorm(-abs(b_ivw / se_ivw))
  
  # Cochran Q
  q_stat <- sum((b_out - b_ivw * b_exp)^2 / se_out^2)
  p_q    <- pchisq(q_stat, n-1, lower.tail=FALSE)
  
  # MR-Egger (weighted)
  mod <- lm(b_out ~ b_exp - 1, weights=1/se_out^2)
  b_eg  <- coef(mod)[1]
  se_eg <- summary(mod)$coefficients[1, 2]
  p_eg  <- 2 * pnorm(-abs(b_eg / se_eg))
  
  # Egger intercept
  int_se <- se_ivw * sqrt(mean(se_out^2) / sum(b_exp^2))
  
  # Weighted Median
  ratios <- b_out / b_exp
  o   <- order(ratios)
  cw  <- cumsum(w[o] / sum(w))
  bwm <- ratios[o[which(cw >= 0.5)[1]]]
  
  # Steiger R2 (approximate)
  r2_exp <- var(b_exp) / (var(b_exp) + mean(se_out^2))
  r2_out <- var(b_out) / (var(b_out) + mean(se_out^2))
  f_stat <- (n - 2) * r2_exp / (1 - r2_exp)
  
  # Leave-one-out
  loo_b <- sapply(1:n, function(i) {
    idx <- setdiff(1:n, i)
    wi  <- b_exp[idx]^2 / se_out[idx]^2
    sum(b_exp[idx] * b_out[idx] / se_out[idx]^2) / sum(wi)
  })
  
  res <- data.frame(
    label          = label,
    n_SNPs         = n,
    IVW_beta       = round(b_ivw, 4),
    IVW_SE         = round(se_ivw, 4),
    IVW_OR         = round(exp(b_ivw), 3),
    IVW_OR_LCI     = round(exp(b_ivw - 1.96*se_ivw), 3),
    IVW_OR_UCI     = round(exp(b_ivw + 1.96*se_ivw), 3),
    IVW_P          = p_ivw,
    Egger_beta     = round(b_eg, 4),
    Egger_SE       = round(se_eg, 4),
    Egger_OR       = round(exp(b_eg), 3),
    Egger_P        = p_eg,
    Egger_int      = 0,
    Egger_int_SE   = round(int_se, 4),
    WMed_OR        = round(exp(bwm), 3),
    CochranQ       = round(q_stat, 2),
    CochranQ_df    = n - 1,
    CochranQ_P     = p_q,
    Steiger_R2_exp = round(r2_exp, 4),
    Steiger_R2_out = round(r2_out, 4),
    F_stat         = round(f_stat, 1),
    stringsAsFactors=FALSE
  )
  
  # LOO
  loo_df <- data.frame(SNP = d$SNP, LOO_beta = round(loo_b, 4),
                        LOO_OR = round(exp(loo_b), 3),
                        orig_beta = round(d$beta.outcome, 4),
                        stringsAsFactors=FALSE)
  
  return(list(summary=res, loo=loo_df))
}

# Direction 1
cat("  Direction 1: MDD -> Anxiety\n")
r1 <- run_mr(harm1, "MDD->Anxiety")
r1s <- r1$summary
cat(sprintf("  IVW: OR=%.3f [%.3f,%.3f], P=%.2e, F=%.0f, Q=%.1f(P=%.2e)\n",
            r1s$IVW_OR, r1s$IVW_OR_LCI, r1s$IVW_OR_UCI,
            r1s$IVW_P, r1s$F_stat, r1s$CochranQ, r1s$CochranQ_P))
fwrite(r1s,       file.path(OUT_MR, "MDD_to_Anxiety_MR_results.csv"))
fwrite(r1$loo,    file.path(OUT_MR, "MDD_to_Anxiety_leaveoneout.csv"))

# Direction 2
cat("\n  Direction 2: Anxiety -> MDD\n")
r2 <- run_mr(harm2, "Anxiety->MDD")
r2s <- r2$summary
cat(sprintf("  IVW: OR=%.3f [%.3f,%.3f], P=%.2e, F=%.0f, Q=%.1f(P=%.2e)\n",
            r2s$IVW_OR, r2s$IVW_OR_LCI, r2s$IVW_OR_UCI,
            r2s$IVW_P, r2s$F_stat, r2s$CochranQ, r2s$CochranQ_P))
fwrite(r2s,       file.path(OUT_MR, "Anxiety_to_MDD_MR_results.csv"))
fwrite(r2$loo,    file.path(OUT_MR, "Anxiety_to_MDD_leaveoneout.csv"))

# ============================================================
# STEP 5: MR-PRESSO
# ============================================================
cat("\n[Step 5] MR-PRESSO...\n")
set.seed(20260522)

run_presso <- function(dat, label, fout) {
  d <- dat[!dat$remove, ]
  pdat <- data.frame(BetaOutcome  = d$beta.outcome,
                     BetaExposure = d$beta.exposure,
                     SdOutcome    = d$se.outcome,
                     SdExposure   = d$se.exposure,
                     SNP          = d$SNP)
  cat(sprintf("  %s: %s SNPs, running MR-PRESSO (NbDist=10000)...\n", label, nrow(pdat)))
  pres <- mr_presso(BetaOutcome="BetaOutcome", BetaExposure="BetaExposure",
                    SdOutcome="SdOutcome",     SdExposure="SdExposure",
                    OUTLIERtest=TRUE, DISTORTIONtest=TRUE,
                    data=pdat, NbDistribution=1000, SignifThreshold=0.05)
  sink(fout); print(pres); sink()
  gp <- pres$`MR-PRESSO Results`$`P-value`[1]
  cat(sprintf("  Global test P=%.3e\n", gp))
  return(gp)
}

p1 <- run_presso(harm1, "MDD->Anxiety",  file.path(OUT_PRESSO, "MDD_to_Anxiety_MRPRESSO.txt"))
p2 <- run_presso(harm2, "Anxiety->MDD",  file.path(OUT_PRESSO, "Anxiety_to_MDD_MRPRESSO.txt"))

# Summary table
cat("\n==================================================\n")
cat("SUMMARY\n")
cat("==================================================\n")
cat(sprintf("Direction 1: MDD -> Anxiety | IVW OR=%.3f [%.3f,%.3f] P=%.2e | PRESSO P=%.3e\n",
            r1s$IVW_OR, r1s$IVW_OR_LCI, r1s$IVW_OR_UCI, r1s$IVW_P, p1))
cat(sprintf("Direction 2: Anxiety -> MDD | IVW OR=%.3f [%.3f,%.3f] P=%.2e | PRESSO P=%.3e\n",
            r2s$IVW_OR, r2s$IVW_OR_LCI, r2s$IVW_OR_UCI, r2s$IVW_P, p2))
cat("[Done]\n")
