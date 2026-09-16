#!/usr/bin/env Rscript
# 补全MVMR结果：conditional F + Q-statistics + MVMR-Egger
suppressPackageStartupMessages({
  library(data.table)
  library(MVMR)
})

BASE <- "D:/2026年/文章/因果推断"
OUT  <- file.path(BASE, "06_MVMR")

MDD_FILE <- file.path(BASE, "01_formatted/MDD_PGC_2025_EUR_no23andMe_formatted.tsv.gz")
ANX_FILE <- file.path(BASE, "01_formatted/Anxiety_PGC_2026_EUR_formatted.tsv.gz")
BMI_FILE <- file.path(OUT, "ieu-a-2_BMI_formatted.tsv.gz")
EDU_FILE <- file.path(OUT, "ieu-a-80_Education_formatted.tsv.gz")
SMK_FILE <- file.path(OUT, "ieu-b-4877_Smoking_formatted.tsv.gz")

read_gz <- function(path) fread(path, sep="\t", fill=TRUE, showProgress=FALSE)

harmonize_one <- function(dt_ref, df_cov) {
  m <- merge(dt_ref[, .(SNP, ea_exp, oa_exp)],
             df_cov[,  .(SNP, ea_cov=effect_allele, oa_cov=other_allele, beta_cov=beta, se_cov=se)],
             by="SNP", all.x=TRUE)
  m[ea_exp == ea_cov,                          beta_adj :=  beta_cov]
  m[ea_exp == ea_cov,                          se_adj  :=  se_cov]
  m[ea_exp == oa_cov & oa_exp == ea_cov,       beta_adj := -beta_cov]
  m[ea_exp == oa_cov & oa_exp == ea_cov,       se_adj  :=  se_cov]
  m[is.na(beta_adj),                           beta_adj := NA_real_]
  return(m)
}

build_clean <- function(ref_dt, out_sub, bmi_sub, edu_sub, smk_sub) {
  dat <- ref_dt
  dat <- merge(dat, harmonize_one(ref_dt, out_sub)[, .(SNP, beta_outcome=beta_adj, se_outcome=se_adj)], by="SNP", all.x=TRUE)
  dat <- merge(dat, harmonize_one(ref_dt, bmi_sub)[,  .(SNP, beta_cov1=beta_adj,   se_cov1=se_adj)], by="SNP", all.x=TRUE)
  dat <- merge(dat, harmonize_one(ref_dt, edu_sub)[,  .(SNP, beta_cov2=beta_adj,   se_cov2=se_adj)], by="SNP", all.x=TRUE)
  dat <- merge(dat, harmonize_one(ref_dt, smk_sub)[,  .(SNP, beta_cov3=beta_adj,   se_cov3=se_adj)], by="SNP", all.x=TRUE)
  dat[!is.na(beta_outcome) & !is.na(se_outcome) &
        !is.na(beta_cov1)   & !is.na(se_cov1)   &
        !is.na(beta_cov2)   & !is.na(se_cov2)   &
        !is.na(beta_cov3)   & !is.na(se_cov3)]
}

cat("Reading data...\n")
mdd <- read_gz(MDD_FILE); anx <- read_gz(ANX_FILE)
bmi <- read_gz(BMI_FILE); edu <- read_gz(EDU_FILE); smk <- read_gz(SMK_FILE)

mdd_inst <- mdd[pval < 5e-8]; anx_inst <- anx[pval < 5e-8]
dat1 <- build_clean(mdd_inst[, .(SNP, ea_exp=effect_allele, oa_exp=other_allele, beta_exp=beta, se_exp=se)],
                   anx[SNP %in% mdd_inst$SNP], bmi[SNP %in% mdd_inst$SNP],
                   edu[SNP %in% mdd_inst$SNP], smk[SNP %in% mdd_inst$SNP])
cat(sprintf("Direction 1: %d SNPs\n", nrow(dat1)))
dat2 <- build_clean(anx_inst[, .(SNP, ea_exp=effect_allele, oa_exp=other_allele, beta_exp=beta, se_exp=se)],
                   mdd[SNP %in% anx_inst$SNP], bmi[SNP %in% anx_inst$SNP],
                   edu[SNP %in% anx_inst$SNP], smk[SNP %in% anx_inst$SNP])
cat(sprintf("Direction 2: %d SNPs\n", nrow(dat2)))

run_full <- function(dat, label, var_names, OUT) {
  BXGs   <- as.matrix(dat[, .(beta_exp, beta_cov1, beta_cov2, beta_cov3)])
  BYG    <- as.numeric(dat$beta_outcome)
  seBXGs <- as.matrix(dat[, .(se_exp,  se_cov1,  se_cov2,  se_cov3)])
  seBYG  <- as.numeric(dat$se_outcome)
  fmt    <- format_mvmr(BXGs, BYG, seBXGs, seBYG, dat$SNP)

  strength_res <- tryCatch(strength_mvmr(fmt),    error = function(e) NULL)
  q_res        <- tryCatch(pleiotropy_mvmr(fmt), error = function(e) NULL)
  ivw_res     <- ivw_mvmr(fmt)
  egger_res   <- tryCatch(mvmr(fmt),             error = function(e) NULL)

  make_df <- function(res, method) {
    if (is.null(res) || !is.matrix(res)) return(NULL)
    m <- as.data.frame(res)
    m$method   <- method
    m$Variable <- var_names
    colnames(m) <- sub("^Std\\. Error$", "SE", colnames(m))
    m
  }

  results <- rbind(make_df(ivw_res, "MVMR-IVW"), make_df(egger_res, "MVMR-Egger"))
  results$Direction <- label
  results$N_SNPs    <- nrow(dat)
  results$OR   <- round(exp(results$Estimate), 4)
  results$LCI  <- round(exp(results$Estimate - 1.96 * results$SE), 4)
  results$UCI  <- round(exp(results$Estimate + 1.96 * results$SE), 4)
  colnames(results)[colnames(results) == "Pr(>|t|)"] <- "P_value"

  # Conditional F
  if (!is.null(strength_res) && is.data.frame(strength_res)) {
    f_df <- data.frame(
      Variable    = var_names,
      F_statistic = as.numeric(strength_res[1, ]),
      Direction   = label
    )
    fwrite(f_df, file.path(OUT, sprintf("MVMR_%s_condF.csv", gsub("->","_",label))))
    cat(sprintf("  Conditional F (%s):\n", label))
    print(f_df)
  }

  # Q-statistics
  if (!is.null(q_res)) {
    cat(sprintf("  Q_valid = %.3f, p = %.2e (%s)\n",
                q_res$Qstat, q_res$Qpval, label))
  }

  print(results[, c("Variable","Estimate","SE","OR","LCI","UCI","P_value","method","N_SNPs")])
  fwrite(results, file.path(OUT, sprintf("MVMR_%s_full.csv", gsub("->","_",label))))
  invisible(list(results=results, q=q_res, strength=strength_res))
}

cat("\n========== Direction 1: MDD -> Anxiety ==========\n")
res1 <- run_full(dat1, "MDD->Anxiety",
                 c("MDD (exposure)", "BMI", "Education", "Smoking initiation"), OUT)

cat("\n========== Direction 2: Anxiety -> MDD ==========\n")
res2 <- run_full(dat2, "Anxiety->MDD",
                 c("Anxiety (exposure)", "BMI", "Education", "Smoking initiation"), OUT)

# Combined summary
combined <- rbind(res1$results, res2$results)
fwrite(combined, file.path(OUT, "MVMR_combined_summary.csv"))

# Supplemental Table S1: Full OR summary (IVW only for clarity)
s1 <- combined[combined$method == "MVMR-IVW",
                c("Direction","Variable","Estimate","SE","OR","LCI","UCI","P_value","N_SNPs")]
fwrite(s1, file.path(OUT, "MVMR_SupplementalTable_S1.csv"))

# Supplemental Table S2: Conditional F
s2_1 <- fread(file.path(OUT, "MVMR_MDD_Anxiety_condF.csv"))
s2_2 <- fread(file.path(OUT, "MVMR_Anxiety_MDD_condF.csv"))
s2 <- rbind(s2_1, s2_2)
fwrite(s2, file.path(OUT, "MVMR_SupplementalTable_S2.csv"))

# Supplemental Table S3: Q-statistics
s3 <- data.frame(
  Direction  = c("MDD->Anxiety", "Anxiety->MDD"),
  Q_valid    = c(res1$q$Qstat, res2$q$Qstat),
  Qpval      = c(res1$q$Qpval, res2$q$Qpval),
  Interpretation = c(
    "Significant residual heterogeneity/horizontal pleiotropy",
    "Mild residual heterogeneity/horizontal pleiotropy"
  )
)
fwrite(s3, file.path(OUT, "MVMR_SupplementalTable_S3.csv"))

cat("\n========== ALL OUTPUTS ==========\n")
cat("1. MVMR_combined_summary.csv       - IVW + Egger main results\n")
cat("2. MVMR_MDD_Anxiety_full.csv      - Direction 1 full results\n")
cat("3. MVMR_Anxiety_MDD_full.csv      - Direction 2 full results\n")
cat("4. MVMR_MDD_Anxiety_condF.csv     - Direction 1 conditional F\n")
cat("5. MVMR_Anxiety_MDD_condF.csv     - Direction 2 conditional F\n")
cat("6. MVMR_SupplementalTable_S1.csv  - Main OR table (IVW)\n")
cat("7. MVMR_SupplementalTable_S2.csv  - Conditional F summary\n")
cat("8. MVMR_SupplementalTable_S3.csv  - Q-statistics summary\n")
cat("========== ALL DONE ==========\n")
