#!/usr/bin/env Rscript
# ============================================================
# 本地MVMR分析：MDD/Anxiety + BMI/Education/Smoking
# 使用预先格式化好的本地文件，无需OpenGWAS token
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(MVMR)
  library(TwoSampleMR)
})

# ---------------- 配置 ----------------
BASE  <- "D:/2026年/文章/因果推断"
OUT   <- file.path(BASE, "06_MVMR")
dir.create(OUT, showWarnings = FALSE)

# 格式化好的文件（TSV.gz）
MDD_FILE  <- file.path(BASE, "01_formatted/MDD_PGC_2025_EUR_no23andMe_formatted.tsv.gz")
ANX_FILE  <- file.path(BASE, "01_formatted/Anxiety_PGC_2026_EUR_formatted.tsv.gz")
BMI_FILE  <- file.path(OUT, "ieu-a-2_BMI_formatted.tsv.gz")
EDU_FILE  <- file.path(OUT, "ieu-a-80_Education_formatted.tsv.gz")
SMK_FILE  <- file.path(OUT, "ieu-b-4877_Smoking_formatted.tsv.gz")

# ---------------- 工具函数 ----------------
read_gz <- function(path) {
  # data.table::fread 原生支持 .gz 文件（Windows/Linux 均可用）
  fread(path, sep = "\t", fill = TRUE, showProgress = TRUE)
}

harmonize_one <- function(dt_ref, df_cov, label) {
  # dt_ref: 参考表，含 SNP, ea_exp, oa_exp, beta_exp, se_exp
  # df_cov:  协变量 data.table，含 SNP, effect_allele, other_allele, beta, se
  # 返回：对齐后的 data.table，beta_cov 已对齐到 ea_exp
  
  # 先把协变量等位基因列名统一
  m <- merge(
    dt_ref[, .(SNP, ea_exp, oa_exp)],
    df_cov[,  .(SNP,
                 ea_cov = effect_allele,
                 oa_cov = other_allele,
                 beta_cov = beta,
                 se_cov   = se)],
    by = "SNP", all.x = TRUE
  )
  
  # 情况1：effect_allele 完全一致 → 直接用
  m[ea_exp == ea_cov, beta_adj := beta_cov]
  m[ea_exp == ea_cov, se_adj  := se_cov]
  
  # 情况2：effect_allele 互换 → beta 取反
  m[ea_exp == oa_cov & oa_exp == ea_cov, beta_adj := -beta_cov]
  m[ea_exp == oa_cov & oa_exp == ea_cov, se_adj  := se_cov]
  
  # 情况3：无法对齐 → NA
  m[is.na(beta_adj), beta_adj := NA_real_]
  
  n_align  <- sum(m$ea_exp == m$ea_cov, na.rm = TRUE)
  n_flip  <- sum(m$ea_exp == m$oa_cov & m$oa_exp == m$ea_cov, na.rm = TRUE)
  n_miss  <- sum(is.na(m$beta_adj))
  cat(sprintf("  [%s] Total=%d, Aligned=%d, Flipped=%d, Missing=%d\n",
              label, nrow(m), n_align, n_flip, n_miss))
  return(m)
}

# ============================================================
# 方向1：MDD → Anxiety，调整 BMI + Education + Smoking
# ============================================================
cat("========== Direction 1: MDD -> Anxiety (MVMR) ==========\n")

cat("[1/5] Reading MDD instruments (p < 5e-8)...\n")
mdd <- read_gz(MDD_FILE)
mdd_inst <- mdd[pval < 5e-8]
cat(sprintf("  MDD instruments: %d SNPs\n", nrow(mdd_inst)))

cat("[2/5] Reading Anxiety outcome...\n")
anx <- read_gz(ANX_FILE)
anx_sub <- anx[SNP %in% mdd_inst$SNP]
cat(sprintf("  Anxiety overlap: %d SNPs\n", nrow(anx_sub)))

cat("[3/5] Reading covariates (BMI, Education, Smoking)...\n")
bmi <- read_gz(BMI_FILE)
edu <- read_gz(EDU_FILE)
smk <- read_gz(SMK_FILE)

bmi_sub <- bmi[SNP %in% mdd_inst$SNP]
edu_sub <- edu[SNP %in% mdd_inst$SNP]
smk_sub <- smk[SNP %in% mdd_inst$SNP]
cat(sprintf("  BMI overlap: %d, Education: %d, Smoking: %d\n",
            nrow(bmi_sub), nrow(edu_sub), nrow(smk_sub)))

cat("[4/5] Harmonizing alleles...\n")
# 先构建暴露（MDD）的参考表
ref <- mdd_inst[, .(SNP, ea_exp = effect_allele, oa_exp = other_allele,
                     beta_exp = beta, se_exp = se)]

# Harmonize Anxiety
h_anx <- harmonize_one(ref, anx_sub, "Anxiety")
# Harmonize BMI
h_bmi <- harmonize_one(ref, bmi_sub, "BMI")
# Harmonize Education
h_edu <- harmonize_one(ref, edu_sub, "Education")
# Harmonize Smoking
h_smk <- harmonize_one(ref, smk_sub, "Smoking")

cat("[5/5] Merging and running MVMR...\n")
dat <- ref
dat <- merge(dat, h_anx[, .(SNP, beta_outcome = beta_adj, se_outcome = se_adj)], by = "SNP", all.x = TRUE)
dat <- merge(dat, h_bmi[, .(SNP, beta_cov1   = beta_adj, se_cov1   = se_adj)], by = "SNP", all.x = TRUE)
dat <- merge(dat, h_edu[, .(SNP, beta_cov2   = beta_adj, se_cov2   = se_adj)], by = "SNP", all.x = TRUE)
dat <- merge(dat, h_smk[, .(SNP, beta_cov3   = beta_adj, se_cov3   = se_adj)], by = "SNP", all.x = TRUE)

# 去掉任何含NA的行
dat_clean <- dat[!is.na(beta_outcome) & !is.na(se_outcome) &
                 !is.na(beta_cov1)    & !is.na(se_cov1)    &
                 !is.na(beta_cov2)    & !is.na(se_cov2)    &
                 !is.na(beta_cov3)    & !is.na(se_cov3)]
cat(sprintf("  Complete cases: %d / %d SNPs\n", nrow(dat_clean), nrow(dat)))

if (nrow(dat_clean) < 10) {
  cat("  ERROR: Too few complete SNPs for MVMR!\n")
  quit(status = 1)
}

# 准备MVMR包所需格式
mvdat <- data.frame(
  SNP       = dat_clean$SNP,
  beta.exposure = as.numeric(dat_clean$beta_exp),
  beta.outcome  = as.numeric(dat_clean$beta_outcome),
  se.exposure   = as.numeric(dat_clean$se_exp),
  se.outcome    = as.numeric(dat_clean$se_outcome),
  beta.cov1     = as.numeric(dat_clean$beta_cov1),
  se.cov1       = as.numeric(dat_clean$se_cov1),
  beta.cov2     = as.numeric(dat_clean$beta_cov2),
  se.cov2       = as.numeric(dat_clean$se_cov2),
  beta.cov3     = as.numeric(dat_clean$beta_cov3),
  se.cov3       = as.numeric(dat_clean$se_cov3),
  exposure = "MDD",
  outcome  = "Anxiety"
)
rownames(mvdat) <- mvdat$SNP

cat("  Running mvMR_estimation...\n")
mv_res <- tryCatch({
  MVMR::mvMR_estimation(mvdat,
                         beta.outcome  = "beta.outcome",
                         beta.exposure = "beta.exposure",
                         se.outcome    = "se.outcome",
                         se.exposure   = "se.exposure",
                         beta.covariate = c("beta.cov1","beta.cov2","beta.cov3"),
                         se.covariate   = c("se.cov1","se.cov2","se.cov3"))
}, error = function(e) {
  cat("  MVMR error:", e$message, "\n")
  NULL
})

if (!is.null(mv_res)) {
  cat("\n--- MVMR Results (MDD -> Anxiety) ---\n")
  print(mv_res)
  
  # 计算OR和CI
  if (!is.null(mv_res$Estimate)) {
    or_tab <- data.frame(
      Variable = rownames(mv_res$Estimate),
      Beta     = mv_res$Estimate[,1],
      SE       = mv_res$Estimate[,2],
      OR       = round(exp(mv_res$Estimate[,1]), 3),
      CI_lower = round(exp(mv_res$Estimate[,1] - 1.96 * mv_res$Estimate[,2]), 3),
      CI_upper = round(exp(mv_res$Estimate[,1] + 1.96 * mv_res$Estimate[,2]), 3),
      P_value  = mv_res$Estimate[,4]
    )
    print(or_tab)
    fwrite(or_tab, file.path(OUT, "MVMR_MDD_Anxiety_result.csv"))
    cat(sprintf("  Saved to %s\n", file.path(OUT, "MVMR_MDD_Anxiety_result.csv")))
  }
}

# ============================================================
# 方向2：Anxiety → MDD，调整 BMI + Education + Smoking
# ============================================================
cat("\n========== Direction 2: Anxiety -> MDD (MVMR) ==========\n")

cat("[1/5] Extracting Anxiety instruments (p < 5e-8)...\n")
anx_inst <- anx[pval < 5e-8]
cat(sprintf("  Anxiety instruments: %d SNPs\n", nrow(anx_inst)))

cat("[2/5] Reading MDD outcome...\n")
mdd_sub <- mdd[SNP %in% anx_inst$SNP]
cat(sprintf("  MDD overlap: %d SNPs\n", nrow(mdd_sub)))

cat("[3/5] Reading covariates for Anxiety instruments...\n")
bmi_sub2 <- bmi[SNP %in% anx_inst$SNP]
edu_sub2 <- edu[SNP %in% anx_inst$SNP]
smk_sub2 <- smk[SNP %in% anx_inst$SNP]
cat(sprintf("  BMI: %d, Education: %d, Smoking: %d\n",
            nrow(bmi_sub2), nrow(edu_sub2), nrow(smk_sub2)))

cat("[4/5] Harmonizing...\n")
ref2 <- anx_inst[, .(SNP, ea_exp = effect_allele, oa_exp = other_allele,
                     beta_exp = beta, se_exp = se)]
h_mdd2 <- harmonize_one(ref2, mdd_sub, "MDD")
h_bmi2 <- harmonize_one(ref2, bmi_sub2, "BMI")
h_edu2 <- harmonize_one(ref2, edu_sub2, "Education")
h_smk2 <- harmonize_one(ref2, smk_sub2, "Smoking")

cat("[5/5] Merging and running MVMR...\n")
dat2 <- ref2
dat2 <- merge(dat2, h_mdd2[, .(SNP, beta_outcome = beta_adj, se_outcome = se_adj)], by = "SNP", all.x = TRUE)
dat2 <- merge(dat2, h_bmi2[, .(SNP, beta_cov1 = beta_adj, se_cov1 = se_adj)],  by = "SNP", all.x = TRUE)
dat2 <- merge(dat2, h_edu2[, .(SNP, beta_cov2 = beta_adj, se_cov2 = se_adj)],  by = "SNP", all.x = TRUE)
dat2 <- merge(dat2, h_smk2[, .(SNP, beta_cov3 = beta_adj, se_cov3 = se_adj)],  by = "SNP", all.x = TRUE)

dat2_clean <- dat2[!is.na(beta_outcome) & !is.na(se_outcome) &
                   !is.na(beta_cov1)    & !is.na(se_cov1)    &
                   !is.na(beta_cov2)    & !is.na(se_cov2)    &
                   !is.na(beta_cov3)    & !is.na(se_cov3)]
cat(sprintf("  Complete cases: %d / %d SNPs\n", nrow(dat2_clean), nrow(dat2)))

if (nrow(dat2_clean) >= 10) {
  mvdat2 <- data.frame(
    SNP            = dat2_clean$SNP,
    beta.exposure  = as.numeric(dat2_clean$beta_exp),
    beta.outcome   = as.numeric(dat2_clean$beta_outcome),
    se.exposure    = as.numeric(dat2_clean$se_exp),
    se.outcome     = as.numeric(dat2_clean$se_outcome),
    beta.cov1      = as.numeric(dat2_clean$beta_cov1),
    se.cov1        = as.numeric(dat2_clean$se_cov1),
    beta.cov2      = as.numeric(dat2_clean$beta_cov2),
    se.cov2        = as.numeric(dat2_clean$se_cov2),
    beta.cov3      = as.numeric(dat2_clean$beta_cov3),
    se.cov3        = as.numeric(dat2_clean$se_cov3),
    exposure = "Anxiety",
    outcome  = "MDD"
  )
  rownames(mvdat2) <- mvdat2$SNP
  
  mv_res2 <- tryCatch({
    MVMR::mvMR_estimation(mvdat2,
                           beta.outcome  = "beta.outcome",
                           beta.exposure = "beta.exposure",
                           se.outcome    = "se.outcome",
                           se.exposure   = "se.exposure",
                           beta.covariate = c("beta.cov1","beta.cov2","beta.cov3"),
                           se.covariate   = c("se.cov1","se.cov2","se.cov3"))
  }, error = function(e) {
    cat("  MVMR error:", e$message, "\n")
    NULL
  })
  
  if (!is.null(mv_res2)) {
    cat("\n--- MVMR Results (Anxiety -> MDD) ---\n")
    print(mv_res2)
    
    if (!is.null(mv_res2$Estimate)) {
      or_tab2 <- data.frame(
        Variable = rownames(mv_res2$Estimate),
        Beta     = mv_res2$Estimate[,1],
        SE       = mv_res2$Estimate[,2],
        OR       = round(exp(mv_res2$Estimate[,1]), 3),
        CI_lower = round(exp(mv_res2$Estimate[,1] - 1.96 * mv_res2$Estimate[,2]), 3),
        CI_upper = round(exp(mv_res2$Estimate[,1] + 1.96 * mv_res2$Estimate[,2]), 3),
        P_value  = mv_res2$Estimate[,4]
      )
      print(or_tab2)
      fwrite(or_tab2, file.path(OUT, "MVMR_Anxiety_MDD_result.csv"))
    }
  }
} else {
  cat("  ERROR: Too few complete SNPs for MVMR!\n")
}

cat("\n========== ALL DONE ==========\n")
