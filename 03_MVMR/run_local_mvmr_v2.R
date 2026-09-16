#!/usr/bin/env Rscript
# ============================================================
# 本地MVMR分析：MDD/Anxiety + BMI/Education/Smoking
# 使用预先格式化好的本地文件，无需OpenGWAS token
# 正确调用 MVMR 包：format_mvmr() + mvmr() / ivw_mvmr()
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(MVMR)
})

# ---------------- 配置 ----------------
BASE <- "D:/2026年/文章/因果推断"
OUT  <- file.path(BASE, "06_MVMR")
dir.create(OUT, showWarnings = FALSE)

MDD_FILE <- file.path(BASE, "01_formatted/MDD_PGC_2025_EUR_no23andMe_formatted.tsv.gz")
ANX_FILE <- file.path(BASE, "01_formatted/Anxiety_PGC_2026_EUR_formatted.tsv.gz")
BMI_FILE <- file.path(OUT, "ieu-a-2_BMI_formatted.tsv.gz")
EDU_FILE <- file.path(OUT, "ieu-a-80_Education_formatted.tsv.gz")
SMK_FILE <- file.path(OUT, "ieu-b-4877_Smoking_formatted.tsv.gz")

# ---------------- 工具函数 ----------------
read_gz <- function(path) {
  fread(path, sep = "\t", fill = TRUE, showProgress = TRUE)
}

harmonize_one <- function(dt_ref, df_cov, label) {
  # dt_ref: 参考表，含 SNP, ea_exp, oa_exp
  # df_cov:  协变量 DT，含 SNP, effect_allele, other_allele, beta, se
  # 返回：对齐后的DT，beta_adj/se_adj 已对齐到 ea_exp
  m <- merge(
    dt_ref[, .(SNP, ea_exp, oa_exp)],
    df_cov[, .(SNP,
                 ea_cov = effect_allele,
                 oa_cov = other_allele,
                 beta_cov = beta,
                 se_cov   = se)],
    by = "SNP", all.x = TRUE
  )
  # 情况1：effect_allele 完全一致
  m[ea_exp == ea_cov, beta_adj := beta_cov]
  m[ea_exp == ea_cov, se_adj  := se_cov]
  # 情况2：effect_allele 互换 → beta 取反
  m[ea_exp == oa_cov & oa_exp == ea_cov, beta_adj := -beta_cov]
  m[ea_exp == oa_cov & oa_exp == ea_cov, se_adj  := se_cov]
  # 情况3：无法对齐 → NA
  m[is.na(beta_adj), beta_adj := NA_real_]
  
  cat(sprintf("  [%s] Total=%d, Aligned=%d, Flipped=%d, Missing=%d\n",
              label, nrow(m),
              sum(m$ea_exp == m$ea_cov, na.rm = TRUE),
              sum(m$ea_exp == m$oa_cov & m$oa_exp == m$ea_cov, na.rm = TRUE),
              sum(is.na(m$beta_adj))))
  return(m)
}

format_and_run_mvmr <- function(dat_clean, direction_label, OUT) {
  # dat_clean: data.table，含 beta_exp, se_exp, beta_outcome, se_outcome,
  #                               beta_cov1/2/3, se_cov1/2/3, SNP
  # 返回：结果 data.frame
  
  cat(sprintf("\n  [%s] Formatting data for MVMR...\n", direction_label))
  
  # 构建 format_mvmr 所需矩阵
  # BXGs: n_SNP x 4 矩阵，列 = exposure + 3个协变量
  BXGs  <- as.matrix(dat_clean[, .(beta_exp, beta_cov1, beta_cov2, beta_cov3)])
  BYG   <- as.numeric(dat_clean$beta_outcome)
  seBXGs <- as.matrix(dat_clean[, .(se_exp,  se_cov1,  se_cov2,  se_cov3)])
  seBYG  <- as.numeric(dat_clean$se_outcome)
  RSID   <- dat_clean$SNP
  
  # 检查是否有NA
  na_rows <- unique(c(which(is.na(BXGs), arr.ind = TRUE)[,1],
                     which(is.na(seBXGs), arr.ind = TRUE)[,1],
                     which(is.na(BYG)),
                     which(is.na(seBYG))))
  if (length(na_rows) > 0) {
    cat(sprintf("  WARNING: %d rows still have NA, removing...\n", length(na_rows)))
    keep <- setdiff(seq_len(nrow(dat_clean)), na_rows)
    BXGs   <- BXGs[keep, , drop = FALSE]
    BYG     <- BYG[keep]
    seBXGs  <- seBXGs[keep, , drop = FALSE]
    seBYG   <- seBYG[keep]
    RSID    <- RSID[keep]
  }
  
  cat(sprintf("  Final SNPs for MVMR: %d\n", length(RSID)))
  
  # format_mvmr: BXGs/BYG/seBXGs/seBYG 均为矩阵/向量
  fmt <- tryCatch({
    format_mvmr(BXGs = BXGs, BYG = BYG,
                 seBXGs = seBXGs, seBYG = seBYG,
                 RSID = RSID)
  }, error = function(e) {
    cat("  format_mvmr ERROR:", e$message, "\n")
    NULL
  })
  
  if (is.null(fmt)) return(NULL)
  
  # 跑 MVMR-Egger
  cat("  Running MVMR-Egger...\n")
  res_egger <- tryCatch({
    mvmr(fmt)
  }, error = function(e) {
    cat("  mvmr() ERROR:", e$message, "\n")
    NULL
  })
  
  # 跑 MVMR-IVW
  cat("  Running MVMR-IVW...\n")
  res_ivw <- tryCatch({
    ivw_mvmr(fmt)
  }, error = function(e) {
    cat("  ivw_mvmr() ERROR:", e$message, "\n")
    NULL
  })
  
  # 整理结果
  cat("\n--- MVMR Results ---\n")
  if (!is.null(res_egger)) {
    cat("=== MVMR-Egger ===\n")
    print(res_egger)
  }
  if (!is.null(res_ivw)) {
    cat("\n=== MVMR-IVW (inverse variance weighted) ===\n")
    print(res_ivw)
  }
  
  # 保存
  result_df <- data.frame()
  if (!is.null(res_egger) && is.matrix(res_egger)) {
    dg <- as.data.frame(res_egger)
    dg$method <- "MVMR-Egger"
    result_df <- rbind(result_df, cbind(Variable = rownames(res_egger), dg))
  }
  if (!is.null(res_ivw) && is.matrix(res_ivw)) {
    dg <- as.data.frame(res_ivw)
    dg$method <- "MVMR-IVW"
    result_df <- rbind(result_df, cbind(Variable = rownames(res_ivw), dg))
  }
  
  if (nrow(result_df) > 0) {
    out_csv <- file.path(OUT, sprintf("MVMR_%s_result.csv", gsub("->", "_", gsub(" ", "_", direction_label))))
    fwrite(result_df, out_csv)
    cat(sprintf("\n  Saved: %s\n", out_csv))
  }
  
  return(invisible(result_df))
}

# ============================================================
# 方向1：MDD → Anxiety，调整 BMI + Education + Smoking
# ============================================================
cat("========== Direction 1: MDD -> Anxiety (MVMR) ==========\n")

cat("[1/5] Reading MDD...\n")
mdd <- read_gz(MDD_FILE)
mdd_inst <- mdd[pval < 5e-8]
cat(sprintf("  MDD instruments: %d SNPs\n", nrow(mdd_inst)))

cat("[2/5] Reading Anxiety...\n")
anx <- read_gz(ANX_FILE)
anx_sub <- anx[SNP %in% mdd_inst$SNP]
cat(sprintf("  Anxiety overlap: %d\n", nrow(anx_sub)))

cat("[3/5] Reading covariates...\n")
bmi <- read_gz(BMI_FILE); bmi_sub <- bmi[SNP %in% mdd_inst$SNP]
edu <- read_gz(EDU_FILE); edu_sub <- edu[SNP %in% mdd_inst$SNP]
smk <- read_gz(SMK_FILE); smk_sub <- smk[SNP %in% mdd_inst$SNP]
cat(sprintf("  BMI=%d, Edu=%d, Smk=%d\n", nrow(bmi_sub), nrow(edu_sub), nrow(smk_sub)))

cat("[4/5] Harmonizing...\n")
ref <- mdd_inst[, .(SNP, ea_exp = effect_allele, oa_exp = other_allele,
                     beta_exp = beta, se_exp = se)]
h_anx <- harmonize_one(ref, anx_sub, "Anxiety")
h_bmi <- harmonize_one(ref, bmi_sub, "BMI")
h_edu <- harmonize_one(ref, edu_sub, "Education")
h_smk <- harmonize_one(ref, smk_sub, "Smoking")

cat("[5/5] Merging...\n")
dat <- ref
dat <- merge(dat, h_anx[, .(SNP, beta_outcome = beta_adj, se_outcome = se_adj)], by = "SNP", all.x = TRUE)
dat <- merge(dat, h_bmi[, .(SNP, beta_cov1   = beta_adj, se_cov1   = se_adj)], by = "SNP", all.x = TRUE)
dat <- merge(dat, h_edu[, .(SNP, beta_cov2   = beta_adj, se_cov2   = se_adj)], by = "SNP", all.x = TRUE)
dat <- merge(dat, h_smk[, .(SNP, beta_cov3   = beta_adj, se_cov3   = se_adj)], by = "SNP", all.x = TRUE)

dat_clean <- dat[!is.na(beta_outcome) & !is.na(se_outcome) &
                 !is.na(beta_cov1)    & !is.na(se_cov1)    &
                 !is.na(beta_cov2)    & !is.na(se_cov2)    &
                 !is.na(beta_cov3)    & !is.na(se_cov3)]
cat(sprintf("  Complete cases: %d / %d\n", nrow(dat_clean), nrow(dat)))

if (nrow(dat_clean) >= 10) {
  format_and_run_mvmr(dat_clean, "MDD->Anxiety", OUT)
} else {
  cat("  SKIP: too few complete SNPs\n")
}

# ============================================================
# 方向2：Anxiety → MDD，调整 BMI + Education + Smoking
# ============================================================
cat("\n========== Direction 2: Anxiety -> MDD (MVMR) ==========\n")

cat("[1/5] Extracting Anxiety instruments...\n")
anx_inst <- anx[pval < 5e-8]
cat(sprintf("  Anxiety instruments: %d SNPs\n", nrow(anx_inst)))

cat("[2/5] Reading MDD as outcome...\n")
mdd_sub2 <- mdd[SNP %in% anx_inst$SNP]
cat(sprintf("  MDD overlap: %d\n", nrow(mdd_sub2)))

cat("[3/5] Subsetting covariates...\n")
bmi_sub2 <- bmi[SNP %in% anx_inst$SNP]
edu_sub2 <- edu[SNP %in% anx_inst$SNP]
smk_sub2 <- smk[SNP %in% anx_inst$SNP]
cat(sprintf("  BMI=%d, Edu=%d, Smk=%d\n", nrow(bmi_sub2), nrow(edu_sub2), nrow(smk_sub2)))

cat("[4/5] Harmonizing...\n")
ref2 <- anx_inst[, .(SNP, ea_exp = effect_allele, oa_exp = other_allele,
                      beta_exp = beta, se_exp = se)]
h_mdd2 <- harmonize_one(ref2, mdd_sub2, "MDD")
h_bmi2 <- harmonize_one(ref2, bmi_sub2, "BMI")
h_edu2 <- harmonize_one(ref2, edu_sub2, "Education")
h_smk2 <- harmonize_one(ref2, smk_sub2, "Smoking")

cat("[5/5] Merging...\n")
dat2 <- ref2
dat2 <- merge(dat2, h_mdd2[, .(SNP, beta_outcome = beta_adj, se_outcome = se_adj)], by = "SNP", all.x = TRUE)
dat2 <- merge(dat2, h_bmi2[, .(SNP, beta_cov1   = beta_adj, se_cov1   = se_adj)], by = "SNP", all.x = TRUE)
dat2 <- merge(dat2, h_edu2[, .(SNP, beta_cov2   = beta_adj, se_cov2   = se_adj)], by = "SNP", all.x = TRUE)
dat2 <- merge(dat2, h_smk2[, .(SNP, beta_cov3   = beta_adj, se_cov3   = se_adj)], by = "SNP", all.x = TRUE)

dat2_clean <- dat2[!is.na(beta_outcome) & !is.na(se_outcome) &
                    !is.na(beta_cov1)    & !is.na(se_cov1)    &
                    !is.na(beta_cov2)    & !is.na(se_cov2)    &
                    !is.na(beta_cov3)    & !is.na(se_cov3)]
cat(sprintf("  Complete cases: %d / %d\n", nrow(dat2_clean), nrow(dat2)))

if (nrow(dat2_clean) >= 10) {
  format_and_run_mvmr(dat2_clean, "Anxiety->MDD", OUT)
} else {
  cat("  SKIP: too few complete SNPs\n")
}

cat("\n========== ALL DONE ==========\n")
