#!/usr/bin/env Rscript
# 使用真实LD clumping后的SNP运行MVMR（调整BMI+教育+吸烟）
# 策略：直接传向量/矩阵给 ivw_mvmr() 和 mvmr()，绕过 format_mvmr()

suppressPackageStartupMessages({
  library(data.table)
  library(MVMR)
})

BASE <- "D:/2026年/文章/因果推断"
OUT  <- file.path(BASE, "06_MVMR_realLD")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ===== 1. 读取真实LD clumping后的SNP =====
cat("=== Reading real LD clumped SNPs ===\n")

clumped_mdd <- fread(file.path(BASE, "02_instruments_realLD/MDD_realLD_clumped.tsv"))$SNP
clumped_anx <- fread(file.path(BASE, "02_instruments_realLD/Anxiety_realLD_clumped.tsv"))$SNP

cat("MDD clumped SNPs:", length(clumped_mdd), "\n")
cat("Anxiety clumped SNPs:", length(clumped_anx), "\n")

# ===== 2. 读取GWAS数据 =====
cat("\n=== Reading GWAS data ===\n")

mdd_formatted <- fread(file.path(BASE, "01_formatted/MDD_PGC_2025_EUR_no23andMe_formatted.tsv.gz"))
anx_formatted <- fread(file.path(BASE, "01_formatted/Anxiety_PGC_2026_EUR_formatted.tsv.gz"))
bmi     <- fread(file.path(BASE, "06_MVMR/ieu-a-2_BMI_formatted.tsv.gz"))
edu     <- fread(file.path(BASE, "06_MVMR/ieu-a-80_Education_formatted.tsv.gz"))
smoking <- fread(file.path(BASE, "06_MVMR/ieu-b-4877_Smoking_formatted.tsv.gz"))

# ===== 3. MDD -> Anxiety (MVMR) =====
cat("\n=== MDD -> Anxiety MVMR (real LD) ===\n")

# 提取暴露（MDD）
mdd_exp <- mdd_formatted[SNP %in% clumped_mdd, .(SNP, beta.exposure = beta, se.exposure = se)]

# 提取结局（Anxiety）
anx_out <- anx_formatted[SNP %in% clumped_mdd, .(SNP, beta.outcome = beta, se.outcome = se)]

# 提取协变量
bmi_cov     <- bmi[SNP %in% clumped_mdd, .(SNP, beta.bmi = beta, se.bmi = se)]
edu_cov     <- edu[SNP %in% clumped_mdd, .(SNP, beta.edu = beta, se.edu = se)]
smoking_cov <- smoking[SNP %in% clumped_mdd, .(SNP, beta.smoking = beta, se.smoking = se)]

# 按SNP合并
dat <- merge(mdd_exp, anx_out, by = "SNP", all = FALSE)
dat <- merge(dat, bmi_cov, by = "SNP", all = FALSE)
dat <- merge(dat, edu_cov, by = "SNP", all = FALSE)
dat <- merge(dat, smoking_cov, by = "SNP", all = FALSE)

cat("Complete cases for MVMR (MDD->Anxiety):", nrow(dat), "\n")

# 运行MVMR
if (nrow(dat) >= 5) {
  # 构造MVMR输入
  # r1 = 暴露的beta向量
  # r2 = 协变量的beta矩阵 (nSNP x nCovariates)
  # r1se = 暴露的se向量
  # r2se = 协变量的se矩阵 (nSNP x nCovariates)
  
  r1    <- dat$beta.exposure
  r2    <- as.matrix(dat[, .(beta.bmi, beta.edu, beta.smoking)])
  r1se  <- dat$se.exposure
  r2se  <- as.matrix(dat[, .(se.bmi, se.edu, se.smoking)])
  
  # 计算条件F统计量（手动计算，因为strength_mvmr参数搞不定）
  cat("Calculating conditional F-statistics manually...\n")
  # 简化计算：F = (R2/(1-R2)) * ((n-k-1)/k)
  # 这里我们用 MVMR 包的 strength_mvmr 函数，如果它还不工作，我们就报告 NA
  condF <- tryCatch({
    strength_mvmr(r1, r2, r1se, r2se)$condF
  }, error = function(e) {
    cat("WARNING: strength_mvmr() failed, skipping condF calculation.\n")
    return(NA)
  })
  cat("Conditional F-statistics (MDD):", condF, "\n")
  write.csv(data.frame(exposure="MDD", condF=condF), 
            file.path(OUT, "MVMR_MDD_Anxiety_condF.csv"), row.names=FALSE)
  
  # MVMR-IVW
  cat("Running MVMR-IVW...\n")
  mv_res <- ivw_mvmr(r1, r2, r1se, r2se)
  cat("MVMR-IVW (MDD->Anxiety): OR =", exp(mv_res$bhat), "\n")
  write.csv(data.frame(method="IVW", OR=exp(mv_res$bhat), se=mv_res$se, pval=mv_res$pval),
            file.path(OUT, "MVMR_MDD_Anxiety_result.csv"), row.names=FALSE)
  
  # MVMR-Egger
  cat("Running MVMR-Egger...\n")
  mv_egger <- mvmr(r1, r2, r1se, r2se, method="egger")
  cat("MVMR-Egger (MDD->Anxiety): OR =", exp(mv_egger$bhat), "\n")
  write.csv(data.frame(method="Egger", OR=exp(mv_egger$bhat), se=mv_egger$se, pval=mv_egger$pval),
            file.path(OUT, "MVMR_MDD_Anxiety_egger.csv"), row.names=FALSE)
} else {
  cat("WARNING: Not enough complete cases for MVMR (MDD->Anxiety)\n")
}

# ===== 4. Anxiety -> MDD (MVMR) =====
cat("\n=== Anxiety -> MDD MVMR (real LD) ===\n")

# 提取暴露（Anxiety）
anx_exp <- anx_formatted[SNP %in% clumped_anx, .(SNP, beta.exposure = beta, se.exposure = se)]

# 提取结局（MDD）
mdd_out <- mdd_formatted[SNP %in% clumped_anx, .(SNP, beta.outcome = beta, se.outcome = se)]

# 提取协变量
bmi_cov2     <- bmi[SNP %in% clumped_anx, .(SNP, beta.bmi = beta, se.bmi = se)]
edu_cov2     <- edu[SNP %in% clumped_anx, .(SNP, beta.edu = beta, se.edu = se)]
smoking_cov2 <- smoking[SNP %in% clumped_anx, .(SNP, beta.smoking = beta, se.smoking = se)]

# 按SNP合并
dat2 <- merge(anx_exp, mdd_out, by = "SNP", all = FALSE)
dat2 <- merge(dat2, bmi_cov2, by = "SNP", all = FALSE)
dat2 <- merge(dat2, edu_cov2, by = "SNP", all = FALSE)
dat2 <- merge(dat2, smoking_cov2, by = "SNP", all = FALSE)

cat("Complete cases for MVMR (Anxiety->MDD):", nrow(dat2), "\n")

# 运行MVMR
if (nrow(dat2) >= 5) {
  # 构造MVMR输入
  r1    <- dat2$beta.exposure
  r2    <- as.matrix(dat2[, .(beta.bmi, beta.edu, beta.smoking)])
  r1se  <- dat2$se.exposure
  r2se  <- as.matrix(dat2[, .(se.bmi, se.edu, se.smoking)])
  
  # 计算条件F统计量
  cat("Calculating conditional F-statistics manually...\n")
  condF2 <- tryCatch({
    strength_mvmr(r1, r2, r1se, r2se)$condF
  }, error = function(e) {
    cat("WARNING: strength_mvmr() failed, skipping condF calculation.\n")
    return(NA)
  })
  cat("Conditional F-statistics (Anxiety):", condF2, "\n")
  write.csv(data.frame(exposure="Anxiety", condF=condF2), 
            file.path(OUT, "MVMR_Anxiety_MDD_condF.csv"), row.names=FALSE)
  
  # MVMR-IVW
  cat("Running MVMR-IVW...\n")
  mv_res2 <- ivw_mvmr(r1, r2, r1se, r2se)
  cat("MVMR-IVW (Anxiety->MDD): OR =", exp(mv_res2$bhat), "\n")
  write.csv(data.frame(method="IVW", OR=exp(mv_res2$bhat), se=mv_res2$se, pval=mv_res2$pval),
            file.path(OUT, "MVMR_Anxiety_MDD_result.csv"), row.names=FALSE)
  
  # MVMR-Egger
  cat("Running MVMR-Egger...\n")
  mv_egger2 <- mvmr(r1, r2, r1se, r2se, method="egger")
  cat("MVMR-Egger (Anxiety->MDD): OR =", exp(mv_egger2$bhat), "\n")
  write.csv(data.frame(method="Egger", OR=exp(mv_egger2$bhat), se=mv_egger2$se, pval=mv_egger2$pval),
            file.path(OUT, "MVMR_Anxiety_MDD_egger.csv"), row.names=FALSE)
} else {
  cat("WARNING: Not enough complete cases for MVMR (Anxiety->MDD)\n")
}

cat("\n=== Done ===\n")
cat("Results saved to:", OUT, "\n")
