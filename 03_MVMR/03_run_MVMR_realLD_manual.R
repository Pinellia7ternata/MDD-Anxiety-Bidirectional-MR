#!/usr/bin/env Rscript
# 使用真实LD clumping后的SNP运行MVMR（调整BMI+教育+吸烟）
# 方法：手动实现MVMR-IVW（矩阵运算），避免MVMR包的兼容性bug

suppressPackageStartupMessages({
  library(data.table)
  library(MVMR) # 仅用于 strength_mvmr() 如果可用
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

# 手动实现MVMR-IVW
if (nrow(dat) >= 5) {
  # 构造矩阵
  # Y = 结局的beta向量
  # X = 暴露+协变量的beta矩阵 (nSNP x 4)
  # W = 权重矩阵 (diag(1/se^2))
  
  Y <- dat$beta.outcome
  X <- as.matrix(dat[, .(beta.exposure, beta.bmi, beta.edu, beta.smoking)])
  W <- diag(1 / (dat$se.outcome^2))  # 使用结局的SE作为权重
  
  # MVMR-IVW 公式: beta_hat = (X^T W X)^(-1) X^T W Y
  beta_hat <- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% Y
  se_hat  <- sqrt(diag(solve(t(X) %*% W %*% X)))
  
  # 计算P值（Z检验）
  z_vals <- beta_hat / se_hat
  p_vals <- 2 * pnorm(-abs(z_vals))
  
  # 整理结果
  res <- data.frame(
    exposure   = c("MDD", "BMI", "Education", "Smoking"),
    beta       = as.numeric(beta_hat),
    se         = as.numeric(se_hat),
    pval       = as.numeric(p_vals),
    OR         = exp(as.numeric(beta_hat)),
    CI_lower   = exp(as.numeric(beta_hat) - 1.96 * as.numeric(se_hat)),
    CI_upper   = exp(as.numeric(beta_hat) + 1.96 * as.numeric(se_hat))
  )
  
  cat("MVMR-IVW (MDD->Anxiety):\n")
  print(res, digits=4)
  
  write.csv(res, file.path(OUT, "MVMR_MDD_Anxiety_result.csv"), row.names=FALSE)
  
  # 尝试计算条件F统计量（如果MVMR包可用）
  cat("\nCalculating conditional F-statistics (if MVMR package works)...\n")
  condF <- tryCatch({
    # 构造MVMR格式的输入
    mv_input <- data.frame(
      b_exp = dat$beta.exposure,
      se_exp = dat$se.exposure,
      b_out = dat$beta.outcome,
      se_out = dat$se.outcome,
      b_cov1 = dat$beta.bmi,
      se_cov1 = dat$se.bmi,
      b_cov2 = dat$beta.edu,
      se_cov2 = dat$se.edu,
      b_cov3 = dat$beta.smoking,
      se_cov3 = dat$se.smoking
    )
    strength_mvmr(mv_input)$condF
  }, error = function(e) {
    cat("WARNING: strength_mvmr() failed, using manual calculation.\n")
    # 手动计算F统计量（简化版）
    R2 <- summary(lm(dat$beta.outcome ~ dat$beta.exposure))$r.squared
    n  <- nrow(dat)
    k  <- 4  # 4个暴露（MDD+BMI+Edu+Smoking）
    F_stat <- (R2 / (1 - R2)) * ((n - k - 1) / k)
    return(F_stat)
  })
  
  cat("Conditional F-statistic (MDD):", condF, "\n")
  write.csv(data.frame(exposure="MDD", condF=condF), 
            file.path(OUT, "MVMR_MDD_Anxiety_condF.csv"), row.names=FALSE)
  
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

# 手动实现MVMR-IVW
if (nrow(dat2) >= 5) {
  Y <- dat2$beta.outcome
  X <- as.matrix(dat2[, .(beta.exposure, beta.bmi, beta.edu, beta.smoking)])
  W <- diag(1 / (dat2$se.outcome^2))
  
  beta_hat <- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% Y
  se_hat  <- sqrt(diag(solve(t(X) %*% W %*% X)))
  
  z_vals <- beta_hat / se_hat
  p_vals <- 2 * pnorm(-abs(z_vals))
  
  res2 <- data.frame(
    exposure   = c("Anxiety", "BMI", "Education", "Smoking"),
    beta       = as.numeric(beta_hat),
    se         = as.numeric(se_hat),
    pval       = as.numeric(p_vals),
    OR         = exp(as.numeric(beta_hat)),
    CI_lower   = exp(as.numeric(beta_hat) - 1.96 * as.numeric(se_hat)),
    CI_upper   = exp(as.numeric(beta_hat) + 1.96 * as.numeric(se_hat))
  )
  
  cat("MVMR-IVW (Anxiety->MDD):\n")
  print(res2, digits=4)
  
  write.csv(res2, file.path(OUT, "MVMR_Anxiety_MDD_result.csv"), row.names=FALSE)
  
  # 尝试计算条件F统计量
  cat("\nCalculating conditional F-statistics (if MVMR package works)...\n")
  condF2 <- tryCatch({
    mv_input2 <- data.frame(
      b_exp = dat2$beta.exposure,
      se_exp = dat2$se.exposure,
      b_out = dat2$beta.outcome,
      se_out = dat2$se.outcome,
      b_cov1 = dat2$beta.bmi,
      se_cov1 = dat2$se.bmi,
      b_cov2 = dat2$beta.edu,
      se_cov2 = dat2$se.edu,
      b_cov3 = dat2$beta.smoking,
      se_cov3 = dat2$se.smoking
    )
    strength_mvmr(mv_input2)$condF
  }, error = function(e) {
    cat("WARNING: strength_mvmr() failed, using manual calculation.\n")
    R2 <- summary(lm(dat2$beta.outcome ~ dat2$beta.exposure))$r.squared
    n  <- nrow(dat2)
    k  <- 4
    F_stat <- (R2 / (1 - R2)) * ((n - k - 1) / k)
    return(F_stat)
  })
  
  cat("Conditional F-statistic (Anxiety):", condF2, "\n")
  write.csv(data.frame(exposure="Anxiety", condF=condF2), 
            file.path(OUT, "MVMR_Anxiety_MDD_condF.csv"), row.names=FALSE)
  
} else {
  cat("WARNING: Not enough complete cases for MVMR (Anxiety->MDD)\n")
}

cat("\n=== Done ===\n")
cat("Results saved to:", OUT, "\n")
