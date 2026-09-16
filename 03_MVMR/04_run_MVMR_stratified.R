#!/usr/bin/env Rscript
# 分层MVMR：不同协变量组合，评估条件F统计量变化
# Model 1: 主暴露 + BMI
# Model 2: 主暴露 + Education
# Model 3: 主暴露 + Smoking
# Model 4: 主暴露 + BMI + Smoking
# Model 5: 完整模型 (BMI + Education + Smoking)

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- "D:/2026年/文章/因果推断"
OUT  <- file.path(BASE, "07_MVMR_stratified")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ===== 读取数据 =====
cat("=== Reading data ===\n")

clumped_mdd <- fread(file.path(BASE, "02_instruments_realLD/MDD_realLD_clumped.tsv"))$SNP
clumped_anx <- fread(file.path(BASE, "02_instruments_realLD/Anxiety_realLD_clumped.tsv"))$SNP

mdd_fmt     <- fread(file.path(BASE, "01_formatted/MDD_PGC_2025_EUR_no23andMe_formatted.tsv.gz"))
anx_fmt     <- fread(file.path(BASE, "01_formatted/Anxiety_PGC_2026_EUR_formatted.tsv.gz"))
bmi         <- fread(file.path(BASE, "06_MVMR/ieu-a-2_BMI_formatted.tsv.gz"))
edu         <- fread(file.path(BASE, "06_MVMR/ieu-a-80_Education_formatted.tsv.gz"))
smoking     <- fread(file.path(BASE, "06_MVMR/ieu-b-4877_Smoking_formatted.tsv.gz"))

cat("MDD SNPs:", length(clumped_mdd), "| Anxiety SNPs:", length(clumped_anx), "\n")

# ===== 定义分层MVMR函数 =====
run_stratified_mvmr <- function(exp_data, out_data, cov_list, direction_label, model_name) {
  # exp_data: data.frame(SNP, beta.exposure, se.exposure)
  # out_data: data.frame(SNP, beta.outcome, se.outcome)
  # cov_list: list of data.frames, each with (SNP, beta.XX, se.XX)
  
  # 合并所有数据
  dat <- merge(exp_data, out_data, by = "SNP", all = FALSE)
  for (cov in cov_list) {
    dat <- merge(dat, cov, by = "SNP", all = FALSE)
  }
  
  n_complete <- nrow(dat)
  n_cov <- length(cov_list)
  
  if (n_complete < 5) {
    cat(sprintf("  %s: SKIP (only %d complete cases)\n", model_name, n_complete))
    return(NULL)
  }
  
  # 构造矩阵：第一列是主暴露，其余是协变量
  beta_cols <- c("beta.exposure", sapply(cov_list, function(x) names(x)[2]))
  X <- as.matrix(dat[, ..beta_cols])
  Y <- dat$beta.outcome
  W <- diag(1 / (dat$se.outcome^2))
  
  # MVMR-IVW
  XtWX <- t(X) %*% W %*% X
  beta_hat <- solve(XtWX) %*% t(X) %*% W %*% Y
  se_hat  <- sqrt(diag(solve(XtWX)))
  
  z_vals  <- beta_hat / se_hat
  p_vals  <- 2 * pnorm(-abs(z_vals))
  
  # 暴露名
  exp_names <- c(names(exp_data)[2],  # "beta.exposure" -> 需要映射
                 sapply(cov_list, function(x) gsub("beta\\.", "", names(x)[2])))
  # 更好的命名
  name_map <- c("beta.exposure" = ifelse(grepl("MDD", direction_label), "MDD", "Anxiety"),
                "beta.bmi" = "BMI", "beta.edu" = "Education", "beta.smoking" = "Smoking")
  exp_names_clean <- name_map[beta_cols]
  
  res <- data.frame(
    Direction   = direction_label,
    Model       = model_name,
    N_SNPs      = n_complete,
    N_covariates = n_cov,
    Exposure    = exp_names_clean,
    Beta        = as.numeric(beta_hat),
    SE          = as.numeric(se_hat),
    P           = as.numeric(p_vals),
    OR          = exp(as.numeric(beta_hat)),
    OR_95CI_LCI = exp(as.numeric(beta_hat) - 1.96 * as.numeric(se_hat)),
    OR_95CI_UCI = exp(as.numeric(beta_hat) + 1.96 * as.numeric(se_hat)),
    stringsAsFactors = FALSE
  )
  
  # 条件F统计量（手动计算：用主暴露对结局的偏F）
  # 方法：从(X'WX)^{-1}的对角元素计算
  # 简化版：用主暴露的t^2作为条件F的近似
  cond_F_main <- (beta_hat[1] / se_hat[1])^2
  
  cat(sprintf("  %s: N=%d, MainExp OR=%.3f [%.3f-%.3f], P=%.2e, CondF=%.2f\n",
              model_name, n_complete, res$OR[1], res$OR_95CI_LCI[1], res$OR_95CI_UCI[1],
              res$P[1], cond_F_main))
  
  return(list(results = res, condF_main = cond_F_main, n_complete = n_complete))
}

# ===== 定义协变量数据提取函数 =====
extract_cov <- function(gwas_df, snp_vec, beta_name, se_name) {
  gwas_df[SNP %in% snp_vec, .(SNP, beta = get(beta_name), se = get(se_name))]
}
setnames(extract_cov(bmi, clumped_mdd, "beta", "se"), c("SNP", "beta.bmi", "se.bmi"))

# ===== MDD -> Anxiety 方向 =====
cat("\n========== MDD -> Anxiety ==========\n")

mdd_exp <- mdd_fmt[SNP %in% clumped_mdd, .(SNP, beta.exposure = beta, se.exposure = se)]
anx_out <- anx_fmt[SNP %in% clumped_mdd, .(SNP, beta.outcome = beta, se.outcome = se)]

bmi_mdd     <- bmi[SNP %in% clumped_mdd, .(SNP, beta.bmi = beta, se.bmi = se)]
edu_mdd     <- edu[SNP %in% clumped_mdd, .(SNP, beta.edu = beta, se.edu = se)]
smoking_mdd <- smoking[SNP %in% clumped_mdd, .(SNP, beta.smoking = beta, se.smoking = se)]

results_ma <- list()
condF_ma  <- c()

# Model 1: MDD + BMI
r1 <- run_stratified_mvmr(mdd_exp, anx_out, list(bmi_mdd), "MDD->Anxiety", "Model1_MDD+BMI")
if (!is.null(r1)) { results_ma[[length(results_ma)+1]] <- r1$results; condF_ma["Model1_BMI"] <- r1$condF_main }

# Model 2: MDD + Education
r2 <- run_stratified_mvmr(mdd_exp, anx_out, list(edu_mdd), "MDD->Anxiety", "Model2_MDD+Edu")
if (!is.null(r2)) { results_ma[[length(results_ma)+1]] <- r2$results; condF_ma["Model2_Edu"] <- r2$condF_main }

# Model 3: MDD + Smoking
r3 <- run_stratified_mvmr(mdd_exp, anx_out, list(smoking_mdd), "MDD->Anxiety", "Model3_MDD+Smoking")
if (!is.null(r3)) { results_ma[[length(results_ma)+1]] <- r3$results; condF_ma["Model3_Smoking"] <- r3$condF_main }

# Model 4: MDD + BMI + Smoking
r4 <- run_stratified_mvmr(mdd_exp, anx_out, list(bmi_mdd, smoking_mdd), "MDD->Anxiety", "Model4_MDD+BMI+Smoking")
if (!is.null(r4)) { results_ma[[length(results_ma)+1]] <- r4$results; condF_ma["Model4_BMI+Smoke"] <- r4$condF_main }

# Model 5: Full model (MDD + BMI + Education + Smoking)
r5 <- run_stratified_mvmr(mdd_exp, anx_out, list(bmi_mdd, edu_mdd, smoking_mdd), "MDD->Anxiety", "Model5_Full")
if (!is.null(r5)) { results_ma[[length(results_ma)+1]] <- r5$results; condF_ma["Model5_Full"] <- r5$condF_main }

# ===== Anxiety -> MDD 方向 =====
cat("\n========== Anxiety -> MDD ==========\n")

anx_exp <- anx_fmt[SNP %in% clumped_anx, .(SNP, beta.exposure = beta, se.exposure = se)]
mdd_out <- mdd_fmt[SNP %in% clumped_anx, .(SNP, beta.outcome = beta, se.outcome = se)]

bmi_anx     <- bmi[SNP %in% clumped_anx, .(SNP, beta.bmi = beta, se.bmi = se)]
edu_anx     <- edu[SNP %in% clumped_anx, .(SNP, beta.edu = beta, se.edu = se)]
smoking_anx <- smoking[SNP %in% clumped_anx, .(SNP, beta.smoking = beta, se.smoking = se)]

results_am <- list()
condF_am  <- c()

# Model 1: Anxiety + BMI
r1a <- run_stratified_mvmr(anx_exp, mdd_out, list(bmi_anx), "Anxiety->MDD", "Model1_Anxiety+BMI")
if (!is.null(r1a)) { results_am[[length(results_am)+1]] <- r1a$results; condF_am["Model1_BMI"] <- r1a$condF_main }

# Model 2: Anxiety + Education
r2a <- run_stratified_mvmr(anx_exp, mdd_out, list(edu_anx), "Anxiety->MDD", "Model2_Anxiety+Edu")
if (!is.null(r2a)) { results_am[[length(results_am)+1]] <- r2a$results; condF_am["Model2_Edu"] <- r2a$condF_main }

# Model 3: Anxiety + Smoking
r3a <- run_stratified_mvmr(anx_exp, mdd_out, list(smoking_anx), "Anxiety->MDD", "Model3_Anxiety+Smoking")
if (!is.null(r3a)) { results_am[[length(results_am)+1]] <- r3a$results; condF_am["Model3_Smoking"] <- r3a$condF_main }

# Model 4: Anxiety + BMI + Smoking
r4a <- run_stratified_mvmr(anx_exp, mdd_out, list(bmi_anx, smoking_anx), "Anxiety->MDD", "Model4_Anxiety+BMI+Smoking")
if (!is.null(r4a)) { results_am[[length(results_am)+1]] <- r4a$results; condF_am["Model4_BMI+Smoke"] <- r4a$condF_main }

# Model 5: Full model
r5a <- run_stratified_mvmr(anx_exp, mdd_out, list(bmi_anx, edu_anx, smoking_anx), "Anxiety->MDD", "Model5_Full")
if (!is.null(r5a)) { results_am[[length(results_am)+1]] <- r5a$results; condF_am["Model5_Full"] <- r5a$condF_main }

# ===== 合并保存结果 =====
cat("\n=== Saving results ===\n")

all_results <- do.call(rbind, c(results_ma, results_am))
write.csv(all_results, file.path(OUT, "MVMR_stratified_all_models.csv"), row.names=FALSE)

# 条件F汇总
condF_summary <- data.frame(
  Direction = c(rep("MDD->Anxiety", length(condF_ma)), rep("Anxiety->MDD", length(condF_am))),
  Model     = c(names(condF_ma), names(condF_am)),
  Conditional_F = c(condF_ma, condF_am),
  stringsAsFactors = FALSE
)
write.csv(condF_summary, file.path(OUT, "MVMR_stratified_conditional_F.csv"), row.names=FALSE)

cat("\n=== Stratified MVMR Summary ===\n")
cat("\n--- Conditional F Statistics ---\n")
print(condF_summary)

cat("\n--- All Results ---\n")
print(all_results[, .(Direction, Model, Exposure, OR, OR_95CI_LCI, OR_95CI_UCI, P, N_SNPs)])

cat("\nDone! Results saved to:", OUT, "\n")
