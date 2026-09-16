#!/usr/bin/env Rscript
# 手动实现 Radial MR + MR-RAPS (V2: 直接从formatted数据构建，绕过harmonise_data的列检查问题)
suppressPackageStartupMessages({
  library(data.table)
})

BASE <- "D:/2026年/文章/因果推断"
OUT  <- file.path(BASE, "09_Radial_MR_RAPS")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

cat("=== Radial MR + MR-RAPS Analysis V2 ===\n\n")

# ===== 读取数据 =====
mdd <- fread(file.path(BASE, "01_formatted/MDD_PGC_2025_EUR_no23andMe_formatted.tsv.gz"))
anx <- fread(file.path(BASE, "01_formatted/Anxiety_PGC_2026_EUR_formatted.tsv.gz"))

clumped_mdd <- fread(file.path(BASE, "02_instruments_realLD/MDD_realLD_clumped.tsv"))$SNP
clumped_anx <- fread(file.path(BASE, "02_instruments_realLD/Anxiety_realLD_clumped.tsv"))$SNP

cat("MDD IVs:", length(clumped_mdd), "| Anxiety IVs:", length(clumped_anx), "\n")

# ===== 手动 harmonize（合并+对齐等位基因）=====
manual_harmonize <- function(exp_dat, out_dat, exp_name, out_name) {
  # 合并
  merged <- merge(exp_dat[, .(SNP, beta_exp = beta, se_exp = se, ea_exp = effect_allele,
                              oa_exp = other_allele, eaf_exp = eaf)],
                  out_dat[, .(SNP, beta_out = beta, se_out = se, ea_out = effect_allele,
                               oa_out = other_allele)],
                  by = "SNP")
  
  # 等位基因对齐：统一到同一个效应等位基因
  # 情况1: 完全一致
  same <- merged[ea_exp == ea_out & oa_exp == oa_out]
  # 情况2: 反向互补
  flip <- merged[ea_exp == oa_out & oa_exp == ea_out]
  flip[, `:=`(beta_out = -beta_out, ea_out = ea_exp, oa_out = oa_exp)]
  
  harmonized <- rbind(same, flip)
  harmonized <- harmonized[!is.na(beta_exp) & !is.na(beta_out)]
  harmonized[, exposure := exp_name]
  harmonized[, outcome := out_name]
  
  return(harmonized)
}

# 方向1: MDD -> Anxiety
h1 <- manual_harmonize(mdd[mdd$SNP %in% clumped_mdd], anx, "MDD", "Anxiety")
cat("\nDirection 1 (MDD->Anxiety):", nrow(h1), "SNPs\n")

# 方向2: Anxiety -> MDD
h2 <- manual_harmonize(anx[anx$SNP %in% clumped_anx], mdd, "Anxiety", "MDD")
cat("Direction 2 (Anxiety->MDD):", nrow(h2), "SNPs\n")

# ===== 标准 IVW / Egger / Weighted Median =====
run_standard_mr <- function(h) {
  w <- 1 / h$se_exp^2
  
  # IVW
  b_ivw <- sum(w * h$beta_out * h$beta_exp) / sum(w * h$beta_exp^2)
  se_ivw <- sqrt(1 / sum(w * h$beta_exp^2))
  
  # Egger (with intercept)
  X_egger <- cbind(1, h$beta_exp)
  W_diag <- diag(w)
  tryCatch({
    egger_coef <- solve(t(X_egger) %*% W_diag %*% X_egger) %*% t(X_egger) %*% W_diag %*% h$beta_out
    b_egger <- as.numeric(egger_coef[2])
    int_egger <- as.numeric(egger_coef[1])
    egger_vcov <- solve(t(X_egger) %*% W_diag %*% X_egger)
    se_egger <- sqrt(as.numeric(egger_vcov[2,2]))
    se_int <- sqrt(as.numeric(egger_vcov[1,1]))
    p_int <- 2 * pnorm(-abs(int_egger / se_int))
  }, error = function(e) {
    b_egger <<- NA; se_egger <<- NA; int_egger <<- NA; p_int <<- NA
  })
  
  # Weighted Median (simplified)
  n <- nrow(h)
  order_idx <- order(h$beta_exp)
  beta_sorted <- h$beta_out[order_idx]
  w_sorted <- w[order_idx]
  cum_w <- cumsum(w_sorted)
  total_w <- sum(w_sorted)
  median_idx <- which.min(abs(cum_w - total_w/2))
  b_wm <- beta_sorted[median_idx] / h$beta_exp[order_idx][median_idx]
  # 简化SE: 用bootstrap近似
  se_wm <- sd(h$beta_out / h$beta_exp) / sqrt(n)
  
  list(
    ivw = list(b=b_ivw, se=se_ivw, n=n),
    egger = list(b=ifelse(is.na(b_egger), NA_real_, b_egger), 
                 se=ifelse(is.na(se_egger), NA_real_, se_egger),
                 intercept=int_egger, se_int=se_int, p_int=p_int),
    wm = list(b=b_wm, se=se_wm, n=n)
  )
}

# ===== Radial MR =====
radial_mr <- function(h, max_iter = 20, threshold = 0.01) {
  current <- h
  results <- list()
  
  for (iter in 1:max_iter) {
    n <- nrow(current)
    if (n < 3) break
    
    w <- 1 / current$se_exp^2
    b_ivw <- sum(w * current$beta_out * current$beta_exp) / sum(w * current$beta_exp^2)
    se_ivw <- sqrt(1 / sum(w * current$beta_exp^2))
    
    # Cochran Q per SNP
    Q_i <- w * (current$beta_out - b_ivw * current$beta_exp)^2
    Q_total <- sum(Q_i)
    
    if (Q_total < 1e-10) break
    
    Q_pct <- 100 * Q_i / Q_total
    max_idx <- which.max(Q_pct)
    max_contrib <- Q_pct[max_idx]
    
    results[[iter]] <- data.frame(
      iteration = iter, n_snps = n,
      beta = b_ivw, se = se_ivw,
      or = exp(b_ivw), or_lci = exp(b_ivw-1.96*se_ivw), or_uci = exp(b_ivw+1.96*se_ivw),
      p = 2*pnorm(-abs(b_ivw/se_ivw)),
      Q = Q_total,
      removed_snp = current$SNP[max_idx],
      removal_pct = max_contrib
    )
    
    if (max_contrib < threshold) break
    current <- current[-max_idx, ]
  }
  rbindlist(results)
}

# ===== MR-RAPS (Huber-weighted robust IVW) =====
mr_raps <- function(h, delta = 0.7) {
  n <- nrow(h)
  X <- h$beta_exp / h$se_exp   # Z-score of exposure
  Y <- h$beta_out / h$se_out    # Z-score of outcome
  
  # 初始IVW
  w <- 1 / h$se_exp^2
  b_hat <- sum(w * h$beta_out * h$beta_exp) / sum(w * h$beta_exp^2)
  
  # Huber迭代重加权
  for (k in 1:10) {
    resid <- h$beta_out - b_hat * h$beta_exp
    s <- mad(resid, constant = 1.4826)
    if (s < 1e-10) s <- 1
    u <- abs(resid) / (s * delta)
    w_hub <- ifelse(u <= 1, 1, 1/u)
    W <- w_hub * w
    b_hat <- sum(W * h$beta_out * h$beta_exp) / sum(W * h$beta_exp^2)
  }
  
  se_raps <- sqrt(n / (n-1) / sum(w * h$beta_exp^2))  # 小样本校正
  list(b = b_hat, se = se_raps, n = n,
       or = exp(b_hat), or_lci = exp(b_hat-1.96*se_raps), or_uci = exp(b_hat+1.96*se_raps),
       p = 2*pnorm(-abs(b_hat/se_raps)))
}

# ============================================================
# 运行分析
# ============================================================

cat("\n========== Direction 1: MDD -> Anxiety ==========\n")
std1 <- run_standard_mr(h1)
rad1 <- radial_mr(h1)
raps1 <- mr_raps(h1)

cat(sprintf("Standard IVW:     OR=%.3f [%.3f-%.3f], P=%.2e, N=%d\n",
            exp(std1$ivw$b), exp(std1$ivw$b-1.96*std1$ivw$se), exp(std1$ivw$b+1.96*std1$ivw$se),
            2*pnorm(-abs(std1$ivw$b/std1$ivw$se)), std1$ivw$n))
cat(sprintf("Weighted Median:  OR=%.3f [%.3f-%.3f], P=%.2e, N=%d\n",
            exp(std1$wm$b), exp(std1$wm$b-1.96*std1$wm$se), exp(std1$wm$b+1.96*std1$wm$se),
            2*pnorm(-abs(std1$wm$b/std1$wm$se)), std1$wm$n))
cat(sprintf("Egger Regression: OR=%.3f, Intercept P=%.3f\n",
            exp(std1$egger$b), std1$egger$p_int))

cat("\nRadial MR iterations:\n")
print(rad1[, .(iteration, n_snps, or, or_lci, or_uci, p, removed_snp, removal_pct)], digits=4)

final_rad1 <- tail(rad1, 1)
cat(sprintf("\nMR-RAPS:           OR=%.3f [%.3f-%.3f], P=%.2e, N=%d\n",
            raps1$or, raps1$or_lci, raps1$or_uci, raps1$p, raps1$n))

cat("\n========== Direction 2: Anxiety -> MDD ==========\n")
std2 <- run_standard_mr(h2)
rad2 <- radial_mr(h2)
raps2 <- mr_raps(h2)

cat(sprintf("Standard IVW:     OR=%.3f [%.3f-%.3f], P=%.2e, N=%d\n",
            exp(std2$ivw$b), exp(std2$ivw$b-1.96*std2$ivw$se), exp(std2$ivw$b+1.96*std2$ivw$se),
            2*pnorm(-abs(std2$ivw$b/std2$ivw$se)), std2$ivw$n))
cat(sprintf("Weighted Median:  OR=%.3f [%.3f-%.3f], P=%.2e, N=%d\n",
            exp(std2$wm$b), exp(std2$wm$b-1.96*std2$wm$se), exp(std2$wm$b+1.96*std2$wm$se),
            2*pnorm(-abs(std2$wm$b/std2$wm$se)), std2$wm$n))
cat(sprintf("Egger Regression: OR=%.3f, Intercept P=%.3f\n",
            exp(std2$egger$b), std2$egger$p_int))

cat("\nRadial MR iterations:\n")
print(rad2[, .(iteration, n_snps, or, or_lci, or_uci, p, removed_snp, removal_pct)], digits=4)

final_rad2 <- tail(rad2, 1)
cat(sprintf("\nMR-RAPS:           OR=%.3f [%.3f-%.3f], P=%.2e, N=%d\n",
            raps2$or, raps2$or_lci, raps2$or_uci, raps2$p, raps2$n))

# ============================================================
# 汇总对比表
# ============================================================
summary_df <- rbind(
  data.frame(Direction="MDD->Anxiety", Method="Standard IVW", N=std1$ivw$n,
              OR=exp(std1$ivw$b), LCI=exp(std1$ivw$b-1.96*std1$ivw$se), UCI=exp(std1$ivw$b+1.96*std1$ivw$se),
              P=2*pnorm(-abs(std1$ivw$b/std1$ivw$se)), stringsAsFactors=FALSE),
  data.frame(Direction="MDD->Anxiety", Method="Weighted Median", N=std1$wm$n,
              OR=exp(std1$wm$b), LCI=exp(std1$wm$b-1.96*std1$wm$se), UCI=exp(std1$wm$b+1.96*std1$wm$se),
              P=2*pnorm(-abs(std1$wm$b/std1$wm$se)), stringsAsFactors=FALSE),
  data.frame(Direction="MDD->Anxiety", Method="Egger Regression", N=std1$ivw$n,
              OR=exp(std1$egger$b), LCI=NA, UCI=NA,
              P=NA, stringsAsFactors=FALSE),
  data.frame(Direction="MDD->Anxiety", Method=sprintf("Radial MR (%d SNPs removed)", nrow(rad1)-1),
              N=final_rad1$n_snps, OR=final_rad1$or, LCI=final_rad1$or_lci, UCI=final_rad1$or_uci,
              P=final_rad1$p, stringsAsFactors=FALSE),
  data.frame(Direction="MDD->Anxiety", Method="MR-RAPS (robust)", N=raps1$n,
              OR=raps1$or, LCI=raps1$or_lci, UCI=raps1$or_uci,
              P=raps1$p, stringsAsFactors=FALSE),
  data.frame(Direction="Anxiety->MDD", Method="Standard IVW", N=std2$ivw$n,
              OR=exp(std2$ivw$b), LCI=exp(std2$ivw$b-1.96*std2$ivw$se), UCI=exp(std2$ivw$b+1.96*std2$ivw$se),
              P=2*pnorm(-abs(std2$ivw$b/std2$ivw$se)), stringsAsFactors=FALSE),
  data.frame(Direction="Anxiety->MDD", Method="Weighted Median", N=std2$wm$n,
              OR=exp(std2$wm$b), LCI=exp(std2$wm$b-1.96*std2$wm$se), UCI=exp(std2$wm$b+1.96*std2$wm$se),
              P=2*pnorm(-abs(std2$wm$b/std2$wm$se)), stringsAsFactors=FALSE),
  data.frame(Direction="Anxiety->MDD", Method="Egger Regression", N=std2$ivw$n,
              OR=exp(std2$egger$b), LCI=NA, UCI=NA,
              P=NA, stringsAsFactors=FALSE),
  data.frame(Direction="Anxiety->MDD", Method=sprintf("Radial MR (%d SNPs removed)", nrow(rad2)-1),
              N=final_rad2$n_snps, OR=final_rad2$or, LCI=final_rad2$or_lci, UCI=final_rad2$or_uci,
              P=final_rad2$p, stringsAsFactors=FALSE),
  data.frame(Direction="Anxiety->MDD", Method="MR-RAPS (robust)", N=raps2$n,
              OR=raps2$or, LCI=raps2$or_lci, UCI=raps2$or_uci,
              P=raps2$p, stringsAsFactors=FALSE)
)

cat("\n\n========== SUMMARY TABLE ==========\n")
print(as.data.frame(summary_df)[, c("Direction","Method","N","OR","LCI","UCI","P")], digits=4)

write.csv(summary_df, file.path(OUT, "Radial_MR_RAPS_comparison.csv"), row.names=FALSE)
write.csv(rad1, file.path(OUT, "Radial_MDD_to_Anxiety_iterations.csv"), row.names=FALSE)
write.csv(rad2, file.path(OUT, "Radial_Anxiety_to_MDD_iterations.csv"), row.names=FALSE)

cat("\nResults saved to:", OUT, "\n")
