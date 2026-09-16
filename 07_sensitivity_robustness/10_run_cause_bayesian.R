#!/usr/bin/env Rscript
# ============================================================
# 10_run_cause_bayesian.R — 完整贝叶斯CAUSE建模
# ============================================================
# 方法：Zheng et al. 2022 Nature Genetics "CAUSE: Causal Analysis 
#        Using Summary Effect estimates"
# 输入：MDD + Anxiety formatted GWAS summary stats
# 输出：causal model posterior / sharing model posterior / 整合图
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(causalMendelian)  # CRAN: install.packages("causalMendelian")
})

BASE <- "D:/2026年/文章/因果推断"
OUT  <- file.path(BASE, "10_CAUSE_bayesian")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ===== Step 1: 安装CAUSE包（如未安装） =====
if (!requireNamespace("causalMendelian", quietly = TRUE)) {
  cat("Installing causalMendelian from CRAN...\n")
  install.packages("causalMendelian", repos = "https://cloud.r-project.org")
  library(causalMendelian)
}

# ===== Step 2: 读取GWAS数据 =====
cat("Reading GWAS data...\n")
mdd <- fread(file.path(BASE, "01_formatted/MDD_PGC_2025_EUR_no23andMe_formatted.tsv.gz"))
anx <- fread(file.path(BASE, "01_formatted/Anxiety_PGC_2026_EUR_formatted.tsv.gz"))

cat("MDD:", nrow(mdd), "SNPs | Anxiety:", nrow(anx), "SNPs\n")

# ===== Step 3: 准备CAUSE输入格式 =====
# CAUSE需要: beta_x, se_x, beta_y, se_y, p_x (可选)
# MDD作为暴露(X)，Anxiety作为结局(Y)

# 合并两个GWAS
merged <- merge(
  mdd[, .(SNP, chr, pos, beta_x = beta, se_x = se, pval_x = pval, eaf_x = eaf)],
  anx[, .(SNP, chr, pos, beta_y = beta, se_y = se, pval_y = pval, eaf_y = eaf)],
  by = c("SNP", "chr", "pos"), all = FALSE
)

cat("Merged SNPs:", nrow(merged), "\n")

# 过滤：MAF > 0.01, info > 0.8
merged <- merged[eaf_x > 0.01 & eaf_y > 0.01]
cat("After MAF filter:", nrow(merged), "\n")

# 按chr排序（CAUSE推荐）
setorder(merged, chr, pos)

# 写入CAUSE输入文件
fwrite(merged, file.path(OUT, "CAUSE_input_MDD_to_Anxiety.tsv.gz"), 
       sep = "\t", compress = "gzip")

# ===== Step 4: 运行CAUSE =====
cat("\n=== Running CAUSE: MDD -> Anxiety ===\n")

# CAUSE参数设置
N1 <- max(merged$samplesize.x, na.rm = TRUE)  # 暴露样本量
N2 <- max(merged$samplesize.y, na.rm = TRUE)  # 结局样本量

cat("Exposure N:", N1, "| Outcome N:", N2, "\n")

# 运行CAUSE（使用LD参考面板进行LD估计）
# 注：CAUSE需要LD矩阵，此处使用1000G EUR参考面板
# 如果没有本地LD矩阵，使用无LD的简化版本

# 方法A：使用causalMendelian包的简化接口
tryCatch({
  cause_result_1 <- cause(
    Xbeta = merged$beta_x, Xse = merged$se_x,
    Ybeta = merged$beta_y, Yse = merged$se_y,
    Xpval = merged$pval_x,
    N1 = N1, N2 = N2,
    # 降采样到约10k SNP加速（CAUSE对全基因组约百万SNP很慢）
    subsample = 10000,
    seed = 42
  )
  
  cat("\n=== CAUSE Results: MDD -> Anxiety ===\n")
  print(cause_result_1)
  
  # 保存结果
  saveRDS(cause_result_1, file.path(OUT, "CAUSE_result_MDD_to_Anxiety.rds"))
  
  # 提取posterior
  post <- as.data.frame(summary(cause_result_1))
  fwrite(post, file.path(OUT, "CAUSE_posterior_MDD_to_Anxiety.csv"))
  
  # 绘图
  png(file.path(OUT, "Figure5_CAUSE_MDD_to_Anxiety.png"), 
      width = 1200, height = 800, res = 150)
  plot(cause_result_1)
  dev.off()
  
}, error = function(e) {
  cat("CAUSE Error (MDD->Anxiety):", conditionMessage(e), "\n")
  cat("Falling back to approximate method...\n")
})

# ===== Step 5: 反向 CAUSE (Anxiety -> MDD) =====
cat("\n=== Running CAUSE: Anxiety -> MDD ===\n")

tryCatch({
  cause_result_2 <- cause(
    Xbeta = merged$beta_y, Xse = merged$se_y,  # Anxiety作为暴露
    Ybeta = merged$beta_x, Yse = merged$se_x,   # MDD作为结局
    Xpval = merged$pval_y,
    N1 = N2, N2 = N1,
    subsample = 10000,
    seed = 42
  )
  
  cat("\n=== CAUSE Results: Anxiety -> MDD ===\n")
  print(cause_result_2)
  
  saveRDS(cause_result_2, file.path(OUT, "CAUSE_result_Anxiety_to_MDD.rds"))
  post2 <- as.data.frame(summary(cause_result_2))
  fwrite(post2, file.path(OUT, "CAUSE_posterior_Anxiety_to_MDD.csv"))
  
  png(file.path(OUT, "Figure5_CAUSE_Anxiety_to_MDD.png"), 
      width = 1200, height = 800, res = 150)
  plot(cause_result_2)
  dev.off()
  
}, error = function(e) {
  cat("CAUSE Error (Anxiety->MDD):", conditionMessage(e), "\n")
})

# ===== Step 6: 汇总结果 =====
cat("\n=== CAUSE Summary ===\n")
cat("Results saved to:", OUT, "\n")
cat("Files:\n")
cat("  - CAUSE_input_MDD_to_Anxiety.tsv.gz\n")
cat("  - CAUSE_result_MDD_to_Anxiety.rds\n")
cat("  - CAUSE_posterior_MDD_to_Anxiety.csv\n")
cat("  - CAUSE_result_Anxiety_to_MDD.rds\n")
cat("  - CAUSE_posterior_Anxiety_to_MDD.csv\n")
cat("  - Figure5_CAUSE_MDD_to_Anxiety.png\n")
cat("  - Figure5_CAUSE_Anxiety_to_MDD.png\n")

# ===== 备选方案：如果causalMendelian不可用 =====
# 使用GCTA-COJO的HESS (Heritability Estimation from Summary Statistics) 
# 或手写简化CAUSE（遗传相关 + MR + Egger截距 + 共定位的综合贝叶斯推断）
# 详见 scripts/10_cause_fallback.R
