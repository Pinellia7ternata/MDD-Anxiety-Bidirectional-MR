#!/usr/bin/env Rscript
# ============================================================
# 10_run_lcv_causal_inference.R - LCV因果推断分析
# ============================================================
# 方法：O'Connor & Price 2020 Nat Genet "Latent Causal Variable"
#       区分遗传相关与因果效应
# 替代CAUSE（因CAUSE需要Rtools编译，LCV在CRAN上可直接安装）
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# 安装LCV包
if (!requireNamespace("LCV", quietly = TRUE)) {
  cat("Installing LCV package from CRAN...\n")
  install.packages("LCV", repos = "https://cloud.r-project.org")
  library(LCV)
}

BASE <- "D:/2026年/文章/因果推断"
OUT  <- file.path(BASE, "10_LCV_causal")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ===== Step 1: 读取GWAS数据 =====
cat("Step 1: Reading GWAS data...\n")

mdd <- fread(file.path(BASE, "01_formatted/MDD_PGC_2025_EUR_no23andMe_formatted.tsv.gz"))
anx <- fread(file.path(BASE, "01_formatted/Anxiety_PGC_2026_EUR_formatted.tsv.gz"))

cat("  MDD:", nrow(mdd), "SNPs\n")
cat("  Anxiety:", nrow(anx), "SNPs\n")

# ===== Step 2: 合并GWAS =====
cat("\nStep 2: Merging GWAS data...\n")

merged <- merge(
  mdd[, .(SNP, chr, pos, beta_x = beta, se_x = se, pval_x = pval, eaf_x = eaf, n_x = samplesize)],
  anx[, .(SNP, chr, pos, beta_y = beta, se_y = se, pval_y = pval, eaf_y = eaf, n_y = samplesize)],
  by = c("SNP", "chr", "pos"), all = FALSE
)

cat("  Overlapping SNPs:", nrow(merged), "\n")

# 过滤MAF > 0.01
merged <- merged[eaf_x > 0.01 & eaf_y < 0.99 & eaf_y > 0.01 & eaf_y < 0.99]
cat("  After MAF filter:", nrow(merged), "\n")

# ===== Step 3: 准备LCV输入 =====
cat("\nStep 3: Preparing LCV input...\n")

# LCV需要Z-scores
merged[, z_x := beta_x / se_x]
merged[, z_y := beta_y / se_y]

# 按染色体排序
setorder(merged, chr, pos)

# 写入LCV格式
lcv_input <- merged[, .(SNP, chr, pos, z_x, z_y, n_x, n_y)]
fwrite(lcv_input, file.path(OUT, "LCV_input.txt"), sep = "\t", quote = FALSE)

# ===== Step 4: 运行LCV =====
cat("\nStep 4: Running LCV analysis...\n")
cat("⚠️ LCV运行时间较长，请耐心等待...\n")

tryCatch({
  # LCV主函数
  lcv_result <- LCV::lcv(
    z1 = lcv_input$z_x,
    z2 = lcv_input$z_y,
    n1 = max(lcv_input$n_x, na.rm = TRUE),
    n2 = max(lcv_input$n_y, na.rm = TRUE),
    # 使用block jackknife估计标准误
    nblock = 200
  )
  
  cat("\n=== LCV Results: MDD -> Anxiety ===\n")
  cat("Genetic correlation (rg):", lcv_result$rg, "\n")
  cat("Causal proportion (p):", lcv_result$p, "\n")
  cat("P-value for p≠0:", lcv_result$p_pvalue, "\n")
  cat("Direction:", ifelse(lcv_result$p > 0, "MDD -> Anxiety", "Anxiety -> MDD"), "\n")
  
  # 保存结果
  saveRDS(lcv_result, file.path(OUT, "LCV_result_MDD_Anxiety.rds"))
  
  result_df <- data.table(
    Analysis = "MDD_Anxiety",
    Genetic_Correlation_rg = lcv_result$rg,
    Causal_Proportion_p = lcv_result$p,
    p_pvalue = lcv_result$p_pvalue,
    Direction = ifelse(lcv_result$p > 0, "MDD -> Anxiety", "Anxiety -> MDD")
  )
  
  fwrite(result_df, file.path(OUT, "LCV_summary.csv"))
  
}, error = function(e) {
  cat("LCV Error:", conditionMessage(e), "\n")
  cat("⚠️ LCV可能需要较长时间运行，或需要更多SNP\n")
})

# ===== Step 5: 反向分析（可选） =====
cat("\nStep 5: Reverse direction (Anxiety -> MDD)...\n")

tryCatch({
  lcv_result_rev <- LCV::lcv(
    z1 = lcv_input$z_y,  # Anxiety作为暴露
    z2 = lcv_input$z_x,  # MDD作为结局
    n1 = max(lcv_input$n_y, na.rm = TRUE),
    n2 = max(lcv_input$n_x, na.rm = TRUE),
    nblock = 200
  )
  
  cat("\n=== LCV Results: Anxiety -> MDD ===\n")
  cat("Genetic correlation (rg):", lcv_result_rev$rg, "\n")
  cat("Causal proportion (p):", lcv_result_rev$p, "\n")
  cat("P-value for p≠0:", lcv_result_rev$p_pvalue, "\n")
  
  saveRDS(lcv_result_rev, file.path(OUT, "LCV_result_Anxiety_MDD.rds"))
  
}, error = function(e) {
  cat("LCV Reverse Error:", conditionMessage(e), "\n")
})

cat("\n✅ LCV分析完成！\n")
cat("结果保存在:", OUT, "\n")
