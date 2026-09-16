#!/usr/bin/env Rscript
# 使用真实LD clumping后的SNP重新跑MR分析
# 输入：02_instruments_realLD/MDD_realLD_clumped.tsv, Anxiety_realLD_clumped.tsv
# 输出：04_univariable_MR_realLD/

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
})

BASE <- "D:/2026年/文章/因果推断"
OUT  <- file.path(BASE, "04_univariable_MR_realLD")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ===== 1. 读取真实LD clumping后的SNP =====
cat("=== Reading real LD clumped SNPs ===\n")

clumped_mdd   <- fread(file.path(BASE, "02_instruments_realLD/MDD_realLD_clumped.tsv"))
clumped_anx   <- fread(file.path(BASE, "02_instruments_realLD/Anxiety_realLD_clumped.tsv"))

cat("MDD clumped SNPs:", nrow(clumped_mdd), "\n")
cat("Anxiety clumped SNPs:", nrow(clumped_anx), "\n")

# ===== 2. 读取暴露和结局的格式化数据 =====
cat("\n=== Reading formatted GWAS ===\n")

# 注意：列名是 beta, se, pval, eaf, effect_allele, other_allele, chr, pos
mdd_formatted <- fread(file.path(BASE, "01_formatted/MDD_PGC_2025_EUR_no23andMe_formatted.tsv.gz"))
anx_formatted <- fread(file.path(BASE, "01_formatted/Anxiety_PGC_2026_EUR_formatted.tsv.gz"))

# ===== 3. 提取clumped SNPs的数据 =====
cat("\n=== Extracting clumped SNPs ===\n")

# MDD -> Anxiety
mdd_exp <- mdd_formatted[SNP %in% clumped_mdd$SNP]
anx_out  <- anx_formatted[SNP %in% clumped_mdd$SNP]

cat("MDD->Anxiety: ", nrow(mdd_exp), "SNPs extracted\n")

# Anxiety -> MDD
anx_exp <- anx_formatted[SNP %in% clumped_anx$SNP]
mdd_out  <- mdd_formatted[SNP %in% clumped_anx$SNP]

cat("Anxiety->MDD:", nrow(anx_exp), "SNPs extracted\n")

# ===== 4. 格式化为TwoSampleMR对象 =====
cat("\n=== Formatting for TwoSampleMR ===\n")

# 转换为data.frame（format_data()不支持data.table）
mdd_exp_df <- as.data.frame(mdd_exp)
anx_out_df  <- as.data.frame(anx_out)
anx_exp_df <- as.data.frame(anx_exp)
mdd_out_df  <- as.data.frame(mdd_out)

# MDD -> Anxiety
# 注意：列名是 beta, se, pval, eaf, effect_allele, other_allele, chr, pos
exp_dat <- format_data(mdd_exp_df, type="exposure", 
                      snp_col="SNP", beta_col="beta", se_col="se", pval_col="pval",
                      effect_allele_col="effect_allele", other_allele_col="other_allele",
                      eaf_col="eaf", chr_col="chr", pos_col="pos")
out_dat <- format_data(anx_out_df, type="outcome", 
                      snp_col="SNP", beta_col="beta", se_col="se", pval_col="pval",
                      effect_allele_col="effect_allele", other_allele_col="other_allele",
                      eaf_col="eaf", chr_col="chr", pos_col="pos")

# Harmonise
cat("Harmonising...\n")
dat <- harmonise_data(exp_dat, out_dat)

# Run MR (所有12种方法)
cat("Running MR...\n")
res_all <- mr(dat, method_list=c("mr_egger_regression", "mr_two_sample_ml", "mr_weighted_median",
                                 "mr_penalised_weighted_median", "mr_ivw", "mr_simple_mode",
                                 "mr_weighted_mode", "mr_simple_median"))
write.csv(res_all, file.path(OUT, "MDD_to_Anxiety_realLD_MR_results.csv"), row.names=FALSE)
cat("MDD->Anxiety MR results saved.\n")

# Heterogeneity
het <- mr_heterogeneity(dat)
write.csv(het, file.path(OUT, "MDD_to_Anxiety_realLD_heterogeneity.csv"), row.names=FALSE)

# Pleiotropy
pleio <- mr_pleiotropy_test(dat)
write.csv(pleio, file.path(OUT, "MDD_to_Anxiety_realLD_pleiotropy.csv"), row.names=FALSE)

# Leave-one-out
loo <- mr_leaveoneout(dat)
write.csv(loo, file.path(OUT, "MDD_to_Anxiety_realLD_leaveoneout.csv"), row.names=FALSE)

# ===== 5. Anxiety -> MDD =====
cat("\n=== Running Anxiety -> MDD MR ===\n")

exp_dat2 <- format_data(anx_exp_df, type="exposure", 
                       snp_col="SNP", beta_col="beta", se_col="se", pval_col="pval",
                       effect_allele_col="effect_allele", other_allele_col="other_allele",
                       eaf_col="eaf", chr_col="chr", pos_col="pos")
out_dat2 <- format_data(mdd_out_df, type="outcome", 
                       snp_col="SNP", beta_col="beta", se_col="se", pval_col="pval",
                       effect_allele_col="effect_allele", other_allele_col="other_allele",
                       eaf_col="eaf", chr_col="chr", pos_col="pos")

# Harmonise
cat("Harmonising...\n")
dat2 <- harmonise_data(exp_dat2, out_dat2)

# Run MR
cat("Running MR...\n")
res2 <- mr(dat2, method_list=c("mr_egger_regression", "mr_two_sample_ml", "mr_weighted_median",
                                "mr_penalised_weighted_median", "mr_ivw", "mr_simple_mode",
                                "mr_weighted_mode", "mr_simple_median"))
write.csv(res2, file.path(OUT, "Anxiety_to_MDD_realLD_MR_results.csv"), row.names=FALSE)
cat("Anxiety->MDD MR results saved.\n")

# Heterogeneity
het2 <- mr_heterogeneity(dat2)
write.csv(het2, file.path(OUT, "Anxiety_to_MDD_realLD_heterogeneity.csv"), row.names=FALSE)

# Pleiotropy
pleio2 <- mr_pleiotropy_test(dat2)
write.csv(pleio2, file.path(OUT, "Anxiety_to_MDD_realLD_pleiotropy.csv"), row.names=FALSE)

# Leave-one-out
loo2 <- mr_leaveoneout(dat2)
write.csv(loo2, file.path(OUT, "Anxiety_to_MDD_realLD_leaveoneout.csv"), row.names=FALSE)

cat("\n=== Done ===\n")
cat("Results saved to:", OUT, "\n")
