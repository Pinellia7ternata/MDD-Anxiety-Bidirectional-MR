#!/usr/bin/env Rscript
# 真实LD clumping - 使用1000G v3 EUR参考面板
# 依赖：PLINK + 1000G EUR .bed/.bim/.fam

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
})

# ===== 路径设置 =====
BASE     <- "D:/2026年/文章/因果推断"
PLINK    <- "D:/2026年/文章/因果推断/00_reference/plink/plink.exe"
BFILE    <- "D:/2026年/文章/因果推断/00_reference/1000G_EUR/EUR"
OUT      <- file.path(BASE, "02_instruments_realLD")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ===== 读取并格式化暴露数据 =====
# data.table::fread() 可直接读取 .gz 文件，无需外部 gzip 命令
read_formatted <- function(filepath) {
  fread(filepath, sep="\t")
}

# ===== 真实LD clumping函数 =====
run_real_clump <- function(infile, outfile, trait_id) {
  cat("\n=== Clumping:", trait_id, "===\n")
  d <- read_formatted(infile)
  
  # 确保列名符合TwoSampleMR要求
  if (!"chr_name" %in% names(d) && "chr" %in% names(d)) setnames(d, "chr", "chr_name")
  if (!"chrom_start" %in% names(d) && "pos" %in% names(d)) setnames(d, "pos", "chrom_start")
  if (!"pval.exposure" %in% names(d) && "pval" %in% names(d)) setnames(d, "pval", "pval.exposure")
  if (!"beta.exposure" %in% names(d) && "beta" %in% names(d)) setnames(d, "beta", "beta.exposure")
  if (!"se.exposure" %in% names(d) && "se" %in% names(d)) setnames(d, "se", "se.exposure")
  if (!"eaf.exposure" %in% names(d) && "eaf" %in% names(d)) setnames(d, "eaf", "eaf.exposure")
  
  d$id.exposure <- trait_id
  
  # 筛选基因组显著SNP (p < 5e-8)
  sig <- d[pval.exposure < 5e-8]
  cat("Significant SNPs (p<5e-8):", nrow(sig), "\n")
  
  if (nrow(sig) == 0) {
    cat("WARNING: No significant SNPs found!\n")
    return(NULL)
  }
  
  # 转换为data.frame（TwoSampleMR需要）
  sig_df <- as.data.frame(sig)
  
  # 真实LD clumping（使用1000G EUR参考面板）
  # 注意：clump_data()不支持clump_p2参数（这是PLINK命令行参数）
  cat("Running PLINK LD clumping with 1000G EUR reference panel...\n")
  cl <- clump_data(
    sig_df,
    clump_kb    = 10000,
    clump_r2    = 0.001,
    clump_p1    = 5e-8,
    pop         = "EUR",
    bfile       = BFILE,
    plink_bin   = PLINK
  )
  
  if (is.null(cl) || nrow(cl) == 0) {
    cat("WARNING: Clumping returned 0 SNPs. Check PLINK path and BFILE path.\n")
    return(NULL)
  }
  
  cat("Clumped SNPs:", nrow(cl), "\n")
  fwrite(cl, outfile, sep='\t')
  return(cl)
}

# ===== 执行clumping =====
cat("=== Real LD Clumping with 1000G v3 EUR Reference Panel ===\n")
cat("PLINK:", PLINK, "\n")
cat("BFILE:", BFILE, "\n")

# MDD
cl_mdd <- run_real_clump(
  file.path(BASE, "01_formatted/MDD_PGC_2025_EUR_no23andMe_formatted.tsv.gz"),
  file.path(OUT, "MDD_realLD_clumped.tsv"),
  "MDD"
)

# Anxiety
cl_anx <- run_real_clump(
  file.path(BASE, "01_formatted/Anxiety_PGC_2026_EUR_formatted.tsv.gz"),
  file.path(OUT, "Anxiety_realLD_clumped.tsv"),
  "Anxiety"
)

cat("\n=== Done ===\n")
if (!is.null(cl_mdd)) cat("MDD clumped SNPs:", nrow(cl_mdd), "\n")
if (!is.null(cl_anx)) cat("Anxiety clumped SNPs:", nrow(cl_anx), "\n")
