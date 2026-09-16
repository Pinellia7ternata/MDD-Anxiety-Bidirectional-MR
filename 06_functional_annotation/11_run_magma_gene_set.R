#!/usr/bin/env Rscript
# ============================================================
# 11_run_magma_gene_set.R — MAGMA基因集/通路富集分析
# ============================================================
# 方法：de Leeuw et al. 2015 Nat Genet "MAGMA: gene-set analysis 
#       of GWAS data"
# 输入：GWAS summary stats (MDD + Anxiety)
# 输出：基因水平P值、基因集富集、通路分析
# ============================================================
# 前置依赖：
#   1. MAGMA软件（https://ctg.cncr.nl/software/magma）
#   2. 1000G EUR LD参考面板（可从MAGMA官网下载）
#   3. 基因注释文件（NCBI或Ensembl）
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE <- "D:/2026年/文章/因果推断"
OUT  <- file.path(BASE, "11_MAGMA")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ===== Step 0: 检查MAGMA是否已安装 =====
MAGMA_PATH <- "D:/tools/magma/magma.exe"  # 修改为实际路径

if (!file.exists(MAGMA_PATH)) {
  cat("⚠️ MAGMA未安装或路径不正确\n")
  cat("请从以下地址下载并解压：\n")
  cat("  https://ctg.cncr.nl/software/magma\n")
  cat("并修改脚本中的MAGMA_PATH变量\n")
  cat("\n备选方案：使用FUMA在线平台（见12_run_fuma_annotation.R）\n")
  stop("MAGMA not found")
}

# ===== Step 1: 准备MAGMA输入文件 =====
cat("Step 1: Preparing MAGMA input files...\n")

# 读取GWAS数据
mdd <- fread(file.path(BASE, "01_formatted/MDD_PGC_2025_EUR_no23andMe_formatted.tsv.gz"))
anx <- fread(file.path(BASE, "01_formatted/Anxiety_PGC_2026_EUR_formatted.tsv.gz"))

# MAGMA需要的格式：SNP, CHR, POS, P, N, REF, ALT（可选）
magma_input_mdd <- mdd[, .(SNP, CHR = chr, BP = pos, P = pval, N = samplesize)]
magma_input_anx <- anx[, .(SNP, CHR = chr, BP = pos, P = pval, N = samplesize)]

# 写入MAGMA格式
fwrite(magma_input_mdd, file.path(OUT, "MDD_for_magma.txt"), 
       sep = "\t", quote = FALSE)
fwrite(magma_input_anx, file.path(OUT, "Anxiety_for_magma.txt"), 
       sep = "\t", quote = FALSE)

cat("  MDD SNPs:", nrow(magma_input_mdd), "\n")
cat("  Anxiety SNPs:", nrow(magma_input_anx), "\n")

# ===== Step 2: 基因注释（需要基因位置文件） =====
cat("\nStep 2: Gene annotation...\n")

# 基因位置文件路径（需下载）
GENE_LOC <- "D:/tools/magma/NCBI37.3.gene.loc"  # 或Ensembl

if (!file.exists(GENE_LOC)) {
  cat("⚠️ 基因位置文件未找到，请下载：\n")
  cat("  NCBI: https://ctg.cncr.nl/software/magma\n")
  cat("  或从Ensembl生成：biomaRt::getBM()\n")
  stop("Gene location file not found")
}

# 运行MAGMA annotate
cmd_annotate <- sprintf(
  '%s --annotate --snp-loc %s --gene-loc %s --out %s/MDD_annotated',
  MAGMA_PATH,
  file.path(OUT, "MDD_for_magma.txt"),
  GENE_LOC,
  OUT
)

cat("Running:", cmd_annotate, "\n")
system(cmd_annotate)

# ===== Step 3: 基因水平分析（需要LD参考面板） =====
cat("\nStep 3: Gene-level analysis...\n")

LD_REF <- "D:/tools/magma/1000G_EUR.ld_reference"  # 1000G EUR LD面板

if (!file.exists(LD_REF)) {
  cat("⚠️ LD参考面板未找到，请从MAGMA官网下载：\n")
  cat("  g1000_eur.zip (约300MB)\n")
  stop("LD reference panel not found")
}

# 运行MAGMA基因分析
cmd_gene <- sprintf(
  '%s --bfile %s --pval %s N=%d --gene-annot %s --out %s/MDD_gene',
  MAGMA_PATH,
  LD_REF,
  file.path(OUT, "MDD_for_magma.txt"),
  max(mdd$samplesize, na.rm = TRUE),
  file.path(OUT, "MDD_annotated.genes.annot"),
  OUT
)

cat("Running:", cmd_gene, "\n")
system(cmd_gene)

# 对Anxiety重复
cmd_gene_anx <- sprintf(
  '%s --bfile %s --pval %s N=%d --gene-annot %s --out %s/Anxiety_gene',
  MAGMA_PATH,
  LD_REF,
  file.path(OUT, "Anxiety_for_magma.txt"),
  max(anx$samplesize, na.rm = TRUE),
  file.path(OUT, "Anxiety_annotated.genes.annot"),  # 需先生成
  OUT
)

cat("Running:", cmd_gene_anx, "\n")
system(cmd_gene_anx)

# ===== Step 4: 基因集分析 =====
cat("\nStep 4: Gene-set analysis...\n")

# 使用MAGMA内置的GO/KEGG基因集，或自定义
GENE_SETS <- "D:/tools/magma/gene_sets.gmt"  # 基因集文件

if (file.exists(GENE_SETS)) {
  cmd_set <- sprintf(
    '%s --gene-results %s/MDD_gene.genes.raw --set-annot %s --out %s/MDD_geneset',
    MAGMA_PATH,
    OUT,
    GENE_SETS,
    OUT
  )
  
  cat("Running:", cmd_set, "\n")
  system(cmd_set)
}

# ===== Step 5: 解析结果 =====
cat("\nStep 5: Parsing MAGMA results...\n")

# 读取基因水平结果
gene_results_mdd <- fread(file.path(OUT, "MDD_gene.genes.out"))
gene_results_anx <- fread(file.path(OUT, "Anxiety_gene.genes.out"))

# 显著基因（FDR < 0.05）
sig_genes_mdd <- gene_results_mdd[p < 0.05 / nrow(gene_results_mdd)]
sig_genes_anx <- gene_results_anx[p < 0.05 / nrow(gene_results_anx)]

cat("Significant genes (Bonferroni):\n")
cat("  MDD:", nrow(sig_genes_mdd), "\n")
cat("  Anxiety:", nrow(sig_genes_anx), "\n")

# 共享显著基因
shared_sig <- intersect(sig_genes_mdd$GENE, sig_genes_anx$GENE)
cat("  Shared:", length(shared_sig), "\n")

# 保存结果
fwrite(gene_results_mdd, file.path(OUT, "MDD_gene_results.csv"))
fwrite(gene_results_anx, file.path(OUT, "Anxiety_gene_results.csv"))
fwrite(data.table(GENE = shared_sig), file.path(OUT, "shared_significant_genes.csv"))

# ===== Step 6: 可视化 =====
cat("\nStep 6: Visualization...\n")

# Manhattan plot for genes
if (nrow(sig_genes_mdd) > 0) {
  p1 <- ggplot(gene_results_mdd, aes(x = POS, y = -log10(p))) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05 / nrow(gene_results_mdd)), 
               color = "red", linetype = "dashed") +
    theme_minimal() +
    labs(title = "MAGMA Gene-level Association: MDD",
         x = "Genomic Position", y = "-log10(P)")
  
  ggsave(file.path(OUT, "Figure6_MAGMA_MDD_manhattan.png"), 
         p1, width = 12, height = 6, dpi = 150)
}

cat("\n✅ MAGMA分析完成！\n")
cat("结果保存在:", OUT, "\n")
