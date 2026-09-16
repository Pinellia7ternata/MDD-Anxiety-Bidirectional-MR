#!/usr/bin/env Rscript
# ============================================================
# 13_prepare_fuma_input.R — FUMA在线平台输入文件准备
# ============================================================
# FUMA (Functional Mapping and Annotation) 是在线GWAS注释平台
# 网址：https://fuma.ctglab.nl/
# 优势：无需本地安装，自动完成gene mapping、chromatin state、
#        eQTL、MAGMA、通路富集等
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
})

BASE <- "D:/2026年/文章/因果推断"
OUT  <- file.path(BASE, "13_FUMA_input")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ===== Step 1: 准备FUMA输入文件 =====
cat("Step 1: Preparing FUMA input files...\n\n")

# FUMA接受的输入格式：
# 1. GWAS summary statistics（必须包含：SNP, CHR, POS, P）
# 2. 可选：beta, SE, N, EAF, Effect allele, Other allele

mdd <- fread(file.path(BASE, "01_formatted/MDD_PGC_2025_EUR_no23andMe_formatted.tsv.gz"))
anx <- fread(file.path(BASE, "01_formatted/Anxiety_PGC_2026_EUR_formatted.tsv.gz"))

# FUMA标准格式
fuma_cols <- c("SNP", "chr", "pos", "pval", "beta", "se", "eaf", "samplesize")

mdd_fuma <- mdd[, .SD, .SDcols = intersect(names(mdd), fuma_cols)]
anx_fuma <- anx[, .SD, .SDcols = intersect(names(anx), fuma_cols)]

# 重命名列以匹配FUMA要求
setnames(mdd_fuma, 
         c("chr", "pos", "pval", "samplesize"),
         c("CHR", "BP", "P", "N"),
         skip_absent = TRUE)
setnames(anx_fuma,
         c("chr", "pos", "pval", "samplesize"),
         c("CHR", "BP", "P", "N"),
         skip_absent = TRUE)

# 写入FUMA格式文件
fwrite(mdd_fuma, file.path(OUT, "FUMA_input_MDD.txt"), 
       sep = "\t", quote = FALSE)
fwrite(anx_fuma, file.path(OUT, "FUMA_input_Anxiety.txt"),
       sep = "\t", quote = FALSE)

cat("✅ FUMA输入文件已生成：\n")
cat("  MDD: ", file.path(OUT, "FUMA_input_MDD.txt"), " (", nrow(mdd_fuma), " SNPs)\n")
cat("  Anxiety: ", file.path(OUT, "FUMA_input_Anxiety.txt"), " (", nrow(anx_fuma), " SNPs)\n\n")

# ===== Step 2: FUMA上传说明 =====
cat("\n========================================\n")
cat("FUMA在线分析步骤：\n")
cat("========================================\n\n")

cat("1. 访问 FUMA 平台：\n")
cat("   https://fuma.ctglab.nl/\n\n")

cat("2. 创建新Job：\n")
cat("   - Job name: MDD_Anxiety_MR_study\n")
cat("   - Population: European\n")
cat("   - Genome build: hg38（或hg19，根据你的数据）\n\n")

cat("3. 上传GWAS数据：\n")
cat("   - 上传 FUMA_input_MDD.txt\n")
cat("   - 上传 FUMA_input_Anxiety.txt\n")
cat("   - 选择列映射（FUMA会自动识别）\n\n")

cat("4. 配置参数：\n")
cat("   - Lead SNP P-value: 5e-8\n")
cat("   - LD window: 250kb\n")
cat("   - LD reference: 1000G Phase 3 EUR\n")
cat("   - MAF threshold: 0.01\n\n")

cat("5. 选择分析模块：\n")
cat("   ✅ Gene mapping (positional, eQTL, chromatin interaction)\n")
cat("   ✅ MAGMA gene analysis\n")
cat("   ✅ GSEA gene-set analysis\n")
cat("   ✅ Tissue expression (GTEx v8)\n")
cat("   ✅ Pathway enrichment (Reactome, KEGG)\n\n")

cat("6. 运行并等待结果（通常需要2-24小时）\n\n")

cat("7. 下载结果文件：\n")
cat("   - 基因注释表\n")
cat("   - MAGMA结果\n")
cat("   - 通路富集结果\n")
cat("   - eQTL和染色质交互图\n\n")

# ===== Step 3: 生成位置映射文件（可选） =====
cat("\nStep 3: Creating SNP position file for FUMA...\n")

# 如果FUMA需要独立的SNP位置文件
snp_pos <- rbind(
  mdd[, .(SNP, CHR = chr, BP = pos)],
  anx[, .(SNP, CHR = chr, BP = pos)]
)
snp_pos <- unique(snp_pos)

fwrite(snp_pos, file.path(OUT, "SNP_positions.txt"),
       sep = "\t", quote = FALSE)

cat("  Unique SNPs:", nrow(snp_pos), "\n")

# ===== Step 4: 生成lead SNPs列表 =====
cat("\nStep 4: Creating lead SNPs list...\n")

lead_snps <- fread(file.path(BASE, "02_instruments_realLD/MDD_realLD_clumped.tsv"))
lead_anx <- fread(file.path(BASE, "02_instruments_realLD/Anxiety_realLD_clumped.tsv"))

lead_all <- rbind(
  lead_snps[, .(SNP, Trait = "MDD")],
  lead_anx[, .(SNP, Trait = "Anxiety")]
)

fwrite(lead_all, file.path(OUT, "lead_SNPs_for_FUMA.txt"),
       sep = "\t", quote = FALSE)

cat("  MDD lead SNPs:", nrow(lead_snps), "\n")
cat("  Anxiety lead SNPs:", nrow(lead_anx), "\n")

cat("\n✅ FUMA输入文件准备完成！\n")
cat("文件保存在:", OUT, "\n")
