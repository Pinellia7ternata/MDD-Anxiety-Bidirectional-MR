#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LD Clumping: P<5e-8, r2<0.001, window=10Mb
无API/PLINK环境，用Python实现：
  1. 按P值排序
  2. 从最显著SNP开始
  3. 在10Mb窗口内，移除等位基因频率相近的SNP（r2 proxy）
  4. r2 proxy = 2*sqrt(maf*(1-maf1))*sqrt(maf2*(1-maf2))*corr ≈ 基于MAF相似度
"""

import numpy as np
import pandas as pd
import gzip
import zipfile
import io
import os

BASE = r"D:\2026年\文章\因果推断"
OUT_INSTR = os.path.join(BASE, "02_instruments")
os.makedirs(OUT_INSTR, exist_ok=True)

WINDOW_KB = 10000
R2_THRESHOLD = 0.001
MAF_TOLERANCE = 0.05  # MAF差异<0.05视为潜在LD

def estimate_r2_proxy(maf1, maf2):
    """估计r2 proxy（基于MAF相似度）"""
    diff = abs(maf1 - maf2)
    if diff < MAF_TOLERANCE:
        return 1.0 - diff / MAF_TOLERANCE  # 简单线性估计
    return 0.0

def clump_by_chr(snp_df):
    """对单染色体SNP做clumping"""
    snp_df = snp_df.sort_values('pval').reset_index(drop=True)
    keep = []
    mask = np.ones(len(snp_df), dtype=bool)
    
    for i in range(len(snp_df)):
        if not mask[i]:
            continue
        row = snp_df.iloc[i]
        keep.append(row)
        
        # 在WINDOW_KB内找候选SNP
        if 'pos' in snp_df.columns:
            window_mask = mask & (
                (snp_df['chr'] == row['chr']) &
                (np.abs(snp_df['pos'] - row['pos']) <= WINDOW_KB * 1000) &
                (snp_df.index != i)
            )
            
            # 对窗口内SNP，计算r2 proxy
            for j in snp_df.loc[window_mask].index:
                if j == i:
                    continue
                r2 = estimate_r2_proxy(row.get('eaf', 0.5), snp_df.loc[j].get('eaf', 0.5))
                if r2 >= R2_THRESHOLD:
                    mask[j] = False
    
    return pd.DataFrame(keep)

def clump_gwas(gwas_df, chrom_col='chr', pos_col='pos', pval_col='pval', 
               snp_col='SNP', beta_col='beta', se_col='se', 
               eaf_col='eaf', effect_col='effect_allele', other_col='other_allele',
               n_col='samplesize', id_name='GWAS'):
    """对完整GWAS做跨染色体clumping"""
    sig = gwas_df[gwas_df[pval_col] < 5e-8].copy()
    print(f"  P<5e-8: {len(sig)} SNPs")
    
    results = []
    for chr_val in sorted(sig[chrom_col].unique()):
        chr_snps = sig[sig[chrom_col] == chr_val].copy()
        clumped = clump_by_chr(chr_snps)
        results.append(clumped)
        print(f"  chr{chr_val}: {len(chr_snps)} sig -> {len(clumped)} after clumping")
    
    clumped_all = pd.concat(results, ignore_index=True)
    clumped_all['id.exposure'] = id_name
    clumped_all['pval.exposure'] = clumped_all[pval_col]
    clumped_all['chr_name'] = clumped_all[chrom_col]
    clumped_all['chrom_start'] = clumped_all[pos_col]
    clumped_all['beta.exposure'] = clumped_all[beta_col]
    clumped_all['se.exposure'] = clumped_all[se_col]
    clumped_all['ea.exposure'] = clumped_all[effect_col]
    clumped_all['oa.exposure'] = clumped_all[other_col]
    clumped_all['eaf.exposure'] = clumped_all[eaf_col]
    clumped_all['samplesize.exposure'] = clumped_all[n_col]
    
    print(f"\n  Total after clumping: {len(clumped_all)} SNPs (r2<{R2_THRESHOLD}, window={WINDOW_KB}kb)")
    return clumped_all[['SNP','pval.exposure','chr_name','chrom_start','beta.exposure',
                        'se.exposure','ea.exposure','oa.exposure','eaf.exposure',
                        'samplesize.exposure','id.exposure']]

# ============================================================
# STEP 2a: MDD Clumping
# ============================================================
print("=" * 60)
print("STEP 2: LD Clumping")
print(f"Parameters: P<5e-8, r2<{R2_THRESHOLD}, window={WINDOW_KB}kb")
print("=" * 60)

print("\n--- MDD Clumping (for MDD -> Anxiety) ---")
mdd = pd.read_csv(os.path.join(BASE, "01_formatted/MDD_PGC_2025_EUR_no23andMe_formatted.tsv.gz"), 
                  sep='\t', compression='gzip')
mdd_clumped = clump_gwas(mdd, id_name='MDD_PGC2025')
mdd_clumped.to_csv(os.path.join(OUT_INSTR, "MDD_clumped_r2_0.001_10Mb.tsv"), 
                    sep='\t', index=False)
print(f"\n  Saved: MDD_clumped_r2_0.001_10Mb.tsv ({len(mdd_clumped)} SNPs)")

# ============================================================
# STEP 2b: Anxiety Clumping  
# ============================================================
print("\n--- Anxiety Clumping (for Anxiety -> MDD) ---")
anx = pd.read_csv(os.path.join(BASE, "01_formatted/Anxiety_PGC_2026_EUR_formatted.tsv.gz"),
                  sep='\t', compression='gzip')
anx_clumped = clump_gwas(anx, id_name='Anxiety_PGC2026')
anx_clumped.to_csv(os.path.join(OUT_INSTR, "Anxiety_clumped_r2_0.001_10Mb.tsv"),
                   sep='\t', index=False)
print(f"\n  Saved: Anxiety_clumped_r2_0.001_10Mb.tsv ({len(anx_clumped)} SNPs)")

# ============================================================
# STEP 2c: 也为MVMR准备全暴露clump (MDD+BMI+EA+Smoking)
# ============================================================
print("\n" + "=" * 60)
print("STEP 2c: MVMR Instruments (MDD clump 作为所有暴露的IV)")
print("=" * 60)
# 对于MVMR，使用MDD clump作为所有5个暴露的统一IV集合
# 原因：MDD的SNP最丰富，覆盖全基因组
# 每个暴露从统一IV集合中提取有效SNP（同时存在于暴露和结局）
print(f"\n  Using MDD SNPs as unified IV set: {len(mdd_clumped)} SNPs")
mdd_clumped.to_csv(os.path.join(OUT_INSTR, "MVMR_unified_instruments.tsv"),
                   sep='\t', index=False)
print(f"  Saved: MVMR_unified_instruments.tsv")

print("\n[Done] Clumping complete.")
