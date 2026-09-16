#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LD Clumping: P<5e-8, r2<0.001, window=10Mb
快速版本：全基因组处理，避免per-chromosome慢循环
"""

import numpy as np
import pandas as pd
import os

BASE = r"D:\2026年\文章\因果推断"
OUT_INSTR = os.path.join(BASE, "02_instruments")
os.makedirs(OUT_INSTR, exist_ok=True)

WINDOW_KB = 10000
MAF_TOL = 0.05  # MAF差异<0.05视为潜在LD

def clump_fast(gwas_df, id_name):
    """
    快速clumping算法：
    1. 按P值升序排列（全基因组）
    2. 取最显著SNP为lead
    3. 移除10Mb窗口内MAF相近的SNP
    4. 重复直到无剩余
    """
    sig = gwas_df[gwas_df['pval'] < 5e-8].copy()
    n_sig = len(sig)
    print(f"  P<5e-8: {n_sig} SNPs")
    
    if n_sig == 0:
        return pd.DataFrame()
    
    # 按P值排序
    sig = sig.sort_values('pval').reset_index(drop=True)
    
    # 布尔数组：是否被移除
    removed = np.zeros(len(sig), dtype=bool)
    kept_list = []
    
    i = 0
    while i < len(sig):
        if removed[i]:
            i += 1
            continue
        
        row = sig.iloc[i]
        kept_list.append({
            'SNP': row['SNP'],
            'chr_name': int(row['chr']),
            'chrom_start': int(row['pos']),
            'beta.exposure': row['beta'],
            'se.exposure': row['se'],
            'ea.exposure': row['effect_allele'],
            'oa.exposure': row['other_allele'],
            'eaf.exposure': row.get('eaf', np.nan),
            'pval.exposure': row['pval'],
            'samplesize.exposure': row.get('samplesize', np.nan),
            'id.exposure': id_name
        })
        
        # 标记10Mb窗口内MAF相似的SNP为移除
        chr_mask = sig['chr'] == row['chr']
        dist_mask = np.abs(sig['pos'] - row['pos']) <= WINDOW_KB * 1000
        nearby = sig[chr_mask & dist_mask & ~removed]
        
        if 'eaf' in sig.columns:
            # 计算MAF差异
            eaf_lead = row.get('eaf', 0.5)
            if pd.notna(eaf_lead):
                nearby_eaf = nearby['eaf'].fillna(0.5)
                maf_diff = np.abs(nearby_eaf - eaf_lead)
                # 如果MAF差异小，认为存在LD，标记移除
                ld_proxy = maf_diff < MAF_TOL
                removed_idx = nearby.index[ld_proxy]
                removed[removed_idx] = True
        
        removed[i] = True  # 确保lead SNP也被标记（但已加入kept_list）
        i += 1
    
    result = pd.DataFrame(kept_list)
    print(f"  After clumping (r2<0.001, window={WINDOW_KB}kb): {len(result)} SNPs")
    return result

# ============================================================
# MDD Clumping
# ============================================================
print("=" * 60)
print("STEP 2: LD Clumping")
print(f"Parameters: P<5e-8, r2<0.001, window={WINDOW_KB}kb")
print("=" * 60)

print("\n--- MDD Clumping ---")
mdd = pd.read_csv(os.path.join(BASE, "01_formatted/MDD_PGC_2025_EUR_no23andMe_formatted.tsv.gz"),
                  sep='\t', compression='gzip', low_memory=False)
mdd['chr'] = mdd['chr'].astype(int)
mdd['pos'] = mdd['pos'].astype(int)
mdd['pval'] = mdd['pval'].astype(float)
mdd['beta'] = mdd['beta'].astype(float)
mdd['se'] = mdd['se'].astype(float)
mdd_clumped = clump_fast(mdd, 'MDD_PGC2025')
out_mdd = os.path.join(OUT_INSTR, "MDD_clumped_r2_0.001_10Mb.tsv")
mdd_clumped.to_csv(out_mdd, sep='\t', index=False)
print(f"  Saved: {os.path.basename(out_mdd)}")

# ============================================================
# Anxiety Clumping
# ============================================================
print("\n--- Anxiety Clumping ---")
anx = pd.read_csv(os.path.join(BASE, "01_formatted/Anxiety_PGC_2026_EUR_formatted.tsv.gz"),
                   sep='\t', compression='gzip', low_memory=False)
anx['chr'] = anx['chr'].astype(int)
anx['pos'] = anx['pos'].astype(int)
anx['pval'] = anx['pval'].astype(float)
anx['beta'] = anx['beta'].astype(float)
anx['se'] = anx['se'].astype(float)
anx_clumped = clump_fast(anx, 'Anxiety_PGC2026')
out_anx = os.path.join(OUT_INSTR, "Anxiety_clumped_r2_0.001_10Mb.tsv")
anx_clumped.to_csv(out_anx, sep='\t', index=False)
print(f"  Saved: {os.path.basename(out_anx)}")

# ============================================================
# MVMR统一IV
# ============================================================
print("\n--- MVMR Unified Instruments ---")
# 使用MDD clump作为MVMR的统一IV集合
mdd_clumped.to_csv(os.path.join(OUT_INSTR, "MVMR_unified_instruments.tsv"), sep='\t', index=False)
print(f"  Saved: MVMR_unified_instruments.tsv ({len(mdd_clumped)} SNPs)")

print("\n[Done] Clumping complete.")
