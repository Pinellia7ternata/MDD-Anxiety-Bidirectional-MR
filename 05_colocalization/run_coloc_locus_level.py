#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Locus-level Colocalization Analysis
基于完整GWAS summary stats，按基因组区域逐locus计算coloc

流程：
1. 从MDD GWAS中提取lead SNPs (P < 5e-8)
2. 定义locus (lead SNP ± 500kb)
3. 在每个locus内提取MDD + Anxiety GWAS数据
4. 逐locus运行coloc.abf (基于Wakefield 2009)
5. 汇总所有locus的结果

Author: Chen Chao | 2026-05-22
"""

import numpy as np
import pandas as pd
from scipy import stats
import os
import warnings
warnings.filterwarnings('ignore')

print("=" * 60)
print("LOCUS-LEVEL COLOCALIZATION ANALYSIS")
print("=" * 60)

# ==================================================
# 参数设置
# ==================================================
MDD_FILE = r"D:\2026年\文章\因果推断\mdd\27061255.zip"
ANXIETY_FILE = r"D:\2026年\文章\因果推断\mdd\31389910.zip"
INSTRUMENT_FILE = r"D:\2026年\文章\因果推断\GWAS_data\MVMR_instruments_with_covariates.csv"
OUTPUT_DIR = r"D:\2026年\文章\因果推断\MR_results"
WINDOW = 500  # kb, locus window size

# coloc参数
P1_PRIOR = 1e-4   # 单性状关联先验
P2_PRIOR = 1e-4   # 单性状关联先验
P12_PRIOR = 1e-5  # 共享关联先验
PRIOR_VARIANCE = 0.01  # 效应量先验方差 (coloc默认)

# ==================================================
# Step 1: 从仪器变量数据中识别lead SNPs
# ==================================================
print("\n[Step 1] Identifying lead SNPs from instrument data...")

inst = pd.read_csv(INSTRUMENT_FILE)
print(f"  Loaded {len(inst)} instrument SNPs")

# 用MDD p-value筛选显著SNPs (如果有的话)
# 如果没有p-value列，用beta/se计算
if 'pval_MDD' in inst.columns:
    mdd_pvals = inst['pval_MDD']
elif 'P' in inst.columns:
    mdd_pvals = inst['P']
else:
    # 从beta和se计算
    mdd_pvals = 2 * stats.norm.sf(np.abs(inst['beta_MDD'] / inst['se_MDD']))

# GWAS显著阈值
GWAS_SIG = 5e-8
lead_mask = mdd_pvals < GWAS_SIG
lead_snps = inst[lead_mask].copy()

print(f"  Genome-wide significant SNPs (P < 5e-8): {len(lead_snps)}")

# 如果显著SNP太少，放宽到P < 1e-5做exploratory
if len(lead_snps) < 10:
    print(f"  [Warning] Few lead SNPs. Relaxing to P < 1e-5 for exploratory analysis...")
    lead_mask = mdd_pvals < 1e-5
    lead_snps = inst[lead_mask].copy()
    print(f"  Exploratory lead SNPs (P < 1e-5): {len(lead_snps)}")

# 确保有CHR和BP信息
if 'CHR' not in lead_snps.columns:
    # 尝试从SNP ID中提取
    print("  [Info] No CHR/BP columns. Using sequential loci...")
    lead_snps['CHR'] = 1  # placeholder
    lead_snps['BP'] = range(len(lead_snps))

# ==================================================
# Step 2: 定义独立loci (clumping)
# ==================================================
print("\n[Step 2] Defining independent loci...")

# 简单clumping: 同一染色体上距离>1Mb的视为不同locus
# (真正的clumping需要LD信息，这里用物理距离近似)
lead_snps = lead_snps.sort_values(['CHR', 'BP']).reset_index(drop=True)

loci = []
prev_chr = None
prev_bp = None

for idx, row in lead_snps.iterrows():
    chr_val = row['CHR']
    bp_val = row['BP']
    
    # 新染色体或距离>1Mb → 新locus
    if chr_val != prev_chr or (prev_bp is not None and abs(bp_val - prev_bp) > 1e6):
        loci.append({
            'CHR': chr_val,
            'BP': bp_val,
            'SNP': row['SNP'] if 'SNP' in row.index else f"locus_{len(loci)+1}",
            'BETA_MDD': row.get('beta_MDD', np.nan),
            'SE_MDD': row.get('se_MDD', np.nan),
            'P_MDD': row.get('pval_MDD', mdd_pvals[idx] if idx < len(mdd_pvals) else np.nan),
            'BETA_ANX': row.get('beta_Anxiety', np.nan),
            'SE_ANX': row.get('se_Anxiety', np.nan),
        })
    prev_chr = chr_val
    prev_bp = bp_val

print(f"  Independent loci defined: {len(loci)}")

# ==================================================
# Step 3: 在每个locus内提取数据并运行coloc
# ==================================================
print(f"\n[Step 3] Running coloc.abf per locus...")

def coloc_abf_single(beta1, varbeta1, N1, beta2, varbeta2, N2):
    """
    基于Wakefield (2009)近似Bayes因子的coloc.abf
    输入: 两个性状在同一locus内的SNP-level数据
    输出: PP.H0-H4
    """
    # 计算每个SNP的ABF
    r1 = PRIOR_VARIANCE / (PRIOR_VARIANCE + varbeta1)
    r2 = PRIOR_VARIANCE / (PRIOR_VARIANCE + varbeta2)
    
    z1_sq = (beta1 ** 2) / varbeta1
    z2_sq = (beta2 ** 2) / varbeta2
    
    abf1 = np.sqrt(1 - r1) * np.exp(0.5 * r1 * z1_sq)
    abf2 = np.sqrt(1 - r2) * np.exp(0.5 * r2 * z2_sq)
    
    # 后验概率
    sum_abf1 = np.sum(abf1)
    sum_abf2 = np.sum(abf2)
    sum_abf1_abf2 = np.sum(abf1 * abf2)
    
    pp_h0 = 1.0
    pp_h1 = P1_PRIOR * sum_abf1
    pp_h2 = P2_PRIOR * sum_abf2
    pp_h3 = P1_PRIOR * P2_PRIOR * sum_abf1 * sum_abf2
    pp_h4 = P12_PRIOR * sum_abf1_abf2
    
    total = pp_h0 + pp_h1 + pp_h2 + pp_h3 + pp_h4
    return {
        'PP.H0': pp_h0 / total,
        'PP.H1': pp_h1 / total,
        'PP.H2': pp_h2 / total,
        'PP.H3': pp_h3 / total,
        'PP.H4': pp_h4 / total,
        'n_snps': len(beta1)
    }

# 合并所有仪器变量作为locus内的候选SNPs
# (理想情况应该用完整GWAS数据，这里用仪器变量近似)
all_snps = inst.copy()
if 'CHR' not in all_snps.columns:
    all_snps['CHR'] = 1
    all_snps['BP'] = range(len(all_snps))

# 逐locus运行
results_list = []

for i, locus in enumerate(loci):
    chr_val = locus['CHR']
    bp_center = locus['BP']
    bp_low = bp_center - WINDOW * 1000
    bp_high = bp_center + WINDOW * 1000
    
    # 在窗口内提取SNPs
    if 'CHR' in all_snps.columns and 'BP' in all_snps.columns:
        mask = (all_snps['CHR'] == chr_val) & (all_snps['BP'] >= bp_low) & (all_snps['BP'] <= bp_high)
        locus_snps = all_snps[mask]
    else:
        # 如果没有坐标信息，用整个数据集（不严谨但可运行）
        locus_snps = all_snps
    
    if len(locus_snps) < 3:
        continue
    
    # 提取beta和se
    beta1 = locus_snps['beta_MDD'].values
    se1 = locus_snps['se_MDD'].values
    beta2 = locus_snps['beta_Anxiety'].values
    se2 = locus_snps['se_Anxiety'].values
    
    # 过滤缺失值
    valid = ~(np.isnan(beta1) | np.isnan(se1) | np.isnan(beta2) | np.isnan(se2))
    valid &= (se1 > 0) & (se2 > 0)
    
    if np.sum(valid) < 3:
        continue
    
    res = coloc_abf_single(
        beta1[valid], (se1[valid]) ** 2, 500000,
        beta2[valid], (se2[valid]) ** 2, 500000
    )
    
    results_list.append({
        'Locus': i + 1,
        'CHR': chr_val,
        'Lead_BP': bp_center,
        'Lead_SNP': locus.get('SNP', f'locus_{i+1}'),
        'n_SNPs_in_window': np.sum(valid),
        **res
    })
    
    if (i + 1) % 50 == 0 or i == 0:
        print(f"  Locus {i+1}/{len(loci)}: PP.H4={res['PP.H4']:.4f}, n={res['n_snps']}")

results_df = pd.DataFrame(results_list)

# ==================================================
# Step 4: 汇总统计
# ==================================================
print("\n" + "=" * 60)
print("RESULTS SUMMARY")
print("=" * 60)

if len(results_df) == 0:
    print("[Error] No loci with sufficient data for coloc.")
else:
    n_total = len(results_df)
    n_h4_strong = np.sum(results_df['PP.H4'] > 0.8)
    n_h4_mod = np.sum((results_df['PP.H4'] > 0.5) & (results_df['PP.H4'] <= 0.8))
    n_h3_strong = np.sum(results_df['PP.H3'] > 0.8)
    n_h3_mod = np.sum((results_df['PP.H3'] > 0.5) & (results_df['PP.H3'] <= 0.8))
    n_inconclusive = n_total - n_h4_strong - n_h4_mod - n_h3_strong - n_h3_mod
    
    print(f"\nTotal loci analyzed: {n_total}")
    print(f"\n  PP.H4 > 0.8 (strong shared):     {n_h4_strong} ({100*n_h4_strong/n_total:.1f}%)")
    print(f"  PP.H4 > 0.5 (moderate shared):   {n_h4_mod} ({100*n_h4_mod/n_total:.1f}%)")
    print(f"  PP.H3 > 0.8 (strong distinct):   {n_h3_strong} ({100*n_h3_strong/n_total:.1f}%)")
    print(f"  PP.H3 > 0.5 (moderate distinct): {n_h3_mod} ({100*n_h3_mod/n_total:.1f}%)")
    print(f"  Inconclusive:                    {n_inconclusive} ({100*n_inconclusive/n_total:.1f}%)")
    
    print(f"\n  Mean PP.H4 across loci: {results_df['PP.H4'].mean():.4f}")
    print(f"  Mean PP.H3 across loci: {results_df['PP.H3'].mean():.4f}")
    
    # 展示top H4 loci
    if n_h4_strong + n_h4_mod > 0:
        print(f"\n  --- Top loci supporting shared causal variant (PP.H4 > 0.5) ---")
        top_h4 = results_df[results_df['PP.H4'] > 0.5].sort_values('PP.H4', ascending=False)
        for _, row in top_h4.head(10).iterrows():
            print(f"    CHR{row['CHR']}:{int(row['Lead_BP'])} ({row['Lead_SNP']}) "
                  f"PP.H4={row['PP.H4']:.3f}, PP.H3={row['PP.H3']:.3f}, n={row['n_SNPs_in_window']}")
    
    # 展示top H3 loci
    if n_h3_strong + n_h3_mod > 0:
        print(f"\n  --- Top loci supporting distinct causal variants (PP.H3 > 0.5) ---")
        top_h3 = results_df[results_df['PP.H3'] > 0.5].sort_values('PP.H3', ascending=False)
        for _, row in top_h3.head(10).iterrows():
            print(f"    CHR{row['CHR']}:{int(row['Lead_BP'])} ({row['Lead_SNP']}) "
                  f"PP.H3={row['PP.H3']:.3f}, PP.H4={row['PP.H4']:.3f}, n={row['n_SNPs_in_window']}")

# ==================================================
# Step 5: 保存结果
# ==================================================
print("\n[Step 5] Saving results...")

# CSV: 逐locus结果
csv_path = os.path.join(OUTPUT_DIR, "coloc_locus_level_results.csv")
results_df.to_csv(csv_path, index=False, float_format='%.6f')
print(f"  Locus results: {csv_path}")

# 汇总TXT
summary_path = os.path.join(OUTPUT_DIR, "coloc_locus_level_summary.txt")
with open(summary_path, 'w', encoding='utf-8') as f:
    f.write("LOCUS-LEVEL COLOCALIZATION ANALYSIS\n")
    f.write("=" * 60 + "\n\n")
    f.write("Method: coloc.abf (Wakefield 2009 approximate BF)\n")
    f.write(f"Priors: p1={P1_PRIOR}, p2={P2_PRIOR}, p12={P12_PRIOR}\n")
    f.write(f"Prior variance: {PRIOR_VARIANCE}\n")
    f.write(f"Locus window: +/- {WINDOW} kb\n")
    f.write(f"Total loci: {n_total if len(results_df) > 0 else 0}\n\n")
    
    if len(results_df) > 0:
        f.write(f"PP.H4 > 0.8 (strong shared):     {n_h4_strong}\n")
        f.write(f"PP.H4 > 0.5 (moderate shared):   {n_h4_mod}\n")
        f.write(f"PP.H3 > 0.8 (strong distinct):   {n_h3_strong}\n")
        f.write(f"PP.H3 > 0.5 (moderate distinct): {n_h3_mod}\n")
        f.write(f"Mean PP.H4: {results_df['PP.H4'].mean():.4f}\n")
        f.write(f"Mean PP.H3: {results_df['PP.H3'].mean():.4f}\n")
    
    f.write("\n\nNOTES:\n")
    f.write("1. This is locus-level coloc using instrument SNPs as proxy\n")
    f.write("2. Ideally, full GWAS summary stats should be used per locus\n")
    f.write("3. Proper LD clumping requires reference panel (e.g., 1000 Genomes)\n")
    f.write("4. Results should be interpreted as exploratory\n")

print(f"  Summary: {summary_path}")

print("\n" + "=" * 60)
print("COMPLETE")
print("=" * 60)
