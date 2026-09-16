#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Strict Locus-level Colocalization Analysis
使用完整GWAS summary stats，按基因组区域逐locus计算coloc

数据源:
  MDD:  pgc-mdd2025 (27061255.zip, EUR, ~928k Neff)
  Anxiety: PGC anx2026 (31389910.zip, ~845k Neff)

流程:
  1. 读取MDD GWAS → 识别lead SNPs (P<5e-8)
  2. Clumping (物理距离 > 250kb = 独立locus)
  3. 定义每个locus: lead SNP +/- 500kb
  4. 在每个locus内提取MDD + Anxiety GWAS数据
  5. 逐locus运行coloc.abf
  6. 汇总

Author: Chen Chao | 2026-05-22
"""

import numpy as np
import pandas as pd
import gzip
import zipfile
import io
import os
import warnings
warnings.filterwarnings('ignore')

print("=" * 60)
print("STRICT LOCUS-LEVEL COLOCALIZATION")
print("=" * 60)

# ==================================================
# 参数
# ==================================================
MDD_ZIP = r"D:\2026年\文章\因果推断\mdd\27061255.zip"
MDD_ENTRY = "daner/daner_pgc_mdd_no23andMe_eur_hg19_v3.49.24.11.neff.gz"
ANX_ZIP = r"D:\2026年\文章\因果推断\mdd\31389910.zip"
ANX_ENTRY = "ANX_2026_daner_fullANX_v12_woUTAH_11022026.gz"
OUTPUT_DIR = r"D:\2026年\文章\因果推断\MR_results"

CLUMP_DIST = 250000  # 250kb clumping distance
WINDOW = 500  # 500kb locus window

P1 = 1e-4; P2 = 1e-4; P12 = 1e-5
SDY = 0.01  # prior variance

# ==================================================
# Step 1: 读取MDD GWAS (仅显著SNPs)
# ==================================================
print("\n[Step 1] Loading MDD GWAS - extracting genome-wide significant SNPs...")

# 读取MDD GWAS，只保留显著SNPs以节省内存
print(f"  Reading: {MDD_ENTRY}")
mdd_sig = []
with zipfile.ZipFile(MDD_ZIP) as zf:
    with gzip.open(io.BytesIO(zf.read(MDD_ENTRY)), 'rt') as f:
        header = f.readline().strip().split('\t')
        for line in f:
            parts = line.strip().split('\t')
            p_val = float(parts[header.index('P')])
            if p_val < 5e-8:
                mdd_sig.append({
                    'CHR': int(parts[header.index('CHR')]),
                    'BP': int(parts[header.index('BP')]),
                    'SNP': parts[header.index('SNP')],
                    'A1': parts[header.index('A1')],
                    'A2': parts[header.index('A2')],
                    'BETA': np.log(float(parts[header.index('OR')])),
                    'SE': float(parts[header.index('SE')]),
                    'P': p_val,
                    'Neff': float(parts[header.index('Neff')])
                })

mdd_sig = pd.DataFrame(mdd_sig)
mdd_sig = mdd_sig.sort_values(['CHR', 'BP']).reset_index(drop=True)
print(f"  MDD genome-wide significant SNPs: {len(mdd_sig)}")

# ==================================================
# Step 2: Clumping (物理距离)
# ==================================================
print("\n[Step 2] Clumping MDD lead SNPs (distance > 250kb)...")

lead_snps = []
prev_chr = None
prev_bp = -1e9

for _, row in mdd_sig.iterrows():
    if row['CHR'] != prev_chr or row['BP'] - prev_bp > CLUMP_DIST:
        lead_snps.append({
            'CHR': row['CHR'],
            'BP': row['BP'],
            'SNP': row['SNP'],
            'A1': row['A1'],
            'A2': row['A2'],
            'P_MDD': row['P']
        })
        prev_chr = row['CHR']
        prev_bp = row['BP']

lead_df = pd.DataFrame(lead_snps)
print(f"  Independent lead SNPs: {len(lead_df)}")

# ==================================================
# Step 3: 读取Anxiety GWAS (逐染色体)
# ==================================================
print("\n[Step 3] Loading Anxiety GWAS...")
print(f"  Reading: {ANX_ENTRY}")

anx_data = []
with zipfile.ZipFile(ANX_ZIP) as zf:
    with gzip.open(io.BytesIO(zf.read(ANX_ENTRY)), 'rt') as f:
        header = f.readline().strip().split('\t')
        col_idx = {h: header.index(h) for h in header}
        for line in f:
            parts = line.strip().split('\t')
            try:
                anx_data.append({
                    'CHR': int(parts[col_idx['CHR']]),
                    'BP': int(parts[col_idx['BP']]),
                    'SNP': parts[col_idx['SNP']],
                    'A1': parts[col_idx['A1']],
                    'A2': parts[col_idx['A2']],
                    'BETA': np.log(float(parts[col_idx['OR']])),
                    'SE': float(parts[col_idx['SE']]),
                    'P': float(parts[col_idx['P']])
                })
            except (ValueError, IndexError):
                continue

anx_df = pd.DataFrame(anx_data)
del anx_data
print(f"  Anxiety SNPs loaded: {len(anx_df)}")

# ==================================================
# Step 4: 对每个locus读取MDD数据 + 运行coloc
# ==================================================
print(f"\n[Step 4] Running locus-level coloc (window: +/- {WINDOW}kb)...")

def coloc_abf(b1, vb1, b2, vb2):
    """coloc.abf using Wakefield (2009) approximate BF"""
    r1 = SDY / (SDY + vb1)
    r2 = SDY / (SDY + vb2)
    z1sq = b1**2 / vb1
    z2sq = b2**2 / vb2
    abf1 = np.sqrt(1 - r1) * np.exp(0.5 * r1 * z1sq)
    abf2 = np.sqrt(1 - r2) * np.exp(0.5 * r2 * z2sq)
    s1 = abf1.sum()
    s2 = abf2.sum()
    s12 = (abf1 * abf2).sum()
    h0 = 1.0
    h1 = P1 * s1
    h2 = P2 * s2
    h3 = P1 * P2 * s1 * s2
    h4 = P12 * s12
    t = h0 + h1 + h2 + h3 + h4
    return h0/t, h1/t, h2/t, h3/t, h4/t

# 读取完整MDD GWAS (逐染色体读取以节省内存)
# 先把MDD按染色体存入dict
print("  Loading MDD GWAS by chromosome...")
mdd_by_chr = {}
with zipfile.ZipFile(MDD_ZIP) as zf:
    with gzip.open(io.BytesIO(zf.read(MDD_ENTRY)), 'rt') as f:
        header = f.readline().strip().split('\t')
        col_idx = {h: header.index(h) for h in header}
        for line in f:
            parts = line.strip().split('\t')
            try:
                chr_val = int(parts[col_idx['CHR']])
                if chr_val not in mdd_by_chr:
                    mdd_by_chr[chr_val] = []
                mdd_by_chr[chr_val].append({
                    'CHR': chr_val,
                    'BP': int(parts[col_idx['BP']]),
                    'SNP': parts[col_idx['SNP']],
                    'A1': parts[col_idx['A1']],
                    'A2': parts[col_idx['A2']],
                    'BETA': np.log(float(parts[col_idx['OR']])),
                    'SE': float(parts[col_idx['SE']]),
                    'P': float(parts[col_idx['P']])
                })
            except (ValueError, IndexError):
                continue

for c in mdd_by_chr:
    mdd_by_chr[c] = pd.DataFrame(mdd_by_chr[c])

print(f"  MDD chromosomes loaded: {len(mdd_by_chr)}")

# 逐locus运行coloc
results = []
n_processed = 0

for _, lead in lead_df.iterrows():
    chr_val = int(lead['CHR'])
    bp_center = int(lead['BP'])
    bp_low = bp_center - WINDOW * 1000
    bp_high = bp_center + WINDOW * 1000
    
    # 提取locus内MDD数据
    if chr_val not in mdd_by_chr:
        continue
    mdd_locus = mdd_by_chr[chr_val]
    mdd_locus = mdd_locus[(mdd_locus['BP'] >= bp_low) & (mdd_locus['BP'] <= bp_high)]
    
    if len(mdd_locus) < 5:
        continue
    
    # 提取locus内Anxiety数据
    anx_locus = anx_df[(anx_df['CHR'] == chr_val) & (anx_df['BP'] >= bp_low) & (anx_df['BP'] <= bp_high)]
    
    if len(anx_locus) < 5:
        continue
    
    # 通过SNP匹配，确保allele方向一致
    merged = mdd_locus.merge(anx_locus, on='SNP', suffixes=('_mdd', '_anx'))
    
    # allele alignment: 要求A1/A2一致 (简化处理，忽略flips)
    same_dir = (merged['A1_mdd'] == merged['A1_anx']) & (merged['A2_mdd'] == merged['A2_anx'])
    flip_dir = (merged['A1_mdd'] == merged['A2_anx']) & (merged['A2_mdd'] == merged['A1_anx'])
    
    # 使用方向一致的SNPs
    valid = same_dir | flip_dir
    merged = merged[valid].copy()
    
    # 翻转allele不一致的SNPs
    flip_mask = merged.index.isin(merged[flip_dir].index)
    merged.loc[flip_mask, 'BETA_anx'] = -merged.loc[flip_mask, 'BETA_anx']
    
    # 过滤
    mask = (merged['SE_mdd'] > 0) & (merged['SE_anx'] > 0) & ~merged['BETA_mdd'].isna() & ~merged['BETA_anx'].isna()
    merged = merged[mask]
    
    if len(merged) < 5:
        continue
    
    # 运行coloc
    h0, h1, h2, h3, h4 = coloc_abf(
        merged['BETA_mdd'].values, merged['SE_mdd'].values**2,
        merged['BETA_anx'].values, merged['SE_anx'].values**2
    )
    
    results.append({
        'Locus': len(results) + 1,
        'CHR': chr_val,
        'Lead_BP': bp_center,
        'Lead_SNP': lead['SNP'],
        'P_MDD': lead['P_MDD'],
        'n_MDD': len(mdd_locus),
        'n_Anxiety': len(anx_locus),
        'n_shared': len(merged),
        'PP.H0': h0, 'PP.H1': h1, 'PP.H2': h2, 'PP.H3': h3, 'PP.H4': h4
    })
    
    n_processed += 1
    if n_processed % 20 == 0:
        print(f"  Processed {n_processed}/{len(lead_df)} loci...")

print(f"  Total loci with coloc results: {len(results)}")

# ==================================================
# Step 5: 汇总
# ==================================================
print("\n" + "=" * 60)
print("RESULTS SUMMARY")
print("=" * 60)

res_df = pd.DataFrame(results)
n = len(res_df)
if n == 0:
    print("[Error] No loci had sufficient data.")
else:
    h4_80 = (res_df['PP.H4'] > 0.8).sum()
    h4_50 = ((res_df['PP.H4'] > 0.5) & (res_df['PP.H4'] <= 0.8)).sum()
    h3_80 = (res_df['PP.H3'] > 0.8).sum()
    h3_50 = ((res_df['PP.H3'] > 0.5) & (res_df['PP.H3'] <= 0.8)).sum()
    incon = n - h4_80 - h4_50 - h3_80 - h3_50
    
    print(f"\nTotal loci: {n}")
    print(f"\n  PP.H4 > 0.8 (strong shared causal):     {h4_80:>4} ({100*h4_80/n:.1f}%)")
    print(f"  0.5 < PP.H4 <= 0.8 (moderate shared):    {h4_50:>4} ({100*h4_50/n:.1f}%)")
    print(f"  PP.H3 > 0.8 (strong distinct):           {h3_80:>4} ({100*h3_80/n:.1f}%)")
    print(f"  0.5 < PP.H3 <= 0.8 (moderate distinct):  {h3_50:>4} ({100*h3_50/n:.1f}%)")
    print(f"  Inconclusive:                             {incon:>4} ({100*incon/n:.1f}%)")
    
    print(f"\n  Median PP.H4: {res_df['PP.H4'].median():.4f}")
    print(f"  Median PP.H3: {res_df['PP.H3'].median():.4f}")
    
    # Top H4 loci
    top_h4 = res_df[res_df['PP.H4'] > 0.5].nlargest(10, 'PP.H4')
    if len(top_h4) > 0:
        print(f"\n  --- Shared causal variants (PP.H4 > 0.5) ---")
        for _, r in top_h4.iterrows():
            print(f"    Locus {r['Locus']:>3}: chr{r['CHR']}:{int(r['Lead_BP']):>12} ({r['Lead_SNP']}) "
                  f"H4={r['PP.H4']:.3f} H3={r['PP.H3']:.3f} n_shared={r['n_shared']}")
    
    # Top H3 loci  
    top_h3 = res_df[res_df['PP.H3'] > 0.5].nlargest(10, 'PP.H3')
    if len(top_h3) > 0:
        print(f"\n  --- Distinct causal variants (PP.H3 > 0.5) ---")
        for _, r in top_h3.iterrows():
            print(f"    Locus {r['Locus']:>3}: chr{r['CHR']}:{int(r['Lead_BP']):>12} ({r['Lead_SNP']}) "
                  f"H3={r['PP.H3']:.3f} H4={r['PP.H4']:.3f} n_shared={r['n_shared']}")

# ==================================================
# Step 6: 保存
# ==================================================
print("\n[Step 6] Saving results...")

csv_path = os.path.join(OUTPUT_DIR, "coloc_locus_level_strict.csv")
res_df.to_csv(csv_path, index=False, float_format='%.6f')
print(f"  Results: {csv_path}")

summary_path = os.path.join(OUTPUT_DIR, "coloc_locus_level_strict_summary.txt")
with open(summary_path, 'w', encoding='utf-8') as f:
    f.write("STRICT LOCUS-LEVEL COLOCALIZATION\n")
    f.write("=" * 50 + "\n\n")
    f.write("Data:\n")
    f.write(f"  MDD: PGC mdd2025 EUR (no23andMe)\n")
    f.write(f"  Anxiety: PGC anx2026 (fullANX woUTAH)\n\n")
    f.write("Method:\n")
    f.write(f"  coloc.abf (Wakefield 2009)\n")
    f.write(f"  Priors: p1={P1}, p2={P2}, p12={P12}, sdy={SDY}\n")
    f.write(f"  Clumping: {CLUMP_DIST/1000:.0f}kb\n")
    f.write(f"  Locus window: +/- {WINDOW}kb\n\n")
    f.write(f"Total loci: {n}\n")
    if n > 0:
        f.write(f"PP.H4>0.8: {h4_80}\n")
        f.write(f"PP.H4>0.5: {h4_80+h4_50}\n")
        f.write(f"PP.H3>0.8: {h3_80}\n")
        f.write(f"PP.H3>0.5: {h3_80+h3_50}\n")
        f.write(f"Median PP.H4: {res_df['PP.H4'].median():.4f}\n")
        f.write(f"Median PP.H3: {res_df['PP.H3'].median():.4f}\n")

print(f"  Summary: {summary_path}")
print("\n" + "=" * 60)
print("COMPLETE")
print("=" * 60)
