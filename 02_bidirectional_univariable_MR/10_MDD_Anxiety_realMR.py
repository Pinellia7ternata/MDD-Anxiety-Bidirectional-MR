# -*- coding: utf-8 -*-
"""
MR分析 - MDD -> Anxiety (真实PGC GWAS数据)
暴露: MDD (PGC mdd2025, Cell 2025)
结局: Anxiety (PGC anx2026, Nature Genetics 2026)
"""
import zipfile, gzip, io, os, time
import pandas as pd
import numpy as np
from scipy import stats

DATA_DIR = r"D:\2026年\文章\因果推断\GWAS_data"
OUTPUT_DIR = r"D:\2026年\文章\因果推断\MR_results"
os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(OUTPUT_DIR, exist_ok=True)

MDD_ZIP = r"D:\2026年\文章\因果推断\mdd\27061255.zip"
ANX_ZIP = r"D:\2026年\文章\因果推断\mdd\31389910.zip"

# ============================================================
# 1. 读取MDD数据 (mdd2025)
# ============================================================
print("=" * 70)
print("  MR Analysis: MDD (mdd2025) -> Anxiety (anx2026)")
print("=" * 70)

print("\n[1] Loading MDD exposure data...")
t0 = time.time()
zf_mdd = zipfile.ZipFile(MDD_ZIP)
mdd_file = None
for n in zf_mdd.namelist():
    if 'no23andMe_eur' in n and n.endswith('.tsv.gz'):
        mdd_file = n
        break
if not mdd_file:
    # fallback to any eur file
    for n in zf_mdd.namelist():
        if 'eur' in n and n.endswith('.tsv.gz') and 'top10k' not in n:
            mdd_file = n
            break

print(f"  File: {mdd_file}")
data_mdd = zf_mdd.read(mdd_file)
gz_mdd = gzip.GzipFile(fileobj=io.BytesIO(data_mdd))

# Parse line by line (skip ## comments)
rows = []
header = None
for line in gz_mdd:
    raw = line.decode('utf-8', errors='replace').strip()
    if raw.startswith('##') or not raw:
        continue
    if header is None:
        header = raw.split('\t')
        continue
    vals = raw.split('\t')
    if len(vals) >= len(header):
        rows.append(vals)

gz_mdd.close()
zf_mdd.close()

df_exp = pd.DataFrame(rows, columns=header)
print(f"  Loaded {len(df_exp):,} SNPs in {time.time()-t0:.1f}s")
print(f"  Columns: {list(df_exp.columns)}")

# ============================================================
# 2. 读取Anxiety数据 (anx2026)
# ============================================================
print("\n[2] Loading Anxiety outcome data...")
t0 = time.time()
zf_anx = zipfile.ZipFile(ANX_ZIP)
anx_file = None
for n in zf_anx.namelist():
    if n.endswith('.gz') and not n.endswith(('.pdf', '.xls')):
        anx_file = n
        break

print(f"  File: {anx_file}")
data_anx = zf_anx.read(anx_file)
gz_anx = gzip.GzipFile(fileobj=io.BytesIO(data_anx))

rows2 = []
header2 = None
for line in gz_anx:
    raw = line.decode('utf-8', errors='replace').strip()
    if raw.startswith('##') or not raw:
        continue
    if header2 is None:
        header2 = raw.split('\t')
        continue
    vals = raw.split('\t')
    if len(vals) >= len(header2):
        rows2.append(vals)

gz_anx.close()
zf_anx.close()

df_out = pd.DataFrame(rows2, columns=header2)
print(f"  Loaded {len(df_out):,} SNPs in {time.time()-t0:.1f}s")
print(f"  Columns: {list(df_out.columns)}")

# ============================================================
# 3. 数据预处理
# ============================================================
print("\n[3] Data preprocessing...")

# MDD columns: CHROM, POS, ID, EA, NEA, BETA, SE, PVAL, ...
# Strip # from column names
df_exp.columns = [c.lstrip('#') for c in df_exp.columns]
df_exp = df_exp.rename(columns={'ID': 'SNP', 'EA': 'A1', 'NEA': 'A2',
                                 'BETA': 'beta_exp', 'SE': 'se_exp', 'PVAL': 'pval_exp'})
df_exp['SNP'] = df_exp['SNP'].str.strip()
df_exp = df_exp[['SNP', 'CHROM', 'POS', 'A1', 'A2', 'beta_exp', 'se_exp', 'pval_exp']].copy()
df_exp[['beta_exp', 'se_exp', 'pval_exp']] = df_exp[['beta_exp', 'se_exp', 'pval_exp']].astype(float)
df_exp['CHROM'] = df_exp['CHROM'].astype(str).str.replace('chr', '', case=False)

# Anxiety columns: CHR, SNP, BP, A1, A2, ..., OR, SE, P, ...
df_out = df_out.rename(columns={'CHR': 'CHR_out', 'SNP': 'SNP', 'BP': 'BP_out',
                                 'OR': 'or_out', 'SE': 'se_out', 'P': 'pval_out'})
df_out['SNP'] = df_out['SNP'].str.strip()
df_out = df_out[['SNP', 'CHR_out', 'BP_out', 'A1', 'A2', 'or_out', 'se_out', 'pval_out']].copy()
df_out[['or_out', 'se_out', 'pval_out']] = df_out[['or_out', 'se_out', 'pval_out']].astype(float)
df_out['CHR_out'] = df_out['CHR_out'].astype(str).str.replace('chr', '', case=False)

# Merge on SNP
print(f"  MDD SNPs: {len(df_exp):,}")
print(f"  ANX SNPs: {len(df_out):,}")

df_merged = df_exp.merge(df_out, on='SNP', how='inner')
print(f"  Overlapping SNPs: {len(df_merged):,}")

# ============================================================
# 4. Harmonize alleles
# ============================================================
print("\n[4] Harmonizing alleles...")

def harmonize(row):
    """Check allele alignment between exposure and outcome"""
    a1_exp = str(row.get('A1_x', '')).upper()
    a2_exp = str(row.get('A2_x', '')).upper()
    a1_out = str(row.get('A1_y', '')).upper()
    a2_out = str(row.get('A2_y', '')).upper()
    
    # Same orientation
    if a1_exp == a1_out and a2_exp == a2_out:
        return 'same'
    elif a1_exp == a2_out and a2_exp == a1_out:
        return 'flipped'
    else:
        return 'mismatch'

df_merged['orientation'] = df_merged.apply(harmonize, axis=1)
print(f"  Same orientation: {(df_merged.orientation == 'same').sum():,}")
print(f"  Flipped: {(df_merged.orientation == 'flipped').sum():,}")
print(f"  Mismatched: {(df_merged.orientation == 'mismatch').sum():,}")

# Keep only valid pairs
df_harm = df_merged[df_merged.orientation != 'mismatch'].copy()
print(f"  After removing mismatches: {len(df_harm):,}")

# Flip outcome for flipped pairs
flip_mask = df_harm.orientation == 'flipped'
# For outcome with OR: flip means 1/OR, log(OR) becomes -log(OR)
df_harm.loc[flip_mask, 'or_out'] = 1 / df_harm.loc[flip_mask, 'or_out']
# se stays approximately same for log scale

# Convert outcome OR to log(OR) for MR
df_harm['beta_out'] = np.log(df_harm['or_out'])
df_harm['se_out_log'] = df_harm['se_out'] / df_harm['or_out']  # delta method: SE(logOR) ≈ SE(OR)/OR

# ============================================================
# 5. Select instruments (clumping + p-value threshold)
# ============================================================
print("\n[5] Selecting instruments...")

# P-value threshold for significance
P_THRESH = 5e-8
sig = df_harm[df_harm['pval_exp'] < P_THRESH].copy()
print(f"  Genome-wide significant (p<{P_THRESH:.0e}): {len(sig):,}")

if len(sig) < 3:
    print("  WARNING: Too few genome-wide significant SNPs! Relaxing threshold...")
    P_THRESH = 5e-6
    sig = df_harm[df_harm['pval_exp'] < P_THRESH].copy()
    print(f"  Relaxed (p<{P_THRESH:.0e}): {len(sig):,}")

if len(sig) < 3:
    print("  ERROR: Still too few SNPs!")
    exit(1)

# Simple clumping: keep the most significant SNP per LD block (approximate by position)
# Since we don't have LD reference, use proxy clumping by keeping top SNP per 500kb window
sig = sig.sort_values('pval_exp')
sig['abs_pos'] = sig['POS'].astype(int)
kept = []
last_chrom = None
last_pos = -1e9
CLUMP_WINDOW = 500000  # 500kb

for _, row in sig.iterrows():
    chrom = row['CHROM']
    pos = row['abs_pos']
    if chrom != last_chrom or pos - last_pos > CLUMP_WINDOW:
        kept.append(row.name)
        last_chrom = chrom
        last_pos = pos

df_inst = df_harm.loc[kept].copy()
print(f"  After clumping (500kb window): {len(df_inst):,}")

# F-statistic calculation
df_inst['F'] = (df_inst['beta_exp'] / df_inst['se_exp']) ** 2
print(f"  Mean F-statistic: {df_inst['F'].mean():.1f}")
print(f"  Min F-statistic: {df_inst['F'].min():.1f}")
weak_instr = (df_inst['F'] < 10).sum()
print(f"  Weak instruments (F<10): {weak_instr}")

# ============================================================
# 6. MR Analysis
# ============================================================
print("\n" + "=" * 70)
print("  MR RESULTS")
print("=" * 70)

bx = df_inst['beta_exp'].values.astype(float)
by = df_inst['beta_out'].values.astype(float)
se_by = df_inst['se_out_log'].values.astype(float)
rsids = df_inst['SNP'].values
n = len(bx)

def mr_ivw(bx, by, se):
    w = 1 / se**2
    b = np.sum(w * by) / np.sum(w * bx)
    s = np.sqrt(1 / np.sum(w * bx**2))
    z = b / s
    p = 2 * (1 - stats.norm.cdf(abs(z)))
    return b, s, p

def mr_egger(bx, by, se):
    X = np.column_stack([np.ones(n), bx])
    W = np.diag(1 / se**2)
    XtWX = X.T @ W @ X
    XtWy = X.T @ W @ by
    try:
        beta = np.linalg.solve(XtWX, XtWy)
    except:
        beta = np.linalg.lstsq(XtWX, XtWy, rcond=None)[0]
    resid = by - X @ beta
    s2 = np.sum(resid**2) / max(n - 2, 1)
    V = s2 * np.linalg.pinv(XtWX)
    se_s = np.sqrt(max(V[1, 1], 1e-12))
    se_i = np.sqrt(max(V[0, 0], 1e-12))
    p_s = 2 * (1 - stats.norm.cdf(abs(beta[1] / se_s)))
    p_i = 2 * (1 - stats.norm.cdf(abs(beta[0] / se_i)))
    return beta[1], se_s, p_s, beta[0], se_i, p_i

def mr_wmedian(bx, by, se):
    ratios = by / bx
    w = 1 / se**2
    order = np.argsort(ratios)
    ratios, w = ratios[order], w[order]
    cum = np.cumsum(w)
    idx = np.searchsorted(cum, cum[-1] / 2)
    b = ratios[min(idx, len(ratios)-1)]
    s = np.std(ratios) / np.sqrt(len(ratios))
    p = 2 * (1 - stats.norm.cdf(abs(b / max(s, 1e-12))))
    return b, s, p

def cochran_q(bx, by, se):
    w = 1 / se**2
    b_ivw = np.sum(w * by) / np.sum(w * bx)
    Q = np.sum(w * (by - b_ivw * bx)**2)
    df = n - 1
    p_Q = 1 - stats.chi2.cdf(Q, df) if Q > 0 else 1.0
    I2 = max(0, (Q - df) / Q * 100) if Q > 0 else 0
    return Q, df, p_Q, I2

def leave_one_out(bx, by, se, labels):
    results = []
    for i in range(n):
        bx_loo = np.concatenate([bx[:i], bx[i+1:]])
        by_loo = np.concatenate([by[:i], by[i+1:]])
        se_loo = np.concatenate([se[:i], se[i+1:]])
        if len(bx_loo) >= 3:
            b, s, p = mr_ivw(bx_loo, by_loo, se_loo)
            results.append({
                'removed': labels[i],
                'beta': b, 'se': s, 'p': p,
                'or': np.exp(b),
                'or_lo': np.exp(b-1.96*s),
                'or_hi': np.exp(b+1.96*s)
            })
    return pd.DataFrame(results)

# Run analyses
print(f"\n  Exposure: MDD (PGC mdd2025, Cell 2025)")
print(f"  Outcome: Anxiety (PGC anx2026, Nat Genet 2026)")
print(f"  Instruments: {n} SNPs (mean F={df_inst['F'].mean():.1f})")
print()

# IVW
b_ivw, s_ivw, p_ivw = mr_ivw(bx, by, se_by)
or_ivw = np.exp(b_ivw)
ci_lo = np.exp(b_ivw - 1.96*s_ivw)
ci_hi = np.exp(b_ivw + 1.96*s_ivw)
print(f"  {'Method':20s} {'Beta':>10s} {'SE':>10s} {'P-value':>14s} {'OR':>8s} {'95%CI':>22s}")
print(f"  {'-'*74}")
print(f"  {'IVW':20s} {b_ivw:+10.4f} {s_ivw:10.4f} {p_ivw:14.2e} {or_ivw:8.3f} [{ci_lo:.3f}, {ci_hi:.3f}]")

# MR-Egger
b_e, s_e, p_e, int_e, s_i, p_i = mr_egger(bx, by, se_by)
or_e = np.exp(b_e)
print(f"  {'MR-Egger':20s} {b_e:+10.4f} {s_e:10.4f} {p_e:14.2e} {or_e:8.3f}")
print(f"  {'  Intercept':20s} {int_e:+10.4f} {s_i:10.4f} {p_i:14.2e}")

# Weighted Median
b_w, s_w, p_w = mr_wmedian(bx, by, se_by)
or_w = np.exp(b_w)
print(f"  {'Weighted Median':20s} {b_w:+10.4f} {s_w:10.4f} {p_w:14.2e} {or_w:8.3f}")

# Heterogeneity
Q, df_Q, p_Q, I2 = cochran_q(bx, by, se_by)
print(f"\n  Heterogeneity: Q={Q:.2f}, df={df_Q}, P={p_Q:.2e}, I2={I2:.1f}%")

# Directionality test (Egger intercept)
if p_i < 0.05:
    direction = "WARNING: Potential pleiotropy detected (Egger intercept P<0.05)"
else:
    direction = "No evidence of directional pleiotropy"
print(f"  Pleiotropy test: {direction}")

# ============================================================
# 7. Save results
# ============================================================
print("\n[7] Saving results...")

# Main results table
results_df = pd.DataFrame([{
    'method': 'IVW', 'n_snp': n,
    'beta': b_ivw, 'se': s_ivw, 'p': p_ivw,
    'or': or_ivw, 'or_lo': ci_lo, 'or_hi': ci_hi
}, {
    'method': 'MR-Egger', 'n_snp': n,
    'beta': b_e, 'se': s_e, 'p': p_e,
    'or': or_e, 'or_lo': np.nan, 'or_hi': np.nan
}, {
    'method': 'Weighted Median', 'n_snp': n,
    'beta': b_w, 'se': s_w, 'p': p_w,
    'or': or_w, 'or_lo': np.nan, 'or_hi': np.nan
}])
results_file = os.path.join(OUTPUT_DIR, "MDD_Anxiety_MR_results.csv")
results_df.to_csv(results_file, index=False, encoding='utf-8-sig')
print(f"  Results: {results_file}")

# Leave-one-out
loo_df = leave_one_out(bx, by, se_by, rsids)
loo_file = os.path.join(OUTPUT_DIR, "MDD_Anxiety_LOO.csv")
loo_df.to_csv(loo_file, index=False, encoding='utf-8-sig')
print(f"  Leave-one-out: {loo_file}")

# Instrument details
inst_df = df_inst[['SNP', 'CHROM', 'POS', 'A1_x', 'A2_x', 'beta_exp', 'se_exp',
                    'pval_exp', 'or_out', 'se_out', 'pval_out', 'beta_out', 'F']].copy()
inst_df.columns = ['SNP', 'CHR', 'POS', 'EA', 'NEA', 'beta_exp', 'se_exp',
                   'pval_exp', 'or_out', 'se_out', 'pval_out', 'beta_out', 'F_stat']
inst_file = os.path.join(DATA_DIR, "MDD_Anxiety_instruments.csv")
inst_df.to_csv(inst_file, index=False, encoding='utf-8-sig')
print(f"  Instruments: {inst_file}")

# Scatter plot
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    colors = df_inst['CHROM'].astype(int) % 20
    cmap = plt.cm.tab20
    scatter = ax.scatter(bx, by, c=[cmap(c) for c in colors], s=60, alpha=0.7, edgecolors='black', linewidths=0.5)
    
    # IVW regression line
    x_line = np.array([bx.min() - 0.01, bx.max() + 0.01])
    y_line = b_ivw * x_line
    ax.plot(x_line, y_line, 'r-', linewidth=2, label=f'IVW (OR={or_ivw:.2f})')
    
    ax.set_xlabel('MDD Effect (beta)', fontsize=12)
    ax.set_ylabel('Anxiety Effect (log OR)', fontsize=12)
    ax.set_title(f'MDD -> Anxiety Two-Sample MR\n({n} instruments, mean F={df_inst["F"].mean():.1f})', fontsize=14)
    ax.legend(fontsize=11)
    ax.axhline(y=0, color='gray', linestyle='--', linewidth=0.8)
    ax.axvline(x=0, color='gray', linestyle='--', linewidth=0.8)
    ax.grid(True, alpha=0.3)
    
    plot_file = os.path.join(OUTPUT_DIR, "MDD_Anxiety_scatter.png")
    plt.tight_layout()
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Scatter plot: {plot_file}")
except Exception as e:
    print(f"  Plot error: {e}")

# Forest plot (leave-one-out)
try:
    fig2, ax2 = plt.subplots(figsize=(12, max(4, n*0.35)))
    
    loo_sorted = loo_df.sort_values('beta')
    y_pos = range(len(loo_sorted))
    
    # Error bars (use CI width from beta and se)
    loo_err_lo = loo_sorted['beta'] - 1.96 * loo_sorted['se']
    loo_err_hi = 1.96 * loo_sorted['se']
    ax2.errorbar(loo_sorted['beta'].values, y_pos,
                 xerr=[loo_err_lo, loo_err_hi],
                 fmt='o', color='steelblue', ecolor='gray', capsize=3, markersize=6)
    
    # Overall estimate
    ax2.axvline(x=b_ivw, color='red', linestyle='--', linewidth=1.5, label=f'IVW overall (OR={or_ivw:.2f})')
    
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(loo_sorted['removed'].values, fontsize=8)
    ax2.set_xlabel('log(Odds Ratio)', fontsize=11)
    ax2.set_title('Leave-One-Out Sensitivity Analysis', fontsize=13)
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3, axis='x')
    
    forest_file = os.path.join(OUTPUT_DIR, "MDD_Anxiety_forest.png")
    plt.tight_layout()
    plt.savefig(forest_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Forest plot: {forest_file}")
except Exception as e:
    print(f"  Forest plot error: {e}")

print(f"\n{'=' * 70}")
print("  DONE!")
print(f"{'=' * 70}")
