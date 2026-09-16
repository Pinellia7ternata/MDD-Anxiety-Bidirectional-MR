# -*- coding: utf-8 -*-
"""
双向MR分析
1. 正向MR: MDD → Anxiety
2. 反向MR: Anxiety → MDD
3. Steiger方向性检验
"""
import zipfile, gzip, io, os, time
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUTPUT_DIR = r"D:\2026年\文章\因果推断\MR_results"
MDD_ZIP = r"D:\2026年\文章\因果推断\mdd\27061255.zip"
ANX_ZIP = r"D:\2026年\文章\因果推断\mdd\31389910.zip"

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================
# 工具函数
# ============================================================
def load_gwas(zip_path, file_pattern, rename_map, numeric_cols):
    """从zip加载GWAS"""
    print(f"  Loading {os.path.basename(zip_path)}...")
    t0 = time.time()
    zf = zipfile.ZipFile(zip_path)
    target = None
    for n in zf.namelist():
        if file_pattern in n and n.endswith('.gz'):
            target = n
            break
    if not target:
        raise FileNotFoundError(f"No {file_pattern} in zip")
    
    data = zf.read(target)
    gz = gzip.GzipFile(fileobj=io.BytesIO(data))
    
    rows, header = [], None
    for line in gz:
        raw = line.decode('utf-8', errors='replace').strip()
        if raw.startswith('##') or not raw:
            continue
        if header is None:
            header = [c.lstrip('#') for c in raw.split('\t')]
            continue
        vals = raw.split('\t')
        if len(vals) >= len(header):
            rows.append(vals)
    gz.close()
    zf.close()
    
    df = pd.DataFrame(rows, columns=header)
    df = df.rename(columns=rename_map)
    for c in numeric_cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors='coerce')
    print(f"    {len(df):,} SNPs in {time.time()-t0:.1f}s")
    return df

def harmonize(df_exp, df_out):
    """Harmonize alleles"""
    df = df_exp.merge(df_out, on='SNP', how='inner', suffixes=('_exp', '_out'))
    
    def orient(row):
        e1, e2 = str(row.get('A1_exp','')).upper(), str(row.get('A2_exp','')).upper()
        o1, o2 = str(row.get('A1_out','')).upper(), str(row.get('A2_out','')).upper()
        if e1 == o1 and e2 == o2: return 'same'
        if e1 == o2 and e2 == o1: return 'flipped'
        return 'mismatch'
    
    df['orient'] = df.apply(orient, axis=1)
    df = df[df.orient != 'mismatch'].copy()
    
    flip = df.orient == 'flipped'
    if 'or_out' in df.columns:
        df.loc[flip, 'or_out'] = 1 / df.loc[flip, 'or_out']
        df['beta_out'] = np.log(df['or_out'])
        df['se_out_log'] = df['se_out'] / df['or_out']
    elif 'beta_out' in df.columns:
        df.loc[flip, 'beta_out'] = -df.loc[flip, 'beta_out']
        df['se_out_log'] = df['se_out']
    
    return df

def select_iv(df, p_thresh=5e-8, window=500000):
    """Select instruments with clumping"""
    sig = df[df['pval_exp'] < p_thresh].sort_values('pval_exp').copy()
    if len(sig) < 3:
        return None
    
    sig['pos'] = sig['BP_exp'].astype(int) if 'BP_exp' in sig.columns else sig['POS_exp'].astype(int)
    kept, last_chr, last_pos = [], None, -1e9
    
    for _, r in sig.iterrows():
        chr = str(r.get('CHR_exp', r.get('CHROM_exp', '')))
        pos = r['pos']
        if chr != last_chr or pos - last_pos > window:
            kept.append(r.name)
            last_chr, last_pos = chr, pos
    
    df_iv = df.loc[kept].copy()
    df_iv['F'] = (df_iv['beta_exp'] / df_iv['se_exp']) ** 2
    return df_iv

def mr_ivw(bx, by, se):
    w = 1 / se**2
    b = np.sum(w * by) / np.sum(w * bx)
    s = np.sqrt(1 / np.sum(w * bx**2))
    p = 2 * (1 - stats.norm.cdf(abs(b/s)))
    return b, s, p

def mr_egger(bx, by, se):
    n = len(bx)
    X = np.column_stack([np.ones(n), bx])
    W = np.diag(1 / se**2)
    b = np.linalg.lstsq(X.T @ W @ X, X.T @ W @ by, rcond=None)[0]
    resid = by - X @ b
    s2 = np.sum(resid**2) / max(n-2, 1)
    V = s2 * np.linalg.pinv(X.T @ W @ X)
    se_s, se_i = np.sqrt(max(V[1,1], 1e-12)), np.sqrt(max(V[0,0], 1e-12))
    p_s = 2 * (1 - stats.norm.cdf(abs(b[1]/se_s)))
    p_i = 2 * (1 - stats.norm.cdf(abs(b[0]/se_i)))
    return b[1], se_s, p_s, b[0], se_i, p_i

def mr_wm(bx, by, se):
    r = by / bx
    w = 1 / se**2
    order = np.argsort(r)
    r, w = r[order], w[order]
    idx = np.searchsorted(np.cumsum(w), np.cumsum(w)[-1]/2)
    b = r[min(idx, len(r)-1)]
    s = np.std(r) / np.sqrt(len(r))
    p = 2 * (1 - stats.norm.cdf(abs(b/max(s,1e-12))))
    return b, s, p

def run_mr(df, exp, out):
    bx = df['beta_exp'].values.astype(float)
    by = df['beta_out'].values.astype(float)
    se = df['se_out_log'].values if 'se_out_log' in df.columns else df['se_out'].values.astype(float)
    
    b1, s1, p1 = mr_ivw(bx, by, se)
    b2, s2, p2, i2, si2, pi2 = mr_egger(bx, by, se)
    b3, s3, p3 = mr_wm(bx, by, se)
    
    w = 1 / se**2
    Q = np.sum(w * (by - b1 * bx)**2)
    I2 = max(0, (Q - len(bx) + 1) / Q * 100) if Q > 0 else 0
    
    return {
        'exposure': exp, 'outcome': out, 'n': len(bx), 'F': df['F'].mean(),
        'IVW_b': b1, 'IVW_se': s1, 'IVW_p': p1, 'IVW_OR': np.exp(b1),
        'Egger_b': b2, 'Egger_se': s2, 'Egger_p': p2, 'Egger_OR': np.exp(b2),
        'Egger_int': i2, 'Egger_int_p': pi2,
        'WM_b': b3, 'WM_se': s3, 'WM_p': p3, 'WM_OR': np.exp(b3),
        'I2': I2
    }

# ============================================================
# 主程序
# ============================================================
print("=" * 70)
print("  Bidirectional MR: MDD <-> Anxiety")
print("=" * 70)

# Load data
print("\n[1] Loading GWAS...")
df_mdd = load_gwas(MDD_ZIP, 'no23andMe_eur',
                   {'ID':'SNP','EA':'A1','NEA':'A2','BETA':'beta','SE':'se','PVAL':'pval'},
                   ['beta','se','pval'])
df_mdd = df_mdd[['SNP','CHROM','POS','A1','A2','beta','se','pval']].copy()
df_mdd.columns = ['SNP','CHR','BP','A1','A2','beta_mdd','se_mdd','pval_mdd']

df_anx = load_gwas(ANX_ZIP, 'daner_fullANX',
                   {'SNP':'SNP','A1':'A1','A2':'A2','OR':'or','SE':'se','P':'pval'},
                   ['or','se','pval'])
df_anx = df_anx[['SNP','CHR','BP','A1','A2','or','se','pval']].copy()
df_anx.columns = ['SNP','CHR','BP','A1','A2','or_anx','se_anx','pval_anx']
df_anx['beta_anx'] = np.log(df_anx['or_anx'])
df_anx['se_anx_log'] = df_anx['se_anx'] / df_anx['or_anx']

# Forward: MDD -> Anxiety
print("\n[2] Forward MR: MDD -> Anxiety...")
df_exp = df_mdd.rename(columns={'beta_mdd':'beta_exp','se_mdd':'se_exp','pval_mdd':'pval_exp',
                                 'BP':'BP_exp','A1':'A1_exp','A2':'A2_exp'})
df_out = df_anx.rename(columns={'or_anx':'or_out','se_anx_log':'se_out_log','beta_anx':'beta_out',
                                 'A1':'A1_out','A2':'A2_out','se_anx':'se_out'})
df_harm_fwd = harmonize(df_exp, df_out)
print(f"  Harmonized: {len(df_harm_fwd):,}")

df_iv_fwd = select_iv(df_harm_fwd)
if df_iv_fwd is not None:
    res_fwd = run_mr(df_iv_fwd, 'MDD', 'Anxiety')
    print(f"  IVW: OR={res_fwd['IVW_OR']:.3f}, P={res_fwd['IVW_p']:.2e}")
    print(f"  Egger: OR={res_fwd['Egger_OR']:.3f}, P={res_fwd['Egger_p']:.2e}")
    print(f"  WM: OR={res_fwd['WM_OR']:.3f}, P={res_fwd['WM_p']:.2e}")
else:
    res_fwd = None
    print("  ERROR: No instruments!")

# Reverse: Anxiety -> MDD
print("\n[3] Reverse MR: Anxiety -> MDD...")
df_exp2 = df_anx.rename(columns={'beta_anx':'beta_exp','se_anx_log':'se_exp','pval_anx':'pval_exp',
                                  'BP':'BP_exp','A1':'A1_exp','A2':'A2_exp'})
df_out2 = df_mdd.rename(columns={'beta_mdd':'beta_out','se_mdd':'se_out',
                                  'A1':'A1_out','A2':'A2_out'})
df_out2['se_out_log'] = df_out2['se_out']
df_harm_rev = harmonize(df_exp2, df_out2)
print(f"  Harmonized: {len(df_harm_rev):,}")

df_iv_rev = select_iv(df_harm_rev)
if df_iv_rev is not None:
    res_rev = run_mr(df_iv_rev, 'Anxiety', 'MDD')
    print(f"  IVW: OR={res_rev['IVW_OR']:.3f}, P={res_rev['IVW_p']:.2e}")
    print(f"  Egger: OR={res_rev['Egger_OR']:.3f}, P={res_rev['Egger_p']:.2e}")
    print(f"  WM: OR={res_rev['WM_OR']:.3f}, P={res_rev['WM_p']:.2e}")
else:
    res_rev = None
    print("  ERROR: No instruments!")

# Summary
print("\n" + "=" * 70)
print("  RESULTS SUMMARY")
print("=" * 70)

if res_fwd and res_rev:
    print(f"\n  Forward (MDD→Anxiety):")
    print(f"    n={res_fwd['n']}, F={res_fwd['F']:.1f}, I2={res_fwd['I2']:.1f}%")
    print(f"    IVW OR={res_fwd['IVW_OR']:.3f} [{np.exp(res_fwd['IVW_b']-1.96*res_fwd['IVW_se']):.3f}, {np.exp(res_fwd['IVW_b']+1.96*res_fwd['IVW_se']):.3f}]")
    
    print(f"\n  Reverse (Anxiety→MDD):")
    print(f"    n={res_rev['n']}, F={res_rev['F']:.1f}, I2={res_rev['I2']:.1f}%")
    print(f"    IVW OR={res_rev['IVW_OR']:.3f} [{np.exp(res_rev['IVW_b']-1.96*res_rev['IVW_se']):.3f}, {np.exp(res_rev['IVW_b']+1.96*res_rev['IVW_se']):.3f}]")
    
    # Save
    df_res = pd.DataFrame([res_fwd, res_rev])
    df_res.to_csv(os.path.join(OUTPUT_DIR, "Bidirectional_MR_results.csv"), index=False, encoding='utf-8-sig')
    
    # Plot
    fig, ax = plt.subplots(figsize=(12, 6))
    methods = ['IVW', 'Egger', 'WM']
    
    for i, m in enumerate(methods):
        # Forward
        b, s = res_fwd[f'{m}_b'], res_fwd[f'{m}_se']
        ax.errorbar(np.exp(b), 3-i-0.15, xerr=[[np.exp(b)-np.exp(b-1.96*s)], [np.exp(b+1.96*s)-np.exp(b)]],
                    fmt='o', color='steelblue', markersize=10, capsize=4, label='MDD→Anxiety' if i==0 else '')
        # Reverse
        b2, s2 = res_rev[f'{m}_b'], res_rev[f'{m}_se']
        ax.errorbar(np.exp(b2), 3-i+0.15, xerr=[[np.exp(b2)-np.exp(b2-1.96*s2)], [np.exp(b2+1.96*s2)-np.exp(b2)]],
                    fmt='s', color='coral', markersize=10, capsize=4, label='Anxiety→MDD' if i==0 else '')
    
    ax.axvline(1, color='gray', ls='--', lw=1.5)
    ax.set_yticks([0.85, 1.85, 2.85])
    ax.set_yticklabels(['Weighted Median', 'MR-Egger', 'IVW'], fontsize=11)
    ax.set_xlabel('Odds Ratio', fontsize=12)
    ax.set_title('Bidirectional MR: MDD ↔ Anxiety', fontsize=14, fontweight='bold')
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3, axis='x')
    
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "Bidirectional_MR_forest.png"), dpi=300, bbox_inches='tight')
    plt.close()
    print(f"\n  Plot saved: {OUTPUT_DIR}/Bidirectional_MR_forest.png")

print("\n" + "=" * 70)
print("  DONE")
print("=" * 70)
