"""Stage 6: the matched comparison the reviewer asked for, plus the diagnostics behind it.

For each instrument set that appears in the new Table 3 this computes
  (a) the univariable multiplicative random-effects IVW estimate on exactly those variants, so the direct
      estimate is compared with the right reference rather than with the 198- or 50-variant primary set;
  (b) the primary multivariable estimate as fitted by stage 5 (MVMR 0.4.8, dispersion-scaled WLS);
  (c) the Sanderson-Windmeijer measurement-error-corrected form, with the per-SNP covariance of the exposure
      estimates taken as diagonal, which is the estimator the manuscript cites;
  (d) a sensitivity that also removes every variant reaching |z| > 5.45 on the outcome side, i.e. not just
      the variants that were formally selected for the outcome;
  (e) the number of strand-ambiguous A/T and C/G variants in the set and the estimate with them removed,
      because the covariate matrix was aligned to a single reference allele.

It also records why the forward model has 84 rather than 198 instruments: coverage of each covariate GWAS at
the primary instrument list, which is what the methods section now has to state.

Writes mvmr_independent/mvmr_matched_comparison.csv and prints the coverage table.
"""
import os
import numpy as np
import pandas as pd

DIR = os.path.dirname(os.path.abspath(__file__))
M = pd.read_csv(f'{DIR}/mvmr_matrix_harmonised.csv')
ST5 = pd.read_csv(f'{DIR}/mvmr_exposure_only.csv')
S1 = pd.read_csv('data/S1_MDD_to_Anxiety_harmonised_198_instruments.csv')
S2 = pd.read_csv('data/S2_Anxiety_to_MDD_harmonised_50_instruments.csv')
COV = ['beta_BMI', 'beta_Education', 'beta_SmokingInit']
CSE = ['se_BMI', 'se_Education', 'se_SmokingInit']
PAL = {frozenset(x) for x in [('A', 'T'), ('C', 'G')]}

al = {}
for f in ['BMI_ieu_a_2_on_union_instruments.csv', 'Education_ieu_a_80_on_union_instruments.csv',
          'SmokingInit_ieu_b_4877_on_union_instruments.csv']:
    for _, r in pd.read_csv(f'{DIR}/{f}').iterrows():
        al[r['snp']] = (str(r['ref']).upper(), str(r['alt']).upper())
for s in (S1, S2):
    for _, r in s.iterrows():
        al.setdefault(r['SNP'], (str(r['EA']).upper(), str(r['OA']).upper()))
M['palindromic'] = M['snp'].map(lambda x: frozenset(al[x]) in PAL if x in al else np.nan)


def re_ivw(bx, by, sey):
    w = 1 / sey ** 2
    b = np.sum(w * bx * by) / np.sum(w * bx ** 2)
    q = np.sum(w * (by - b * bx) ** 2)
    k = len(bx)
    se = np.sqrt(1 / np.sum(w * bx ** 2)) * max(1.0, np.sqrt(q / (k - 1)))
    return b, se, q / (k - 1)


def sw_corrected(rows, xb, yb, ys):
    X = rows[[xb] + COV].values
    seX = rows[[xb.replace('beta', 'se')] + CSE].values
    y = rows[yb].values
    sey = rows[ys].values
    w = 1 / sey ** 2
    XtWX = X.T @ (X * w[:, None])
    A = np.zeros((X.shape[1], X.shape[1]))
    for j in range(len(rows)):
        A += w[j] * np.diag(seX[j] ** 2)
    C = np.linalg.inv(XtWX - A)
    b = C @ (X.T @ (y * w))
    se = np.sqrt(np.diag(C))
    r = y - X @ b
    disp = np.sqrt(np.sum(w * r ** 2) / (len(rows) - X.shape[1]))
    return b, se, disp


SETS = {
    'forward | exposure-only union': dict(xb='beta_MDD', yb='beta_Anxiety', ys='se_Anxiety',
                                          mask=lambda d: d.is_Anx_iv == 0),
    'forward | covariates from exposure instruments': dict(xb='beta_MDD', yb='beta_Anxiety', ys='se_Anxiety',
                                                           mask=lambda d: (d.is_MDD_iv == 1) & (d.indep_BMI == 0) &
                                                                          (d.indep_Education == 0) & (d.indep_SmokingInit == 0)),
    'reverse | exposure-only union': dict(xb='beta_Anxiety', yb='beta_MDD', ys='se_MDD',
                                          mask=lambda d: d.is_MDD_iv == 0),
    'reverse | covariates from exposure instruments': dict(xb='beta_Anxiety', yb='beta_MDD', ys='se_MDD',
                                                           mask=lambda d: (d.is_Anx_iv == 1) & (d.indep_BMI == 0) &
                                                                          (d.indep_Education == 0) & (d.indep_SmokingInit == 0)),
}

rows_out = []
for name, spec in SETS.items():
    rows = M[spec['mask'](M)]
    xb, yb, ys = spec['xb'], spec['yb'], spec['ys']
    ub, use, uq = re_ivw(rows[xb].values, rows[yb].values, rows[ys].values)
    p = ST5[ST5['set'] == name].iloc[0]
    cb, cse, cdisp = sw_corrected(rows, xb, yb, ys)
    z_out = np.abs(rows[yb].values / rows[ys].values)
    strict = rows[z_out < 5.45]
    sb, sse, _ = sw_corrected(strict, xb, yb, ys)
    pb, pse, _ = re_ivw(strict[xb].values, strict[yb].values, strict[ys].values)
    npal = int((rows.palindromic == True).sum())
    nop = rows[rows.palindromic != True]
    ob, ose, _ = sw_corrected(nop, xb, yb, ys)
    rows_out.append(dict(
        set=name, n=len(rows), palindromic=npal,
        uni_beta=round(ub, 4), uni_se=round(use, 4), uni_or=round(np.exp(ub), 3),
        uni_lo=round(np.exp(ub - 1.96 * use), 3), uni_hi=round(np.exp(ub + 1.96 * use), 3), uni_q_per_df=round(uq, 2),
        pkg_beta=round(p['beta'], 4), pkg_se=round(p['se'], 4), pkg_or=round(p['or_'], 3),
        pkg_lo=round(p['lo'], 3), pkg_hi=round(p['hi'], 3), pkg_dispersion=round(p['scale'], 3),
        sw_beta=round(cb[0], 4), sw_se=round(cse[0], 4), sw_or=round(np.exp(cb[0]), 3),
        sw_lo=round(np.exp(cb[0] - 1.96 * cse[0]), 3), sw_hi=round(np.exp(cb[0] + 1.96 * cse[0]), 3),
        sw_dispersion=round(cdisp, 3),
        ratio_sw_over_uni=round(cb[0] / ub, 3), ratio_pkg_over_uni=round(p['beta'] / ub, 3),
        strict_n=len(strict), strict_sw_or=round(np.exp(sb[0]), 3), strict_sw_se=round(sse[0], 4),
        strict_uni_or=round(np.exp(pb), 3), strict_uni_se=round(pse, 4),
        no_palindromic_n=len(nop), no_pal_sw_or=round(np.exp(ob[0]), 3),
        condF_primary=p['condF_primary'], condF_BMI=p['condF_BMI'], condF_Education=p['condF_Education'],
        condF_Smoking=p['condF_Smoking'],
        bmi_coeff=round(sw_corrected(rows, xb, yb, ys)[0][1], 4),
        edu_coeff=round(cb[2], 4), smoke_coeff=round(cb[3], 4),
    ))

res = pd.DataFrame(rows_out)
res.to_csv(f'{DIR}/mvmr_matched_comparison.csv', index=False)
pd.set_option('display.width', 250)
print(res[['set', 'n', 'uni_or', 'pkg_or', 'sw_or', 'sw_lo', 'sw_hi', 'ratio_sw_over_uni',
           'uni_q_per_df', 'strict_n', 'strict_sw_or', 'no_pal_sw_or', 'palindromic']].to_string(index=False))
print()
print(res[['set', 'condF_primary', 'condF_BMI', 'condF_Education', 'condF_Smoking',
           'bmi_coeff', 'edu_coeff', 'smoke_coeff', 'sw_dispersion', 'pkg_dispersion']].to_string(index=False))

print('\n' + '=' * 100)
print('why the forward model has 84 rather than 198 instruments: covariate GWAS coverage at the primary lists')
cov = {}
for f, lab in [('BMI_ieu_a_2', 'BMI'), ('Education_ieu_a_80', 'education'), ('SmokingInit_ieu_b_4877', 'smoking')]:
    d = pd.read_csv(f'{DIR}/{f}_on_primary_instruments.csv')
    cov[lab] = d
mdd_iv = set(S1['SNP']); anx_iv = set(S2['SNP'])

print('\ncoverage of each covariate GWAS at the primary instrument lists (the n column of these VCFs holds an')
print('identifier, not a sample size, so no N is quoted here; only beta and se are used by the fits):')
for lab, d in cov.items():
    have = set(d.loc[d['beta'].notna(), 'snp'])
    print(f'  {lab:9s} has data for {len(have & mdd_iv):3d} of 198 MDD instruments and {len(have & anx_iv):2d} of 50 anxiety instruments')
allthree = set.intersection(*[set(d.loc[d['beta'].notna(), 'snp']) for d in cov.values()])
print(f'  all three    {len(allthree & mdd_iv):3d} of 198 MDD and {len(allthree & anx_iv):2d} of 50 anxiety instruments')
missing_by = {lab: len(mdd_iv) - len(set(d.loc[d['beta'].notna(), 'snp']) & mdd_iv) for lab, d in cov.items()}
print('  MDD instruments missing each covariate GWAS:', missing_by)
