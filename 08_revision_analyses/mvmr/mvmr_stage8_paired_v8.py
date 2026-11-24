"""Stage 8: the fifth review round's numbers.

Five things are computed here, all of them from the shipped instrument tables, all of them feeding Table 3,
Supplementary Table S7 and Figure 3 of v8.

1. PAIRED COMPARISON UNDER ONE ESTIMATOR. In v7 the direct column was the measurement-error-corrected fit and
   the univariable column was plain multiplicative random-effects IVW, so the ratio mixed two specifications;
   the correction alone inflates a beta by 2.4-13.8%. Here both sides are computed under the corrected
   estimator and, as a check, both under the plain one.
2. JOINT TEST OF THE COVARIATES in each model, from the corrected covariance (chi-square on 3 df, and the
   dispersion-scaled F). MVMR 0.4.8 exports no difference test between the univariable and multivariable
   estimate of the same exposure, so that part is reported as not done rather than approximated.
3. LD OVERLEFT BY POSITION CLUMPING: covariate instruments lying within 1 Mb, and within 100 kb, of a primary
   instrument, and the refit with them removed.
4. THE 114 MDD INSTRUMENTS dropped for lack of covariate data: their univariable estimate, which is the
   representativeness check the reviewer asked for.
5. WHAT THE PER-SD LIABILITY CONVERSION WOULD NEED: the observed-scale figures the paper currently uses, and
   the range implied by plausible liability-scale h2_SNP, to justify removing the sentence.

Writes mvmr_independent/mvmr_paired_v8.csv.
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

DIR = os.path.dirname(os.path.abspath(__file__))
OUT = f'{DIR}/mvmr_paired_v8.csv'
M = pd.read_csv(os.path.join(DIR, 'mvmr_matrix_harmonised.csv'))
S1 = pd.read_csv(os.path.join('data', 'S1_MDD_to_Anxiety_harmonised_198_instruments.csv'))
S2 = pd.read_csv('data/S2_Anxiety_to_MDD_harmonised_50_instruments.csv')
COV = ['beta_BMI', 'beta_Education', 'beta_SmokingInit']
CSE = ['se_BMI', 'se_Education', 'se_SmokingInit']

SETS = {
    'forward | exposure-only union': dict(xb='beta_MDD', yb='beta_Anxiety', ys='se_Anxiety',
                                          mask=lambda d: d.is_Anx_iv == 0),
    'forward | covariates from exposure instruments': dict(
        xb='beta_MDD', yb='beta_Anxiety', ys='se_Anxiety',
        mask=lambda d: (d.is_MDD_iv == 1) & (d.indep_BMI == 0) & (d.indep_Education == 0) & (d.indep_SmokingInit == 0)),
    'reverse | exposure-only union': dict(xb='beta_Anxiety', yb='beta_MDD', ys='se_MDD',
                                          mask=lambda d: d.is_MDD_iv == 0),
    'reverse | covariates from exposure instruments': dict(
        xb='beta_Anxiety', yb='beta_MDD', ys='se_MDD',
        mask=lambda d: (d.is_Anx_iv == 1) & (d.indep_BMI == 0) & (d.indep_Education == 0) & (d.indep_SmokingInit == 0)),
}


def mvmr(rows, xb, yb, ys, correct=True):
    X = rows[[xb] + COV].values
    seX = rows[[xb.replace('beta', 'se')] + CSE].values
    y = rows[yb].values
    sey = rows[ys].values
    w = 1 / sey ** 2
    XtWX = X.T @ (X * w[:, None])
    A = np.zeros((4, 4))
    if correct:
        for j in range(len(rows)):
            A += w[j] * np.diag(seX[j] ** 2)
    Cinv = np.linalg.inv(XtWX - A)
    b = Cinv @ (X.T @ (y * w))
    r = y - X @ b
    disp = np.sqrt(np.sum(w * r ** 2) / (len(rows) - 4))
    se = np.sqrt(np.diag(Cinv)) * (max(1.0, np.sqrt(disp)) if correct else 1.0)
    se_plain = np.sqrt(np.diag(np.linalg.inv(XtWX))) * disp
    cov_scaled = Cinv * max(1.0, disp)
    bb = b[1:4]
    Cb = cov_scaled[np.ix_(range(1, 4), range(1, 4))]
    chi2 = float(bb @ np.linalg.solve(Cb, bb))
    Fstat = chi2 / 3.0
    return dict(beta=b[0], se=se[0], disp=disp, se_pkg=se_plain[0], beta_pkg=(
        np.linalg.solve(XtWX, X.T @ (y * w)))[0],
        wald_chi2=chi2, wald_p=stats.chi2.sf(chi2, 3), wald_F=Fstat,
        wald_p_F=stats.f.sf(Fstat, 3, len(rows) - 4))


def uvmr(rows, xb, yb, ys, correct=True):
    x = rows[xb].values
    sex = rows[xb.replace('beta', 'se')].values
    y = rows[yb].values
    sey = rows[ys].values
    w = 1 / sey ** 2
    den = np.sum(w * x ** 2) - (np.sum(w * sex ** 2) if correct else 0.0)
    b = np.sum(w * x * y) / den
    r = y - b * x
    q = np.sum(w * r ** 2)
    se = np.sqrt(1 / den) * max(1.0, np.sqrt(q / (len(rows) - 1)))
    return b, se, q / (len(rows) - 1)


rows_out = []
print('=' * 112)
print('1. paired comparison, both sides under the same estimator')
for name, spec in SETS.items():
    rows = M[spec['mask'](M)]
    uc, usc, qc = uvmr(rows, spec['xb'], spec['yb'], spec['ys'], True)
    up, usp, qp = uvmr(rows, spec['xb'], spec['yb'], spec['ys'], False)
    mc = mvmr(rows, spec['xb'], spec['yb'], spec['ys'], True)
    mp = mvmr(rows, spec['xb'], spec['yb'], spec['ys'], False)
    rows_out.append(dict(model=name, n=len(rows),
                         uni_corr_beta=uc, uni_corr_se=usc, uni_corr_or=np.exp(uc),
                         uni_corr_lo=np.exp(uc - 1.96 * usc), uni_corr_hi=np.exp(uc + 1.96 * usc),
                         uni_plain_beta=up, uni_plain_se=usp, uni_plain_or=np.exp(up),
                         uni_plain_lo=np.exp(up - 1.96 * usp), uni_plain_hi=np.exp(up + 1.96 * usp),
                         dir_corr_beta=mc['beta'], dir_corr_se=mc['se'], dir_corr_or=np.exp(mc['beta']),
                         dir_corr_lo=np.exp(mc['beta'] - 1.96 * mc['se']), dir_corr_hi=np.exp(mc['beta'] + 1.96 * mc['se']),
                         dir_plain_beta=mp['beta'], dir_plain_se=mp['se_pkg'], dir_plain_or=np.exp(mp['beta']),
                         dir_plain_lo=np.exp(mp['beta'] - 1.96 * mp['se_pkg']), dir_plain_hi=np.exp(mp['beta'] + 1.96 * mp['se_pkg']),
                         ratio_corr=mc['beta'] / uc, ratio_plain=mp['beta'] / up,
                         ratio_v7=mc['beta'] / up, dispersion=mc['disp'],
                         wald_chi2_3df=mc['wald_chi2'], wald_p=mc['wald_p'],
                         wald_F=mc['wald_F'], wald_p_F=mc['wald_p_F'],
                         uni_q_per_df=qc))
    print(f'  {name:48s} n={len(rows):3d}')
    print(f'      corrected : univariable {uc:+.4f} ({usc:.4f}) OR {np.exp(uc):.3f} | direct {mc["beta"]:+.4f} '
          f'({mc["se"]:.4f}) OR {np.exp(mc["beta"]):.3f} | ratio {mc["beta"]/uc:.3f}')
    print(f'      plain     : univariable {up:+.4f} ({usp:.4f}) OR {np.exp(up):.3f} | direct {mp["beta"]:+.4f} '
          f'({mp["se_pkg"]:.4f}) OR {np.exp(mp["beta"]):.3f} | ratio {mp["beta"]/up:.3f}')
    print(f'      v7 mixed  : ratio {mc["beta"]/up:.3f}   covariate joint test: chi2(3)={mc["wald_chi2"]:.2f} '
          f'P={mc["wald_p"]:.3g}  F(3,{len(rows)-4})={mc["wald_F"]:.2f} P={mc["wald_p_F"]:.3g}')

print('\n' + '=' * 112)
print('3. covariate instruments in LD range of a primary instrument, and the refit without them')
prim = M[(M.is_MDD_iv == 1) | (M.is_Anx_iv == 1)][['snp', 'chrom', 'pos']]
covr = M[(M.indep_BMI == 1) | (M.indep_Education == 1) | (M.indep_SmokingInit == 1)][['snp', 'chrom', 'pos']].copy()
mind = []
for _, c in covr.iterrows():
    p = prim[prim.chrom == c.chrom]['pos'].values
    mind.append(min(np.abs(p - c['pos'])) if len(p) else np.inf)
covr['dist'] = mind
near1 = set(covr.loc[covr.dist < 1_000_000, 'snp'])
near01 = set(covr.loc[covr.dist < 100_000, 'snp'])
print(f'  covariate instruments within 1 Mb of a primary instrument: {len(near1)} of {len(covr)}; within 100 kb: {len(near01)}')
print(f'  closest: {sorted(zip(covr.snp, covr.dist.astype(int)), key=lambda x: x[1])[:4]}')
for lab, drop in [('drop the 100 kb pairs', near01), ('drop all within 1 Mb', near1)]:
    for name, spec in [('forward | exposure-only union', SETS['forward | exposure-only union']),
                       ('reverse | exposure-only union', SETS['reverse | exposure-only union'])]:
        rows = M[spec['mask'](M) & ~M.snp.isin(drop)]
        uc, usc, _ = uvmr(rows, spec['xb'], spec['yb'], spec['ys'], True)
        mc = mvmr(rows, spec['xb'], spec['yb'], spec['ys'], True)
        print(f'  {lab:24s} {name:32s} n={len(rows):3d} univariable {np.exp(uc):.3f} direct {np.exp(mc["beta"]):.3f} '
              f'ratio {mc["beta"]/uc:.3f}')
        rows_out.append(dict(model=f'{name} | {lab}', n=len(rows), uni_corr_beta=uc, uni_corr_se=usc,
                             uni_corr_or=np.exp(uc), uni_corr_lo=np.exp(uc - 1.96 * usc),
                             uni_corr_hi=np.exp(uc + 1.96 * usc), dir_corr_beta=mc['beta'], dir_corr_se=mc['se'],
                             dir_corr_or=np.exp(mc['beta']), dir_corr_lo=np.exp(mc['beta'] - 1.96 * mc['se']),
                             dir_corr_hi=np.exp(mc['beta'] + 1.96 * mc['se']), ratio_corr=mc['beta'] / uc,
                             dispersion=mc['disp'], wald_chi2_3df=mc['wald_chi2'], wald_p=mc['wald_p']))

print('\n' + '=' * 112)
print('4. the 114 MDD instruments excluded for lack of covariate data')
kept = set(M[M.is_MDD_iv == 1]['snp'])
for lab, sub in [('retained 84', S1[S1.SNP.isin(kept)]), ('excluded 114', S1[~S1.SNP.isin(kept)]),
                 ('all 198', S1)]:
    b, se, q = uvmr(sub.assign(beta_MDD=sub.Beta_exposure, se_MDD=sub.SE_exposure,
                               beta_Anxiety=sub.Beta_outcome, se_Anxiety=sub.SE_outcome),
                    'beta_MDD', 'beta_Anxiety', 'se_Anxiety', False)
    bc, sec, _ = uvmr(sub.assign(beta_MDD=sub.Beta_exposure, se_MDD=sub.SE_exposure,
                                 beta_Anxiety=sub.Beta_outcome, se_Anxiety=sub.SE_outcome),
                      'beta_MDD', 'beta_Anxiety', 'se_Anxiety', True)
    eaf = pd.to_numeric(sub['EAF_exposure'], errors='coerce')
    print(f'  {lab:12s} n={len(sub):3d}  plain OR {np.exp(b):.3f} ({np.exp(b-1.96*se):.3f}-{np.exp(b+1.96*se):.3f})'
          f'  corrected OR {np.exp(bc):.3f}  EAF median {eaf.median():.3f}  Q/df {q:.2f}')
    rows_out.append(dict(model=f'MDD univariable on the {lab} instruments', n=len(sub),
                         uni_plain_or=np.exp(b), uni_plain_lo=np.exp(b - 1.96 * se), uni_plain_hi=np.exp(b + 1.96 * se),
                         uni_corr_or=np.exp(bc), uni_q_per_df=q))

print('\n' + '=' * 112)
print('5. the per-SD-of-genetic-liability conversion')
for lab, h2f, h2r in [('observed-scale R2 as used in v7', 0.0044, 0.0022),
                      ('liability h2_SNP 0.02 / 0.02', 0.02, 0.02),
                      ('liability h2_SNP 0.05 / 0.03', 0.05, 0.03),
                      ('liability h2_SNP 0.05 / 0.05', 0.05, 0.05)]:
    vf, vr = 0.8516 * np.sqrt(h2f), 0.5170 * np.sqrt(h2r)
    print(f'  {lab:34s} forward {vf:.4f}  reverse {vr:.4f}  ratio {vf/vr:.2f}')
print('  no liability-scale h2_SNP was computed in this project, so the figure cannot be made defensible;')
print('  the sentence is removed from the abstract and Section 3.1 rather than restated with an assumption.')

out = pd.DataFrame(rows_out)
out.to_csv(OUT, index=False)
print(f'\nwrote {OUT} ({len(out)} rows)')
