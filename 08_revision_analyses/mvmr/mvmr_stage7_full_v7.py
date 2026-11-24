"""Stage 7: every coefficient of the four v7 multivariable models, under both estimators.

Table 3 of the manuscript reports the direct estimate; Supplementary Table S7 has to carry the whole model. For
each of the four fits (forward and reverse, exposure-only instrument set and covariates-from-the-exposure
instrument set) this writes the coefficient, standard error and P value for all four exposures under

  pkg  - the dispersion-scaled weighted least squares that ivw_mvmr in MVMR 0.4.8 returns, and
  sw   - the Sanderson-Windmeijer measurement-error-corrected estimator, with and without the residual
         dispersion applied to its standard errors (the manuscript reports the scaled form).

Conditional F-statistics come from stage 5 (strength_mvmr). The palindromic count and the estimate with those
variants removed are recorded here too, because the reviewer asked how strand ambiguity was handled in a matrix
that was aligned to a single reference allele.

Writes mvmr_independent/mvmr_full_models_v7.csv.
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

DIR = os.path.dirname(os.path.abspath(__file__))
M = pd.read_csv(f'{DIR}/mvmr_matrix_harmonised.csv')
ST5 = pd.read_csv(f'{DIR}/mvmr_exposure_only.csv')
CMP = pd.read_csv(f'{DIR}/mvmr_matched_comparison.csv')
COV = ['beta_BMI', 'beta_Education', 'beta_SmokingInit']
CSE = ['se_BMI', 'se_Education', 'se_SmokingInit']
PAL = {frozenset(x) for x in [('A', 'T'), ('C', 'G')]}

al = {}
for f in ['BMI_ieu_a_2_on_union_instruments.csv', 'Education_ieu_a_80_on_union_instruments.csv',
          'SmokingInit_ieu_b_4877_on_union_instruments.csv']:
    for _, r in pd.read_csv(f'{DIR}/{f}').iterrows():
        al[r['snp']] = (str(r['ref']).upper(), str(r['alt']).upper())
for f, lab in [('BMI_ieu_a_2', 'BMI'), ('Education_ieu_a_80', 'education'), ('SmokingInit_ieu_b_4877', 'smoking')]:
    d = pd.read_csv(f'{DIR}/{f}_on_primary_instruments.csv')
    for _, r in d.iterrows():
        al.setdefault(r['snp'], (str(r['ref']).upper(), str(r['alt']).upper()))
M['palindromic'] = M['snp'].map(lambda x: frozenset(al[x]) in PAL if x in al else np.nan)

SETS = {
    'forward | exposure-only union': dict(xb='beta_MDD', yb='beta_Anxiety', ys='se_Anxiety',
                                          mask=lambda d: d.is_Anx_iv == 0, names=['MDD', 'BMI', 'education', 'smoking initiation']),
    'forward | covariates from exposure instruments': dict(
        xb='beta_MDD', yb='beta_Anxiety', ys='se_Anxiety',
        mask=lambda d: (d.is_MDD_iv == 1) & (d.indep_BMI == 0) & (d.indep_Education == 0) & (d.indep_SmokingInit == 0),
        names=['MDD', 'BMI', 'education', 'smoking initiation']),
    'reverse | exposure-only union': dict(xb='beta_Anxiety', yb='beta_MDD', ys='se_MDD',
                                          mask=lambda d: d.is_MDD_iv == 0,
                                          names=['anxiety disorders', 'BMI', 'education', 'smoking initiation']),
    'reverse | covariates from exposure instruments': dict(
        xb='beta_Anxiety', yb='beta_MDD', ys='se_MDD',
        mask=lambda d: (d.is_Anx_iv == 1) & (d.indep_BMI == 0) & (d.indep_Education == 0) & (d.indep_SmokingInit == 0),
        names=['anxiety disorders', 'BMI', 'education', 'smoking initiation']),
}


def fits(rows, xb, yb, ys):
    X = rows[[xb] + COV].values
    seX = rows[[xb.replace('beta', 'se')] + CSE].values
    y = rows[yb].values
    sey = rows[ys].values
    w = 1 / sey ** 2
    XtWX = X.T @ (X * w[:, None])
    b_pkg = np.linalg.solve(XtWX, X.T @ (y * w))
    r = y - X @ b_pkg
    disp = np.sqrt(np.sum(w * r ** 2) / (len(rows) - X.shape[1]))
    se_pkg = np.sqrt(np.diag(np.linalg.inv(XtWX))) * disp
    A = np.zeros((X.shape[1], X.shape[1]))
    for j in range(len(rows)):
        A += w[j] * np.diag(seX[j] ** 2)
    C = np.linalg.inv(XtWX - A)
    b_sw = C @ (X.T @ (y * w))
    se_sw = np.sqrt(np.diag(C))
    r2 = y - X @ b_sw
    disp2 = np.sqrt(np.sum(w * r2 ** 2) / (len(rows) - X.shape[1]))
    return dict(b_pkg=b_pkg, se_pkg=se_pkg, disp_pkg=disp,
                b_sw=b_sw, se_sw=se_sw, se_sw_scaled=se_sw * max(1.0, np.sqrt(disp2)), disp_sw=disp2)


rows_out = []
for name, spec in SETS.items():
    rows = M[spec['mask'](M)]
    f = fits(rows, spec['xb'], spec['yb'], spec['ys'])
    st5 = ST5[ST5['set'] == name].iloc[0]
    cf = [st5['condF_primary'], st5['condF_BMI'], st5['condF_Education'], st5['condF_Smoking']]
    nop = rows[rows.palindromic != True]
    fp = fits(nop, spec['xb'], spec['yb'], spec['ys'])
    for i, exp in enumerate(spec['names']):
        t = f['b_sw'][i] / f['se_sw_scaled'][i]
        rows_out.append(dict(
            model=name, exposure=exp, n=len(rows),
            beta_pkg=round(f['b_pkg'][i], 4), se_pkg=round(f['se_pkg'][i], 4),
            beta_sw=round(f['b_sw'][i], 4), se_sw=round(f['se_sw'][i], 4),
            se_sw_scaled=round(f['se_sw_scaled'][i], 4),
            or_sw=round(np.exp(f['b_sw'][i]), 4),
            lo=round(np.exp(f['b_sw'][i] - 1.96 * f['se_sw_scaled'][i]), 4),
            hi=round(np.exp(f['b_sw'][i] + 1.96 * f['se_sw_scaled'][i]), 4),
            p_sw=round(2 * stats.t.sf(abs(t), len(rows) - 4), 8),
            conditional_F=round(cf[i], 2),
            dispersion_pkg=round(f['disp_pkg'], 3), dispersion_sw=round(f['disp_sw'], 3),
            palindromic=int((rows.palindromic == True).sum()),
            beta_sw_no_palindromic=round(fp['b_sw'][i], 4),
            or_sw_no_palindromic=round(np.exp(fp['b_sw'][i]), 4),
        ))

out = pd.DataFrame(rows_out)
out.to_csv(f'{DIR}/mvmr_full_models_v7.csv', index=False)
pd.set_option('display.width', 260)
print(out[['model', 'exposure', 'n', 'beta_sw', 'se_sw_scaled', 'or_sw', 'lo', 'hi', 'p_sw',
           'conditional_F', 'beta_pkg', 'se_pkg', 'palindromic', 'or_sw_no_palindromic']].to_string(index=False))
print('\ndispersion (pkg / sw):')
print(out.groupby('model')[['dispersion_pkg', 'dispersion_sw']].first().to_string())
print('\nstage-5 package primary estimates reproduced here:',
      np.allclose(out[out.exposure.isin(['MDD', 'anxiety disorders'])]['beta_pkg'].values,
                  ST5['beta'].values, atol=2e-4))
