"""Corrected re-estimation of the bidirectional MR statistics.

Audits and replaces the numbers that failed review in the rejected JAD version:
  * variance explained was reported as sum(2*pq*beta^2) with the denominator dropped
  * "MR-RAPS" and "Radial MR" came from a script using SE = 1/sqrt(sum F)
  * coloc-filtered (locus-specific) instruments were never analysed

Only estimators that reproduce the manuscript's verified source values are used here:
IVW under a multiplicative random-effects model (reproduces Table 2 exactly), Cochran Q,
and a per-instrument Q screen. Egger / weighted-median / mode rows are left as published
because this script's own implementations do not reproduce TwoSampleMR's values and must
not silently replace them.

Run from the project root:
    python JAD_Submission_Package_MDD_Anxiety_v2/08_Revision_20261102/recompute_corrected_estimates.py
"""
import csv, math
from _locate import locate, output

CSV = {
    'fwd': locate('S1_MDD_to_Anxiety_harmonised_198_instruments.csv'),
    'rev': locate('S2_Anxiety_to_MDD_harmonised_50_instruments.csv'),
}
COLOC = locate('coloc_locus_level_strict.csv')
OUT = output('coloc_filtered_MR.csv')

N = {'MDD': (412305, 1588397), 'ANX': (122083, 729602)}
DIR = {'fwd': ('MDD', 'ANX'), 'rev': ('ANX', 'MDD')}
WINDOW = 500000
CHISQ_1_DF_95 = 3.841459


def log10_p(z):
    """log10 of a two-sided normal p-value, stable for large |z|."""
    z = abs(z)
    if z < 8:
        return math.log10(max(2 * 0.5 * (1 + math.erf(-z / math.sqrt(2))), 1e-320))
    return (math.log10(math.sqrt(2 / math.pi)) - math.log10(z) - z * z / (2 * math.log(10)) + math.log10(2))


def load(path):
    out = []
    with open(path, encoding='utf-8-sig') as fh:
        for r in csv.DictReader(fh):
            if r.get('Included_in_MR', 'True').strip().lower() != 'true':
                continue
            out.append(dict(snp=r['SNP'], chr=int(r['Chr']), pos=int(r['Position']),
                            bx=float(r['Beta_exposure']), sx=float(r['SE_exposure']),
                            by=float(r['Beta_outcome']), sy=float(r['SE_outcome']),
                            eaf=float(r['EAF_exposure'])))
    return out


def ivw(d):
    n = len(d)
    den = sum(x['bx'] ** 2 / x['sy'] ** 2 for x in d)
    b = sum(x['bx'] * x['by'] / x['sy'] ** 2 for x in d) / den
    se_fe = math.sqrt(1 / den)
    qi = [(x['by'] - b * x['bx']) ** 2 / x['sy'] ** 2 for x in d]
    q = sum(qi)
    df = n - 1
    se_re = se_fe * math.sqrt(max(q / df, 1.0))
    return dict(n=n, b=b, se=se_re, se_fe=se_fe, q=q, df=df, infl=se_re / se_fe,
                or_=math.exp(b), lo=math.exp(b - 1.96 * se_re), hi=math.exp(b + 1.96 * se_re),
                log10p=log10_p(b / se_re), n_outlier_q95=sum(1 for v in qi if v > CHISQ_1_DF_95))


def r2(beta, se, eaf, n_total):
    """Approximate variance explained: 2pq*b^2 / (2pq*b^2 + N*2pq*se^2)  (quantitative-scale approximation)."""
    p = 2 * eaf * (1 - eaf)
    vx = p * beta ** 2
    ve = p * se ** 2
    den = vx + n_total * ve
    return vx / den if den > 0 else 0.0


def variance_and_steiger(d, ex, ou):
    nex, nou = N[ex][0] + N[ex][1], N[ou][0] + N[ou][1]
    rx = [r2(x['bx'], x['sx'], x['eaf'], nex) for x in d]
    ry = [r2(x['by'], x['sy'], x['eaf'], nou) for x in d]
    zs = []
    for a, bb in zip(rx, ry):
        r1 = math.sqrt(min(a, 0.999999))
        r2v = math.sqrt(min(bb, 0.999999))
        f = 0.5 * math.log((1 + r1) / (1 - r1)) - 0.5 * math.log((1 + r2v) / (1 - r2v))
        n = min(nex, nou)
        zs.append(f / math.sqrt(2 / (n - 3)))
    k = len(zs)
    mz = sum(zs) / k
    sd = math.sqrt(sum((z - mz) ** 2 for z in zs) / (k - 1))
    z = mz / (sd / math.sqrt(k)) if sd > 0 else float('inf')
    return dict(sum_r2_ex=100 * sum(rx), sum_r2_ou=100 * sum(ry), mean_r2_ex=1000 * sum(rx) / k,
                n_correct=sum(1 for v in zs if v > 0), k=k, z=z, log10p=log10_p(z))


def fmt(name, i, s=None):
    txt = (f'{name}: n={i["n"]:3d}  beta={i["b"]:.4f} (SE {i["se"]:.4f})  '
           f'OR={i["or_"]:.3f} ({i["lo"]:.3f}-{i["hi"]:.3f})  log10P={i["log10p"]:.1f}  '
           f'Q={i["q"]:.1f}/{i["df"]} = {i["q"]/i["df"]:.2f}/SNP  Q_i>3.84: {i["n_outlier_q95"]}')
    if s:
        txt += f'\n     R2 exposure={s["sum_r2_ex"]:.3f}%  outcome={s["sum_r2_ou"]:.3f}%  ' \
               f'per-SNP R2x>R2y: {s["n_correct"]}/{s["k"]}  Z={s["z"]:.2f} log10P={s["log10p"]:.1f}'
    return txt


def main():
    data = {t: load(p) for t, p in CSV.items()}
    loci = list(csv.DictReader(open(COLOC, encoding='utf-8-sig')))
    print('=== 0) coloc locus counts (source: coloc_locus_level_strict.csv) ===')
    for th in (0.5, 0.9, 0.95, 0.99, 0.995):
        print(f'  PP.H4 > {th}: {sum(1 for l in loci if float(l["PP.H4"]) > th)} of {len(loci)} loci')

    print('\n=== 1) primary estimates (must reproduce Table 2) + per-instrument Q screen ===')
    rows = []
    prim = {}
    for t, d in data.items():
        ex, ou = DIR[t]
        i = ivw(d)
        s = variance_and_steiger(d, ex, ou)
        prim[t] = i
        print(fmt(f'{t} ALL', i, s))
        rows.append(dict(coloc_thresh='none', direction=t, instrument_set='all',
                         n_instruments=i['n'], n_total=len(d),
                         ivw_beta=round(i['b'], 4), ivw_se=round(i['se'], 4),
                         ivw_or=round(i['or_'], 4), ivw_lo=round(i['lo'], 4), ivw_hi=round(i['hi'], 4),
                         ivw_p_log10=round(i['log10p'], 1),
                         Q=round(i['q'], 1), df=i['df'], Q_per_snp=round(i['q'] / i['df'], 2),
                         n_q_outliers=i['n_outlier_q95'],
                         sum_r2_exposure_pct=round(s['sum_r2_ex'], 3),
                         sum_r2_outcome_pct=round(s['sum_r2_ou'], 3),
                         steiger_z=round(s['z'], 2), steiger_log10p=round(s['log10p'], 1)))

    print('\n=== 2) locus-specific instrument sets (shared coloc loci removed) ===')
    for thr in (0.5, 0.95):
        L = [(int(l['CHR']), int(l['Lead_BP'])) for l in loci if float(l['PP.H4']) > thr]
        for t, d in data.items():
            ex, ou = DIR[t]
            specific = [x for x in d if not any(c == x['chr'] and abs(p - x['pos']) < WINDOW for c, p in L)]
            shared = [x for x in d if any(c == x['chr'] and abs(p - x['pos']) < WINDOW for c, p in L)]
            for lab, sub in (('non-shared', specific), ('shared-only', shared)):
                if len(sub) < 5:
                    print(f'  PP.H4>{thr} {t} {lab}: n={len(sub)} - not analysable')
                    continue
                i = ivw(sub)
                s = variance_and_steiger(sub, ex, ou)
                print(fmt(f'  PP.H4>{thr} {t} [{lab}]', i, s))
                rows.append(dict(coloc_thresh=thr, direction=t, instrument_set=lab,
                                 n_instruments=i['n'], n_total=len(d),
                                 ivw_beta=round(i['b'], 4), ivw_se=round(i['se'], 4),
                                 ivw_or=round(i['or_'], 4), ivw_lo=round(i['lo'], 4), ivw_hi=round(i['hi'], 4),
                                 ivw_p_log10=round(i['log10p'], 1),
                                 Q=round(i['q'], 1), df=i['df'], Q_per_snp=round(i['q'] / i['df'], 2),
                                 n_q_outliers=i['n_outlier_q95'],
                                 sum_r2_exposure_pct=round(s['sum_r2_ex'], 3),
                                 sum_r2_outcome_pct=round(s['sum_r2_ou'], 3),
                                 steiger_z=round(s['z'], 2), steiger_log10p=round(s['log10p'], 1)))
    with open(OUT, 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f'\n-> {OUT}')


if __name__ == '__main__':
    main()
