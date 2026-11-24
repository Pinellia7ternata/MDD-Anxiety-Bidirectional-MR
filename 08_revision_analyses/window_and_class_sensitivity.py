"""Reviewer response 2: window sensitivity, coloc hypothesis-class decomposition, and
scale-normalised directional estimates.

Answers reviewer M3 (the +-500 kb window is arbitrary; loci with two distinct causal
variants (H3) are retained as "locus-specific" although they are the configuration where
MR is most easily biased) and M4 (forward and reverse ORs are on different exposure scales
and cannot be compared as printed).

Hypothesis classes follow Giambartolomei et al.: H1 exposure-only signal, H2 outcome-only,
H3 both traits associated through distinct causal variants, H4 shared causal variant.
"""
import csv, math
from _locate import locate, output

CSV = {
    'fwd': (locate('S1_MDD_to_Anxiety_harmonised_198_instruments.csv'), 'MDD', 'ANX'),
    'rev': (locate('S2_Anxiety_to_MDD_harmonised_50_instruments.csv'), 'ANX', 'MDD'),
}
COLOC = locate('coloc_locus_level_strict.csv')
N = {'MDD': (412305, 1588397), 'ANX': (122083, 729602)}
OUT = output('window_and_class_sensitivity.csv')


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
    den = sum(x['bx'] ** 2 / x['sy'] ** 2 for x in d)
    b = sum(x['bx'] * x['by'] / x['sy'] ** 2 for x in d) / den
    se_fe = math.sqrt(1 / den)
    q = sum((x['by'] - b * x['bx']) ** 2 / x['sy'] ** 2 for x in d)
    df = len(d) - 1
    se = se_fe * math.sqrt(max(q / df, 1.0))
    return dict(n=len(d), b=b, se=se, or_=math.exp(b), lo=math.exp(b - 1.96 * se),
                hi=math.exp(b + 1.96 * se), q_per=q / df)


def r2(x, trait):
    nc, no = N[trait]
    Nt = nc + no
    p = 2 * x['eaf'] * (1 - x['eaf'])
    vx = p * x['bx'] ** 2
    return vx / (vx + Nt * p * x['sx'] ** 2)


def main():
    loci = list(csv.DictReader(open(COLOC, encoding='utf-8-sig')))
    for l in loci:
        pp = {k: float(l[k]) for k in ('PP.H0', 'PP.H1', 'PP.H2', 'PP.H3', 'PP.H4')}
        l['class'] = max(pp, key=pp.get).replace('PP.', '')
        l['pp_h4'] = pp['PP.H4']
        l['chr'] = int(l['CHR'])
        l['bp'] = int(l['Lead_BP'])
    rows = []
    for tag, (path, ex, ou) in CSV.items():
        d = load(path)
        print(f'\n===== {tag} (exposure {ex}, {len(d)} instruments) =====')
        for win in (250000, 500000, 1000000):
            for thr in (0.5, 0.95):
                shared = [(l['chr'], l['bp']) for l in loci if l['pp_h4'] > thr]
                keep = [x for x in d if not any(c == x['chr'] and abs(p - x['pos']) < win for c, p in shared)]
                if len(keep) < 5:
                    continue
                r = ivw(keep)
                sdg = math.sqrt(sum(r2(x, ex) for x in keep))
                print(f'  window +-{win/1000:.0f}kb  PP.H4>{thr}: n={r["n"]:3d} OR={r["or_"]:.3f} '
                      f'({r["lo"]:.3f}-{r["hi"]:.3f}) Q/instr={r["q_per"]:.2f}  '
                      f'beta per SD of genetic liability={r["b"]*sdg:.4f}')
                rows.append(dict(analysis=f'window_+-{win//1000}kb_H4>{thr}', direction=tag, n=r['n'],
                                 beta=round(r['b'], 4), se=round(r['se'], 4), or_=round(r['or_'], 3),
                                 lo=round(r['lo'], 3), hi=round(r['hi'], 3), q_per_instrument=round(r['q_per'], 3),
                                 sum_r2_exposure_pct=round(100 * sdg ** 2, 3),
                                 beta_per_SD_genetic_liability=round(r['b'] * sdg, 4)))
        # hypothesis-class decomposition at +-500 kb
        def cls_of(x):
            cand = [l for l in loci if l['chr'] == x['chr'] and abs(l['bp'] - x['pos']) < 500000]
            return max(cand, key=lambda l: l['pp_h4'])['class'] if cand else 'no_locus'
        groups = {}
        for x in d:
            groups.setdefault(cls_of(x), []).append(x)
        print('  coloc hypothesis class of each instrument\'s locus (+-500 kb):')
        for k in sorted(groups, key=lambda g: -len(groups[g])):
            g = groups[k]
            if len(g) < 5:
                print(f'    {k:9s} n={len(g):3d}  (too few to analyse)')
                continue
            r = ivw(g)
            sdg = math.sqrt(sum(r2(x, ex) for x in g))
            print(f'    {k:9s} n={len(g):3d}  OR={r["or_"]:.3f} ({r["lo"]:.3f}-{r["hi"]:.3f}) '
                  f'Q/instr={r["q_per"]:.2f}  beta per SD={r["b"]*sdg:.4f}')
            rows.append(dict(analysis=f'class_{k}', direction=tag, n=r['n'], beta=round(r['b'], 4),
                             se=round(r['se'], 4), or_=round(r['or_'], 3), lo=round(r['lo'], 3),
                             hi=round(r['hi'], 3), q_per_instrument=round(r['q_per'], 3),
                             sum_r2_exposure_pct=round(100 * sdg ** 2, 3),
                             beta_per_SD_genetic_liability=round(r['b'] * sdg, 4)))
        # the reviewer-relevant strict set: exposure-specific only.
        # coloc was run MDD as trait 1 and anxiety as trait 2 (columns P_MDD / n_Anxiety),
        # so H1 = MDD-only signal and H2 = anxiety-only signal; the exposure-specific class
        # therefore depends on the direction.
        exp_only = 'H1' if ex == 'MDD' else 'H2'
        strict = [x for x in d if cls_of(x) in (exp_only, 'no_locus')]
        untested = [x for x in d if cls_of(x) == 'no_locus']
        h3 = [x for x in d if cls_of(x) == 'H3']
        print(f'  instruments with no coloc locus within +-500 kb (region untested): {len(untested)}')
        if len(strict) >= 5:
            r = ivw(strict)
            print(f'  STRICT exposure-specific ({exp_only} or no coloc locus): n={r["n"]} OR={r["or_"]:.3f} '
                  f'({r["lo"]:.3f}-{r["hi"]:.3f}) Q/instr={r["q_per"]:.2f}')
            rows.append(dict(analysis='strict_exposure_specific', direction=tag, n=r['n'], beta=round(r['b'], 4),
                             se=round(r['se'], 4), or_=round(r['or_'], 3), lo=round(r['lo'], 3),
                             hi=round(r['hi'], 3), q_per_instrument=round(r['q_per'], 3)))
        if len(h3) >= 5:
            r = ivw(h3)
            print(f'  H3 loci (both traits, DISTINCT causal variants): n={r["n"]} OR={r["or_"]:.3f} '
                  f'({r["lo"]:.3f}-{r["hi"]:.3f}) Q/instr={r["q_per"]:.2f}   <-- retained by the current filter')
            rows.append(dict(analysis='H3_distinct_variants', direction=tag, n=r['n'], beta=round(r['b'], 4),
                             se=round(r['se'], 4), or_=round(r['or_'], 3), lo=round(r['lo'], 3),
                             hi=round(r['hi'], 3), q_per_instrument=round(r['q_per'], 3)))
    with open(OUT, 'w', newline='', encoding='utf-8') as fh:
        keys = ['analysis', 'direction', 'n', 'beta', 'se', 'or_', 'lo', 'hi', 'q_per_instrument',
                'sum_r2_exposure_pct', 'beta_per_SD_genetic_liability']
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, '') for k in keys})
    print(f'\n-> {OUT}')


if __name__ == '__main__':
    main()
