"""Reviewer response 1: null controls for the colocalization-filtered analysis.

Reviewer M2: "removing the 90 instruments that colocalize removes exactly the instruments
with the largest Q_i, so any such deletion lowers residual heterogeneity. Provide a null
control: does random or Q-ordered deletion of the same number of instruments also drive
Q per instrument to 1.18?"

Four comparators per direction and per coloc threshold, all deleting the same number of
instruments as the colocalization filter:
  random          - 2000 uniformly random subsets (seeded)
  F-matched       - random subsets drawn within deciles of exposure-side F, so the
                    removed set matches the strength profile of the colocalizing set
  Q-greedy        - iteratively delete the single instrument with the largest Q_i
                    (the most aggressive heterogeneity-reducing deletion possible)
  coloc (observed)- the analysis reported in Table 4

Outputs both IVW OR and Q per instrument so the reviewer can see whether the effect
attenuation, the heterogeneity drop, or both are non-generic.
"""
import csv, math
from _locate import locate, output, random, statistics

CSV = {
    'fwd': locate('S1_MDD_to_Anxiety_harmonised_198_instruments.csv'),
    'rev': locate('S2_Anxiety_to_MDD_harmonised_50_instruments.csv'),
}
COLOC = locate('coloc_locus_level_strict.csv')
OUT = output('locus_specific_null_controls.csv')
WINDOW = 500000
NREP = 2000
SEED = 20261102


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
    qi = [(x['by'] - b * x['bx']) ** 2 / x['sy'] ** 2 for x in d]
    q = sum(qi)
    df = len(d) - 1
    se = se_fe * math.sqrt(max(q / df, 1.0))
    return dict(b=b, se=se, or_=math.exp(b), lo=math.exp(b - 1.96 * se), hi=math.exp(b + 1.96 * se),
                q_per=df and q / df, qi=qi, n=len(d))


def q_greedy(d, k):
    cur = list(d)
    for _ in range(k):
        if len(cur) < 5:
            break
        r = ivw(cur)
        j = max(range(len(cur)), key=lambda i: r['qi'][i])
        cur = cur[:j] + cur[j + 1:]
    return cur


def stratified_sample(d, k, rng):
    """Sample k instruments for deletion, stratified by decile of exposure-side F."""
    ranked = sorted(d, key=lambda x: (x['bx'] / x['sx']) ** 2)
    bins = [[] for _ in range(10)]
    for i, x in enumerate(ranked):
        bins[min(i * 10 // len(ranked), 9)].append(x)
    drop = []
    quota = {b: 0 for b in range(10)}
    counts = {b: len(bins[b]) for b in range(10)}
    total = len(d)
    for b in range(10):
        quota[b] = round(k * counts[b] / total)
    # fix rounding drift
    drift = k - sum(quota.values())
    order = sorted(range(10), key=lambda b: -counts[b])
    i = 0
    while drift != 0 and i < 1000:
        b = order[i % 10]
        if drift > 0 and quota[b] < counts[b]:
            quota[b] += 1; drift -= 1
        elif drift < 0 and quota[b] > 0:
            quota[b] -= 1; drift += 1
        i += 1
    for b in range(10):
        drop.extend(rng.sample(bins[b], min(quota[b], counts[b])))
    ids = {x['snp'] for x in drop}
    return [x for x in d if x['snp'] not in ids]


def summarize(label, results, observed, stat):
    vals = [r for r in results if r is not None]
    if not vals:
        return None
    mean = statistics.mean(vals)
    sd = statistics.pstdev(vals) if len(vals) > 1 else 0
    if len(vals) == 1:  # deterministic comparator: no sampling distribution to compare against
        print(f'    {label:26s} value={mean:.3f}  observed={observed:.3f}  '
              f'({"observed is lower" if observed < mean else "observed is higher"} by {abs(mean-observed):.3f})')
        return dict(comparator=label, statistic=stat, mean=round(mean, 4), sd='', observed=round(observed, 4),
                    z='', empirical_p='', n_rep=1)
    below = sum(1 for v in vals if v <= observed)
    above = sum(1 for v in vals if v >= observed)
    p = min(2 * min(below, above) / len(vals), 1.0)
    z = (observed - mean) / sd if sd else float('nan')
    print(f'    {label:26s} mean={mean:.3f} sd={sd:.3f}  observed={observed:.3f}  '
          f'z={z:+.2f}  two-sided empirical P={p:.4g}')
    return dict(comparator=label, statistic=stat, mean=round(mean, 4), sd=round(sd, 4),
                observed=round(observed, 4), z=round(z, 2), empirical_p=round(p, 6), n_rep=len(vals))


def main():
    rng = random.Random(SEED)
    loci = list(csv.DictReader(open(COLOC, encoding='utf-8-sig')))
    rows = []
    for thr in (0.5, 0.95):
        shared = [(int(l['CHR']), int(l['Lead_BP'])) for l in loci if float(l['PP.H4']) > thr]
        for tag, path in CSV.items():
            d = load(path)
            is_shared = lambda x: any(c == x['chr'] and abs(p - x['pos']) < WINDOW for c, p in shared)
            coloc_keep = [x for x in d if not is_shared(x)]
            k = len(d) - len(coloc_keep)
            if k < 1 or len(coloc_keep) < 5:
                continue
            obs = ivw(coloc_keep)
            base = ivw(d)
            print(f'\n=== PP.H4>{thr} {tag}: delete {k} of {len(d)} (Q/instr {base["q_per"]:.2f} -> {obs["q_per"]:.2f}; '
                  f'OR {base["or_"]:.3f} -> {obs["or_"]:.3f}) ===')
            res_q, res_or, res_fq, res_for = [], [], [], []
            for _ in range(NREP):
                sub = rng.sample(d, len(coloc_keep))
                r = ivw(sub)
                res_q.append(r['q_per']); res_or.append(r['or_'])
            for _ in range(NREP):
                sub = stratified_sample(d, k, rng)
                r = ivw(sub)
                res_fq.append(r['q_per']); res_for.append(r['or_'])
            greedy = ivw(q_greedy(d, k))
            print('    Q per instrument:')
            a = summarize('random subsets', res_q, obs['q_per'], 'Q_per_instrument')
            b = summarize('F-stratified random', res_fq, obs['q_per'], 'Q_per_instrument')
            c = summarize('Q-greedy deletion', [greedy['q_per']], obs['q_per'], 'Q_per_instrument')
            print('    IVW OR:')
            e = summarize('random subsets', res_or, obs['or_'], 'IVW_OR')
            f = summarize('F-stratified random', res_for, obs['or_'], 'IVW_OR')
            g = summarize('Q-greedy deletion', [greedy['or_']], obs['or_'], 'IVW_OR')
            for r in (a, b, c, e, f, g):
                if r:
                    rows.append(dict(coloc_thresh=thr, direction=tag, n_total=len(d), n_deleted=k,
                                     n_retained=len(coloc_keep), baseline_Q=round(base['q_per'], 3),
                                     baseline_OR=round(base['or_'], 3), **r))
    with open(OUT, 'w', newline='', encoding='utf-8') as fh:
        keys = ['coloc_thresh', 'direction', 'n_total', 'n_deleted', 'n_retained', 'baseline_Q',
                'baseline_OR', 'statistic', 'comparator', 'mean', 'sd', 'observed', 'z',
                'empirical_p', 'n_rep']
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, '') for k in keys})
    print(f'\n-> {OUT}')


if __name__ == '__main__':
    main()
