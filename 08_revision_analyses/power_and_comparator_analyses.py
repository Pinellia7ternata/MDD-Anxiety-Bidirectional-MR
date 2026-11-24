"""Reviewer response 1, 3 and 4 in one runnable script.

Produces three CSVs used by the revised manuscript:

  coloc_power_experiment.csv   - how coloc's class assignment moves with outcome-GWAS power
  table5_composition_matched.csv - the deletion comparator matched to the RETAINED set's F
                                   composition, replacing the previous proportional-decile column
  winners_curse_terciles.csv   - unfiltered versus filtered estimates within exposure-F terciles

The power experiment is the one that changes the manuscript's central claim: with the anxiety
standard errors halved, 68% of the instruments currently at H1 loci are reclassified as H4, so
H1 cannot be read as evidence of trait-specificity.
"""
import csv, math, os, random, statistics, sys
from _locate import locate, output

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from coloc_power_simulation import coloc_single, load, S1, COLOC  # noqa: E402

REPO = os.environ.get('MR_OUT_DIR') or '.'
S2 = locate('S2_Anxiety_to_MDD_harmonised_50_instruments.csv')
POWER_CSV = f'{REPO}/coloc_power_experiment.csv'
T5_CSV = f'{REPO}/table5_composition_matched.csv'
WC_CSV = f'{REPO}/winners_curse_terciles.csv'
NREP = 2000
SEED = 20261102


def ivw(d):
    den = sum(x['bx'] ** 2 / x['sy'] ** 2 for x in d)
    b = sum(x['bx'] * x['by'] / x['sy'] ** 2 for x in d) / den
    se_fe = math.sqrt(1 / den)
    q = sum((x['by'] - b * x['bx']) ** 2 / x['sy'] ** 2 for x in d)
    df = len(d) - 1
    se = se_fe * math.sqrt(max(q / df, 1.0))
    return dict(b=b, se=se, or_=math.exp(b), lo=math.exp(b - 1.96 * se), hi=math.exp(b + 1.96 * se),
                q_per=q / df, n=len(d))


def published_classes():
    loci = list(csv.DictReader(open(COLOC, encoding='utf-8-sig')))
    for l in loci:
        pp = {k: float(l[k]) for k in ('PP.H0', 'PP.H1', 'PP.H2', 'PP.H3', 'PP.H4')}
        l['cls'] = max(pp, key=pp.get).replace('PP.', '')
        l['chr'] = int(l['CHR'])
        l['bp'] = int(l['Lead_BP'])

    def f(x):
        c = [l for l in loci if l['chr'] == x['chr'] and abs(l['bp'] - x['pos']) < 500000]
        return max(c, key=lambda l: float(l['PP.H4']))['cls'] if c else 'none'
    return f


def power_experiment(d, cls_of):
    rows = []
    for src, target in (('H4', 'H1'), ('H1', 'H4')):
        g = [x for x in d if cls_of(x) == src]
        factors = (sorted([1.0, 1.5, 2.0, 2.4, 3.0], reverse=True) if src == 'H4'
                   else sorted([1.0, 0.7, 0.5, 0.35, 0.25], reverse=True))
        for f in factors:
            moved = sum(1 for x in g if coloc_single(x['bx'], x['sx'], x['by'], x['sy'] * f)[0] == target)
            rows.append(dict(source_class=src, n=len(g), outcome_se_multiplier=f,
                             reclassified_to=target, n_reclassified=moved,
                             percent=round(100 * moved / len(g), 1)))
            print(f'  {src} (n={len(g)}) with outcome SE x{f:4.2f} -> {moved:3d} reclassified as {target} '
                  f'({100*moved/len(g):5.1f}%)')
    with open(POWER_CSV, 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(f'-> {POWER_CSV}')


def composition_matched(d, shared, direction, thr):
    inside = lambda x: any(c == x['chr'] and abs(p - x['pos']) < 500000 for c, p in shared)
    keep = [x for x in d if not inside(x)]
    k_del = len(d) - len(keep)
    obs = ivw(keep)
    base = ivw(d)
    meanF = lambda g: statistics.mean((x['bx'] / x['sx']) ** 2 for x in g)
    print(f'\n  {direction}: delete {k_del} of {len(d)} | mean F full {meanF(d):.1f} '
          f'retained {meanF(keep):.1f} deleted {meanF([x for x in d if inside(x)]):.1f}')

    def stratified(pool, n_target, rng):
        ranked = sorted(pool, key=lambda x: (x['bx'] / x['sx']) ** 2)
        bins = [[] for _ in range(10)]
        for i, x in enumerate(ranked):
            bins[min(i * 10 // len(ranked), 9)].append(x)
        tgt = sorted(keep, key=lambda x: (x['bx'] / x['sx']) ** 2)
        quota = [0] * 10
        for i, x in enumerate(tgt):
            quota[min(i * 10 // len(tgt), 9)] += 1
        sub = []
        for b in range(10):
            if bins[b]:
                sub.extend(rng.sample(bins[b], min(quota[b], len(bins[b]))))
        while len(sub) < n_target:
            extra = [x for x in pool if x not in sub]
            if not extra:
                break
            sub.append(rng.choice(extra))
        return sub[:n_target]

    rng = random.Random(SEED)
    qs, ors, qr, orr = [], [], [], []
    for _ in range(NREP):
        r = ivw(rng.sample(d, len(keep)))
        qs.append(r['q_per']); ors.append(r['or_'])
        r2 = ivw(stratified(d, len(keep), rng))
        qr.append(r2['q_per']); orr.append(r2['or_'])

    def summ(vals, obs_v):
        m, sd = statistics.mean(vals), statistics.pstdev(vals)
        pct = 100 * sum(1 for v in vals if v <= obs_v) / len(vals)
        return m, sd, (obs_v - m) / sd, pct

    out = []
    for stat, vals, obs_v in (('Q_per_instrument', qs, obs['q_per']), ('IVW_OR', ors, obs['or_']),
                              ('Q_per_instrument', qr, obs['q_per']), ('IVW_OR', orr, obs['or_'])):
        comp = 'uniform random' if vals is qs or vals is ors else 'retained-composition matched'
        m, sd, z, pct = summ(vals, obs_v)
        print(f'    {comp:28s} {stat:17s} null {m:.3f} (SD {sd:.3f})  observed {obs_v:.3f}  z={z:+.2f}  percentile {pct:.1f}%')
        out.append(dict(coloc_thresh=thr, direction=direction, comparator=comp, statistic=stat, n_deleted=k_del,
                        n_retained=len(keep), null_mean=round(m, 4), null_sd=round(sd, 4),
                        observed=round(obs_v, 4), z=round(z, 2), percentile_of_observed=round(pct, 2),
                        n_rep=len(vals), baseline_value=round(base['q_per'] if stat.startswith('Q') else base['or_'], 4)))
    return out


def winners_curse(d, shared, direction):
    inside = lambda x: any(c == x['chr'] and abs(p - x['pos']) < 500000 for c, p in shared)
    ranked = sorted(d, key=lambda x: (x['bx'] / x['sx']) ** 2)
    terc = [ranked[:len(ranked) // 3], ranked[len(ranked) // 3:2 * len(ranked) // 3], ranked[2 * len(ranked) // 3:]]
    rows = []
    print(f'  {direction}: estimates within terciles of exposure-side F')
    for i, t in enumerate(terc, 1):
        a = ivw(t)
        b = ivw([x for x in t if not inside(x)])
        fmean = statistics.mean((x['bx'] / x['sx']) ** 2 for x in t)
        print(f'    tercile {i} (mean F {fmean:6.1f}): all n={a["n"]:3d} OR {a["or_"]:.3f} | '
              f'non-coloc n={b["n"]:3d} OR {b["or_"]:.3f}  ratio {b["or_"]/a["or_"]:.3f}')
        rows.append(dict(direction=direction, f_tercile=i, mean_F=round(fmean, 1),
                         n_all=a['n'], or_all=round(a['or_'], 3), ci_all=f"{a['lo']:.3f}-{a['hi']:.3f}",
                         q_all=round(a['q_per'], 2), n_noncoloc=b['n'], or_noncoloc=round(b['or_'], 3),
                         ci_noncoloc=f"{b['lo']:.3f}-{b['hi']:.3f}", q_noncoloc=round(b['q_per'], 2),
                         or_ratio_noncoloc_over_all=round(b['or_'] / a['or_'], 3)))
    return rows


def main():
    fwd = load(S1)
    rev = load(S2)
    loci = list(csv.DictReader(open(COLOC, encoding='utf-8-sig')))
    shared = [(int(l['CHR']), int(l['Lead_BP'])) for l in loci if float(l['PP.H4']) > 0.5]

    print('== power experiment (reviewer point 1) ==')
    power_experiment(fwd, published_classes())

    print('\n== composition-matched deletion comparators (reviewer point 3) ==')
    t5 = []
    for thr in (0.5, 0.95):
        sh = [(int(l['CHR']), int(l['Lead_BP'])) for l in loci if float(l['PP.H4']) > thr]
        print(f'  --- PP.H4 > {thr} ---')
        t5 += composition_matched(fwd, sh, 'MDD -> anxiety', thr) + composition_matched(rev, sh, 'Anxiety -> MDD', thr)
    with open(T5_CSV, 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=list(t5[0].keys()))
        w.writeheader(); w.writerows(t5)
    print(f'-> {T5_CSV}')

    print('\n== winner\'s-curse tercile analysis (reviewer point 4) ==')
    wc = winners_curse(fwd, shared, 'MDD -> anxiety') + winners_curse(rev, shared, 'Anxiety -> MDD')
    with open(WC_CSV, 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=list(wc[0].keys()))
        w.writeheader(); w.writerows(wc)
    print(f'-> {WC_CSV}')


if __name__ == '__main__':
    main()
