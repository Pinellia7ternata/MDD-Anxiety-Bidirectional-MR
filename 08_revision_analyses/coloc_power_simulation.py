"""Reviewer response 1: is coloc's H1 class evidence of trait-specificity, or of power?

The revised manuscript builds its central asymmetry claim on the 75 MDD instruments sitting at
H1 loci ("no independent anxiety signal"). But the anxiety GWAS is 2.4x smaller than the MDD
GWAS, so a genuinely shared variant can be assigned H1 simply because its anxiety effect is
not detected. This script quantifies that.

Step 1  calibration - recompute the single-SNP coloc.abf posterior from the published lead-SNP
        statistics and measure how well it reproduces the classes the manuscript reports.
Step 2  null of complete sharing - simulate anxiety effects under the hypothesis that EVERY
        instrument acts through one shared causal variant with the fitted forward slope
        (beta_y = 0.8516 x beta_x), using each instrument's actual anxiety standard error,
        and record how often coloc.abf returns H1.
Step 3  the number the manuscript needs - expected H1 count under complete sharing versus the
        observed 75, i.e. how many H1 loci can plausibly be called trait-specific.

Also reports the analytic check: the anxiety-side z that a fully shared model predicts at H1
loci versus the z actually observed there.
"""
import csv, math, random, statistics
from _locate import locate, output

S1 = locate('S1_MDD_to_Anxiety_harmonised_198_instruments.csv')
COLOC = locate('coloc_locus_level_strict.csv')
OUT = output('coloc_power_simulation.csv')
P1, P2, P12 = 1e-4, 1e-4, 1e-5
W = 0.2 ** 2          # coloc default sd.prior for case-control traits
NSIM = 200
SEED = 20261102


def logabf(b, se):
    V = se ** 2
    z = b / se
    return 0.5 * math.log(V / (V + W)) + 0.5 * z ** 2 * W / (V + W)


def coloc_single(bx, sx, by, sy):
    l1, l2 = logabf(bx, sx), logabf(by, sy)
    # H3 uses the product of the two single-trait sums; H4 the joint sum
    terms = {'H0': 0.0, 'H1': math.log(P1) + l1, 'H2': math.log(P2) + l2,
             'H3': math.log(P1 * P2) + l1 + l2, 'H4': math.log(P12) + l1 + l2}
    m = max(terms.values())
    tot = sum(math.exp(v - m) for v in terms.values())
    pp = {k: math.exp(v - m) / tot for k, v in terms.items()}
    return max(pp, key=pp.get), pp


def load(path):
    out = []
    for r in csv.DictReader(open(path, encoding='utf-8-sig')):
        if r.get('Included_in_MR', 'True').strip().lower() != 'true':
            continue
        out.append(dict(snp=r['SNP'], chr=int(r['Chr']), pos=int(r['Position']),
                        bx=float(r['Beta_exposure']), sx=float(r['SE_exposure']),
                        by=float(r['Beta_outcome']), sy=float(r['SE_outcome'])))
    return out


def main():
    d = load(S1)
    loci = list(csv.DictReader(open(COLOC, encoding='utf-8-sig')))
    for l in loci:
        pp = {k: float(l[k]) for k in ('PP.H0', 'PP.H1', 'PP.H2', 'PP.H3', 'PP.H4')}
        l['cls'] = max(pp, key=pp.get).replace('PP.', '')
        l['chr'] = int(l['CHR'])
        l['bp'] = int(l['Lead_BP'])

    def published(x):
        c = [l for l in loci if l['chr'] == x['chr'] and abs(l['bp'] - x['pos']) < 500000]
        return max(c, key=lambda l: float(l['PP.H4']))['cls'] if c else 'none'

    # forward IVW slope = the shared-model prediction coefficient
    den = sum(x['bx'] ** 2 / x['sy'] ** 2 for x in d)
    slope = sum(x['bx'] * x['by'] / x['sy'] ** 2 for x in d) / den
    print(f'forward IVW slope (shared-model predictor) = {slope:.4f}\n')

    print('== step 1: calibration of the single-SNP coloc.abf against the published classes ==')
    agree = {}
    for x in d:
        pub = published(x)
        sim, _ = coloc_single(x['bx'], x['sx'], x['by'], x['sy'])
        agree.setdefault(pub, []).append(sim == pub)
    overall = []
    for pub, flags in sorted(agree.items()):
        overall.extend(flags)
        print(f'  published {pub:8s} n={len(flags):3d}  reproduced by single-SNP coloc.abf: {100*sum(flags)/len(flags):5.1f}%')
    print(f'  overall agreement: {100*sum(overall)/len(overall):.1f}%')

    obs = {}
    for x in d:
        obs[published(x)] = obs.get(published(x), 0) + 1
    print(f'  observed class counts: {dict(sorted(obs.items()))}')

    print('\n== step 2: expected classes if EVERY variant is shared (complete-sharing null) ==')
    rng = random.Random(SEED)
    counts = {}
    h1_given_shared = []
    for _ in range(NSIM):
        c = {'H1': 0, 'H2': 0, 'H3': 0, 'H4': 0, 'H0': 0}
        for x in d:
            by_sim = rng.gauss(slope * x['bx'], x['sy'])
            cls, _ = coloc_single(x['bx'], x['sx'], by_sim, x['sy'])
            c[cls] += 1
        for k in c:
            counts.setdefault(k, []).append(c[k])
    print(f'  simulated {NSIM} replicates, 198 instruments each')
    rows = []
    for k in ('H1', 'H3', 'H4', 'H2', 'H0'):
        if k not in counts:
            continue
        m = statistics.mean(counts[k])
        sd = statistics.pstdev(counts[k])
        lo, hi = sorted(counts[k])[int(0.025 * NSIM)], sorted(counts[k])[int(0.975 * NSIM)]
        o = obs.get(k, 0)
        print(f'  {k}: expected under complete sharing {m:6.1f} (SD {sd:4.1f}, 95% range {lo}-{hi})   observed {o:3d}')
        rows.append(dict(cls=k, expected_shared=round(m, 1), sd=round(sd, 2), lo95=lo, hi95=hi, observed=o))
    h1mean = statistics.mean(counts['H1'])
    print(f'\n  => of the {obs.get("H1",0)} instruments classed H1, {h1mean:.0f} are expected '
          f'({100*h1mean/obs.get("H1",1):.0f}% of them) even under COMPLETE sharing.')
    print(f'     instruments at H1 loci not expected under sharing: {max(obs.get("H1",0) - h1mean, 0):.0f} of 198')

    print('\n== step 3: analytic check on the H1 instruments themselves ==')
    h1 = [x for x in d if published(x) == 'H1']
    h4 = [x for x in d if published(x) == 'H4']
    for name, g in (('H1', h1), ('H4', h4)):
        zo = [abs(x['by'] / x['sy']) for x in g]
        zp = [abs(slope * x['bx'] / x['sy']) for x in g]
        nom = sum(1 for x in g if abs(x['by'] / x['sy']) > 1.96) / len(g)
        print(f'  {name}: n={len(g):3d} observed anxiety-side |z| mean {statistics.mean(zo):.2f} '
              f'| shared-model predicted |z| mean {statistics.mean(zp):.2f} '
              f'| ratio observed/predicted {statistics.mean(zo)/statistics.mean(zp):.2f} '
              f'| nominally significant {100*nom:.0f}%')
    print('  (a fully specific instrument would show observed |z| ~1.0 and ~5% nominal significance)')

    with open(OUT, 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=['cls', 'expected_shared', 'sd', 'lo95', 'hi95', 'observed'])
        w.writeheader()
        w.writerows(rows)
    print(f'\n-> {OUT}')


if __name__ == '__main__':
    main()
