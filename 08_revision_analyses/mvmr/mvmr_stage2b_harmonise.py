"""Stage 2b: harmonised MVMR matrix, all five traits aligned to ONE reference allele.

Stage 2 stored each trait's beta in that trait's own effect allele, which is not a harmonised
multivariable input - it is why the old-design control run failed to reproduce Table 3.
Here every beta is flipped onto the effect allele of the MDD GWAS, the outcome of the reverse
model and one of the exposures of the forward model; variants whose alleles cannot be matched
are dropped and counted rather than guessed at.
"""
import csv, gzip, os

BASE = os.environ.get('MR_SUMSTATS_DIR', '')
MD = f'{BASE}/数据/01_formatted/MDD_PGC_2025_EUR_no23andMe_formatted.tsv.gz'
AX = f'{BASE}/数据/01_formatted/Anxiety_PGC_2026_EUR_formatted.tsv.gz'
DIR = os.path.dirname(os.path.abspath(__file__))
S1 = os.path.join('data', 'S1_MDD_to_Anxiety_harmonised_198_instruments.csv')
S2 = os.path.join('data', 'S2_Anxiety_to_MDD_harmonised_50_instruments.csv')
COVS = {'BMI': 'BMI_ieu_a_2', 'Education': 'Education_ieu_a_80', 'SmokingInit': 'SmokingInit_ieu_b_4877'}



if not BASE or not os.path.isdir(BASE):
    raise SystemExit('set MR_SUMSTATS_DIR to the directory holding the PGC and IEU OpenGWAS summary statistics; they are not redistributed')

def read_gwas(path, wanted):
    out = {}
    with gzip.open(path, 'rt', encoding='utf-8', errors='replace') as fh:
        hdr = fh.readline().rstrip('\n').split('\t')
        idx = {c: i for i, c in enumerate(hdr)}
        for line in fh:
            f = line.rstrip('\n').split('\t')
            snp = f[idx['SNP']]
            if snp in wanted:
                out[snp] = dict(ea=f[idx['effect_allele']].upper(), oa=f[idx['other_allele']].upper(),
                                beta=float(f[idx['beta']]), se=float(f[idx['se']]),
                                chrom=f[idx['chr']], pos=f[idx['pos']])
    return out


def main():
    prim = set()
    for f in (S1, S2):
        for r in csv.DictReader(open(f, encoding='utf-8-sig')):
            prim.add(r['SNP'])
    mdd_iv = {r['SNP'] for r in csv.DictReader(open(S1, encoding='utf-8-sig'))}
    anx_iv = {r['SNP'] for r in csv.DictReader(open(S2, encoding='utf-8-sig'))}

    cov, indep = {}, {}
    for short, name in COVS.items():
        with open(f'{DIR}/{name}_on_union_instruments.csv', encoding='utf-8') as fh:
            cov[short] = {r['snp']: r for r in csv.DictReader(fh)}
        with open(f'{DIR}/{name}_gws_clumped.csv', encoding='utf-8') as fh:
            indep[short] = {r['snp'] for r in csv.DictReader(fh)}
    wanted = prim | {s for d in indep.values() for s in d}

    mdd = read_gwas(MD, wanted)
    anx = read_gwas(AX, wanted)
    ref = {s: mdd[s]['ea'] for s in mdd}          # align everything to the MDD effect allele

    rows = []
    tally = dict(no_mdd=0, no_anx=0, allele_mismatch=0, no_cov=0, ok=0)
    for snp in sorted(wanted):
        m = mdd.get(snp)
        if m is None:
            tally['no_mdd'] += 1
            continue
        a = anx.get(snp)
        if a is None:
            tally['no_anx'] += 1
            continue
        R = ref[snp]

        def align(beta, ea, oa):
            """return beta on allele R, or None if the alleles cannot be matched"""
            if ea == R and oa != R:
                return beta
            if oa == R and ea != R:
                return -beta
            return None

        av = align(a['beta'], a['ea'], a['oa'])
        if av is None:
            tally['allele_mismatch'] += 1
            continue
        rec = dict(snp=snp, chrom=m['chrom'], pos=m['pos'], ref_allele=R,
                   beta_MDD=m['beta'], se_MDD=m['se'],
                   beta_Anxiety=av, se_Anxiety=a['se'],
                   is_MDD_iv=int(snp in mdd_iv), is_Anx_iv=int(snp in anx_iv))
        bad = False
        for short, name in COVS.items():
            r = cov[short].get(snp)
            if r is None:
                tally['no_cov'] += 1
                bad = True
                break
            alt, other = r['alt'].upper(), r['ref'].upper()
            cv = align(float(r['beta']), alt, other)
            if cv is None:
                tally['allele_mismatch'] += 1
                bad = True
                break
            rec[f'beta_{short}'] = cv
            rec[f'se_{short}'] = float(r['se'])
            rec[f'indep_{short}'] = int(snp in indep[short])
        if bad:
            continue
        tally['ok'] += 1
        rows.append(rec)

    cols = ['snp', 'chrom', 'pos', 'ref_allele', 'beta_MDD', 'se_MDD', 'beta_Anxiety', 'se_Anxiety',
            'is_MDD_iv', 'is_Anx_iv'] + [k for s in COVS for k in (f'beta_{s}', f'se_{s}', f'indep_{s}')]
    with open(f'{DIR}/mvmr_matrix_harmonised.csv', 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    print('usable instruments:', len(rows), '| dropped:', tally)
    for s in COVS:
        print(f'  {s}: {sum(r[f"indep_{s}"] for r in rows)} independently selected instruments retained')
    print(f'  MDD instruments present: {sum(r["is_MDD_iv"] for r in rows)} / 198')
    print(f'  anxiety instruments present: {sum(r["is_Anx_iv"] for r in rows)} / 50')
    print(f'-> {DIR}/mvmr_matrix_harmonised.csv')


if __name__ == '__main__':
    main()
