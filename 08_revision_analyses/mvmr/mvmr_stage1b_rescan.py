"""Stage 1b: re-scan the covariate GWASs for the FULL union instrument set.

Stage 1 only recorded covariate statistics for the 246 primary instruments plus each
covariate's own hits, so a BMI-selected instrument had no education or smoking coefficient
and multivariable IVW, which needs a complete exposure matrix per instrument, discarded all
of them. This pass reads each file again against the union list.
"""
import csv, gzip, os

BASE = os.environ.get('MR_SUMSTATS_DIR', '')
DIR = os.path.dirname(os.path.abspath(__file__))
COVS = {'BMI_ieu_a_2': 'ieu-a-2', 'Education_ieu_a_80': 'ieu-a-80', 'SmokingInit_ieu_b_4877': 'ieu-b-4877'}



if not BASE or not os.path.isdir(BASE):
    raise SystemExit('set MR_SUMSTATS_DIR to the directory holding the PGC and IEU OpenGWAS summary statistics; they are not redistributed')

def union_snps():
    s = set()
    for c in COVS:
        with open(f'{DIR}/{c}_gws_clumped.csv', encoding='utf-8') as fh:
            s |= {r['snp'] for r in csv.DictReader(fh)}
        with open(f'{DIR}/{c}_on_primary_instruments.csv', encoding='utf-8') as fh:
            s |= {r['snp'] for r in csv.DictReader(fh)}
    return s


def scan(path, wanted):
    out = {}
    n = 0
    with gzip.open(path, 'rt', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            n += 1
            if len(f) < 10 or f[2] not in wanted:
                continue
            keys = f[8].split(':')
            if keys[:3] != ['ES', 'SE', 'LP']:
                continue
            v = f[9].split(':')
            try:
                es, se = float(v[0]), float(v[1])
            except ValueError:
                continue
            out[f[2]] = dict(snp=f[2], chrom=f[0], pos=f[1], ref=f[3], alt=f[4],
                             beta=es, se=se, lp=v[2], af=v[3] if len(v) > 3 else '',
                             n=v[4] if len(v) > 4 else '')
    return n, out


def main():
    wanted = union_snps()
    print(f'union set to look up in every covariate GWAS: {len(wanted)} SNPs')
    for c, ds in COVS.items():
        n, hits = scan(f'{BASE}/影响因素/{ds}.vcf.gz', wanted)
        with open(f'{DIR}/{c}_on_union_instruments.csv', 'w', newline='', encoding='utf-8') as fh:
            w = csv.DictWriter(fh, fieldnames=['snp', 'chrom', 'pos', 'ref', 'alt', 'beta', 'se', 'lp', 'af', 'n'])
            w.writeheader()
            w.writerows([hits[k] for k in sorted(hits)])
        print(f'  {c:24s} scanned {n:>10,} | found {len(hits):>4}/{len(wanted)}')


if __name__ == '__main__':
    main()
