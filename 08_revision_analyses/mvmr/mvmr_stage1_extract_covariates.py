"""Stage 1 of the independent-instrument MVMR requested by the reviewer.

The current MVMR instruments the covariate exposures from the primary-exposure instrument set.
This stage builds the alternative: instruments selected independently in each covariate GWAS.

Covariate summary statistics on disk (影响因素/):
  ieu-a-2      BMI (GIANT 2018)
  ieu-a-80     years of education (Lee 2018 / Ramin 2017)
  ieu-b-4877   smoking initiation (GSCAN)
  Smit_GIANT_BMI_2025_PRSCS   a newer BMI GWAS, plain text

GWAS Catalog VCFs pack the statistics in the FORMAT field ES:SE:LP:AF:SS:ID, where LP is
-log10(P). Everything is written to CSV so each step can be inspected before the next runs.
Clumping here is POSITION based (>= 1 Mb apart, lowest P first): no PLINK binary exists on
this machine, so the 1000G panel in 数据/00_reference cannot be used. That limitation is
recorded in the output and must be stated in the manuscript.
"""
import csv, gzip, json, math, os

BASE = os.environ.get('MR_SUMSTATS_DIR', '')
OUTDIR = os.path.dirname(os.path.abspath(__file__))
S1 = os.path.join('data', 'S1_MDD_to_Anxiety_harmonised_198_instruments.csv')
S2 = os.path.join('data', 'S2_Anxiety_to_MDD_harmonised_50_instruments.csv')

COVS = {
    'BMI_ieu_a_2': f'{BASE}/影响因素/ieu-a-2.vcf.gz',
    'Education_ieu_a_80': f'{BASE}/影响因素/ieu-a-80.vcf.gz',
    'SmokingInit_ieu_b_4877': f'{BASE}/影响因素/ieu-b-4877.vcf.gz',
}
P_THRESHOLD_LP = -math.log10(5e-8)          # 7.301
MIN_SEP = 1_000_000



if not BASE or not os.path.isdir(BASE):
    raise SystemExit('set MR_SUMSTATS_DIR to the directory holding the PGC and IEU OpenGWAS summary statistics; they are not redistributed')

def targets():
    out = {}
    for f, lab in ((S1, 'MDD'), (S2, 'Anxiety')):
        for r in csv.DictReader(open(f, encoding='utf-8-sig')):
            out[r['SNP']] = lab
    return out


def scan_vcf(path, wanted):
    """Return (stats for wanted SNPs, all genome-wide-significant rows)."""
    hits, gws = {}, []
    n = 0
    with gzip.open(path, 'rt', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 10:
                continue
            n += 1
            chrom, pos, rsid, ref, alt = f[0], f[1], f[2], f[3], f[4]
            fmt, vals = f[8], f[9]
            keys = fmt.split(':')
            if keys[:4] != ['ES', 'SE', 'LP', 'AF']:
                continue
            v = vals.split(':')
            try:
                es, se, lp = float(v[0]), float(v[1]), float(v[2])
            except (ValueError, IndexError):
                continue
            rec = dict(snp=rsid, chrom=chrom, pos=pos, ref=ref, alt=alt, beta=es, se=se,
                       lp=lp, af=(float(v[3]) if len(v) > 3 else ''), n=(v[4] if len(v) > 4 else ''))
            if rsid in wanted:
                hits[rsid] = rec
            if lp >= P_THRESHOLD_LP:
                gws.append(rec)
    return n, hits, gws


def clump_by_position(rows):
    rows = sorted(rows, key=lambda r: -float(r['lp']))
    by_chr = {}
    kept = []
    for r in rows:
        c = r['chrom']
        try:
            p = int(r['pos'])
        except ValueError:
            continue
        if any(abs(p - q) < MIN_SEP for q in by_chr.setdefault(c, [])):
            continue
        by_chr[c].append(p)
        kept.append(r)
    return kept


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    want = targets()
    print(f'target instrument SNPs: {len(want)} (MDD {sum(1 for v in want.values() if v=="MDD")}, '
          f'anxiety {sum(1 for v in want.values() if v=="Anxiety")})')
    summary = {}
    for name, path in COVS.items():
        if not os.path.exists(path):
            print('MISSING', path)
            continue
        n, hits, gws = scan_vcf(path, want)
        kept = clump_by_position(gws)
        with open(f'{OUTDIR}/{name}_gws_clumped.csv', 'w', newline='', encoding='utf-8') as fh:
            w = csv.DictWriter(fh, fieldnames=['snp', 'chrom', 'pos', 'ref', 'alt', 'beta', 'se', 'lp', 'af', 'n'])
            w.writeheader()
            w.writerows(kept)
        with open(f'{OUTDIR}/{name}_on_primary_instruments.csv', 'w', newline='', encoding='utf-8') as fh:
            w = csv.DictWriter(fh, fieldnames=['snp', 'chrom', 'pos', 'ref', 'alt', 'beta', 'se', 'lp', 'af', 'n'])
            w.writeheader()
            w.writerows([hits[k] for k in sorted(hits)])
        summary[name] = dict(variants_scanned=n, gws=len(gws), gws_clumped=len(kept),
                             primary_instruments_with_data=len(hits))
        print(f'  {name:24s} scanned {n:>9,} | P<5e-8: {len(gws):>6,} | position-clumped: {len(kept):>4,} '
              f'| stats found for {len(hits):>3}/248 primary instruments')
    with open(f'{OUTDIR}/stage1_summary.json', 'w', encoding='utf-8') as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)
    print(f'\n-> {OUTDIR}/stage1_summary.json')
    print('NOTE: clumping is position based (no PLINK on this machine); record this in the manuscript.')


if __name__ == '__main__':
    main()
