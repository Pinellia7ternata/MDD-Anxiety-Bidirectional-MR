# Multivariable MR with covariate instruments selected independently

Round-4 revision analyses. The reviewer comment that started this: the MVMR covariates (BMI, education years,
smoking initiation) were instrumented by the primary-exposure instruments, so their conditional F-statistics were
3.3, 1.0 and 3.5 and the direct estimates could not be read as effects of separate pathways. The fix is to
instrument each covariate in its own GWAS.

## Two corrections came out of doing it

**The previous Table 3 was not reproducible.** On the same 84 instruments multivariable IVW gives a forward
direct beta of 0.829 (SE 0.057), not the 0.144 (SE 0.029) reported. Three variance treatments were tested - the
plain weighted least squares that `ivw_mvmr` in MVMR 0.4.8 performs (weights 1 / se_Y^2, standard errors scaled by
residual dispersion), the same fit with the per-SNP exposure covariance taken as diagonal, and with it taken as the
outer product of the standard errors - giving 0.829, 0.768 and 1.339. None is 0.144, and the covariate statistics
match the source files variant by variant. The table was withdrawn rather than reinterpreted.

**A union instrument set that includes outcome-selected variants inflates the direct estimate.** The first
independent-instrument design used the union of the exposure and covariate instruments, 264 variants, and that set
contains the 27 variants selected because they are genome-wide significant for anxiety disorders and the 84
selected for MDD. Each of those appears as an instrument in the model where the trait that selected it is the
outcome, and outcome-side winner's curse pulls a multivariable slope up. The forward direct estimate came out at
2.82 against a univariable 2.34, which was read as BMI masking the association; on the same 264 variants the
univariable estimate is already 2.79, so the increase was the sample. Stages 5-7 exclude outcome-selected variants
and pair every direct estimate with the univariable estimate on exactly the same instruments.

## Result (the numbers now in the manuscript)

| model | n | univariable IVW on the same variants | direct (corrected, dispersion-scaled) | ratio | conditional F: primary / BMI / educ / smoke |
|---|---|---|---|---|---|
| forward, exposure-only union | 237 | 2.41 (2.24-2.58) | 2.57 (2.40-2.76) | 1.08 | 12.1 / 16.6 / 8.1 / 8.7 |
| forward, covariates from exposure instruments | 84 | 2.44 (2.26-2.64) | 2.16 (1.82-2.56) | 0.86 | 5.0 / 3.3 / 1.0 / 3.5 |
| reverse, exposure-only union | 180 | 1.63 (1.55-1.72) | 1.71 (1.64-1.79) | 1.10 | 7.6 / 18.9 / 10.3 / 11.2 |
| reverse, covariates from exposure instruments | 27 | 1.65 (1.55-1.75) | 1.64 (1.53-1.76) | 1.00 | 12.9 / 2.5 / 1.5 / 6.6 |

Conditioning on BMI, education and smoking initiation moves neither estimate by more than a tenth on the
instruments it uses. The forward BMI coefficient is negative in both designs (-0.067 and -0.081), which is why the
direct estimate is marginally the larger one forward. Manuscript locations: Table 3, Figure 3, Supplementary
Table S7.

Diagnostics recorded in `mvmr_matched_comparison.csv`: requiring data on all three covariate GWAS retains 84 of
the 198 MDD instruments because BMI and education each cover 89 of them (smoking initiation covers 188); 16 of the
237 forward variants and 12 of the 180 reverse variants are strand-ambiguous A/T or C/G pairs, and excluding them
moves the forward direct estimate from 2.57 to 2.55 and the reverse from 1.71 to 1.72. In the
exposure-instrumented models the covariate coefficients are not identified - education has conditional F 1.0
forward, its corrected coefficient is OR 5.65 (0.83-38.6), and removing five strand-ambiguous variants changes it
by three orders of magnitude - which is why they are reported but not interpreted.

## Stages

| stage | script | input | output |
|---|---|---|---|
| 1 | `mvmr_stage1_extract_covariates.py` | ieu-a-2, ieu-a-80, ieu-b-4877 VCFs | `*_gws_clumped.csv`, `*_on_primary_instruments.csv`, `stage1_summary.json` |
| 1b | `mvmr_stage1b_rescan.py` | the same VCFs | `*_on_union_instruments.csv` |
| 2b | `mvmr_stage2b_harmonise.py` | the two PGC sumstats plus stage-1/1b tables | `mvmr_matrix_harmonised.csv` |
| 3 | `mvmr_stage3_run.R` | `mvmr_matrix_harmonised.csv` | `mvmr_old_vs_new_design.csv` (the 264-variant union design, superseded) |
| 4 | `mvmr_stage4_full_models.R` | `mvmr_matrix_harmonised.csv` | `mvmr_full_models.csv` (same, superseded) |
| 5 | `mvmr_stage5_exposure_only.R` | `mvmr_matrix_harmonised.csv` | `mvmr_exposure_only.csv` (package fit and conditional F, four models) |
| 6 | `mvmr_stage6_matched.py` | stage-5 output plus the matrix | `mvmr_matched_comparison.csv` (matched univariable, both estimators, coverage and strand diagnostics) |
| 7 | `mvmr_stage7_full_v7.py` | matrix plus stage-5 output | `mvmr_full_models_v7.csv` (all sixteen coefficients, both estimators) |

Stages 3 and 4 are kept because they document the design that was withdrawn; the manuscript does not report them.
Stage 2 (`mvmr_stage2_matrix.py`) is not here: it aligned the traits inconsistently and retained zero independent
covariate instruments, and 2b supersedes it.

## Figure 3

`make_figure3.py` regenerates the submitted figure from `mvmr_matched_comparison.csv`. It supersedes
`figures/create_figure3.py`, which still draws the withdrawn estimates and must not be used. Run from the
repository root:

    python -X utf8 08_revision_analyses/mvmr/make_figure3.py

## Running it

Stages 3-7 run from the repository alone:

    Rscript 08_revision_analyses/mvmr/mvmr_stage5_exposure_only.R
    python -X utf8 08_revision_analyses/mvmr/mvmr_stage6_matched.py
    python -X utf8 08_revision_analyses/mvmr/mvmr_stage7_full_v7.py

Stages 1, 1b and 2 need the summary statistics, which are not redistributed. Point `MR_SUMSTATS_DIR` at a
directory holding `影响因素/ieu-a-2.vcf.gz`, `ieu-a-80.vcf.gz`, `ieu-b-4877.vcf.gz` and
`数据/01_formatted/MDD_PGC_2025_EUR_no23andMe_formatted.tsv.gz`, `Anxiety_PGC_2026_EUR_formatted.tsv.gz`,
then run from the repository root. Covariate instruments were clumped by position at 1 Mb rather than by LD,
because no LD-clumping tool was available on the analysis machine; that is a stated limitation in the manuscript.
## The second correction: pairing a corrected estimate with an uncorrected reference

Fixing the instrument sets exposed a second mismatch. v7 of the manuscript put a measurement-error-CORRECTED
direct estimate next to a PLAIN random-effects IVW reference and reported ratios of 1.08 and 1.10, which read as
conditioning raising the estimate. The correction alone inflates a beta by 2.4-13.8% depending on instrument
strength, so part of that was the estimator. With both sides under the same estimator:

| model | n | both plain | both corrected | v7 (mixed) |
|---|---|---|---|---|
| forward, exposure-only | 237 | 1.00 | 1.02 | 1.08 |
| forward, exposure-instrumented | 84 | 0.93 | 0.84 | 0.86 |
| reverse, exposure-only | 180 | 0.96 | 0.96 | 1.10 |
| reverse, exposure-instrumented | 27 | 0.97 | 0.97 | 1.00 |

The manuscript now reports the corrected column and gives the plain one as a check. The honest statement is that
conditioning moves the forward estimate by 2% and the reverse by 4%, not that it raises them.

## Checks added in the same round

`mvmr_paired_v8.csv` also carries the joint Wald test that the three covariate coefficients are zero (chi-square
on 3 df: 8.78, P = 0.032 forward and 36.5, P = 6 x 10^-8 reverse in the exposure-only models), and the estimate
for the 114 MDD instruments dropped by the covariate-coverage filter (univariable OR 2.27, 95% CI 2.13-2.42,
against 2.44, 2.26-2.64 for the 84 retained). `mvmr_ld_sensitivity.csv` refits after removing the covariate
instruments that position clumping leaves near a primary instrument - 14 within 1 Mb, 4 within 100 kb, closest
pair 1.06 kb - and finds little: the forward direct estimate moves from 2.42 to 2.41 and the BMI conditional F
from 16.6 to 15.9. The difference test of Sanderson and Windmeijer between the univariable and multivariable
estimate of one exposure is not implemented in MVMR 0.4.8 and was not implemented here; the joint Wald test is
reported instead and the manuscript says so.

## Running it

    Rscript 08_revision_analyses/mvmr/mvmr_stage8b_ld.R
    python -X utf8 08_revision_analyses/mvmr/mvmr_stage8_paired_v8.py
    python -X utf8 08_revision_analyses/mvmr/make_figure3.py
