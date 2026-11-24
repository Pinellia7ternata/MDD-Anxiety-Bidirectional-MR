# Bidirectional Mendelian Randomization of Major Depressive Disorder and Anxiety Disorders

R/Python analysis codebase accompanying the manuscript
*"Shared genetic architecture and the limits of directional inference in major depressive disorder and anxiety disorders comorbidity: a bidirectional Mendelian randomization study"*
(currently under peer review).

This repository provides the full analysis pipeline: instrument
selection, harmonisation, bidirectional univariable MR, multivariable MR (MVMR),
LDSC genetic correlation, colocalization, functional annotation, and the
sensitivity / robustness and revision-stage analyses reported in the paper.

## Repository structure

| Folder | Contents |
|--------|----------|
| `01_instrument_selection/` | LD-clumping (real-LD via IEU OpenGWAS), strict MR instrument-selection pipeline (clumping → harmonisation → export). |
| `02_bidirectional_univariable_MR/` | Bidirectional two-sample MR (MDD → Anxiety and Anxiety → MDD) under real LD. |
| `03_MVMR/` | Multivariable MR (original design with BMI, education, smoking) and MVMR supplementary analyses. Superseded by the refits in `08_revision_analyses/mvmr/`; see its README. |
| `04_LDSC/` | Comprehensive supplementary analyses script (includes LDSC genetic correlation, `rg`). |
| `05_colocalization/` | Colocalization (`coloc.abf` locus-level + strict) and GTEx eQTL colocalization. |
| `06_functional_annotation/` | MAGMA gene-set analysis, FUMA input preparation, and SuSiE fine-mapping. |
| `07_sensitivity_robustness/` | MR-PRESSO scripts (residual heterogeneity checks). |
| `08_revision_analyses/` | Revision-stage analyses reported in the current manuscript: corrected variance-explained and Cochran Q statistics, the colocalization filter, deletion comparators (random / composition-matched / Q-greedy), the colocalization power experiment, window-and-class sensitivity, and winner's-curse stratification. `08_revision_analyses/mvmr/` holds the multivariable refits behind Table 3 and Supplementary Tables S7 and S15 (harmonised instrument matrix, covariate coverage diagnostics, paired-estimator comparison, covariate joint tests, LD-overlap sensitivity). |
| `SUPERSEDED_FAILED/` | Scripts for analyses that were attempted but are not reported in the manuscript (a Radial MR / "MR-RAPS" prototype with known defects, CAUSE, LCV, UKB external validation), with the reason for each recorded in `SUPERSEDED_FAILED/README.md`. Kept for audit transparency. |
| `figures/` | Python scripts that regenerate the manuscript figures (Forest plots, MVMR, supplementary figures). Note: the Figure 3 source of record is `08_revision_analyses/mvmr/make_figure3.py`. |
| `data/` | Harmonised instrument lists used in the primary analyses (`S1`: MDD→Anxiety, 198 instruments; `S2`: Anxiety→MDD, 50 instruments), the strict locus-level colocalization table (`coloc_locus_level_strict.csv`), and revision-stage analysis records (`S10`–`S12`). |

## Dependencies

**R** (tested on R 4.5.x): `TwoSampleMR`, `MRInstruments`, `ieugwasr`,
`MRPRESSO`, `MVMR`, `RadialMR`, `coloc`, `susieR`, `CAUSE`, `LCV`,
`data.table`, `dplyr`.

**Python** (3.12+): `pandas`, `numpy`, `scipy`, `matplotlib`, plus the
[`ldsc`](https://github.com/bulik/ldsc) package for the LDSC step.

An IEU OpenGWAS API token is required for the instrument-selection and MR
scripts (`ieugwasr::get_access_token()`); set it before running
`01_instrument_selection/`.

## Recommended run order

1. `01_instrument_selection/` – select and clump instruments (real LD), harmonise, export.
2. `02_bidirectional_univariable_MR/` – primary bidirectional MR.
3. `04_LDSC/` – genetic correlation.
4. `05_colocalization/` – colocalization with eQTL/pQTL resources.
5. `06_functional_annotation/` – MAGMA / FUMA / SuSiE.
6. `07_sensitivity_robustness/` – MR-PRESSO.
7. `08_revision_analyses/` – revision-stage analyses (colocalization filter, deletion comparators, power experiment, winner's-curse stratification); then `08_revision_analyses/mvmr/` for the multivariable refits. These reproduce from the repository alone; see `08_revision_analyses/README.md` and `08_revision_analyses/mvmr/README.md` for the verified run order.
8. `figures/` – regenerate manuscript figures from the outputs above.

## Reproducibility note

The multivariable refits in `08_revision_analyses/mvmr/` reproduce from the
repository alone and have been verified end-to-end from a fresh checkout. The
older pipeline scripts (steps 1–6) reference upstream GWAS summary statistics
that are not redistributed, so those steps document the analysis rather than
reproducing it unattended. Two corrections made during revision are documented
in `08_revision_analyses/mvmr/README.md`: the withdrawal of the original
Table 3 (which could not be reproduced) and the replacement of a mixed
corrected/uncorrected estimator pairing with a same-estimator comparison.

## Data availability

Exposure and outcome GWAS summary statistics are drawn from public resources
(Psychiatric Genomics Consortium, IEU OpenGWAS, and the originating consortia).
The harmonised instrument lists underpinning the primary results are provided in
`data/`. Full summary statistics are available from the originating consortia
under their respective access terms.

## License

Code is released under the MIT License (see `LICENSE`). Manuscript text and
figures are © the authors.

## Citation

If you use this code, please cite the accompanying manuscript; citation details
will be added here on publication.
