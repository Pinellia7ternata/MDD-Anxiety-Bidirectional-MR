# Bidirectional Mendelian Randomization of Major Depressive Disorder and Anxiety Disorders

R/Python analysis codebase accompanying the manuscript
*"Bidirectional causal relationship between major depressive disorder and anxiety disorders: a two-sample and multivariable Mendelian randomization study"*
(submitted to the *Journal of Affective Disorders*).

This repository provides the full, reproducible analysis pipeline: instrument
selection, harmonisation, bidirectional univariable MR, multivariable MR (MVMR),
LDSC genetic correlation, colocalization, functional annotation, and the
sensitivity / robustness analyses reported in the paper.

## Repository structure

| Folder | Contents |
|--------|----------|
| `01_instrument_selection/` | LD-clumping (real-LD via IEU OpenGWAS), strict MR instrument-selection pipeline (clumping → harmonisation → export). |
| `02_bidirectional_univariable_MR/` | Bidirectional two-sample MR (MDD → Anxiety and Anxiety → MDD) under real LD. |
| `03_MVMR/` | Multivariable MR (full model + stratified/sensitivity models with BMI, education, smoking) and MVMR supplementary analyses. |
| `04_LDSC/` | Comprehensive supplementary analyses script (includes LDSC genetic correlation, `rg`). |
| `05_colocalization/` | Colocalization (`coloc.abf` locus-level + strict) and GTEx eQTL colocalization. |
| `06_functional_annotation/` | MAGMA gene-set analysis, FUMA input preparation, and SuSiE fine-mapping. |
| `07_sensitivity_robustness/` | Radial MR + MR-RAPS, MR-PRESSO, CAUSE, LCV, external/UKB validation, and the causal-inference summary orchestrator. |
| `figures/` | Python scripts that regenerate the manuscript figures (Forest plots, MVMR, funnel/Supt. figures). |
| `data/` | Harmonised instrument lists used in the primary analyses (`S1`: MDD→Anxiety, 198 instruments; `S2`: Anxiety→MDD, 50 instruments). |

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
3. `03_MVMR/` – multivariable MR (full + stratified).
4. `04_LDSC/` – genetic correlation.
5. `05_colocalization/` – colocalization with eQTL/pQTL resources.
6. `06_functional_annotation/` – MAGMA / FUMA / SuSiE.
7. `07_sensitivity_robustness/` – Radial MR, MR-PRESSO, CAUSE, LCV, validation.
8. `figures/` – regenerate manuscript figures from the outputs above.

## Data availability

Exposure and outcome GWAS summary statistics are drawn from public resources
(IEU OpenGWAS / EBI GWAS Catalog / UK Biobank and published MDD/anxiety GWAS).
The harmonised instrument lists underpinning the primary results are provided in
`data/`. Full summary statistics are available from the originating consortia
under their respective access terms.

## License

Code is released under the MIT License (see `LICENSE`). Manuscript text and
figures are © the authors.

## Citation

If you use this code, please cite the accompanying *Journal of Affective Disorders* manuscript.
