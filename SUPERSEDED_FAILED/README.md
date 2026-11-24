# SUPERSEDED / FAILED analyses

Moved here on 2026-11-02 while responding to peer review. Nothing in
this directory is reported in the manuscript. Each file is kept with the reason it was set
aside, so the audit trail stays intact.

| file | why it is here |
|---|---|
| `08_Radial_MR_RAPS.R` | Named "MR-RAPS" and "Radial MR" but implements neither. se = sqrt(1/sum(beta_exp^2/se_exp^2)) = 1/sqrt(sum F) uses the exposure standard errors as IVW weights, so every confidence interval it emits is roughly 40% too narrow; the outlier loop stopped at max_iter rather than on convergence, and residual Q per instrument rose after the "outliers" were removed. The manuscript no longer reports either method. Superseded by 08_revision_analyses/. |
| `06_cause_analysis.R` | CAUSE was never run to completion (the companion LCV script header records that CAUSE needed Rtools and was abandoned). Not reported in the manuscript. |
| `10_run_cause_bayesian.R` | Bayesian CAUSE modelling; no model posterior output exists on disk and it is not reported. Retained here because a CAUSE or LCV decomposition is the natural next analysis for the question this paper asks, and it should be run before this script is reused. |
| `10_run_lcv_causal_inference.R` | Latent causal variable analysis written as a CAUSE substitute; no saved results, not reported in the manuscript. |
| `10_ukb_validation.R` | FAILED. ukb_match_statistics.csv records N_Matched = 0 and Match_Rate = 0 for both UKB phenotypes (20002_1287, 20002_1286). The downstream estimates are therefore coordinate-matched at best and carry a scale error: "UKB_Depression -> PGC_MDD" gives beta = 15.18, OR = 3.9e6. |
| `10_ukb_external_validation.R` | FAILED for the same reason as 10_ukb_validation.R; its companion ukb_comparison_with_primary.csv marks a primary OR of 2.343 against an external OR of 1.005 as Consistent = TRUE. |
| `10_external_validation.R` | External-validation driver script; results not reported and not reproducible from the repository contents. |
| `15_causal_inference_summary.R` | Emits CAUSE_style_summary.csv, whose "conditional F = 173.3 / 144.2" and "genetic correlation r = 0.263" contradict the manuscript values (15.6 / 12.9 and rg = 0.90). Stale summary from a superseded model specification; must not be used. |

## Known limitation of the remaining pipeline

Scripts still outside this folder hard-code the author's local layout
(`BASE <- "D:/2026年/文章/因果推断"`) and read upstream GWAS summary statistics that are not
redistributed here. They document the analysis but will not run as-is on another machine;
paths must be parameterised before the repository can be described as a runnable pipeline.
