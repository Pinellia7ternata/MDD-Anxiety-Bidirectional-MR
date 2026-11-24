# Revision analyses

Scripts and tabulated outputs produced while responding to peer review. Run from the
repository root; the `mvmr/` subdirectory holds the multivariable rerun and its own README, which explains why the previous Table 3 was withdrawn rather than revised; inputs resolve through `_locate.py` (set `MR_DATA_DIR` to override, or
`MR_OUT_DIR` for where results are written). The `data/` directory holds the harmonised
instrument tables and the locus-level colocalization output these scripts read.

| script | produces | manuscript location |
|---|---|---|
| `recompute_corrected_estimates.py` | random-effects IVW, Cochran Q, variance explained, directionality, coloc-filtered estimates | Tables 2 and 4 |
| `null_controls_locus_specific.py` | uniform-random, proportional-decile and Q-greedy deletion comparators | Table 5 (round-1 values) |
| `window_and_class_sensitivity.py` | window widths and coloc hypothesis-class decomposition | Sections 3.4, Tables 4; S11, S12 |
| `coloc_power_simulation.py` | single-SNP posterior recomputation and the complete-sharing null | Section 3.4; S13 |
| `power_and_comparator_analyses.py` | power experiment, retained-composition matched comparators, winner's-curse terciles | Sections 3.5-3.6, Table 5; S13, S14 |
| `mvmr/mvmr_stage5_exposure_only.R`, `mvmr/mvmr_stage6_matched.py`, `mvmr/mvmr_stage7_full_v7.py` | multivariable IVW on exposure-only instrument sets, each direct estimate paired with the univariable estimate on the same variants | Table 3, Figure 3; S7 |

The comparators in `power_and_comparator_analyses.py` supersede the proportional-decile
column shipped in the first revision pass: they resample within the F-decile profile of
the RETAINED instrument set, because the colocalizing instruments are systematically
stronger (mean F 49.8 deleted against 39.1 retained). `table5_composition_matched.csv`
is the file the manuscript's Table 5 is built from.

Independent-implementation check: `recompute_corrected_estimates.py` reproduces the
published random-effects IVW exactly (beta 0.8516, SE 0.0248, Q 363.2 forward; 0.5170,
0.0238, 151.3 reverse) and the reported F-statistics, which is why these scripts are
trusted for the filtering and comparator analyses and why Table 2 remains as produced by
TwoSampleMR.
