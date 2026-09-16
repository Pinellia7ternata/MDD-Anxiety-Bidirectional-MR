#!/usr/bin/env Rscript
# ============================================================
# MDD-Anxiety MVMR pipeline: local PGC MDD/Anxiety + OpenGWAS covariates
# Author: ChatGPT-assisted reproducible template
# Purpose: Complete MVMR adjustment for BMI, education, and smoking initiation.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(TwoSampleMR)
  library(ieugwasr)
  library(MVMR)
})

# ---------------- USER CONFIG ----------------
BASE <- "D:/2026年/文章/因果推断"
OUTDIR <- file.path(BASE, "06_MVMR")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# Local formatted GWAS files. Required columns can use common aliases.
MDD_FILE <- file.path(BASE, "01_formatted/MDD_PGC_2025_EUR_no23andMe_formatted.tsv.gz")
ANX_FILE <- file.path(BASE, "01_formatted/Anxiety_PGC_2026_EUR_formatted.tsv.gz")

# OpenGWAS IDs for covariates. Replace these if you have better/current local datasets.
ID_BMI <- "ieu-a-2"          # Body mass index
ID_EDU <- "ieu-a-1239"       # Years of schooling; fallback: ieu-a-755
ID_SMK <- "ieu-b-4877"       # Smoking initiation
COV_IDS <- c(BMI = ID_BMI, Education = ID_EDU, Smoking = ID_SMK)

P_THRESHOLD <- 5e-8
CLUMP_R2 <- 0.001
CLUMP_KB <- 10000
HARMONISE_ACTION <- 2

# Optional: local PLINK LD clumping. Preferred for reproducibility.
# If left blank, the script uses OpenGWAS API clumping, which may require a JWT token.
PLINK_BIN <- Sys.getenv("PLINK_BIN", unset = "")
BFILE_EUR <- Sys.getenv("BFILE_EUR", unset = "")  # e.g. D:/ref/1000G_EUR/1000G.EUR.QC

# ---------------- HELPERS ----------------
find_col <- function(dt, candidates, required = TRUE) {
  nms <- names(dt)
  low <- tolower(nms)
  cand <- tolower(candidates)
  idx <- match(cand, low)
  idx <- idx[!is.na(idx)]
  if (length(idx) > 0) return(nms[idx[1]])
  if (required) stop("Missing required column. Tried: ", paste(candidates, collapse = ", "))
  return(NA_character_)
}

read_local_gwas <- function(file, trait_name) {
  message("Reading ", trait_name, ": ", file)
  dt <- fread(file)
  col_snp <- find_col(dt, c("SNP", "rsid", "variant_id", "ID"))
  col_beta <- find_col(dt, c("beta", "BETA", "beta.exposure", "beta.outcome"))
  col_se <- find_col(dt, c("se", "SE", "standard_error", "se.exposure", "se.outcome"))
  col_p <- find_col(dt, c("p", "pval", "P", "PVAL", "p_value", "pval.exposure", "pval.outcome"))
  col_ea <- find_col(dt, c("effect_allele", "ea", "EA", "A1", "ALT", "effect_allele.exposure"))
  col_oa <- find_col(dt, c("other_allele", "oa", "OA", "A2", "REF", "other_allele.exposure"))
  col_eaf <- find_col(dt, c("eaf", "EAF", "af", "freq", "effect_allele_frequency"), required = FALSE)
  col_n <- find_col(dt, c("samplesize", "n", "N", "sample_size"), required = FALSE)
  out <- data.table(
    SNP = as.character(dt[[col_snp]]),
    beta = as.numeric(dt[[col_beta]]),
    se = as.numeric(dt[[col_se]]),
    pval = as.numeric(dt[[col_p]]),
    effect_allele = toupper(as.character(dt[[col_ea]])),
    other_allele = toupper(as.character(dt[[col_oa]]))
  )
  out$eaf <- if (!is.na(col_eaf)) as.numeric(dt[[col_eaf]]) else NA_real_
  out$samplesize <- if (!is.na(col_n)) as.numeric(dt[[col_n]]) else NA_real_
  out <- out[!is.na(SNP) & !is.na(beta) & !is.na(se) & !is.na(pval)]
  out <- unique(out, by = "SNP")
  out[, trait := trait_name]
  out
}

as_exposure_dat <- function(dt, trait_name) {
  data.frame(
    SNP = dt$SNP,
    exposure = trait_name,
    id.exposure = trait_name,
    beta.exposure = dt$beta,
    se.exposure = dt$se,
    pval.exposure = dt$pval,
    effect_allele.exposure = dt$effect_allele,
    other_allele.exposure = dt$other_allele,
    eaf.exposure = dt$eaf,
    samplesize.exposure = dt$samplesize,
    stringsAsFactors = FALSE
  )
}

as_outcome_dat <- function(dt, trait_name) {
  data.frame(
    SNP = dt$SNP,
    outcome = trait_name,
    id.outcome = trait_name,
    beta.outcome = dt$beta,
    se.outcome = dt$se,
    pval.outcome = dt$pval,
    effect_allele.outcome = dt$effect_allele,
    other_allele.outcome = dt$other_allele,
    eaf.outcome = dt$eaf,
    samplesize.outcome = dt$samplesize,
    stringsAsFactors = FALSE
  )
}

extract_local_for_snps <- function(gwas_dt, snps, trait_name, role = c("exposure", "outcome")) {
  role <- match.arg(role)
  sub <- gwas_dt[SNP %in% snps]
  if (role == "exposure") as_exposure_dat(sub, trait_name) else as_outcome_dat(sub, trait_name)
}

api_extract_instruments_safe <- function(id, trait_label) {
  message("Extracting instruments from OpenGWAS: ", trait_label, " / ", id)
  x <- extract_instruments(outcomes = id, p1 = P_THRESHOLD, clump = FALSE)
  if (nrow(x) == 0) stop("No instruments returned for ", trait_label, " (", id, ")")
  data.table(SNP = x$SNP, pval = x$pval.exposure, trait = trait_label)
}

clump_union_snps <- function(inst_dt) {
  inst_min <- inst_dt[, .(pval = min(pval, na.rm = TRUE)), by = SNP]
  clump_df <- data.frame(SNP = inst_min$SNP, pval = inst_min$pval, id = "MVMR_union")
  if (nzchar(PLINK_BIN) && nzchar(BFILE_EUR) && file.exists(PLINK_BIN)) {
    message("Local PLINK clumping requested, but TwoSampleMR::clump_data does not directly accept PLINK path in all versions.")
    message("Using TwoSampleMR::clump_data / ieugwasr backend unless you implement system PLINK clumping separately.")
  }
  clumped <- clump_data(clump_df, clump_kb = CLUMP_KB, clump_r2 = CLUMP_R2, pop = "EUR")
  unique(clumped$SNP)
}

harmonise_one_exposure_to_outcome <- function(exp_dat, out_dat, exposure_name, outcome_name) {
  h <- harmonise_data(exp_dat, out_dat, action = HARMONISE_ACTION)
  h <- h[h$mr_keep == TRUE, ]
  data.table(
    SNP = h$SNP,
    beta_exp = h$beta.exposure,
    se_exp = h$se.exposure,
    pval_exp = h$pval.exposure,
    beta_out = h$beta.outcome,
    se_out = h$se.outcome,
    exposure = exposure_name,
    outcome = outcome_name
  )
}

run_mvmr_direction <- function(main_exposure = c("MDD", "Anxiety"), outcome = c("Anxiety", "MDD"), mdd_gwas, anx_gwas) {
  main_exposure <- match.arg(main_exposure)
  outcome <- match.arg(outcome)
  message("\n==============================")
  message("Running MVMR: ", main_exposure, " + BMI + Education + Smoking -> ", outcome)
  message("==============================")

  # instruments: union of main exposure and covariate genome-wide significant SNPs
  main_gwas <- if (main_exposure == "MDD") mdd_gwas else anx_gwas
  outcome_gwas <- if (outcome == "MDD") mdd_gwas else anx_gwas
  main_inst <- main_gwas[pval < P_THRESHOLD, .(SNP, pval, trait = main_exposure)]
  cov_inst <- rbindlist(lapply(names(COV_IDS), function(nm) api_extract_instruments_safe(COV_IDS[[nm]], nm)), fill = TRUE)
  union_inst <- rbind(main_inst, cov_inst, fill = TRUE)
  fwrite(union_inst, file.path(OUTDIR, paste0(main_exposure, "_to_", outcome, "_union_instruments_preclump.csv")))

  snps <- clump_union_snps(union_inst)
  message("Clumped union SNPs: ", length(snps))
  fwrite(data.table(SNP = snps), file.path(OUTDIR, paste0(main_exposure, "_to_", outcome, "_union_instruments_clumped.csv")))

  # outcome local data
  out_dat <- extract_local_for_snps(outcome_gwas, snps, outcome, role = "outcome")

  # main exposure local data
  exp_main <- extract_local_for_snps(main_gwas, snps, main_exposure, role = "exposure")
  h_main <- harmonise_one_exposure_to_outcome(exp_main, out_dat, main_exposure, outcome)
  setnames(h_main, c("beta_exp", "se_exp", "pval_exp"), c(paste0("beta_", main_exposure), paste0("se_", main_exposure), paste0("pval_", main_exposure)))
  base <- h_main[, .(SNP, beta_out, se_out, get(paste0("beta_", main_exposure)), get(paste0("se_", main_exposure)), get(paste0("pval_", main_exposure)))]
  setnames(base, c("V4", "V5", "V6"), c(paste0("beta_", main_exposure), paste0("se_", main_exposure), paste0("pval_", main_exposure)))

  # covariates from OpenGWAS as exposures, harmonised to same local outcome
  for (covnm in names(COV_IDS)) {
    message("Extracting covariate associations for ", covnm)
    od <- extract_outcome_data(snps = snps, outcomes = COV_IDS[[covnm]], proxies = FALSE)
    exp_cov <- data.frame(
      SNP = od$SNP,
      exposure = covnm,
      id.exposure = covnm,
      beta.exposure = od$beta.outcome,
      se.exposure = od$se.outcome,
      pval.exposure = od$pval.outcome,
      effect_allele.exposure = od$effect_allele.outcome,
      other_allele.exposure = od$other_allele.outcome,
      eaf.exposure = od$eaf.outcome,
      samplesize.exposure = od$samplesize.outcome,
      stringsAsFactors = FALSE
    )
    h_cov <- harmonise_one_exposure_to_outcome(exp_cov, out_dat, covnm, outcome)
    h_cov <- h_cov[, .(SNP, beta_exp, se_exp, pval_exp)]
    setnames(h_cov, c("beta_exp", "se_exp", "pval_exp"), c(paste0("beta_", covnm), paste0("se_", covnm), paste0("pval_", covnm)))
    base <- merge(base, h_cov, by = "SNP", all = FALSE)
  }

  fwrite(base, file.path(OUTDIR, paste0(main_exposure, "_to_", outcome, "_MVMR_harmonised_matrix.csv")))

  exposure_names <- c(main_exposure, "BMI", "Education", "Smoking")
  BX <- as.matrix(base[, paste0("beta_", exposure_names), with = FALSE])
  seBX <- as.matrix(base[, paste0("se_", exposure_names), with = FALSE])
  BY <- base$beta_out
  seBY <- base$se_out

  r_input <- format_mvmr(BXGs = BX, BYG = BY, seBXGs = seBX, seBYG = seBY, RSID = base$SNP)

  ivw <- ivw_mvmr(r_input, model = "random")
  strength <- tryCatch(strength_mvmr(r_input, gencov = 0), error = function(e) e)
  pleio <- tryCatch(pleiotropy_mvmr(r_input, gencov = 0), error = function(e) e)

  sink(file.path(OUTDIR, paste0(main_exposure, "_to_", outcome, "_MVMR_raw_output.txt")))
  cat("MVMR direction:", main_exposure, "+ covariates ->", outcome, "\n\n")
  cat("Number of SNPs:", nrow(base), "\n\n")
  cat("--- IVW MVMR random effects ---\n")
  print(ivw)
  cat("\n--- Conditional strength ---\n")
  print(strength)
  cat("\n--- Pleiotropy / heterogeneity ---\n")
  print(pleio)
  sink()

  # Try to coerce core estimates into a CSV. MVMR package output format may vary across versions.
  est <- as.data.frame(ivw)
  fwrite(est, file.path(OUTDIR, paste0(main_exposure, "_to_", outcome, "_MVMR_IVW_results.csv")))
  invisible(list(harmonised = base, ivw = ivw, strength = strength, pleio = pleio))
}

# ---------------- RUN ----------------
setwd(BASE)
mdd <- read_local_gwas(MDD_FILE, "MDD")
anx <- read_local_gwas(ANX_FILE, "Anxiety")

# Main MVMR direction and reverse direction.
res1 <- run_mvmr_direction("MDD", "Anxiety", mdd, anx)
res2 <- run_mvmr_direction("Anxiety", "MDD", mdd, anx)

message("\nDone. Outputs written to: ", OUTDIR)
