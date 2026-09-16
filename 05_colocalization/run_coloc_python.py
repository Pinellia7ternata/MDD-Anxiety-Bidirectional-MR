#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Colocalization Analysis using coloc.abf formula
Based on Giambartolomei et al. PLoS Genet 2014
Author: Chen Chao | 2026-05-22
"""

import numpy as np
import pandas as pd
from scipy import stats

print("=== COLOCALIZATION ANALYSIS (Python Implementation) ===\n")

# ==================================================
# Step 1: Load instrument data
# ==================================================
print("[Step 1] Loading instrument data...")

instrument_file = r"D:\2026年\文章\因果推断\GWAS_data\MVMR_instruments_with_covariates.csv"
data = pd.read_csv(instrument_file)

print(f"Loaded {len(data)} SNPs")

# ==================================================
# Step 2: Calculate coloc posterior probabilities
# ==================================================
print("\n[Step 2] Calculating coloc posterior probabilities...")

# Extract MDD and Anxiety data
beta1 = data['beta_MDD'].values  # MDD
se1 = data['se_MDD'].values
beta2 = data['beta_Anxiety'].values  # Anxiety
se2 = data['se_Anxiety'].values

# Calculate Z-scores and p-values
z1 = beta1 / se1
z2 = beta2 / se2
p1 = 2 * stats.norm.sf(np.abs(z1))
p2 = 2 * stats.norm.sf(np.abs(z2))

# Calculate approximate Bayes factors for each hypothesis
# Following coloc.abf formula from Giambartolomei et al.

# Prior probabilities (default from coloc)
p1_prior = 1e-4  # Prior on SNP being associated with trait 1
p2_prior = 1e-4  # Prior on SNP being associated with trait 2
p12_prior = 1e-5  # Prior on SNP being associated with both traits

# Calculate Bayes factors for each SNP
# Using approximate BF from Wakefield (2009)
def calculate_abf(beta, se, N=500000):
    """Calculate approximate Bayes factor"""
    # Assume prior variance = 0.01 (default in coloc)
    omega = 0.01
    # Variance of beta
    V = se ** 2
    # Bayes factor
    r = omega / (omega + V)
    z = beta / se
    abf = np.sqrt(1 - r) * np.exp(0.5 * r * z ** 2)
    return abf

abf1 = calculate_abf(beta1, se1)
abf2 = calculate_abf(beta2, se2)

# Calculate posterior probabilities for each hypothesis
# H0: no association
# H1: associated with trait 1 only
# H2: associated with trait 2 only
# H3: associated with both traits (different causal variants)
# H4: associated with both traits (shared causal variant)

# Sum over all SNPs
sum_abf1 = np.sum(abf1)
sum_abf2 = np.sum(abf2)
sum_abf1_abf2 = np.sum(abf1 * abf2)

# Posterior probabilities
pp_h0 = 1.0
pp_h1 = p1_prior * sum_abf1
pp_h2 = p2_prior * sum_abf2
pp_h3 = p1_prior * p2_prior * sum_abf1 * sum_abf2
pp_h4 = p12_prior * sum_abf1_abf2

# Normalize
total = pp_h0 + pp_h1 + pp_h2 + pp_h3 + pp_h4
pp_h0 /= total
pp_h1 /= total
pp_h2 /= total
pp_h3 /= total
pp_h4 /= total

# ==================================================
# Step 3: Print results
# ==================================================
print("\n[Step 3] Colocalization Results:\n")

print("Posterior Probabilities:")
print(f"  PP.H0 (no association):    {pp_h0:.4f}")
print(f"  PP.H1 (MDD only):          {pp_h1:.4f}")
print(f"  PP.H2 (Anxiety only):      {pp_h2:.4f}")
print(f"  PP.H3 (two distinct):      {pp_h3:.4f}")
print(f"  PP.H4 (shared causal):     {pp_h4:.4f}")

print("\n=== INTERPRETATION ===")
if pp_h4 > 0.5:
    print("[STRONG] PP.H4 > 0.5: Strong evidence for shared causal variant")
    print("   MDD and Anxiety share a common causal SNP at this locus.")
elif pp_h4 > 0.2:
    print("[MODERATE] PP.H4 > 0.2: Moderate evidence for shared causal variant")
    print("   Results suggest possible shared causality, but not conclusive.")
elif pp_h3 > 0.5:
    print("[DISTINCT] PP.H3 > 0.5: Evidence for two distinct causal variants")
    print("   MDD and Anxiety have different causal SNPs at this locus.")
else:
    print("[INCONCLUSIVE] No strong evidence for colocalization.")
    print("   Results are inconclusive or suggest no shared causality.")

# ==================================================
# Step 4: Save results
# ==================================================
print("\n[Step 4] Saving results...")

results = pd.DataFrame({
    'Hypothesis': ['PP.H0', 'PP.H1', 'PP.H2', 'PP.H3', 'PP.H4'],
    'Description': [
        'No association',
        'MDD only',
        'Anxiety only',
        'Two distinct causal variants',
        'Shared causal variant'
    ],
    'Posterior_Probability': [pp_h0, pp_h1, pp_h2, pp_h3, pp_h4]
})

output_file = r"D:\2026年\文章\因果推断\MR_results\coloc_python_results.csv"
results.to_csv(output_file, index=False)
print(f"Results saved to: {output_file}")

# Save summary
summary_file = r"D:\2026年\文章\因果推断\MR_results\coloc_python_summary.txt"
with open(summary_file, 'w', encoding='utf-8') as f:
    f.write("=== COLOCALIZATION RESULTS (Python Implementation) ===\n\n")
    f.write("Data source:\n")
    f.write(f"  MDD: PGC 2025 (5,419 SNPs)\n")
    f.write(f"  Anxiety: PGC 2026 (5,419 SNPs)\n\n")
    f.write("Posterior Probabilities:\n")
    f.write(f"  PP.H0 (no association):    {pp_h0:.4f}\n")
    f.write(f"  PP.H1 (MDD only):          {pp_h1:.4f}\n")
    f.write(f"  PP.H2 (Anxiety only):      {pp_h2:.4f}\n")
    f.write(f"  PP.H3 (two distinct):      {pp_h3:.4f}\n")
    f.write(f"  PP.H4 (shared causal):     {pp_h4:.4f}\n\n")
    f.write("=== INTERPRETATION ===\n")
    if pp_h4 > 0.5:
        f.write("✅ PP.H4 > 0.5: STRONG evidence for shared causal variant\n")
        f.write("   MDD and Anxiety share a common causal SNP at this locus.\n")
    elif pp_h4 > 0.2:
        f.write("⚠️  PP.H4 > 0.2: MODERATE evidence for shared causal variant\n")
        f.write("   Results suggest possible shared causality, but not conclusive.\n")
    else:
        f.write("⚠️  No strong evidence for colocalization.\n")

print(f"Summary saved to: {summary_file}")

print("\n=== COMPLETE ===\n")
