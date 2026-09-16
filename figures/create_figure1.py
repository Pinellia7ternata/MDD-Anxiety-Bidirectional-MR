#!/usr/bin/env python3
"""
Figure 1: Study design and analysis workflow
MAXIMUM SPACING - Clean SVG with editable text
"""
import matplotlib
matplotlib.use('SVG')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import numpy as np

# Configure for editable text
plt.rcParams['svg.fonttype'] = 'none'  # Keep text as text elements

# Enlarged figure with MAXIMUM spacing
fig, ax = plt.subplots(figsize=(26, 36))
ax.set_xlim(0, 26)
ax.set_ylim(0, 36)
ax.axis('off')
ax.set_facecolor('white')
fig.patch.set_facecolor('white')

# Color scheme
BLUE_DARK = '#1A5276'
BLUE_MED = '#2471A3'
PURPLE = '#76448A'
RED = '#A93226'
ORANGE = '#CA6F1E'
GREEN = '#1E8449'
GRAY = '#717D7E'

def draw_box(ax, x, y, w, h, text, color, fontsize=14, bold=True, text_color='white'):
    """Draw a rounded box with text - MAXIMUM spacing"""
    box = FancyBboxPatch((x - w/2, y - h/2), w, h,
                          boxstyle="round,pad=0.15,rounding_size=0.4",
                          facecolor=color, edgecolor='white', linewidth=4)
    ax.add_patch(box)
    weight = 'bold' if bold else 'normal'
    ax.text(x, y, text, ha='center', va='center', fontsize=fontsize,
            color=text_color, fontweight=weight, wrap=True,
            multialignment='center', linespacing=1.8)

def draw_arrow(ax, x1, y1, x2, y2, color='#2C3E50'):
    """Draw an arrow with more space"""
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='->', color=color, lw=4))

# Title - LARGER
ax.text(13, 34.5, 'Figure 1', fontsize=24, fontweight='bold', color='#1A5276', ha='center')
ax.text(13, 33.8, 'Study design and workflow of the bidirectional MR, MVMR, LDSC,\ncolocalization, and functional annotation analyses',
        fontsize=16, color='#2C3E50', ha='center', style='italic')

# Level 1: GWAS Data Sources - MORE SPACE ABOVE
draw_box(ax, 13, 31.0, 22, 2.5, 
         'GWAS Data Sources\n\nPGC 2025 MDD: 688,808 cases + 4,364,225 controls  |  PGC 2026 Anxiety: 122,341 cases\nCovariates: BMI (ieu-a-2), Education (ieu-a-80), Smoking initiation (ieu-b-4877)',
         BLUE_DARK, fontsize=15)

# Level 2: Quality Control - MORE SPACE
draw_box(ax, 13, 27.5, 16, 1.8, 'Quality Control and Harmonization', BLUE_MED, fontsize=16)
draw_arrow(ax, 13, 29.8, 13, 28.6)

# Level 3: LD Clumping - MORE SPACE
draw_box(ax, 13, 24.2, 18, 1.8, 
         'LD Clumping against the 1000 Genomes European Reference Panel\nP < 5 x 10^-8, r2 < 0.001, 10,000 kb window', 
         PURPLE, fontsize=15)
draw_arrow(ax, 13, 26.3, 13, 25.3)

# Level 4: Final Instruments - MUCH MORE SPACE
draw_box(ax, 6.5, 21.0, 10, 1.8, 'MDD Instruments\n198 SNPs', BLUE_MED, fontsize=16)
draw_box(ax, 19.5, 21.0, 10, 1.8, 'Anxiety Instruments\n57 SNPs\n(50 after harmonization)', BLUE_MED, fontsize=16)
draw_arrow(ax, 9, 23.4, 7, 22.2)
draw_arrow(ax, 17, 23.4, 18.5, 22.2)

# Level 5: Bidirectional MR Methods - MUCH MORE SPACE
ax.text(6.5, 18.5, 'MDD to Anxiety MR\n(198 SNPs)', ha='center', fontsize=16, 
        fontweight='bold', color=BLUE_DARK)
ax.text(19.5, 18.5, 'Anxiety to MDD MR\n(50 SNPs)', ha='center', fontsize=16, 
        fontweight='bold', color=BLUE_DARK)

draw_box(ax, 6.5, 15.5, 11, 3.2, 
         'IVW, MR-Egger,\nWeighted Median,\nWeighted Mode, Simple Mode,\nMR-PRESSO, Radial MR,\nHeterogeneity Test,\nEgger Intercept',
         RED, fontsize=14, bold=False)

draw_box(ax, 19.5, 15.5, 11, 3.2, 
         'IVW, MR-Egger,\nWeighted Median,\nWeighted Mode, Simple Mode,\nMR-PRESSO, Radial MR,\nHeterogeneity Test,\nEgger Intercept',
         RED, fontsize=14, bold=False)

draw_arrow(ax, 6.5, 17.3, 6.5, 16.0)
draw_arrow(ax, 19.5, 17.3, 19.5, 16.0)

# Level 6: MVMR - MORE SPACE
draw_box(ax, 13, 11.5, 20, 2.2, 
         'Multivariable MR (MVMR)\nAdjusting for BMI, Education, and Smoking Initiation\n\nMDD to Anxiety: conditional F = 15.6  |  Anxiety to MDD: conditional F = 12.9',
         ORANGE, fontsize=15)

draw_arrow(ax, 6.5, 13.7, 6.5, 12.5, color=ORANGE)
draw_arrow(ax, 19.5, 13.7, 19.5, 12.5, color=ORANGE)

# Level 7: Shared Genetic Architecture - MORE SPACE
draw_box(ax, 13, 8.2, 18, 1.8, 
         'Shared Genetic Architecture Analyses\nLD Score Regression (LDSC)  |  Colocalization (coloc)',
         GREEN, fontsize=15)
draw_arrow(ax, 13, 10.3, 13, 9.1)

# Level 8: Exploratory Functional Annotation - MORE SPACE
draw_box(ax, 13, 5.2, 20, 1.8, 
         'Exploratory Functional Annotation (Hypothesis-Generating)\nFUMA  |  MAGMA Gene Analysis  |  GSEA  |  Brain eQTL Enrichment',
         GRAY, fontsize=15)
draw_arrow(ax, 13, 7.2, 13, 6.1)

# Level 9: Results - MUCH MORE SPACE
ax.text(6.5, 2.8, 'Results:', ha='center', fontsize=18, fontweight='bold', color=BLUE_DARK)
ax.text(19.5, 2.8, 'Results:', ha='center', fontsize=18, fontweight='bold', color=BLUE_DARK)

draw_box(ax, 6.5, 0.5, 11, 4.0, 
         'MDD to Anxiety\n\nIVW OR = 2.34\n(95% CI: 2.23-2.46)\nP = 1.66 x 10^-257\n\nNo directional pleiotropy\n(Egger P = 0.864)',
         GREEN, fontsize=15)

draw_box(ax, 19.5, 0.5, 11, 4.0, 
         'Anxiety to MDD\n\nIVW OR = 1.68\n(95% CI: 1.60-1.76)\nP = 1.97 x 10^-104\n\nDirectional pleiotropy\n(Egger P = 0.006)',
         ORANGE, fontsize=15)

draw_arrow(ax, 6.5, 4.8, 6.5, 2.5, color=GREEN)
draw_arrow(ax, 19.5, 4.8, 19.5, 2.5, color=ORANGE)

plt.tight_layout()
plt.savefig(r'D:\2026年\文章\因果推断\图表\Figure_1_Study_Design_Workflow.svg', 
            format='svg', facecolor='white')
print("[SUCCESS] Figure 1 saved")
plt.close()
