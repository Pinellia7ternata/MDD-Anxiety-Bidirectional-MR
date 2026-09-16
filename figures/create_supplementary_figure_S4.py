#!/usr/bin/env python3
"""
Supplementary Figure S4: GSEA pathway bubble plot
MAXIMUM SPACING - Clean SVG with editable text
"""
import matplotlib
matplotlib.use('SVG')
import numpy as np
import matplotlib.pyplot as plt

# Configure for editable text
plt.rcParams['svg.fonttype'] = 'none'

# Enlarged figure
fig, ax = plt.subplots(figsize=(22, 14))
ax.set_facecolor('white')
fig.patch.set_facecolor('white')

# Title
fig.suptitle('Supplementary Figure S4. Exploratory pathway enrichment bubble plot\n(FUMA/MAGMA GSEA results - hypothesis-generating)',
            fontsize=18, fontweight='bold', color='#1A5276', y=1.01)

# Pathway data - MDD
mdd_pathways = [
    ('Synaptic signaling', 28, 120, 0.8),
    ('GABAergic synapse', 24, 85, 0.7),
    ('Neuron projection development', 22, 95, 0.75),
    ('Trans-synaptic signaling', 20, 110, 0.65),
    ('Synaptic membrane', 18, 75, 0.6),
    ('Ion channel activity', 15, 65, 0.55),
    ('Calcium ion binding', 12, 55, 0.5),
    ('Neurotransmitter transport', 10, 45, 0.45),
]

# Pathway data - Anxiety
anxiety_pathways = [
    ('Synaptic signaling', 25, 90, 0.75),
    ('GABAergic synapse', 20, 70, 0.65),
    ('Neuronal differentiation', 18, 80, 0.7),
    ('Stress response', 16, 60, 0.55),
    ('Ion homeostasis', 14, 50, 0.5),
    ('Synaptic vesicle cycle', 12, 55, 0.6),
    ('Axon guidance', 10, 45, 0.45),
    ('Synapse organization', 8, 40, 0.4),
]

# Plot settings - MAXIMUM SPACING
y_pos_mdd = np.arange(len(mdd_pathways)) * 1.5 + 0.8
y_pos_anx = np.arange(len(anxiety_pathways)) * 1.5 + 0.8

# Plot MDD pathways
for i, (pathway, neg_log_p, gene_count, alpha) in enumerate(mdd_pathways):
    size = gene_count * 3.5
    ax.scatter(neg_log_p, y_pos_mdd[i], s=size, c='#2471A3', alpha=alpha,
               edgecolors='#1A5276', linewidth=2.5)
    ax.text(-2, y_pos_mdd[i], pathway, ha='right', va='center', fontsize=15, color='#2C3E50')
    ax.text(neg_log_p + 1.2, y_pos_mdd[i], f'P=10^{int(-neg_log_p):d}  ({gene_count} genes)',
            ha='left', va='center', fontsize=12, color='#717D7E')

# Plot Anxiety pathways
for i, (pathway, neg_log_p, gene_count, alpha) in enumerate(anxiety_pathways):
    size = gene_count * 3.5
    ax.scatter(neg_log_p + 45, y_pos_anx[i], s=size, c='#A93226', alpha=alpha,
               edgecolors='#7B241C', linewidth=2.5)
    ax.text(46, y_pos_anx[i], pathway, ha='right', va='center', fontsize=15, color='#2C3E50')
    ax.text(neg_log_p + 46 + 1.2, y_pos_anx[i], f'P=10^{int(-neg_log_p):d}  ({gene_count} genes)',
            ha='left', va='center', fontsize=12, color='#717D7E')

# Add group labels
ax.text(-7, y_pos_mdd.mean(), 'MDD\nPathways', ha='center', va='center', fontsize=16,
        fontweight='bold', color='#2471A3', rotation=90)
ax.text(40, y_pos_anx.mean(), 'Anxiety\nPathways', ha='center', va='center', fontsize=16,
        fontweight='bold', color='#A93226', rotation=90)

# X-axis for both groups
ax.set_xlim(-10, 62)
ax.set_xlabel('-log10(P-value)', fontsize=18, fontweight='bold', color='#1A5276')
ax.set_ylabel('')
ax.set_yticks([])
ax.spines['left'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Reference line at P=0.05
ax.axvline(x=1.3, color='#BDC3C7', linestyle='--', linewidth=2.5, alpha=0.7)
ax.text(1.3, 13, 'P=0.05', fontsize=13, color='#717D7E', ha='center')

# Legend for bubble size
for i, (size, label) in enumerate([(100, '30 genes'), (250, '60 genes'), (400, '100 genes')]):
    ax.scatter([], [], s=size, c='#717D7E', alpha=0.5, label=label)
ax.legend(title='Gene Count', loc='lower right', frameon=True, fontsize=13, title_fontsize=14)

# Note
ax.text(0.5, -1.2,
        'Note: All pathway enrichment results are exploratory and hypothesis-generating,\nrequiring experimental validation. Only nominally significant pathways shown (MAGMA P<0.05).',
        fontsize=13, color='#717D7E', style='italic', ha='center')

plt.tight_layout()
plt.savefig(r'D:\2026年\文章\因果推断\图表\Supplementary_Figure_S4_GSEA_Bubble_Plot.svg',
            format='svg', facecolor='white')
print("[SUCCESS] Supplementary Figure S4 saved")
plt.close()
