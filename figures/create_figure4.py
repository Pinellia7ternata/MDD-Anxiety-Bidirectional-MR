#!/usr/bin/env python3
"""
Figure 4: Shared genetic architecture and functional annotation summary
MAXIMUM SPACING - Clean SVG with editable text
"""
import matplotlib
matplotlib.use('SVG')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import matplotlib.gridspec as gridspec

# Configure for editable text
plt.rcParams['svg.fonttype'] = 'none'

# Enlarged figure - MAXIMUM SPACING
fig = plt.figure(figsize=(24, 20))
gs = gridspec.GridSpec(2, 2, figure=fig, height_ratios=[1.3, 1.8], wspace=0.4, hspace=0.5)
fig.patch.set_facecolor('white')

# Main title
fig.suptitle('Figure 4. Summary of shared genetic architecture and exploratory functional annotation\nbetween MDD and anxiety disorders',
            fontsize=20, fontweight='bold', color='#1A5276', y=0.995)

# Panel A: Colocalization donut chart
ax_a = fig.add_subplot(gs[0, 0])
ax_a.set_facecolor('white')

sizes = [48, 52]
colors = ['#1E8449', '#BDC3C7']
wedges, texts = ax_a.pie(sizes, colors=colors, startangle=90,
                        wedgeprops=dict(width=0.7, edgecolor='white', linewidth=4))
ax_a.text(0, 0, '375\nLoci', ha='center', va='center', fontsize=22, fontweight='bold', color='#2C3E50')
ax_a.set_title('Panel A: Colocalization\n(375 Loci Tested)', fontsize=18, fontweight='bold',
               color='#1A5276', pad=20)

# Legend
ax_a.legend(['PP.H4 > 0.5 (48%)', 'PP.H4 <= 0.5 (52%)'], loc='lower center',
            fontsize=13, frameon=True)

# Panel B: LDSC genetic correlation
ax_b = fig.add_subplot(gs[0, 1])
ax_b.set_xlim(0, 16)
ax_b.set_ylim(0, 8)
ax_b.axis('off')

# LDSC box
rect = FancyBboxPatch((1.5, 2.0), 13, 4.5, boxstyle="round,pad=0.2,rounding_size=0.5",
                       facecolor='#2471A3', edgecolor='#1A5276', linewidth=5)
ax_b.add_patch(rect)
ax_b.text(8, 5.0, 'rg = 0.263', ha='center', va='center', fontsize=34,
          fontweight='bold', color='white')

# SE and P-value
ax_b.text(8, 1.4, 'SE = 0.031', ha='center', va='center', fontsize=16, color='#717D7E')
ax_b.text(8, 0.3, 'P = 2.1 x 10^-14', ha='center', va='center', fontsize=16,
          color='#A93226', fontweight='bold')

ax_b.set_title('Panel B: LDSC Genetic Correlation', fontsize=18, fontweight='bold',
               color='#1A5276', pad=20)

# Panel C: High-confidence loci summary
ax_c = fig.add_subplot(gs[1, 0])
ax_c.set_xlim(0, 16)
ax_c.set_ylim(0, 7)
ax_c.axis('off')

# Summary table title
ax_c.text(8, 6.2, 'Panel C: Colocalization Summary', fontsize=18, fontweight='bold',
          ha='center', color='#1A5276')

# Table headers
headers = ['PP.H4 Threshold', 'Number of Loci', 'Percentage']
col_widths = [6, 5.5, 5]
row_height = 1.2
start_x = 0.5
start_y = 4.8

x = start_x
for j, (header, width) in enumerate(zip(headers, col_widths)):
    rect = FancyBboxPatch((x, start_y), width, row_height,
                           boxstyle="round,pad=0.03,rounding_size=0.1",
                           facecolor='#1A5276', edgecolor='#1A5276')
    ax_c.add_patch(rect)
    ax_c.text(x + width/2, start_y + row_height/2, header, ha='center', va='center',
             fontsize=15, fontweight='bold', color='white')
    x += width

# Data rows
rows = [
    ('PP.H4 > 0.95', '11', '2.9%'),
    ('PP.H4 > 0.7', '85', '22.7%'),
    ('PP.H4 > 0.5', '180', '48.0%'),
]
for i, row_data in enumerate(rows):
    y = start_y - (i+1) * row_height
    x = start_x
    row_color = '#EBF5FB' if i % 2 == 0 else 'white'
    for j, (cell, width) in enumerate(zip(row_data, col_widths)):
        rect = FancyBboxPatch((x, y), width, row_height,
                               boxstyle="round,pad=0.03,rounding_size=0.05",
                               facecolor=row_color, edgecolor='#BDC3C7', linewidth=1)
        ax_c.add_patch(rect)
        weight = 'bold' if j == 0 else 'normal'
        color = '#1E8449' if '0.95' in cell else '#2C3E50'
        ax_c.text(x + width/2, y + row_height/2, cell, ha='center', va='center',
                 fontsize=14, fontweight=weight, color=color)
        x += width

# Panel D: Functional annotation summary table
ax_d = fig.add_subplot(gs[1, 1])
ax_d.set_xlim(0, 16)
ax_d.set_ylim(0, 7)
ax_d.axis('off')

# Table title
ax_d.text(8, 6.2, 'Panel D: Functional Annotation Summary', fontsize=18, fontweight='bold',
          ha='center', color='#1A5276')

# Table headers
headers2 = ['Category', 'MDD', 'Anxiety', 'Shared']
col_widths2 = [5.5, 3.2, 3.2, 4]
row_height2 = 1.1
start_x2 = 0.3
start_y2 = 4.8

x = start_x2
for j, (header, width) in enumerate(zip(headers2, col_widths2)):
    rect = FancyBboxPatch((x, start_y2), width, row_height2,
                           boxstyle="round,pad=0.03,rounding_size=0.1",
                           facecolor='#1A5276', edgecolor='#1A5276')
    ax_d.add_patch(rect)
    ax_d.text(x + width/2, start_y2 + row_height2/2, header, ha='center', va='center',
             fontsize=14, fontweight='bold', color='white')
    x += width

# Data rows
rows2 = [
    ('Lead SNPs', '327', '68', '-'),
    ('Significant Genes', '5,625', '3,627', '2,128 (58.7%)'),
    ('Pathways', '64', '23', '18'),
    ('Brain eQTL Genes', '252', '85', '72 (85%)'),
]
for i, row_data in enumerate(rows2):
    y = start_y2 - (i+1) * row_height2
    x = start_x2
    row_color = '#EBF5FB' if i % 2 == 0 else 'white'
    for j, (cell, width) in enumerate(zip(row_data, col_widths2)):
        rect = FancyBboxPatch((x, y), width, row_height2,
                               boxstyle="round,pad=0.03,rounding_size=0.05",
                               facecolor=row_color, edgecolor='#BDC3C7', linewidth=1)
        ax_d.add_patch(rect)
        weight = 'bold' if j == 3 else 'normal'
        color = '#1E8449' if j == 3 else '#2C3E50'
        ax_d.text(x + width/2, y + row_height2/2, cell, ha='center', va='center',
                 fontsize=13, fontweight=weight, color=color)
        x += width

# Note at bottom
fig.text(0.5, 0.04, 'Note: All functional annotation results are exploratory and hypothesis-generating, requiring experimental validation.',
         fontsize=13, color='#717D7E', style='italic', ha='center')

# Legend
fig.text(0.5, 0.09, 'Blue: MDD results    Red: Anxiety results    Green: Shared/Overlap',
         fontsize=13, color='#2C3E50', ha='center')

plt.savefig(r'D:\2026年\文章\因果推断\图表\Figure_4_Shared_Genetic_Architecture.svg',
            format='svg', facecolor='white')
print("[SUCCESS] Figure 4 saved")
plt.close()
