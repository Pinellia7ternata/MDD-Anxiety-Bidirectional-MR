#!/usr/bin/env python3
"""
Figure 2: Forest plot of bidirectional MR results
MAXIMUM SPACING - Clean SVG with editable text
"""
import matplotlib
matplotlib.use('SVG')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Configure for editable text
plt.rcParams['svg.fonttype'] = 'none'

# Enlarged figure with MAXIMUM spacing
fig, ax = plt.subplots(figsize=(22, 16))
ax.set_facecolor('white')
fig.patch.set_facecolor('white')

# Clear and set up
ax.set_xlim(0, 22)
ax.set_ylim(-1, 16)
ax.axis('off')

# Title - LARGER
ax.text(11, 15.2, 'Figure 2', fontsize=20, fontweight='bold', color='#1A5276', ha='center')
ax.text(11, 14.5, 'Forest plot of bidirectional Mendelian randomization estimates\nbetween MDD and anxiety disorders',
        fontsize=15, color='#2C3E50', ha='center', style='italic')

# Data
mdd_anxiety = [
    ('IVW', 2.34, 2.23, 2.46, 1.66e-257, '#2471A3'),
    ('MR-Egger', 2.30, 1.84, 2.87, 7.85e-12, '#5DADE2'),
    ('Weighted Median', 2.28, 2.14, 2.42, 9.69e-155, '#85929E'),
]

anxiety_mdd = [
    ('IVW', 1.68, 1.60, 1.76, 1.97e-104, '#A93226'),
    ('MR-Egger', 1.28, 1.06, 1.55, 1.34e-2, '#EC7063'),
    ('Weighted Median', 1.63, 1.55, 1.72, 3.96e-72, '#B7950B'),
]

# Function to draw forest plot with MAXIMUM SPACING
def draw_forest_group_fixed(ax, data, start_y, group_label, color_main, x_offset=0):
    # Group label
    ax.text(x_offset + 0.5, start_y + 1.0, group_label, fontsize=18, fontweight='bold', 
            color=color_main, ha='left')
    
    # Column headers
    ax.text(x_offset + 4.5, start_y + 0.4, 'OR', fontsize=14, fontweight='bold', ha='center', color='#2C3E50')
    ax.text(x_offset + 6.5, start_y + 0.4, '95% CI', fontsize=14, fontweight='bold', ha='center', color='#2C3E50')
    ax.text(x_offset + 9.0, start_y + 0.4, 'P-value', fontsize=14, fontweight='bold', ha='center', color='#2C3E50')
    
    for i, (method, or_val, ci_low, ci_high, pval, color) in enumerate(data):
        y = start_y - i * 2.8  # MAXIMUM SPACING
        
        # Method name
        ax.text(x_offset + 0.5, y, method, fontsize=15, va='center', ha='left', 
                color='#2C3E50')
        
        # OR value
        ax.text(x_offset + 4.5, y, f'{or_val:.2f}', fontsize=15, va='center', ha='center', 
                fontweight='bold', color='#2C3E50')
        
        # CI in parentheses
        ax.text(x_offset + 6.5, y, f'({ci_low:.2f}-{ci_high:.2f})', fontsize=14, 
                va='center', ha='center', color='#717D7E')
        
        # P-value
        if pval < 1e-100:
            ptext = 'P < 1e-100'
        else:
            ptext = f'P = {pval:.2e}'
        ax.text(x_offset + 9.0, y, ptext, fontsize=13, va='center', ha='center', 
                color='#717D7E')
        
        # Draw forest plot bar on right side
        log_or = np.log(or_val)
        log_low = np.log(ci_low)
        log_high = np.log(ci_high)
        
        # Map to display coordinates
        log_range = np.log(4.0)
        display_x = 11 + (log_or / log_range) * 10
        display_low = 11 + (log_low / log_range) * 10
        display_high = 11 + (log_high / log_range) * 10
        
        # Draw CI bar
        ax.plot([display_low, display_high], [y, y], color=color, linewidth=5)
        ax.plot([display_low, display_low], [y-0.3, y+0.3], color=color, linewidth=3)
        ax.plot([display_high, display_high], [y-0.3, y+0.3], color=color, linewidth=3)
        ax.plot(display_x, y, 'o', color=color, markersize=16, 
                markeredgecolor='white', markeredgewidth=2)
    
    return start_y - len(data) * 2.8 - 1.5

# Draw reference line
ax.axvline(x=11, color='#717D7E', linestyle='--', linewidth=2, alpha=0.7)
ax.text(11.3, 13.0, 'OR=1', fontsize=14, color='#717D7E', va='bottom', fontweight='bold')

# Draw MDD -> Anxiety
y_end = draw_forest_group_fixed(ax, mdd_anxiety, 12.5, 'MDD to Anxiety', '#2471A3', x_offset=0)

# Draw Anxiety -> MDD
y_end2 = draw_forest_group_fixed(ax, anxiety_mdd, 6.0, 'Anxiety to MDD', '#A93226', x_offset=0)

# X-axis
ax.text(16, -0.5, 'Odds Ratio (OR)', fontsize=16, ha='center', color='#2C3E50', fontweight='bold')
for val, label in [(0.5, '0.5'), (1, '1'), (1.5, '1.5'), (2, '2'), (2.5, '2.5'), (3, '3'), (3.5, '3.5')]:
    log_val = np.log(val)
    display_x = 11 + (log_val / np.log(4.0)) * 10
    ax.text(display_x, -0.9, label, fontsize=13, ha='center', color='#717D7E')
    ax.plot([display_x, display_x], [-0.6, -1.1], color='#BDC3C7', linewidth=1)

# Legend
legend_elements = [
    mpatches.Patch(facecolor='#2471A3', label='MDD to Anxiety'),
    mpatches.Patch(facecolor='#A93226', label='Anxiety to MDD'),
]
ax.legend(handles=legend_elements, loc='upper right', frameon=True, fontsize=14)

plt.tight_layout()
plt.savefig(r'D:\2026年\文章\因果推断\图表\Figure_2_Bidirectional_MR_Forest_Plot.svg', 
            format='svg', facecolor='white')
print("[SUCCESS] Figure 2 saved")
plt.close()
