#!/usr/bin/env python3
"""
Figure 3: Forest plot of full MVMR results
MAXIMUM SPACING - Clean SVG with editable text
"""
import matplotlib
matplotlib.use('SVG')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch

# Configure for editable text
plt.rcParams['svg.fonttype'] = 'none'

# Enlarged figure with two panels - MAXIMUM SPACING
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(24, 12))
fig.patch.set_facecolor('white')
fig.suptitle('Figure 3. Full multivariable Mendelian randomization estimates\nadjusted for BMI, education, and smoking initiation',
            fontsize=18, fontweight='bold', color='#1A5276', y=1.06)

# Data for Panel A: MDD -> Anxiety
mdd_mvmr = [
    ('MDD (Primary)', 1.15, 1.09, 1.22, 6.71e-7, '#2471A3'),
    ('BMI', 1.32, 1.13, 1.55, 4.49e-4, '#5DADE2'),
    ('Education', 1.26, 1.04, 1.53, 0.0205, '#1ABC9C'),
    ('Smoking Initiation', 3.26, 2.93, 3.63, 4.18e-104, '#A93226'),
]

# Data for Panel B: Anxiety -> MDD
anx_mvmr = [
    ('Anxiety (Primary)', 1.10, 1.05, 1.15, 7.90e-6, '#A93226'),
    ('BMI', 1.00, 0.85, 1.19, 0.983, '#EC7063'),
    ('Education', 1.77, 1.48, 2.11, 1.89e-10, '#F39C12'),
    ('Smoking Initiation', 1.02, 0.91, 1.15, 0.742, '#717D7E'),
]

def draw_mvmr_panel_fixed(ax, data, title, cond_f, x_start=0.5, x_end=5.5):
    # Clear and setup
    ax.cla()
    ax.set_facecolor('white')
    
    # Set limits with MAXIMUM MARGINS
    ax.set_xlim(x_start - 0.8, x_end + 3.0)
    ax.set_ylim(-1.5, 6.0)
    
    # Reference line at OR = 1
    ax.axvline(x=1.0, color='#717D7E', linestyle='--', linewidth=3, alpha=0.7)
    ax.text(1.15, 5.2, 'OR=1', fontsize=14, color='#717D7E', va='bottom', fontweight='bold')
    
    # Panel title
    ax.text((x_start + x_end) / 2, 5.5, title, fontsize=18, fontweight='bold', 
            ha='center', color='#1A5276')
    ax.text((x_start + x_end) / 2, 5.0, f'Conditional F = {cond_f}', 
            fontsize=15, color='#2471A3', ha='center', style='italic')
    
    for i, (exposure, or_val, ci_low, ci_high, pval, color) in enumerate(data):
        y = 4.2 - i * 1.3  # MAXIMUM SPACING
        
        # Exposure name
        ax.text(x_start, y, exposure, fontsize=15, va='center', ha='left', 
                color='#2C3E50', fontweight='bold' if 'Primary' in exposure else 'normal')
        
        # OR value
        ax.text(x_start + 2.0, y, f'{or_val:.2f}', fontsize=16, va='center', ha='center', 
                fontweight='bold', color='#2C3E50')
        
        # CI
        ax.text(x_start + 3.0, y, f'({ci_low:.2f}-{ci_high:.2f})', fontsize=14, 
                va='center', ha='left', color='#717D7E')
        
        # P-value
        ax.text(x_start + 4.5, y, f'P={pval:.2e}', fontsize=14, va='center', 
                ha='left', color='#717D7E')
        
        # Forest plot
        log_or = np.log(or_val)
        log_low = np.log(ci_low)
        log_high = np.log(ci_high)
        
        # Map to plot range
        log_range = np.log(x_end)
        display_x = x_start + (log_or / log_range) * (x_end - x_start)
        display_low = x_start + (log_low / log_range) * (x_end - x_start)
        display_high = x_start + (log_high / log_range) * (x_end - x_start)
        
        # Draw CI
        ax.plot([display_low, display_high], [y, y], color=color, linewidth=5)
        ax.plot([display_low, display_low], [y-0.2, y+0.2], color=color, linewidth=3)
        ax.plot([display_high, display_high], [y-0.2, y+0.2], color=color, linewidth=3)
        ax.plot(display_x, y, 'o', color=color, markersize=16, 
               markeredgecolor='white', markeredgewidth=2)
        
        # Significance marker
        if pval < 0.05:
            ax.text(x_end + 0.5, y, '*', fontsize=26, color='red', va='center', ha='left')
    
    # X-axis
    ax.set_xlabel('Odds Ratio (OR)', fontsize=16, fontweight='bold', color='#1A5276')
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)

# Draw panels
draw_mvmr_panel_fixed(ax1, mdd_mvmr, 'Panel A: MDD to Anxiety', 15.6, x_start=0.5, x_end=5.5)
draw_mvmr_panel_fixed(ax2, anx_mvmr, 'Panel B: Anxiety to MDD', 12.9, x_start=0.7, x_end=3.0)

# Add note
fig.text(0.5, -0.08, '* Significant at P < 0.05. CI = confidence interval.',
        fontsize=13, color='#717D7E', style='italic', ha='center')

plt.tight_layout()
plt.savefig(r'D:\2026年\文章\因果推断\图表\Figure_3_MVMR_Forest_Plot.svg', 
            format='svg', facecolor='white')
print("[SUCCESS] Figure 3 saved")
plt.close()
