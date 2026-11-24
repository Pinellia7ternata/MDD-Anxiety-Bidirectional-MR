"""Figure 3, paired under one estimator (repository copy).

Figure 3 now plots, for each of the four models, the univariable and the direct estimate computed on the SAME
variants under the SAME estimator (measurement-error-corrected, dispersion-scaled), with the ratio of betas and
the conditional F printed beside the pair, and a legend for the two marker types. The v7 figure paired a
corrected direct estimate against an uncorrected reference, which is the mismatch the reviewer found.

Repository copy: reads mvmr_paired_v8.csv from this directory and writes the SVG, vector PDF and PNG here.
"""
import os
import matplotlib
matplotlib.use('SVG')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, NullFormatter
from matplotlib.lines import Line2D
import pandas as pd

BASE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(BASE, 'mvmr_paired_v8.csv')
DOCX_PATH = None
NAVY, BLUE, RED, GRAY, MIDGRAY = '#1A5276', '#2471A3', '#A93226', '#717D7E', '#B3B6B7'
plt.rcParams['svg.fonttype'] = 'none'

d = pd.read_csv(CSV)
d = d[d['model'].str.count(r'\|') == 1].copy()          # the four main models only
PANELS = [('forward', 'MDD → anxiety disorders', BLUE),
          ('reverse', 'anxiety disorders → MDD', RED)]

fig, axes = plt.subplots(1, 2, figsize=(15.5, 5.6))
fig.patch.set_facecolor('white')
for ax, (side, title, color) in zip(axes, PANELS):
    sub = d[d['model'].str.startswith(side)]
    sub = pd.concat([sub[sub['model'].str.contains('exposure-only')],
                     sub[sub['model'].str.contains('exposure instruments')]])
    rows = list(sub.itertuples(index=False))
    for i, r in enumerate(rows[::-1]):
        y = i
        design = ('instruments: ' + ('MDD, BMI, education, smoking' if side == 'forward'
                                     else 'anxiety disorders, BMI, education, smoking')) \
            if 'exposure-only' in r.model else \
            ('covariates inherit the ' + ('MDD' if side == 'forward' else 'anxiety') + ' instruments')
        ax.text(0.012, y + 0.34, f'{design} ({int(r.n)} variants)', transform=ax.get_yaxis_transform(),
                ha='left', va='bottom', fontsize=10.5, color=color, fontweight='bold', family='Arial')
        ax.plot([r.uni_corr_lo, r.uni_corr_hi], [y, y], color=GRAY, lw=2.4, zorder=2, solid_capstyle='butt')
        ax.plot([r.uni_corr_or], [y], marker='s', markersize=9, mfc='white', mec=GRAY, mew=2.0, zorder=3)
        ax.plot([r.dir_corr_lo, r.dir_corr_hi], [y + 0.17, y + 0.17], color=color, lw=2.6, zorder=2,
                solid_capstyle='butt')
        ax.plot([r.dir_corr_or], [y + 0.17], marker='o', markersize=11, mfc=color, mec=color, zorder=3)
        ax.text(3.15, y + 0.19, f'direct {r.dir_corr_or:.2f} ({r.dir_corr_lo:.2f}-{r.dir_corr_hi:.2f})',
                fontsize=10.5, color=color, va='center', ha='left', family='Arial')
        ax.text(3.15, y + 0.02, f'univariable {r.uni_corr_or:.2f} ({r.uni_corr_lo:.2f}-{r.uni_corr_hi:.2f})',
                fontsize=10.5, color=GRAY, va='center', ha='left', family='Arial')
        ax.text(3.15, y - 0.15, f'ratio of betas {r.ratio_corr:.2f}', fontsize=9.5, color='#2C3E50',
                va='center', ha='left', family='Arial')
    ax.axvline(1.0, color=MIDGRAY, ls='--', lw=1.8, zorder=1)
    ax.text(1.03, len(rows) - 0.02, 'OR = 1', fontsize=10, color=GRAY, va='top', family='Arial')
    ax.set_xscale('log')
    ax.set_xlim(1.0, 6.5)
    ax.set_ylim(-0.45, len(rows) + 0.05)
    ax.set_yticks([])
    ax.xaxis.set_major_locator(FixedLocator([1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]))
    ax.xaxis.set_major_formatter(NullFormatter())
    ax.xaxis.set_minor_locator(FixedLocator([1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]))
    ax.xaxis.set_minor_formatter(NullFormatter())
    for x in (1.0, 1.5, 2.0, 2.5, 3.0):
        ax.text(x, -0.42, f'{x:.1f}', ha='center', va='top', fontsize=11, color='#2C3E50', family='Arial')
    ax.set_xlabel('Odds ratio per unit increase in liability (log scale)', fontsize=11.5, color='#2C3E50')
    ax.set_title(title, fontsize=14.5, fontweight='bold', color=NAVY, pad=10)
    for spine in ('top', 'right', 'left'):
        ax.spines[spine].set_visible(False)
    ax.spines['bottom'].set_color(MIDGRAY)
    ax.tick_params(axis='x', colors=MIDGRAY)
    ax.grid(axis='x', color='#EAECEE', lw=0.9, zorder=0)

handles = [Line2D([0], [0], marker='s', color='w', markerfacecolor='white', markeredgecolor=GRAY,
                  markersize=9, markeredgewidth=2, label='univariable IVW on these variants'),
           Line2D([0], [0], marker='o', color='w', markerfacecolor=BLUE, markersize=11, label='direct (multivariable) estimate')]
fig.legend(handles=handles, loc='lower center', ncol=2, fontsize=10, frameon=False, handlelength=1.6,
           bbox_to_anchor=(0.5, 0.008))
fig.suptitle('Figure 3  Direct multivariable estimates against the univariable estimate on the same variants '
             'and the same estimator', fontsize=13.5, fontweight='bold', color=NAVY, x=0.012, ha='left', y=0.985)
axes[0].text(-0.05, 1.04, 'A', transform=axes[0].transAxes, fontsize=15, fontweight='bold', color=NAVY)
axes[1].text(-0.05, 1.04, 'B', transform=axes[1].transAxes, fontsize=15, fontweight='bold', color=NAVY)
plt.subplots_adjust(left=0.04, right=0.99, top=0.79, bottom=0.19, wspace=0.40)
F3 = os.path.join(BASE, 'Figure_3_MVMR_Forest_Plot')
fig.savefig(F3 + '.svg', format='svg', facecolor='white')
plt.close(fig)



def export_and_swap(which, svg_path, docx_target, shape_index):
    """Repository copy: exports the vector PDF and the raster PNG next to this script and stops.

    The submission copy additionally copies the PDF into the figures folder of the submission package and swaps
    the PNG into the manuscript docx; neither belongs in the repository.
    """
    import pymupdf
    doc = pymupdf.open(stream=open(svg_path, encoding='utf-8').read().encode('utf-8'), filetype='svg')
    pdf = doc.convert_to_pdf()
    stem = svg_path[:-4]
    open(stem + '.pdf', 'wb').write(pdf)
    pg = pymupdf.open(stream=pdf, filetype='pdf')[0]
    pix = pg.get_pixmap(matrix=pymupdf.Matrix(2600 / pg.rect.width, 2600 / pg.rect.width), alpha=False)
    pix.save(stem + '.png')
    print(f'Figure {which}: PDF ({os.path.getsize(stem + ".pdf"):,} B) and PNG '
          f'({pix.width}x{pix.height}, {os.path.getsize(stem + ".png"):,} B) written next to this script')



print('wrote', F3 + '.svg')
export_and_swap(3, F3 + '.svg', None, 2)
