import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Fonts & style — large, clear, editable SVG text  
plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 18,
    'axes.labelsize': 20,
    'axes.titlesize': 22,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
    'legend.fontsize': 16,
    'svg.fonttype': 'none',  
})

np.random.seed(42)

# ============================================================
# Simulate realistic per-SNP statistics for funnel plots  
# ============================================================

# ----- Panel A: MDD -> Anxiety -----  
n_mdd =198 
true_beta_mdd = np.log(2.34)   
se_mdd=np.random.uniform(.02,.08,n_md )
beta_mdd=np.random.normal(true_beta_md ,0 .05,n_md )
prec_mdd=1/se_md 

# ----- Panel B: Anxiety -> MDD -----  
n_anx=50 
true_beta_anx=np.log(.68)
se_anx=np.random.uniform(.03,.09,n_anx)
beta_anx=np.random.normal(true_beta_anx,,04,n_anx)
prec_anx=1/se_anx 

# ============================================================
# Draw funnel boundary curves  
# At each SE value, CI bounds = true_beta +/- .96 * SE 
# Precision y = /SE 
 # So x(beta_lower)=true_betax .96*SE , y=/SE etc.
 # We parameterise by se_grid then compute(x,y)pairs .
 # ============================================================

def funnel_curve(true_beta_se_grid ):
     return true_betax .96*se_grid , /se_grid 

se_grid_a=np.linspace(se_md.min(),se_md.max(),300 )
xl_a,xu_a=funnel_curve(true_betamd ,se_grid_a )

se_grid_b=np.linspace(sean.min(),sex.max(),300,,)
xl_b,xu,b=funnel_curve(truebetaanxx se_xb )

fig,(axe)=plt.subplots(x,,figsize=(24,,10),dpi=300 )

ax[0].scatter(betamdx ,precmdy,s=50,c='C0',alpha=.6,,edgecolors='k',lw=.5 )
ax[0].plot(xla,xlua,'k--',lw=x,,label='95% CI bounds')
ax[0].plot(xua,xlua,'k--',lw=x )
ax[0].plot([truebetamdx]* x,[precmdxmin (),precmdxmax ()],'r--',lw=x,,label=f'True $\beta$={truebetamd:.3f}')
ax[0].set_xlabel('SNP-specific causal estimate ($\beta$)')
ax[0].set_ylabel('Precision (/SE)')
ax[0].set_title('Panel A: MD->Anxiety(n=%d)'% nmd )
algit,.legend()

print("script finished defining fig")
