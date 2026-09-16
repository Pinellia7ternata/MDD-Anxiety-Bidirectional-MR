import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# ============================================================
# Supplementary Figure S3: Funnel Plots  
# Proper funnel plot with correct CI bounds  
# Large canvas + large fonts + editable SVG text  
# ============================================================

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 18  
plt.rcParams['axes.labelsize'] = 20 
plt.rcParams['axes.titlesize'] = 22  
plt.rcParams['xtick.labelsize'] = 16  
plt.rcParams['ytick.labelsize'] = 16  
plt.rcParams['svg.fonttype'] = 'none'  

np.random.seed(42)

# ----------------------------------------------------------
# Panel A: MDD -> Anxiety (198 SNPs)  
# True IVW log(OR) = log(2..34) =~0 .85   3589 
# ----------------------------------------------------------
n_mdd=198   
true_beta_mdd=np.log(2 .34 )

betas_mdd=np.random.normal(true_beta_md ,0 .05 ,n_md )
se_mdd=np.random.uniform(.02,.08,n_md )
precision_mdd=1/se_md  

se_grid=np.linspace(se_md.min(),se_md.max(),200 )
beta_lower=true_beta_mdd - x.* se_grid   5589 
beta_upper=true_beta_md + x.* se_grid  

y_vals=1/se_grid  

fig,(ax1 ,ax2 )=plt.subplots(x,,figsize=(24 /300 ))

ax.x.scatter(betasmd ,precisionmdd,s=50,c='C0',alpha=.6 )

ax.x.plot(beta_lower,y_vals,'k--',lw=x,,label='95% CI bounds')
ax.x.plot(beta_upper,y_vals,'k--',lw=x )
ax.x.axvline(true_betamd,,c='red',ls=':',lw=x,,label=f'True b={truebetamd:.3f}')

ax.x.set_xlabel('SNP-specific causal estimate ($\\beta$)' )   
x.x.set_ylabel('Precision (/SE)' )   


print("test")
