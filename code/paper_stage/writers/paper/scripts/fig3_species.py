# Fig. 3: flux-complete J=1-0 N/S contrasts of line/HCO+ ratios, 2.5" beam (Table 1).
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig,ax=plt.subplots(figsize=(4.8,4))
lab=['HCN/HCO$^+$','HNC/HCO$^+$','CN/HCO$^+$']; v=[1.56,1.90,2.10]; e=[.03,.11,.04]
ax.bar(range(3),v,.5,yerr=e,color='#0072B2',capsize=3); ax.axhline(1,color='gray',ls=':')
ax.set_xticks(range(3)); ax.set_xticklabels(lab); ax.set_ylim(0,2.6)
ax.set_ylabel('N/S contrast of the ratio (statistical errors)')
ax.set_title('J=1-0, 2.5", flux-complete 12m+7m',fontsize=9)
fig.tight_layout(); fig.savefig('../figs/fig3_species.png',dpi=160)
