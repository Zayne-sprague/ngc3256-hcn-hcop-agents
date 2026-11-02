# Fig. 4: flux-complete comparisons only (values from Sec. 3.2 / Table 1, Table A2)
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, numpy as np
c=dict(b='#0072B2',o='#E69F00')
plt.rcParams.update({'font.size':11})
fig,ax=plt.subplots(1,2,figsize=(10,4.4),gridspec_kw={'width_ratios':[1.2,1]})
a=ax[0]; x=np.arange(3); w=.36
a.bar(x-w/2,[1.46,1.62,1.49],w,yerr=[.02,.06,.07],color=c['b'],capsize=3,label='J = 1-0')
a.bar(x+w/2,[1.48,1.60,1.64],w,yerr=[.04,.07,.13],color=c['o'],capsize=3,label='J = 3-2')
a.axhline(1,color='gray',ls=':'); a.set_xticks(x); a.set_xticklabels(['2.2"\n(430 pc)','1.5"\n(290 pc)','2.5"\n(490 pc)'])
a.set_xlabel('matched circular beam, flux-complete 12m+7m'); a.set_ylabel('N/S contrast of HCN/HCO$^+$ ratio')
a.set_ylim(0,2.6); a.legend(loc='upper left'); a.set_title('(a) 1-0 and 3-2 (300-500 pc blend)')
b=ax[1]; v=[1.56,1.90,2.10]; e=[.03,.11,.04]
b.bar(range(3),v,.5,yerr=e,color=c['b'],capsize=3); b.axhline(1,color='gray',ls=':')
b.set_xticks(range(3)); b.set_xticklabels(['HCN/HCO$^+$','HNC/HCO$^+$','CN/HCO$^+$']); b.set_ylim(0,2.6)
b.set_ylabel('N/S contrast (statistical errors)'); b.set_title('(b) species, J=1-0, 2.5", flux-complete')
plt.tight_layout(); plt.savefig('figs/fig4_contrast_J_species.png',dpi=150)
