import matplotlib.ticker, json, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
r=json.load(open('/workspace/paper_answers/q01/q01_values.json'))
plt.rcParams.update({'font.size':11})
fig,ax=plt.subplots(figsize=(6.4,6.2))
bs=[0.3,0.6,1.05,1.5,2.2]
for est,c,m,lab,dx in (('ann025','#0072B2','o',r'Xb4b convolved, core $r<0.25^{\prime\prime}$ (Table C1)',0.97),
                       ('ann1','#009E73','s',r'Xb4b convolved, $r<1^{\prime\prime}$ (Table C1)',1.0),
                       ('ap','#D55E00','^',r'Xb4b convolved, aperture $r$=beam FWHM (uvtaper estimator)',1.03)):
    ax.errorbar(np.array(bs)*dx,[r[str(b)][est]['NS'] for b in bs],[r[str(b)][est]['err'] for b in bs],color=c,marker=m,ms=6,capsize=2,lw=1.3,label=lab)
ax.errorbar([1.5,2.2],[1.532,1.441],[0.026,0.023],color='k',marker='D',ls='none',ms=7,capsize=2,label='12m+7m uvtaper, 3-2 (aperture)')
ax.errorbar([1.5,2.2],[1.578,1.451],[0.027,0.027],mfc='none',color='k',marker='D',ls='none',ms=7,capsize=2,label='12m-only uvtaper, 3-2')
ax.errorbar([1.5,2.2,2.5],[1.616,1.467,1.56],[0.059,0.022,0.03],color='#CC79A7',marker='v',ls='none',ms=7,capsize=2,label='HCN/HCO$^+$ 1-0 (Band 3)')
ax.plot([1.5],[2.03],marker='*',color='gray',ls='none',ms=11,label=r'team-a 3-2 Xb4b, 1.5$^{\prime\prime}$ (p0065)')
ax.axvline(0.246,color='gray',ls=':',lw=1); ax.text(0.25,1.28,'native\n0.25"',fontsize=8,color='gray',va='bottom',ha='left')
ax.set_xscale('log'); ax.set_xticks([0.3,0.6,1,1.5,2.2]); ax.set_xticklabels(['0.3','0.6','1','1.5','2.2'])
ax.set_xlim(0.2,3); ax.set_ylim(1.3,2.5); ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
ax.set_xlabel('Circular beam FWHM (arcsec)'); ax.set_ylabel(r'N/S ratio of HCN/HCO$^+$ (moment-0)')
sec=ax.secondary_xaxis('top',functions=(lambda x:x*194,lambda x:x/194)); sec.set_xticks([60,100,200,300,430]); sec.xaxis.set_minor_locator(matplotlib.ticker.NullLocator()); sec.set_xticklabels(['60','100','200','300','430']); sec.set_xlabel('Linear scale (pc, 194 pc/arcsec)')
ax.legend(fontsize=7.5,loc='upper center',bbox_to_anchor=(0.5,-0.15),frameon=False,ncol=1); fig.tight_layout()
fig.savefig('/workspace/paper_answers/q01/dilution_NS_vs_beam.png',dpi=180,bbox_inches='tight')
