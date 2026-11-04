"""r09 figure from r09_profile.json (made by r09_profile.py). usage: python3 -I r09_fig.py"""
import json, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10,'axes.labelsize':11})
R=json.load(open('/workspace/paper_answers/r09/r09_profile.json'))
ANN=['0-0.25','0.25-0.5','0.5-1.0','1.0-2.0','2.0-3.0']
def rad(a): lo,hi=map(float,a.split('-')); return 0.5*(lo+hi),0.5*(hi-lo)
SER=[('Xb4b_0.26','Xb4b, 0.26" beam','#0072B2','o',-1),('Xb4b_1.0','Xb4b, 1.0" beam','#E69F00','s',0),('Xb4d_1.0','Xb4d, 1.0" beam','#009E73','^',1)]
fig,ax=plt.subplots(1,3,figsize=(13.5,4.4))
for k,lab,c,mk,sh in SER:
    for ni,(rn,ls) in enumerate((('N','-'),('S','--'))):
        x=np.array([rad(a)[0] for a in ANN])*(1+0.04*sh); y=[R[k][rn]['rows'][a]['R'] for a in ANN]; e=[R[k][rn]['rows'][a]['eR'] for a in ANN]
        ax[0].errorbar(x,y,e,color=c,marker=mk,ls=ls,mfc=c if rn=='N' else 'w',ms=5,lw=1,capsize=2,label=f'{rn}, {lab}')
for i,rn in enumerate('NS'):
    a=ax[i+1]
    for k,lab,c,mk,sh in SER:
        x=np.array([rad(q)[0] for q in ANN])*(1+0.04*sh)
        for ln,ls,fc in (('HCN','-',c),('HCO+',':',c)):
            y=[R[k][rn]['rows'][q][ln] for q in ANN]; e=[R[k][rn]['rows'][q]['e'+ln] for q in ANN]
            a.errorbar(x,y,e,color=c,marker=mk,ls=ls,mfc=c if ln=='HCN' else 'w',ms=5,lw=1,capsize=2,label=f'{ln}, {lab}')
    a.set_yscale('log'); a.set_title(f'{rn} nucleus: $T_b$ profiles (solid HCN, dotted HCO$^+$)',fontsize=10); a.set_ylabel(r'annulus-mean $\int T_b\,dv$ (K km s$^{-1}$, $\pm$350 km/s)')
ax[0].set_ylabel(r'HCN/HCO$^+$ 3-2 $T_b$ ratio'); ax[0].set_ylim(0.25,1.1); ax[0].set_title(r'Ratio (filled/solid N; open/dashed S)',fontsize=10)
for a in ax:
    a.set_xscale('log'); a.set_xlabel('annulus radius from HCO$^+$ peak (arcsec; 1" = 194 pc)'); a.set_xlim(0.09,3.4)
    for b,t in ((0.26,':'),(1.0,'--')): a.axvline(b,color='0.7',ls=t,lw=0.8,zorder=0)
    a.grid(alpha=0.2)
    s=a.secondary_xaxis('top',functions=(lambda x:x*194,lambda x:x/194)); s.set_xlabel('pc',fontsize=9); s.set_xticks([20,50,100,200,500]); s.set_xticklabels(['20','50','100','200','500']); s.minorticks_off()
ax[0].legend(fontsize=6.5,ncol=1,loc='upper right'); ax[1].legend(fontsize=6,loc='lower left')
ax[0].text(0.27,0.27,'grey lines: 0.26" and 1.0" beam FWHM',fontsize=7,color='0.4')
fig.tight_layout(); fig.savefig('/workspace/paper_answers/r09/fig_ratio_profile.png',dpi=170)
