# Overview figure per observation set: tot_r8 / N / S spectra vs rest frequency (at VSYS) with catalogue ticks (tier-1 ground-state lines)
import sys; sys.path.insert(0,'/workspace/scripts')
from common import *
import pandas as pd, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from s3c_scan import LL
T=LL[(LL.tier==1)&(LL.upper_state_energy_K<300)]
T=T[T.sijmu2.fillna(30)*np.exp(-T.upper_state_energy_K/40)>3]
sets={'Xb55':['spw22','spw20','spw16','spw18'],'Xb47':['spw22','spw20','spw16','spw18'],'Xb4f':['spw22','spw20','spw18','spw16'],'Xb4d':['spw31','spw29','spw27','spw25']}
for ous,sp in sets.items():
    fig,ax=plt.subplots(len(sp),1,figsize=(16,3.2*len(sp)))
    for a,s in zip(ax,sp):
        z=np.load(f'/workspace/out/spectra/{ous}.{s}.cube.pbcor.npz'); fr=rest_of(z['nu'])/1e9
        for k,c in [('tot_r8','k'),('N_ap','C0'),('S_ap','C3')]:
            y=z[k]; sc=1e3 if k!='tot_r8' else 1e3
            a.plot(fr,y*sc/(10 if k=='tot_r8' else 1),color=c,lw=.7,label=k+(' /10 [mJy]' if k=='tot_r8' else ' [mJy/b]'))
        ym=a.get_ylim()[1]
        for _,r in T[(T.freq_GHz>fr.min())&(T.freq_GHz<fr.max())].iterrows():
            a.axvline(r.freq_GHz,color='g',alpha=.3,lw=.7); a.text(r.freq_GHz,ym*0.9,r['name'].split()[0]+' '+str(r.resolved_QNs).strip()[:8],rotation=90,fontsize=5,va='top')
        a.set_title(f'{ous} {s}',fontsize=8); a.legend(fontsize=6); a.set_ylim(-0.1*ym,ym)
    plt.tight_layout(); plt.savefig(f'/workspace/out/overview_{ous}.png',dpi=55); plt.close()
