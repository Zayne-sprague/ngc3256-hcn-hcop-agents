# Zoom plots of candidate weak features with catalogue ticks (all tiers<=2, Eu<300 K, proxy>2)
import sys; sys.path.insert(0,'/workspace/scripts')
from common import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from s3c_scan import LL
T=LL[(LL.tier<=2)&(LL.upper_state_energy_K<300)&(LL.proxy>2)]
Z=[('Xb4f.spw20',253.5,254.5),('Xb47.spw22',342.4,343.6),('Xb55.spw22',164.8,166.7),('Xb47.spw20',343.8,345.55),('Xb4f.spw22',251.3,252.4),('Xb4d.spw29',253.5,254.9),('Xb55.spw16',176.2,177.0),('Xb47.spw16',353.8,354.35)]
fig,ax=plt.subplots(len(Z),1,figsize=(15,3.4*len(Z)))
for a,(t,lo,hi) in zip(ax,Z):
    z=np.load(f'/workspace/out/spectra/{t}.cube.pbcor.npz'); fr=rest_of(z['nu'])/1e9; m=(fr>lo)&(fr<hi)
    for k,c in [('N_ap','C0'),('S_ap','C3')]: a.plot(fr[m],z[k][m]*1e3,color=c,lw=.8,label=k)
    a.plot(fr[m],z['tot_r8'][m]*1e2,color='k',lw=.8,label='tot_r8 x0.1')
    ym=a.get_ylim()[1]
    for _,r in T[(T.freq_GHz>lo)&(T.freq_GHz<hi)].iterrows():
        a.axvline(r.freq_GHz,color='g',alpha=.25,lw=.6); a.text(r.freq_GHz,ym,r['name'].split()[0]+str(r.resolved_QNs).strip()[:14],rotation=90,fontsize=5,va='top')
    a.set_title(t,fontsize=8); a.legend(fontsize=6)
plt.tight_layout(); plt.savefig('/workspace/out/zoom_candidates.png',dpi=60)
