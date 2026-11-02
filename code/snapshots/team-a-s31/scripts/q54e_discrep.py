"""Q54 r2 fix1: why HNC/HCN(1-0) differs between 2015.1.00993 (X961) and 2015.1.00412 (X145/X149)?
Same-line integrals at identical circular beams, several windows/baselines; spectra overlaid. out/q54/discrep.csv, discrep.png"""
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd, q4_lib, q54_chem as c, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from q4_lib import spectrum,kb,vel
pairs=[('HCN1-0_00993','HCN1-0_00412'),('HNC1-0_00993','HNC1-0_00412')]
rows=[];fig,ax=plt.subplots(2,2,figsize=(11,7))
for ip,(a,b) in enumerate(pairs):
  for ir,reg in enumerate('NS'):
    for key,col in ((a,'k'),(b,'r')):
        path,f0=c.L[key]; q4_lib.f=lambda *x,**k: path
        bm=2.5
        fr,s,_,nb=spectrum('x','x',reg,tb=(bm,bm,0.)); v=vel(fr,f0); y=s*kb(fr.mean(),bm,bm)
        ax[ip,ir].plot(v,y,col,lw=.8,label=f'{key} nb={nb[0]:.2f}x{nb[1]:.2f} dv={abs(np.median(np.diff(v))):.1f}')
        vc=c.VC[reg]
        for hw in (250,350,450):
          for bl in ((450,950),(600,1200)):
            c.BL=bl
            I,e=c.win(v,y,vc,vc-hw,vc+hw)
            rows.append(dict(line=key,reg=reg,hw=hw,bl=str(bl),I=I,e=e,peak=y[np.abs(v-vc)<150].max(),vmin=v.min(),vmax=v.max()))
    ax[ip,ir].set_title(f'{pairs[ip][0][:5]} {reg} 2.5"'); ax[ip,ir].set_xlim(vc-1200,vc+1200); ax[ip,ir].legend(fontsize=6)
c.BL=(450,950)
plt.savefig('/workspace/out/q54/discrep.png',dpi=90)
d=pd.DataFrame(rows).round(2); d.to_csv('/workspace/out/q54/discrep.csv',index=False)
print(d.to_string())
