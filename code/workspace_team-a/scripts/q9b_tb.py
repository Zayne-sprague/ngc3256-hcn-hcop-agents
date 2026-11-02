# Q9b fix 4b: per-channel Tb(HCN) vs Tb(HCO+) (aperture-mean brightness temperature, r<0.25" about the HCO+ peak, Xb4b convolved to 0.30"; also Xb4d@1.05" r<0.5")
# with Tb = 1.222e6 S[Jy/beam]/(nu_GHz^2 theta_maj theta_min) (Rayleigh-Jeans); linear fit through origin & ratio of Tb per channel; peak Tb vs Tkin ceilings. python3 -I q9b_tb.py
import sys; sys.path.insert(0,'/workspace/scripts')
from q9b_lib import *
import pandas as pd, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
NU={'HCN':265.886434*(1-2780/C),'HCO+':267.557626*(1-2780/C)}    # sky GHz at vsys
rows=[]; fig,ax=plt.subplots(1,2,figsize=(10,4.5))
for (ous,tb,ap,k) in (('Xb4b',(0.30,)*2+(0.,),0.25,0),('Xb4d',(1.05,)*2+(0.,),0.5,1)):
    for reg,col in (('N','#1f77b4'),('S','#d95f02')):
        x0,y0=centre(ous,tb,reg); T={}
        for ln in ('HCN','HCO+'):
            c=cube(ous,ln); s=smcube(ous,ln,tb); m=apmask(c,x0,y0,0,ap); I=s[:,m].mean(1)      # mean Jy/beam(tb)
            T[ln]=(c.v,1.222e6*I/(NU[ln]**2*tb[0]*tb[1]))
            pk=np.argmax(I); print(f'{ous}@{tb[0]} {reg} {ln}: peak Tb(mean r<{ap}") {T[ln][1][pk]:.1f} K at v={c.v[pk]:.0f}; channel rms {1.4826*np.median(np.abs(T[ln][1][c.bm]-np.median(T[ln][1][c.bm]))):.2f} K')
        vO=np.interp(T['HCN'][0],T['HCO+'][0][::-1],T['HCO+'][1][::-1]) if T['HCO+'][0][0]>T['HCO+'][0][-1] else np.interp(T['HCN'][0],*T['HCO+'])
        v=T['HCN'][0]; h=T['HCN'][1]; o=vO; sel=(np.abs(v)<250)&(o>0.1*o.max()); sl=(h[sel]*o[sel]).sum()/(o[sel]**2).sum()
        print(f'  Tb(HCN)=slope*Tb(HCO+) through origin, channels with Tb(HCO+)>10% peak: slope {sl:.3f}; ratio in channel Tb(HCO+)>50% peak {np.median(h[o>0.5*o.max()]/o[o>0.5*o.max()]):.3f}')
        if k==0: ax[0].scatter(o[np.abs(v)<250],h[np.abs(v)<250],s=8,color=col,label=f'{reg} (Xb4b 0.30", r<0.25")')
        else: ax[1].scatter(o[np.abs(v)<250],h[np.abs(v)<250],s=8,color=col,label=f'{reg} (Xb4d 1.05", r<0.5")')
        rows.append(dict(ous=ous,beam=tb[0],reg=reg,TbHCN_peak=h.max(),TbHCO_peak=o.max(),slope=sl))
for a,t in zip(ax,('Xb4b 0.30"','Xb4d 1.05"')):
    a.plot([0,40],[0,40],'k:',lw=1); a.set_xlabel('Tb HCO+ 3-2 (K)'); a.set_ylabel('Tb HCN 3-2 (K)'); a.legend(); a.set_title(t+', per 4.4 km/s channel, |v|<250')
plt.tight_layout(); plt.savefig('/workspace/out/q9b/tb_scatter.png',dpi=110); pd.DataFrame(rows).to_csv('/workspace/out/q9b/tb_peaks.csv',index=False)
