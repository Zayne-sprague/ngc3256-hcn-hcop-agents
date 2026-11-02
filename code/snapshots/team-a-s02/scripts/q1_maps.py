# Q1: HCN/HCO+ 3-2 ratio maps (brightness-temperature ratio) at matched beams, S/N masks, pixel statistics per zone.
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
from q1_lib import *
import pandas as pd, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
CASES=[(1.05,'Xb4d'),(1.05,'Xb4b'),(6.5,'Xb4f'),(6.5,'Xb4d'),(6.5,'Xb4b')]
ZONES={'N':(0,2.55),'S':(-0.15,-2.55),'Bridge':(0.3,0.15),'Wclump':(-1.8,-2.4),'Soutflow':(-0.15,-4.5)}
ZR={1.05:0.8,6.5:2.0}   # zone radius (arcsec)
K=lambda key,fw: Tb_factor(LINES[key][2]/(1+VSYS/C),fw)
rows=[];fig,ax=plt.subplots(1,5,figsize=(26,5.4)); out={}
for i,(fw,ous) in enumerate(CASES):
    kn,kd=f'HCN 3-2@{ous}',f'HCO+ 3-2@{ous}'
    n,dn=[load(k,crop=11.0 if ous=='Xb4b' else None) for k in (kn,kd)]
    sn,sd=smooth(n,fw),smooth(dn,fw)
    pn,pd_=noise_profile(n,fw),noise_profile(dn,fw)
    r=np.hypot(n['dx'],n['dy']); en=np.vectorize(lambda q:sigma_at(pn,q))(r); ed=np.vectorize(lambda q:sigma_at(pd_,q))(r)
    # both beam-matched: ratio of Jy/beam km/s maps; Tb ratio = x (nu_HCO+/nu_HCN)^2
    cT=K(kn,fw)/K(kd,fw)
    for thr in (3,5):
        m=(sn/en>thr)&(sd/ed>thr)&(r<10)
        R=np.where(m,sn/sd*cT,np.nan); eR=np.abs(R)*np.hypot(en/sn,ed/sd)
        if thr==5: out[(fw,ous)]=(R,eR,n['dx'],n['dy'],sd)
        stride=max(1,int(round(fw/2/n['pix'])))   # ~2 samples per beam for the statistics
        for z,(x0,y0) in ZONES.items():
            if fw>2: continue   # 6.5" zones hold <3 independent samples
            zm=m&(np.hypot(n["dx"]-x0,n["dy"]-y0)<ZR[fw]); zs=np.zeros_like(zm); zs[::stride,::stride]=True; zm=zm&zs
            if zm.sum()<3: continue
            v=R[zm]; q=np.percentile(v,[16,50,84])
            rows.append(dict(beam=fw,ous=ous,snr_cut=thr,zone=z,n_samples=int(zm.sum()),median=q[1],p16=q[0],p84=q[2],IQR=np.subtract(*np.percentile(v,[75,25])),
                             sum_ratio=float(np.sum(sn[zm])/np.sum(sd[zm])*cT),min=v.min(),max=v.max()))
    R,eR,dx,dy,sdm=out[(fw,ous)]
    a=ax[i]; im=a.pcolormesh(dx,dy,R,vmin=0,vmax=1.2,cmap='viridis',shading='auto')
    a.contour(dx,dy,sdm,levels=np.nanmax(sdm)*np.array([.05,.1,.2,.4,.8]),colors='w',linewidths=.5)
    a.set_xlim(8,-8); a.set_ylim(-8,8); a.set_aspect('equal'); a.set_title(f'HCN/HCO+ 3-2 (Tb), {fw}" beam, {ous}, S/N>5'); a.plot([0,-.15],[2.55,-2.55],'rx'); plt.colorbar(im,ax=a,fraction=.046)
    np.savez_compressed(f'/workspace/out/q1/ratio_map_{ous}_{fw}.npz',R=R,eR=eR,dx=dx,dy=dy)
plt.tight_layout(); plt.savefig('/workspace/out/q1/ratio_maps.png',dpi=90)
S=pd.DataFrame(rows); S.to_csv('/workspace/out/q1/q1_ratio_map_stats.csv',index=False)
pd.set_option('display.width',250,'display.max_rows',400); print(S[S.snr_cut==5].round(3).to_string())
