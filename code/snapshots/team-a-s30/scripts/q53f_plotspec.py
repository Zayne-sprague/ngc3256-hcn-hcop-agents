import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
d=np.load('/workspace/out/q53/spec43.npz'); VC={'N':2785.,'S':2750.}
fig,ax=plt.subplots(2,2,figsize=(11,7))
for i,bm in enumerate((1.0,1.9)):
  for j,reg in enumerate('NS'):
    a=ax[i,j]
    for ln,c in (('HCN','C0'),('HCO+','C3')):
        v=d[f'{ln}_{reg}_{bm}_v']; y=d[f'{ln}_{reg}_{bm}_y']
        # rebin to ~10 km/s
        n=max(1,int(round(10/abs(np.median(np.diff(v)))))); m=len(v)//n*n
        vv=v[:m].reshape(-1,n).mean(1); yy=y[:m].reshape(-1,n).mean(1)
        a.step(vv,yy,c,label=ln,where='mid')
    for s in (-1,1): a.axvline(VC[reg]+s*350,color='k',ls=':'); a.axvline(VC[reg]+s*450,color='gray',ls=':')
    a.set_xlim(1900,3700); a.set_title(f'{reg} {bm}" beam'); a.set_xlabel('v_radio km/s'); a.set_ylabel('T_B (K)'); a.legend()
plt.tight_layout(); plt.savefig('/workspace/out/q53/spec43_overlay.png',dpi=90)
