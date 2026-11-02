# Q2 step a: where is the S-nucleus high-velocity (wing) emission? Xb4d HCN/HCO+ 3-2 (1.05" circular beam) and Xb4b smoothed to 0.5".
# (1) wing moment maps (blue -405..-150, red +150..+405 km/s about VSYS=2780 radio; baseline per pixel median v in -900..-450),
# (2) azimuthal distribution of wing flux in an annulus 1.0-3.5" around S (-> outflow axis PA, empirical),
# (3) noise: random line-free apertures.  Output out/q2/s74a.txt/json, s74a_maps.npz, s74a_maps.png
import sys, json; sys.path.insert(0,'/workspace/scripts')
from s72_lib import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
OUT='/workspace/out/q2'; txt=[]; P=lambda s:(txt.append(s),print(s,flush=True))
WIN={'blue':(-405,-150),'red':(150,405)}
res={}; maps={}
for ous,fw,half in (('Xb4d',1.05,24.0),('Xb4b',0.5,10.0)):
    for mol,spw in (('HCN','spw27'),('HCO+','spw25')):
        c,v,ps,wc,x0,y0,b0=load_cube(ous,spw,F0[mol][3],half=half)
        c=conv(c,ps,b0,(fw,fw,0.0)); c=c-np.median(c[(v>-900)&(v<-450)],axis=0)
        dv=abs(np.median(np.diff(v))); n=c.shape[1]; pxy=pxy_fn(wc,x0,y0); ba=1.1331*fw**2/ps**2
        sx,sy=pxy(*REGIONS['S'][:2]); Y,X=np.mgrid[:n,:n]
        rr=np.hypot(X-sx,Y-sy)*ps; dE=-(X-sx)*ps; dN=(Y-sy)*ps; pa=np.rad2deg(np.arctan2(dE,dN))%360
        for w,(a,b) in WIN.items():
            k=(v>=a)&(v<=b); M=c[k].sum(0)*dv      # Jy/beam km/s (beam-averaged brightness, pbcor)
            # noise: robust sigma in line-free field (rr 4..9" for Xb4b, 4..12 for Xb4d), bright masked
            m0=c[np.abs(v)<350].sum(0)*dv; sg0=np.median(np.abs(m0-np.median(m0)))*1.4826
            br=binary_dilation(m0>4*sg0,iterations=int(1.5/ps)); rcen=np.hypot(X-pxy(*CEN)[0],Y-pxy(*CEN)[1])*ps
            fld=(~br)&(rcen>3.5)&(rcen<(half*0.9)); sig=float(np.std(M[fld]))
            maps[f'{ous}_{mol}_{w}']=M
            ann=(rr>1.0)&(rr<3.5); bins=np.arange(0,361,30); prof=[]
            for i in range(12):
                m=ann&(pa>=bins[i])&(pa<bins[i+1]); prof.append(float(M[m].mean()) if m.any() else np.nan)
            res[f'{ous}_{mol}_{w}']=dict(sig_pix=sig,az_bins=bins[:-1].tolist(),az_mean=prof,ps=ps,beam=fw)
            P(f'{ous} {mol} {w} {a}..{b}: pixel sigma {sig:.4f} Jy/b km/s; annulus 1-3.5" azimuth (PA E of N, 30 deg bins) mean: '+' '.join(f'{p:+.3f}' for p in prof))
        maps[f'{ous}_{mol}_ps']=ps; maps[f'{ous}_{mol}_S']=np.array([sx,sy]); maps[f'{ous}_{mol}_m0']=c[np.abs(v)<350].sum(0)*dv
        del c
json.dump(res,open(OUT+'/s74a.json','w'),indent=1); np.savez(OUT+'/s74a_maps.npz',**maps)
open(OUT+'/s74a.txt','w').write('\n'.join(txt))
fig,ax=plt.subplots(2,4,figsize=(18,9))
for r,ous in enumerate(('Xb4d','Xb4b')):
  for j,(mol,w) in enumerate((('HCN','blue'),('HCO+','blue'),('HCN','red'),('HCO+','red'))):
    M=maps[f'{ous}_{mol}_{w}']; ps=float(maps[f'{ous}_{mol}_ps']); sx,sy=maps[f'{ous}_{mol}_S']; s=res[f'{ous}_{mol}_{w}']['sig_pix']
    a=ax[r,j]; a.imshow(M,origin='lower',cmap='RdBu_r',vmin=-8*s,vmax=8*s,extent=[0,M.shape[1]*ps,0,M.shape[0]*ps])
    m0=maps[f'{ous}_{mol}_m0']; a.contour((np.arange(M.shape[1])+.5)*ps,(np.arange(M.shape[0])+.5)*ps,m0,levels=m0.max()*np.array([.05,.2,.5]),colors='k',linewidths=.6)
    a.plot(sx*ps,sy*ps,'g+',ms=12); a.set_xlim((sx-8/ps)*ps,(sx+8/ps)*ps); a.set_ylim((sy-8/ps)*ps,(sy+8/ps)*ps); a.invert_xaxis(); a.set_title(f'{ous} {mol} {w} wing (rms {s*1e3:.0f} mJy/b km/s)',fontsize=9)
plt.tight_layout(); plt.savefig(OUT+'/s74a_maps.png',dpi=80)
