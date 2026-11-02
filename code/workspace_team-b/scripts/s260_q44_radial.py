"""Q44 (a): radial test of H13CO+/HCO+ and H13CN/HCN (1-0, +-150 km/s, 2.2" circular beam, X961/X95d 12m+ACA-like cubes) at N, S, Bridge.
Reuses the s185 machinery (imported: it re-runs s185, deterministic, and rewrites out/q29/s185*.json identically). Apertures are centred on the +-150 km/s HCO+ peak (s185 positions),
cumulative r<R and annuli, errors = max(formal*null-calibration, null std) from the s185 null apertures of the SAME shape; +8% cross-project scale on the ratio (as s185).
Also prints Harada 1.7" peak (9.4/205.25 = 0.0458+-0.0049).   usage: python3 -I s260_q44_radial.py -> out/q44/s260.json, s260.txt"""
import sys, json, os, numpy as np
sys.path.insert(0,'/workspace/scripts')
import s185_q29_weaklines as S
O='/workspace/out/q44/'; os.makedirs(O,exist_ok=True); log=[]
def P(s=''): print(s,flush=True); log.append(s)
BOX=(-150,150)
def ap_ann(Tb,pix,cen,r0,r1):
    ny,nx=Tb.shape[1:]; Y,X=np.mgrid[:ny,:nx]; d=np.hypot(X-cen[0],Y-cen[1])*pix; m=(d<=r1)&(d>=r0) if r0>0 else (d<=r1); return Tb[:,m].mean(axis=1)
def spectra_r(pos,r0,r1):
    out={}
    for k,(Tb,nu,pix,bm,pxy) in S.DATA.items():
        x,y=pxy(pos); out[k]=ap_ann(Tb,pix,(x,y),r0,r1)
    return out
APS=[(0,0.3),(0,0.7),(0,1.1),(0,1.5),(0,2.2),(0,3.0),(0.7,1.5),(1.5,2.5),(2.5,3.5)]
res={}
for (r0,r1) in APS:
    nul={}
    for (rr,ang,p,mi) in S.nulls:
        for kk,(a,e) in S.boxfit(spectra_r(p,r0,r1),BOX).items(): nul.setdefault(kk,[]).append((a,e))
    for rn in ('N','S','Bridge'):
        pos=S.positions(rn); sp=spectra_r(pos,r0,r1); R=S.boxfit(sp,BOX)
        IO=S.refint(sp['X961_spw0'],S.F0['HCOp'],BOX); IH=S.refint(sp['X961_spw0'],S.F0['HCN'],BOX)
        row={'HCOp_int':IO,'HCN_int':IH}
        for kk,(a,e) in R.items():
            l=kk[0]
            if l not in ('H13CN','H13COp') or (l=='H13COp' and kk[1]!='comb') or (l=='H13CN' and kk[1]!='X95d_spw0'): continue
            L=np.array(nul[kk]); c=max(1.0,float(np.sqrt(np.mean((L[:,0]/L[:,1])**2)))); e2=max(e*c,float(np.std(L[:,0])))
            ref=IH if l=='H13CN' else IO; rat=a/ref; er=np.hypot(e2/ref,0.08*abs(rat)); row[l]=dict(amp=a,err=e2,ratio=rat,ratio_err=er)
        res[f'{rn}|{r0}-{r1}']=row
        h=row.get('H13COp'); n=row.get('H13CN')
        P(f' {rn:6s} r={r0:.1f}-{r1:.1f}"  HCO+ {IO:6.1f} HCN {IH:6.1f} K km/s | H13CO+/HCO+ {h["ratio"]:.4f}+-{h["ratio_err"]:.4f} ({h["amp"]/h["err"]:.1f} sig) | H13CN/HCN {n["ratio"]:.4f}+-{n["ratio_err"]:.4f} ({n["amp"]/n["err"]:.1f} sig)')
P('Harada+2018 Table 13 (1.7" peak pixel) N: H13CO+/HCO+ 9.4/205.25 = 0.0458+-0.0049; C 3.44/139.3? (see lit table)')
json.dump(res,open(O+'s260.json','w'),indent=1,default=float); open(O+'s260.txt','w').write('\n'.join(log))
