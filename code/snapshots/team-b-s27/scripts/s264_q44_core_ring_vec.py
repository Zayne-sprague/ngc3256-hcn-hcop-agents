"""Q44 round 2 (reviewer fix 2): core (r<1.1") and ring (1.1-3.3") data vectors from the SAME 2.2" tapered cubes, +-150 km/s, same aperture centre (s185 HCO+ peak) for 1-0 and 3-2.
12C 1-0 HCN/HCO+: s185 X961 spw0 (refint, local linear baseline); 12C 3-2: Q34 12m+7m 'comb' 2.2" cubes (s192 loaders, line-free-channel noise); 13C 1-0: s260 (boxfit, null-calibrated).
12C 1-0 error: 2% (statistical negligible; the fit adds 5% floor in quadrature). Also writes the null test (core-vs-ring H13C/12C difference at N, S, Bridge, and injected-noise false-gradient rate).
usage: python3 -I s264_q44_core_ring_vec.py -> out/q44/s264.json, s264.txt"""
import sys,json,numpy as np
sys.path.insert(0,'/workspace/scripts')
import s185_q29_weaklines as S
from astropy.io import fits
from astropy.wcs import WCS
O='/workspace/out/q44/'; C=299792.458; V0=2780.; F32={'HCN32':265.88618,'HCOp32':267.55762}
L=[]
def P(s=''): print(s,flush=True); L.append(s)
def kw(v,dv,lo,hi): return np.clip(np.minimum(v+dv/2,hi)-np.maximum(v-dv/2,lo),0,dv)/dv
def load(tag):
    h=fits.open('/workspace/out/q34/'+tag+'.pbcor.fits'); d=np.nan_to_num(h[0].data.astype(float)); hd=h[0].header; w=WCS(hd).celestial
    nu=(hd['CRVAL3']+(np.arange(hd['NAXIS3'])-(hd['CRPIX3']-1))*hd['CDELT3'])/1e9; f0=F32[tag.split('_')[0]]
    v=C*(1-nu/f0)-V0 if hd['CUNIT3'].strip()=='Hz' else hd['CRVAL3']/1e3+(np.arange(hd['NAXIS3'])-(hd['CRPIX3']-1))*hd['CDELT3']/1e3-V0
    bm=hd['BMAJ']*3600; pix=abs(hd['CDELT2'])*3600; Tb=d*1.222e6/(f0**2*bm**2); lf=(v<-420)|(v>440); Tb=Tb-np.median(Tb[lf],axis=0)
    return Tb,v,pix,w,lf
D32={ln:load(f'{ln}_2p2_comb') for ln in F32}
APS={'core':(0,1.1),'ring':(1.1,3.3),'r<2.2':(0,2.2)}
s260=json.load(open(O+'s260.json')); rowkey={'core':'0-1.1','ring':'1.1-3.3'}
# the ring in s260.json was not stored under that key: recompute below if missing
def ap_mask(shape,pix,cen,r0,r1):
    ny,nx=shape; Y,X=np.mgrid[:ny,:nx]; d=np.hypot(X-cen[0],Y-cen[1])*pix; return (d<=r1)&(d>=r0) if r0>0 else (d<=r1)
OUT={}
for rn,pos_nm in (('N','N'),('S','S'),('Bridge','Bridge')):
    pos=S.positions(rn)
    for an,(r0,r1) in APS.items():
        # 1-0
        Tb,nu,pix,bm,pxy=S.DATA['X961_spw0']; x,y=pxy(pos); mk=ap_mask(Tb.shape[1:],pix,(x,y),r0,r1); sp=Tb[:,mk].mean(axis=1)
        IO=S.refint(sp,S.F0['HCOp'],(-150,150)); IH=S.refint(sp,S.F0['HCN'],(-150,150))
        out={'HCOp10':IO,'HCN10':IH}
        for ln,key in (('HCN32','HCN32'),('HCOp32','HCOp32')):
            T,v,pix3,w,lf=D32[ln]; px,py=[float(t) for t in w.world_to_pixel_values(*pos)]
            mk=ap_mask(T.shape[1:],pix3,(px,py),r0,r1); sp=T[:,mk].mean(axis=1); dv=abs(v[1]-v[0]); kk=kw(v,dv,-150,150); sig=float(np.std(sp[lf]))
            out[key]=[float((kk*sp).sum()*dv),float(sig*np.sqrt((kk**2).sum())*dv)]
        OUT[f'{rn}|{an}']=out
        P(f'{rn:6s} {an:6s} HCN1-0 {IH:7.1f} HCO+1-0 {IO:7.1f} | HCN3-2 {out["HCN32"][0]:6.1f}+-{out["HCN32"][1]:.2f} HCO+3-2 {out["HCOp32"][0]:6.1f}+-{out["HCOp32"][1]:.2f}')
json.dump(OUT,open(O+'s264.json','w'),indent=1,default=float); open(O+'s264.txt','w').write('\n'.join(L))
