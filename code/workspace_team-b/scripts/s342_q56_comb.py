"""Q56: R32/10 at N,S with the uv-complete (12m+7m) Band 6 3-2 cubes (s341, circular 2.5" and 3.5") vs the Xb4d image-plane 3-2 (s340), same pixels (s320 POS),
windows (+-100,+-150,+-350 km/s), X961 1-0 integrals via s340 functions. Comb cubes: radio LSRK 20 km/s, baseline = median of channels |v|>380 (cube spans ~-600..+550).
usage: python3 -I s342_q56_comb.py -> out/q56/s342.txt/json"""
import sys, json, numpy as np
sys.path.insert(0,'/workspace/scripts')
from astropy.io import fits
from astropy.wcs import WCS
import s320_q52_chem_ratios as m
import s340_q56_sweep_lib as L
O='/workspace/out/q34/'; C=299792.458
F0={'HCN':265.88618,'HCO+':267.55762}; TAG={'HCN':'HCN32','HCO+':'HCOp32'}
def comb(ln,th,rn):
    tg=f'{TAG[ln]}_{str(th).replace(".","p")}_comb'
    h=fits.open(O+tg+'.pbcor.fits'); d=np.nan_to_num(h[0].data.astype(float)); hd=h[0].header; w=WCS(hd).celestial
    pos,v0=m.POS[rn]; x,y=[float(t) for t in w.world_to_pixel_values(*pos)]
    f0=F0[ln]; bm=hd['BMAJ']*3600
    nu=(hd['CRVAL3']+(np.arange(hd['NAXIS3'])-(hd['CRPIX3']-1))*hd['CDELT3'])/1e9
    v=C*(1-nu/f0)-v0   # frequency axis (LSRK), radio velocity w.r.t. the line rest frequency
    Tb=d*1.222e6/(f0**2*bm**2); sp=Tb[:,int(round(y)),int(round(x))]
    free=abs(v)>380; sp=sp-np.median(sp[free]); sig=np.std(sp[free]); dv=np.full_like(v,abs(np.mean(np.diff(v))))
    return v,sp,sig,dv,bm
out=[]; res={}
def P(t): print(t,flush=True); out.append(t)
for w in (100,150,350):
    P(f'\n== window |v|<{w}: R32/10 uv-complete (12m+7m) vs Xb4d-image 3-2; 1-0 = X961 (12m+ACA), Tb K km/s')
    for th in (2.5,3.5):
        for ln in ('HCN','HCO+'):
            row=f'{th}" {ln}: '; R={}
            for rn in 'NS':
                v,sp,sig,dv,bm=comb(ln,th,rn)
                a=(L.integ(v,sp,dv,-w,w),sig*np.sqrt(2*w*dv[0]))
                b=L.I('10',ln,th,rn,w); x=L.I('32',ln,th,rn,w)
                R[rn]=dict(comb32=a[0],xb4d32=x[0],I10=b[0],r=a[0]/b[0],rx=x[0]/b[0],er=abs(a[0]/b[0])*np.hypot(a[1]/a[0],b[1]/b[0]),fl=a[0]/x[0])
                res[f'{w}|{th}|{ln}|{rn}']=dict(R[rn],beam=bm)
            row+=' | '.join(f'{rn}: R32/10 comb {R[rn]["r"]:.3f}+-{R[rn]["er"]:.3f} (Xb4d {R[rn]["rx"]:.3f}; comb/Xb4d {R[rn]["fl"]:.2f}; I10 {R[rn]["I10"]:.1f}, I32 {R[rn]["comb32"]:.1f})' for rn in 'NS')
            row+=f' | S/N comb {R["S"]["r"]/R["N"]["r"]:.2f} (Xb4d {R["S"]["rx"]/R["N"]["rx"]:.2f})'
            P(row)
        # N/S contrast of HCN/HCO+
    for th in (2.5,3.5):
        g=lambda ln,rn,k: res[f'{w}|{th}|{ln}|{rn}'][k]
        for lab,kk in (('3-2 comb','comb32'),('3-2 Xb4d','xb4d32'),('1-0','I10')):
            c=(g('HCN','N',kk)/g('HCO+','N',kk))/(g('HCN','S',kk)/g('HCO+','S',kk)); P(f' {th}" N/S HCN/HCO+ contrast {lab}: {c:.2f}')
json.dump(res,open('/workspace/out/q56/s342.json','w'),indent=1,default=float); open('/workspace/out/q56/s342.txt','w').write('\n'.join(out))
