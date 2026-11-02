"""Q38: robustness of the disk (E1/E2/W clump) HCN/HCO+ decline with J (0.52/0.34/0.14 for 2-1/3-2/4-3 at 10", Q3/Q13) against PB, zero level, nuclear leakage, short spacings.
A: 10" position values (s30 maps) with additive zero-level offsets (fractions of the far-field mean, same sign both lines) and 2.5" apertures; PB model effect on the ratio.
B: decontaminated clump totals (s58d NNLS, free-sign constant) -> R(J) with sign-flip + variant errors; HCN 4-3 as a limit.
C: uv-complete matched 2.2" clump ratios: 1-0 (X961 12m), 3-2 (Q34 12m+7m 'comb' and 12m 'c12'); r(3-2/1-0).
D: 1-0 at 10" (X961, same recipe as s30 positions) to add J=1-0 to the 10" ladder.
usage: python3 -I s240_q38_disk_clumps.py -> out/q38/s240.json/.txt"""
import sys,json,numpy as np
sys.path.insert(0,'/workspace/scripts')
from s180_q29_load import *
from regions import REGIONS
from common import files
from astropy.io import fits
O='/workspace/out/'; log=[]; res={}
def P(s=''): print(s,flush=True); log.append(s)
POS={k:REGIONS[k][:2] for k in ('E1','E2','W','N')}
KF={2:(178.37506/177.26111)**2,3:(267.55763/265.88618)**2,4:(356.73424/354.50547)**2}
FW=2*np.sqrt(2*np.log(2))
# ---------------- A: 10" maps
P('## A: 10" beam position values with zero-level offsets (Jy/beam km/s maps from s30; offset d = f*far-field mean of each line, f in -1,-0.5,0,+0.5,+1; value-at-position bilinear, 0.7" pix interpolation)')
from scipy.ndimage import map_coordinates
from astropy.wcs import WCS
FAR={'Xb55':(30,45),'Xb4f':(26,40),'Xb4d':(17,19),'Xb47':(17,24)}
import s30_matchedJ_ratios as S
CEN=S.CEN
TAG={'Xb55':('Xb55_J2_10as',2),'Xb4f':('Xb4f_J3_10as',3),'Xb4d':('Xb4d_J3_10as',3),'Xb47':('Xb47_J4_10as',4)}
SPWH={'Xb55':'spw16','Xb4f':'spw18','Xb4d':'spw27','Xb47':'spw16'}
A={}
for ous,(tag,J) in TAG.items():
    d=np.load(O+f's30_maps_{tag}.npz'); H,Oo,eH,eO=d['HCN'],d['HCOp'],d['eHCN'],d['eHCOp']; rr=d['rr']; ps=float(d['ps']); cx,cy=float(d['cx']),float(d['cy'])
    lo,hi=FAR[ous]; mk=(rr>=lo)&(rr<hi)&(eH>0); zH,zO=H[mk].mean(),Oo[mk].mean()
    fwpb=1.13*299792458/({'Xb55':177.3,'Xb4f':265.9,'Xb4d':265.9,'Xb47':354.5}[ous]*1e9)/12*206265
    for pn,(ra,dec) in POS.items():
        # pixel of position in map frame: use wcs of same cube geometry
        fn=files(f'*{ous}.ngc_3256_sci.{SPWH[ous]}.cube.I.pbcor.fits')[0]
        hd=fits.getheader(fn); w=WCS(hd).celestial
        # map origin offset: s30 maps are cropped about CEN: pixel (cx,cy)=CEN
        pcx,pcy=w.world_to_pixel_values(*CEN); px,py=w.world_to_pixel_values(ra,dec); x=cx+(px-pcx); y=cy+(py-pcy)
        rpos=np.hypot(px-pcx,py-pcy)*ps
        row={}
        for f in (-1,-0.5,0,0.5,1):
            h=map_coordinates(H,[[y],[x]],order=1)[0]-f*zH; o=map_coordinates(Oo,[[y],[x]],order=1)[0]-f*zO
            row[f]=(h,o,h/o*KF[J] if o>0 and h>0 else float('nan'))
        eh=map_coordinates(eH,[[y],[x]],order=1)[0]; eo=map_coordinates(eO,[[y],[x]],order=1)[0]
        pbg=np.exp(-4*np.log(2)*(rpos/fwpb)**2)
        # PB-error effect on the ratio: HCN and HCO+ pbcor use PB(nu): ratio of PB(HCN)/PB(HCO+) for a PB FWHM in [0.5,2]x nominal
        nuH,nuO=(F0H:={2:177.26111,3:265.88618,4:354.50547}[J]),{2:178.37506,3:267.55763,4:356.73424}[J]
        dpb=[float(np.exp(-4*np.log(2)*(rpos/(fw))**2*((nuO/nuH)**2-1)*(-1)))-1 for fw in (fwpb*0.7,fwpb,fwpb*1.4)]
        A[(ous,pn)]=dict(r=rpos,PB12=pbg,eh=eh,eo=eo,row=row,dpb=dpb,J=J)
        P(f' {ous} J={J} {pn:3s} r={rpos:4.1f}" PB(12m model)={pbg:.2f} HCN {row[0][0]:.2f}+-{eh:.2f} HCO+ {row[0][1]:.2f}+-{eo:.2f} | R_K(offset f=-1..+1): '+' '.join(f'{row[f][2]:.3f}' for f in (-1,-0.5,0,0.5,1))+f' | far-field mean HCN {zH:.2f} HCO+ {zO:.2f} | PB(nu)-ratio bias at FWHM x0.7/1/1.4: '+'/'.join(f'{100*x:.1f}%' for x in dpb))
res['A']={f'{k[0]}|{k[1]}':dict(r=v['r'],PB12=v['PB12'],R_K={str(f):v['row'][f][2] for f in v['row']},HCN=v['row'][0][0],HCOp=v['row'][0][1],eHCN=v['eh'],eHCOp=v['eo']) for k,v in A.items()}
# ---------------- B: clump totals
P('\n## B: de-contaminated clump totals (s58d NNLS with free constant; error = sign-flip stat (+) variant spread (sys), Jy km/s): R(J)=HCN/HCO+ (K scale)')
cf=json.load(open(O+'q12b/clumpfit.json')); txt=open(O+'q12b/clumpfit.txt').read().splitlines()
import re
tot={}
for l in txt:
    m=re.match(r'(HCO\+|HCN)\s+(Xb\w+) beam .*clump total flux \[Jy km/s\] E1 ([\d.]+)\+-([\d.]+)\(sys ([\d.]+)\)\s+E2 ([\d.]+)\+-([\d.]+)\(sys ([\d.]+)\)\s+W ([\d.]+)\+-([\d.]+)\(sys ([\d.]+)\)',l)
    if m:
        g=m.groups(); tot[(g[0],g[1])]={'E1':tuple(map(float,g[2:5])),'E2':tuple(map(float,g[5:8])),'W':tuple(map(float,g[8:11]))}
JX={'Xb55':2,'Xb4f':3,'Xb4d':3,'Xb47':4,'Xb4b':3}
B={}
for ous in ('Xb55','Xb4f','Xb4d','Xb47','Xb4b'):
    for c in ('E1','E2','W'):
        h=tot[('HCN',ous)][c]; o=tot[('HCO+',ous)][c]; e_h=np.hypot(h[1],h[2]); e_o=np.hypot(o[1],o[2])
        if h[0]>2*e_h:
            r=h[0]/o[0]*KF[JX[ous]]; er=r*np.hypot(e_h/h[0],e_o/o[0]); s=f'{r:.3f}+-{er:.3f}'; B[(ous,c)]=(r,er,False)
        else:
            lim=(h[0]+2*e_h)/o[0]*KF[JX[ous]]; s=f'<{lim:.3f} (2 sigma; HCN {h[0]:.2f}+-{e_h:.2f})'; B[(ous,c)]=(lim,np.nan,True)
        P(f' {ous} J={JX[ous]} {c}: HCN {h[0]:.2f}+-{e_h:.2f} HCO+ {o[0]:.2f}+-{e_o:.2f} Jy km/s -> R_K {s}')
res['B']={f'{k[0]}|{k[1]}':v for k,v in B.items()}
for c in ('E1','E2','W'):
    r21,e21,_=B[('Xb55',c)]
    for a in ('Xb4f','Xb4d'):
        r32,e32,lim=B[(a,c)]
        if not lim: P(f' {c}: R32({a})/R21 = {r32/r21:.2f}+-{r32/r21*np.hypot(e32/r32,e21/r21):.2f}')
        else: P(f' {c}: R32({a})/R21 < {r32/r21:.2f}')
    r43,e43,l43=B[('Xb47',c)]; P(f' {c}: R43/R21 {"<" if l43 else "="} {r43/r21:.2f}')
# ---------------- C: matched 2.2" uv-complete
P('\n## C: matched 2.2" circular beam, clump aperture r=2.2" at the HCO+ peak within 1" of nominal; K km/s (RJ), +-350 km/s; 1-0 = X961 12m, 3-2 = Q34 comb (12m+7m) and c12 (12m)')
def kw(v,dv,lo,hi): return np.clip(np.minimum(v+dv/2,hi)-np.maximum(v-dv/2,lo),0,dv)/dv
def vel(nu,f0): return C*(1-nu/f0)-V0
C_={}
X961='/archive/files/2015.1.00993.S/member.uid___A001_X2fb_X961.ari_l.ngc3256_sci.spw0_87827MHz.12m.cube.I.pbcor.fits'
def meas(Tb,v_h,v_o,pix,pxy,label,fwhm):
    ny,nx=Tb[0].shape[1:] if isinstance(Tb,tuple) else Tb.shape[1:]
out={}
def proc(label,TH,TO,vH,vO,pix,pxy,rad,free_h,free_o):
    dvh=abs(vH[1]-vH[0]); dvo=abs(vO[1]-vO[0]); ny,nx=TO.shape[1:]; Y,X=np.mgrid[:ny,:nx]
    k=kw(vO,dvo,-350,350); m0=np.tensordot(k,TO,axes=(0,0))*dvo
    for pn in ('E1','E2','W','N'):
        px,py=pxy(POS[pn]); mm=np.where(np.hypot(X-px,Y-py)*pix<1.0,m0,-np.inf); j=np.argmax(mm); mp=(float(X.flat[j]),float(Y.flat[j]))
        ap=np.hypot(X-mp[0],Y-mp[1])*pix<=rad; r={}
        for ln,T,v,fr,dv in (('HCN',TH,vH,free_h,dvh),('HCOp',TO,vO,free_o,dvo)):
            sp=T[:,ap].mean(axis=1); sig=float(np.std(sp[fr])); kk=kw(v,dv,-350,350)
            r[ln]=((kk*sp).sum()*dv,sig*np.sqrt((kk**2).sum())*dv)
        out[f'{label}|{pn}']=r
        a,b=r['HCN'],r['HCOp']; P(f' {label:12s} {pn:3s} HCN {a[0]:6.2f}+-{a[1]:.2f} HCO+ {b[0]:6.2f}+-{b[1]:.2f} K km/s  R={a[0]/b[0]:.3f}+-{a[0]/b[0]*np.hypot(a[1]/a[0],b[1]/b[0]):.3f}  (offset from nominal {np.hypot(mp[0]-px,mp[1]-py)*pix:.2f}")')
# 1-0
for th,rad in ((2.2,2.2),(10.0,5.0)):
    cut,nu,pix,bm,pxy=load(X961,half=34)
    Tb=to_tb_conv(cut,nu,pix,bm,th,np.mean(nu)); vH=vel(nu,88.6318); vO=vel(nu,89.1885); fr=(np.abs(vH)>600)&(np.abs(vO)>600)
    Tb=Tb-np.median(Tb[fr],axis=0)
    # use the same cube for both lines
    proc(f'1-0 X961 {th}"',Tb,Tb,vH,vO,pix,pxy,rad,fr,fr)
# 3-2 comb / c12 at 2.2"
def l32(tag):
    h=fits.open(O+f'q34/{tag}.pbcor.fits'); d=np.nan_to_num(h[0].data.astype(float)); hd=h[0].header; w=WCS(hd).celestial
    nu=(hd['CRVAL3']+(np.arange(hd['NAXIS3'])-(hd['CRPIX3']-1))*hd['CDELT3'])/1e9; f0=265.88618 if 'HCN' in tag else 267.55762
    v=vel(nu,f0) if hd['CUNIT3'].strip()=='Hz' else (hd['CRVAL3']/1e3+(np.arange(hd['NAXIS3'])-(hd['CRPIX3']-1))*hd['CDELT3']/1e3-V0)
    bm=hd['BMAJ']*3600; pix=abs(hd['CDELT2'])*3600; T=d*1.222e6/(f0**2*bm**2); fr=(v<-420)|(v>440); T=T-np.median(T[fr],axis=0)
    return T,v,pix,w,fr
for vn in ('comb','c12'):
    TH,vH,pix,w,fh=l32(f'HCN32_2p2_{vn}'); TO,vO,_,_,fo=l32(f'HCOp32_2p2_{vn}')
    pxy=lambda pos:tuple(float(t) for t in w.world_to_pixel_values(*pos))
    proc(f'3-2 {vn} 2.2"',TH,TO,vH,vO,pix,pxy,2.2,fh,fo)
res['C']={k:v for k,v in out.items()}
P('\n## C2: 3-2/1-0 ratio of HCN/HCO+ (double ratio) at matched 2.2" per clump')
for pn in ('E1','E2','W','N'):
    r10=out[f'1-0 X961 2.2"|{pn}']; r10v=r10['HCN'][0]/r10['HCOp'][0]; e10=r10v*np.hypot(r10['HCN'][1]/r10['HCN'][0],r10['HCOp'][1]/r10['HCOp'][0])
    for vn in ('comb','c12'):
        r=out[f'3-2 {vn} 2.2"|{pn}']; rv=r['HCN'][0]/r['HCOp'][0]; e=rv*np.hypot(r['HCN'][1]/r['HCN'][0],r['HCOp'][1]/r['HCOp'][0])
        P(f' {pn} {vn}: R10={r10v:.3f}+-{e10:.3f} R32={rv:.3f}+-{e:.3f}  R32/R10={rv/r10v:.3f}+-{rv/r10v*np.hypot(e/rv,e10/r10v):.3f}')

P('\n## C3: r31 = T_b(3-2)/T_b(1-0) per line at matched 2.2" (comb) and excitation contrast clump/N (double ratio of r31)')
def v(k,l): return out[k][l]
for pn in ('E1','E2','W','N'):
    row=f' {pn}:'
    for ln in ('HCN','HCOp'):
        a=v(f'3-2 comb 2.2"|{pn}',ln); b=v(f'1-0 X961 2.2"|{pn}',ln); r=a[0]/b[0]; e=r*np.hypot(a[1]/a[0],b[1]/b[0])
        n_a=v('3-2 comb 2.2"|N',ln); n_b=v('1-0 X961 2.2"|N',ln); rn_=n_a[0]/n_b[0]
        row+=f' r31({ln}) {r:.3f}+-{e:.3f} (clump/N {r/rn_:.2f}+-{e/rn_:.2f});'
        res.setdefault('C3',{})[f'{pn}|{ln}']=[r,e,r/rn_]
    P(row)
P('\n## D: PB mask: positions with PB>=cut (12-m Gaussian model | 7-m-like FWHM x1.7); R_K at the position')
for ous in ('Xb55','Xb4f','Xb4d','Xb47'):
    for pn in ('E1','E2','W'):
        a=A[(ous,pn)]; fw12=a['r']/np.sqrt(-np.log(a['PB12'])/(4*np.log(2))); pb7=np.exp(-4*np.log(2)*(a['r']/(1.7*fw12))**2)
        ok=[f'{c}:{"Y" if a["PB12"]>=c else "n"}|{"Y" if pb7>=c else "n"}' for c in (0.3,0.5,0.7)]
        P(f' {ous} {pn} PB12={a["PB12"]:.2f} PB7={pb7:.2f} R_K={a["row"][0][2]:.3f}  '+'  '.join(ok))
json.dump(res,open(O+'q38/s240.json','w'),indent=1,default=float); open(O+'q38/s240.txt','w').write('\n'.join(log))
