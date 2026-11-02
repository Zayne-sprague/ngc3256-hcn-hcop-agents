"""Q34 s192: measure uv-matched Band 6 3-2 cubes (s191: 2017.1.00379.S Xb4d 12m EB Xe07 [+ Xb4f 7m EBs], taper + restoringbeam 1.5"/2.2" circular, pbcor, 20 km/s)
at N, S, Bridge with p0034 conventions (aperture r = beam FWHM centred on HCO+ 3-2 moment peak within 0.4" of nominal, +-350 and +-150 km/s about 2780, RJ T_b, K km/s,
noise from line-free channels of the aperture-mean spectrum).  Compares with the Band 3 1-0 values of s183 (X145 1.5", X961 2.2") and with the image-plane Xb4d convolved to the same beam (s183 xb4d_32).
Variants: comb (12m+7m, uv >~7.7 kl) and c12 (12m only).  Also a 6" aperture (disk) flux fraction comb vs Xb4d.  usage: python3 -I s192_q34_measure.py"""
import sys, json, numpy as np
sys.path.insert(0,'/workspace/scripts')
from astropy.io import fits
from astropy.wcs import WCS
from regions import REGIONS
O='/workspace/out/q34/'; C=299792.458; V0=2780.
F0={'HCN32':265.88618,'HCOp32':267.55762}
REG={k:REGIONS[k][:2] for k in ('N','S','Bridge')}
r183=json.load(open('/workspace/out/q29/s183.json'))
def load(tag):
    h=fits.open(O+tag+'.pbcor.fits'); d=h[0].data.astype(float); hd=h[0].header; w=WCS(hd).celestial
    nu=(hd['CRVAL3']+(np.arange(hd['NAXIS3'])-(hd['CRPIX3']-1))*hd['CDELT3'])/1e9
    f0=F0[tag.split('_')[0]]; v=C*(1-nu/f0)-V0 if hd['CUNIT3'].strip()=='Hz' else None
    if v is None: v=hd['CRVAL3']/1e3+(np.arange(hd['NAXIS3'])-(hd['CRPIX3']-1))*hd['CDELT3']/1e3-V0   # radio velocity axis (m/s)
    bm=hd['BMAJ']*3600; pix=abs(hd['CDELT2'])*3600
    Tb=d*1.222e6/(f0**2*bm**2)           # K (circular beam)
    return Tb,v,pix,w
def kw(v,dv,lo,hi): return np.clip(np.minimum(v+dv/2,hi)-np.maximum(v-dv/2,lo),0,dv)/dv
res={}; log=[]
def P(s=''): print(s,flush=True); log.append(s)
for bn,th in (('1p5',1.5),('2p2',2.2)):
  for vn in ('comb','c12'):
    D={}
    for ln in ('HCN32','HCOp32'):
        Tb,v,pix,w=load(f'{ln}_{bn}_{vn}'); Tb=np.nan_to_num(Tb); lf=(v<-420)|(v>440); Tb=Tb-np.median(Tb[lf],axis=0); D[ln]=(Tb,v,pix,w,lf)
    Tb,v,pix,w,lf=D['HCOp32']; ny,nx=Tb.shape[1:]; Y,X=np.mgrid[:ny,:nx]; dv=abs(v[1]-v[0])
    k=kw(v,dv,-350,350); m0=np.tensordot(k,Tb,axes=(0,0))*dv
    for rn,pos in REG.items():
        px,py=[float(t) for t in w.world_to_pixel_values(*pos)]
        mm=np.where(np.hypot(X-px,Y-py)*pix<0.4,m0,-np.inf); j=np.argmax(mm); mp=(float(X.flat[j]),float(Y.flat[j]))
        out={}
        for rad,rl in ((th,'ap'),(6.0,'disk6')):
            apm=np.hypot(X-mp[0],Y-mp[1])*pix<=rad
            for ln in ('HCN32','HCOp32'):
                T,v_,_,_,lf_=D[ln]; sp=T[:,apm].mean(axis=1); sig=float(np.std(sp[lf_]))
                for wn,(lo,hi) in (('w350',(-350,350)),('h150',(-150,150))):
                    kk=kw(v_,dv,lo,hi); out[f'{ln}|{rl}|{wn}']=[float((kk*sp).sum()*dv),float(sig*np.sqrt((kk**2).sum())*dv)]
        res[f'{bn}|{vn}|{rn}']=out
    P(f'beam {th}" {vn}: peaks done')
# comparison
def ctr(a,b,key,i):  # N/S contrast of HCN/HCO+
    return (a[f'HCN32|ap|{key}'][0]/a[f'HCOp32|ap|{key}'][0])/(b[f'HCN32|ap|{key}'][0]/b[f'HCOp32|ap|{key}'][0])
summ={}
for bn,th,k10 in (('1p5',1.5,'2015.1.00412 X145|1.5'),('2p2',2.2,'2015.1.00993 X961|2.2')):
  thk='%.2f'%th
  P(f'\n################ beam {th}"  (1-0: {k10}; Xb4d image-plane convolved 3-2 for reference)')
  for vn in ('comb','c12'):
    for wn,w10,w32 in (('w350','w350_ap','win'),('h150','h150_ap','h150')):
        d={};
        for rn in ('N','S','Bridge'):
            m=res[f'{bn}|{vn}|{rn}']; a=r183[f'{k10}|{rn}']; x=r183['xb4d_32'][thk][rn]
            d[rn]=dict(H32=m[f'HCN32|ap|{wn}'],O32=m[f'HCOp32|ap|{wn}'],H10=a['HCN'][w10],O10=a['HCOp'][w10],Hd=x[w32]['HCN32'],Od=x[w32]['HCOp32'])
        def con(a,b,h,o): return (d[a][h][0]/d[a][o][0])/(d[b][h][0]/d[b][o][0])
        def cerr(a,b,h,o): return con(a,b,h,o)*np.sqrt(sum((d[q][t][1]/d[q][t][0])**2 for q in (a,b) for t in (h,o)))
        for a,b in (('N','S'),('N','Bridge'),('S','Bridge')):
            c10=con(a,b,'H10','O10'); e10=cerr(a,b,'H10','O10'); c32=con(a,b,'H32','O32'); e32=cerr(a,b,'H32','O32'); cd=con(a,b,'Hd','Od')
            Dv=c32/c10; eD=Dv*np.hypot(e10/c10,e32/c32)
            P(f' {vn} {wn} {a}/{b}: 1-0 {c10:.3f}+-{e10:.3f} | 3-2 matched-uv {c32:.3f}+-{e32:.3f} (Xb4d image {cd:.3f}) | D={Dv:.3f}+-{eD:.3f} ({(Dv-1)/eD:.1f} sig)')
            summ[f'{bn}|{vn}|{wn}|{a}/{b}']=dict(c10=c10,e10=e10,c32=c32,e32=e32,cd=cd,D=Dv,eD=eD)
        for ln,(h1,h3,hd) in (('HCN',('H10','H32','Hd')),('HCO+',('O10','O32','Od'))):
            row=f' {vn} {wn} r31 {ln}: '
            for rn in ('N','S','Bridge'):
                q=d[rn]; r=q[h3][0]/q[h1][0]; e=r*np.hypot(q[h3][1]/q[h3][0],q[h1][1]/q[h1][0]); rd=q[hd][0]/q[h1][0]
                f=q[h3][0]/q[hd][0]; summ[f'{bn}|{vn}|{wn}|r31|{ln}|{rn}']=dict(r31=r,e=e,r31_xb4d=rd,frac_vs_xb4d=f); row+=f'{rn} {r:.3f}+-{e:.3f} (Xb4d {rd:.3f}; flux/Xb4d {f:.3f}) | '
            nsr=summ[f'{bn}|{vn}|{wn}|r31|{ln}|N']['r31']/summ[f'{bn}|{vn}|{wn}|r31|{ln}|S']['r31']; nsd=summ[f'{bn}|{vn}|{wn}|r31|{ln}|N']['r31_xb4d']/summ[f'{bn}|{vn}|{wn}|r31|{ln}|S']['r31_xb4d']
            en=nsr*np.hypot(summ[f'{bn}|{vn}|{wn}|r31|{ln}|N']['e']/summ[f'{bn}|{vn}|{wn}|r31|{ln}|N']['r31'],summ[f'{bn}|{vn}|{wn}|r31|{ln}|S']['e']/summ[f'{bn}|{vn}|{wn}|r31|{ln}|S']['r31'])
            summ[f'{bn}|{vn}|{wn}|r31NS|{ln}']=dict(NS=nsr,e=en,NS_xb4d=nsd); P(row+f'N/S {nsr:.3f}+-{en:.3f} (Xb4d {nsd:.3f})')
  # disk6 flux fractions vs c12
  for rn in ('N','S','Bridge'):
    for ln in ('HCN32','HCOp32'):
        a=res[f'{bn}|comb|{rn}'][f'{ln}|disk6|w350'][0]; b=res[f'{bn}|c12|{rn}'][f'{ln}|disk6|w350'][0]
        P(f' r=6" aperture {rn} {ln}: comb {a:.0f} K km/s, c12 {b:.0f} -> c12/comb {b/a:.3f}; r=beam aperture c12/comb {res[f"{bn}|c12|{rn}"][f"{ln}|ap|w350"][0]/res[f"{bn}|comb|{rn}"][f"{ln}|ap|w350"][0]:.3f}')
        summ[f'{bn}|disk6|{rn}|{ln}']=dict(comb=a,c12=b)
json.dump(dict(res=res,summ=summ),open(O+'s192.json','w'),indent=1,default=float); open(O+'s192.txt','w').write('\n'.join(log))
