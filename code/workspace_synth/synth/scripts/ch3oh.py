# CH3OH: LTE template fits (all catalogue CH3OH lines, Eu<250 K, relative strengths nu*S*mu2*exp(-Eu/Trot), common v and FWHM) to the region spectra of every OUS.
# Free: A (K km/s scale), v, FWHM, linear baseline per window; Trot scanned on a grid (profile chi2).  Also free-amplitude flux of Band-5 Q-branch groups.
# Output: out/ch3oh_fits.csv, out/ch3oh_groups.csv, out/ch3oh_fits.png
import sys,warnings,json; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd
from scipy.optimize import minimize
from astropy.io import fits
from ident import cat
C=299792.458
d=np.load('/workspace/out/region_spectra.npz')
c=cat(300); c=c[c.sp.str.startswith('CH3OH')&(c.logA>-6.5)].copy()
BEAM={}
def beam(ous,spw):
    k=(ous,spw)
    if k not in BEAM:
        import glob
        p=glob.glob(f'/data/member.uid___A001_X1273_{ous}.ngc_3256_sci.{spw}.cube.I.pbcor.fits')[0]
        h=fits.getheader(p); BEAM[k]=(h['BMAJ']*3600,h['BMIN']*3600)  # nominal (median) beam
    return BEAM[k]
def toK(s,fGHz,bm): return 1.222e6*s/(fGHz**2*bm[0]*bm[1])
# windows: (ous, spw) -> list of (flo,fhi) GHz sky containing CH3OH lines for the model (avoid strong non-CH3OH lines)
WIN={'Xb55':{'spw22':[(163.30,164.45)],'spw20':[(164.55,165.60)]},
     'Xb4f':{'spw20':[(251.45,251.95)],'spw22':[(249.10,249.80)],'spw18':[(262.90,264.60)],'spw16':[(264.70,265.10)]},
     'Xb4d':{'spw29':[(251.45,251.95)],'spw31':[(249.10,249.80)],'spw27':[(262.90,264.60)],'spw25':[(264.70,265.10)]},
     'Xb4b':{'spw29':[(251.45,251.95)],'spw31':[(249.10,249.80)],'spw27':[(262.90,264.60)],'spw25':[(264.70,265.10)]}}
def template(fr,T,v,fw,fref,lines=None):
    """K (relative A=1 K km/s per unit weight) LTE optically thin Tb(fr); weights w=nu*Smu2*exp(-Eu/T)/wmax"""
    l=c if lines is None else lines
    w=l.f0*l.sijmu2*np.exp(-l.Eu/T); w=w/w.max()
    fl=l.f0.values*(1-v/C); sg=fw/2.3548/C*fl
    out=np.zeros_like(fr)
    for fi,wi,si in zip(fl,w,sg):
        m=np.abs(fr-fi)<6*si
        if m.any(): out[m]+=wi*np.exp(-0.5*((fr[m]-fi)/si)**2)/(si/fi*C*2.5066)   # K per (K km/s) in velocity units
    return out
def getdata(ous,reg):
    segs=[]
    for spw,ws in WIN[ous].items():
        kf=f'{ous}__{spw}__freq'; kk=f'{ous}__{spw}__{reg}__pk'
        if kf not in d.files or kk not in d.files: continue
        fr=d[kf]; s=d[kk]; bm=beam(ous,spw)
        for lo,hi in ws:
            m=(fr>=lo)&(fr<=hi)&np.isfinite(s)
            if m.sum()<8: continue
            f=fr[m]; y=toK(s[m],f,bm)
            # local rms: robust, from the whole spw with 3-sigma clipping on high-pass residuals
            ya=toK(s[np.isfinite(s)],fr[np.isfinite(s)],bm); r=np.diff(ya)/np.sqrt(2); sg=1.4826*np.median(np.abs(r-np.median(r)))
            segs.append(dict(spw=spw,f=f,y=y,sg=sg))
    return segs
def fit(segs,T,v,fw,fref=None):
    """linear solve for A and a linear baseline per segment at given T,v,fw; returns chi2, A, eA, model"""
    cols=[];Y=[];S=[]
    nseg=len(segs); n=sum(len(s['f']) for s in segs)
    M=np.zeros((n,1+2*nseg)); Yv=np.zeros(n); W=np.zeros(n); i=0
    for k,s in enumerate(segs):
        m=len(s['f']); t=template(s['f'],T,v,fw,None)
        M[i:i+m,0]=t; M[i:i+m,1+2*k]=1; M[i:i+m,2+2*k]=(s['f']-s['f'].mean())/0.1
        Yv[i:i+m]=s['y']; W[i:i+m]=1/s['sg']; i+=m
    Mw=M*W[:,None]; Yw=Yv*W
    p,res,rk,sv=np.linalg.lstsq(Mw,Yw,rcond=None)
    if p[0]<0:   # A>=0: refit without the template
        M[:,0]=0; Mw=M*W[:,None]; p=np.linalg.lstsq(Mw,Yw,rcond=None)[0]; p[0]=0.
    cov=np.linalg.pinv(Mw.T@Mw+1e-12*np.eye(Mw.shape[1])); chi2=np.sum((Mw@p-Yw)**2)
    return chi2,p[0],np.sqrt(cov[0,0]),M@p,n
TGRID=[4,6,8,10,13,17,22,30,45,70,120]
def run(ous,reg,v0,vwin=40,fwlist=(120,150,185,230)):
    segs=getdata(ous,reg)
    if not segs: return None
    best=None; prof={}
    for T in TGRID:
        bT=None
        for v in np.arange(v0-vwin,v0+vwin+1,10):
            for fw in fwlist:
                chi2,A,eA,mod,n=fit(segs,T,v,fw)
                if bT is None or chi2<bT[0]: bT=(chi2,A,eA,v,fw,n)
        prof[T]=bT
        if best is None or bT[0]<best[0][0]: best=(bT,T)
    # null chi2 (A=0)
    chi0=0; 
    for s in segs:
        M=np.vstack([np.ones_like(s['f']),(s['f']-s['f'].mean())/0.1]).T; p=np.linalg.lstsq(M/s['sg'],s['y']/s['sg'],rcond=None)[0]; chi0+=np.sum(((M@p-s['y'])/s['sg'])**2)
    return best,prof,chi0,segs
if __name__=='__main__':
    VLOC={('Xb55','N'):2777,('Xb4f','N'):2779,('Xb4f','S'):2764,('Xb4f','Bridge'):2778,('Xb4d','N'):2790,('Xb4d','S'):2722,('Xb4d','Bridge'):2783,('Xb4b','N'):2817,('Xb4b','S'):2710,('Xb4b','Bridge'):2801,
          ('Xb4f','NEclump'):2840,('Xb4f','Eclump'):2862,('Xb4f','Soutflow'):2774,('Xb4d','NEclump'):2833,('Xb4d','Eclump'):2865,('Xb4d','Soutflow'):2730}
    rows=[]
    for (ous,reg),v0 in VLOC.items():
        r=run(ous,reg,v0)
        if r is None: continue
        (bT,T),prof,chi0,segs=r
        chi2,A,eA,v,fw,n=bT
        prof_d={k:(p[0]-bT[0]) for k,p in prof.items()}
        # A at each T (K km/s of strongest line) and 3sig limit style
        rows.append(dict(ous=ous,region=reg,Tbest=T,A_Kkms=A,eA=eA,snr=A/eA,v=v,fwhm=fw,chi2=chi2,chi2_null=chi0,dchi2=chi0-chi2,npts=n,
              **{f'dchi2_T{k}':round(prof[k][0]-chi2,1) for k in TGRID},**{f'A_T{k}':prof[k][1] for k in TGRID},**{f'eA_T{k}':prof[k][2] for k in TGRID}))
        print(ous,reg,f'Tbest={T} A={A:.3f}+-{eA:.3f} ({A/eA:.1f}s) v={v:.0f} fw={fw} dchi2={chi0-chi2:.0f}',flush=True)
    pd.DataFrame(rows).to_csv('/workspace/out/ch3oh_fits.csv',index=False)
