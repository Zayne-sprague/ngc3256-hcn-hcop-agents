"""Q60 step 2 (CASA host; run: rawdata.py casa q60c_uvspec.py). For every 4-3 split MS (Xb47 ACA x9, 00493 X1df3 ACA x3 [HCN],
00223 X85a ACA x2 [HCO+], 01246 X398b 12m x1 [HCN]) fit, CHANNEL BY CHANNEL and directly in the uv plane (no imaging), the visibilities
(continuum NOT subtracted) with a linear model  V = F_N g_N + F_S g_S (+ F_E g_E)  [g = circular Gaussian of FWHM 0.6" at the N / S positions; g_E = 8" Gaussian at the
mid-point, 'extended' term], in several uv sets (all >=8 klam, and annuli). Output: spectra F_N(ch),F_S(ch),F_E(ch) [Jy, Jy/channel-equivalent flux of the component] + errors, per MS/line/uvset/model.
Everything downstream (baseline, window, curve of growth, N/S ratios) is done on these spectra in numpy (q60d_analyse.py). Output out/q60/vis/<eb>_<line>.npz"""
import numpy as np, glob, os
from casatools import table as tbtool, ms as mstool
C=299792.458; REST={'HCN':354.505477e9,'HCO+':356.734223e9}
RA0,DE0=156.96346,-43.90461
POS={'N':(0.0,2.55),'S':(-0.15,-2.55)}   # arcsec offsets E,N from nominal centre (out/regions.csv)
SETS={'all':(8e3,1e9),'A':(8e3,15e3),'B':(15e3,25e3),'C':(25e3,50e3),'D':(50e3,1e9)}
JOBS=[('20261029T170032-7502',0,'HCN'),('20261029T170032-7502',1,'HCO+'),('20261029T170032-e34b',2,'HCN'),('20261029T170032-8fc3',1,'HCO+'),('20261029T170058-3ae5',0,'HCN')]
def gauss(u,v,fwhm):  # FT of unit-flux circular Gaussian, u,v in wavelengths
    s=fwhm/2.3548/206265.; return np.exp(-2*np.pi**2*s**2*(u**2+v**2))
def comp(u,v,dx,dy,fwhm,sg):
    l=dx/206265.; m=dy/206265.; return gauss(u,v,fwhm)*np.exp(-2j*np.pi*sg*(u*l+v*m))
def fitspec(V,W,B):    # V,W [nch,nr]; B [nr,k] complex basis -> F [nch,k] real, err [nch,k]
    k=B.shape[1]; BH=B.conj()
    # real-valued normal equations for real flux parameters: Re(sum w conj(B_i) B_j), Re(sum w conj(B_i) V)
    M=np.einsum('cr,ri,rj->cij',W,BH,B).real; b=np.einsum('cr,ri,cr->ci',W,BH,V).real
    Mi=np.linalg.inv(M+1e-30*np.eye(k)); F=np.einsum('cij,cj->ci',Mi,b); E=np.sqrt(np.einsum('cii->ci',Mi)); return F,E
for d,spw,ln in JOBS:
  for m in sorted(glob.glob(f'/archive/derived/team-a/{d}/*.split.ms')):
    eb=os.path.basename(m).split('.')[0]; out=f'/workspace/out/q60/vis/{eb}_{ln}.npz'
    if os.path.exists(out): continue
    s=mstool(); s.open(m)
    try: fr=np.array(s.cvelfreqs(spwids=[spw],mode='channel',nchan=-1,start=0,width=1,outframe='LSRK'))
    except Exception as e: fr=None
    s.close()
    t=tbtool(); t.open(m+'/SPECTRAL_WINDOW'); frt=t.getcell('CHAN_FREQ',spw); t.close()
    if fr is None or len(fr)!=len(frt): fr=frt
    t=tbtool(); t.open(m+'/FIELD'); ph=t.getcol('PHASE_DIR')[:,0,:]; t.close()
    ra,de=np.degrees(ph[0,0]),np.degrees(ph[1,0])
    t=tbtool(); t.open(m); cols=t.colnames(); dc='CORRECTED_DATA' if 'CORRECTED_DATA' in cols else 'DATA'
    q=t.query('DATA_DESC_ID==%d'%spw if False else 'DATA_DESC_ID in [%s]'%spw); 
    # DDID != spw in general: map via DATA_DESCRIPTION
    t.close(); t=tbtool(); t.open(m+'/DATA_DESCRIPTION'); sp=t.getcol('SPECTRAL_WINDOW_ID'); t.close(); dd=int(np.where(sp==spw)[0][0])
    t=tbtool(); t.open(m); q=t.query('DATA_DESC_ID==%d'%dd)
    D=q.getcol(dc); F=q.getcol('FLAG'); uvw=q.getcol('UVW'); a1=q.getcol('ANTENNA1'); a2=q.getcol('ANTENNA2'); tm=q.getcol('TIME'); q.close(); t.close()
    ok=a1!=a2; D=D[:,:,ok]; F=F[:,:,ok]; uvw=uvw[:,ok]
    W=(~F).astype(float); Dm=np.where(F,0,D); w2=W.sum(0); V=Dm.sum(0)/np.maximum(w2,1)   # pol-average [nch,nr]
    Wt=(w2>0).astype(float)
    # per-row noise from line-free channels is not used; uniform weights (equal-sensitivity baselines within an EB)
    f0=np.mean(fr); lam=C*1e3/f0; u=uvw[0]/lam; v=uvw[1]/lam; rho=np.hypot(u,v)
    # field offset of nominal positions relative to this MS phase centre
    cd=np.cos(np.radians(DE0)); res={}
    # sign convention test
    best=None
    for sg in (+1,-1):
        dxs={k:((RA0-ra)*3600*cd+p[0],(DE0-de)*3600+p[1]) for k,p in POS.items()}
        B=np.stack([comp(u,v,*dxs['N'],0.6,sg),comp(u,v,*dxs['S'],0.6,sg)],1)
        sel=(rho>8e3)
        # use channel-integrated amplitude of whole spw as quick test: total power coherent with model
        Fq,Eq=fitspec(V[:,sel],Wt[:,sel],B[sel]); score=(Fq/Eq).sum() if False else np.nansum(np.abs(Fq/Eq)**2)
        if best is None or score>best[0]: best=(score,sg)
    sg=best[1]; dxs={k:((RA0-ra)*3600*cd+p[0],(DE0-de)*3600+p[1]) for k,p in POS.items()}
    dxm=((RA0-ra)*3600*cd,(DE0-de)*3600)
    BN=comp(u,v,*dxs['N'],0.6,sg); BS=comp(u,v,*dxs['S'],0.6,sg); BE=comp(u,v,dxm[0],dxm[1],8.0,sg)
    for sname,(lo,hi) in SETS.items():
        sel=(rho>=lo)&(rho<hi)
        if sel.sum()<50: continue
        for mod,Bs in (('2',[BN,BS]),('3',[BN,BS,BE])):
            B=np.stack(Bs,1)[sel]; Fx,Ex=fitspec(V[:,sel],Wt[:,sel],B)
            res[f'F_{sname}_{mod}']=Fx; res[f'E_{sname}_{mod}']=Ex
    kl=rho/1e3
    np.savez(out,freq=fr,rest=REST[ln],sign=sg,nrow=len(u),rmin=kl.min(),rmax=kl.max(),r50=np.median(kl),ph=np.array([ra,de]),
             uvhist=np.histogram(kl,bins=[0,8,15,25,50,100,300,1e4])[0],tmid=tm.mean(),**res)
    print(eb,ln,'rows',len(u),'uv klam %.1f-%.1f med %.1f'%(kl.min(),kl.max(),np.median(kl)),'sign',sg,'nch',len(fr),flush=True)
