"""Q60 step 2b (CASA host; rawdata.py casa q60f_mosaic.py): redo q60c for MULTI-FIELD MS (2018.1.00493 ACA X1df3 is a 3-pointing mosaic;
q60c used one phase centre for all fields -> decorrelation). Per field: own phase centre, own Gaussian primary-beam gain (7m FWHM=1.13 lam/D)
at the N/S/extended positions (V_row = sum_k F_k PB_f(k) g_k(u,v;phase centre of f)). Fluxes F_k are sky fluxes (PB-corrected). Same outputs/format as q60c.
Also reruns the 12m EB (single field) with the PB factor (12m FWHM 1.13 lam/12m)."""
import numpy as np, os
from casatools import table as tbtool, ms as mstool
C=299792.458; REST={'HCN':354.505477e9,'HCO+':356.734223e9}
RA0,DE0=156.96346,-43.90461
POS={'N':(0.0,2.55),'S':(-0.15,-2.55)}
SETS={'all':(8e3,1e9),'A':(8e3,15e3),'B':(15e3,25e3),'C':(25e3,50e3),'D':(50e3,1e9)}
JOBS=[('20261029T170032-e34b',2,'HCN',7.0),('20261029T170058-3ae5',0,'HCN',12.0)]
import glob
def gauss(u,v,fwhm):
    s=fwhm/2.3548/206265.; return np.exp(-2*np.pi**2*s**2*(u**2+v**2))
def comp(u,v,dx,dy,fwhm,sg,pb):
    return pb*gauss(u,v,fwhm)*np.exp(-2j*np.pi*sg*(u*dx/206265.+v*dy/206265.))
def fitspec(V,W,B):
    k=B.shape[1]; BH=B.conj()
    M=np.einsum('cr,ri,rj->cij',W,BH,B).real; b=np.einsum('cr,ri,cr->ci',W,BH,V).real
    Mi=np.linalg.inv(M+1e-30*np.eye(k)); return np.einsum('cij,cj->ci',Mi,b),np.sqrt(np.einsum('cii->ci',Mi))
for d,spw,ln,D in JOBS:
  for m in sorted(glob.glob(f'/archive/derived/team-a/{d}/*.split.ms')):
    eb=os.path.basename(m).split('.')[0]; out=f'/workspace/out/q60/vis/{eb}_{ln}.npz'
    s=mstool(); s.open(m)
    try: fr=np.array(s.cvelfreqs(spwids=[spw],mode='channel',nchan=-1,start=0,width=1,outframe='LSRK'))
    except Exception: fr=None
    s.close()
    t=tbtool(); t.open(m+'/SPECTRAL_WINDOW'); frt=t.getcell('CHAN_FREQ',spw); t.close()
    if fr is None or len(fr)!=len(frt): fr=frt
    t=tbtool(); t.open(m+'/FIELD'); ph=t.getcol('PHASE_DIR')[:,0,:]; t.close()
    t=tbtool(); t.open(m+'/DATA_DESCRIPTION'); sp=t.getcol('SPECTRAL_WINDOW_ID'); t.close(); dd=int(np.where(sp==spw)[0][0])
    t=tbtool(); t.open(m); cols=t.colnames(); dc='CORRECTED_DATA' if 'CORRECTED_DATA' in cols else 'DATA'
    q=t.query('DATA_DESC_ID==%d'%dd)
    Dd=q.getcol(dc); F=q.getcol('FLAG'); uvw=q.getcol('UVW'); a1=q.getcol('ANTENNA1'); a2=q.getcol('ANTENNA2'); fid=q.getcol('FIELD_ID'); tm=q.getcol('TIME'); q.close(); t.close()
    ok=a1!=a2; Dd=Dd[:,:,ok]; F=F[:,:,ok]; uvw=uvw[:,ok]; fid=fid[ok]
    W=(~F).astype(float); Dm=np.where(F,0,Dd); w2=W.sum(0); V=Dm.sum(0)/np.maximum(w2,1); Wt=(w2>0).astype(float)
    f0=np.mean(fr); lam=C*1e3/f0; u=uvw[0]/lam; v=uvw[1]/lam; rho=np.hypot(u,v)
    pbfw=1.13*lam/D*206265.; cd=np.cos(np.radians(DE0))
    def basis(sg):
        cols=[]
        for key in ('N','S','E'):
            col=np.zeros(len(u),complex)
            for f in np.unique(fid):
                ra,de=np.degrees(ph[0,f]),np.degrees(ph[1,f]); r=fid==f
                if key=='E': dx=(RA0-ra)*3600*cd; dy=(DE0-de)*3600; fw=8.0
                else: dx=(RA0-ra)*3600*cd+POS[key][0]; dy=(DE0-de)*3600+POS[key][1]; fw=0.6
                # sky position relative to pointing centre (for PB): same as (dx,dy) up to the sign convention of the offsets
                pb=np.exp(-4*np.log(2)*(dx**2+dy**2)/pbfw**2)
                col[r]=comp(u[r],v[r],dx,dy,fw,sg,pb)
            cols.append(col)
        return cols
    best=None
    for sg in (+1,-1):
        B=np.stack(basis(sg)[:2],1); sel=rho>8e3
        Fq,Eq=fitspec(V[:,sel],Wt[:,sel],B[sel]); sc=np.nansum(np.abs(Fq/Eq)**2)
        if best is None or sc>best[0]: best=(sc,sg)
    sg=best[1]; BN,BS,BE=basis(sg); res={}
    for sname,(lo,hi) in SETS.items():
        sel=(rho>=lo)&(rho<hi)
        if sel.sum()<50: continue
        for mod,Bs in (('2',[BN,BS]),('3',[BN,BS,BE])):
            B=np.stack(Bs,1)[sel]; Fx,Ex=fitspec(V[:,sel],Wt[:,sel],B); res[f'F_{sname}_{mod}']=Fx; res[f'E_{sname}_{mod}']=Ex
    kl=rho/1e3
    np.savez(out,freq=fr,rest=REST[ln],sign=sg,nrow=len(u),rmin=kl.min(),rmax=kl.max(),r50=np.median(kl),ph=ph[:,0],
             uvhist=np.histogram(kl,bins=[0,8,15,25,50,100,300,1e4])[0],tmid=tm.mean(),**res)
    print(eb,ln,'fields',np.unique(fid),'rows',len(u),'sign',sg,flush=True)
