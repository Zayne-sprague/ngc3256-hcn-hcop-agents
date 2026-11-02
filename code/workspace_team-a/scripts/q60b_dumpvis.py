"""Q60 step 2 (CASA host): per MS and line, build line-integrated, continuum-subtracted visibilities
V_L(u,v) = sum_{|v-VC|<HW} (V - linear baseline fit to channels 470<|v-VC|<1100 km/s) * dv   [Jy km/s, LSRK radio velocity], pol XX,YY averaged.
Writes /workspace/out/q60/vis/<eb>_<line>.npz (u,v in m, field, time, ant1, ant2, VL, sig, freq0, field phase centres). Run: rawdata.py casa q60b_dumpvis.py"""
import numpy as np, glob, os
from casatools import ms as mstool, table as tbtool
C=299792.458; REST={'HCN':354.505477e9,'HCO+':356.734223e9}; VC=2767.5; HW=370.; BL=(470.,1100.)
J=[]
for d,sp,ln in [('20261029T170032-7502',(0,'HCN'),None),('20261029T170032-7502',(1,'HCO+'),None),
    ('20261029T170032-e34b',(2,'HCN'),None),('20261029T170032-8fc3',(1,'HCO+'),None),('20261029T170058-3ae5',(0,'HCN'),None)]:
    for m in sorted(glob.glob(f'/archive/derived/team-a/{d}/*.split.ms')): J.append((m,sp[0],sp[1]))
for m,spw,ln in J:
    eb=os.path.basename(m).split('.')[0]; out=f'/workspace/out/q60/vis/{eb}_{ln}.npz'
    if os.path.exists(out): continue
    s=mstool(); s.open(m); 
    try: fr=np.array(s.cvelfreqs(spwids=[spw],mode='channel',nchan=-1,start=0,width=1,outframe='LSRK'))
    except Exception as e:
        print('cvelfreqs failed',e); t=tbtool(); t.open(m+'/SPECTRAL_WINDOW'); fr=t.getcell('CHAN_FREQ',spw); t.close()
    s.selectinit(datadescid=spw)
    r=s.getdata(['data','flag','u','v','field_id','time','antenna1','antenna2','axis_info'],ifraxis=False)
    s.close()
    D=r['data']; F=r['flag']; v=(1-fr/REST[ln])*C; dv=abs(np.median(np.diff(v)))
    D=np.where(F,np.nan,D); D=np.nanmean(D,axis=0) if False else D
    with np.errstate(all='ignore'):
        w=(~F).astype(float); Dm=np.where(F,0,D); V=(Dm.sum(0))/np.maximum(w.sum(0),1)  # [nchan,nrow]
        good=(w.sum(0)>0)
    wb=(np.abs(v-VC)>BL[0])&(np.abs(v-VC)<BL[1]); wl=np.abs(v-VC)<HW
    x=(v-VC)/1000.
    A=np.vstack([np.ones(wb.sum()),x[wb]]).T; P=np.linalg.pinv(A)  # [2,nb]
    coef=P@V[wb]           # [2,nrow]
    base=coef[0][None,:]+coef[1][None,:]*x[:,None]
    res=V-base
    sig_ch=np.sqrt(0.5*(np.var(res[wb].real,axis=0)+np.var(res[wb].imag,axis=0)))  # per-channel rms from the line-free channels
    VL=res[wl].sum(0)*dv; sig=sig_ch*dv*np.sqrt(wl.sum())
    # fully-covered line window only?
    cover=wl.sum()
    t=tbtool(); t.open(m+'/FIELD'); ph=t.getcol('PHASE_DIR')[:,0,:]; t.close()
    ok=good&np.isfinite(sig)&(sig>0)
    np.savez(out,u=r['u'][ok],v=r['v'][ok],field=r['field_id'][ok],time=r['time'][ok],a1=r['antenna1'][ok],a2=r['antenna2'][ok],
        VL=VL[ok],sig=sig[ok],freq0=REST[ln]/(1+VC/C) if False else fr[np.argmin(abs(v-VC))],phase_dir=ph,nwin=cover,nbase=wb.sum(),
        vmin=v.min(),vmax=v.max(),spw=spw)
    print(eb,ln,'rows',ok.sum(),'nwin',cover,'nbase',wb.sum(),'v',v.min(),v.max(),flush=True)
