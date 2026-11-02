# Q14b: CH3OH J_-1-J_0 E Q-branch blend (J=1-3, 165.05-165.10 GHz rest, Xb55 spw22 low-frequency edge) with an explicit edge-ramp baseline,
# versus the 254.015 GHz 2(0)-1(-1)E line (Xb4f smoothed to the Xb55 beam): thin-LTE Trot with baseline systematic; column; CH3OH/HCO+ 2-1.
# Re-implements (independently of team-b p0020) the test that the negative spw22 edge ramp scales with the continuum, then three baselines:
#   L linear, Q quadratic, R linear + exp(-ch/tau) ramp (tau grid), all fitted outside line masks; integral over +-200 km/s of v0=2777.
# Output out/q14/ramp_continuum.csv, blend_integrals.csv, trot.csv
import sys,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'/workspace/scripts')
import numpy as np, pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import map_coordinates
from scipy.optimize import curve_fit
from q4_lib import spectrum,kb,vel,robust_sigma,REG,C
from common import f as fpath
REST=[165.050175,165.061130,165.099240]; SMU=[None]*3
JM=[165.190475,165.369341,165.678649]   # Q J=4,5,6 (rest) masked
V0=2777.; HWIN=200.; INF=1.3
REGS=['N','S','Bridge','Wclump']
cat=pd.read_csv('/workspace/out/inband_lines_curated.csv')
def line(f0): return cat[(cat.f0-f0).abs()<2e-4].iloc[0]
def cont(reg):
    p=fpath('Xb55','spw16_18_20_22',kind='cont.I.pbcor',sep='.') if False else '/data/member.uid___A001_X1273_Xb55.ngc_3256_sci.spw16_18_20_22.cont.I.pbcor.fits'
    h=fits.open(p); hd=h[0].header; im=np.squeeze(h[0].data); w=WCS(hd).celestial
    px,py=w.world_to_pixel_values(REG.loc[reg,'RA_deg'],REG.loc[reg,'Dec_deg']); return float(map_coordinates(np.nan_to_num(im),[[py],[px]],order=1)[0])
def masks(fr,v):
    m=np.zeros(len(fr),bool)
    for f0 in REST: m|=np.abs(vel(fr,f0)-V0)<HWIN+120
    for f0 in JM: m|=np.abs(vel(fr,f0)-V0)<250
    return m
def baselines(fr,y,v):
    ch=np.arange(len(fr),dtype=float); m=masks(fr,v); x=(fr-fr[len(fr)//2])/(fr[-1]-fr[0])*2
    sel=~m; out={}
    # local anchors for L/Q: within +-0.5 GHz of the blend
    loc=sel&(np.abs(fr-163.53)<0.45)
    for name,deg in (('L',1),('Q',2)):
        p=np.polyfit(x[loc],y[loc],deg); out[name]=np.polyval(p,x)
    best=None
    for tau in (15,25,40,60,90,130):
        A=np.vstack([np.ones_like(x),x,np.exp(-ch/tau)]).T; use=sel&(np.abs(fr-163.53)<0.9)
        c,*_=np.linalg.lstsq(A[use],y[use],rcond=None); r=y[use]-A[use]@c; chi=np.sum(r**2)
        if best is None or chi<best[0]: best=(chi,tau,A@c,c)
    out['R']=best[2]; out['tau']=best[1]; out['rampamp']=best[3][2]
    return out,m
def fitg(v,y):
    m=np.abs(v-2780)<700
    g=lambda v,A,v0,fw,a: A/(fw/2.3548*2.5066)*np.exp(-0.5*((v-v0)/(fw/2.3548))**2)+a
    p,cv=curve_fit(g,v[m],y[m],p0=[y.max()*300,2780,150,0]); return p,np.sqrt(abs(cv[0,0]))
def _lam(p):
    L=open(p).read().split('\n'); i=[k for k,l in enumerate(L) if 'NUMBER OF ENERGY' in l][0]; n=int(L[i+1]); lv=[l.split() for l in L[i+3:i+3+n]]
    return np.array([[float(a[1])*1.4388,float(a[2])] for a in lv]),[a[3] for a in lv]
LE,_q=_lam('/opt/lamda/e-ch3oh.dat'); LA,_=_lam('/opt/lamda/ch3oh_a.dat')
EU=float(LE[_q.index('2_0')][0]) if '2_0' in _q else 20.09
print('E-CH3OH 2_0 level energy (K, LAMDA) =',EU)
QE=lambda T: float(np.sum(LE[:,1]*np.exp(-LE[:,0]/T))); QA=lambda T: float(np.sum(LA[:,1]*np.exp(-LA[:,0]/T)))
rows=[];ramp=[];trot=[]
# thin LTE line strengths for E-type lines: W propto nu*S mu^2*exp(-Eu/T)*bose
def bose(nu,T): x=4.799e-11*nu*1e9/T; xb=4.799e-11*nu*1e9/2.73; return 1-(np.exp(x)-1)/(np.exp(xb)-1)
def wl(f0,T): l=line(f0); return f0*l.sijmu2*np.exp(-l.Eu/T)*bose(f0,T)
Ts=np.arange(4.,60.1,0.5); rng=np.random.default_rng(1)
for reg in REGS:
    fr,s,beam,nb=spectrum('Xb55','spw22',reg); v=vel(fr,REST[0]); y=s*1e3  # mJy/beam
    sig=robust_sigma(y[np.abs(y-np.median(y))<10])  # per channel
    c=cont(reg)*1e3
    B,m=baselines(fr,y,v)
    kk=kb(fr.mean(),*nb[:2])/1e3     # K per mJy/beam
    # ramp/continuum from first 30 channels (mean) vs local model baseline far
    ramp.append(dict(region=reg,cont_mJy=c,edge0_30_mean_mJy=float(np.mean(y[:30])),ramp_over_cont=float(np.mean(y[:30])/c),R_amp_fit_mJy=float(B['rampamp']),tau_ch=B['tau'],rms_ch_mJy=sig))
    dv=abs(np.median(np.diff(v))); sel=np.abs(vel(fr,REST[0])-V0)<HWIN   # window around the J=1 component; blend spans J=1-3 within ~+-0.04 GHz
    sel=(fr>163.53*(1-0.0)-0.0)&False
    fc=np.mean(REST)*(1-V0/C); sel=np.abs(C*(1-fr/fc))<HWIN+0   # centred on mean rest freq (blend components within 0.049 GHz=90 km/s)
    sel=np.abs(C*(1-fr/(np.mean(REST)*(1-V0/C))))<HWIN+45
    res={}
    for k in 'LQR':
        res[k]=np.sum((y-B[k])[sel])*dv*kk   # K km/s
    ls=np.array(list(res.values())); e_stat=sig*kk*dv*np.sqrt(sel.sum())*INF
    rows.append(dict(region=reg,I_L=res['L'],I_Q=res['Q'],I_R=res['R'],I_mean=ls.mean(),e_stat=e_stat,sys_halfrange=(ls.max()-ls.min())/2,I_blend_Jykms=ls.mean()/kk/1e3,beam='%.2fx%.2f'%nb[:2],tau_ramp_ch=B['tau']))
    # 254 line at Xb55 beam (Xb4f, spw20) and HCO+ 2-1 (Xb55 spw18) at the same position
    tb=(nb[0],nb[1],nb[2])
    f4,s4,_,nb4=spectrum('Xb4f','spw20',reg,tb=tb); v4=vel(f4,254.015377); y4=s4*kb(f4.mean(),*tb[:2])
    fh,sh,_,nbh=spectrum('Xb55','spw18',reg); vh=vel(fh,178.375056); yh=sh*kb(fh.mean(),*nbh[:2])
    ph,eh=fitg(vh,yh)
    # 254: template Gaussian at v0=ph[1], fw=ph[2], linear baseline
    mm=np.abs(v4-ph[1])<450; x=v4[mm]-ph[1]; g=np.exp(-0.5*(x/(ph[2]/2.3548))**2)/(ph[2]/2.3548*2.5066)
    M=np.vstack([g,np.ones_like(x),x/100.]).T; use=np.ones(mm.sum(),bool)
    for _ in range(4):
        p=np.linalg.lstsq(M[use],y4[mm][use],rcond=None)[0]; r=y4[mm]-M@p; far=np.abs(x)>ph[2]; sg=1.4826*np.median(np.abs(r[far&use]-np.median(r[far&use]))); use=np.abs(r)<4*sg
    cv=np.linalg.pinv(M[use].T@M[use])*sg**2; I254,e254=p[0],np.sqrt(cv[0,0])*INF
    Ib=np.array([res[k] for k in 'LQR']); Imean=Ib.mean(); sysb=(Ib.max()-Ib.min())/2
    rows[-1].update(I254=I254,e254=e254,I_HCOp21=ph[0],blend_over_HCOp21=Imean/ph[0],c254_over_HCOp21=I254/ph[0])
    # Trot posterior: baseline uniform in [min,max] of methods (widened by 1 stat), 254 Gaussian error; flat prior in T 4-60 K
    N=20000; Ibs=rng.uniform(Ib.min(),Ib.max(),N)+rng.normal(0,e_stat,N); I2=I254+rng.normal(0,e254,N)
    Rm=I2/Ibs; Rt=np.array([wl(254.015377,T)/sum(wl(f0,T) for f0 in REST) for T in Ts])
    # map each draw to T (monotonic increasing Rt); draws outside range -> piled at edges, flagged
    ok=(Rm>0)
    Tdraw=np.interp(Rm[ok],Rt[::-1],Ts[::-1]); fr_ok=1-np.mean((Rm[ok]>=Rt[0])|(Rm[ok]<=Rt[-1]))   # clamps to 4/60 K at the grid edges; frac_valid = fraction strictly inside the grid
    q=np.percentile(Tdraw,[2.5,16,50,84,97.5]); pe=np.mean(Rm[ok]<=Rt[-1])
    # 3-sigma upper limit on R if I254<3 sigma
    trot.append(dict(region=reg,R_K=I254/Imean,R_K_lo=np.percentile(Rm[ok],16),R_K_hi=np.percentile(Rm[ok],84),T16=q[1],T50=q[2],T84=q[3],T2_5=q[0],T97_5=q[4],frac_valid=fr_ok,frac_T_ge_60=pe,snr254=I254/e254,
        Tgrid_R=str(np.round(np.interp([6,10,15,20,30],Ts,Rt),3).tolist())))
    # column: N_u=8 pi k nu^2 W/(h c^3 A) per RJ thin; N(E)=N_u*Q_E(T)/g_u/exp(-Eu/T) with LAMDA-level Q ; use T=10,20,30 K assumptions
    A=10**line(254.015377).logA; nu=254.015377e9; W=I254*1e5  # K cm/s
    Nu=8*np.pi*1.380649e-16*nu**2*W/(6.62607e-27*2.998e10**3*A)
    rows[-1]['Nu_254_cm-2']=Nu
    for T in (10.,15.,20.,30.): rows[-1][f'N_E_T{int(T)}']=Nu*QE(T)/5.*np.exp(EU/T); rows[-1][f'N_AE_T{int(T)}']=rows[-1][f'N_E_T{int(T)}']*(1+QA(T)/QE(T))
rd=pd.DataFrame(rows); td=pd.DataFrame(trot); rp=pd.DataFrame(ramp)
rd.to_csv('/workspace/out/q14/blend_integrals.csv',index=False); td.to_csv('/workspace/out/q14/trot.csv',index=False); rp.to_csv('/workspace/out/q14/ramp_continuum.csv',index=False)
pd.set_option('display.width',250); print(rp.round(3).to_string()); print(rd.round(3).to_string()); print(td.round(3).to_string())
